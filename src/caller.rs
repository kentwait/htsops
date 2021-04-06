use std::collections::HashMap;
use std::fmt;

use bitflags::bitflags;
use fishers_exact::fishers_exact;

use crate::pileup::{SitePileup, SpatialSitePileup};
use crate::pileup::{FullBaseCount, AlleleCount, AlleleSet};
use crate::filter::{FilterParameters, ControlFilterScore, TargetFilterScore, ControlFilterResult, TargetFilterResult};
use crate::filter::{SNVParameters, SNVScore};


#[derive(Debug)]
pub struct PooledSNVResult {
    pub minor_allele: char,
    pub fr_ratio: f32,
    pub major_fr_ratio: f32,
    pub minor_fr_ratio: f32,
    pub full_base_count: FullBaseCount,
    pub bias_pval: f64,
}
impl fmt::Display for PooledSNVResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f, 
            "{pooled_fr_ratio:.6}\t{pooled_major_fr_ratio:.6}\t{pooled_minor_fr_ratio:.6}\t{f_a}:{r_a}:{f_c}:{r_c}:{f_g}:{r_g}:{f_t}:{r_t}\t{pooled_bias_pval:.6}",
            pooled_fr_ratio=self.fr_ratio, 
            pooled_major_fr_ratio=self.major_fr_ratio,
            pooled_minor_fr_ratio=self.minor_fr_ratio,
            f_a=self.full_base_count.f_a(),
            r_a=self.full_base_count.r_a(),
            f_c=self.full_base_count.f_c(),
            r_c=self.full_base_count.r_c(),
            f_g=self.full_base_count.f_g(),
            r_g=self.full_base_count.r_g(),
            f_t=self.full_base_count.f_t(),
            r_t=self.full_base_count.r_t(),
            pooled_bias_pval=self.bias_pval,
        )
    }
}

#[derive(Debug)]
pub struct SampleSNVResult {
    pub sample_name: String,
    pub filt_bitscore: TargetFilterScore,
    pub snv_bitscore: SNVScore,
    pub fr_ratio: f32,
    pub major_fr_ratio: f32,
    pub minor_fr_ratio: f32,
    pub full_base_count: FullBaseCount,
    pub bias_pval: f64,
}
impl fmt::Display for SampleSNVResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f, 
            "{sample_name}\t{filt_bitscore}\t{snv_bitscore}\t{fr_ratio:.6}\t{major_fr_ratio:.6}\t{minor_fr_ratio:.6}\t{f_a}:{r_a}:{f_c}:{r_c}:{f_g}:{r_g}:{f_t}:{r_t}\t{bias_pval:.6}",
            sample_name=self.sample_name,
            filt_bitscore=self.filt_bitscore.to_bits(),
            snv_bitscore=self.snv_bitscore.to_bits(),
            fr_ratio=self.fr_ratio,
            major_fr_ratio=self.major_fr_ratio,
            minor_fr_ratio=self.minor_fr_ratio,
            f_a=self.full_base_count.f_a(),
            r_a=self.full_base_count.r_a(),
            f_c=self.full_base_count.f_c(),
            r_c=self.full_base_count.r_c(),
            f_g=self.full_base_count.f_g(),
            r_g=self.full_base_count.r_g(),
            f_t=self.full_base_count.f_t(),
            r_t=self.full_base_count.r_t(),
            bias_pval=self.bias_pval,
        )
    }
}

#[derive(Debug)]
pub struct RegionwideSNVResult {
    pub pooled_result: PooledSNVResult,
    pub sample_results: Vec<SampleSNVResult>,
}
impl RegionwideSNVResult {
    pub fn minor_allele(&self) -> &char { &self.pooled_result.minor_allele }
    pub fn pooled_fr_ratio(&self) -> f32 { self.pooled_result.fr_ratio }
    pub fn pooled_major_fr_ratio(&self) -> f32 { self.pooled_result.major_fr_ratio }
    pub fn pooled_minor_fr_ratio(&self) -> f32 { self.pooled_result.minor_fr_ratio }
    pub fn pooled_full_base_count(&self) -> &FullBaseCount { &self.pooled_result.full_base_count }
    pub fn pooled_bias_pval(&self) -> f64 { self.pooled_result.bias_pval }

    pub fn passed_samples_bitstring(&self, filt_flags: TargetFilterScore, snv_flags: SNVScore) -> String {
        self.sample_results.iter()
            .map(|result| {
                let filt_result = result.filt_bitscore.contains(filt_flags);
                let snv_result = result.snv_bitscore.contains(snv_flags);
                match filt_result & snv_result {
                    true => '1',
                    false => '0',
                }
            })
            .collect::<String>()
    }
}

pub fn filter_control(pileup: &mut SitePileup, params: &FilterParameters) -> ControlFilterResult {
    let mut score: ControlFilterScore = Default::default();
    // Drop Ns and/or drop deletions
    if params.drop_n || params.drop_del {
        pileup.cleanup(params.drop_n, params.drop_del);
    }
    // Quality filter
    // Minimum coverage post filter
    let cov: usize = if let Some(c) = pileup.quality_filter(params.min_bq, params.min_mq) {
        // if c > 0 { score.insert(ControlFilterScore::PassedQualFilter) }
        if c >= params.min_cov { score.insert(ControlFilterScore::PassedMinCov) }
        c
    } else {
        0
    };
    // Check if site in control has a variant
    // Currently works only for homozygous sites
    let full_base_count = pileup.full_base_count();
    let allele_set = AlleleSet::from_basecount(&full_base_count.to_basecount());
    let total_count = allele_set.total_count() as f32;
    match allele_set.len() {
        // 0 => panic!("No alleles present. Possibly empty pileup?"),
        0 => (),
        1 => {
            score.insert(ControlFilterScore::InvariantSite);
            // score.insert(ControlFilterScore::PassedMaxVariantCount);
        },
        _ => {
            if (allele_set.len() == 2) && (allele_set.alleles[1].count() < params.max_minor_ac) {
                score.insert(ControlFilterScore::InvariantSite);
                // score.insert(ControlFilterScore::PassedMaxVariantCount);
            }
        },
    };
    // Check forward/reverse balance
    let fr_ratio = pileup.fr_ratio();
    if fr_ratio > params.fr_ratio_range.0 && fr_ratio < params.fr_ratio_range.1 {
        score.insert(ControlFilterScore::PassedFRRatio);
    }
    // Return result
    ControlFilterResult {
        score,
        cov,
        fr_ratio,
        full_base_count,
        alleles: allele_set.alleles.iter().map(|a| AlleleCount::new(a.base(), a.count())).collect(),
    }
}

pub fn quality_clean_pileups(regionwide_pileups: &mut SpatialSitePileup, params: &FilterParameters) {
    for pileup in regionwide_pileups.pileups.iter_mut() {
        if params.drop_n || params.drop_del { 
            pileup.cleanup(params.drop_n, params.drop_del);
        }
        pileup.quality_filter(params.min_bq, params.min_mq);
    }
}

pub fn call_regionwide_minor_allele(multi_pileup: &mut SpatialSitePileup, major_allele: &char) -> Option<PooledSNVResult> {
    // Pool base counts for all target samples regionwide
    let mut full_base_count: FullBaseCount = {
        let mut bc = FullBaseCount::empty();
        for p in multi_pileup.pileups.iter() {
            bc = bc + p.full_base_count();
        }
        bc
    };
    // Determine region-wide major and minor allele
    // Major is based on control allele, 
    // Must also be the top in the pooled
    let allele_set: AlleleSet = AlleleSet::from_basecount(&full_base_count.to_basecount());
    let total = allele_set.total_count();
    if allele_set.len() == 0 {
        panic!("No alleles present. Possibly empty pileup?")
    }
    if major_allele != allele_set.alleles[0].base() {
        panic!("Major allele in control [{}] and target pool [{}] do not match", 
            major_allele, allele_set.alleles[0].base());
    }
    // Minor is second highest overall
    // Major + minor should be greater than 98% of total coverage
    let minor_allele = match allele_set.len() {
        1 => None,
        2 => {
            if (allele_set.alleles[1].count() as f32 / total as f32) < 0.02 {
                None
            } else {
                Some(allele_set.alleles[1].base().to_owned())
            }
        },
        _ => {
            let total_count = allele_set.total_count() as f32;
            let major_minor_sum = (allele_set.alleles[0].count() + allele_set.alleles[1].count()) as f32;
            if (allele_set.alleles[1].count() as f32 / total as f32) < 0.02 {
                None
            } else if major_minor_sum / total_count >= 0.98 {
                Some(allele_set.alleles[1].base().to_owned())
            } else {
                None
            }
        },
    };
    if let None = minor_allele {
        return None
    }
    let minor_allele = minor_allele.unwrap();

    // Compute stats
    let total_count = full_base_count.total();
    let fr_ratio = full_base_count.forward() as f32 / total as f32;
    let (major_f, major_r) = match major_allele {
        'a' => (full_base_count.f_a(), full_base_count.r_a()),
        'c' => (full_base_count.f_c(), full_base_count.r_c()),
        'g' => (full_base_count.f_g(), full_base_count.r_g()),
        't' => (full_base_count.f_t(), full_base_count.r_t()),
        _ => panic!("")
    };
    let major_allele_count = major_f + major_r;
    let (minor_f, minor_r) = match minor_allele {
        'a' => (full_base_count.f_a(), full_base_count.r_a()),
        'c' => (full_base_count.f_c(), full_base_count.r_c()),
        'g' => (full_base_count.f_g(), full_base_count.r_g()),
        't' => (full_base_count.f_t(), full_base_count.r_t()),
        _ => panic!("")
    };
    let minor_allele_count = minor_f + minor_r;
    let major_fr_ratio = major_f as f32 / major_allele_count as f32;
    let minor_fr_ratio = minor_f as f32 / minor_allele_count as f32;
    let bias_pval = fishers_exact(&[
            major_f as u32, major_r as u32,
            minor_f as u32, minor_r as u32
        ]).unwrap().two_tail_pvalue;

    Some(PooledSNVResult {
        minor_allele,
        fr_ratio,
        major_fr_ratio,
        minor_fr_ratio,
        full_base_count,
        bias_pval,
    })
}

pub fn check_sample_snv(pileup: &mut SitePileup, major_allele: &char, minor_allele: &char, filt_params: &FilterParameters, snv_params: &SNVParameters) -> SampleSNVResult {
    let mut snv_bitscore: SNVScore = Default::default();
    let mut filt_bitscore: TargetFilterScore = Default::default();

    // Quality filter
    if pileup.cov() > 0 {
        filt_bitscore.insert(TargetFilterScore::PassedQualFilter);
    }
    if pileup.cov() >= filt_params.min_minor_ac {
        filt_bitscore.insert(TargetFilterScore::PassedMinCov);
    }
    // FR balance filter
    let fr_ratio = pileup.fr_ratio();
    if fr_ratio > filt_params.fr_ratio_range.0 && fr_ratio < filt_params.fr_ratio_range.1 {
        filt_bitscore.insert(TargetFilterScore::PassedFRRatio);
    }
    // Call if variant site or not
    let full_base_count = pileup.full_base_count();
    let allele_set = AlleleSet::from_basecount(&full_base_count.to_basecount());
    let total_count = allele_set.total_count() as f32;
    match allele_set.len() {
        0 => (),
        1 => (),
        2 => {
            filt_bitscore.insert(TargetFilterScore::VariantSite);
            filt_bitscore.insert(TargetFilterScore::PassedMaxOtherCount);
        },
        _ => {
            filt_bitscore.insert(TargetFilterScore::VariantSite);
            // Check if the total of other variants is less than the maximum
            if allele_set.alleles.iter().skip(2).map(|a| a.count()).sum::<usize>() < filt_params.max_minor_ac {
                filt_bitscore.insert(TargetFilterScore::PassedMaxOtherCount);
            }
        },
    }; 

    // Call alleles
    let (major_f, major_r) = match major_allele {
        'a' => (full_base_count.f_a(), full_base_count.r_a()),
        'c' => (full_base_count.f_c(), full_base_count.r_c()),
        'g' => (full_base_count.f_g(), full_base_count.r_g()),
        't' => (full_base_count.f_t(), full_base_count.r_t()),
        _ => panic!("")
    };
    let major_allele_count = major_f + major_r;
    let (minor_f, minor_r) = match minor_allele {
        'a' => (full_base_count.f_a(), full_base_count.r_a()),
        'c' => (full_base_count.f_c(), full_base_count.r_c()),
        'g' => (full_base_count.f_g(), full_base_count.r_g()),
        't' => (full_base_count.f_t(), full_base_count.r_t()),
        _ => panic!("")
    };
    let minor_allele_count = minor_f + minor_r;
    let minor_allele_freq = (minor_allele_count as f32) / (total_count as f32);

    // Minimum minor allele count
    if minor_allele_count >= snv_params.min_mac {
        snv_bitscore.insert(SNVScore::PassedMinMac)
    }

    // Minimum minor allele frequency
    if minor_allele_freq >= snv_params.min_maf {
        snv_bitscore.insert(SNVScore::PassedMinMaf)
    }

    let other_allele_count: usize = match (major_allele, minor_allele) {
        ('a', 'c') => full_base_count.g() + full_base_count.t(),
        ('a', 'g') => full_base_count.c() + full_base_count.t(),
        ('a', 't') => full_base_count.c() + full_base_count.g(),

        ('c', 'a') => full_base_count.g() + full_base_count.t(),
        ('c', 'g') => full_base_count.a() + full_base_count.t(),
        ('c', 't') => full_base_count.a() + full_base_count.g(),

        ('g', 'a') => full_base_count.c() + full_base_count.t(),
        ('g', 'c') => full_base_count.a() + full_base_count.t(),
        ('g', 't') => full_base_count.a() + full_base_count.c(),

        ('t', 'a') => full_base_count.c() + full_base_count.g(),
        ('t', 'c') => full_base_count.a() + full_base_count.g(),
        ('t', 'g') => full_base_count.a() + full_base_count.c(),

        _ => panic!("")
    };
    // Maximum other allele count
    if other_allele_count < snv_params.max_oac { snv_bitscore.insert(SNVScore::PassedMaxOac) }

    // Maximum other allele frequency
    if (other_allele_count as f32 / snv_params.max_oaf) < snv_params.max_oaf { snv_bitscore.insert(SNVScore::PassedMaxOaf) }

    // Within forward/reverse ratio
    // FR balance filter
    let fr_ratio = pileup.fr_ratio();
    let major_fr_ratio: f32 = major_f as f32 / major_allele_count as f32;
    let minor_fr_ratio: f32 = minor_f as f32 / minor_allele_count as f32;
    if minor_fr_ratio > snv_params.minor_fr_ratio_range.0 && minor_fr_ratio < snv_params.minor_fr_ratio_range.1 {
        snv_bitscore.insert(SNVScore::PassedMinorFRRatio);
    }
    let bias_pval = fishers_exact(&[
            major_f as u32, major_r as u32,
            minor_f as u32, minor_r as u32
        ]).unwrap().two_tail_pvalue;

    SampleSNVResult {
        sample_name: pileup.sample_name.clone(),
        filt_bitscore,
        snv_bitscore,
        fr_ratio,
        major_fr_ratio,
        minor_fr_ratio,
        full_base_count,
        bias_pval,
    }
}

pub fn check_snv_regionwide(regionwide_pileups: &mut SpatialSitePileup, control_filter_result: &ControlFilterResult, target_filt_params: &FilterParameters, target_snv_params: &SNVParameters) -> Option<RegionwideSNVResult> {
    // Determine region-wide major and minor allele
    let major_allele = control_filter_result.alleles[0].base();
    let pooled_result = match call_regionwide_minor_allele(regionwide_pileups, major_allele) {
        Some(c) => c,
        None => return None,
    };

    let minor_allele = &pooled_result.minor_allele;
    let mut sample_results: Vec<SampleSNVResult> = Vec::new();
    for pileup in regionwide_pileups.pileups.iter_mut() {
        let snv_result = check_sample_snv(pileup, major_allele, minor_allele, target_filt_params, target_snv_params);
        sample_results.push(snv_result);
    }
    Some(RegionwideSNVResult{ pooled_result, sample_results })
}
