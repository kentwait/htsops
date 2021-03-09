#[macro_use]
use std::path::Path;
use std::collections::HashMap;
use rust_htslib::tbx::{self, Read as TbxRead};

use htsops::pileup::{SitePileup, SpatialSitePileup};
use htsops::pileup::{BaseCount, FullBaseCount, AlleleCount, AlleleFreq, AlleleSet};
use htsops::filter::{FilterParameters, ControlFilterScore, TargetFilterScore, ControlFilterResult, TargetFilterResult};
use htsops::filter::{SNVParameters, SNVScore, SiteStatus, SNVResult};


// # Pass 1 columns
// chr
// pos
// major_char
// minor_char
// control_filt_bitscore
// control_fr_ratio
// control_basecount A:a:C:c:G:g:T:t
// pooled_fr_ratio
// pooled_major_fr_ratio
// pooled_minor_fr_ratio
// pooled_basecount A:a:C:c:G:g:T:t
// pooled_bias_fisher
// sample_name
// sample_filt_bitscore
// sample_snv_bitscore
// sample_fr_ratio
// sample_major_fr_ratio
// sample_minor_fr_ratio
// sample_basecount A:a:C:c:G:g:T:t
// sample_bias_fisher
// passed_samples_bitscore
// # Header
// colum names
// control_filt_score guide
// sample_filt_score guide
// sample_snv_score guide
// control_filt params
// sample_filt params
// sample_snv params
// # Footer
// for each sample: for each chromosome: coverage depth mean, variance, 95%CI
// for each sample: binned coverage depth 0 to 100 percentile, in 10 percentile increments
// for each sample: binned coverage depth vs fr_ratio
// for each sample: binned coverage depth vs ma_fr_ratio
// for each sample: binned coverage depth vs sample_bias_fisher
// major->minor vs count and freq
// passed_samples_bitscore vs count and frequency

// # Pass 2 columns
// chr
// pos
// major_char
// minor_char
// control_filt_bitscore
// control_fr_ratio
// control_basecount A:a:C:c:G:g:T:t
// pooled_fr_ratio
// pooled_major_fr_ratio
// pooled_minor_fr_ratio
// pooled_basecount A:a:C:c:G:g:T:t
// pooled_bias_fisher
// sample_name
// sample_filt_bitscore
// sample_snv_bitscore
// sample_fr_ratio
// sample_major_fr_ratio
// sample_minor_fr_ratio
// sample_basecount A:a:C:c:G:g:T:t
// sample_bias_fisher
// passed_samples_bitscore
// passed_site_bitscores


fn filter_control(pileup: &mut SitePileup, params: &FilterParameters) -> ControlFilterResult {
    let mut score: ControlFilterScore = Default::default();
    // Drop Ns and/or drop deletions
    if params.drop_n || params.drop_del {
        pileup.cleanup(params.drop_n, params.drop_del);
    }
    // Quality filter
    // Minimum coverage post filter
    let cov: usize = if let Some(c) = pileup.quality_filter(params.min_bq, params.min_mq) {
        if c > 0 { score.insert(ControlFilterScore::PassedQualFilter) }
        if c >= params.min_cov { score.insert(ControlFilterScore::PassedMinCov) }
        c
    } else {
        0
    };
    // Check if site in control has a variant
    // Currently works only for homozygous sites
    let allele_set = AlleleSet::from_basecount(pileup.base_count());
    let total_count = allele_set.total_count() as f32;
    match allele_set.len() {
        // 0 => panic!("No alleles present. Possibly empty pileup?"),
        0 => (),
        1 => {
            score.insert(ControlFilterScore::InvariantSite);
            score.insert(ControlFilterScore::PassedMaxVariantCount);
        },
        _ => {
            if (allele_set.len() == 2) && (allele_set.alleles[1].count() < params.max_minor_ac) {
                score.insert(ControlFilterScore::InvariantSite);
                score.insert(ControlFilterScore::PassedMaxVariantCount);
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
        alleles: allele_set.alleles.iter().map(|a| AlleleCount::new(a.base(), a.count())).collect(),
    }
}

fn quality_clean_pileup(pileup: &mut SitePileup, params: &FilterParameters) {
    if params.drop_n || params.drop_del { 
        pileup.cleanup(params.drop_n, params.drop_del);
    }
    pileup.quality_filter(params.min_bq, params.min_mq);
}

enum MinorAlleleCall {
    Some(char),
    None,
    Indeterminate,
}

fn call_regionwide_minor_allele(multi_pileup: &mut SpatialSitePileup, sample_names: &Vec<&str>, major_allele: &char) -> MinorAlleleCall {
    let mut base_count: BaseCount = BaseCount::empty();
    for sample_name in sample_names.iter() {
        let pileup = multi_pileup.pileups.get_mut(*sample_name).unwrap();
        base_count = base_count + pileup.base_count();
    }
    // Determine region-wide major and minor allele
    // Major is based on control allele, 
    // this must also be the top in the pooled
    // Minor is second highest overall
    // Major + minor should be greater than 98% of total coverage
    let allele_set: AlleleSet = AlleleSet::from_basecount(base_count);
    let total = allele_set.total_count();
    if allele_set.len() == 0 {
        panic!("No alleles present. Possibly empty pileup?")
    }
    if major_allele != allele_set.alleles[0].base() {
        panic!("Major allele in control [{}] and target pool [{}] do not match", 
            major_allele, allele_set.alleles[0].base());
    }
    match allele_set.len() {
        1 => MinorAlleleCall::None,
        2 => {
            if (allele_set.alleles[1].count() as f32 / total as f32) < 0.02 {
                MinorAlleleCall::None
            } else {
                MinorAlleleCall::Some(allele_set.alleles[1].base().to_owned())
            }
        },
        _ => {
            let total_count = allele_set.total_count() as f32;
            let major_minor_sum = (allele_set.alleles[0].count() + allele_set.alleles[1].count()) as f32;
            if (allele_set.alleles[1].count() as f32 / total as f32) < 0.02 {
                MinorAlleleCall::Indeterminate
            } else if major_minor_sum / total_count >= 0.98 {
                MinorAlleleCall::Some(allele_set.alleles[1].base().to_owned())
            } else {
                MinorAlleleCall::Indeterminate
            }
        },
    }
}

fn filter_target(pileup: &mut SitePileup, params: &FilterParameters) -> TargetFilterScore {
    let mut score: TargetFilterScore = Default::default();
    // Quality filter
    if pileup.cov() > 0 {
        score.insert(TargetFilterScore::PassedQualFilter);
    }
    if pileup.cov() >= params.min_minor_ac {
        score.insert(TargetFilterScore::PassedMinCov);
    }
    // FR balance filter
    let fr_ratio = pileup.fr_ratio();
    if fr_ratio > params.fr_ratio_range.0 && fr_ratio < params.fr_ratio_range.1 {
        score.insert(TargetFilterScore::PassedFRRatio);
    }
    // Call alleles
    let allele_set = AlleleSet::from_basecount(pileup.base_count());
    let total_count = allele_set.total_count() as f32;
    match allele_set.len() {
        0 => (),
        1 => (),
        2 => {
            score.insert(TargetFilterScore::VariantSite);
            score.insert(TargetFilterScore::PassedMaxOtherCount);
        },
        _ => {
            score.insert(TargetFilterScore::VariantSite);
            // Check if the total of other variants is less than the maximum
            if allele_set.alleles.iter().skip(2).map(|a| a.count()).sum::<usize>() < params.max_minor_ac {
                score.insert(TargetFilterScore::PassedMaxOtherCount);
            }
        },
    }; 
    score
}

fn call_target_snv(pileup: &mut SitePileup, params: &SNVParameters) -> Option<SNVResult> {
    let mut score: SNVScore = Default::default();

    // Call alleles
    let full_base_count = pileup.full_base_count();
    let total = full_base_count.total();
    let major_allele_count: usize = match params.major_allele {
        'a' => full_base_count.a(),
        'c' => full_base_count.c(),
        'g' => full_base_count.g(),
        't' => full_base_count.t(),
        _ => panic!("")
    };
    let minor_allele_count: usize = match params.minor_allele {
        'a' => full_base_count.a(),
        'c' => full_base_count.c(),
        'g' => full_base_count.g(),
        't' => full_base_count.t(),
        _ => panic!("")
    };
    let minor_allele_freq = (minor_allele_count as f32) / (total  as f32);

    // Minimum minor allele count
    if minor_allele_count >= params.min_mac { score.insert(SNVScore::PassedMinMac) }

    // Minimum minor allele frequency
    if minor_allele_freq >= params.min_maf { score.insert(SNVScore::PassedMinMaf) }

    let other_allele_count: usize = match (params.major_allele, params.minor_allele) {
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
    if other_allele_count < params.max_oac { score.insert(SNVScore::PassedMaxOac) }

    // Maximum other allele frequency
    if (other_allele_count as f32 / params.max_oaf) < params.max_oaf { score.insert(SNVScore::PassedMaxOaf) }

    // Within forward/reverse ratio
    // FR balance filter
    let fr_ratio = pileup.fr_ratio();
    let minor_fr_ratio: f32 = match params.major_allele {
        'a' => (full_base_count.f_a() as f32) / (full_base_count.a() as f32),
        'c' => (full_base_count.f_c() as f32) / (full_base_count.c() as f32),
        'g' => (full_base_count.f_g() as f32) / (full_base_count.g() as f32),
        't' => (full_base_count.f_t() as f32) / (full_base_count.t() as f32),
        _ => panic!("")
    };
    if minor_fr_ratio > params.minor_fr_ratio_range.0 && minor_fr_ratio < params.minor_fr_ratio_range.1 {
        score.insert(SNVScore::PassedMinorFRRatio);
    }

    // If all passed then call
    if score.is_all() {
        return Some(SNVResult{
            major_allele: AlleleCount::new(&params.major_allele, major_allele_count),
            minor_allele: Some(AlleleCount::new(&params.minor_allele, minor_allele_count)),
            cov: total,
            fr_ratio,
            minor_fr_ratio: minor_fr_ratio,
            snv_score: score,
            status: SiteStatus::SomaticMutation,
        })
    }

    None
}

fn main() {
    let pileup_db_path = Path::new("/Users/kent/Desktop/UPN59.pileup.sled");
    let mpileup_path = Path::new("/Volumes/NVME4TB/local_data/yokoyama/UPN59.pileup.bq0.mq0.adjmq50.bgz");
    let sample_list: Vec<&str> = vec![
        "UPN59_B",
        "UPN59_1_1",
        "UPN59_1_3",
        "UPN59_1_4",
        "UPN59_1_5",
        "UPN59_1_6",
        "UPN59_1_7",
        "UPN59_1_8",
        "UPN59_1_9",
        "UPN59_1_10",
        "UPN59_1_11",
        "UPN59_1_12",
        "UPN59_1_13",
        "UPN59_1_14",
        "UPN59_1_15",
        "UPN59_1_16",
        "UPN59_1_17",
        "UPN59_1_18",
        "UPN59_1_19",
        "UPN59_1_21",
        "UPN59_1_22",
        "UPN59_1_23",
        "UPN59_1_24",
        "UPN59_1_25",
    ];
    let control_name: &str = sample_list[0];

    // Program settings
    let NUM_THREADS = 4;
    let FETCH_CHUNKSIZE: u64 = 1_000_000;
    // Filtering parameters
    let MIN_BQ = 20;
    let MIN_MQ = 40;
    let MIN_COV = 20;
    let MAX_CONTROL_MINAC = 1;
    let MIN_SAMPLE_AC = 3;
    let (FR_RATIO_LB, FR_RATIO_UB) = (0.3, 0.7);

    // Declare filters
    let control_filter_params = FilterParameters {
        drop_n: true,
        drop_del: true,
        min_bq: 20,
        min_mq: 40,
        min_cov: 20,
        fr_ratio_range: (0.3, 0.7),
        min_minor_ac: 0,
        max_minor_ac: 2,
    };
    let target_filter_params = FilterParameters {
        drop_n: true,
        drop_del: true,
        min_bq: 20,
        min_mq: 40,
        min_cov: 20,
        fr_ratio_range: (0.3, 0.7),
        min_minor_ac: 0,
        max_minor_ac: 1,
    };

    // DB configuration
    // Store pileups hierarchically
    // key is chromosome, value is a BtreeMap whose key is (start, end) tuple
    // value Vec<PileupRecord>
    // use sled embedded db
    // let db: sled::Db = sled::Config::new()
    //     .path(pileup_db_path)
    //     .cache_capacity(10_000_000_000)
    //     .flush_every_ms(Some(1000))
    //     // .create_new(true)
    //     .open()
    //     .unwrap();

    // Open tabix-indexed bgzipped mpileup file
    let mut tbx_reader = tbx::Reader::from_path(&mpileup_path)
        .expect(&format!("Could not open {}", mpileup_path.display()));
    tbx_reader.set_threads(NUM_THREADS);

    // Resolve chromosome name to numeric ID.
    let target_chromosomes = vec![
        "chr1", "chr2", "chr3", "chr4", "chr5", 
        "chr6", "chr7", "chr8", "chr9", "chr10",
        "chr11", "chr12", "chr13", "chr14", "chr15", 
        "chr16", "chr17", "chr18", "chr19", "chr20",
        "chr21", "chr22", "chrX", "chrY"];
    let chrom_size: HashMap<&str, usize> = vec![
            ("chr1", 249250621),
            ("chr2", 243199373),
            ("chr3", 198022430),
            ("chr4", 191154276),
            ("chr5", 180915260), 
            ("chr6", 171115067),
            ("chr7", 159138663),
            ("chrX", 155270560),
            ("chr8", 146364022),
            ("chr9", 141213431),
            ("chr10", 135534747),
            ("chr11", 135006516),
            ("chr12", 133851895),
            ("chr13", 115169878),
            ("chr14", 107349540),
            ("chr15", 102531392), 
            ("chr16", 90354753),
            ("chr17", 81195210),
            ("chr18", 78077248),
            ("chr20", 63025520),
            ("chrY", 59373566),
            ("chr19", 59128983),
            ("chr22", 51304566),
            ("chr21", 48129895),
        ].into_iter().collect();
    let chrom_tid_lookup: HashMap<&str, u64> = target_chromosomes.iter().map(|c| {
            match tbx_reader.tid(c) {
                Ok(tid) => (c.to_owned(), tid),
                Err(_) => panic!("Could not resolve contig ID"),
            }
        }).collect();

    // Read through all records in all chromosomes
    // Look at control first and apply filter
    // If control passed, set major allele
    // Then pool targets and apply filter as a pool
    // If pool passed, call minor allele
    // Then, apply filter to each sample and call SNV independently
    // Return SiteStatus per site per sample

    // FIRST PASS
    // Read through all records in all chromosomes
    // Compute statistics:
    // - cov
    // - fr ratio
    // - median 
    for chrom in target_chromosomes.iter() {
        let ul: u64 = (chrom_size.get(chrom).unwrap().to_owned() as u64 / FETCH_CHUNKSIZE) + 1;
        for i in 0..ul {
            tbx_reader.fetch(*chrom_tid_lookup.get(chrom).unwrap(), i*FETCH_CHUNKSIZE, (i+1)*FETCH_CHUNKSIZE).unwrap();
            for record in tbx_reader.records() {
                // Convert binary mpileup line and parse
                let record_string = String::from_utf8(record.unwrap()).unwrap();
                let mut multi_pileups = SpatialSitePileup::parse_mpileup_row(
                    &record_string, &sample_list);

                // Process control pileup first
                let control_pileup = multi_pileups.pileups.get_mut(control_name).unwrap();
                let control_filter_result = filter_control(control_pileup, &control_filter_params);

                // println!("{}\t{}\t{:?}", multi_pileups.chrom, multi_pileups.pos, control_filter_result.score);

                if !control_filter_result.score.is_all() {
                    continue
                }
                // Set major allele from control
                let major_allele = control_filter_result.alleles[0].base();

                // Process non-control pileups at the same site
                for (sample_name, pileup) in multi_pileups.pileups.iter_mut() {
                    if sample_name == control_name { continue }
                    quality_clean_pileup(pileup, &control_filter_params);
                }
                // Determine region-wide major and minor allele
                let minor_allele_call: MinorAlleleCall = call_regionwide_minor_allele(&mut multi_pileups, &sample_list, major_allele);

                for (sample_name, pileup) in multi_pileups.pileups.iter_mut() {
                    if sample_name == control_name { continue }
                    if let MinorAlleleCall::Some(b) = minor_allele_call {
                        let target_score = filter_target(pileup, &target_filter_params);
                        if target_score.is_all() {
                            let snv_params = SNVParameters {
                                major_allele: major_allele.to_owned(),
                                minor_allele: b,
                                min_maf: 0.01,
                                min_mac: 1,
                                max_oaf: 0.01,
                                max_oac: 1,
                                minor_fr_ratio_range: (0.3, 0.7),
                            };                        
                            if let Some(result) = call_target_snv(pileup, &snv_params) {
                                // println!("\t{:?}", result)
                                println!("{}:{}\t{}\t{:?}", 
                                    multi_pileups.chrom, multi_pileups.pos,
                                    sample_name, result)
                            }
                        }
                    }
                }
            }
        }
    }
}
