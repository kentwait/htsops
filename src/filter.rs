use bitflags::bitflags;

// use crate::pileup::{AlleleCount, FullBaseCount};
use crate::Base;
use crate::pileup::PileupStats;


bitflags! {
    #[derive(Default)]
    pub struct ControlBitscore: usize {
        // Greater than or equal to the min coverage post quality filtering
        const PassedMinCov          = 0b00000001;
        // Ratio between forward and reverse reads are within acceptable bounds
        const PassedFRRatio         = 0b00000010;
        // Site is solely composed of a single base, ignoring singletons
        const InvariantSite         = 0b00000100;
        // Site has a different major allele compared to the hg19 reference
        const RefseqVariant         = 0b00001000;
        // Homozygous
        const Homozygous            = 0b00010000;
        // Heterozygote
        const Heterozygous          = 0b00100000;
        // Homozygous
        const IndeterminateGenotype = 0b01000000;
        // Elevated coverage depth if > 95th percentile of the bootstrap 95th percentile  
        const ElevatedCov           = 0b10000000;
    }
}
impl ControlBitscore {
    // constructor
    pub fn compute_score_from_stats(stats: &PileupStats, ref_base: Base, params: ControlScoreParams) -> Self {
        let mut bitscore = ControlBitscore::empty();
        // check if site passes minimum coverage
        bitscore.score_min_cov(stats, params.min_cov());
        // check if coverage if greater than threshold
        bitscore.score_cov_threshold(stats, params.cov_threshold());
        // check if site passes minimum forward ratio
        bitscore.score_stand_bias(stats, params.min_fratio(), params.max_fratio());
        // check if site's major allele matches the refseq base
        bitscore.score_refseq_variant(stats, ref_base);
        // check if site is invariant
        // check if hom, het or indeterminate
        bitscore.score_genotype(stats, params.max_hom_freq(), params.min_het_freq());
        bitscore
    }
    // scorers
    pub fn score_min_cov(&mut self, stats: &PileupStats, min_cov: usize) -> bool {
        if stats.cov() >= min_cov {
            self.insert(ControlBitscore::PassedMinCov);
            return true
        }
        false
    }
    pub fn score_cov_threshold(&mut self, stats: &PileupStats, cov_threshold: usize) -> bool {
        if stats.cov() > cov_threshold {
            self.insert(ControlBitscore::ElevatedCov);
            return true
        }
        false
    }
    pub fn score_stand_bias(&mut self, stats: &PileupStats, min_fratio: f64, max_fratio: f64) -> bool {
        let f_ratio = stats.base_count().forward() as f64 / (stats.cov() as f64);
        if (f_ratio >= min_fratio) && (f_ratio < max_fratio) {
            self.insert(ControlBitscore::PassedFRRatio);
            return true
        }
        false
    }
    pub fn score_refseq_variant(&mut self, stats: &PileupStats, ref_base: Base) -> bool {
        if stats.major_allele() != ref_base {
            self.insert(ControlBitscore::RefseqVariant);
            return true
        }
        false
    }
    pub fn score_genotype(&mut self, stats: &PileupStats, max_hom_freq: f64, min_het_freq: f64) -> (bool, bool, bool, bool) {
        match stats.num_alleles() {
            0 => panic!("No alleles found"),
            1 => {
                self.insert(ControlBitscore::InvariantSite);
                return (true, false, false, false)
            },
            _ => {
                // check if site is homo or heterozygous 
                let minor_freq = stats.minor_freq();
                // naive 30-70% interval filtering
                if (minor_freq >= min_het_freq) && (minor_freq < max_hom_freq) {
                    self.insert(ControlBitscore::IndeterminateGenotype);
                    return (false, false, false, true)
                } else {
                    if minor_freq >= min_het_freq {
                        self.insert(ControlBitscore::Heterozygous);
                        return (false, false, true, false)
                    }
                    self.insert(ControlBitscore::Homozygous);
                    return (false, true, false, false)
                }
            }
        }
    }

    // transformations
    pub fn to_bits(&self) -> usize { self.bits }
}
impl std::fmt::Display for ControlBitscore {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.to_bits())
    }
}

#[derive(Debug, Clone)]
pub struct ControlScoreParams {
    min_cov: usize,
    min_fratio: f64,
    max_fratio: f64,
    max_hom_freq: f64,
    min_het_freq: f64,
    cov_threshold: usize,
}
impl ControlScoreParams {
    // constructor
    pub fn new(min_cov: usize, cov_threshold: usize, min_max_fratio: (f64, f64), max_hom_freq: f64, min_het_freq: f64) -> Self {
        Self { 
            min_cov,
            min_fratio: min_max_fratio.0,
            max_fratio: min_max_fratio.1,
            max_hom_freq,
            min_het_freq,
            cov_threshold
        }
    }
    // getters
    pub fn min_cov(&self) -> usize { self.min_cov }
    pub fn min_fratio(&self) -> f64 { self.min_fratio }
    pub fn max_fratio(&self) -> f64 { self.max_fratio }
    pub fn min_het_freq(&self) -> f64 { self.min_het_freq }
    pub fn max_hom_freq(&self) -> f64 { self.max_hom_freq }
    pub fn cov_threshold(&self) -> usize { self.cov_threshold }

    // setters
    pub fn set_min_cov(&mut self, v: usize) -> Option<usize> { 
        let old_v = self.min_cov;
        self.min_cov = v;
        Some(old_v)
    }
    pub fn set_min_fratio(&mut self, v: f64) -> Option<f64> {
        let old_v = self.min_fratio;
        self.min_fratio = v;
        Some(old_v)
    }
    pub fn set_max_fratio(&mut self, v: f64) -> Option<f64> {
        let old_v = self.max_fratio;
        self.max_fratio = v;
        Some(old_v)
    }
    pub fn set_min_het_freq(&mut self, v: f64) -> Option<f64> {
        let old_v = self.min_het_freq;
        self.min_het_freq = v;
        Some(old_v)
    }
    pub fn set_max_hom_freq(&mut self, v: f64) -> Option<f64> {
        let old_v = self.max_hom_freq;
        self.max_hom_freq = v;
        Some(old_v)
    }
    pub fn set_cov_threshold(&mut self, v: usize) -> Option<usize> {
        let old_v = self.cov_threshold;
        self.cov_threshold = v;
        Some(old_v)
    }
}


// #[derive(Debug)]
// pub struct FilterParameters {
//     // Drop N and non ACGT bases
//     pub drop_n: bool,
//     // Drop deleted bases
//     pub drop_del: bool,
//     // Minimum base quality
//     pub min_bq: u8,
//     // Minimum mapping quality
//     pub min_mq: u8,
//     // Minimum coverage post quality filtering
//     pub min_cov: usize,
//     // Lower and upper bounds of acceptable forward/reverse read ratio
//     pub fr_ratio_range: (f32, f32),
//     // Minimum total count of variant bases
//     // Set to 0 if not needed (always pass)
//     pub min_minor_ac: usize,
//     // Maximum total count of variant bases
//     // Set to 10_000 if not needed (always pass)
//     pub max_minor_ac: usize,
// }

// #[derive(Debug)]
// pub struct ControlFilterResult {
//     // Bit flags indicating which filters passed
//     pub score: ControlFilterScore,
//     // Covarage post quality filtering and drops
//     pub cov: usize,
//     // Forward/reverse read ratio
//     pub fr_ratio: f32,
//     pub full_base_count: FullBaseCount,
//     // List of alleles post filtering
//     pub alleles: Vec<AlleleCount>,
// }
// impl ControlFilterResult {
//     pub fn naive_major_allele(&self) -> &char {
//         &self.alleles[0].base()
//     }
// }


// bitflags! {
//     #[derive(Default)]
//     pub struct TargetFilterScore: usize {
//         // At least one read after applying min bq and mq thresholds and drops
//         const PassedQualFilter    = 0b00000001;
//         // Greater than or equal to the min coverage post quality filtering
//         const PassedMinCov        = 0b00000010;
//         // Ratio between forward and reverse reads are within acceptable bounds
//         const PassedFRRatio       = 0b00000100;
//         // Site is solely composed of a single base
//         const VariantSite         = 0b00001000;
//         // Site has variants but less than or equal to the max allowable number 
//         const PassedMaxOtherCount = 0b00010000;
//     }
// }
// impl TargetFilterScore {
//     pub fn to_bits(&self) -> usize { self.bits }
// }

// #[derive(Debug)]
// pub struct TargetFilterResult {
//     // Bit flags indicating which filters passed
//     pub score: TargetFilterScore,
//     // Covarage post quality filtering and drops
//     pub cov: usize,
//     // Forward/reverse read ratio
//     pub fr_ratio: f32,
//     // List of alleles post filtering
//     pub alleles: Vec<AlleleCount>,
// }

// #[derive(Debug)]
// pub struct SNVParameters {
//     // Minimum minor allele freq
//     pub min_maf: f32,
//     // Minimum minor allele count
//     pub min_mac: usize,
//     // Maximum total freq of other bases
//     pub max_oaf: f32,
//     // Maximum total count of other bases
//     pub max_oac: usize,
//     // Lower and upper bounds of acceptable 
//     // forward/reverse read ratio of minor allele
//     pub minor_fr_ratio_range: (f32, f32),
// }

// bitflags! {
//     #[derive(Default)]
//     pub struct SNVScore: usize {
//         // At least one read after applying min bq and mq thresholds and drops
//         const PassedMinMaf       = 0b00000001;
//         // Greater than or equal to the min coverage post quality filtering
//         const PassedMinMac       = 0b00000010;
//         // Ratio between forward and reverse reads are within acceptable bounds
//         const PassedMaxOaf       = 0b00000100;
//         // Site is solely composed of a single base
//         const PassedMaxOac       = 0b00001000;
//         // Site has variants that is more than or equal to the minimum number 
//         const PassedMinorFRRatio = 0b00010000;
//     }
// }
// impl SNVScore {
//     pub fn to_bits(&self) -> usize { self.bits }
// }

// bitflags! {
//     #[derive(Default)]
//     pub struct SiteStatus: usize { 
//         const Heterozygous     = 0b00000001;
//         const GermlineMutation = 0b00000010;
//         const SomaticMutation  = 0b00000100;
//         const Deletion         = 0b00001000;
//         const Insertion        = 0b00010000;
//         const CopyNumber       = 0b00100000;
//         const MultipleAlleles  = 0b01000000;
//         const Indeterminate    = 0b10000000;
//     }
// }
// impl SiteStatus {
//     pub fn to_bits(&self) -> usize { self.bits }
// }

// #[derive(Debug)]
// pub struct SNVResult {
//     pub major_allele: AlleleCount,
//     pub minor_allele: Option<AlleleCount>,
//     pub cov: usize,
//     pub fr_ratio: f32,
//     pub minor_fr_ratio: f32,
//     pub snv_score: SNVScore,
//     pub status: SiteStatus,
// }
