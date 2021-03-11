use bitflags::bitflags;

use crate::pileup::{AlleleCount, FullBaseCount};

#[derive(Debug)]
pub struct FilterParameters {
    // Drop N and non ACGT bases
    pub drop_n: bool,
    // Drop deleted bases
    pub drop_del: bool,
    // Minimum base quality
    pub min_bq: u8,
    // Minimum mapping quality
    pub min_mq: u8,
    // Minimum coverage post quality filtering
    pub min_cov: usize,
    // Lower and upper bounds of acceptable forward/reverse read ratio
    pub fr_ratio_range: (f32, f32),
    // Minimum total count of variant bases
    // Set to 0 if not needed (always pass)
    pub min_minor_ac: usize,
    // Maximum total count of variant bases
    // Set to 10_000 if not needed (always pass)
    pub max_minor_ac: usize,
}


bitflags! {
    #[derive(Default)]
    pub struct ControlFilterScore: usize {
        // At least one read after applying min bq and mq thresholds and drops
        const PassedQualFilter      = 0b00000001;
        // Greater than or equal to the min coverage post quality filtering
        const PassedMinCov          = 0b00000010;
        // Ratio between forward and reverse reads are within acceptable bounds
        const PassedFRRatio         = 0b00000100;
        // Site is solely composed of a single base
        const InvariantSite         = 0b00001000;
        // Site has variants but less than or equal to the max allowable number 
        const PassedMaxVariantCount = 0b00010000;
    }
}
impl ControlFilterScore {
    pub fn to_bits(&self) -> usize { self.bits }
}

#[derive(Debug)]
pub struct ControlFilterResult {
    // Bit flags indicating which filters passed
    pub score: ControlFilterScore,
    // Covarage post quality filtering and drops
    pub cov: usize,
    // Forward/reverse read ratio
    pub fr_ratio: f32,
    pub full_base_count: FullBaseCount,
    // List of alleles post filtering
    pub alleles: Vec<AlleleCount>,
}
impl ControlFilterResult {
    pub fn naive_major_allele(&self) -> &char {
        &self.alleles[0].base()
    }
}


bitflags! {
    #[derive(Default)]
    pub struct TargetFilterScore: usize {
        // At least one read after applying min bq and mq thresholds and drops
        const PassedQualFilter    = 0b00000001;
        // Greater than or equal to the min coverage post quality filtering
        const PassedMinCov        = 0b00000010;
        // Ratio between forward and reverse reads are within acceptable bounds
        const PassedFRRatio       = 0b00000100;
        // Site is solely composed of a single base
        const VariantSite         = 0b00001000;
        // Site has variants but less than or equal to the max allowable number 
        const PassedMaxOtherCount = 0b00010000;
    }
}
impl TargetFilterScore {
    pub fn to_bits(&self) -> usize { self.bits }
}

#[derive(Debug)]
pub struct TargetFilterResult {
    // Bit flags indicating which filters passed
    pub score: TargetFilterScore,
    // Covarage post quality filtering and drops
    pub cov: usize,
    // Forward/reverse read ratio
    pub fr_ratio: f32,
    // List of alleles post filtering
    pub alleles: Vec<AlleleCount>,
}

#[derive(Debug)]
pub struct SNVParameters {
    // Minimum minor allele freq
    pub min_maf: f32,
    // Minimum minor allele count
    pub min_mac: usize,
    // Maximum total freq of other bases
    pub max_oaf: f32,
    // Maximum total count of other bases
    pub max_oac: usize,
    // Lower and upper bounds of acceptable 
    // forward/reverse read ratio of minor allele
    pub minor_fr_ratio_range: (f32, f32),
}

bitflags! {
    #[derive(Default)]
    pub struct SNVScore: usize {
        // At least one read after applying min bq and mq thresholds and drops
        const PassedMinMaf       = 0b00000001;
        // Greater than or equal to the min coverage post quality filtering
        const PassedMinMac       = 0b00000010;
        // Ratio between forward and reverse reads are within acceptable bounds
        const PassedMaxOaf       = 0b00000100;
        // Site is solely composed of a single base
        const PassedMaxOac       = 0b00001000;
        // Site has variants that is more than or equal to the minimum number 
        const PassedMinorFRRatio = 0b00010000;
    }
}
impl SNVScore {
    pub fn to_bits(&self) -> usize { self.bits }
}

bitflags! {
    #[derive(Default)]
    pub struct SiteStatus: usize { 
        const Heterozygous     = 0b00000001;
        const GermlineMutation = 0b00000010;
        const SomaticMutation  = 0b00000100;
        const Deletion         = 0b00001000;
        const Insertion        = 0b00010000;
        const CopyNumber       = 0b00100000;
        const MultipleAlleles  = 0b01000000;
        const Indeterminate    = 0b10000000;
    }
}
impl SiteStatus {
    pub fn to_bits(&self) -> usize { self.bits }
}

#[derive(Debug)]
pub struct SNVResult {
    pub major_allele: AlleleCount,
    pub minor_allele: Option<AlleleCount>,
    pub cov: usize,
    pub fr_ratio: f32,
    pub minor_fr_ratio: f32,
    pub snv_score: SNVScore,
    pub status: SiteStatus,
}
