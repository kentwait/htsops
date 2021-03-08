#[macro_use]
use bitflags;


pub struct FilterParameters {
    // Drop N and non ACGT bases
    pub drop_n: bool
    // Drop deleted bases
    pub drop_del: bool
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
    struct FilterScore: u32 {
        // At least one read after applying min bq and mq thresholds and drops
        const PassedQualFilter    = 0b00000001;
        // Greater than or equal to the min coverage post quality filtering
        const PassedMinCov        = 0b00000010;
        // Ratio between forward and reverse reads are within acceptable bounds
        const PassedFRRatio       = 0b00000100;
        // Site is solely composed of a single base
        const InvariantSite       = 0b00001000;
        // Site has variants that is more than or equal to the minimum number 
        const PassedMinMinorCount = 0b00010000;
        // Site has variants but less than or equal to the max allowable number 
        const PassedMaxMinorCount = 0b00100000;
    }
}


pub struct FilterResult {
    // Bit flags indicating which filters passed
    pub score: FilterScore,
    // Covarage post quality filtering and drops
    pub cov: usize,
    // Forward/reverse read ratio
    pub fr_ratio: f32,
    // List of alleles post filtering
    pub alleles: Vec<AlleleFreq>,
}


pub struct SNVParameters {
    // Major allele based on control
    pub major_allele: char,
    // Minor allele based on pooling
    pub minor_allele: char,
    // Minimum minor allele freq
    pub min_maf: f64,
    // Minimum minor allele count
    pub min_mac: usize,
    // Maximum total freq of other bases
    pub max_oaf: f64,
    // Maximum total count of other bases
    pub min_oac: usize,
    // Lower and upper bounds of acceptable 
    // forward/reverse read ratio of minor allele
    pub minor_fr_ratio_range: (f32, f32),
}


pub enum SNVCall{
    NoMutation,
    GermlineMutation,
    SomaticMutation,
    Deletion,
    CopyNumberVariation,
    Indeterminate,
}

bitflags! {
    struct SiteStatus: u32 { 
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


pub struct SNVResult {
    pub major_allele: AlleleFreq,
    pub minor_allele: Option<AlleleFreq>,
    pub cov: usize,
    pub fr_ratio: f32,
    pub minor_fr_ratio: f32,
    pub status: SiteStatus,
}


pub struct AlleleCount(char, usize);
impl AlleleCount {
    pub fn base(&self) -> &char { &self.0 }
    pub fn base_index(&self) -> usize { 
        match self.0.to_ascii_lowercase() {
            'a' => 0,
            'c' => 1,
            'g' => 2,
            't' => 3,
            'n' => 4,
            b => panic!("Invalid base [{}]", b),
        }    
    }
    pub fn count(&self) -> usize { self.1}
}


pub struct AlleleFreq(char, f64);
impl AlleleFreq {
    pub fn base(&self) -> &char { &self.0 }
    pub fn freq(&self) -> f64 { self.1}
}