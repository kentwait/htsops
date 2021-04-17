use std::collections::HashMap;

use indexmap::IndexMap;
use lazy_static::lazy_static;


pub const FETCH_CHUNKSIZE: u64 = 1_000_000;
pub const FREQ_DIVISOR: f64 = 100.0;

lazy_static! {
    /// hg19 chromosome list
    pub static ref HG19_CHROMS: Vec<&'static str> = vec![
        "chr1", "chr2", "chr3", "chr4", "chr5", 
        "chr6", "chr7", "chr8", "chr9", "chr10",
        "chr11", "chr12", "chr13", "chr14", "chr15", 
        "chr16", "chr17", "chr18", "chr19", "chr20",
        "chr21", "chr22", "chrX", "chrY"];
    /// hg19 chromosome lengths
    pub static ref HG19_CHROM_LEN: HashMap<&'static str, u64> = vec![
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
    ]
    .into_iter()
    .collect();
    /// hg19 chromosome cummulative length
    pub static ref HG19_CHROM_CUMM_LEN: IndexMap<&'static str, u64> = vec![
        ("chr1", 249250621),
        ("chr2", 492449994),
        ("chr3", 690472424),
        ("chr4", 881626700),
        ("chr5", 1062541960),
        ("chr6", 1233657027),
        ("chr7", 1392795690),
        ("chr8", 1539159712),
        ("chr9", 1680373143),
        ("chr10", 1815907890),
        ("chr11", 1950914406),
        ("chr12", 2084766301),
        ("chr13", 2199936179),
        ("chr14", 2307285719),
        ("chr15", 2409817111),
        ("chr16", 2500171864),
        ("chr17", 2581367074),
        ("chr18", 2659444322),
        ("chr19", 2718573305),
        ("chr20", 2781598825),
        ("chr21", 2829728720),
        ("chr22", 2881033286),
        ("chrX", 3036303846),
        ("chrY", 3095677412),
    ]
    .into_iter()
    .collect();

}