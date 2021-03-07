use std::path::Path;
use std::collections::HashMap;
use rust_htslib::tbx::{self, Read as TbxRead};

use htsops::pileup::{SitePileup, SpatialSitePileup, ControlFilterResult};

pub struct FilterParameters {
    pub min_bq: u8,
    pub min_mq: u8,
    pub min_cov: usize,
    pub min_ac: usize,
    pub fr_ratio_range: (f32, f32),
    pub drop_n: bool
}

pub struct CallSNVParameters {
    pub min_bq: u8,
    pub min_mq: u8,
    pub min_cov: usize,
    pub min_ac: usize,
    pub fr_ratio_range: (f32, f32),
    pub minor_fr_ratio_range: (f32, f32),
    pub drop_n: bool
}

pub struct CallResult {
    pub score: u8,
    pub fr_ratio: f32,
    pub cov: usize,
    pub major_allele: Allele,
    pub minor_allele: Option<Allele>,
}

pub struct Allele(char, usize);
impl Allele {
    pub fn base(&self) -> &char { &self.0 }
    pub fn count(&self) -> usize { self.1}
}

fn filter_normal(pileup: &mut SitePileup, params: &FilterParameters) -> CallResult {
    let mut score = 0;
    // Quality filter
    let cov = match pileup.quality_filter(params.min_bq, params.min_mq, params.drop_n) {
        Some(c) if c >= params.min_cov => { score += 3; c },
        Some(c) if c < params.min_cov => { score += 1; c },
        _ => 0,
    };
    // Major allele and count
    let (major_allele, minor_allele) = {
        let alleles = pileup.allele_count();
        // TODO: Tier perfect and imperfect
        // Doesnt consider hetero
        match alleles.len() {
            0 => panic!("No alleles present. Possibly empty pileup?"),
            1 => {
                score += 4;
                (Allele(alleles[0].0, alleles[0].1), None)
            },
            2 => if alleles[1].1 == params.min_ac {
                    score += 4;
                    (Allele(alleles[0].0, alleles[0].1), Some(Allele(alleles[1].0, alleles[1].1)))
                 } else {
                    (Allele(alleles[0].0, alleles[0].1), None)
                 }
            _ => (Allele(alleles[0].0, alleles[0].1), None),
        }
    };
    // FR balance filter
    let fr_ratio = pileup.fr_ratio();
    if fr_ratio > params.fr_ratio_range.0 && fr_ratio < params.fr_ratio_range.1 {
        score += 8;
    }
    CallResult{
        score,
        fr_ratio: 0.,
        cov: 0,
        major_allele,
        minor_allele,
    }
}

fn call_regionwide_minor_alleles(multi_pileup: &mut SpatialSitePileup, sample_names: &Vec<&str>) -> Option<Allele> {
    let (mut tot_a, mut tot_c, mut tot_g, mut tot_t) = (0, 0, 0, 0);
    for sample_name in sample_names.iter() {
        let mut pileup = multi_pileup.pileups.get_mut(*sample_name).unwrap();
        let (_a, _c, _g, _t, _) = pileup.base_count();
        tot_a += _a;
        tot_c += _c;
        tot_g += _g;
        tot_t += _t;
    }
    // Determine region-wide major and minor allele
    let mut alleles: Vec<Allele> = vec![
        ('a', tot_a), ('c', tot_c), ('g', tot_g), ('t', tot_t)]
    .into_iter()
    .filter_map(|(k, v)| match v {
        0 => None,
        _ => Some(Allele(k, v))
    }).collect();
    alleles.sort_by(|Allele(_, a), Allele(_, b)| b.cmp(a));

    match alleles.len() {
        2 => {
            Some(alleles[1])
        },
        3 => {
            if alleles[1].count() > 1 && alleles[2].count() == 1 {
                Some(alleles[1])
            } else {
                None
            }
        },
        4 => {
            if alleles[1].count() > 1 && alleles[3].count() == 1 {
                Some(alleles[1])
            } else {
                None
            }
        },
        _ => None,
    }
}

fn call_snv_on_target(pileup: &mut SitePileup, params: &CallSNVParameters) -> CallResult {
    let mut score = 0;
    // Quality filter
    let cov = match pileup.quality_filter(params.min_bq, params.min_mq, params.drop_n) {
        Some(c) if c >= params.min_cov => { score += 3; c },
        Some(c) if c < params.min_cov => { score += 1; c },
        _ => 0,
    };
    // Major allele and count
    let (major_allele, major_cnt, minor_allele, minor_cnt) = {
        let alleles = pileup.allele_count();
        // TODO: Tier perfect and imperfect
        // Doesnt consider hetero
        match alleles.len() {
            0 => panic!("No alleles present. Possibly empty pileup?"),
            1 => {
                score += 4;
                (alleles[0].0, alleles[0].1, None, None)
            },
            2 => if alleles[1].1 == params.min_ac {
                    score += 4;
                    (alleles[0].0, alleles[0].1, Some(alleles[1].0), Some(alleles[1].1))
                 } else {
                    (alleles[0].0, alleles[0].1, None, None)
                 }
            _ => (alleles[0].0, alleles[0].1, None, None),
        }
    };
    // FR balance filter
    let fr_ratio = pileup.fr_ratio();
    if fr_ratio > params.fr_ratio_lb && fr_ratio < params.fr_ratio_ub {
        score += 8;
    }
    CallResult{
        score,
        fr_ratio: 0.,
        cov: 0,
        major_allele,
        minor_allele,
    }
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
    let control_filter_params = FilterParameters {
        min_bq: 20,
        min_mq: 40,
        min_cov: 20,
        min_ac: 1,
        fr_ratio_lb: 0.3,
        fr_ratio_ub: 0.7,
        drop_n: true,
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
                Err(_) => panic!("Could not resolve 'chr1' to contig ID"),
            }
        }).collect();

    
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
                let control_filter_result = filter_normal(control_pileup, &control_filter_params);

                if control_filter_result.score < 15 {
                    continue
                }

                // Process non-control pileups at the same site
                let (mut tot_a, mut tot_c, mut tot_g, mut tot_t) = (0, 0, 0, 0);
                let mut sample_base_counts: HashMap<&str, [usize; 4]> = HashMap::with_capacity(multi_pileups.pileups.len());
                let mut sample_filt_covs: HashMap<&str, usize> = HashMap::with_capacity(multi_pileups.pileups.len());
                let mut sample_fr_ratio: HashMap<&str, f64> = HashMap::with_capacity(multi_pileups.pileups.len());
                for (sample_name, pileup) in multi_pileups.pileups.iter_mut() {
                    if sample_name == control_name { continue }
                    if let Some(filt_cov) = pileup.quality_filter(MIN_BQ, MIN_MQ, true) {
                        // Add to total
                        let (_a, _c, _g, _t, _) = pileup.base_count();
                        tot_a += _a;
                        tot_c += _c;
                        tot_g += _g;
                        tot_t += _t;
                        sample_base_counts.insert(sample_name, [_a, _c, _g, _t]);
                        sample_filt_covs.insert(sample_name, filt_cov);
                        sample_fr_ratio.insert(sample_name, pileup.fr_ratio());
                    };
                }
                // Determine region-wide major and minor allele
                let mut alleles: Vec<(char, usize)> = vec![
                        ('a', tot_a), ('c', tot_c), ('g', tot_g), ('t', tot_t)]
                    .into_iter()
                    .filter_map(|(k, v)| match v {
                        0 => None,
                        _ => Some((k, v))
                    }).collect();
                alleles.sort_by(|(_, a), (_, b)| b.cmp(a));
                if alleles.len() < 2 { continue }
                // Make sure that major+minor make up more than 95% of depth
                let major_minor_percent = (alleles[0].1 + alleles[1].1) as f64 / ((tot_a + tot_c + tot_g + tot_t) as f64);
                if major_minor_percent < 0.95 { continue }
                if alleles[1].1 < MIN_SAMPLE_AC { continue }
                let (major, minor) = (alleles[0].0, alleles[1].0);
                let major_idx = match alleles[0].0 {
                    'a' => 0,
                    'c' => 1,
                    'g' => 2,
                    't' => 3,
                    _ => 4,
                };
                let minor_idx = match alleles[1].0 {
                    'a' => 0,
                    'c' => 1,
                    'g' => 2,
                    't' => 3,
                    _ => 4,
                };
                if control_major != major { continue }

                print!("{chrom}\t{pos}\t{cov}:{major}:{minor}:{fr}:{score}:{passed}", 
                    chrom=chrom,
                    pos=pos,
                    cov=cov,
                    major=major,
                    minor=minor,
                    fr=fr_ratio,
                    score=filter_score,
                    passed=passed
                );
                // Call SNV
                for sample_name in sample_list.iter().skip(1) {
                    let base_count = sample_base_counts.get(sample_name).unwrap();
                    let major_cnt = base_count[major_idx];
                    let minor_cnt = base_count[minor_idx];

                    let mut filter_score: usize = 0;
                    if minor_cnt > 0 {
                        filter_score = 1;
                        // Quality filter
                        let filt_cov = sample_filt_covs.get(sample_name).unwrap();
                        if *filt_cov > 0 {
                            filter_score += 2;
                            // Coverage filter
                            if *filt_cov >= MIN_COV { 
                                filter_score += 4;
                                // Min count
                                // if minor_cnt >= MIN_SAMPLE_AC {
                                //     filter_score += 4;
                                // }
                            }
                        }
                        // FR balance filter
                        if filter_score == 7 {
                            let fr_ratio = sample_fr_ratio.get(sample_name).unwrap();
                            if *fr_ratio > FR_RATIO_LB && *fr_ratio < FR_RATIO_UB {
                                filter_score += 8;
                            }
                        }
                    }
                    let passed = match filter_score {
                        0 => "-",
                        15 => "PASS",
                        _ => "FAIL"
                    };
                    print!("\t{sample_name}={cov}:{major}:{minor}:{fr}:{score}:{passed}", 
                        sample_name=sample_name,
                        cov=cov,
                        major=major,
                        minor=minor,
                        fr=fr_ratio,
                        score=filter_score,
                        passed=passed
                    );
                }
            }
        }
    }

}
