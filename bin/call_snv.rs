#[macro_use]
use std::path::Path;
use std::collections::HashMap;
use rust_htslib::tbx::{self, Read as TbxRead};

use htsops::pileup::{SitePileup, SpatialSitePileup, ControlFilterResult};
use htsops::filter::{FilterParameters, FilterResult, SNVCall, SNVResult, AlleleCount, AlleleFreq};

fn filter_normal(pileup: &mut SitePileup, params: &FilterParameters) -> CallResult {
    // bit score
    // 1 - passed quality filter
    // 2 - passed min cov post quality filter
    // 4 - has 1 or 2 alleles // consider splitting
    // 8 - has a forward-reverse ratio within bounds
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

fn call_regionwide_minor_allele(multi_pileup: &mut SpatialSitePileup, 
        sample_names: &Vec<&str>, major_allele: &char) -> Option<AlleleCount> {
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
    // Major is based on control allele, 
    // this must also be the top in the pooled
    // Minor is second highest overall
    // Major + minor should be greater than 95% of total coverage
    let mut alleles: Vec<AlleleCount> = vec![
        ('a', tot_a), ('c', tot_c), ('g', tot_g), ('t', tot_t)]
    .into_iter()
    .filter_map(|(k, v)| match v {
        0 => None,
        _ => Some(AlleleCount(k, v))
    }).collect();
    alleles.sort_by(|AlleleCount(_, a), AlleleCount(_, b)| b.cmp(a));
    if alleles.len() > 0 {
        if major_allele != alleles[0].base() {
            panic!("Major allele in control [{}] and target pool [{}] do not match", 
                major_allele, alleles[0].base());
        }
        if alleles.len() == 1 { return None }
        else if alleles.len() == 2 { return Some(alleles[1]) }
        else if alleles.len() > 2 {
            let total_count = (tot_a + tot_c + tot_g + tot_t) as f64;
            let total_majmin = (alleles[0].count() + alleles[0].count()) as f64;
            if total_majmin / total_count >= 0.95 {
                return Some(alleles[1])
            }
        }
    }
    return None
}

fn call_snv_on_target(pileup: &mut SitePileup, params: &CallSNVParameters) -> CallResult {
    // bit score
    // 1 - passed quality filter
    // 2 - passed min cov post quality filter
    // 4 - has 1 or 2 alleles // consider splitting
    // 8 - has a forward-reverse ratio within bounds
    let mut score = 0;
    // Quality filter
    let cov = match pileup.quality_filter(params.min_bq, params.min_mq, params.drop_n) {
        Some(c) if c >= params.min_cov => { score += 3; c },
        Some(c) if c < params.min_cov => { score += 1; c },
        _ => 0,
    };
    // Major allele and count
    let (a, c, g, t, n) = pileup.base_count();
    let major_count = match params.major_allele.to_ascii_lowercase() {
        'a' => a,
        'c' => c,
        'g' => g,
        't' => t,
        'n' => n,
        b => panic!("Invalid base [{}]", b),
    };
    let minor_count = match params.minor_allele.to_ascii_lowercase() {
        'a' => a,
        'c' => c,
        'g' => g,
        't' => t,
        'n' => n,
        b => panic!("Invalid base [{}]", b),
    };
    // 

    
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
