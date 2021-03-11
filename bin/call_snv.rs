use std::path::Path;
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{self, BufRead};

use rust_htslib::tbx::{self, Read as TbxRead};

use htsops::caller::{filter_control, quality_clean_pileups, check_snv_regionwide};
use htsops::filter::{FilterParameters, SNVParameters, TargetFilterScore, SNVScore};
use htsops::pileup::SpatialSitePileup;

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




fn main() {
    let args: Vec<String> = env::args().collect();
    let mpileup_path = Path::new(&args[1]);
    let sample_list_path = &args[2];

    // let pileup_db_path = Path::new("/Users/kent/Desktop/UPN39.pileup.sled");
    // let mpileup_path = Path::new("/Volumes/NVME4TB/local_data/yokoyama/UPN39.pileup.bq0.mq0.adjmq50.bgz");
    let sample_list: Vec<&str> = vec![
        // "UPN59_B",
        // "UPN59_1_1",
        // "UPN59_1_3",
        // "UPN59_1_4",
        // "UPN59_1_5",
        // "UPN59_1_6",
        // "UPN59_1_7",
        // "UPN59_1_8",
        // "UPN59_1_9",
        // "UPN59_1_10",
        // "UPN59_1_11",
        // "UPN59_1_12",
        // "UPN59_1_13",
        // "UPN59_1_14",
        // "UPN59_1_15",
        // "UPN59_1_16",
        // "UPN59_1_17",
        // "UPN59_1_18",
        // "UPN59_1_19",
        // "UPN59_1_21",
        // "UPN59_1_22",
        // "UPN59_1_23",
        // "UPN59_1_24",
        // "UPN59_1_25",
        "UPN39_B",
        "UPN39_2_1",
        "UPN39_2_2",
        "UPN39_2_3",
        "UPN39_2_4",
        "UPN39_2_5",
        "UPN39_2_6",
        "UPN39_2_7",
        "UPN39_2_8",
        "UPN39_2_9",
        "UPN39_2_10",
        "UPN39_2_11",
        "UPN39_2_12",
        "UPN39_2_13",
        "UPN39_2_14",
        "UPN39_2_15",
        "UPN39_2_16",
        "UPN39_2_17",
        "UPN39_2_18",
        "UPN39_2_20",
        "UPN39_2_21",
        "UPN39_2_22",
        "UPN39_2_23",
        "UPN39_2_24",
        "UPN39_2_25",
        "UPN39_2_26",
        "UPN39_2_27",
        "UPN39_2_28",
        "UPN39_2_29",
        "UPN39_2_30",
        "UPN39_2_31",
    ];
    // TODO: read sample_list_path and parse line by line for sample names
    // expectation is control name is first
    let file = File::open(sample_list_path).unwrap();
    let sample_list: Vec<&str> = io::BufReader::new(file).lines()
        .filter_map(|l| if let Ok(line) = l { Some(line) } else { None })
        .map(|line| {
            line.split('.').collect::<Vec<&str>>()[0]
        } )
        .collect();
    let control_name: &str = sample_list[0];

    // Program settings
    const NUM_THREADS: usize = 4;
    const FETCH_CHUNKSIZE: u64 = 1_000_000;

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
    let snv_params = SNVParameters {
        min_maf: 0.01,
        min_mac: 2,
        max_oaf: 0.01,
        max_oac: 1,
        minor_fr_ratio_range: (0.3, 0.7),
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
    tbx_reader.set_threads(NUM_THREADS).unwrap();

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
                let mut spatial = SpatialSitePileup::parse_mpileup_row(&record_string, &sample_list, control_name);

                // Process control pileup first
                // Continue only if control passed all filters
                let control_filter_result = filter_control(&mut spatial.control_pileup, &control_filter_params);
                if !control_filter_result.score.is_all() {
                    continue
                }
                // Process non-control pileups at the same site
                quality_clean_pileups(&mut spatial, &control_filter_params);
    
                // Determine region-wide major and minor allele
                // Call SNV for each sample
                let target_results = match check_snv_regionwide(&mut spatial, &control_filter_result, &target_filter_params, &snv_params) {
                    Some(r) => r,
                    None => continue,
                };

                // print results
                let site_info = format!(
                    "{chrom}\t{pos}\t{ref_char}\t{major_char}\t{minor_char}",
                    chrom=spatial.chrom, 
                    pos=spatial.pos,
                    ref_char=spatial.ref_char,
                    major_char=control_filter_result.naive_major_allele(),
                    minor_char=target_results.minor_allele(),
                );
                let site_info_numcols = 5;

                let control_info = format!(
                    "{ctrl_filt_bitscore}\t{ctrl_fr_ratio:.6}\t{f_a}:{r_a}:{f_c}:{r_c}:{f_g}:{r_g}:{f_t}:{r_t}",
                    ctrl_filt_bitscore=control_filter_result.score.to_bits(), 
                    ctrl_fr_ratio=control_filter_result.fr_ratio,
                    f_a=control_filter_result.full_base_count.f_a(),
                    r_a=control_filter_result.full_base_count.r_a(),
                    f_c=control_filter_result.full_base_count.f_c(),
                    r_c=control_filter_result.full_base_count.r_c(),
                    f_g=control_filter_result.full_base_count.f_g(),
                    r_g=control_filter_result.full_base_count.r_g(),
                    f_t=control_filter_result.full_base_count.f_t(),
                    r_t=control_filter_result.full_base_count.r_t(),
                );
                let control_info_numcols = 3;

                let pooled_info = format!(
                    "{pooled_fr_ratio:.6}\t{pooled_major_fr_ratio:.6}\t{pooled_minor_fr_ratio:.6}\t{f_a}:{r_a}:{f_c}:{r_c}:{f_g}:{r_g}:{f_t}:{r_t}\t{pooled_bias_pval:.6}",
                    pooled_fr_ratio=target_results.pooled_fr_ratio(), 
                    pooled_major_fr_ratio=target_results.pooled_major_fr_ratio(),
                    pooled_minor_fr_ratio=target_results.pooled_minor_fr_ratio(),
                    f_a=target_results.pooled_full_base_count().f_a(),
                    r_a=target_results.pooled_full_base_count().r_a(),
                    f_c=target_results.pooled_full_base_count().f_c(),
                    r_c=target_results.pooled_full_base_count().r_c(),
                    f_g=target_results.pooled_full_base_count().f_g(),
                    r_g=target_results.pooled_full_base_count().r_g(),
                    f_t=target_results.pooled_full_base_count().f_t(),
                    r_t=target_results.pooled_full_base_count().r_t(),
                    pooled_bias_pval=target_results.pooled_bias_pval(),
                );
                let pooled_info_numcols = 5;

                let sample_section: String = target_results.sample_results.iter()
                    .map(|result| {
                        format!(
                            "{sample_name}\t{filt_bitscore}\t{snv_bitscore}\t{fr_ratio:.6}\t{major_fr_ratio:.6}\t{minor_fr_ratio:.6}\t{f_a}:{r_a}:{f_c}:{r_c}:{f_g}:{r_g}:{f_t}:{r_t}\t{bias_pval:.6}",
                            sample_name=result.sample_name,
                            filt_bitscore=result.filt_bitscore.to_bits(),
                            snv_bitscore=result.snv_bitscore.to_bits(),
                            fr_ratio=result.fr_ratio,
                            major_fr_ratio=result.major_fr_ratio,
                            minor_fr_ratio=result.minor_fr_ratio,
                            f_a=result.full_base_count.f_a(),
                            r_a=result.full_base_count.r_a(),
                            f_c=result.full_base_count.f_c(),
                            r_c=result.full_base_count.r_c(),
                            f_g=result.full_base_count.f_g(),
                            r_g=result.full_base_count.r_g(),
                            f_t=result.full_base_count.f_t(),
                            r_t=result.full_base_count.r_t(),
                            bias_pval=result.bias_pval,
                        )
                    }).collect::<Vec<String>>().join("\t");
                let pooled_info_numcols = 8;

                let passed_samples_bitstring = target_results.passed_samples_bitstring(TargetFilterScore::all(), SNVScore::all());
                println!(
                    "{}\t{}\t{}\t{}\t{}",
                    site_info,
                    passed_samples_bitstring,
                    control_info,
                    pooled_info,
                    sample_section,
                );
            }
        }
    }
}
