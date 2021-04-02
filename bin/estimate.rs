// estimate.rs
//
// First-pass program to estimate and record various statistics about each pile up site.
// This program writes data in tabular form to stdout.
// The output can be indexed by the first two columns (chrom and pos) using tabix.
// Header:
// line 0 - path of pileup being read
// line 1 - comma-sep sample names
// line 2 - comma-sep column names
// Table format:
//  0 - chrom
//  1 - pos
//  2 - ref_char (of hg19) in big caps
//  3 - estimation bitflags for the site over the control and all regions
//  4 - Fisher exact test for I(control, pool), J(control major count, not major), 2-tailed pval
//  5 - CMH text for I(control, sample), J(control major count, not major), over samples, 2-tailed 
//  6 - control_bq_histogram in intervals of 5 as counts [0]:[1-5]:[6-10]:[11-15]:[16-20]:[21-25]:[26-30]:[31-35]:[36-40]:[41+]
//  7 - control_mq_histogram in intervals of 10 as counts [0]:[1-10]:[11-20]:[21-30]:[31-40]:[41-40]:[51-60]:[61+]
//  8 - control_raw_basecount as A:a:C:c:G:g:T:t
//  9 - control_qf_basecount as A:a:C:c:G:g:T:t for bases q>20 and reads q>40
// 10 - control_qf_cov
// 11 - pooled_bq_histogram in intervals of 5 as counts [0]:[1-5]:[6-10]:[11-15]:[16-20]:[21-25]:[26-30]:[31-35]:[36-40]:[41+]
// 12 - pooled_mq_histogram in intervals of 10 as counts [0]:[1-10]:[11-20]:[21-30]:[31-40]:[41-40]:[51-60]:[61+]
// 13 - pooled_raw_basecount as A:a:C:c:G:g:T:t
// 14 - pooled_qf_basecount as A:a:C:c:G:g:T:t for bases q>20 and reads q>40
// Sample data repeats from colum 15 to 19
// 15 - sample_bq_histogram in intervals of 5 as counts [0]:[1-5]:[6-10]:[11-15]:[16-20]:[21-25]:[26-30]:[31-35]:[36-40]:[41+]
// 16 - sample_mq_histogram in intervals of 10 as counts [0]:[1-10]:[11-20]:[21-30]:[31-40]:[41-40]:[51-60]:[61+]
// 17 - sample_raw_basecount as A:a:C:c:G:g:T:t
// 18 - sample_qf_basecount as A:a:C:c:G:g:T:t for bases q>20 and reads q>40
// 19 - sample_qf_cov

use std::path::Path;
use std::collections::{HashMap, BTreeMap};
use std::fs::File;
use std::io::{self, BufRead};

use clap::{Arg, App};
use rust_htslib::tbx::{self, Read as TbxRead};
use fishers_exact::fishers_exact;

use htsops::pileup::{SpatialSitePileup, SitePileup, FullBaseCount, AlleleSet};

fn validate_qual(v: String) -> Result<(), String> {
    match v.parse::<u8>() {
        Ok(_) => Ok(()),
        Err(_) => Err(String::from("Value is not a positive integer")),
    }
}

fn validate_path(v: String) -> Result<(), String> {
    let mpileup_path = Path::new(&v);
    match mpileup_path.is_file() {
        true => Ok(()),
        false => Err(String::from("Input file does not exist")),
    }
}

fn get_sample_names(path_str: &str) -> (String, Vec<String>) {
    let file = File::open(path_str).unwrap();
    let sample_names: Vec<String> = io::BufReader::new(file).lines()
        .filter_map(|l| {
            if let Ok(line) = l { 
                if line.len() > 0 { Some(line) } else { None }
            } 
            else { None }
        })
        .map(|line| {
            // TODO: Change to more robust solution for paths
            let split: Vec<&str> = line.split('/').collect();
            let filename = split[split.len() - 1];
            let sample_name = filename.split('.').collect::<Vec<&str>>()[0];
            sample_name.to_owned()
        })
        .collect();
    let control_name: String = sample_names[0].to_owned();
    (control_name, sample_names)
}

fn flatten<T>(nested: Vec<Vec<T>>) -> Vec<T> {
    nested.into_iter().flatten().collect()
}

fn site_bq_histogram_str(scores: &Vec<u8>) -> String {
    let mut hist: HashMap<usize, usize> = vec![
            (0, 0),
            (5, 0),
            (10, 0),
            (15, 0),
            (20, 0),
            (25, 0),
            (30, 0),
            (35, 0),
            (40, 0),
            (41, 0),
        ].into_iter().collect();
    scores.iter().for_each(|s| {
        match s {
            0 => *hist.entry(0).or_insert(0) += 1,
            v if *v >  0 && *v <=  5 => *hist.entry( 5).or_insert(0) += 1,
            v if *v >  5 && *v <= 10 => *hist.entry(10).or_insert(0) += 1,
            v if *v > 10 && *v <= 15 => *hist.entry(15).or_insert(0) += 1,
            v if *v > 15 && *v <= 20 => *hist.entry(20).or_insert(0) += 1,
            v if *v > 20 && *v <= 25 => *hist.entry(25).or_insert(0) += 1,
            v if *v > 25 && *v <= 30 => *hist.entry(30).or_insert(0) += 1,
            v if *v > 30 && *v <= 35 => *hist.entry(35).or_insert(0) += 1,
            v if *v > 35 && *v <= 40 => *hist.entry(40).or_insert(0) += 1,
            _                      => *hist.entry(41).or_insert(0) += 1,
        }
    });
    hist.into_iter().map(|(_,v)| format!("{}", v)).collect::<Vec<String>>().join(":")
}

fn site_mq_histogram_str(scores: &Vec<u8>) -> String {
    let mut hist: HashMap<usize, usize> = vec![
            (0, 0),
            (10, 0),
            (20, 0),
            (30, 0),
            (40, 0),
            (50, 0),
            (60, 0),
            (61, 0),
        ].into_iter().collect();
    scores.iter().for_each(|s| {
        match s {
            0 => *hist.entry(0).or_insert(0) += 1,
            v if *v >  0 && *v <= 10 => *hist.entry(10).or_insert(0) += 1,
            v if *v > 10 && *v <= 20 => *hist.entry(20).or_insert(0) += 1,
            v if *v > 20 && *v <= 30 => *hist.entry(30).or_insert(0) += 1,
            v if *v > 30 && *v <= 40 => *hist.entry(40).or_insert(0) += 1,
            v if *v > 40 && *v <= 50 => *hist.entry(50).or_insert(0) += 1,
            v if *v > 50 && *v <= 60 => *hist.entry(60).or_insert(0) += 1,
            _                      => *hist.entry(61).or_insert(0) += 1,
        }
    });
    hist.into_iter().map(|(_,v)| format!("{}", v)).collect::<Vec<String>>().join(":")
}

fn qf_base_count(pileup: &SitePileup, min_bq: u8, min_mq: u8) -> FullBaseCount {
    if pileup.bases.len() != pileup.bqs.len() {
        panic!("bases len ({}) != bqs len ({})", pileup.bases.len(), pileup.bqs.len())
    }
    if pileup.bases.len() != pileup.mqs.len() {
        panic!("bases len ({}) != mqs len ({})", pileup.bases.len(), pileup.mqs.len())
    }
    let passed_bases: Vec<char> = pileup.bases.iter().enumerate()
        .filter_map(|(i, &b)| {
            if (b == 'N') || (b == 'n') { return None }
            if (b == 'D') || (b == 'd') { return None }
            let bq = pileup.bqs[i];
            let mq = pileup.mqs[i];
            if (bq >= min_bq) && (mq >= min_mq) { return Some(b) }
            return None
        }).collect();
    FullBaseCount::from_vec(&passed_bases)
}

fn main() {
    // Argument parser to get input path and minimum quality params
    let matches = App::new("htsops estimate.rs")
        .version("0.0.1")
        .author("Kent Kawashima <kentkawashima@gmail.com>")
        .about("Outputs per-site statistics for variant calling")
        .arg(Arg::with_name("min_bq")
            .short("b")
            .long("min_bq")
            .default_value("20")
            .value_name("INT")
            .help("Minimum base quality. Bases whose base quality is less than this value will not be counted.")
            .takes_value(true)
            .validator(validate_qual))
        .arg(Arg::with_name("min_mq")
            .short("m")
            .long("min_mq")
            .default_value("40")
            .value_name("INT")
            .help("Minimum mapping quality. Bases belonging to reads with mapping quality less than this value will not be counted.")
            .takes_value(true)
            .validator(validate_qual))
            .arg(Arg::with_name("min_cov")
            .short("c")
            .long("min_cov")
            .default_value("20")
            .value_name("INT")
            .help("Minimum control coverage depth. Sites where the control coverage depth is less than this value will not be printed.")
            .takes_value(true)
            .validator(validate_qual))
        .arg(Arg::with_name("threads")
            .short("t")
            .long("threads")
            .default_value("1")
            .value_name("INT")
            .help("Number of threads to use")
            .takes_value(true))
        .arg(Arg::with_name("INPUT")
            .help("Sets the input file to use")
            .required(true)
            .index(1)
            .validator(validate_path))
        .arg(Arg::with_name("SAMPLE_LIST")
            .help("Sets the sample list to use for parsing the input")
            .required(true)
            .index(2)
            .validator(validate_path))
        .get_matches();

    let min_bq = matches.value_of("min_bq").unwrap_or("20").parse::<u8>().unwrap();
    let min_mq = matches.value_of("min_mq").unwrap_or("40").parse::<u8>().unwrap();
    let min_cov = matches.value_of("min_cov").unwrap_or("20").parse::<usize>().unwrap();
    let threads = matches.value_of("threads").unwrap_or("1").parse::<usize>().unwrap();
    let mpileup_path = matches.value_of("INPUT").unwrap();
    let sample_list_path = matches.value_of("SAMPLE_LIST").unwrap();

    // Set sample_names and control name
    let (control_name, sample_names) = get_sample_names(sample_list_path);
    let sample_names = sample_names.iter().map(|s| &**s).collect::<Vec<&str>>();

    // Open tabix-indexed bgzipped mpileup file
    let mut tbx_reader = tbx::Reader::from_path(mpileup_path)
        .expect(&format!("Tabix reader could not open {}", mpileup_path));
    tbx_reader.set_threads(threads).unwrap();

    // Resolve chromosome name to numeric ID.
    const FETCH_CHUNKSIZE: u64 = 1_000_000;
    let target_chrom = vec![
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
    let chrom_tid_lookup: HashMap<&str, u64> = target_chrom.iter().map(|c| {
            match tbx_reader.tid(c) {
                Ok(tid) => (c.to_owned(), tid),
                Err(_) => panic!("Could not resolve contig ID"),
            }
        }).collect();

    // Collect stats
    // histogram of coverage per chromosome
    // key1: chrom, key2: cov, value: counter
    let mut ctrl_cov_freq: HashMap<&str,  BTreeMap<usize, usize>> = HashMap::new();
    // histogram of major allele frequency per chromosome
    // key1: chrom, key2: major allele frequency as usize where 1.0 is 1000, value: counter
    let mut ctrl_maf_freq: HashMap<&str,  BTreeMap<usize, usize>> = HashMap::new();
    let mut ctrl_totals: HashMap<&str, usize> = HashMap::new();

    // Print header
    let sample_col_names: Vec<String> = sample_names.iter()
    .map(|s| {
        format!("{s}_bq_hist,{s}_mq_hist,{s}_raw_bc,{s}_bc,{s}_cov", s=s)
    })
    .collect();
    println!("# mpileup: {}", mpileup_path);
    println!("# min_bq: {}", min_bq);
    println!("# min_mq: {}", min_mq);
    println!("# min_cov: {}", min_cov);
    println!("# sample_names: {}", sample_names.join(","));
    println!("# column_names: chrom,pos,ref_char,summary_flags,reads_or,exact_test,cmh_test,ctrl_bq_hist,ctrl_mq_hist,ctrl_raw_bc,ctrl_bc,ctrl_cov,pool_bq_hist,pool_mq_hist,pool_raw_bc,pool_bc,{}", sample_col_names.join(","));
    println!("# bq_histogram_fmt: [0]:[1,5]:[6,10]:[11,15]:[16,20]:[21,25]:[26,30]:[31,35]:[36,40]:[41,+)");
    println!("# mq_histogram_fmt: [0]:[1,10]:[11,20]:[21,30]:[31,40]:[41,40]:[51,60]:[61,+)");
    println!("# raw_base_count_fmt: A:a:C:c:G:g:T:t");
    println!("# base_count_fmt: A:a:C:c:G:g:T:t");

    // TODO: Make a mpileup iterator
    let mut total_count = 0;
    for &chrom in target_chrom.iter() {
        let ul: u64 = (chrom_size.get(chrom).unwrap().to_owned() as u64 / FETCH_CHUNKSIZE) + 1;
        // Initialize counters for chromosome
        ctrl_cov_freq.insert(chrom, BTreeMap::new());
        let chrom_ctrl_cov_freq = ctrl_cov_freq.get_mut(chrom).unwrap();
        ctrl_maf_freq.insert(chrom, BTreeMap::new());
        let chrom_ctrl_maf_freq = ctrl_maf_freq.get_mut(chrom).unwrap();
        let mut chrom_total = 0;

        for i in 0..ul {
            tbx_reader.fetch(*chrom_tid_lookup.get(chrom).unwrap(), i*FETCH_CHUNKSIZE, (i+1)*FETCH_CHUNKSIZE).unwrap();
            for record in tbx_reader.records() {

                // Convert binary mpileup line and parse
                let record_string = String::from_utf8(record.unwrap()).unwrap();
                let sp_data = SpatialSitePileup::parse_mpileup_row(&record_string, &sample_names, &control_name);

                let ctrl_bc = qf_base_count(&sp_data.control_pileup, min_bq, min_mq);
                let ctrl_cov = ctrl_bc.total();
                if ctrl_cov < min_cov {
                    continue
                }

                let site_info = format!(
                    "{chrom}\t{pos}\t{ref_char}\t{bitflag}",
                    chrom=sp_data.chrom, 
                    pos=sp_data.pos,
                    ref_char=sp_data.ref_char.to_ascii_uppercase(),
                    bitflag=0,  // change later
                );

                let ctrl_info = format!(
                    "{ctrl_bq_hist}\t{ctrl_mq_hist}\t{ctrl_raw_bc}\t{ctrl_bc}\t{ctrl_cov}",
                    ctrl_bq_hist=site_bq_histogram_str(&sp_data.control_pileup.bqs),
                    ctrl_mq_hist=site_mq_histogram_str(&sp_data.control_pileup.mqs),
                    ctrl_raw_bc=sp_data.control_pileup.full_base_count(),
                    ctrl_bc=ctrl_bc,
                    ctrl_cov=ctrl_cov,
                );
                let site_info_list: Vec<_> = sp_data.pileups.iter().map(|p| {
                    let fbc = qf_base_count(p, min_bq, min_mq);
                    let cov = fbc.total();
                    (
                        site_bq_histogram_str(&p.bqs),
                        site_mq_histogram_str(&p.mqs),
                        p.full_base_count(),
                        fbc,
                        cov,
                    )
                }).collect();
                let mut pool_raw_bc: FullBaseCount = FullBaseCount::empty();
                let mut pool_bc: FullBaseCount = FullBaseCount::empty();
                let mut pool_cov: usize = 0;
                for (_, _, raw_bc, bc, cov) in site_info_list.iter() {
                    pool_raw_bc = pool_raw_bc + raw_bc;
                    pool_bc = pool_bc + bc;
                    pool_cov += cov;
                }

                let pool_info = format!(
                    "{pool_bq_hist}\t{pool_mq_hist}\t{pool_raw_bc}\t{pool_bc}\t{pool_cov}",
                    pool_bq_hist=site_bq_histogram_str(&flatten(sp_data.pileups.iter().map(|p| p.bqs.clone()).collect())),
                    pool_mq_hist=site_mq_histogram_str(&flatten(sp_data.pileups.iter().map(|p| p.mqs.clone()).collect())),
                    pool_raw_bc=pool_raw_bc,
                    pool_bc=pool_bc,
                    pool_cov=pool_cov,
                );

                let sample_info = site_info_list.iter()
                    .map(|(bq_hist, mq_hist, raw_bc, bc, cov)| {
                        format!(
                            "{bq_hist}\t{mq_hist}\t{raw_bc}\t{bc}\t{cov}",
                            bq_hist=bq_hist,
                            mq_hist=mq_hist,
                            raw_bc=raw_bc,
                            bc=bc,
                            cov=cov,
                        )
                    })
                    .collect::<Vec<String>>()
                    .join("\t");
                
                // Update control cov depth counter
                *chrom_ctrl_cov_freq.entry(ctrl_cov).or_insert(0) += 1;

                // Update control maf counter
                let ctrl_allele_set = AlleleSet::from_fullbasecount(&ctrl_bc);
                // let major_char = ctrl_allele_set.alleles[0].base();
                let ctrl_major_base = ctrl_allele_set.alleles[0].base();
                let ctrl_major_count = ctrl_allele_set.alleles[0].count();
                let ctrl_major_freq = ((ctrl_major_count as f64 / ctrl_cov as f64) * 1000.0).round() as usize;
                *chrom_ctrl_maf_freq.entry(ctrl_major_freq).or_insert(0) += 1;
                total_count += 1;
                chrom_total += 1;

                let pool_major_count = pool_bc.count_of(ctrl_major_base);
                let contingency_tbl: [u32; 4] = [
                    ctrl_major_count as u32, (ctrl_cov-ctrl_major_count) as u32,
                    pool_major_count as u32, (pool_cov-pool_major_count) as u32,
                ];
                let exact_test = fishers_exact::fishers_exact(&contingency_tbl).unwrap();
                let ctrl_reads_or: f64 = {
                    let reverse = ctrl_bc.reverse();
                    if reverse == 0 {
                        -1.0
                    } else {
                        (ctrl_bc.forward() as f64) / (reverse as f64)
                    }
                };
                let site_stats = format!(
                    "{reads_or:>7.4}\t{exact_test:.4}\t{cmh_test:.4}",
                    reads_or=ctrl_reads_or,
                    exact_test=exact_test.two_tail_pvalue,
                    cmh_test="",
                );

                // Output data for control + all locations
                println!(
                    "{site_info}\t{site_stats}\t{ctrl_info}\t{pool_info}\t{sample_info}",
                    site_info=site_info,
                    site_stats=site_stats,
                    ctrl_info=ctrl_info,
                    pool_info=pool_info,
                    sample_info=sample_info,
                );
            }
            break
        }
        ctrl_totals.insert(chrom, chrom_total);
    }

    // Print stats
    println!("# mpileup: {}", mpileup_path);
    println!("# Total non-zero cov sites: {}", total_count);
    // Print coverage histogram
    println!("# Control per-chromosome coverage histogram");
    for &chrom in target_chrom.iter() {
        if let Some(chrom_ctrl_cov_freq) = ctrl_cov_freq.get(chrom) {
            if chrom_ctrl_cov_freq.len() == 0 { continue }
            println!("# {:<10} total:{:>8}", chrom, ctrl_totals.get(chrom).unwrap());
            for (cov, count) in chrom_ctrl_cov_freq.iter() {
                let freq: f64 = (*count as f64) / (total_count as f64);
                println!("# \t{:>5}:\t{:>7}\t{:.6}", cov, count, freq);
            }
        }
    }
    // Print major allele frequency
    println!("# mpileup: {}", mpileup_path);
    println!("# Total non-zero cov sites: {}", total_count);
    println!("# Control per-chromosome major allele freq histogram");
    for chrom in target_chrom.iter() {
        if let Some(chrom_ctrl_maf_freq) = ctrl_maf_freq.get(chrom) {
            if chrom_ctrl_maf_freq.len() == 0 { continue }
            println!("# {:<10} total:{:>8}", chrom, ctrl_totals.get(chrom).unwrap());
            for (maf, count) in chrom_ctrl_maf_freq.iter() {
                let maf: f64 = (*maf as f64) / 1000.0;
                let freq: f64 = (*count as f64) / (total_count as f64);
                println!("# \t\t{:.3}:\t{:>7}\t{:.6}", maf, count, freq);
            }
        }
    }
}