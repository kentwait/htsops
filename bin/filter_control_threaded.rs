// spatial_pileup.rs
// Inputs:
// - Filelist of mpileups to process, 1 pileup path per line
//   Assumes control is at the top 
// - Coords list
//   Format is "{row}\t{col} per line, corresponds to the order of the pileup filelist
//   Note that origin (0,0) is top-left position
// Parameters:
// --min-control-cov : Minimum coverage depth for the control calculated after bq and mq minimums
// --min-control-bq  : Minimum base quality to consider in the control
// --min-control-mq  : Minimum mapping quality to consider in the control
// --min-target-cov  : Minimum coverage depth for each sample
// --min-target-bq   : Minimum base quality to consider for each sample
// --min-target-mq   : Minimum mapping quality to consider for each sample
// --threads         : Number of threads to use
use std::collections::{HashMap, HashSet, BTreeMap, VecDeque};
use std::fs::File;
use std::io::{BufWriter, Write};

use clap::{Arg, App};
use rust_htslib::tbx::{self, Read as TbxRead};
use indexmap::IndexMap;
use std::time::Duration;
use crossbeam::channel::{self, select, bounded};
use crossbeam::thread;


use htsops::util::{validate_path, validate_qual, read_tabix};
use htsops::pileup::{SitePileup, FullBaseCount, AlleleSet};
use htsops::filter::ControlFilterScore;
use htsops::constant::*;

fn site_bq_histogram_str(scores: &Vec<u8>) -> String {
    let mut hist: BTreeMap<usize, usize> = vec![
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
    // format!("{:?}", scores)
    scores.into_iter().for_each(|s| {
        match *s {
            v if v == 0            => *hist.entry( 0).or_insert(0) += 1,
            v if v >  0 && v <=  5 => *hist.entry( 5).or_insert(0) += 1,
            v if v >  5 && v <= 10 => *hist.entry(10).or_insert(0) += 1,
            v if v > 10 && v <= 15 => *hist.entry(15).or_insert(0) += 1,
            v if v > 15 && v <= 20 => *hist.entry(20).or_insert(0) += 1,
            v if v > 20 && v <= 25 => *hist.entry(25).or_insert(0) += 1,
            v if v > 25 && v <= 30 => *hist.entry(30).or_insert(0) += 1,
            v if v > 30 && v <= 35 => *hist.entry(35).or_insert(0) += 1,
            v if v > 35 && v <= 40 => *hist.entry(40).or_insert(0) += 1,
            _                      => *hist.entry(41).or_insert(0) += 1,
        }
    });
    hist.into_iter().map(|(_,v)| format!("{}", v)).collect::<Vec<String>>().join(":")
}

fn site_mq_histogram_str(scores: &Vec<u8>) -> String {
    let mut hist: BTreeMap<usize, usize> = vec![
            (0, 0),
            (10, 0),
            (20, 0),
            (30, 0),
            (40, 0),
            (50, 0),
            (60, 0),
            (61, 0),
        ].into_iter().collect();
    // format!("{:?}", scores)
    scores.into_iter().for_each(|s| {
        match *s {
            v if v == 0            => *hist.entry( 0).or_insert(0) += 1,
            v if v >  0 && v <= 10 => *hist.entry(10).or_insert(0) += 1,
            v if v > 10 && v <= 20 => *hist.entry(20).or_insert(0) += 1,
            v if v > 20 && v <= 30 => *hist.entry(30).or_insert(0) += 1,
            v if v > 30 && v <= 40 => *hist.entry(40).or_insert(0) += 1,
            v if v > 40 && v <= 50 => *hist.entry(50).or_insert(0) += 1,
            v if v > 50 && v <= 60 => *hist.entry(60).or_insert(0) += 1,
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

fn hist_mean(hist: &BTreeMap<usize, usize>, total: usize) -> f64 {
    let mean: f64 = hist.iter().map(|(v, i)| (v*i) as f64 ).sum::<f64>() / (total as f64);
    mean
}

fn hist_median(hist: &BTreeMap<usize, usize>, total: usize) -> f64 {
    let midpoint = ((total as f64) / 2.0).round() as usize;
    let median: f64 = {
        let mut ci = 0;
        let mut median = 0;
        for (v, i) in hist.iter() {
            ci += i;
            if ci >= midpoint { 
                median = *v;
                break
            }
        }
        median as f64
    };
    median
}

fn hist_to_file(per_chrom_hist: &IndexMap<&str, BTreeMap<usize, usize>>, per_chrom_total: &IndexMap<&str, usize>, mut f: BufWriter<File>) {
    let total_sites = per_chrom_total.iter().map(|(_,c)| c ).sum::<usize>();
    let total_counted = per_chrom_hist.iter().map(|(_,hist)| {
            hist.iter().map(|(_,c)| c ).sum::<usize>()
        })
        .sum::<usize>();
    // Overall
    let mut all_hist: BTreeMap<usize, usize> = BTreeMap::new();
    for (_, hist) in per_chrom_hist.iter() {
        for (&val, &count) in hist.iter() {
            *all_hist.entry(val).or_insert(0) += count;
        }
    }
    // Print overall stats
    f.write_fmt(format_args!("# Overall histogram\n"))
        .expect("Unable to write data");
    f.write_fmt(format_args!("# total_sites:{}\n", total_sites))
        .expect("Unable to write data");
    f.write_fmt(format_args!("# total_counted:{}\n", total_counted))
        .expect("Unable to write data");
    let mean: f64 = hist_mean(&all_hist, total_counted);
    let median: f64 = hist_median(&all_hist, total_counted);
    f.write_fmt(format_args!("# mean:{mean:.2}\tmedian:{median:.1}\n",
        mean=mean, median=median)).expect("Unable to write data");
    let mut cum_freq = 0.0;
    for (val, count) in all_hist.iter() {
        let freq: f64 = (*count as f64) / (total_counted as f64);
        cum_freq += freq;
        f.write_fmt(format_args!("\t{:>5}:\t{:>9}\t{:.6}\t{:.6}\n",
            val, count, freq, cum_freq)).expect("Unable to write data");
    }        
    f.flush().expect("Unable to flush data");
    // Print coverage histogram
    f.write_fmt(format_args!("# Per-chromosome histogram\n")).expect("Unable to write data");
    for (&chrom, hist) in per_chrom_hist.iter() {
        if hist.len() == 0 { continue }
        // calculate mean
        let this_chrom_total: usize = per_chrom_total.get(chrom).unwrap().to_owned();
        let this_chrom_counted: usize = per_chrom_hist.get(chrom).unwrap().iter().map(|(_,c)| c ).sum::<usize>();
        let mean: f64 = hist_mean(&hist, this_chrom_counted);
        let median: f64 = hist_median(&hist, this_chrom_counted);
        f.write_fmt(format_args!("{:<10}\ttotal:{:>9}\tcounted:{:>9}\tmean:{mean:.2}\tmedian:{median:.1}\n", 
            chrom, this_chrom_total, this_chrom_counted,
            mean=mean, median=median)).expect("Unable to write data");
        let mut cum_freq = 0.0;
        for (val, count) in hist.iter() {
            let freq: f64 = (*count as f64) / (this_chrom_counted as f64);
            cum_freq += freq;
            f.write_fmt(format_args!("\t{:>5}:\t{:>9}\t{:.6}\t{:.6}\n",
                val, count, freq, cum_freq)).expect("Unable to write data");
        }        
    }
    f.flush().expect("Unable to flush data");
    drop(f);
}

fn fhist_to_file(per_chrom_hist: &IndexMap<&str, BTreeMap<usize, usize>>, per_chrom_total: &IndexMap<&str, usize>, mut f: BufWriter<File>, divisor: f64) {
    let total_sites = per_chrom_total.iter().map(|(_,c)| c ).sum::<usize>();
    let total_counted = per_chrom_hist.iter().map(|(_,hist)| {
            hist.iter().map(|(_,c)| c ).sum::<usize>()
        })
        .sum::<usize>();
    // Overall
    let mut all_hist: BTreeMap<usize, usize> = BTreeMap::new();
    for (_, hist) in per_chrom_hist.iter() {
        for (&val, &count) in hist.iter() {
            *all_hist.entry(val).or_insert(0) += count;
        }
    }
    // Print overall stats
    f.write_fmt(format_args!("# Overall histogram\n"))
        .expect("Unable to write data");
    f.write_fmt(format_args!("# total_sites:{}\n", total_sites))
        .expect("Unable to write data");
    f.write_fmt(format_args!("# total_counted:{}\n", total_counted))
        .expect("Unable to write data");
    let mean: f64 = hist_mean(&all_hist, total_counted) / divisor;
    let median: f64 = hist_median(&all_hist, total_counted) / divisor;
    f.write_fmt(format_args!("# mean:{mean:.6}\tmedian:{median:.6}\n",
        mean=mean, median=median)).expect("Unable to write data");
    let mut cum_freq = 0.0;
    for (val, count) in all_hist.iter() {
        let val: f64 = (*val as f64) / divisor;
        let freq: f64 = (*count as f64) / (total_counted as f64);
        cum_freq += freq;
        f.write_fmt(format_args!("\t{:.4}:\t{:>9}\t{:.6}\t{:.6}\n",
            val, count, freq, cum_freq)).expect("Unable to write data"); 
    }        
    f.flush().expect("Unable to flush data");
    // Print coverage histogram
    f.write_fmt(format_args!("# Per-chromosome histogram\n"))
        .expect("Unable to write data");
    for (&chrom, hist) in per_chrom_hist.iter() {
        if hist.len() == 0 { continue }
        // calculate mean
        let this_chrom_total: usize = per_chrom_total.get(chrom).unwrap().to_owned();
        let this_chrom_counted: usize = per_chrom_hist.get(chrom).unwrap().iter()
            .map(|(_,c)| c ).sum::<usize>();
        let mean: f64 = hist_mean(&hist, this_chrom_counted) / divisor;
        let median: f64 = hist_median(&hist, this_chrom_counted) / divisor;
        f.write_fmt(format_args!("{:<10}\ttotal:{:>9}\tcounted:{:>9}\tmean:{mean:.6}\tmedian:{median:.6}\n", 
            chrom, this_chrom_total, this_chrom_counted,
            mean=mean, median=median)).expect("Unable to write data");
        let mut cum_freq = 0.0;
        for (val, count) in hist.iter() {
            let val: f64 = (*val as f64) / divisor;
            let freq: f64 = (*count as f64) / (this_chrom_counted as f64);
            cum_freq += freq;
            f.write_fmt(format_args!("\t{:.4}:\t{:>9}\t{:.6}\t{:.6}\n",
                val, count, freq, cum_freq)).expect("Unable to write data");
        }        
    }
    f.flush().expect("Unable to flush data");
    drop(f);
}

fn pos_to_bed(pos_vec: &Vec<(&str, u64)>, mut f: BufWriter<File>) {
    let (mut start_chrom, mut start_pos) = pos_vec[0];
    let (mut last_chrom, mut last_pos) = pos_vec[0];
    // let mut intervals: Vec<(&str, u64, u64)> = Vec::new();
    // for i in 1..pos_vec.len() {
    //     let (this_chrom, this_pos) = pos_vec[i];
    //     if (last_chrom != this_chrom) || (last_pos+1 != this_pos) {
    //         intervals.push((start_chrom, start_pos, last_pos));
    //         start_chrom = last_chrom;
    //         start_pos = last_pos;
    //     }
    //     last_chrom = this_chrom;
    //     last_pos = this_pos;
    // }
    let intervals: Vec<(&str, u64, u64)> = (1..pos_vec.len()).filter_map(|i| {
            let (this_chrom, this_pos) = pos_vec[i];
            let mut interval = None;
            if (last_chrom != this_chrom) || (last_pos+1 != this_pos) {
                interval = Some((start_chrom, start_pos, last_pos));
                start_chrom = this_chrom;
                start_pos = this_pos;
            }
            last_chrom = this_chrom;
            last_pos = this_pos;
            interval
        })
        .collect();
    for (chrom, start, end) in intervals.iter() {
        f.write_fmt(format_args!("{}\t{}\t{}\n", chrom, start, end)).expect("Unable to write data");
    }
    f.flush().expect("Unable to flush data");
    drop(f);
}

fn bytes_to_pileup(b: Vec<u8>) -> (String, u64, char, SitePileup) {
    // Convert binary mpileup line and parse
    let record_string = String::from_utf8(b).unwrap();
    let cols: Vec<String> = record_string
        .split('\t')
        .map(|s| s.to_owned())
        .collect();
    let chrom: String = cols[0].to_owned();
    let pos: u64 = cols[1].parse::<u64>().unwrap();
    let ref_char: char = cols[2]
        .chars()
        .next()
        .unwrap()
        .to_ascii_uppercase();
    let sample_name = "B";

    (chrom, pos, ref_char, SitePileup::from_str(
        sample_name,
        &ref_char,
        (&cols[3]).parse::<usize>().unwrap(), // cov
        &cols[4], &cols[5], &cols[6],         // base_str, bq_str, mq_str
    ))
}

fn qual_filter_pileup(pileup: &SitePileup, min_bq: u8, min_mq: u8) -> SitePileup {
    if pileup.bases.len() != pileup.bqs.len() {
        panic!("bases len ({}) != bqs len ({})", pileup.bases.len(), pileup.bqs.len())
    }
    if pileup.bases.len() != pileup.mqs.len() {
        panic!("bases len ({}) != mqs len ({})", pileup.bases.len(), pileup.mqs.len())
    }
    let passed_pos: Vec<usize> = pileup.bases.iter().enumerate()
        .filter_map(|(i, &b)| {
            if (b == 'N') || (b == 'n') { return None }
            if (b == 'D') || (b == 'd') { return None }
            let bq = pileup.bqs[i];
            let mq = pileup.mqs[i];
            if (bq >= min_bq) && (mq >= min_mq) { return Some(i) }
            return None
        }).collect();
    SitePileup {
        sample_name: pileup.sample_name.to_owned(),
        bases: passed_pos.iter().map(|&i| pileup.bases[i]).collect(),
        indels: HashMap::new(),
        bqs: passed_pos.iter().map(|&i| pileup.bqs[i]).collect(),
        mqs: passed_pos.iter().map(|&i| pileup.mqs[i]).collect(),
    }
}

fn main() {
    let matches = App::new("htsops spatial_pileup.rs")
        .version("0.0.1")
        .author("Kent Kawashima <kentkawashima@gmail.com>")
        .about("Processes multiple pileups into one")
        .arg(Arg::with_name("min_cov")
            .long("min-cov")
            .default_value("20")
            .value_name("INT")
            .help("Minimum coverage depth for the control sample. Sites where the coverage depth is less than this value will not be used.")
            .takes_value(true)
            .validator(validate_qual))
        .arg(Arg::with_name("min_bq")
            .long("min-bq")
            .default_value("20")
            .value_name("INT")
            .help("Minimum base quality for the control sample. Bases whose base quality is less than this value will not be counted.")
            .takes_value(true)
            .validator(validate_qual))
        .arg(Arg::with_name("min_mq")
            .long("min-mq")
            .default_value("40")
            .value_name("INT")
            .help("Minimum mapping quality for the control sample. Bases belonging to reads with mapping quality less than this value will not be counted.")
            .takes_value(true)
            .validator(validate_qual))
        // TODO: Add validation
        .arg(Arg::with_name("min_fratio")
            .long("min-fratio")
            .default_value("0.33")
            .value_name("FLOAT")
            .help("Minimum acceptable proportion of forward reads over total coverage depth")
            .takes_value(true))
        .arg(Arg::with_name("max_fratio")
            .long("max-fratio")
            .default_value("0.67")
            .value_name("FLOAT")
            .help("Maximum acceptable proportion of forward reads over total coverage depth")
            .takes_value(true))
        .arg(Arg::with_name("CONTROL_PILEUP")
            .help("Control mpileup path")
            .required(true)
            .index(1)
            .validator(validate_path))
        .arg(Arg::with_name("threads")
            .short("t")
            .long("threads")
            .default_value("4")
            .value_name("INT")
            .help("Number of threads to use")
            .takes_value(true))
        // .arg(Arg::with_name("include_bad_cov")
        //     .short("k")
        //     .long("include-bad-cov")
        //     .help("Do not output sites whose coverage depth is less than the minimum"))
        .get_matches();
    let min_cov: usize = matches.value_of("min_cov").unwrap().parse::<usize>().unwrap();
    let min_bq: u8 = matches.value_of("min_bq").unwrap().parse::<u8>().unwrap();
    let min_mq: u8 = matches.value_of("min_mq").unwrap().parse::<u8>().unwrap();
    let min_fratio: f64 = matches.value_of("min_fratio").unwrap().parse::<f64>().unwrap();
    let max_fratio: f64 = matches.value_of("max_fratio").unwrap().parse::<f64>().unwrap();
    let threads: usize = matches.value_of("threads").unwrap().parse::<usize>().unwrap();
    // let include_bad_cov: bool = matches.is_present("include_bad_cov");
    let control_path: &str = matches.value_of("CONTROL_PILEUP").unwrap();

    // Read control bam first and evaluate minimums
    let passed_simple_filter = ControlFilterScore::PassedMinCov | ControlFilterScore::PassedFRRatio | ControlFilterScore::InvariantSite;
    let mut passed_pos: Vec<(&str, u64)> = Vec::new();
    // Open control mpileup
    let (mut tbx_reader, chrom_tid_lookup) = read_tabix(control_path);

    // Collect stats
    // histogram of coverage per chromosome
    let mut chrom_cov_freq: IndexMap<&str,  BTreeMap<usize, usize>> = IndexMap::with_capacity(HG19_CHROMS.len());
    // histogram of major allele frequency per chromosome
    // key1: chrom, key2: major allele frequency as usize where 1.0 is 1000, value: counter
    let mut chrom_maj_freq: IndexMap<&str,  BTreeMap<usize, usize>> = IndexMap::with_capacity(HG19_CHROMS.len());
    // histogram of forward-reverse ratio frequency per chromosome
    // key1: chrom, key2: FR ratio as usize where 1.0 is 1000, value: counter
    let mut chrom_frratio_freq: IndexMap<&str,  BTreeMap<usize, usize>> = IndexMap::with_capacity(HG19_CHROMS.len());
    // Total number of sites visited per chrom
    let mut chrom_totals: IndexMap<&str, usize> = IndexMap::with_capacity(HG19_CHROMS.len());


    // Scoped threading
    // channels
    let (send_rec, recv_rec) = bounded(1);

    let (send_site, recv_site) = bounded(1);
    // let (send_cov, recv_cov) = bounded(1);
    // let (send_fwdratio, recv_fwdratio) = bounded(1);
    // let (send_majfreq, recv_majfreq) = bounded(1);
    // let (send_het, recv_het) = bounded(1);

    thread::scope(|s| {
        // Producer thread
        s.spawn(|_| {
            let mut c: usize = 0;
            for &chrom in HG19_CHROMS.iter() {
                let ul: u64 = (HG19_CHROM_LENS.get(chrom).unwrap().to_owned() as u64 / FETCH_CHUNKSIZE) + 1;
                for i in 0..ul {
                    tbx_reader.fetch(*chrom_tid_lookup.get(chrom).unwrap(), i*FETCH_CHUNKSIZE, (i+1)*FETCH_CHUNKSIZE).unwrap();
                    for record in tbx_reader.records() {
                        send_rec.send((c, record.unwrap())).unwrap();
                        c += 1;
                    }
                }
            }
            drop(send_rec)
        });

        // Parallel processing by n threads
        for _ in 0..threads {
            // Send to sink, receive from source
            let (sendr, recvr) = (send_site.clone(), recv_rec.clone());
            // Spawn workers in separate threads
            s.spawn(move |_| {
                // Receive until channel closes
                for (i, bytes) in recvr.iter() {
                    let (chrom, pos, ref_char, raw_pileup) = bytes_to_pileup(bytes);
                    // quality filter bases based on bq and mq
                    let filt_pileup = qual_filter_pileup(&raw_pileup, min_bq, min_mq);

                    // test if cov after filtering is >= min_cov
                    // skip if below min_cov
                    let filt_cov = filt_pileup.cov();
                    if filt_cov < min_cov { continue }
                    let mut site_flags = ControlFilterScore::PassedMinCov;

                    // Convert to FullBaseCount
                    let filt_fbc = filt_pileup.full_base_count();
                    // compute forward ratio
                    let fwd_ratio = (filt_fbc.forward() as f64) / (filt_cov as f64);
                    if (fwd_ratio >= min_fratio) && (fwd_ratio < max_fratio) {
                        site_flags.insert(ControlFilterScore::PassedFRRatio);
                    }
                    // compute major allele
                    // TODO: Code is naive and cannot handle het sites
                    let allele_set = AlleleSet::from_fullbasecount(&filt_fbc);
                    let maj_base = allele_set.alleles[0].base();
                    let maj_cnt = allele_set.alleles[0].count();
                    let maj_freq = (maj_cnt as f64) / (filt_cov as f64);
                    // check if invariant
                    if maj_freq == 1.0 {
                        site_flags.insert(ControlFilterScore::InvariantSite);
                    }
                    // check if ref_char == maj_base
                    if ref_char == *maj_base {
                        site_flags.insert(ControlFilterScore::ReferenceBase);
                    }
                    sendr.send((i, chrom, pos, filt_pileup, site_flags)).unwrap();
                }
            });
        }
        drop(send_site);

        // Sink
        for (i, chrom, pos, data, site_flags) in recv_site.iter() {
            println!("{} {}\t{}\t{:?}", i, chrom, pos, site_flags);
        }
    }).unwrap();


    // // TODO: Output to cov freq file
    // // Print stats
    // let cov_freq_path = format!("{}{}", control_path.rsplit_once(".").unwrap().0, ".cov.txt");
    // let f = File::create(cov_freq_path).expect("Unable to create file");
    // let mut f = BufWriter::new(f);
    // f.write_fmt(format_args!("# Coverage\n")).expect("Unable to write data");
    // hist_to_file(&chrom_cov_freq, &chrom_totals, f);

    // // TODO: Output to allele freq file
    // // Print stats
    // let allele_freq_path = format!("{}{}", control_path.rsplit_once(".").unwrap().0, ".af.txt");
    // let f = File::create(allele_freq_path).expect("Unable to create file");
    // let mut f = BufWriter::new(f);
    // f.write_fmt(format_args!("# Major allele frequency\n")).expect("Unable to write data");
    // fhist_to_file(&chrom_maj_freq, &chrom_totals, f, FREQ_DIVISOR);

    // // TODO: Output to fr ratio file
    // // Print stats
    // let frratio_path = format!("{}{}", control_path.rsplit_once(".").unwrap().0, ".fr.txt");
    // let f = File::create(frratio_path).expect("Unable to create file");
    // let mut f = BufWriter::new(f);
    // f.write_fmt(format_args!("# FR ratio frequency\n")).expect("Unable to write data");
    // fhist_to_file(&chrom_frratio_freq, &chrom_totals, f, FREQ_DIVISOR);

    // // Output as bed file
    // let bed_path = format!("{}{}", control_path.rsplit_once(".").unwrap().0, ".passed.bed");
    // let f = File::create(bed_path).expect("Unable to create file");
    // pos_to_bed(&passed_pos, BufWriter::new(f));
    
    // Evaluate control first. Do bq and mq filtering and test coverage
    // Evaluate samples. Do bq and mq for each sample and test coverage for each as well
    // If passed control and passed ALL samples, then include in the output

    // Read all the mpileups and get the sites which is found in all
    // Saturate processing from queue based on the number of threads
    // let (tx, rx) = mpsc::channel();
    // for i in 0..threads {
    //     let ttx = tx.clone();
    //     thread::spawn(move || {
    //         ttx.send(i).unwrap();
    //     });
    // }
    // Close send to finish
    // drop(tx);
    // for received in rx {
    //     println!("Got: {}", received);
    // }

}