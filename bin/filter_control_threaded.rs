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
use std::collections::{HashMap, BTreeMap};
use std::fs::File;
use std::io::{self, BufWriter, Write};
// use std::fmt;

use clap::{Arg, App};
use rust_htslib::tbx::Read;
use crossbeam::channel::bounded;
use crossbeam::thread;
use bgzip::BGZFWriter;
use flate2;
use std::process::Command;
use histogram::Histogram;
use rand::prelude::*;
use rand::{SeedableRng, thread_rng};
use rand::rngs::SmallRng;

use htsops::util::{validate_path, validate_qual, read_tabix, bytes_to_pileup};
use htsops::pileup::{SiteInfo, SitePileup, PileupStats};
use htsops::filter::{ControlBitscore, ControlScoreParams};
use htsops::constant::*;


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
            .takes_value(true))
        .arg(Arg::with_name("cov_threshold")
            .long("cov-th")
            .default_value("0")
            .value_name("INT")
            .help("Sites above the coverage threshold will be marked suspicious. If value is 0, the coverage threshold will be the 95th percentile per chromosome coverage distribution")
            .takes_value(true))
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
        .arg(Arg::with_name("max_hom_freq")
            .long("max-hom-freq")
            .default_value("0.3")
            .value_name("FLOAT")
            .help("Maximum minor allele frequency for a site to be considered homozygous")
            .takes_value(true))
        .arg(Arg::with_name("min_het_freq")
            .long("min-het-freq")
            .default_value("0.67")
            .value_name("FLOAT")
            .help("Minimum minor allele frequenct for a site to be considered heterozygous")
            .takes_value(true))
        .arg(Arg::with_name("out")
            .short("o")
            .long("out")
            .value_name("PATH")
            .help("Set output path. If blank or not specified, output is printed to stdout")
            .takes_value(true))
        .arg(Arg::with_name("bgzip")
            .short("z")
            .long("bgzip")
            .help("Compress output using bgzip"))
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
        .arg(Arg::with_name("CONTROL_PILEUP")
            .help("Control mpileup path")
            .required(true)
            .index(1)
            .validator(validate_path))
        .get_matches();

    // Assign to variables
    // control filtering params
    let min_cov: usize = matches.value_of("min_cov").unwrap().parse::<usize>().unwrap();
    let cov_threshold: usize = matches.value_of("cov_threshold").unwrap().parse::<usize>().unwrap();
    let min_bq: u8 = matches.value_of("min_bq").unwrap().parse::<u8>().unwrap();
    let min_mq: u8 = matches.value_of("min_mq").unwrap().parse::<u8>().unwrap();
    let min_fratio: f64 = matches.value_of("min_fratio").unwrap().parse::<f64>().unwrap();
    let max_fratio: f64 = matches.value_of("max_fratio").unwrap().parse::<f64>().unwrap();
    let max_hom_freq: f64 = matches.value_of("max_hom_freq").unwrap().parse::<f64>().unwrap();
    let min_het_freq: f64 = matches.value_of("min_het_freq").unwrap().parse::<f64>().unwrap();
    let score_params = ControlScoreParams::new(
        min_cov,
        cov_threshold,
        (min_fratio, max_fratio),
        max_hom_freq,
        min_het_freq, 
    );

    let output_path: Option<&str> = matches.value_of("out");
    let bgzip: bool = matches.is_present("bgzip");
    let threads: usize = matches.value_of("threads").unwrap().parse::<usize>().unwrap();
    // let include_bad_cov: bool = matches.is_present("include_bad_cov");
    let control_path: &str = matches.value_of("CONTROL_PILEUP").unwrap();

    // Open control mpileup
    let (mut tbx_reader, chrom_tid_lookup) = read_tabix(control_path);

    // Scoped threading
    // channels
    let (send_rec, recv_rec) = bounded(threads);
    let (send_site, recv_site) = bounded(threads);

    let mut collated_data: BTreeMap<usize, (SiteInfo, SitePileup, PileupStats, ControlBitscore)> = BTreeMap::new();
    let mut chrom_covs: HashMap<String, Vec<usize>> = HashMap::new();

    thread::scope(|s| {
        // Producer thread
        s.spawn(|_| {
            let mut c: usize = 0;
            for &chrom in HG19_CHROMS.iter() {
                tbx_reader.fetch(*chrom_tid_lookup.get(chrom).unwrap(), 0, *HG19_CHROM_LEN.get(chrom).unwrap()).unwrap();
                for record in tbx_reader.records() {
                    send_rec.send((c, record.unwrap())).unwrap();
                    c += 1;
                }
            }
            drop(send_rec)
        });

        // Parallel processing by n threads
        for _ in 0..threads {
            // Send to sink, receive from source
            let (sendr, recvr) = (send_site.clone(), recv_rec.clone());
            let score_params = score_params.clone();
            // Spawn workers in separate threads
            s.spawn(move |_| {
                // Receive until channel closes
                for (i, bytes) in recvr.iter() {
                    let (site_info, raw_pileup) = bytes_to_pileup(bytes);
                    // quality filter bases based on bq and mq
                    let filt_pileup = raw_pileup.quality_filter(min_bq, min_mq);
                    // test if cov after filtering is >= min_cov
                    // skip if below min_cov
                    if filt_pileup.cov() < min_cov { continue }
                    // generate stats
                    let mut site_stats = match PileupStats::from_pileup(&filt_pileup) {
                        Some(stats) => stats,
                        None => continue,
                    };
                    // fisher test
                    // odds ratio of minor f/r / major f/r, 0.5/alt_arm_total zero adjustment
                    site_stats.compute_stats();

                    // compute bitscore for min_cov, strand_bias, refseq_variant and genotype
                    let mut bitscore = ControlBitscore::empty();
                    bitscore.score_min_cov(&site_stats, score_params.min_cov());
                    bitscore.score_stand_bias(&site_stats, score_params.min_fratio(), score_params.max_fratio());
                    bitscore.score_refseq_variant(&site_stats, site_info.ref_base());
                    bitscore.score_genotype(&site_stats, score_params.max_hom_freq(), score_params.min_het_freq());
                    // compute bitscore for cov_threshold if not 0
                    if score_params.cov_threshold() != 0 {
                        bitscore.score_cov_threshold(&site_stats, score_params.cov_threshold());
                    }
                    sendr.send((i, site_info, filt_pileup, site_stats, bitscore)).unwrap();
                }
            });
        }
        drop(send_site);

        // Sink
        for (i, info, pileup, stats, bitscore) in recv_site.iter() {
            chrom_covs.entry(info.chrom().to_owned()).or_insert(Vec::new()).push(stats.cov());
            collated_data.insert(i, (info, pileup, stats, bitscore));
        }        
    }).unwrap();

    // compute covearage stats
    let mut chrom_mean_cov: HashMap<&str, f64> = HashMap::new();
    let mut chrom_cumdist: HashMap<&str, BTreeMap<usize, f64>> = HashMap::new();
    let mut chrom_cov_threshold: HashMap<&str, usize> = HashMap::new();
    let (send_cov, recv_cov) = bounded(1);
    let (send_calc, recv_calc) = bounded(1);
    thread::scope(|s| {
        // producer thread
        // send vec of site coverages, 1 for each chromosome
        s.spawn(|_| {
            chrom_covs.iter().for_each(|(chrom, covs)| send_cov.send((chrom, covs)).unwrap() );
            drop(send_cov)
        });
        // parallel processing by n threads
        for _ in 0..threads {
            // send to sink, receive from source
            let (sendr, recvr)= (send_calc.clone(), recv_cov.clone());
            let score_params = score_params.clone();
            // spawn workers in separate threads
            s.spawn(move |_| {
                // receive until channel closes
                // calculate the mean cov per chromosome
                // calculate the cumulative distrivution of cov per chromosome
                for (chrom, covs) in recvr.iter() {
                    let mean = (covs.iter().sum::<usize>() as f64) / (covs.len() as f64);
                    // create an cov freq btree map ordered by coverage val
                    let mut cov_cnt: BTreeMap<usize, usize> = BTreeMap::new();
                    covs.iter().for_each(|cov| *cov_cnt.entry(*cov).or_insert(0) += 1 );
                    // create a cummulative frequency btree map
                    let mut cum_sum = 0.0;
                    let cov_cummfreq: BTreeMap<usize, f64> = cov_cnt.into_iter()
                        .map(|(cov, cnt)| {
                            let pct = (cnt as f64) / (covs.len() as f64);
                            cum_sum += pct;
                            (cov, cum_sum)
                        })
                        .collect();

                    // compute 95th percentile of coverage and use as cov_threshold
                    let per_chrom_cov_threshold: Option<usize> = match score_params.cov_threshold() {
                        0 => {
                            // bootstrap 100 times
                            let mut thread_rng = thread_rng();
                            let mut bs_cov_hist = Histogram::new();
                            for _ in 0..1000 {
                                let mut hist = Histogram::new();
                                let mut rng = SmallRng::from_rng(&mut thread_rng).unwrap();
                                // randomly pick items with replacement from covs
                                let vec_len = covs.len();
                                for _ in 0..vec_len {
                                    hist.increment(covs[rng.gen_range(0..vec_len)] as u64).unwrap();
                                } 
                                // get 99th percentile for this bootstrap
                                let threshold = hist.percentile(99.0).unwrap();
                                bs_cov_hist.increment(threshold).unwrap();
                            }
                            // get 95th percentile of 99th percentile dist
                            Some(bs_cov_hist.percentile(95.0).unwrap() as usize)
                        },
                        _ => None,
                    };
                    
                    sendr.send((chrom, mean, cov_cummfreq, per_chrom_cov_threshold)).unwrap();
                }
            });
        }
        drop(send_calc);
        for (chrom, mean, cumdist, per_chrom_cov_threshold) in recv_calc.iter() {
            chrom_mean_cov.insert(chrom, mean);
            chrom_cumdist.insert(chrom, cumdist);
            if let Some(cov_threshold) = per_chrom_cov_threshold {
                chrom_cov_threshold.insert(chrom, cov_threshold);
            }
        }
    }).unwrap();

    // Output either to stdout or to file
    let stdout = io::stdout();
    let mut write_h: Box<dyn Write> = match output_path {
        Some(path) => {
            let writer = BufWriter::new(File::create(path).expect("Unable to create file"));
            if bgzip {
                Box::new(BGZFWriter::new(writer, flate2::Compression::default()))
            } else {
                Box::new(writer)
            }
        },
        None => {
            Box::new(stdout.lock())
        }
    };
    // Write to writer
    for (_, (info, _, stats, mut bitscore)) in collated_data.iter() {
        let chrom: &str = info.chrom().as_ref();
        let cov_pctile = chrom_cumdist.get(&chrom).unwrap().get(&stats.cov()).unwrap();
        if let Some(cov_threshold) = chrom_cov_threshold.get(&chrom) {
            bitscore.score_cov_threshold(&stats, *cov_threshold);
        }
        // write_h.write_fmt(format_args!("{}\t{}\t{:.4}\n", info, stats, cov_pctile)).unwrap();
        let pval = stats.allele_orientation_fisher_test().unwrap();
        let or = stats.allele_orientation_or().unwrap();
        write_h.write_fmt(format_args!("{}\t\t{}\t{:.4}\t{:.4}\t\t{:.4}\n", info, stats, pval, or, cov_pctile)).unwrap();
    }
    write_h.flush().unwrap();

    // Index if bgzipped
    if let Some(path) = output_path {
        if bgzip {
            let status = Command::new("tabix")
                    .args(&["-s", "1"])
                    .args(&["-b", "2"])
                    .args(&["-e", "2"])
                    .arg(&path)
                    .output()
                    .expect(&format!("error executing tabix for {}", &path))
                    .status;
            if status.success() {
                println!("tabix  : {} {}", &path, status.code().take().unwrap());
            } else {
                match status.code() {
                    Some(code) => println!("ERROR tabix  : {} {}", &path, code),
                    None       => println!("tabix terminated by signal")
                }
            }
        }
    }
}