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
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::collections::BTreeMap;

use clap::{Arg, App};
use rust_htslib::tbx::{self, Read as TbxRead};
use crossbeam::channel::bounded;
use crossbeam::thread;
// use indexmap::IndexMap;

use htsops::{BaseCount, AlleleSet};
use htsops::util::{validate_path, validate_qual, read_tabix, bytes_to_pileup};
use htsops::pileup::{SiteInfo, SitePileup, PileupStats};
use htsops::constant::*;


pub fn bytes_to_pileupstats(b: Vec<u8>) -> Option<(SiteInfo, PileupStats)> {
    // Convert binary mpileup line and parse
    let record_string = String::from_utf8(b).unwrap();
    let data_vec: Vec<&str> = record_string
        .split("\t\t")
        .collect();
    let info = SiteInfo::from_str(data_vec[0], "\t");
    let stats = match PileupStats::from_str(data_vec[1], "\t") {
        Some(stats) => stats,
        None => return None
    };
    Some((info, stats))
}


fn compute_per_site_sample_stats(info_stats_vec: &Vec<(SiteInfo, PileupStats)>, tbx_readers: &mut Vec<tbx::Reader>, tid: u64, num_samples: usize, min_cov: usize, min_bq: u8, min_mq: u8) -> BTreeMap<u64, Vec<PileupStats>> {
    let start_pos = info_stats_vec[0].0.pos();
    let end_pos = info_stats_vec[info_stats_vec.len()-1].0.pos();
    // Map of sites to look at and a vec to hold sample data at that site
    let mut per_site_sample_stats: BTreeMap<u64, Vec<PileupStats>> = info_stats_vec.iter()
        .map(|(info, _)| (info.pos(), Vec::with_capacity(num_samples)))
        .collect();
    for tbx_reader in tbx_readers.iter_mut() {
        tbx_reader.fetch(tid, start_pos, end_pos).unwrap();
        for b in tbx_reader.records() {
            // skip if the current site in this sample is not found in the set of sites to look at
            let (info, pileup) = bytes_to_pileup(b.unwrap());
            let sample_stats_vec = match per_site_sample_stats.get_mut(&info.pos()) {
                Some(vec) => vec,
                None => continue,
            };
            // if !stats.site_flags.contains(filter_threshold) { continue }
            // base and mapping quality filtering
            let filt_pileup = pileup.quality_filter(min_bq, min_mq);
            // skip if current sample is below min_cov after filtering
            if filt_pileup.cov() < min_cov { continue }
            // compute stats 
            let mut stats = match PileupStats::from_pileup(&filt_pileup) {
                Some(stats) => stats,
                None => continue,
            };
            stats.compute_stats();
            sample_stats_vec.push(stats);  
        };
    }
    per_site_sample_stats
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
            .help("Minimum coverage depth for each target sample. Sites where the coverage depth is less than this value will not be used.")
            .takes_value(true)
            .validator(validate_qual))
        .arg(Arg::with_name("min_bq")
            .long("min-bq")
            .default_value("20")
            .value_name("INT")
            .help("Minimum base quality for target samples. Bases whose base quality is less than this value will not be counted.")
            .takes_value(true)
            .validator(validate_qual))
        .arg(Arg::with_name("min_mq")
            .long("min-mq")
            .default_value("40")
            .value_name("INT")
            .help("Minimum mapping quality for target samples. Bases belonging to reads with mapping quality less than this value will not be counted.")
            .takes_value(true)
            .validator(validate_qual))
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
        .arg(Arg::with_name("pthreads")
            .short("p")
            .long("pthreads")
            .default_value("1")
            .value_name("INT")
            .help("Number of producer threads to use")
            .takes_value(true))
        .arg(Arg::with_name("cthreads")
            .short("c")
            .long("cthreads")
            .default_value("4")
            .value_name("INT")
            .help("Number of consumer threads to use")
            .takes_value(true))
        .arg(Arg::with_name("CTRL_FILE")
            .help("Output file from filter_control")
            .required(true)
            .index(1)
            .validator(validate_path))
        .arg(Arg::with_name("PILEUP_LIST")
            .help("Contains list of pileup files to process. Assumes control pileup is the first in the list.")
            .required(true)
            .index(2)
            .validator(validate_path))
        .get_matches();
    let min_cov: usize = matches.value_of("min_cov").unwrap().parse::<usize>().unwrap();
    let min_bq: u8 = matches.value_of("min_bq").unwrap().parse::<u8>().unwrap();
    let min_mq: u8 = matches.value_of("min_mq").unwrap().parse::<u8>().unwrap();
    let min_fratio: f64 = matches.value_of("min_fratio").unwrap().parse::<f64>().unwrap();
    let max_fratio: f64 = matches.value_of("max_fratio").unwrap().parse::<f64>().unwrap();
    let pthreads: usize = matches.value_of("pthreads").unwrap().parse::<usize>().unwrap();
    let cthreads: usize = matches.value_of("cthreads").unwrap().parse::<usize>().unwrap();
    let ctrl_filt_path = matches.value_of("CTRL_FILE").unwrap();
    let pileup_path_file: &str = matches.value_of("PILEUP_LIST").unwrap();

    let target_paths: Vec<String> = {
        let file = File::open(pileup_path_file).unwrap();
        let mut pileup_paths: Vec<String> = BufReader::new(file).lines()
            .filter_map(|l| {
                if let Ok(line) = l { 
                    if line.len() > 0 { Some(line) } else { None }
                } 
                else { None }
            })
            .collect();
        // Ignore control pileup
        pileup_paths.remove(0);
        pileup_paths
    };
    let target_paths: Vec<&str> = target_paths.iter().map(|s| s.as_ref()).collect();

    let (mut ctrl_tbx_reader, chrom_tid_lookup) = read_tabix(ctrl_filt_path);
    let num_samples: usize = target_paths.len();
    // let filter_threshold = ControlFilterScore::PassedMinCov | ControlFilterScore::PassedFRRatio;

    // Read chrrom and pos of control filter file line by line to find intervals
    // Everytime continuity is broken send to channel

    const INTERVAL_SIZE: usize = 1000;

    for &chrom in HG19_CHROMS.iter() {
        // let mut collated_data: BTreeMap<usize, Vec<(SiteInfo, SitePileup)>> = BTreeMap::new();

        let (send_interval, recv_interval) = bounded(cthreads);
        let (send_site, recv_site) = bounded(cthreads);

        // Producer thread
        // Make intervals of N length from filtered control file
        thread::scope(|s| {
            s.spawn(|_| {
                ctrl_tbx_reader.fetch(*chrom_tid_lookup.get(chrom).unwrap(), 0, *HG19_CHROM_LEN.get(chrom).unwrap()).unwrap();
                let mut info_stats_vec: Vec<(SiteInfo, PileupStats)> = Vec::with_capacity(INTERVAL_SIZE);
                for b in ctrl_tbx_reader.records() {
                    let (site_info, pileup_stats) = match bytes_to_pileupstats(b.unwrap()) {
                        Some(data) => data,
                        None => continue,
                    };
                    let current_pos = site_info.pos();
                    info_stats_vec.push((site_info, pileup_stats));
                    // Send if there are N records
                    if info_stats_vec.len() == INTERVAL_SIZE {
                        let mut v: Vec<(SiteInfo, PileupStats)> = Vec::with_capacity(INTERVAL_SIZE);
                        v.append(&mut info_stats_vec);
                        send_interval.send(v).unwrap();
                    }
                }
                drop(send_interval);
            });

            // Parallel processing by n threads of intervals
            for _ in 0..cthreads {
                let (sendr, recvr) = (send_site.clone(), recv_interval.clone());
                let tid: u64 = *chrom_tid_lookup.get(&chrom as &str).unwrap();
                // Open multiple tbx readers
                let mut tbx_readers: Vec<tbx::Reader> = target_paths.iter()
                    .map(|path| {
                        // let split: Vec<&str> = path.split('/').collect();
                        // let filename = split[split.len() - 1];
                        let (tbx_reader, _) = read_tabix(path);
                        tbx_reader
                    })
                    .collect();
                // Spawn workers in separate threads
                s.spawn(move |_| {
                    // Receive position vecs until channel closes
                    for info_stats_vec in recvr.iter() {
                        // Map of sites to look at and a vec to hold sample data at that site
                        let per_site_sample_stats: BTreeMap<u64, Vec<PileupStats>> = compute_per_site_sample_stats(&info_stats_vec, &mut tbx_readers, tid, num_samples, min_cov, min_bq, min_mq);
                        // Convert received vector into a map with position as key
                        let ctrl_info_stats: BTreeMap<u64, (SiteInfo, PileupStats)> = info_stats_vec.into_iter()
                            .map(|(info, pileup)| (info.pos(), (info, pileup)))
                            .collect();
                        // Send only sites present in all samples and has a variant
                        'per_site: for (pos, sample_stats_vec) in per_site_sample_stats.into_iter() {
                            if sample_stats_vec.len() != num_samples { continue }
                            let (info, ctrl_stats) = ctrl_info_stats.get(&pos).unwrap();
                            // major allele is the ctrl base
                            let ctrl_major_allele = ctrl_stats.major_allele();
                            // minor allele based on pooled count
                            let pooled_base_count: BaseCount = sample_stats_vec.iter()
                                .fold(BaseCount::empty(), |acc, current| acc + current.base_count() );
                            let pooled_allele_set = match AlleleSet::from_base_count(pooled_base_count) {
                                Some(allele_set) => allele_set,
                                None => continue 'per_site,
                            };
                            let pooled_major_allele = pooled_allele_set.major_allele();
                            let pooled_minor_allele = pooled_allele_set.minor_allele();
                            // control and pooled major allele should be the same or else skip
                            if ctrl_major_allele != pooled_major_allele { continue 'per_site }
                            // skip invariant sites
                            if pooled_minor_allele.is_none() { continue 'per_site }
                            // Check individual samples if they are consistent with the pooled alleles
                            for stats in sample_stats_vec.iter() {
                                // Skip samples with only 1 allele
                                if stats.num_alleles() < 2 { continue }
                                // Check if pooled minor allele is the minor allele for all samples
                                // if minor allele is present but it is not the minor allele in the current sample
                                // if (pileup_fbc.count_of(&pooled_minor_allele) > 0) && (*sample_allele_set.alleles[1].base() != pooled_minor_allele) {
                                //     continue 'per_site;
                                // }

                                
                                // if (pileup_fbc.full_count_of(&pooled_minor_allele.to_ascii_uppercase()) < 1) || (pileup_fbc.full_count_of(&pooled_minor_allele.to_ascii_lowercase()) < 1) {
                                //     continue 'per_site;
                                // }
                            }
                            let pooled_stats = PileupStats::from_allele_set(pooled_allele_set);
                            // pooled must have at least 0.1% of depth
                            if pooled_stats.minor_freq() < 0.001 { continue 'per_site }
                            sendr.send((pos, info.clone(), ctrl_stats.clone(), pooled_stats, sample_stats_vec)).unwrap();
                        }                        
                    }
                });
            }
            drop(send_site);

            // Sink
            // for (c, site_pileups) in recv_site.iter() {
            //     collated_data.insert(c, site_pileups);
            // }
            s.spawn(move |_| {
                for (_, info, ctrl_stats, pooled_stats, sample_stats_vec) in recv_site.iter() {
                    let concat_sample_stats: String = sample_stats_vec.iter()
                        .map(|stats| format!("{}", stats))
                        .collect::<Vec<String>>()
                        .join("\t\t");
                    // calculate fisher exact test
                    println!("{}\t\t{}\t\t{}\t\t{}", 
                        info, ctrl_stats, pooled_stats, concat_sample_stats);
                }
            });

        }).unwrap();

        
        // Output
        // for (c, site_pileups) in collated_data.iter() {
        //     let (site_info, site_pileup1) = &site_pileups[0];
        //     let (_, site_pileup2) = &site_pileups[1];
        //     let (_, site_pileup3) = &site_pileups[2];
        //     println!("{}\t{}\t{}\t{}\t{}", c, site_info, site_pileup1, site_pileup2, site_pileup3);
        // }
        
    }    
}