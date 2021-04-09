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
use std::io::{self, BufWriter, Write};
use std::fmt;

use clap::{Arg, App};
use rust_htslib::tbx::{self, Read as TbxRead};
use indexmap::IndexMap;
use crossbeam::channel::bounded;
use crossbeam::thread;
use bgzip::BGZFWriter;
use flate2;
use std::process::Command;

use htsops::util::{validate_path, validate_qual, read_tabix};
use htsops::pileup::{SitePileup, FullBaseCount, AlleleSet};
use htsops::filter::ControlFilterScore;
use htsops::constant::*;


/// Convert raw bytes into a Pileup object
fn bytes_to_pileup(b: Vec<u8>) -> (SiteInfo, SitePileup) {
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

    (SiteInfo{ chrom, pos, ref_char },
     SitePileup::from_str(
        sample_name,
        &ref_char,
        (&cols[3]).parse::<usize>().unwrap(), // cov
        &cols[4], &cols[5], &cols[6])         // base_str, bq_str, mq_str
    )
}

/// Filter pileup based on base and mapping quality scores
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

#[derive(Debug)]
pub struct SiteInfo {
    pub chrom: String,
    pub pos: u64,
    pub ref_char: char,
}
impl fmt::Display for SiteInfo {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}\t{}\t{}", self.chrom, self.pos, self.ref_char)
    }
}

#[derive(Debug)]
struct PileupStats {
    full_base_count: FullBaseCount,
    major_allele: char,
    major_freq: f64,
    minor_allele: Option<char>,
    minor_freq: f64,
    site_flags: ControlFilterScore,
    cov: usize,
}
impl PileupStats {
    fn from_pileup(pileup: &SitePileup, ref_char: &char, min_fratio: f64, max_fratio: f64) -> PileupStats {
        let cov = pileup.cov();
        let mut site_flags = ControlFilterScore::PassedMinCov;
        // Convert to FullBaseCount
        let fbc = pileup.full_base_count();
        // compute forward ratio
        let fwd_ratio = (fbc.forward() as f64) / (cov as f64);
        if (fwd_ratio >= min_fratio) && (fwd_ratio < max_fratio) {
            site_flags.insert(ControlFilterScore::PassedFRRatio);
        }
        // compute major allele
        // TODO: Code is naive and cannot handle het sites
        let allele_set = AlleleSet::from_fullbasecount(&fbc);
        let maj_base = allele_set.alleles[0].base();
        let maj_cnt = allele_set.alleles[0].count();
        let maj_freq = (maj_cnt as f64) / (cov as f64);
        // check if invariant
        if maj_freq == 1.0 {
            site_flags.insert(ControlFilterScore::InvariantSite);
        }
        // check if ref_char == maj_base
        if ref_char == maj_base {
            site_flags.insert(ControlFilterScore::ReferenceBase);
        }

        PileupStats {
            cov: fbc.total(),
            full_base_count: fbc,
            major_allele: *maj_base,
            major_freq: maj_freq,
            minor_allele: None,
            minor_freq: 0.0,
            site_flags: site_flags,
        }
    }
    fn major_fratio(&self) -> f64 {
        let f_cnt = self.full_base_count.full_count_of(&self.major_allele.to_ascii_uppercase());
        let fr_cnt = self.full_base_count.count_of(&self.major_allele);
        (f_cnt as f64) / (fr_cnt as f64)
    }
    fn minor_fratio(&self) -> f64 {
        if self.minor_freq == 0.0 { return 0.0 }
        let minor_allele = self.minor_allele.unwrap();
        let f_cnt = self.full_base_count.full_count_of(&minor_allele.to_ascii_uppercase());
        let fr_cnt = self.full_base_count.count_of(&minor_allele);
        (f_cnt as f64) / (fr_cnt as f64)
    }
}
impl fmt::Display for PileupStats {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{site_flags}\t{majorb}:{majorf:.3}\t{minorb}:{minorf:.3}\t{fbc}\t{cov}", 
            site_flags=self.site_flags.to_bits(),
            majorb=self.major_allele,
            majorf=self.major_freq,
            minorb={if let Some(b) = self.minor_allele { b } else { '-' }},
            minorf=self.minor_freq,
            fbc=self.full_base_count,
            cov=self.cov,
        )
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
    let min_cov: usize = matches.value_of("min_cov").unwrap().parse::<usize>().unwrap();
    let min_bq: u8 = matches.value_of("min_bq").unwrap().parse::<u8>().unwrap();
    let min_mq: u8 = matches.value_of("min_mq").unwrap().parse::<u8>().unwrap();
    let min_fratio: f64 = matches.value_of("min_fratio").unwrap().parse::<f64>().unwrap();
    let max_fratio: f64 = matches.value_of("max_fratio").unwrap().parse::<f64>().unwrap();
    let output_path: Option<&str> = matches.value_of("out");
    let bgzip: bool = matches.is_present("bgzip");
    let threads: usize = matches.value_of("threads").unwrap().parse::<usize>().unwrap();
    // let include_bad_cov: bool = matches.is_present("include_bad_cov");
    let control_path: &str = matches.value_of("CONTROL_PILEUP").unwrap();

    // Open control mpileup
    let (mut tbx_reader, chrom_tid_lookup) = read_tabix(control_path);
    let filter_threshold = ControlFilterScore::PassedMinCov | ControlFilterScore::PassedFRRatio | ControlFilterScore::InvariantSite;

    // Scoped threading
    // channels
    let (send_rec, recv_rec) = bounded(1);
    let (send_site, recv_site) = bounded(1);
    // let (send_cov, recv_cov) = bounded(1);
    // let (send_fwdratio, recv_fwdratio) = bounded(1);
    // let (send_majfreq, recv_majfreq) = bounded(1);
    // let (send_het, recv_het) = bounded(1);

    let mut collated_data: BTreeMap<usize, (SiteInfo, SitePileup, PileupStats)> = BTreeMap::new();
    let mut chrom_covs: HashMap<String, Vec<usize>> = HashMap::new();

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
                    let (site_info, raw_pileup) = bytes_to_pileup(bytes);
                    // quality filter bases based on bq and mq
                    let filt_pileup = qual_filter_pileup(&raw_pileup, min_bq, min_mq);
                    // test if cov after filtering is >= min_cov
                    // skip if below min_cov
                    if filt_pileup.cov() < min_cov { continue }
                    // Generate stats
                    let site_stats = PileupStats::from_pileup(&filt_pileup, &site_info.ref_char, min_fratio, max_fratio);
                    // Send to channel
                    if site_stats.site_flags.contains(filter_threshold) {
                        sendr.send((i, site_info, filt_pileup, site_stats)).unwrap();
                    }
                }
            });
        }
        drop(send_site);

        // Sink
        for (i, info, pileup, stats) in recv_site.iter() {
            chrom_covs.entry(info.chrom.clone()).or_insert(Vec::new()).push(stats.full_base_count.total());
            collated_data.insert(i, (info, pileup, stats));
        }
    }).unwrap();

    // compute stats
    let mut chrom_mean_cov: HashMap<&str, f64> = HashMap::new();
    let mut chrom_cumdist: HashMap<&str, BTreeMap<usize, f64>> = HashMap::new();
    let (send_cov, recv_cov) = bounded(1);
    let (send_calc, recv_calc) = bounded(1);
    thread::scope(|s| {
        // Producer thread
        s.spawn(|_| {
            chrom_covs.iter().for_each(|(chrom, covs)| send_cov.send((chrom, covs)).unwrap() );
            drop(send_cov)
        });
        // Parallel processing by n threads
        for _ in 0..threads {
            // Send to sink, receive from source
            let (sendr, recvr)= (send_calc.clone(), recv_cov.clone());
            // Spawn workers in separate threads
            s.spawn(move |_| {
                // Receive until channel closes
                for (chrom, covs) in recvr.iter() {
                    let mean = (covs.iter().sum::<usize>() as f64) / (covs.len() as f64);
                    let mut btm: BTreeMap<usize, usize> = BTreeMap::new();
                    covs.iter().for_each(|cov| *btm.entry(*cov).or_insert(0) += 1 );
                    let pctdist: BTreeMap<usize, f64> = btm.into_iter()
                        .map(|(cov, cnt)| (cov, (cnt as f64) / (covs.len() as f64)))
                        .collect();
                    let mut cum_sum = 0.0;
                    let cumdist: BTreeMap<usize, f64> = pctdist.into_iter()
                        .map(|(cov, pct)| {
                            cum_sum += pct;
                            (cov, cum_sum)
                        })
                        .collect();
                    sendr.send((chrom, mean, cumdist)).unwrap();
                }
            });
        }
        drop(send_calc);
        for (chrom, mean, cumdist) in recv_calc.iter() {
            chrom_mean_cov.insert(chrom, mean);
            chrom_cumdist.insert(chrom, cumdist);
        }
    }).unwrap();

    // Output either to stdout or to file
    let stdout = io::stdout();
    let mut write_h: Box<dyn Write> = match output_path {
        Some(path) => {
            let f = File::create(path).expect("Unable to create file");
            let writer = BufWriter::new(f);
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
    for (_, (info, _, stats)) in collated_data.iter() {
        let chrom: &str = info.chrom.as_ref();
        let cov_pctile = chrom_cumdist.get(&chrom).unwrap().get(&stats.cov).unwrap();
        write_h.write_fmt(format_args!("{}\t{}\t{:.4}\n", info, stats, cov_pctile)).unwrap();
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