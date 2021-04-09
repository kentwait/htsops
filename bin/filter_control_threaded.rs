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
use std::fmt;

use clap::{Arg, App};
use rust_htslib::tbx::{self, Read as TbxRead};
use indexmap::IndexMap;
use crossbeam::channel::{self, select, bounded};
use crossbeam::thread;

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
        write!(f, "{site_flags}\t{majorb}:{majorf}\t{minorb}:{minorf}\t{fbc}", 
            site_flags=self.site_flags.to_bits(),
            majorb=self.major_allele,
            majorf=self.major_freq,
            minorb={if let Some(b) = self.minor_allele { b } else { '-' }},
            minorf=self.minor_freq,
            fbc=self.full_base_count,
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
        let mut collated_data: BTreeMap<usize, (SiteInfo, SitePileup, PileupStats)> = BTreeMap::new();
        for (i, info, pileup, stats) in recv_site.iter() {
            collated_data.insert(i, (info, pileup, stats));
        }
        for (i, (info, pileup, stats)) in collated_data.iter() {
            println!("{}\t{}\t{}", info, stats, pileup);
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