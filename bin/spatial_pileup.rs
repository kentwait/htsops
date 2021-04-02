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
use std::collections::VecDeque;
use std::thread;
use std::sync::mpsc;

use clap::{Arg, App};

use htsops::util::{validate_path, validate_qual, filelist_to_vec};


fn main() {
    let matches = App::new("htsops spatial_pileup.rs")
        .version("0.0.1")
        .author("Kent Kawashima <kentkawashima@gmail.com>")
        .about("Processes multiple pileups into one")
        .arg(Arg::with_name("threads")
            .short("t")
            .long("threads")
            .default_value("4")
            .value_name("INT")
            .help("Number of threads to use")
            .takes_value(true))
        .arg(Arg::with_name("min_control_cov")
            .long("min-control-cov")
            .default_value("20")
            .value_name("INT")
            .help("Minimum coverage depth for the control sample. Sites where the coverage depth is less than this value will not be used.")
            .takes_value(true)
            .validator(validate_qual))
        .arg(Arg::with_name("min_control_bq")
            .long("min-control-bq")
            .default_value("20")
            .value_name("INT")
            .help("Minimum base quality for the control sample. Bases whose base quality is less than this value will not be counted.")
            .takes_value(true)
            .validator(validate_qual))
        .arg(Arg::with_name("min_control_mq")
            .long("min-control-mq")
            .default_value("40")
            .value_name("INT")
            .help("Minimum mapping quality for the control sample. Bases belonging to reads with mapping quality less than this value will not be counted.")
            .takes_value(true)
            .validator(validate_qual))
        .arg(Arg::with_name("min_target_cov")
            .long("min-target-cov")
            .default_value("20")
            .value_name("INT")
            .help("Minimum coverage depth for target samples. Sites where the coverage depth is less than this value will not be used.")
            .takes_value(true)
            .validator(validate_qual))
        .arg(Arg::with_name("min_target_bq")
            .long("min-target-bq")
            .default_value("20")
            .value_name("INT")
            .help("Minimum base quality for target samples. Bases whose base quality is less than this value will not be counted.")
            .takes_value(true)
            .validator(validate_qual))
        .arg(Arg::with_name("min_target_mq")
            .long("min-target-mq")
            .default_value("40")
            .value_name("INT")
            .help("Minimum mapping quality for target samples. Bases belonging to reads with mapping quality less than this value will not be counted.")
            .takes_value(true)
            .validator(validate_qual))
        .arg(Arg::with_name("PILEUP_LIST")
            .help("Contains list of pileup files to process")
            // .required(true)
            .index(1)
            .validator(validate_path))
        .get_matches();
    let min_control_cov: usize = matches.value_of("min_control_cov").unwrap().parse::<usize>().unwrap();
    let min_control_bq: u8 = matches.value_of("min_control_bq").unwrap().parse::<u8>().unwrap();
    let min_control_mq: u8 = matches.value_of("min_control_mq").unwrap().parse::<u8>().unwrap();
    let threads: usize = matches.value_of("threads").unwrap().parse::<usize>().unwrap();
    // let pileup_path_file = matches.value_of("PILEUP_LIST").unwrap();
    // let pileup_path_list: VecDeque<String> = filelist_to_vec(pileup_path_file).into_iter().collect();

    // Read all the mpileups and get the sites which is found in all
    // Saturate processing from queue based on the number of threads
    let (tx, rx) = mpsc::channel();

    for i in 0..threads {
        let ttx = tx.clone();
        thread::spawn(move || {
            ttx.send(i).unwrap();
        });
    }
    // Close send to finish
    drop(tx);

    for received in rx {
        println!("Got: {}", received);
    }

    // Visit each site
    // Evaluate control first. Do bq and mq filtering and test coverage
    // Evaluate samples. Do bq and mq for each sample and test coverage for each as well
    // If passed control and passed ALL samples, then include in the output
}