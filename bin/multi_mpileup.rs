use std::{thread, time};
use std::collections::VecDeque;

use std::process::{Command, Stdio};
use std::fs::File;
use std::os::unix::io::{FromRawFd, IntoRawFd};
use std::io::{self, BufRead};
use std::path::Path;

use itertools;
use clap::{Arg, App};

fn get_paths_for_filelist(path_str: &str) -> Vec<String> {
    let file = File::open(path_str).unwrap();
    let paths: Vec<String> = io::BufReader::new(file).lines()
        .filter_map(|l| {
            if let Ok(line) = l {
                // Skip empty lines
                if line.len() > 0 { Some(line) } else { None }
            } 
            else { None }
        })
        .collect();
    paths
}

fn validate_path(v: String) -> Result<(), String> {
    let path = Path::new(&v);
    match path.is_file() {
        true => Ok(()),
        false => Err(String::from("Input file does not exist")),
    }
}

fn main() {
    // Argument parser to get file list paths
    let matches = App::new("htsops multi_mpileup.rs")
        .version("0.0.1")
        .author("Kent Kawashima <kentkawashima@gmail.com>")
        .about("Runs mpileup for each bam file")
        .arg(Arg::with_name("threads")
            .short("t")
            .long("threads")
            .default_value("4")
            .value_name("INT")
            .help("Number of threads to use")
            .takes_value(true))
        .arg(Arg::with_name("fasta_path")
            .short("f")
            .long("fasta_path")
            .value_name("FILE")
            .help("Reference FASTA file")
            .takes_value(true))
        .arg(Arg::with_name("positions_path")
            .short("p")
            .long("positions_path")
            .value_name("FILE")
            .help("Genomic positions of interest as a BED file")
            .takes_value(true))
        .arg(Arg::with_name("BAM_LIST")
            .help("Contains list of bam files to process")
            .required(true)
            .index(1)
            .validator(validate_path))
        .arg(Arg::with_name("OUTPUT_LIST")
            .help("Contains corresponding list of output paths to save to")
            .required(true)
            .index(2)
            .validator(validate_path))
        .get_matches();

    let bam_paths_file = matches.value_of("BAM_LIST").unwrap();
    let bgz_paths_file = matches.value_of("OUTPUT_LIST").unwrap();
    let ref_fasta_path = matches.value_of("fasta_path").unwrap();
    let positions_path = matches.value_of("positions_path").unwrap();
    let threads = matches.value_of("threads").unwrap_or("1").parse::<usize>().unwrap();

    // List of bam paths
    // List of corresponding bgzip paths
    let bam_paths: Vec<String> = get_paths_for_filelist(bam_paths_file);
    let bgz_paths: Vec<String> = get_paths_for_filelist(bgz_paths_file);
    if bam_paths.len() != bgz_paths.len() {
        panic!("number of BAM files != number of output paths")
    }

    // Create a fifo queue
    let mut bam_bgz_list: VecDeque<(String, String)>= itertools::zip_eq(bam_paths, bgz_paths).collect();
    // Initialize deques to wait on spawned processes
    let mut bam_bgz_queue: VecDeque<(String, String)>= VecDeque::with_capacity(threads);
    let mut bgz_queue: VecDeque<std::process::Child> = VecDeque::with_capacity(threads);

    // Saturate processing from queue depending on the number of threads
    for _ in 0..threads {
        if let Some((bam_path, bgz_path)) = bam_bgz_list.pop_front() {
            let mut mpileup_process = Command::new("samtools")
                .arg("mpileup")
                .args(&["--min-BQ", "0"])
                .args(&["--min-MQ", "0"])
                .args(&["--adjust-MQ", "50"])
                .arg("--ignore-overlaps")
                .arg("--output-MQ")
                .args(&["--fasta-ref", ref_fasta_path])
                .args(&["--positions", positions_path])
                .arg(&bam_path)
                .stdout(Stdio::piped())
                .stderr(Stdio::null())
                .spawn()
                .expect(&format!("failed to execute samtools mpileup for {}", &bam_path));
            let fd = File::create(&bgz_path).unwrap().into_raw_fd();
            let file_out = unsafe {Stdio::from_raw_fd(fd)};
            let bgzip_process = Command::new("bgzip")
                .stdin(mpileup_process.stdout.take().unwrap())
                .stdout(file_out)
                .stderr(Stdio::null())
                .spawn()
                .expect(&format!("failed to execute bgzip for {}", bam_path));
            bam_bgz_queue.push_back((bam_path.to_owned(), bgz_path.to_owned()));
            bgz_queue.push_back(bgzip_process);
            println!("queue  : {}", bam_path);
        }
    }
    // Saturate the queue by replacing finished processes until no more bam files
    while !bam_bgz_list.is_empty() {
        if let Some(mut bgzip_proc) = bgz_queue.pop_front() {
            let (bam_path, bgz_path) = bam_bgz_queue.pop_front().unwrap();

            // check bgzip if done
            match bgzip_proc.try_wait() {
                Ok(Some(_)) => println!("mpileup: {}", bam_path),
                Ok(None) => {
                    thread::sleep(time::Duration::from_millis(1000));
                    // Put back in queue
                    bgz_queue.push_back(bgzip_proc);
                    bam_bgz_queue.push_back((bam_path, bgz_path));
                    continue;
                },
                Err(e) => println!("ERROR bgzip  : {} {}", &bam_path, e),
            }

            // Perform indexing after mpileup and bgzip are finished
            let status = Command::new("tabix")
                    .args(&["-s", "1"])
                    .args(&["-b", "2"])
                    .args(&["-e", "2"])
                    .arg(&bgz_path)
                    .output()
                    .expect(&format!("error executing tabix for {}", &bam_path))
                    .status;
            if status.success() {
                println!("tabix  : {} {}", &bam_path, status.code().take().unwrap());
            } else {
                match status.code() {
                    Some(code) => println!("ERROR tabix  : {} {}", &bam_path, code),
                    None       => println!("tabix terminated by signal")
                }
            }

            // Add a new process to the queue
            if let Some((bam_path, bgz_path)) = bam_bgz_list.pop_front() {
                // TODO: make this into a function because it's repeated twice
                let mut mpileup_process = Command::new("samtools")
                    .arg("mpileup")
                    .args(&["--min-BQ", "0"])
                    .args(&["--min-MQ", "0"])
                    .args(&["--adjust-MQ", "50"])
                    .arg("--ignore-overlaps")
                    .arg("--output-MQ")
                    .args(&["--fasta-ref", ref_fasta_path])
                    .args(&["--positions", positions_path])
                    .arg(&bam_path)
                    .stdout(Stdio::piped())
                    .stderr(Stdio::null())
                    .spawn()
                    .expect(&format!("failed to execute samtools mpileup for {}", &bam_path));
                let fd = File::create(&bgz_path).unwrap().into_raw_fd();
                let file_out = unsafe {Stdio::from_raw_fd(fd)};
                let bgzip_process = Command::new("bgzip")
                    .stdin(mpileup_process.stdout.take().unwrap())
                    .stdout(file_out)
                    .stderr(Stdio::null())
                    .spawn()
                    .expect(&format!("failed to execute bgzip for {}", bam_path));
                bam_bgz_queue.push_back((bam_path.to_owned(), bgz_path.to_owned()));
                bgz_queue.push_back(bgzip_process);
                println!("queue  : {}", bam_path);
            }
        }
    }
    println!("FINISHED mpileup | bgzip && tabix");
}