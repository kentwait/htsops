use std::path::Path;
use std::fs::File;
use std::io::{self, BufRead};
use std::collections::HashMap;

use rust_htslib::tbx::{self, Read as TbxRead};

use crate::constant::*;


pub fn validate_path(v: String) -> Result<(), String> {
    let path = Path::new(&v);
    match path.is_file() {
        true => Ok(()),
        false => Err(String::from("Input file does not exist")),
    }
}

pub fn validate_qual(v: String) -> Result<(), String> {
    match v.parse::<u8>() {
        Ok(_) => Ok(()),
        Err(_) => Err(String::from("Value is not a positive integer")),
    }
}

fn flatten_vec_vec<T>(nested: Vec<Vec<T>>) -> Vec<T> {
    nested.into_iter().flatten().collect()
}

pub fn filelist_to_vec(path_str: &str) -> Vec<String> {
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

pub fn samplelist_to_vec(path_str: &str) -> Vec<String> {
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
    sample_names
}

pub fn read_tabix(path: &str) -> (tbx::Reader, HashMap<&str, u64>) {
    let tbx_reader = tbx::Reader::from_path(&path)
        .expect(&format!("Tabix reader could not open {}", &path));
    let chrom_tid_lookup: HashMap<&str, u64> = HG19_CHROMS.iter().map(|c| {
        match tbx_reader.tid(c) {
            Ok(tid) => (c.to_owned(), tid),
            Err(_) => panic!("Could not resolve contig ID"),
        }
    }).collect();
    (tbx_reader, chrom_tid_lookup)
}