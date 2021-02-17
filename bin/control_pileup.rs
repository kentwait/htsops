use std::path::Path;
use std::fs::File;
use std::io::prelude::*;
use std::io::{BufWriter, Write};
use flate2::GzBuilder;
use flate2::Compression;
use sled;
use std::primitive::u64;

use std::collections::{HashMap, HashSet, BTreeMap, BTreeSet};

use rust_htslib::{bam, bam::Read};
use bio_types::strand::ReqStrand;
use itertools::multizip;
use bio::io::{bed, fasta};

use std::thread;

use bincode::{serialize, deserialize};

use htsops::pileup::PileupRecord;

fn main() {
    let ctrl_bam_path = Path::new("/Volumes/areca42TB/yokoyama_data_tmp_20201021/UPN95_2/UPN95_B.bam");
    let ref_fa_path = Path::new("/Volumes/SSD16TB/reference_seq/hg19.fa");
    let nr_bed_path = Path::new("/Users/kent/Desktop/nr_list.full.bed");
    let output_bin_path = Path::new("/Users/kent/Desktop/UPN95_B.pileup.bc");
    let output_bin_path_gz = Path::new("/Users/kent/Desktop/UPN95_B.pileup.bc.gz");
    let pileup_db_path = Path::new("/Users/kent/Desktop/UPN95_B.pileup.sled");
    let MIN_BQ = 20;
    let MIN_MQ = 20;

    let sled_config = sled::Config::default()
        .path(pileup_db_path)
        .cache_capacity(10_000_000_000)
        .flush_every_ms(Some(1000))
        .create_new(true);

    // process bed record to make a nonredundant list
    let mut bed = bed::Reader::from_file(nr_bed_path).unwrap();
    let fa = fasta::IndexedReader::from_file(&ref_fa_path).unwrap();
    let mut bam = bam::IndexedReader::from_path(ctrl_bam_path).unwrap();

    let mut i = 0;

    // Open bed record one by one
    let mut bed = bed::Reader::from_file(nr_bed_path).unwrap();

    // Store pileups hierarchically
    // key is chromosome, value is a BtreeMap whose key is (start, end) tuple
    // value Vec<PileupRecord>
    // use sled embedded db
    let db: sled::Db = sled::Config::new()
        .path(pileup_db_path)
        .cache_capacity(10_000_000_000)
        .flush_every_ms(Some(1000))
        .create_new(true)
        .open()
        .unwrap();

    // let mut hierarchical_lookup: HashMap<String, BTreeMap<(u64, u64), usize>> = HashMap::new();
    // let mut reverse_hierarchical_lookup: Vec<(String, u64, u64)> = Vec::new();
    // let mut pileups: Vec<Vec<PileupRecord>> = Vec::new();
    
    for (i, record) in bed.records().enumerate() {
        let rec = record.ok().expect("Error reading bed record.");
        let chrom = rec.chrom();
        let start = rec.start();
        let end = rec.end();
        print!("{}:{}-{}\t\tlen/actual\t:\t{}:", chrom, start, end, end-start);

        // fetch interval
        bam.fetch((chrom, start, end)).unwrap();
        let mut pileups: Vec<PileupRecord> = Vec::new();
        for p in bam.pileup() {
            let pileup = p.unwrap(); 
            let pos = pileup.pos() as u64;
            if pos < start { continue }
            if pos >= end { break }

            // Insert pileup
            let pileup = PileupRecord::from_pileup(chrom, pileup).quality_filter(MIN_BQ, MIN_MQ);
            if pileup.total_depth() > 0 {
                pileups.push(pileup);
            }
        }
        println!("{}", pileups.len());
        if pileups.len() > 0 {
            db.insert(format!("{}:{}-{}", chrom, start, end).as_bytes(), serialize(&pileups).unwrap()).unwrap();
        }

        // if i % 999 == 0 {
        //     db.flush().unwrap();
        //     println!("* Flushed to disk *");
        // }
    }
}

// TODO: Make an index for target bams to fetch against


// Make a test that checks whether the base+indel counter and the pileup matches
// for record in bed.records() {
//     let rec = record.ok().expect("Error reading bed record.");
//     let interval_start = rec.start();
//     let interval_end = rec.end();

//     println!("{} {}-{} {:?} {:?} {:?}", rec.chrom(), rec.start(), rec.end(), rec.name(), rec.score(), rec.strand());

//     // fetch interval
//     bam.fetch((rec.chrom(), rec.start(), rec.end())).unwrap();
//     let mut j = 0;
//     let mut prev_pos = interval_start - 1;
//     for p in bam.pileup() {
//         let pileup = p.unwrap();
//         let pos = pileup.pos() as u64;
//         if pos < interval_start || pos > interval_end { continue }
//         if pos != prev_pos + 1 { panic!("current pos:{} previous: {}", pos, prev_pos) }
//         println!("\tC {}:{} depth {}", pileup.tid(), pileup.pos(), pileup.depth());
    
//         // pileup per site in exon as text
//         print!("\t\t");
//         for alignment in pileup.alignments() {
//             if alignment.is_del() || alignment.is_refskip() {
//                 continue;
//             }
//             let base = match (alignment.record().seq()[alignment.qpos().unwrap()], alignment.record().strand()) {
//                 (65, ReqStrand::Forward) => 'A',
//                 (67, ReqStrand::Forward) => 'C',
//                 (71, ReqStrand::Forward) => 'G',
//                 (84, ReqStrand::Forward) => 'T',
//                 (78, ReqStrand::Forward) => 'N',
//                 (65, ReqStrand::Reverse) => 'a',
//                 (67, ReqStrand::Reverse) => 'c',
//                 (71, ReqStrand::Reverse) => 'g',
//                 (84, ReqStrand::Reverse) => 't',
//                 (78, ReqStrand::Reverse) => 'n',
//                 (_, _) => '?',
//             };
//             print!("{}", base);
//         }
//         print!("\n");

//         // pileup per site in exon as map
//         let mut forward_bases: HashMap<char, i64> = [
//             ('a', 0), ('c', 0), ('g', 0), ('t', 0), ('n', 0)]
//             .iter()
//             .cloned()
//             .collect();
//         let mut reverse_bases: HashMap<char, i64> = [
//             ('a', 0), ('c', 0), ('g', 0), ('t', 0), ('n', 0)]
//             .iter()
//             .cloned()
//             .collect();
        
//         let mut forward_depth = 0;
//         let mut reverse_depth = 0;
//         let mut insertion = 0;
//         let mut deletion = 0;
//         for alignment in pileup.alignments() {
//             // mark indel start
//             if alignment.is_refskip() {
//                 insertion += 1;
//                 continue;
//             }
//             if alignment.is_del() {
//                 deletion += 1;
//                 continue;
//             } 
//             //     match alignment.indel() {
//             //         bam::pileup::Indel::Ins(len) => println!("\t\tInsertion of length {} between this and next position.", len),
//             //         bam::pileup::Indel::Del(len) => println!("\t\tDeletion of length {} between this and next position.", len),
//             //         bam::pileup::Indel::None => ()
//             //     }
            
//             let record = alignment.record();
//             match alignment.record().strand() {
//                 ReqStrand::Forward => {
//                     forward_depth += 1;
//                     match record.seq()[alignment.qpos().unwrap()] {
//                         65 => *forward_bases.get_mut(&'a').unwrap() += 1,
//                         67 => *forward_bases.get_mut(&'c').unwrap() += 1,
//                         71 => *forward_bases.get_mut(&'g').unwrap() += 1,
//                         84 => *forward_bases.get_mut(&'t').unwrap() += 1,
//                         78 => *forward_bases.get_mut(&'n').unwrap() += 1,
//                         _ => (),
//                     }
//                 },
//                 ReqStrand::Reverse => {
//                     reverse_depth += 1;
//                     match record.seq()[alignment.qpos().unwrap()] {
//                         65 => *reverse_bases.get_mut(&'a').unwrap() += 1,
//                         67 => *reverse_bases.get_mut(&'c').unwrap() += 1,
//                         71 => *reverse_bases.get_mut(&'g').unwrap() += 1,
//                         84 => *reverse_bases.get_mut(&'t').unwrap() += 1,
//                         78 => *reverse_bases.get_mut(&'n').unwrap() += 1,
//                         _ => (),
//                     }
//                 },
//             };  
//         }
//         println!("\t\t{}:{} A:{} C:{} G:{} T:{} N:{} a:{} c:{} g:{} t:{} n:{}\tins:{} del:{}",
//             forward_depth, reverse_depth,
//             forward_bases.get(&'a').unwrap(),
//             forward_bases.get(&'c').unwrap(),
//             forward_bases.get(&'g').unwrap(),
//             forward_bases.get(&'t').unwrap(),
//             forward_bases.get(&'n').unwrap(),
//             reverse_bases.get(&'a').unwrap(),
//             reverse_bases.get(&'c').unwrap(),
//             reverse_bases.get(&'g').unwrap(),
//             reverse_bases.get(&'t').unwrap(),
//             reverse_bases.get(&'n').unwrap(),
//             insertion, deletion
//         );

//         if j > 10 { break; }
//         j += 1;

//         prev_pos = pos;
//     }

//     if i > 100 { break; }
//     i += 1;
// }

// Example code for pileup
// let mut i = 0;
// for p in bam.pileup() {
//     let pileup = p.unwrap();
//     println!("C {}:{} depth {}", pileup.tid(), pileup.pos(), pileup.depth());

//     // for alignment in pileup.alignments() {
//     //     if !alignment.is_del() && !alignment.is_refskip() {
//     //         println!("Base {}", alignment.record().seq()[alignment.qpos().unwrap()]);
//     //     }
//     //     // mark indel start
//     //     match alignment.indel() {
//     //         bam::pileup::Indel::Ins(len) => println!("Insertion of length {} between this and next position.", len),
//     //         bam::pileup::Indel::Del(len) => println!("Deletion of length {} between this and next position.", len),
//     //         bam::pileup::Indel::None => ()
//     //     }
//     // }
//     if i > 10000 { break; }
//     i += 1;
// }