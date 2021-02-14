use rust_htslib::{bam, bam::Read};
use itertools::multizip;
use bio::io::fasta;
use std::thread;

fn main() {
    let tgt_bam_path = "/Volumes/areca42TB/yokoyama_data/hg19/UPN95/UPN95_2_1.rg.sort.bam";
    let ref_path = "/Volumes/SSD16TB/reference_seq/hg19.fa";

    let ref_record = fasta::Reader::from_file(ref_path);
    let mut tgt_bam = bam::Reader::from_path(tgt_bam_path).unwrap();
    let mut i = 0;

    let handle = thread::spawn(|| {
        let mut i = 0;
        let ctrl_bam_path = "/Volumes/areca42TB/yokoyama_data_tmp_20201021/UPN95_2/UPN95_B.bam";
        let mut bam = bam::Reader::from_path(ctrl_bam_path).unwrap();
        for p in bam.pileup() {
            let pileup = p.unwrap();
            println!("C {}:{} depth {}", pileup.tid(), pileup.pos(), pileup.depth());
        
            // for alignment in pileup.alignments() {
            //     if !alignment.is_del() && !alignment.is_refskip() {
            //         println!("Base {}", alignment.record().seq()[alignment.qpos().unwrap()]);
            //     }
            //     // mark indel start
            //     match alignment.indel() {
            //         bam::pileup::Indel::Ins(len) => println!("Insertion of length {} between this and next position.", len),
            //         bam::pileup::Indel::Del(len) => println!("Deletion of length {} between this and next position.", len),
            //         bam::pileup::Indel::None => ()
            //     }
            // }
            if i > 10000 { break; }
            i += 1;
        }
    });

    // pileup over all covered sites
    for p in tgt_bam.pileup() {
        let pileup = p.unwrap();
        println!("T {}:{} depth {}", pileup.tid(), pileup.pos(), pileup.depth());
    
        // for alignment in pileup.alignments() {
        //     if !alignment.is_del() && !alignment.is_refskip() {
        //         println!("Base {}", alignment.record().seq()[alignment.qpos().unwrap()]);
        //     }
        //     // mark indel start
        //     match alignment.indel() {
        //         bam::pileup::Indel::Ins(len) => println!("Insertion of length {} between this and next position.", len),
        //         bam::pileup::Indel::Del(len) => println!("Deletion of length {} between this and next position.", len),
        //         bam::pileup::Indel::None => ()
        //     }
        // }
        if i > 10000 { break; }
        i += 1;
    }

    handle.join().unwrap();
}
