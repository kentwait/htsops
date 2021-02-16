use std::path::Path;
use std::io::prelude::*;
use std::collections::{HashMap, HashSet, BTreeSet};

use bio::io::bed;

fn main() {
    let exome_bed_path = Path::new("/Volumes/SSD16TB/reference_seq/hg19.ncbirefseqgenes.knowngene.exons.bed");

    // process bed record to make a nonredundant list
    let mut bed = bed::Reader::from_file(exome_bed_path).unwrap();
    let mut target_locs: HashMap<String, BTreeSet<u64>> = HashMap::new();

    let chromosomes: HashSet<String> = vec![
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
        "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
        "chr21", "chr22", "chrX", "chrY"]
            .into_iter()
            .map(|x| x.to_owned())
            .collect();

    // make expanded nonredundant set
    eprintln!("Expanding BED intervals to make a non-redundant list");
    let mut c = 0;
    let mut j = 0;
    for record in bed.records() {
        let rec = record.ok().expect("Error reading bed record.");
        if !chromosomes.contains(rec.chrom()) {
            continue;
        }
        let interval = target_locs.entry(rec.chrom().to_owned()).or_insert(BTreeSet::new());
        (rec.start()..(rec.end()+1)).for_each(|i| {
            interval.insert(i);
        });
        j += 1;
        c += 1;
        if j >= 10000 {
            eprintln!("{}", c);
            j = 0;
        }
    }
    eprintln!("DONE.");
    eprintln!("Processed {} records.\n", c);

    eprintln!("Creating non-redundant chromosome-interval list");
    let mut c = 0;
    let mut j = 0;
    let mut nr_target_intervals: HashMap<String, BTreeSet<(u64, u64)>> = HashMap::new();
    let stdout = std::io::stdout();
    let mut handle = stdout.lock();
    // loop through btreeset for each chromosome and make collapsed intervals
    for chrom in vec![
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
        "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
        "chr21", "chr22", "chrX", "chrY"] {
        if let Some(locs) = target_locs.get(chrom) {
            let chrom_intervals = nr_target_intervals.entry(chrom.to_owned()).or_insert(BTreeSet::new());
            eprint!("{}\t\t", chrom.clone());
            let mut locs_iter = locs.iter();
            let mut start = *locs_iter.next().unwrap();
            let mut last = *locs_iter.next().unwrap();
            for i in locs_iter {
                let i = *i;
                if i != last+1 {
                    chrom_intervals.insert((start, last));
                    handle.write_fmt(format_args!("{}\t{}\t{}\n", chrom, start, last)).unwrap();
                    start = i;
                    j += 1;
                    c += 1;
                    if j >= 10000 {
                        eprint!("{}\n\t\t", c);
                        j = 0;
                    }
                }
                last = i;
            }
            handle.flush().unwrap();
            eprintln!("{}", c);
        }
    }
    eprintln!("DONE.");
    eprintln!("Processed {} records.\n", c);
}