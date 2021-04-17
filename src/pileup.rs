use std::collections::HashMap;
use serde::{Serialize, Deserialize};

use std::fmt;
use crate::{Base, BaseCount, AlleleSet};


#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SitePileup {
    sample_name: String,
    base_chars: Vec<char>,
    // pub bases: Vec<Base>,
    indels: HashMap<String, usize>,
    bqs: Vec<u8>, 
    mqs: Vec<u8>,
}
impl SitePileup {
    // constructor
    pub fn from_str(sample_name: &str, ref_char: &char, cov: usize, base_str: &str, bq_str: &str, mq_str: &str) -> SitePileup {
        // Immediately return if base_str is "*" which means empty
        let sample_name = sample_name.to_owned();
        if base_str == "*" {
            let base_chars = Vec::new();
            let indels = HashMap::new();
            let bqs = Vec::new();
            let mqs = Vec::new();
            return SitePileup{ sample_name, base_chars, indels, bqs, mqs }
        }
        let (mut base_chars, indels) = SitePileup::parse_base_str(base_str, ref_char);
        let mut bqs = SitePileup::parse_qual_str(bq_str);
        let mut mqs = SitePileup::parse_qual_str(mq_str);

        if base_chars.len() != cov {
            panic!(format!("Base length [{}] and coverage [{}] do not match:\n{}", 
            base_chars.len(), cov, base_str));
        }
        else if base_chars.len() != bqs.len() {
            panic!(format!("Base length [{}] and base qualities length [{}] do not match:\n{}\n{}", 
            base_chars.len(), bqs.len(), base_str, bq_str));
        } 
        else if base_chars.len() != mqs.len() {
            panic!(format!("Base length [{}] and mapping qualities length [{}] do not match:\n{}\n{}", 
            base_chars.len(), mqs.len(), base_str, mq_str));
        }
        let keep: Vec<bool> = base_chars.iter().map(|b| {
                if b == &'D' || b == &'d' { false }
                else { true }
            }).collect();
        {
            let mut i = 0;
            base_chars.retain(|_| (keep[i], i += 1).0);
        }
        {
            let mut i = 0;
            bqs.retain(|_| (keep[i], i += 1).0);
        }
        {
            let mut i = 0;
            mqs.retain(|_| (keep[i], i += 1).0);
        }
        SitePileup{ sample_name, base_chars, indels, bqs, mqs }
    }

    // getters
    pub fn sample_name(&self) -> &str { &self.sample_name }
    pub fn base_chars(&self) -> &Vec<char> { &self.base_chars }
    pub fn indels(&self) -> &HashMap<String, usize> { &self.indels }
    pub fn bqs(&self) -> &Vec<u8>  { &self.bqs }
    pub fn mqs(&self) -> &Vec<u8>  { &self.mqs }

    // TODO: transform this into an iterator
    fn parse_qual_str(q_str: &str) -> Vec<u8> {
        q_str.chars().map(|c| (c as u8) - 33).collect()
    }
    
    // TODO: transform this into an iterator that spits out an enum indicating
    // base, indel, del,
    fn parse_base_str(base_str: &str, ref_char: &char) -> (Vec<char>, HashMap<String, usize>) {
        let mut bases: Vec<char> = Vec::new();
        let mut indels: HashMap<String, usize> = HashMap::new();
        
        // temp vars
        let mut i = 0;
        let base_chars: Vec<char> = base_str.chars().map(|c| c.to_owned()).collect();
        let mut indel_str: Vec<char> = Vec::new();
        
        while i < base_chars.len() {
            let c = base_chars[i];
            match c {
                // indel
                ch if (ch == '+') || (ch == '-') => {
                    indel_str.push(ch);
                    i += 1;
                    // search all digits
                    let mut indel_len_vec: Vec<char> = Vec::new();
                    for d in base_chars[i..].iter() {
                        match d.is_ascii_digit() {
                            true => indel_len_vec.push(d.to_owned()),
                            false => break,
                        }
                    }
                    if indel_len_vec.len() == 0 {
                        panic!(format!("Invalid indel_len_vec [{:?}] from base str: {:?}", 
                            indel_len_vec, base_str))
                    }
                    let indel_len: usize = indel_len_vec.iter().cloned().collect::<String>().parse::<usize>().unwrap();
                    base_chars[i..(i+indel_len_vec.len()+indel_len)].iter().for_each(|c| {
                        indel_str.push(c.to_owned());
                    });
                    i += indel_len_vec.len() + indel_len;
                    let indel_string: String = indel_str.iter().cloned().collect();
                    *indels.entry(indel_string).or_insert(0) += 1;
                    indel_str = Vec::new();
                },
                // not indel
                // same as ref_char
                '.' => {
                    bases.push(ref_char.to_ascii_uppercase());
                    i += 1;
                },
                ',' => {
                    bases.push(ref_char.to_ascii_lowercase());
                    i += 1;
                },
                // substitution
                'A' => {
                    bases.push(c);
                    i += 1;
                },
                'C' => {
                    bases.push(c);
                    i += 1;
                },
                'G' => {
                    bases.push(c);
                    i += 1;
                },
                'T' => {
                    bases.push(c);
                    i += 1;
                },
                'N' => {
                    bases.push(c);
                    i += 1;
                },
                'a' => {
                    bases.push(c);
                    i += 1;
                },
                'c' => {
                    bases.push(c);
                    i += 1;
                },
                'g' => {
                    bases.push(c);
                    i += 1;
                },
                't' => {
                    bases.push(c);
                    i += 1;
                },
                'n' => {
                    bases.push(c);
                    i += 1;
                },
                // * or # ref base deletion, CIGAR “D”
                '*' => {
                    bases.push('D');
                    i += 1;
                },
                '#' => {
                    bases.push('d');
                    i += 1;
                },
                // start of read
                '^' => i += 2,
                // others
                // > or < ref skip, CIGAR “N”
                // $ indicating end of read
                _ => i += 1,
            }
        }
        
        (bases, indels)
    }

    // drop N's and deletions
    pub fn cleanup(&mut self, drop_n: bool, drop_del: bool) -> Option<usize> {
        let keep: Vec<bool> = self.base_chars.iter().enumerate()
            .map(|(i, b)| {
                if drop_n == true {
                    if *b == 'N' || *b == 'n' { return false }
                }
                if drop_del == true {
                    if *b == 'D' || *b == 'd' { return false }
                }
                return true
            }).collect();
        {
            let mut i = 0;
            self.base_chars.retain(|_| (keep[i], i += 1).0);
        }
        {
            let mut i = 0;
            self.bqs.retain(|_| (keep[i], i += 1).0);
        }
        {
            let mut i = 0;
            self.mqs.retain(|_| (keep[i], i += 1).0);
        }
        if self.base_chars.len() > 0 {
            return Some(self.base_chars.len())
        }
        return None
    }

    // quality filter and drop bases inplace
    // TODO: Fix algo
    pub fn quality_filter_inplace(&mut self, min_bq: u8, min_mq: u8) -> Option<usize> {
        let keep: Vec<bool> = self.base_chars.iter().enumerate()
            .map(|(i, _)| {
                let bq = self.bqs[i];
                let mq = self.mqs[i];
                if bq >= min_bq && mq >= min_mq { true }
                else { false }
            }).collect();
        {
            let mut i = 0;
            self.base_chars.retain(|_| (keep[i], i += 1).0);
        }
        {
            let mut i = 0;
            self.bqs.retain(|_| (keep[i], i += 1).0);
        }
        {
            let mut i = 0;
            self.mqs.retain(|_| (keep[i], i += 1).0);
        }
        if self.base_chars.len() > 0 {
            return Some(self.base_chars.len())
        }
        return None
    }

    // quality filter and makes a new SitePileup of the result
    pub fn quality_filter(&self, min_bq: u8, min_mq: u8) -> SitePileup {
        let passed_pos: Vec<usize> = self.base_chars.iter().enumerate()
            .filter_map(|(i, &b)| {
                if (b == 'N') || (b == 'n') { return None }
                if (b == 'D') || (b == 'd') { return None }
                let bq = self.bqs[i];
                let mq = self.mqs[i];
                if (bq >= min_bq) && (mq >= min_mq) { return Some(i) }
                return None
            }).collect();
        SitePileup {
            sample_name: self.sample_name.to_owned(),
            base_chars: passed_pos.iter().map(|&i| self.base_chars[i]).collect(),
            indels: HashMap::new(),
            bqs: passed_pos.iter().map(|&i| self.bqs[i]).collect(),
            mqs: passed_pos.iter().map(|&i| self.mqs[i]).collect(),
        }
    }

    pub fn to_base_count(&self) -> BaseCount {
        BaseCount::from_char_vec(&self.base_chars)
    }
    
    pub fn cov(&self) -> usize { self.base_chars.len() }
}
impl fmt::Display for SitePileup {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let bases_str: String = self.base_chars.iter()
            .cloned()
            .collect();
        let bqs_str: String = self.bqs.iter()
            .map(|&u| (u + 33) as char)
            .collect::<Vec<char>>()
            .into_iter()
            .collect();
        let mqs_str: String = self.mqs.iter()
            .map(|&u| (u + 33) as char)
            .collect::<Vec<char>>()
            .into_iter()
            .collect();
        write!(f, "{}\t{}\t{}", bases_str, bqs_str, mqs_str)
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SpatialSitePileup {
    // Chromosome name such as "chr1" in hg19
    pub chrom: String, 
    // Chromosome position (1-based?)
    pub pos: u64, 
    // Reference character in hg19
    pub ref_char: char,
    // Control SitePileup 
    pub control_pileup: SitePileup,
    // Vector of target SitePileups
    pub pileups: Vec<SitePileup>,
}
impl SpatialSitePileup {
    // TODO: Change sample_list to hasmap of sample_list and 2d coords
    pub fn parse_mpileup_row(row: &str, sample_list: &Vec<&str>, control_name: &str) -> SpatialSitePileup {
        let cols: Vec<String> = row.split('\t').map(|s| s.to_owned()).collect();
        let chrom: String = cols[0].to_owned();
        let pos: u64 = cols[1].parse::<u64>().unwrap();
        let ref_char: char = cols[2].chars().next().unwrap();
        
        // ENHANCEMENT: Spawn multiple threads and let each thread work with one group of cols
        let mut control_idx = 0;
        let mut pileups: Vec<SitePileup> = sample_list.iter().enumerate()
            .map(|(i, &sample_name)| {
                if sample_name == control_name {
                    control_idx = i;
                }
                SitePileup::from_str(
                    sample_name,
                    &ref_char, (&cols[3+(i*4)]).parse::<usize>().unwrap(),
                    &cols[4+(i*4)], &cols[5+(i*4)], &cols[6+(i*4)])
            }).collect();
        let control_pileup = pileups.remove(control_idx);
        SpatialSitePileup{ chrom, pos, ref_char, control_pileup, pileups }
    }
}

#[derive(Debug, Clone)]
pub struct SiteInfo {
    chrom: String,
    pos: u64,
    ref_base: Base,
}
impl SiteInfo {
    pub fn new(chrom: &str, pos: u64, ref_char: &char) -> Self {
        Self {
            chrom: chrom.to_owned(),
            pos: pos,
            ref_base: Base::from_char(ref_char).unwrap(),
        }
    }
    pub fn from_str(s: &str, sep: &str) -> Self {
        let data_vec: Vec<&str> = s.split(sep).collect();
        Self {
            chrom: data_vec[0].to_owned(),
            pos: data_vec[1].parse::<u64>().unwrap(),
            ref_base: Base::from_char(&data_vec[2].chars().next().unwrap()).unwrap(),
        }
    }
    pub fn chrom(&self) -> &str {
        self.chrom.as_ref()
    }
    pub fn pos(&self) -> u64 {
        self.pos
    }
    pub fn ref_base(&self) -> Base {
        self.ref_base
    }
}
impl fmt::Display for SiteInfo {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}\t{}\t{}", self.chrom, self.pos, self.ref_base)
    }
}

#[derive(Debug, Clone)]
pub struct PileupStats {
    allele_set: AlleleSet,
    cov: usize,
    allele_orientation_table: [u32; 4],
    allele_orientation_fisher_pval: Option<f64>,
    allele_orientation_odds_ratio: Option<f64>,
}
impl PileupStats {
    // constructors
    // TODO: force a particular base to be the major allele, minor_allele
    pub fn from_pileup(pileup: &SitePileup) -> Option<PileupStats> {
        let base_count = pileup.to_base_count();
        let cov = base_count.total();
        let allele_set = match AlleleSet::from_base_count(base_count) {
            Some(allele_set) => allele_set,
            None => return None
        };
        let major_allele = allele_set.major_allele();
        
        let a: u32 = allele_set.base_count().forward_count_of(major_allele) as u32;
        let b: u32 = allele_set.base_count().reverse_count_of(major_allele) as u32;
        let (c, d) = match allele_set.minor_allele() {
            Some(minor_allele) => (
                allele_set.base_count().forward_count_of(minor_allele) as u32,
                allele_set.base_count().forward_count_of(minor_allele) as u32,
            ),
            None => (0, 0),
        };
        let allele_orientation_table = [
            a, b,
            c, d,
        ];
        Some(PileupStats {
            allele_set,
            cov,
            allele_orientation_table,
            allele_orientation_fisher_pval: None,
            allele_orientation_odds_ratio: None,
        })
    }
    pub fn from_str(s: &str, sep: &str) -> Option<PileupStats> {
        // format
        // {major_base}:{major_count}:{major_freq:.3}\t{minor_base}:{minor_count}:{minor_freq:.3}\t{base_count}\t{cov}
        let data_vec: Vec<&str> = s.split(sep).collect();
        let base_count = BaseCount::from_str(data_vec[data_vec.len()-2], ":");
        let cov: usize = data_vec[data_vec.len()-1].parse::<usize>().unwrap();
        let allele_set = match AlleleSet::from_base_count(base_count) {
            Some(allele_set) => allele_set,
            None => return None
        };
        let major_allele = allele_set.major_allele();
        
        let a: u32 = allele_set.base_count().forward_count_of(major_allele) as u32;
        let b: u32 = allele_set.base_count().reverse_count_of(major_allele) as u32;
        let (c, d) = match allele_set.minor_allele() {
            Some(minor_allele) => (
                allele_set.base_count().forward_count_of(minor_allele) as u32,
                allele_set.base_count().forward_count_of(minor_allele) as u32,
            ),
            None => (0, 0),
        };
        let allele_orientation_table = [
            a, b,
            c, d,
        ];
        Some(PileupStats {
            allele_set,
            cov,
            allele_orientation_table,
            allele_orientation_fisher_pval: None,
            allele_orientation_odds_ratio: None,
        })
    }
    pub fn from_allele_set(allele_set: AlleleSet) -> Self {
        let major_allele = allele_set.major_allele();
        let a: u32 = allele_set.base_count().forward_count_of(major_allele) as u32;
        let b: u32 = allele_set.base_count().reverse_count_of(major_allele) as u32;
        let (c, d) = match allele_set.minor_allele() {
            Some(minor_allele) => (
                allele_set.base_count().forward_count_of(minor_allele) as u32,
                allele_set.base_count().forward_count_of(minor_allele) as u32,
            ),
            None => (0, 0),
        };
        let allele_orientation_table = [
            a, b,
            c, d,
        ];
        Self {
            cov: allele_set.num_bases(),
            allele_set,
            allele_orientation_table,
            allele_orientation_fisher_pval: None,
            allele_orientation_odds_ratio: None,
        }
    }

    // getters
    pub fn major_allele(&self) -> Base { self.allele_set.major_allele() }
    pub fn minor_allele(&self) -> Option<Base> { self.allele_set.minor_allele() }
    pub fn other_alleles(&self) -> Vec<Base> { self.allele_set.other_alleles() }

    pub fn major_count(&self) -> usize { self.allele_set.major_count() }
    pub fn minor_count(&self) -> usize { self.allele_set.minor_count() }
    pub fn others_count(&self) -> usize { self.allele_set.others_count() }

    pub fn major_freq(&self) -> f64 { self.allele_set.major_freq() }
    pub fn minor_freq(&self) -> f64 { self.allele_set.minor_freq() }

    pub fn allele_set(&self) -> &AlleleSet { &self.allele_set }
    pub fn base_count(&self) -> &BaseCount { self.allele_set.base_count() }

    // stats
    pub fn num_alleles(&self) -> usize { self.allele_set.len() }
    pub fn cov(&self) -> usize { self.cov }
    pub fn allele_orientation_fisher_test(&self) -> Option<f64> { self.allele_orientation_fisher_pval }
    pub fn allele_orientation_or(&self) -> Option<f64> { self.allele_orientation_odds_ratio }

    pub fn compute_stats(&mut self) -> (f64, Option<f64>) {
        (
            self.compute_allele_orientation_fisher_pval(),
            self.compute_allele_orientation_odds_ratio(),
        )
    }

    pub fn compute_allele_orientation_fisher_pval(&mut self) -> f64 {
        // Fisher test
        let exact_test = fishers_exact::fishers_exact(&self.allele_orientation_table).unwrap();
        let pval = exact_test.two_tail_pvalue;
        self.allele_orientation_fisher_pval = Some(pval);
        pval
    }
    pub fn compute_allele_orientation_odds_ratio(&mut self) -> Option<f64> {
        if self.minor_count() == 0 { return None }
        // Odds ratio of minor f/r / major f/r, 0.5 zero adjustment
        // k / reciprocal of the opposite treatment arm size
        let major_zero_adj: f64 = 0.5 / (self.minor_count() as f64);
        let minor_zero_adj: f64 = 0.5 / (self.major_count() as f64);
        let minor_odds: f64 = (self.allele_orientation_table[2] as f64 + minor_zero_adj) / (self.allele_orientation_table[3] as f64 + minor_zero_adj);
        let major_odds: f64 = (self.allele_orientation_table[0] as f64 + major_zero_adj) / (self.allele_orientation_table[1] as f64 + major_zero_adj);
        let or = minor_odds / major_odds;
        self.allele_orientation_odds_ratio = Some(or);
        Some(or)
    }

}
impl fmt::Display for PileupStats {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{major_base}:{major_count}:{major_freq:.3}\t{minor_base}:{minor_count}:{minor_freq:.3}\t{base_count}\t{cov}", 
            major_base=self.major_allele(),
            major_count=self.major_count(),
            major_freq=self.major_freq(),
            minor_base={if let Some(base) = self.minor_allele() { base.to_char() } else { '*' }},
            minor_count=self.minor_count(),
            minor_freq=self.minor_freq(),
            base_count=self.base_count(),
            cov=self.cov,
        )
    }
}



// Using rust htslib directly

// use rust_htslib::bam::pileup::{Indel, Pileup};

// #[derive(Serialize, Deserialize, Debug)]
// pub struct PileupRecord {
//     chrom: String,
//     pos: u64,
//     f_depth: u64,
//     r_depth: u64,
//     bases: Vec<char>,
//     bq: Vec<u8>,
//     mq: Vec<u8>,
//     indels: HashMap<String, u64>,
// }
// impl PileupRecord {
//     pub fn chrom(&self) -> &str { &self.chrom }
//     pub fn pos(&self) -> u64 { self.pos }
//     pub fn depth(&self) -> (u64, u64) { (self.f_depth, self.r_depth) }
//     pub fn bases(&self) -> &Vec<char> { &self.bases }
//     pub fn bq(&self) -> &Vec<u8> { &self.bq }
//     pub fn mq(&self) -> &Vec<u8> { &self.mq }
//     pub fn indels(&self) -> &HashMap<String, u64> { &self.indels }

//     pub fn base_counts(&self) -> HashMap<char, u64> {
//         let (mut f_a, mut f_c, mut f_g, mut f_t, mut f_n) = (0, 0, 0, 0, 0);
//         let (mut r_a, mut r_c, mut r_g, mut r_t, mut r_n) = (0, 0, 0, 0, 0);
//         self.bases.iter().for_each(|b| {
//             match b {
//                 'A' => f_a += 1,
//                 'C' => f_c += 1,
//                 'G' => f_g += 1,
//                 'T' => f_t += 1,
//                 'N' => f_n += 1,
//                 'a' => r_a += 1,
//                 'c' => r_c += 1,
//                 'g' => r_g += 1,
//                 't' => r_t += 1,
//                 'n' => r_n += 1,
//                 _ => ()
//             }
//         });
//         let bases: HashMap<char, u64> = [
//                 ('A', f_a), ('C', f_c), ('G', f_g), ('T', f_t), ('N', f_n),
//                 ('a', r_a), ('c', r_c), ('g', r_g), ('t', r_t), ('n', r_n)]
//             .iter()
//             .cloned()
//             .collect();
//         bases
//     }
//     pub fn indel_counts(&self) -> HashMap<String, u64> {
//         self.indels().to_owned()
//     }

//     pub fn total_depth(&self) -> u64 { self.f_depth + self.r_depth }
//     pub fn depth_ratio(&self) -> f32 { self.f_depth as f32 / self.r_depth as f32 }

//     pub fn mean_bq(&self) -> f32 {
//         self.bq.iter().map(|i| *i as i64).sum::<i64>() as f32 / self.bq.len() as f32
//     }
//     pub fn mean_mq(&self) -> f32 {
//         self.mq.iter().map(|i| *i as i64).sum::<i64>() as f32 / self.mq.len() as f32
//     }

//     pub fn new_empty(chrom: &str, pos: u64) -> PileupRecord {
//         PileupRecord {
//             chrom: chrom.to_owned(),
//             pos,
//             f_depth: 0,
//             r_depth: 0,
//             bases: Vec::new(),
//             bq: Vec::new(),
//             mq: Vec::new(),
//             indels: HashMap::new(),
//         }
//     }

//     pub fn from_pileup(chrom: &str, pileup: Pileup) -> PileupRecord {
//         let chrom = chrom.to_owned();
//         let pos = pileup.pos() as u64;
        
//         let mut bases: Vec<char> = Vec::new();
//         let mut bq: Vec<u8> = Vec::new();
//         let mut mq: Vec<u8> = Vec::new();
//         let indels: HashMap<String, u64> = HashMap::new();
        
//         let mut f_depth = 0;
//         let mut r_depth = 0;

//         // Iterate over reads covering this position
//         for aln in pileup.alignments() {
//             let rec = aln.record();

//             // Insertion or deletion
//             match aln.indel() {
//                 Indel::Ins(_) => (),
//                 Indel::Del(_) => (),
//                 Indel::None => ()
//             }
            
//             // Substitutions
//             if let Some(qpos) = aln.qpos() {
//                 // bases
//                 match aln.record().strand() {
//                     ReqStrand::Forward => {
//                         f_depth += 1;
//                         match rec.seq()[qpos] {
//                             65 => bases.push('A'),
//                             67 => bases.push('C'),
//                             71 => bases.push('G'),
//                             84 => bases.push('T'),
//                             78 => bases.push('N'),
//                             _ => (),
//                         }
//                     },
//                     ReqStrand::Reverse => {
//                         r_depth += 1;
//                         match rec.seq()[qpos] {
//                             65 => bases.push('a'),
//                             67 => bases.push('c'),
//                             71 => bases.push('g'),
//                             84 => bases.push('t'),
//                             78 => bases.push('n'),
//                             _ => (),
//                         }
//                     },
//                 };
//                 // bq
//                 bq.push(rec.qual()[qpos]);
//                 // mq
//                 mq.push(rec.mapq());
//             } 
//         }
//         if bases.len() != bq.len() || bases.len() != mq.len() { 
//             panic!("Unequal lengths: {} {} {}", bases.len(), bq.len(), mq.len())
//         }
//         PileupRecord {
//             chrom,
//             pos,
//             f_depth,
//             r_depth,
//             bases,
//             bq,
//             mq,
//             indels,
//         }
//     }

//     pub fn from_pileup_str(pileup: &str) -> PileupRecord {
        
//     }

//     pub fn quality_filter(&self, min_bq: u8, min_mq: u8) -> PileupRecord {
//         // Get positions that passed minimum quality scores
//         let bq_pos: BTreeSet<usize> = self.bq.iter().enumerate().filter(|(_, q)| **q >= min_bq).map(|(i, _)| i as usize).collect();
//         let mq_pos: BTreeSet<usize> = self.mq.iter().enumerate().filter(|(_, q)| **q >= min_mq).map(|(i, _)| i as usize).collect();
//         let target_pos: Vec<usize> = bq_pos.union(&mq_pos).cloned().collect();
//         let bases: Vec<char> = target_pos.iter().map(|&i| self.bases[i].to_owned()).collect();
//         let bq: Vec<u8> = target_pos.iter().map(|&i| self.bq[i].to_owned()).collect();
//         let mq: Vec<u8> = target_pos.iter().map(|&i| self.mq[i].to_owned()).collect();
//         // Recount forward and reverse depth
//         let mut f_depth = 0;
//         let mut r_depth = 0;
//         for b in bases.iter() {
//             match b {
//                 'A' => f_depth += 1,
//                 'C' => f_depth += 1,
//                 'G' => f_depth += 1,
//                 'T' => f_depth += 1,
//                 'N' => f_depth += 1,
//                 'a' => r_depth += 1,
//                 'c' => r_depth += 1,
//                 'g' => r_depth += 1,
//                 't' => r_depth += 1,
//                 'n' => r_depth += 1,
//                 _ => ()
//             }
//         }
//         PileupRecord {
//             chrom: self.chrom.to_owned(),
//             pos: self.pos,
//             f_depth,
//             r_depth,
//             bases,
//             bq,
//             mq,
//             indels: self.indels.clone(),
//         }
//     }    
// }

// impl fmt::Display for PileupRecord {
//     // This trait requires `fmt` with this exact signature.
//     fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
//         // Write strictly the first element into the supplied output
//         // stream: `f`. Returns `fmt::Result` which indicates whether the
//         // operation succeeded or failed. Note that `write!` uses syntax which
//         // is very similar to `println!`.
//         let bases: String = self.bases.iter().collect();
//         let indels: String = self.indels.iter().map(|(k, _)| k.to_string()).collect::<Vec<String>>().join(":");
//         write!(f, "{}:{} depth:{}:{} bases:{} indels:{}", 
//             self.chrom, self.pos, self.f_depth, self.r_depth, bases, indels)
//     }
// }

// pub enum BaseType {
//     A, C, G, T, N
// }

// pub enum Base {
//     A, C, G, T, N, a, c, g, t, n,
// }

// #[derive(Debug)]
// pub struct SiteSummary {
//     depth: usize,
//     ref_char: char,
//     f_a: usize,
//     f_c: usize,
//     f_g: usize,
//     f_t: usize,
//     f_n: usize,
//     r_a: usize,
//     r_c: usize,
//     r_g: usize,
//     r_t: usize,
//     r_n: usize,
//     indels: HashMap<String, usize>
// }

// impl SiteSummary {
//     pub fn new(s: &str, ref_char: char, depth: usize) -> SiteSummary {
//         let mut f_a: usize = 0;
//         let mut r_a: usize = 0;
//         let mut f_c: usize = 0;
//         let mut r_c: usize = 0;
//         let mut f_g: usize = 0;
//         let mut r_g: usize = 0;
//         let mut f_t: usize = 0;
//         let mut r_t: usize = 0;
//         let mut f_n: usize = 0;
//         let mut r_n: usize = 0;

//         let mut indel_dict: HashMap<String, usize> = HashMap::new();

//         let ref_char: char = ref_char.to_ascii_uppercase();

//         let mut i: usize = 0;
//         let s: Vec<char> = s.chars().collect();
//         while i < s.len() {
//             match s[i] {
//                 // Indicates start of read
//                 // Char after ^ indicates MAPQ score
//                 // Skip this char and next
//                 '^' => i += 2,
//                 // Count reference bases if ref is specified
//                 cc if cc == '.' => {
//                     match ref_char {
//                         'A' => {
//                             f_a += 1;
//                             i += 1;
//                         },
//                         'C' => {
//                             f_c += 1;
//                             i += 1;
//                         },
//                         'G' => {
//                             f_g += 1;
//                             i += 1;
//                         },
//                         'T' => {
//                             f_t += 1;
//                             i += 1;
//                         },
//                         'N' => {
//                             f_n += 1;
//                             i += 1;
//                         },
//                         _ => (),
//                     }
//                     i += 1;
//                 },
//                 ',' => {
//                     match ref_char {
//                         'A' => {
//                             r_a += 1;
//                             i += 1;
//                         },
//                         'C' => {
//                             r_c += 1;
//                             i += 1;
//                         },
//                         'G' => {
//                             r_g += 1;
//                             i += 1;
//                         },
//                         'T' => {
//                             r_t += 1;
//                             i += 1;
//                         },
//                         'N' => {
//                             r_n += 1;
//                             i += 1;
//                         },
//                         _ => (),
//                     }
//                     i += 1;
//                 },
//                 // Count bases
//                 'A' => {
//                     f_a += 1;
//                     i += 1;
//                 },
//                 'a' => {
//                     r_a += 1;
//                     i += 1;
//                 },
//                 'C' => {
//                     f_c += 1;
//                     i += 1;
//                 },
//                 'c' => {
//                     r_c += 1;
//                     i += 1;
//                 },
//                 'G' => {
//                     f_g += 1;
//                     i += 1;
//                 },
//                 'g' => {
//                     r_g += 1;
//                     i += 1;
//                 },
//                 'T' => {
//                     f_t += 1;
//                     i += 1;
//                 },
//                 't' => {
//                     r_t += 1;
//                     i += 1;
//                 },
//                 'N' => {
//                     f_n += 1;
//                     i += 1;
//                 },
//                 'n' => {
//                     r_n += 1;
//                     i += 1;
//                 },
//                 // Count indels
//                 cc if cc == '+' || cc == '-' => {
//                     // get indel length
//                     let mut indel_len_str: Vec<&char> = Vec::new();
//                     for c in s[(i+1)..].iter() {
//                         if !c.is_numeric() {
//                             break;
//                         }
//                         indel_len_str.push(c);
//                     }
//                     // get indel characters
//                     let indel_len = indel_len_str.into_iter().collect::<String>().parse::<usize>().unwrap();
//                     let indel_ident = s[(i+2)..(i+2+indel_len)].iter().collect::<String>();
//                     let indel_ident = format!("{}{}{}", cc, indel_len, indel_ident);

//                     // upsert
//                     *indel_dict.entry(indel_ident).or_insert(0) += 1;

//                     i += 2 + indel_len;
//                 },
//                 // skip any other character
//                 _ => i += 1,
//             }
//         }

//         SiteSummary{
//             depth,
//             ref_char,
//             f_a, f_c, f_g, f_t, f_n,
//             r_a, r_c, r_g, r_t, r_n,
//             indels: indel_dict,
//         }
//     }

//     pub fn pool(summaries: &Vec<&SiteSummary>) -> SiteSummary {
//         if summaries.len() == 0 {
//             panic!("Empty summaries vector")
//         }
//         let ref_char: char = summaries[0].ref_char;
//         let mut depth: usize = 0;
//         let mut f_a: usize = 0;
//         let mut f_c: usize = 0;
//         let mut f_g: usize = 0;
//         let mut f_t: usize = 0;
//         let mut f_n: usize = 0;
//         let mut r_a: usize = 0;
//         let mut r_c: usize = 0;
//         let mut r_g: usize = 0;
//         let mut r_t: usize = 0;
//         let mut r_n: usize = 0;
//         let mut indels: HashMap<String, usize> = HashMap::new();

//         summaries.into_iter().for_each(|s| {
//             depth += s.depth;
//             f_a += s.f_a;
//             f_c += s.f_c;
//             f_g += s.f_g;
//             f_t += s.f_t;
//             f_n += s.f_n;
//             r_a += s.r_a;
//             r_c += s.r_c;
//             r_g += s.r_g;
//             r_t += s.r_t;
//             r_n += s.r_n;
//             s.indels.iter().for_each(|(k, v)| {
//                 *indels.entry(k.to_owned()).or_insert(0) += v;
//             });
//         });

//         SiteSummary{
//             depth,
//             ref_char,
//             f_a, f_c, f_g, f_t, f_n,
//             r_a, r_c, r_g, r_t, r_n,
//             indels,
//         }
//     }

//     pub fn to_hashmaps(&self) -> (HashMap<String, usize>, HashMap<String, usize>) {
//         let subs_dict: HashMap<String, usize> = [
//             ("A".to_owned(), self.f_a),
//             ("C".to_owned(), self.f_c),
//             ("G".to_owned(), self.f_g),
//             ("T".to_owned(), self.f_t),
//             ("N".to_owned(), self.f_n),
//             ("a".to_owned(), self.r_a),
//             ("c".to_owned(), self.r_c),
//             ("g".to_owned(), self.r_g),
//             ("t".to_owned(), self.r_t),
//             ("n".to_owned(), self.r_n),
//         ].iter().cloned().collect();
//         let indel_dict: HashMap<String, usize> = self.indels.clone();

//         (subs_dict, indel_dict)
//     }

//     pub fn count_of_any(&self, base: char) -> usize {
//         let base: char = base.to_ascii_uppercase();
//         match base {
//             'A' => self.f_a + self.r_a,
//             'C' => self.f_c + self.r_c,
//             'G' => self.f_g + self.r_g,
//             'T' => self.f_t + self.r_t,
//             'N' => self.f_n + self.r_n,
//             _ => 0,
//         }
//     }

//     pub fn count_of(&self, base: char) -> usize {
//         match base {
//             'A' => self.f_a,
//             'C' => self.f_c,
//             'G' => self.f_g,
//             'T' => self.f_t,
//             'N' => self.f_n,
//             'a' => self.r_a,
//             'c' => self.r_c,
//             'g' => self.r_g,
//             't' => self.r_t,
//             'n' => self.r_n,
//             _ => 0,
//         }
//     }

//     pub fn freq_of_any(&self, base: char) -> f32 {
//         (self.count_of_any(base) as f32) / (self.depth as f32)
//     }

//     pub fn freq_of(&self, base: char) -> f32 {
//         (self.count_of(base) as f32) / (self.depth as f32)
//     }

// }

// #[derive(Debug)]
// pub struct SNV {
//     major_allele: char,
//     major_freq: f32,
//     minor_allele: char,
//     minor_freq: f32,
// }

// impl SNV {
//     pub fn basic_call(s: &SiteSummary) -> Option<SNV> {
//         let sorted_bases: Vec<(char, usize)> = SNV::sort_alleles(s);
//         if sorted_bases[0].1 == s.depth {
//             None
//         } else {
//             Some(SNV{
//                 major_allele: sorted_bases[0].0,
//                 major_freq: (sorted_bases[0].1 as f32) / (s.depth as f32),
//                 minor_allele: sorted_bases[1].0,
//                 minor_freq: (sorted_bases[1].1 as f32) / (s.depth as f32),
//             })
//         }
//     }

//     pub fn threshold_call(s: &SiteSummary, threshold_cnt: Option<usize>, threshold_pct: Option<f32>) -> Option<SNV> {
//         let sorted_bases: Vec<(char, usize)> = SNV::sort_alleles(s);
//         if sorted_bases[0].1 == s.depth {
//             None
//         } else {
//             let minor_cnt = sorted_bases[1].1;
//             let minor_freq = (minor_cnt as f32) / (s.depth as f32);
//             let threshold_cnt = threshold_cnt.unwrap_or(0);
//             let threshold_pct = threshold_pct.unwrap_or(0.);
//             if minor_cnt > threshold_cnt {
//                 if minor_freq > threshold_pct {
//                     return Some(SNV{
//                         major_allele: sorted_bases[0].0,
//                         major_freq: (sorted_bases[0].1 as f32) / (s.depth as f32),
//                         minor_allele: sorted_bases[1].0,
//                         minor_freq: (sorted_bases[1].1 as f32) / (s.depth as f32),
//                     });
//                 } else {
//                     None
//                 }
//             } else {
//                 None
//             }
//         }
//     }

//     pub fn similar_to(&self, other: SNV) -> bool {
//         if (self.major_allele == other.major_allele) &&
//            (self.minor_allele == other.minor_allele) {
//             true
//         } else {
//             false
//         }
//     }

//     fn sort_alleles(s: &SiteSummary) -> Vec<(char, usize)> {
//         let cnt_a: usize = s.count_of_any('A');
//         let cnt_c: usize = s.count_of_any('C');
//         let cnt_g: usize = s.count_of_any('G');
//         let cnt_t: usize = s.count_of_any('T');

//         // sort bases from largest freq
//         // top two will be major and minor respectively
//         let mut sorted_bases: Vec<(char, usize)> = vec![
//             ('A', cnt_a),
//             ('C', cnt_c),
//             ('G', cnt_g),
//             ('T', cnt_t),
//         ];
//         sorted_bases.sort_unstable_by(|(_, v1), (_, v2)| v2.cmp(v1));
//         sorted_bases
//     }
// }

// #[derive(Debug)]
// pub enum SNVCallType {
//     True,
//     ProbablyTrue,
//     Undecided,
//     ProbablyFalse,
//     False,
//     InsufficientData,
// }


// #[derive(Debug)]
// pub struct PairedSNVCall {
//     // target alleles and frequency
//     target_major_allele: char,
//     target_major_freq: f32,
//     target_minor_allele: Option<char>,
//     target_minor_freq: f32,
//     // control alleles and frequency
//     control_major_allele: char,
//     control_major_freq: f32,
//     control_minor_allele: Option<char>,
//     control_minor_freq: f32,
//     call: SNVCallType,
// }

// impl PairedSNVCall {
//     pub fn basic_call(target: &SiteSummary, control: &SiteSummary) -> Option<PairedSNVCall> {
//         let sorted_target: Vec<(char, usize)> = SNV::sort_alleles(target);
//         let sorted_control: Vec<(char, usize)> = SNV::sort_alleles(target);

//         // Target and control have the same major allele
//         if sorted_target[0].0 == sorted_control[0].0 {
//             // Target has no variant
//             if sorted_target[1].1 == 0 {
//                 // Control also has no variants
//                 if sorted_control[1].1 == 0 {
//                     return None
//                 }
//                 // Target has no variants but control has variants
//                 // Since this is in the context of the target, it is false
//                 return Some(PairedSNVCall{
//                     target_major_allele: sorted_target[0].0,
//                     target_major_freq: 1.0,
//                     target_minor_allele: None,
//                     target_minor_freq: 0.0,
//                     // control alleles and frequency
//                     control_major_allele: sorted_control[0].0,
//                     control_major_freq: (sorted_control[0].1 as f32) / (control.depth as f32),
//                     control_minor_allele: Some(sorted_control[1].0),
//                     control_minor_freq: (sorted_control[1].1 as f32) / (control.depth as f32),
//                     call: SNVCallType::False,
//                 })
//             }
//             // Target has a variant
//             if sorted_control[1].1 == 0 {
//                 // Control has 0 variants so target SNV is likely true
//                 return Some(PairedSNVCall{
//                     target_major_allele: sorted_target[0].0,
//                     target_major_freq: (sorted_target[0].1 as f32) / (target.depth as f32),
//                     target_minor_allele: Some(sorted_target[1].0),
//                     target_minor_freq: (sorted_target[1].1 as f32) / (target.depth as f32),
//                     // control alleles and frequency
//                     control_major_allele: sorted_control[0].0,
//                     control_major_freq: 1.0,
//                     control_minor_allele: None,
//                     control_minor_freq: 0.0,
//                     call: SNVCallType::True,
//                 })
//             } else if sorted_target[1].0 == sorted_control[1].0 {
//                 // Control contains at least one variant and 
//                 // it is the same variant then target SNV is probably false
//                 return Some(PairedSNVCall{
//                     target_major_allele: sorted_target[0].0,
//                     target_major_freq: (sorted_target[0].1 as f32) / (target.depth as f32),
//                     target_minor_allele: Some(sorted_target[1].0),
//                     target_minor_freq: (sorted_target[1].1 as f32) / (target.depth as f32),
//                     // control alleles and frequency
//                     control_major_allele: sorted_control[0].0,
//                     control_major_freq: (sorted_control[0].1 as f32) / (control.depth as f32),
//                     control_minor_allele: Some(sorted_control[1].0),
//                     control_minor_freq: (sorted_control[1].1 as f32) / (control.depth as f32),
//                     call: SNVCallType::ProbablyFalse,
//                 })
//             }
//             // Control contains a different variant which indicates probably multiple alleles
//             // or some artifact of sequencing or mapping
//             // then the target SNV is considered false
//             return Some(PairedSNVCall{
//                 target_major_allele: sorted_target[0].0,
//                 target_major_freq: (sorted_target[0].1 as f32) / (target.depth as f32),
//                 target_minor_allele: Some(sorted_target[1].0),
//                 target_minor_freq: (sorted_target[1].1 as f32) / (target.depth as f32),
//                 // control alleles and frequency
//                 control_major_allele: sorted_control[0].0,
//                 control_major_freq: (sorted_control[0].1 as f32) / (control.depth as f32),
//                 control_minor_allele: Some(sorted_control[1].0),
//                 control_minor_freq: (sorted_control[1].1 as f32) / (control.depth as f32),
//                 call: SNVCallType::False,
//             })  
//         }
//         // Target and control have different major alleles
//         // Target has no variant
//         if sorted_target[1].1 == 0 {
//             // Control also has no variants
//             if sorted_control[1].1 == 0 {
//                 return None
//             } 
//             // Target has no variants but control has variants
//             // Since this is in the context of the target, it is false
//             return Some(PairedSNVCall{
//                 target_major_allele: sorted_target[0].0,
//                 target_major_freq: 1.0,
//                 target_minor_allele: None,
//                 target_minor_freq: 0.0,
//                 // control alleles and frequency
//                 control_major_allele: sorted_control[0].0,
//                 control_major_freq: (sorted_control[0].1 as f32) / (control.depth as f32),
//                 control_minor_allele: Some(sorted_control[1].0),
//                 control_minor_freq: (sorted_control[1].1 as f32) / (control.depth as f32),
//                 call: SNVCallType::False,
//             }) 
//         }
//         // Target and control have different major alleles
//         // Target and control have minor alleles
//         return Some(PairedSNVCall{
//             target_major_allele: sorted_target[0].0,
//             target_major_freq: (sorted_target[0].1 as f32) / (target.depth as f32),
//             target_minor_allele: Some(sorted_target[1].0),
//             target_minor_freq: (sorted_target[1].1 as f32) / (target.depth as f32),
//             // control alleles and frequency
//             control_major_allele: sorted_control[0].0,
//             control_major_freq: (sorted_control[0].1 as f32) / (control.depth as f32),
//             control_minor_allele: Some(sorted_control[1].0),
//             control_minor_freq: (sorted_control[1].1 as f32) / (control.depth as f32),
//             call: SNVCallType::False,
//         })
//     }

//     pub fn threshold_call(target: &SiteSummary, control: &SiteSummary,
//             threshold_target_cnt: Option<usize>, threshold_target_freq: Option<f32>,
//             threshold_control_cnt: Option<usize>, threshold_control_freq: Option<f32>) -> Option<PairedSNVCall> {

//         let threshold_target_cnt = threshold_target_cnt.unwrap_or(0);
//         let threshold_target_freq = threshold_target_freq.unwrap_or(0.);
//         let threshold_control_cnt = threshold_control_cnt.unwrap_or(0);
//         let threshold_control_freq = threshold_control_freq.unwrap_or(0.);

//         let sorted_target: Vec<(char, usize)> = SNV::sort_alleles(target);
//         let sorted_control: Vec<(char, usize)> = SNV::sort_alleles(target);

//         let target_maj_allele = sorted_target[0].0;
//         let target_maj_cnt = sorted_target[0].1;
//         let target_maj_freq = (sorted_target[0].1 as f32) / (target.depth as f32);
//         let target_min_allele = sorted_target[1].0;
//         let target_min_cnt = sorted_target[1].1;
//         let target_min_freq = (sorted_target[1].1 as f32) / (target.depth as f32);

//         let control_maj_allele = sorted_control[0].0;
//         let control_maj_cnt = sorted_control[0].1;
//         let control_maj_freq = (sorted_control[0].1 as f32) / (control.depth as f32);
//         let control_min_allele = sorted_control[1].0;
//         let control_min_cnt = sorted_control[1].1;
//         let control_min_freq = (sorted_control[1].1 as f32) / (control.depth as f32);

//         // Target and control have the same major allele
//         if target_maj_allele == control_maj_allele {
//             // Target has effectively no variant
//             // Variant count is 0 or less than threshold count
//             // Variant freq is 0 or less than threshold freq
//             if (target_min_cnt <= threshold_target_cnt) || (target_min_freq <= threshold_target_freq) {
//                 // Control also has effectively no variants
//                 // Variant count is 0 or less than threshold count
//                 // Variant freq is 0 or less than threshold freq
//                 if (control_min_cnt <= threshold_control_cnt) || (control_min_freq <= threshold_control_freq) {
//                     return None
//                 }
//                 // Target has no variants but control has variants greater than the thresholds
//                 // Since this is in the context of the target, it is false
//                 return Some(PairedSNVCall{
//                     target_major_allele: sorted_target[0].0,
//                     target_major_freq: 1.0,
//                     target_minor_allele: None,
//                     target_minor_freq: 0.0,
//                     // control alleles and frequency
//                     control_major_allele: sorted_control[0].0,
//                     control_major_freq: (sorted_control[0].1 as f32) / (control.depth as f32),
//                     control_minor_allele: Some(sorted_control[1].0),
//                     control_minor_freq: (sorted_control[1].1 as f32) / (control.depth as f32),
//                     call: SNVCallType::False,
//                 })
//             }
//             // Target has a variant
//             // Control has effectively no variants so target SNV is likely true
//             if (control_min_cnt <= threshold_control_cnt) || (control_min_freq <= threshold_control_freq) {
//                 return Some(PairedSNVCall{
//                     target_major_allele: sorted_target[0].0,
//                     target_major_freq: (sorted_target[0].1 as f32) / (target.depth as f32),
//                     target_minor_allele: Some(sorted_target[1].0),
//                     target_minor_freq: (sorted_target[1].1 as f32) / (target.depth as f32),
//                     // control alleles and frequency
//                     control_major_allele: sorted_control[0].0,
//                     control_major_freq: 1.0,
//                     control_minor_allele: None,
//                     control_minor_freq: 0.0,
//                     call: SNVCallType::True,
//                 })
//             } else if target_min_allele == control_min_allele {
//                 // Control contains at least one variant and 
//                 // it is the same variant then target SNV is probably false
//                 return Some(PairedSNVCall{
//                     target_major_allele: sorted_target[0].0,
//                     target_major_freq: (sorted_target[0].1 as f32) / (target.depth as f32),
//                     target_minor_allele: Some(sorted_target[1].0),
//                     target_minor_freq: (sorted_target[1].1 as f32) / (target.depth as f32),
//                     // control alleles and frequency
//                     control_major_allele: sorted_control[0].0,
//                     control_major_freq: (sorted_control[0].1 as f32) / (control.depth as f32),
//                     control_minor_allele: Some(sorted_control[1].0),
//                     control_minor_freq: (sorted_control[1].1 as f32) / (control.depth as f32),
//                     call: SNVCallType::ProbablyFalse,
//                 })
//             }
//             // Control contains a different variant which indicates probably multiple alleles
//             // or some artifact of sequencing or mapping
//             // then the target SNV is considered false
//             return Some(PairedSNVCall{
//                 target_major_allele: sorted_target[0].0,
//                 target_major_freq: (sorted_target[0].1 as f32) / (target.depth as f32),
//                 target_minor_allele: Some(sorted_target[1].0),
//                 target_minor_freq: (sorted_target[1].1 as f32) / (target.depth as f32),
//                 // control alleles and frequency
//                 control_major_allele: sorted_control[0].0,
//                 control_major_freq: (sorted_control[0].1 as f32) / (control.depth as f32),
//                 control_minor_allele: Some(sorted_control[1].0),
//                 control_minor_freq: (sorted_control[1].1 as f32) / (control.depth as f32),
//                 call: SNVCallType::False,
//             })  
//         }
//         // Target and control have different major alleles
//         // Target has effectively no variant
//         if (target_min_cnt <= threshold_target_cnt) || (target_min_freq <= threshold_target_freq) {
//             // Control also has effectively no variants
//             if (control_min_cnt <= threshold_control_cnt) || (control_min_freq <= threshold_control_freq) {
//                 return None
//             } 
//             // Target has no variants but control has variants
//             // Since this is in the context of the target, it is false
//             return Some(PairedSNVCall{
//                 target_major_allele: sorted_target[0].0,
//                 target_major_freq: 1.0,
//                 target_minor_allele: None,
//                 target_minor_freq: 0.0,
//                 // control alleles and frequency
//                 control_major_allele: sorted_control[0].0,
//                 control_major_freq: (sorted_control[0].1 as f32) / (control.depth as f32),
//                 control_minor_allele: Some(sorted_control[1].0),
//                 control_minor_freq: (sorted_control[1].1 as f32) / (control.depth as f32),
//                 call: SNVCallType::False,
//             }) 
//         }
//         // Target and control have different major alleles
//         // Target and control have minor alleles
//         return Some(PairedSNVCall{
//             target_major_allele: sorted_target[0].0,
//             target_major_freq: (sorted_target[0].1 as f32) / (target.depth as f32),
//             target_minor_allele: Some(sorted_target[1].0),
//             target_minor_freq: (sorted_target[1].1 as f32) / (target.depth as f32),
//             // control alleles and frequency
//             control_major_allele: sorted_control[0].0,
//             control_major_freq: (sorted_control[0].1 as f32) / (control.depth as f32),
//             control_minor_allele: Some(sorted_control[1].0),
//             control_minor_freq: (sorted_control[1].1 as f32) / (control.depth as f32),
//             call: SNVCallType::False,
//         })
//     }

//     pub fn similar_to(&self, other: &PairedSNVCall) -> bool {
//         if (self.target_major_allele == other.target_major_allele) &&
//            (self.control_major_allele == other.control_major_allele) &&
//            (self.target_minor_allele == other.target_minor_allele) &&
//            (self.control_minor_allele == other.control_minor_allele) {
//             true
//         } else {
//             false
//         }
//     }
// }

// #[derive(Debug)]
// pub struct MultiPairedSNVCall{
//     // Sample values are stored in the same order as the sample names
//     sample: Vec<String>,
//     sample_major_allele: char,
//     sample_major_freq: Vec<f32>,
//     sample_minor_allele: Option<char>,
//     sample_minor_freq: Vec<f32>,
//     // control alleles and frequency
//     control_major_allele: char,
//     control_major_freq: f32,
//     control_minor_allele: Option<char>,
//     control_minor_freq: f32,
//     call: SNVCallType,
// }

// impl MultiPairedSNVCall {
//     pub fn pooled_basic_call(samples: &Vec<&SiteSummary>, control: &SiteSummary) -> Option<PairedSNVCall> {
//         // pool SiteSummary vec into one and call as PairedSNVCall
//         let pooled_target = SiteSummary::pool(samples);
//         PairedSNVCall::basic_call(&pooled_target, control)
//     }

//     pub fn pooled_threshold_call(samples: &Vec<&SiteSummary>, control: &SiteSummary,
//             threshold_target_cnt: Option<usize>, threshold_target_freq: Option<f32>,
//             threshold_control_cnt: Option<usize>, threshold_control_freq: Option<f32>) -> Option<PairedSNVCall> {
//         // pool SiteSummary vec into one and call as PairedSNVCall
//         let pooled_target = SiteSummary::pool(samples);
//         PairedSNVCall::threshold_call(&pooled_target, control, 
//             threshold_target_cnt, threshold_target_freq,
//             threshold_control_cnt, threshold_control_freq)
//     }
    
//     pub fn basic_call(samples: &Vec<&SiteSummary>, control: &SiteSummary) -> Option<MultiPairedSNVCall> {
//         // Make pooled basic call first to establish major-minor alleles
//         if let Some(pc) = MultiPairedSNVCall::pooled_basic_call(samples, control) {
//             // Target and control do not have the same major allele
//             if pc.target_major_allele != pc.control_major_allele {
//                 return Some(MultiPairedSNVCall{
//                     // Sample values are stored in the same order as the sample names
//                     sample: Vec::new(),
//                     sample_major_allele: pc.target_major_allele,
//                     sample_major_freq: samples.iter().map(|s| s.freq_of_any(pc.target_major_allele)).collect(),
//                     sample_minor_allele: pc.target_minor_allele,
//                     sample_minor_freq: match pc.target_minor_allele {
//                         Some(min_allele) => samples.iter().map(|s| s.freq_of_any(min_allele)).collect(),
//                         None => (0..samples.len()).map(|_| 0.0).collect(),
//                     },
//                     // control alleles and frequency
//                     control_major_allele: pc.control_major_allele,
//                     control_major_freq: pc.control_major_freq,
//                     control_minor_allele: pc.control_minor_allele,
//                     control_minor_freq: pc.control_minor_freq,
//                     call: SNVCallType::False,
//                 })
//             }
//             // Target and control major alleles are the same and
//             // Target has at least 1 variant read 
//             if let Some(target_minor_allele) = pc.target_minor_allele {
//                 // Iterate over
//                 // let mut sample: Vec<String> = Vec::with_capacity(samples.len());
//                 let mut sample_major_freq: Vec<f32> = Vec::with_capacity(samples.len());
//                 let mut sample_minor_freq: Vec<f32> = Vec::with_capacity(samples.len());

//                 let mut min_allele_present: Vec<usize> = Vec::new();

//                 samples.iter().enumerate().for_each(|(i, s)| {
//                     let maj_freq = s.freq_of_any(pc.target_major_allele);
//                     let min_freq = s.freq_of_any(target_minor_allele);
//                     if min_freq > 0.0 {
//                         min_allele_present.push(i);
//                     }
//                     sample_major_freq.push(maj_freq);
//                     sample_minor_freq.push(min_freq);
//                 });
//                 let call_type: SNVCallType = match min_allele_present.len() {
//                     1 => SNVCallType::ProbablyFalse,
//                     2 => SNVCallType::ProbablyTrue,
//                     _ => SNVCallType::True,
//                 };
//                 return Some(MultiPairedSNVCall{
//                     // Sample values are stored in the same order as the sample names
//                     sample: Vec::new(),
//                     sample_major_allele: pc.target_major_allele,
//                     sample_major_freq,
//                     sample_minor_allele: Some(target_minor_allele),
//                     sample_minor_freq,
//                     // control alleles and frequency
//                     control_major_allele: pc.control_major_allele,
//                     control_major_freq: pc.control_major_freq,
//                     control_minor_allele: pc.control_minor_allele,
//                     control_minor_freq: pc.control_minor_freq,
//                     call: call_type,
//                 })
//             }
//             // Target has no minor but control has a minor allele
//             return Some(MultiPairedSNVCall{
//                 // Sample values are stored in the same order as the sample names
//                 sample: Vec::new(),
//                 sample_major_allele: pc.target_major_allele,
//                 sample_major_freq: (0..samples.len()).map(|_| 1.0).collect(),
//                 sample_minor_allele: None,
//                 sample_minor_freq: (0..samples.len()).map(|_| 0.0).collect(),
//                 // control alleles and frequency
//                 control_major_allele: pc.control_major_allele,
//                 control_major_freq: pc.control_major_freq,
//                 control_minor_allele: pc.control_minor_allele,
//                 control_minor_freq: pc.control_minor_freq,
//                 call: SNVCallType::False,
//             })
//         }
//         // No SNV call in pooled samples so there is no SNV, no further calling necessary
//         None
//     }

//     pub fn threshold_call(samples: &Vec<&SiteSummary>, control: &SiteSummary,
//             threshold_target_cnt: Option<usize>, threshold_target_freq: Option<f32>,
//             threshold_control_cnt: Option<usize>, threshold_control_freq: Option<f32>) -> Option<MultiPairedSNVCall> {
//         // Make pooled basic call first to establish major-minor alleles
//         if let Some(pc) = MultiPairedSNVCall::pooled_threshold_call(
//             samples, control,
//             threshold_target_cnt, threshold_target_freq,
//             threshold_control_cnt, threshold_control_freq) {
//             // Target and control do not have the same major allele
//             if pc.target_major_allele != pc.control_major_allele {
//                 return Some(MultiPairedSNVCall{
//                     // Sample values are stored in the same order as the sample names
//                     sample: Vec::new(),
//                     sample_major_allele: pc.target_major_allele,
//                     sample_major_freq: samples.iter().map(|s| s.freq_of_any(pc.target_major_allele)).collect(),
//                     sample_minor_allele: pc.target_minor_allele,
//                     sample_minor_freq: match pc.target_minor_allele {
//                         Some(min_allele) => samples.iter().map(|s| s.freq_of_any(min_allele)).collect(),
//                         None => (0..samples.len()).map(|_| 0.0).collect(),
//                     },
//                     // control alleles and frequency
//                     control_major_allele: pc.control_major_allele,
//                     control_major_freq: pc.control_major_freq,
//                     control_minor_allele: pc.control_minor_allele,
//                     control_minor_freq: pc.control_minor_freq,
//                     call: SNVCallType::False,
//                 })
//             }
//             // Target and control major alleles are the same and
//             // Target has at least 1 variant read 
//             if let Some(target_minor_allele) = pc.target_minor_allele {
//                 // Iterate over
//                 // let mut sample: Vec<String> = Vec::with_capacity(samples.len());
//                 let mut sample_major_freq: Vec<f32> = Vec::with_capacity(samples.len());
//                 let mut sample_minor_freq: Vec<f32> = Vec::with_capacity(samples.len());

//                 let mut min_allele_present: Vec<usize> = Vec::new();
//                 let threshold_target_freq = threshold_target_freq.unwrap_or(0.0);
//                 let threshold_target_cnt = threshold_target_cnt.unwrap_or(0);

//                 samples.iter().enumerate().for_each(|(i, s)| {
//                     let maj_freq = s.freq_of_any(pc.target_major_allele);
//                     let min_freq = s.freq_of_any(target_minor_allele);
//                     let min_cnt = s.count_of_any(target_minor_allele);
//                     if (min_freq > threshold_target_freq) && (min_cnt > threshold_target_cnt) {
//                         min_allele_present.push(i);
//                     }
//                     sample_major_freq.push(maj_freq);
//                     sample_minor_freq.push(min_freq);
//                 });
//                 let call_type: SNVCallType = match min_allele_present.len() {
//                     1 => SNVCallType::ProbablyFalse,
//                     2 => SNVCallType::ProbablyTrue,
//                     _ => SNVCallType::True,
//                 };
//                 return Some(MultiPairedSNVCall{
//                     // Sample values are stored in the same order as the sample names
//                     sample: Vec::new(),
//                     sample_major_allele: pc.target_major_allele,
//                     sample_major_freq,
//                     sample_minor_allele: Some(target_minor_allele),
//                     sample_minor_freq,
//                     // control alleles and frequency
//                     control_major_allele: pc.control_major_allele,
//                     control_major_freq: pc.control_major_freq,
//                     control_minor_allele: pc.control_minor_allele,
//                     control_minor_freq: pc.control_minor_freq,
//                     call: call_type,
//                 })
//             }
//             // Target has no minor but control has a minor allele
//             return Some(MultiPairedSNVCall{
//                 // Sample values are stored in the same order as the sample names
//                 sample: Vec::new(),
//                 sample_major_allele: pc.target_major_allele,
//                 sample_major_freq: (0..samples.len()).map(|_| 1.0).collect(),
//                 sample_minor_allele: None,
//                 sample_minor_freq: (0..samples.len()).map(|_| 0.0).collect(),
//                 // control alleles and frequency
//                 control_major_allele: pc.control_major_allele,
//                 control_major_freq: pc.control_major_freq,
//                 control_minor_allele: pc.control_minor_allele,
//                 control_minor_freq: pc.control_minor_freq,
//                 call: SNVCallType::False,
//             })
//         }
//         // No SNV call in pooled samples so there is no SNV, no further calling necessary
//         None
//     }

//     // pub fn similar_to(&self, other: &MultiPairedSNVCall) -> bool {
//     //     false
//     // }

// }