use std::collections::{HashMap, BTreeSet};
use rust_htslib::bam::pileup::{Indel, Pileup};
use bio_types::strand::ReqStrand;


#[derive(Debug)]
pub struct PileupRecord {
    chrom: String,
    pos: u64,
    f_depth: u64,
    r_depth: u64,
    bases: Vec<char>,
    bq: Vec<u8>,
    mq: Vec<u8>,
    indels: HashMap<String, u64>,
}
impl PileupRecord {
    pub fn from_pileup(chrom: &str, pileup: Pileup) -> PileupRecord {
        let chrom = chrom.to_owned();
        let pos = pileup.pos() as u64;
        
        let mut bases: Vec<char> = Vec::new();
        let mut bq: Vec<u8> = Vec::new();
        let mut mq: Vec<u8> = Vec::new();
        let mut indels: HashMap<String, u64> = HashMap::new();
        
        let mut f_depth = 0;
        let mut r_depth = 0;

        // Iterate over reads covering this position
        for aln in pileup.alignments() {
            let rec = aln.record();

            // Insertion or deletion
            match aln.indel() {
                Indel::Ins(len) => (),
                Indel::Del(len) => (),
                Indel::None => ()
            }
            
            // Substitutions
            if let Some(qpos) = aln.qpos() {
                // bases
                match aln.record().strand() {
                    ReqStrand::Forward => {
                        f_depth += 1;
                        match rec.seq()[qpos] {
                            65 => bases.push('A'),
                            67 => bases.push('C'),
                            71 => bases.push('G'),
                            84 => bases.push('T'),
                            78 => bases.push('N'),
                            _ => (),
                        }
                    },
                    ReqStrand::Reverse => {
                        r_depth += 1;
                        match rec.seq()[qpos] {
                            65 => bases.push('a'),
                            67 => bases.push('c'),
                            71 => bases.push('g'),
                            84 => bases.push('t'),
                            78 => bases.push('n'),
                            _ => (),
                        }
                    },
                };
                // bq
                bq.push(rec.qual()[qpos]);
                // mq
                mq.push(rec.mapq());
            } 
        }
        if bases.len() != bq.len() || bases.len() != mq.len() { 
            panic!("Unequal lengths: {} {} {}", bases.len(), bq.len(), mq.len())
        }
        PileupRecord {
            chrom,
            pos,
            f_depth,
            r_depth,
            bases,
            bq,
            mq,
            indels,
        }
    }

    pub fn quality_filter(&self, min_bq: u8, min_mq: u8) -> PileupRecord {
        let bq_pos: BTreeSet<usize> = self.bq.iter().enumerate().filter(|(i, q)| **q >= min_bq).map(|(i, _)| i as usize).collect();
        let mq_pos: BTreeSet<usize> = self.mq.iter().enumerate().filter(|(i, q)| **q >= min_mq).map(|(i, _)| i as usize).collect();
        let target_pos: Vec<usize> = bq_pos.union(&mq_pos).cloned().collect();
        let bases: Vec<char> = target_pos.iter().map(|&i| self.bases[i].to_owned()).collect();
        let bq: Vec<u8> = target_pos.iter().map(|&i| self.bq[i].to_owned()).collect();
        let mq: Vec<u8> = target_pos.iter().map(|&i| self.mq[i].to_owned()).collect();

        let mut f_depth = 0;
        let mut r_depth = 0;
        for b in bases.iter() {
            match b {
                'A' => f_depth += 1,
                'C' => f_depth += 1,
                'G' => f_depth += 1,
                'T' => f_depth += 1,
                'N' => f_depth += 1,
                'a' => r_depth += 1,
                'c' => r_depth += 1,
                'g' => r_depth += 1,
                't' => r_depth += 1,
                'n' => r_depth += 1,
                _ => ()
            }
        }
        PileupRecord {
            chrom: self.chrom.to_owned(),
            pos: self.pos,
            f_depth,
            r_depth,
            bases,
            bq,
            mq,
            indels: self.indels.clone(),
        }
    }

    pub fn depth_ratio(&self) -> f64 {
        self.f_depth as f64 / self.r_depth as f64
    }

    pub fn mean_bq(&self) -> f64 {
        self.bq.iter().map(|i| *i as i64).sum::<i64>() as f64 / self.bq.len() as f64
    }
    pub fn mean_mq(&self) -> f64 {
        self.mq.iter().map(|i| *i as i64).sum::<i64>() as f64 / self.mq.len() as f64
    }
}

pub enum BaseType {
    A, C, G, T, N
}

pub enum Base {
    A, C, G, T, N, a, c, g, t, n,
}

#[derive(Debug)]
pub struct SiteSummary {
    depth: usize,
    ref_char: char,
    f_a: usize,
    f_c: usize,
    f_g: usize,
    f_t: usize,
    f_n: usize,
    r_a: usize,
    r_c: usize,
    r_g: usize,
    r_t: usize,
    r_n: usize,
    indels: HashMap<String, usize>
}

impl SiteSummary {
    pub fn new(s: &str, ref_char: char, depth: usize) -> SiteSummary {
        let mut f_a: usize = 0;
        let mut r_a: usize = 0;
        let mut f_c: usize = 0;
        let mut r_c: usize = 0;
        let mut f_g: usize = 0;
        let mut r_g: usize = 0;
        let mut f_t: usize = 0;
        let mut r_t: usize = 0;
        let mut f_n: usize = 0;
        let mut r_n: usize = 0;

        let mut indel_dict: HashMap<String, usize> = HashMap::new();

        let ref_char: char = ref_char.to_ascii_uppercase();

        let mut i: usize = 0;
        let s: Vec<char> = s.chars().collect();
        while i < s.len() {
            match s[i] {
                // Indicates start of read
                // Char after ^ indicates MAPQ score
                // Skip this char and next
                '^' => i += 2,
                // Count reference bases if ref is specified
                cc if cc == '.' => {
                    match ref_char {
                        'A' => {
                            f_a += 1;
                            i += 1;
                        },
                        'C' => {
                            f_c += 1;
                            i += 1;
                        },
                        'G' => {
                            f_g += 1;
                            i += 1;
                        },
                        'T' => {
                            f_t += 1;
                            i += 1;
                        },
                        'N' => {
                            f_n += 1;
                            i += 1;
                        },
                        _ => (),
                    }
                    i += 1;
                },
                ',' => {
                    match ref_char {
                        'A' => {
                            r_a += 1;
                            i += 1;
                        },
                        'C' => {
                            r_c += 1;
                            i += 1;
                        },
                        'G' => {
                            r_g += 1;
                            i += 1;
                        },
                        'T' => {
                            r_t += 1;
                            i += 1;
                        },
                        'N' => {
                            r_n += 1;
                            i += 1;
                        },
                        _ => (),
                    }
                    i += 1;
                },
                // Count bases
                'A' => {
                    f_a += 1;
                    i += 1;
                },
                'a' => {
                    r_a += 1;
                    i += 1;
                },
                'C' => {
                    f_c += 1;
                    i += 1;
                },
                'c' => {
                    r_c += 1;
                    i += 1;
                },
                'G' => {
                    f_g += 1;
                    i += 1;
                },
                'g' => {
                    r_g += 1;
                    i += 1;
                },
                'T' => {
                    f_t += 1;
                    i += 1;
                },
                't' => {
                    r_t += 1;
                    i += 1;
                },
                'N' => {
                    f_n += 1;
                    i += 1;
                },
                'n' => {
                    r_n += 1;
                    i += 1;
                },
                // Count indels
                cc if cc == '+' || cc == '-' => {
                    // get indel length
                    let mut indel_len_str: Vec<&char> = Vec::new();
                    for c in s[(i+1)..].iter() {
                        if !c.is_numeric() {
                            break;
                        }
                        indel_len_str.push(c);
                    }
                    // get indel characters
                    let indel_len = indel_len_str.into_iter().collect::<String>().parse::<usize>().unwrap();
                    let indel_ident = s[(i+2)..(i+2+indel_len)].iter().collect::<String>();
                    let indel_ident = format!("{}{}{}", cc, indel_len, indel_ident);

                    // upsert
                    *indel_dict.entry(indel_ident).or_insert(0) += 1;

                    i += 2 + indel_len;
                },
                // skip any other character
                _ => i += 1,
            }
        }

        SiteSummary{
            depth,
            ref_char,
            f_a, f_c, f_g, f_t, f_n,
            r_a, r_c, r_g, r_t, r_n,
            indels: indel_dict,
        }
    }

    pub fn pool(summaries: &Vec<&SiteSummary>) -> SiteSummary {
        if summaries.len() == 0 {
            panic!("Empty summaries vector")
        }
        let ref_char: char = summaries[0].ref_char;
        let mut depth: usize = 0;
        let mut f_a: usize = 0;
        let mut f_c: usize = 0;
        let mut f_g: usize = 0;
        let mut f_t: usize = 0;
        let mut f_n: usize = 0;
        let mut r_a: usize = 0;
        let mut r_c: usize = 0;
        let mut r_g: usize = 0;
        let mut r_t: usize = 0;
        let mut r_n: usize = 0;
        let mut indels: HashMap<String, usize> = HashMap::new();

        summaries.into_iter().for_each(|s| {
            depth += s.depth;
            f_a += s.f_a;
            f_c += s.f_c;
            f_g += s.f_g;
            f_t += s.f_t;
            f_n += s.f_n;
            r_a += s.r_a;
            r_c += s.r_c;
            r_g += s.r_g;
            r_t += s.r_t;
            r_n += s.r_n;
            s.indels.iter().for_each(|(k, v)| {
                *indels.entry(k.to_owned()).or_insert(0) += v;
            });
        });

        SiteSummary{
            depth,
            ref_char,
            f_a, f_c, f_g, f_t, f_n,
            r_a, r_c, r_g, r_t, r_n,
            indels,
        }
    }

    pub fn to_hashmaps(&self) -> (HashMap<String, usize>, HashMap<String, usize>) {
        let subs_dict: HashMap<String, usize> = [
            ("A".to_owned(), self.f_a),
            ("C".to_owned(), self.f_c),
            ("G".to_owned(), self.f_g),
            ("T".to_owned(), self.f_t),
            ("N".to_owned(), self.f_n),
            ("a".to_owned(), self.r_a),
            ("c".to_owned(), self.r_c),
            ("g".to_owned(), self.r_g),
            ("t".to_owned(), self.r_t),
            ("n".to_owned(), self.r_n),
        ].iter().cloned().collect();
        let indel_dict: HashMap<String, usize> = self.indels.clone();

        (subs_dict, indel_dict)
    }

    pub fn count_of_any(&self, base: char) -> usize {
        let base: char = base.to_ascii_uppercase();
        match base {
            'A' => self.f_a + self.r_a,
            'C' => self.f_c + self.r_c,
            'G' => self.f_g + self.r_g,
            'T' => self.f_t + self.r_t,
            'N' => self.f_n + self.r_n,
            _ => 0,
        }
    }

    pub fn count_of(&self, base: char) -> usize {
        match base {
            'A' => self.f_a,
            'C' => self.f_c,
            'G' => self.f_g,
            'T' => self.f_t,
            'N' => self.f_n,
            'a' => self.r_a,
            'c' => self.r_c,
            'g' => self.r_g,
            't' => self.r_t,
            'n' => self.r_n,
            _ => 0,
        }
    }

    pub fn freq_of_any(&self, base: char) -> f64 {
        (self.count_of_any(base) as f64) / (self.depth as f64)
    }

    pub fn freq_of(&self, base: char) -> f64 {
        (self.count_of(base) as f64) / (self.depth as f64)
    }

}

#[derive(Debug)]
pub struct SNV {
    major_allele: char,
    major_freq: f64,
    minor_allele: char,
    minor_freq: f64,
}

impl SNV {
    pub fn basic_call(s: &SiteSummary) -> Option<SNV> {
        let sorted_bases: Vec<(char, usize)> = SNV::sort_alleles(s);
        if sorted_bases[0].1 == s.depth {
            None
        } else {
            Some(SNV{
                major_allele: sorted_bases[0].0,
                major_freq: (sorted_bases[0].1 as f64) / (s.depth as f64),
                minor_allele: sorted_bases[1].0,
                minor_freq: (sorted_bases[1].1 as f64) / (s.depth as f64),
            })
        }
    }

    pub fn threshold_call(s: &SiteSummary, threshold_cnt: Option<usize>, threshold_pct: Option<f64>) -> Option<SNV> {
        let sorted_bases: Vec<(char, usize)> = SNV::sort_alleles(s);
        if sorted_bases[0].1 == s.depth {
            None
        } else {
            let minor_cnt = sorted_bases[1].1;
            let minor_freq = (minor_cnt as f64) / (s.depth as f64);
            let threshold_cnt = threshold_cnt.unwrap_or(0);
            let threshold_pct = threshold_pct.unwrap_or(0.);
            if minor_cnt > threshold_cnt {
                if minor_freq > threshold_pct {
                    return Some(SNV{
                        major_allele: sorted_bases[0].0,
                        major_freq: (sorted_bases[0].1 as f64) / (s.depth as f64),
                        minor_allele: sorted_bases[1].0,
                        minor_freq: (sorted_bases[1].1 as f64) / (s.depth as f64),
                    });
                } else {
                    None
                }
            } else {
                None
            }
        }
    }

    pub fn similar_to(&self, other: SNV) -> bool {
        if (self.major_allele == other.major_allele) &&
           (self.minor_allele == other.minor_allele) {
            true
        } else {
            false
        }
    }

    fn sort_alleles(s: &SiteSummary) -> Vec<(char, usize)> {
        let cnt_a: usize = s.count_of_any('A');
        let cnt_c: usize = s.count_of_any('C');
        let cnt_g: usize = s.count_of_any('G');
        let cnt_t: usize = s.count_of_any('T');

        // sort bases from largest freq
        // top two will be major and minor respectively
        let mut sorted_bases: Vec<(char, usize)> = vec![
            ('A', cnt_a),
            ('C', cnt_c),
            ('G', cnt_g),
            ('T', cnt_t),
        ];
        sorted_bases.sort_unstable_by(|(_, v1), (_, v2)| v2.cmp(v1));
        sorted_bases
    }
}

#[derive(Debug)]
pub enum SNVCallType {
    True,
    ProbablyTrue,
    Undecided,
    ProbablyFalse,
    False,
    InsufficientData,
}


#[derive(Debug)]
pub struct PairedSNVCall {
    // target alleles and frequency
    target_major_allele: char,
    target_major_freq: f64,
    target_minor_allele: Option<char>,
    target_minor_freq: f64,
    // control alleles and frequency
    control_major_allele: char,
    control_major_freq: f64,
    control_minor_allele: Option<char>,
    control_minor_freq: f64,
    call: SNVCallType,
}

impl PairedSNVCall {
    pub fn basic_call(target: &SiteSummary, control: &SiteSummary) -> Option<PairedSNVCall> {
        let sorted_target: Vec<(char, usize)> = SNV::sort_alleles(target);
        let sorted_control: Vec<(char, usize)> = SNV::sort_alleles(target);

        // Target and control have the same major allele
        if sorted_target[0].0 == sorted_control[0].0 {
            // Target has no variant
            if sorted_target[1].1 == 0 {
                // Control also has no variants
                if sorted_control[1].1 == 0 {
                    return None
                }
                // Target has no variants but control has variants
                // Since this is in the context of the target, it is false
                return Some(PairedSNVCall{
                    target_major_allele: sorted_target[0].0,
                    target_major_freq: 1.0,
                    target_minor_allele: None,
                    target_minor_freq: 0.0,
                    // control alleles and frequency
                    control_major_allele: sorted_control[0].0,
                    control_major_freq: (sorted_control[0].1 as f64) / (control.depth as f64),
                    control_minor_allele: Some(sorted_control[1].0),
                    control_minor_freq: (sorted_control[1].1 as f64) / (control.depth as f64),
                    call: SNVCallType::False,
                })
            }
            // Target has a variant
            if sorted_control[1].1 == 0 {
                // Control has 0 variants so target SNV is likely true
                return Some(PairedSNVCall{
                    target_major_allele: sorted_target[0].0,
                    target_major_freq: (sorted_target[0].1 as f64) / (target.depth as f64),
                    target_minor_allele: Some(sorted_target[1].0),
                    target_minor_freq: (sorted_target[1].1 as f64) / (target.depth as f64),
                    // control alleles and frequency
                    control_major_allele: sorted_control[0].0,
                    control_major_freq: 1.0,
                    control_minor_allele: None,
                    control_minor_freq: 0.0,
                    call: SNVCallType::True,
                })
            } else if sorted_target[1].0 == sorted_control[1].0 {
                // Control contains at least one variant and 
                // it is the same variant then target SNV is probably false
                return Some(PairedSNVCall{
                    target_major_allele: sorted_target[0].0,
                    target_major_freq: (sorted_target[0].1 as f64) / (target.depth as f64),
                    target_minor_allele: Some(sorted_target[1].0),
                    target_minor_freq: (sorted_target[1].1 as f64) / (target.depth as f64),
                    // control alleles and frequency
                    control_major_allele: sorted_control[0].0,
                    control_major_freq: (sorted_control[0].1 as f64) / (control.depth as f64),
                    control_minor_allele: Some(sorted_control[1].0),
                    control_minor_freq: (sorted_control[1].1 as f64) / (control.depth as f64),
                    call: SNVCallType::ProbablyFalse,
                })
            }
            // Control contains a different variant which indicates probably multiple alleles
            // or some artifact of sequencing or mapping
            // then the target SNV is considered false
            return Some(PairedSNVCall{
                target_major_allele: sorted_target[0].0,
                target_major_freq: (sorted_target[0].1 as f64) / (target.depth as f64),
                target_minor_allele: Some(sorted_target[1].0),
                target_minor_freq: (sorted_target[1].1 as f64) / (target.depth as f64),
                // control alleles and frequency
                control_major_allele: sorted_control[0].0,
                control_major_freq: (sorted_control[0].1 as f64) / (control.depth as f64),
                control_minor_allele: Some(sorted_control[1].0),
                control_minor_freq: (sorted_control[1].1 as f64) / (control.depth as f64),
                call: SNVCallType::False,
            })  
        }
        // Target and control have different major alleles
        // Target has no variant
        if sorted_target[1].1 == 0 {
            // Control also has no variants
            if sorted_control[1].1 == 0 {
                return None
            } 
            // Target has no variants but control has variants
            // Since this is in the context of the target, it is false
            return Some(PairedSNVCall{
                target_major_allele: sorted_target[0].0,
                target_major_freq: 1.0,
                target_minor_allele: None,
                target_minor_freq: 0.0,
                // control alleles and frequency
                control_major_allele: sorted_control[0].0,
                control_major_freq: (sorted_control[0].1 as f64) / (control.depth as f64),
                control_minor_allele: Some(sorted_control[1].0),
                control_minor_freq: (sorted_control[1].1 as f64) / (control.depth as f64),
                call: SNVCallType::False,
            }) 
        }
        // Target and control have different major alleles
        // Target and control have minor alleles
        return Some(PairedSNVCall{
            target_major_allele: sorted_target[0].0,
            target_major_freq: (sorted_target[0].1 as f64) / (target.depth as f64),
            target_minor_allele: Some(sorted_target[1].0),
            target_minor_freq: (sorted_target[1].1 as f64) / (target.depth as f64),
            // control alleles and frequency
            control_major_allele: sorted_control[0].0,
            control_major_freq: (sorted_control[0].1 as f64) / (control.depth as f64),
            control_minor_allele: Some(sorted_control[1].0),
            control_minor_freq: (sorted_control[1].1 as f64) / (control.depth as f64),
            call: SNVCallType::False,
        })
    }

    pub fn threshold_call(target: &SiteSummary, control: &SiteSummary,
            threshold_target_cnt: Option<usize>, threshold_target_freq: Option<f64>,
            threshold_control_cnt: Option<usize>, threshold_control_freq: Option<f64>) -> Option<PairedSNVCall> {

        let threshold_target_cnt = threshold_target_cnt.unwrap_or(0);
        let threshold_target_freq = threshold_target_freq.unwrap_or(0.);
        let threshold_control_cnt = threshold_control_cnt.unwrap_or(0);
        let threshold_control_freq = threshold_control_freq.unwrap_or(0.);

        let sorted_target: Vec<(char, usize)> = SNV::sort_alleles(target);
        let sorted_control: Vec<(char, usize)> = SNV::sort_alleles(target);

        let target_maj_allele = sorted_target[0].0;
        let target_maj_cnt = sorted_target[0].1;
        let target_maj_freq = (sorted_target[0].1 as f64) / (target.depth as f64);
        let target_min_allele = sorted_target[1].0;
        let target_min_cnt = sorted_target[1].1;
        let target_min_freq = (sorted_target[1].1 as f64) / (target.depth as f64);

        let control_maj_allele = sorted_control[0].0;
        let control_maj_cnt = sorted_control[0].1;
        let control_maj_freq = (sorted_control[0].1 as f64) / (control.depth as f64);
        let control_min_allele = sorted_control[1].0;
        let control_min_cnt = sorted_control[1].1;
        let control_min_freq = (sorted_control[1].1 as f64) / (control.depth as f64);

        // Target and control have the same major allele
        if target_maj_allele == control_maj_allele {
            // Target has effectively no variant
            // Variant count is 0 or less than threshold count
            // Variant freq is 0 or less than threshold freq
            if (target_min_cnt <= threshold_target_cnt) || (target_min_freq <= threshold_target_freq) {
                // Control also has effectively no variants
                // Variant count is 0 or less than threshold count
                // Variant freq is 0 or less than threshold freq
                if (control_min_cnt <= threshold_control_cnt) || (control_min_freq <= threshold_control_freq) {
                    return None
                }
                // Target has no variants but control has variants greater than the thresholds
                // Since this is in the context of the target, it is false
                return Some(PairedSNVCall{
                    target_major_allele: sorted_target[0].0,
                    target_major_freq: 1.0,
                    target_minor_allele: None,
                    target_minor_freq: 0.0,
                    // control alleles and frequency
                    control_major_allele: sorted_control[0].0,
                    control_major_freq: (sorted_control[0].1 as f64) / (control.depth as f64),
                    control_minor_allele: Some(sorted_control[1].0),
                    control_minor_freq: (sorted_control[1].1 as f64) / (control.depth as f64),
                    call: SNVCallType::False,
                })
            }
            // Target has a variant
            // Control has effectively no variants so target SNV is likely true
            if (control_min_cnt <= threshold_control_cnt) || (control_min_freq <= threshold_control_freq) {
                return Some(PairedSNVCall{
                    target_major_allele: sorted_target[0].0,
                    target_major_freq: (sorted_target[0].1 as f64) / (target.depth as f64),
                    target_minor_allele: Some(sorted_target[1].0),
                    target_minor_freq: (sorted_target[1].1 as f64) / (target.depth as f64),
                    // control alleles and frequency
                    control_major_allele: sorted_control[0].0,
                    control_major_freq: 1.0,
                    control_minor_allele: None,
                    control_minor_freq: 0.0,
                    call: SNVCallType::True,
                })
            } else if target_min_allele == control_min_allele {
                // Control contains at least one variant and 
                // it is the same variant then target SNV is probably false
                return Some(PairedSNVCall{
                    target_major_allele: sorted_target[0].0,
                    target_major_freq: (sorted_target[0].1 as f64) / (target.depth as f64),
                    target_minor_allele: Some(sorted_target[1].0),
                    target_minor_freq: (sorted_target[1].1 as f64) / (target.depth as f64),
                    // control alleles and frequency
                    control_major_allele: sorted_control[0].0,
                    control_major_freq: (sorted_control[0].1 as f64) / (control.depth as f64),
                    control_minor_allele: Some(sorted_control[1].0),
                    control_minor_freq: (sorted_control[1].1 as f64) / (control.depth as f64),
                    call: SNVCallType::ProbablyFalse,
                })
            }
            // Control contains a different variant which indicates probably multiple alleles
            // or some artifact of sequencing or mapping
            // then the target SNV is considered false
            return Some(PairedSNVCall{
                target_major_allele: sorted_target[0].0,
                target_major_freq: (sorted_target[0].1 as f64) / (target.depth as f64),
                target_minor_allele: Some(sorted_target[1].0),
                target_minor_freq: (sorted_target[1].1 as f64) / (target.depth as f64),
                // control alleles and frequency
                control_major_allele: sorted_control[0].0,
                control_major_freq: (sorted_control[0].1 as f64) / (control.depth as f64),
                control_minor_allele: Some(sorted_control[1].0),
                control_minor_freq: (sorted_control[1].1 as f64) / (control.depth as f64),
                call: SNVCallType::False,
            })  
        }
        // Target and control have different major alleles
        // Target has effectively no variant
        if (target_min_cnt <= threshold_target_cnt) || (target_min_freq <= threshold_target_freq) {
            // Control also has effectively no variants
            if (control_min_cnt <= threshold_control_cnt) || (control_min_freq <= threshold_control_freq) {
                return None
            } 
            // Target has no variants but control has variants
            // Since this is in the context of the target, it is false
            return Some(PairedSNVCall{
                target_major_allele: sorted_target[0].0,
                target_major_freq: 1.0,
                target_minor_allele: None,
                target_minor_freq: 0.0,
                // control alleles and frequency
                control_major_allele: sorted_control[0].0,
                control_major_freq: (sorted_control[0].1 as f64) / (control.depth as f64),
                control_minor_allele: Some(sorted_control[1].0),
                control_minor_freq: (sorted_control[1].1 as f64) / (control.depth as f64),
                call: SNVCallType::False,
            }) 
        }
        // Target and control have different major alleles
        // Target and control have minor alleles
        return Some(PairedSNVCall{
            target_major_allele: sorted_target[0].0,
            target_major_freq: (sorted_target[0].1 as f64) / (target.depth as f64),
            target_minor_allele: Some(sorted_target[1].0),
            target_minor_freq: (sorted_target[1].1 as f64) / (target.depth as f64),
            // control alleles and frequency
            control_major_allele: sorted_control[0].0,
            control_major_freq: (sorted_control[0].1 as f64) / (control.depth as f64),
            control_minor_allele: Some(sorted_control[1].0),
            control_minor_freq: (sorted_control[1].1 as f64) / (control.depth as f64),
            call: SNVCallType::False,
        })
    }

    pub fn similar_to(&self, other: &PairedSNVCall) -> bool {
        if (self.target_major_allele == other.target_major_allele) &&
           (self.control_major_allele == other.control_major_allele) &&
           (self.target_minor_allele == other.target_minor_allele) &&
           (self.control_minor_allele == other.control_minor_allele) {
            true
        } else {
            false
        }
    }
}

#[derive(Debug)]
pub struct MultiPairedSNVCall{
    // Sample values are stored in the same order as the sample names
    sample: Vec<String>,
    sample_major_allele: char,
    sample_major_freq: Vec<f64>,
    sample_minor_allele: Option<char>,
    sample_minor_freq: Vec<f64>,
    // control alleles and frequency
    control_major_allele: char,
    control_major_freq: f64,
    control_minor_allele: Option<char>,
    control_minor_freq: f64,
    call: SNVCallType,
}

impl MultiPairedSNVCall {
    pub fn pooled_basic_call(samples: &Vec<&SiteSummary>, control: &SiteSummary) -> Option<PairedSNVCall> {
        // pool SiteSummary vec into one and call as PairedSNVCall
        let pooled_target = SiteSummary::pool(samples);
        PairedSNVCall::basic_call(&pooled_target, control)
    }

    pub fn pooled_threshold_call(samples: &Vec<&SiteSummary>, control: &SiteSummary,
            threshold_target_cnt: Option<usize>, threshold_target_freq: Option<f64>,
            threshold_control_cnt: Option<usize>, threshold_control_freq: Option<f64>) -> Option<PairedSNVCall> {
        // pool SiteSummary vec into one and call as PairedSNVCall
        let pooled_target = SiteSummary::pool(samples);
        PairedSNVCall::threshold_call(&pooled_target, control, 
            threshold_target_cnt, threshold_target_freq,
            threshold_control_cnt, threshold_control_freq)
    }
    
    pub fn basic_call(samples: &Vec<&SiteSummary>, control: &SiteSummary) -> Option<MultiPairedSNVCall> {
        // Make pooled basic call first to establish major-minor alleles
        if let Some(pc) = MultiPairedSNVCall::pooled_basic_call(samples, control) {
            // Target and control do not have the same major allele
            if pc.target_major_allele != pc.control_major_allele {
                return Some(MultiPairedSNVCall{
                    // Sample values are stored in the same order as the sample names
                    sample: Vec::new(),
                    sample_major_allele: pc.target_major_allele,
                    sample_major_freq: samples.iter().map(|s| s.freq_of_any(pc.target_major_allele)).collect(),
                    sample_minor_allele: pc.target_minor_allele,
                    sample_minor_freq: match pc.target_minor_allele {
                        Some(min_allele) => samples.iter().map(|s| s.freq_of_any(min_allele)).collect(),
                        None => (0..samples.len()).map(|_| 0.0).collect(),
                    },
                    // control alleles and frequency
                    control_major_allele: pc.control_major_allele,
                    control_major_freq: pc.control_major_freq,
                    control_minor_allele: pc.control_minor_allele,
                    control_minor_freq: pc.control_minor_freq,
                    call: SNVCallType::False,
                })
            }
            // Target and control major alleles are the same and
            // Target has at least 1 variant read 
            if let Some(target_minor_allele) = pc.target_minor_allele {
                // Iterate over
                // let mut sample: Vec<String> = Vec::with_capacity(samples.len());
                let mut sample_major_freq: Vec<f64> = Vec::with_capacity(samples.len());
                let mut sample_minor_freq: Vec<f64> = Vec::with_capacity(samples.len());

                let mut min_allele_present: Vec<usize> = Vec::new();

                samples.iter().enumerate().for_each(|(i, s)| {
                    let maj_freq = s.freq_of_any(pc.target_major_allele);
                    let min_freq = s.freq_of_any(target_minor_allele);
                    if min_freq > 0.0 {
                        min_allele_present.push(i);
                    }
                    sample_major_freq.push(maj_freq);
                    sample_minor_freq.push(min_freq);
                });
                let call_type: SNVCallType = match min_allele_present.len() {
                    1 => SNVCallType::ProbablyFalse,
                    2 => SNVCallType::ProbablyTrue,
                    _ => SNVCallType::True,
                };
                return Some(MultiPairedSNVCall{
                    // Sample values are stored in the same order as the sample names
                    sample: Vec::new(),
                    sample_major_allele: pc.target_major_allele,
                    sample_major_freq,
                    sample_minor_allele: Some(target_minor_allele),
                    sample_minor_freq,
                    // control alleles and frequency
                    control_major_allele: pc.control_major_allele,
                    control_major_freq: pc.control_major_freq,
                    control_minor_allele: pc.control_minor_allele,
                    control_minor_freq: pc.control_minor_freq,
                    call: call_type,
                })
            }
            // Target has no minor but control has a minor allele
            return Some(MultiPairedSNVCall{
                // Sample values are stored in the same order as the sample names
                sample: Vec::new(),
                sample_major_allele: pc.target_major_allele,
                sample_major_freq: (0..samples.len()).map(|_| 1.0).collect(),
                sample_minor_allele: None,
                sample_minor_freq: (0..samples.len()).map(|_| 0.0).collect(),
                // control alleles and frequency
                control_major_allele: pc.control_major_allele,
                control_major_freq: pc.control_major_freq,
                control_minor_allele: pc.control_minor_allele,
                control_minor_freq: pc.control_minor_freq,
                call: SNVCallType::False,
            })
        }
        // No SNV call in pooled samples so there is no SNV, no further calling necessary
        None
    }

    pub fn threshold_call(samples: &Vec<&SiteSummary>, control: &SiteSummary,
            threshold_target_cnt: Option<usize>, threshold_target_freq: Option<f64>,
            threshold_control_cnt: Option<usize>, threshold_control_freq: Option<f64>) -> Option<MultiPairedSNVCall> {
        // Make pooled basic call first to establish major-minor alleles
        if let Some(pc) = MultiPairedSNVCall::pooled_threshold_call(
            samples, control,
            threshold_target_cnt, threshold_target_freq,
            threshold_control_cnt, threshold_control_freq) {
            // Target and control do not have the same major allele
            if pc.target_major_allele != pc.control_major_allele {
                return Some(MultiPairedSNVCall{
                    // Sample values are stored in the same order as the sample names
                    sample: Vec::new(),
                    sample_major_allele: pc.target_major_allele,
                    sample_major_freq: samples.iter().map(|s| s.freq_of_any(pc.target_major_allele)).collect(),
                    sample_minor_allele: pc.target_minor_allele,
                    sample_minor_freq: match pc.target_minor_allele {
                        Some(min_allele) => samples.iter().map(|s| s.freq_of_any(min_allele)).collect(),
                        None => (0..samples.len()).map(|_| 0.0).collect(),
                    },
                    // control alleles and frequency
                    control_major_allele: pc.control_major_allele,
                    control_major_freq: pc.control_major_freq,
                    control_minor_allele: pc.control_minor_allele,
                    control_minor_freq: pc.control_minor_freq,
                    call: SNVCallType::False,
                })
            }
            // Target and control major alleles are the same and
            // Target has at least 1 variant read 
            if let Some(target_minor_allele) = pc.target_minor_allele {
                // Iterate over
                // let mut sample: Vec<String> = Vec::with_capacity(samples.len());
                let mut sample_major_freq: Vec<f64> = Vec::with_capacity(samples.len());
                let mut sample_minor_freq: Vec<f64> = Vec::with_capacity(samples.len());

                let mut min_allele_present: Vec<usize> = Vec::new();
                let threshold_target_freq = threshold_target_freq.unwrap_or(0.0);
                let threshold_target_cnt = threshold_target_cnt.unwrap_or(0);

                samples.iter().enumerate().for_each(|(i, s)| {
                    let maj_freq = s.freq_of_any(pc.target_major_allele);
                    let min_freq = s.freq_of_any(target_minor_allele);
                    let min_cnt = s.count_of_any(target_minor_allele);
                    if (min_freq > threshold_target_freq) && (min_cnt > threshold_target_cnt) {
                        min_allele_present.push(i);
                    }
                    sample_major_freq.push(maj_freq);
                    sample_minor_freq.push(min_freq);
                });
                let call_type: SNVCallType = match min_allele_present.len() {
                    1 => SNVCallType::ProbablyFalse,
                    2 => SNVCallType::ProbablyTrue,
                    _ => SNVCallType::True,
                };
                return Some(MultiPairedSNVCall{
                    // Sample values are stored in the same order as the sample names
                    sample: Vec::new(),
                    sample_major_allele: pc.target_major_allele,
                    sample_major_freq,
                    sample_minor_allele: Some(target_minor_allele),
                    sample_minor_freq,
                    // control alleles and frequency
                    control_major_allele: pc.control_major_allele,
                    control_major_freq: pc.control_major_freq,
                    control_minor_allele: pc.control_minor_allele,
                    control_minor_freq: pc.control_minor_freq,
                    call: call_type,
                })
            }
            // Target has no minor but control has a minor allele
            return Some(MultiPairedSNVCall{
                // Sample values are stored in the same order as the sample names
                sample: Vec::new(),
                sample_major_allele: pc.target_major_allele,
                sample_major_freq: (0..samples.len()).map(|_| 1.0).collect(),
                sample_minor_allele: None,
                sample_minor_freq: (0..samples.len()).map(|_| 0.0).collect(),
                // control alleles and frequency
                control_major_allele: pc.control_major_allele,
                control_major_freq: pc.control_major_freq,
                control_minor_allele: pc.control_minor_allele,
                control_minor_freq: pc.control_minor_freq,
                call: SNVCallType::False,
            })
        }
        // No SNV call in pooled samples so there is no SNV, no further calling necessary
        None
    }

    // pub fn similar_to(&self, other: &MultiPairedSNVCall) -> bool {
    //     false
    // }

}