use std::fmt;
use std::ops::Add;

// use pyo3::prelude::*;
// use pyo3::wrap_pyfunction;

// use std::collections::HashMap;
pub mod pileup;
pub mod filter;
pub mod caller;
pub mod constant;
pub mod util;


#[derive(Debug, PartialEq, Copy, Clone)]
pub enum Base {
    A,
    C,
    T,
    G,
}
impl Base {
    // constructor
    pub fn from_char(c: &char) -> Option<Self> {
        match c {
            'A' => Some(Self::A),
            'C' => Some(Self::C),
            'G' => Some(Self::G),
            'T' => Some(Self::T),
            'a' => Some(Self::A),
            'c' => Some(Self::C),
            'g' => Some(Self::G),
            't' => Some(Self::T),
            _ => None,
        }
    }
    pub fn to_char(&self) -> char {
        match self {
            Self::A => 'A',
            Self::C => 'C',
            Self::G => 'G',
            Self::T => 'T',
        }
    }
}
impl fmt::Display for Base {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", match self {
            Self::A => 'A',
            Self::C => 'C',
            Self::G => 'G',
            Self::T => 'T',
        })
    }
}

#[derive(Debug, PartialEq, Copy, Clone)]
pub enum OrientedBase {
    ForwardA,
    ForwardC,
    ForwardG,
    ForwardT,
    ReverseA,
    ReverseC,
    ReverseG,
    ReverseT,
}
impl OrientedBase {
    // constructor
    pub fn from_char(c: &char) -> Option<Self> {
        match c {
            'A' => Some(Self::ForwardA),
            'C' => Some(Self::ForwardC),
            'G' => Some(Self::ForwardG),
            'T' => Some(Self::ForwardT),
            'a' => Some(Self::ReverseA),
            'c' => Some(Self::ReverseC),
            'g' => Some(Self::ReverseG),
            't' => Some(Self::ReverseT),
            _ => None,
        }
    }
    pub fn to_char(&self) -> char {
        match self {
            Self::ForwardA => 'A',
            Self::ForwardC => 'C',
            Self::ForwardG => 'G',
            Self::ForwardT => 'T',
            Self::ReverseA => 'a',
            Self::ReverseC => 'c',
            Self::ReverseG => 'g',
            Self::ReverseT => 't',
        }
    }
}
impl fmt::Display for OrientedBase {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", match self {
            Self::ForwardA => 'A',
            Self::ForwardC => 'C',
            Self::ForwardG => 'G',
            Self::ForwardT => 'T',
            Self::ReverseA => 'a',
            Self::ReverseC => 'c',
            Self::ReverseG => 'g',
            Self::ReverseT => 't',
        })
    }
}

#[derive(Debug, Copy, Clone)]
pub struct BaseCount {
    f_a: usize,
    f_c: usize,
    f_g: usize,
    f_t: usize,
    r_a: usize,
    r_c: usize,
    r_g: usize,
    r_t: usize,
}
// public methods
impl BaseCount {
    // constructors
    pub fn from_char_vec(v: &Vec<char>) -> Self {
        let (mut f_a, mut f_c, mut f_g, mut f_t) = (0, 0, 0, 0);
        let (mut r_a, mut r_c, mut r_g, mut r_t) = (0, 0, 0, 0);
        for base in v {
            match base {
                'A' => f_a += 1,
                'C' => f_c += 1,
                'G' => f_g += 1,
                'T' => f_t += 1,
                'a' => r_a += 1,
                'c' => r_c += 1,
                'g' => r_g += 1,
                't' => r_t += 1,
                _ => ()
            };
        }
        BaseCount {
            f_a, f_c, f_g, f_t,
            r_a, r_c, r_g, r_t,
        }
    }

    pub fn from_enum_vec(v: &Vec<OrientedBase>) -> Self {
        let (mut f_a, mut f_c, mut f_g, mut f_t) = (0, 0, 0, 0);
        let (mut r_a, mut r_c, mut r_g, mut r_t) = (0, 0, 0, 0);
        for oriented_base in v {
            match oriented_base {
                OrientedBase::ForwardA => f_a += 1,
                OrientedBase::ForwardC => f_c += 1,
                OrientedBase::ForwardG => f_g += 1,
                OrientedBase::ForwardT => f_t += 1,
                OrientedBase::ReverseA => r_a += 1,
                OrientedBase::ReverseC => r_c += 1,
                OrientedBase::ReverseG => r_g += 1,
                OrientedBase::ReverseT => r_t += 1,
            };
        }
        BaseCount {
            f_a, f_c, f_g, f_t,
            r_a, r_c, r_g, r_t,
        }
    }

    pub fn from_str(s: &str, sep: &str) -> Self {
        let v: Vec<usize> = s.split(sep)
            .map(|c| c.parse::<usize>().unwrap())
            .collect();
        if v.len() != 8 {
            panic!(format!("cannot convert string into BaseCount using the given separator '{}'", sep))
        }
        BaseCount {
            f_a: v[0],
            r_a: v[1],
            f_c: v[2],
            r_c: v[3],
            f_g: v[4],
            r_g: v[5],
            f_t: v[6],
            r_t: v[7],
        }
    }

    pub fn empty() -> Self {
        BaseCount {
            f_a: 0,
            r_a: 0,
            f_c: 0,
            r_c: 0,
            f_g: 0,
            r_g: 0,
            f_t: 0,
            r_t: 0,
        }
    }

    // base count getters
    // orientation independent
    pub fn a(&self) -> usize { self.f_a + self.r_a }
    pub fn c(&self) -> usize { self.f_c + self.r_c }
    pub fn g(&self) -> usize { self.f_g + self.r_g }
    pub fn t(&self) -> usize { self.f_t + self.r_t }
    pub fn count_of(&self, base: Base) -> usize {
        let (fb, rb) = match base {
            Base::A => (OrientedBase::ForwardA, OrientedBase::ReverseA),
            Base::C => (OrientedBase::ForwardC, OrientedBase::ReverseC),
            Base::G => (OrientedBase::ForwardG, OrientedBase::ReverseG),
            Base::T => (OrientedBase::ForwardT, OrientedBase::ReverseT),
        };
        self.oriented_count_of(fb) + self.oriented_count_of(rb)
    }

    // orientation dependent
    pub fn forward_a(&self) -> usize { self.f_a }
    pub fn forward_c(&self) -> usize { self.f_c }
    pub fn forward_g(&self) -> usize { self.f_g }
    pub fn forward_t(&self) -> usize { self.f_t }
    pub fn reverse_a(&self) -> usize { self.r_a }
    pub fn reverse_c(&self) -> usize { self.r_c }
    pub fn reverse_g(&self) -> usize { self.r_g }
    pub fn reverse_t(&self) -> usize { self.r_t }
    pub fn oriented_count_of(&self, base: OrientedBase) -> usize {
        match base {
            OrientedBase::ForwardA => self.f_a,
            OrientedBase::ForwardC => self.f_c,
            OrientedBase::ForwardG => self.f_g,
            OrientedBase::ForwardT => self.f_t,
            OrientedBase::ReverseA => self.r_a,
            OrientedBase::ReverseC => self.r_c,
            OrientedBase::ReverseG => self.r_g,
            OrientedBase::ReverseT => self.r_t,
        }
    }

    pub fn forward(&self) -> usize { self.f_a + self.f_c + self.f_g + self.f_t }
    pub fn reverse(&self) -> usize { self.r_a + self.r_c + self.r_g + self.r_t }

    pub fn forward_count_of(&self, base: Base) -> usize {
        let oriented_base = match base {
            Base::A => OrientedBase::ForwardA,
            Base::C => OrientedBase::ForwardC,
            Base::G => OrientedBase::ForwardG,
            Base::T => OrientedBase::ForwardT,
        };
        self.oriented_count_of(oriented_base)
    }

    pub fn reverse_count_of(&self, base: Base) -> usize {
        let oriented_base = match base {
            Base::A => OrientedBase::ReverseA,
            Base::C => OrientedBase::ReverseC,
            Base::G => OrientedBase::ReverseG,
            Base::T => OrientedBase::ReverseT,
        };
        self.oriented_count_of(oriented_base)
    }

    // operations
    pub fn add(&self, other: &BaseCount) -> Self {
        BaseCount {
            f_a: self.f_a + other.f_a,
            f_c: self.f_c + other.f_c,
            f_g: self.f_g + other.f_g,
            f_t: self.f_t + other.f_t,
            r_a: self.r_a + other.r_a,
            r_c: self.r_c + other.r_c,
            r_g: self.r_g + other.r_g,
            r_t: self.r_t + other.r_t,
        }
    }
    pub fn total(&self) -> usize {
        self.f_a + self.f_c + self.f_g + self.f_t + self.r_a + self.r_c + self.r_g + self.r_t
    }

    pub fn cov(&self) -> usize {
        self.total()
    }
}

impl Add<BaseCount> for BaseCount {
    type Output = BaseCount;
    fn add(self, other: BaseCount) -> Self::Output {
        self.add(&other)
    }
}
impl Add<BaseCount> for &BaseCount {
    type Output = BaseCount;
    fn add(self, other: BaseCount) -> Self::Output {
        self.add(&other)
    }
}
impl Add<&BaseCount> for BaseCount {
    type Output = BaseCount;
    fn add(self, other_ptr: &BaseCount) -> Self::Output {
        self.add(other_ptr)
    }
}
impl Add<&BaseCount> for &BaseCount {
    type Output = BaseCount;
    fn add(self, other_ptr: &BaseCount) -> Self::Output {
        self.add(other_ptr)
    }
}
impl fmt::Display for BaseCount {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}:{}:{}:{}:{}:{}:{}:{}", 
            self.f_a, self.r_a, 
            self.f_c, self.r_c,
            self.f_g, self.r_g,
            self.f_t, self.r_t,
        )
    }
}

#[derive(Debug, Clone)]
pub struct AlleleSet{
    major: Base,
    minor: Option<Base>,
    others: Vec<Base>,
    base_count: BaseCount,
}
impl AlleleSet {
    // constructor
    pub fn from_base_count(base_count: BaseCount) -> Option<Self> {
        // Get alleles ignoring singletons
        // A base needs to have at least 1 forward and 1 reverse support to be called an allele
        let mut alleles: Vec<(Base, usize)> = vec![
                (Base::A, base_count.forward_a(), base_count.reverse_a()),
                (Base::C, base_count.forward_c(), base_count.reverse_c()),
                (Base::G, base_count.forward_g(), base_count.reverse_g()),
                (Base::T, base_count.forward_t(), base_count.reverse_t()),
            ]
            .into_iter()
            .filter_map(|(base, f, r)| match (f, r) {
                (0, 0) => None,
                (1, 0) => None,
                (0, 1) => None,
                _ => Some((base, f + r))
            })
            .collect();
        if alleles.len() == 0 { return None }
        alleles.sort_by(|a, b| b.1.cmp(&a.1));
        // define major and minor alleles
        let major = alleles[0].0;
        let mut minor = None;
        let mut others = Vec::new();
        if alleles.len() > 1 { 
            minor = Some(alleles[1].0);
            others = alleles.split_off(2).into_iter().map(|(b, _)| b ).collect();
        }
        Some(AlleleSet { major, minor, others, base_count })
    }

    // getters
    pub fn major_allele(&self) -> Base { self.major }
    pub fn minor_allele(&self) -> Option<Base> { self.minor }
    pub fn other_alleles(&self) -> Vec<Base> { self.others.clone() }

    pub fn major_count(&self) -> usize { self.base_count.count_of(self.major) }
    pub fn minor_count(&self) -> usize {
        if self.minor.is_none() { return 0 }
        self.base_count.count_of(self.minor.unwrap()) 
    }
    pub fn others_count(&self) -> usize { 
        self.others.iter()
            .map(|b| self.base_count.count_of(*b))
            .sum()
    }
    pub fn base_count(&self) -> &BaseCount { &self.base_count }

    pub fn major_freq(&self) -> f64 { (self.major_count() as f64) / (self.base_count.total() as f64) }
    pub fn minor_freq(&self) -> f64 { (self.minor_count() as f64) / (self.base_count.total() as f64) }

    pub fn num_bases(&self) -> usize { self.base_count.cov() }
    pub fn num_alleles(&self) -> usize { self.len() }
    pub fn len(&self) -> usize {
        let mut count = 1 + self.others.len();
        if let Some(_) = self.minor {
            count += 1;
        }
        count
    }
}

// #[pyfunction]
// fn parse_pileup_record(record: &str) -> PyResult<String> {
//     // Assumes output_MQ is enabled
//     // Format:
//     // chrom pos refchar depth pileup_str base_quals map_quals 

//     Ok("ok".to_owned())
// }

// /// Counts the number of bases and indels in one site
// #[pyfunction]
// fn parse_pileup_str(s: &str, ref_char: char, depth: usize) -> PyResult<(HashMap<String, usize>, HashMap<String, usize>)> {
//     Ok(pileup::SiteSummary::new(s, ref_char, depth).to_hashmaps())
// }

// /// A Python module implemented in Rust.
// #[pymodule]
// fn htsops(_py: Python, m: &PyModule) -> PyResult<()> {
//     m.add_function(wrap_pyfunction!(parse_pileup_str, m)?)?;

//     Ok(())
// }
