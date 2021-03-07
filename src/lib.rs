// use pyo3::prelude::*;
// use pyo3::wrap_pyfunction;

// use std::collections::HashMap;

pub mod pileup;

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
