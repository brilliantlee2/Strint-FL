use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Parser;
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};
use rust_htslib::bam::{self, Read};
use strint_rust::matrices::{add_unique_umi, is_genomic_placeholder, matrix_axes};

#[derive(Debug, Parser)]
#[command(version, about = "Build gene expression matrix from tagged BAM")]
struct Cli {
    bam: PathBuf,

    #[arg(long = "output", default_value = "gene_expression.tsv")]
    output: PathBuf,

    #[arg(long = "verbosity", default_value_t = 2)]
    _verbosity: u8,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let mut bam = bam::Reader::from_path(&cli.bam).with_context(|| format!("open {}", cli.bam.display()))?;
    let mut matrix = HashMap::default();
    for rec in bam.records() {
        let rec = rec?;
        let gene = require_string_tag(&rec, b"GN", "GN")?;
        let cell = require_string_tag(&rec, b"CB", "CB")?;
        let umi = require_string_tag(&rec, b"UB", "UB")?;
        if is_genomic_placeholder(&gene) {
            continue;
        }
        add_unique_umi(&mut matrix, &gene, &cell, &umi);
    }
    write_matrix(&cli.output, "gene", &matrix)
}

fn write_matrix(
    output: &PathBuf,
    first_col: &str,
    matrix: &HashMap<(String, String), HashSet<String>>,
) -> Result<()> {
    let (rows, cols) = matrix_axes(matrix);
    let mut writer = BufWriter::new(File::create(output).with_context(|| format!("create {}", output.display()))?);
    write!(writer, "{first_col}")?;
    for col in &cols {
        write!(writer, "\t{col}")?;
    }
    writeln!(writer)?;
    for row in &rows {
        write!(writer, "{row}")?;
        for col in &cols {
            let count = matrix
                .get(&(row.clone(), col.clone()))
                .map(|umis| umis.len())
                .unwrap_or(0);
            write!(writer, "\t{count}")?;
        }
        writeln!(writer)?;
    }
    Ok(())
}

fn require_string_tag(rec: &bam::Record, tag: &[u8; 2], label: &str) -> Result<String> {
    match rec.aux(tag).with_context(|| format!("missing {label} tag on read {}", String::from_utf8_lossy(rec.qname())))? {
        bam::record::Aux::String(v) => Ok(v.to_string()),
        _ => anyhow::bail!("non-string {label} tag on read {}", String::from_utf8_lossy(rec.qname())),
    }
}
