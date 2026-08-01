[![English](https://img.shields.io/badge/Language-English-2563eb)](README.md)
[![中文](https://img.shields.io/badge/语言-中文-0f766e)](README_zh-CN.md)

# Strint Rust 0.0.3.5

Strint Rust is an analysis workflow for full-length single-cell RNA sequencing. It combines a Rust implementation of barcode correction, read tagging, gene/transcript assignment, UMI clustering, matrix generation, QC, saturation analysis, and an offline HTML report with Python workflow helpers.

The public release name is `0.0.3.5`. Cargo records it as `0.0.3+5`, because Cargo follows Semantic Versioning and does not accept four numeric version components.

## Repository contents

The GitHub repository must include these files and directories:

```text
Cargo.toml                 Rust package definition
Cargo.lock                 Reproducible Rust dependency versions
src/                       Rust source and command-line binaries
scripts/                   Python fallbacks, QC, plotting, and report code
vendor/rust-htslib/        Patched rust-htslib required by Cargo.toml
run_all.sh                 Single-species workflow
run_all_mixed_species.sh   Mixed-species workflow
main.py                    Python fallback for the barcode workflow
args_parser.py
utils.py
environment.yml            Recommended Conda environment
requirements.txt           Python-only dependency list
tests/                     Release and report tests
```

Do not upload FASTQ/BAM files, reference genomes, analysis output directories, `target/`, `report_new*/`, or `vendor.zip`. These are covered by `.gitignore`.

## System requirements

- Linux x86_64 is recommended. The workflow and examples use Bash and standard HPC/Linux tools.
- Conda, Miniforge, Mambaforge, or Micromamba.
- Internet access is required for the first Conda environment creation and first Cargo build.
- Memory and CPU requirements depend on input size. For large full-length datasets, request enough RAM for alignment and barcode/cell assignment and start with 16-32 threads unless the node has been benchmarked.
- Glycine is an external optional dependency. It is required unless `--skip-glycine` is used.

## Installation

### 1. Download the source

Using Git:

```bash
git clone https://github.com/brilliantlee2/Strint-FL.git
cd Strint-FL
```

Alternatively, download the GitHub source ZIP, extract it, and enter the extracted directory.

Confirm that the patched dependency is present:

```bash
test -f vendor/rust-htslib/Cargo.toml
```

### 2. Create the recommended environment

The Conda `rust` package includes both `rustc` and `cargo`; do not add a separate `cargo` package.

```bash
conda env create -f environment.yml
conda activate strint
```

If Conda is configured to use `libmamba` but that solver plugin is not
available, create the environment once with the classic solver:

```bash
CONDA_SOLVER=classic CONDA_CHANNEL_PRIORITY=strict \
  conda env create -f environment.yml
```

Verify the tools:

```bash
python --version
rustc --version
cargo --version
samtools --version
minimap2 --version
bedtools --version
```

The environment file installs `pysam`, `bioframe`, and `fast-edit-distance`
through pip after resolving the Conda packages. If that pip phase fails on a
restricted cluster, activate the environment and retry it explicitly:

```bash
python -m pip install \
  "pysam>=0.24,<0.25" \
  "bioframe>=0.8,<0.9" \
  "fast-edit-distance>=1.2,<1.3"
```

For an existing CPython 3.11 environment, Python dependencies can
instead be installed with:

```bash
python -m pip install -r requirements.txt
```

RNA Cluster Analysis adds these direct Python dependencies:

```text
scanpy>=1.11,<1.12
anndata>=0.12,<0.13
numpy>=2.4,<2.5
scipy>=1.17,<1.18
numba>=0.66,<0.67
igraph>=1.0,<1.1
leidenalg>=0.12,<0.13
```

With Conda, the `igraph` package is named `python-igraph`:

```bash
conda install -c conda-forge \
  "scanpy>=1.11,<1.12" \
  "anndata>=0.12,<0.13" \
  "numpy>=2.4,<2.5" \
  "scipy>=1.17,<1.18" \
  "numba>=0.66,<0.67" \
  "python-igraph>=1.0,<1.1" \
  "leidenalg>=0.12,<0.13"
```

Verify the clustering dependencies with:

```bash
python -c "import scanpy, anndata, numba, igraph, leidenalg; print(scanpy.__version__)"
```

CPython 3.11 is the supported and tested target. The environment pins Scanpy
to the 1.11 release line, AnnData to 0.12, NumPy to 2.4, and Numba to 0.66.
These ranges reproduce the dependency combination validated for this release
while allowing compatible patch updates.

### 3. Build the Rust binaries

From the repository root:

```bash
cargo build --release
```

For a strictly lock-file-based build:

```bash
cargo build --release --locked
```

The executables are created in `target/release/`, including `strint-rust`, `prepare_read_tags`, `assign_genes`, `cluster_umis_allbam`, `gene_expression`, `rna_qc_metrics`, and other workflow binaries.

Check the installation:

```bash
./target/release/strint-rust --help
bash run_all.sh -h
```

## Glycine configuration

When Glycine should run, install or build Glycine separately and either add its executable directory to `PATH` or pass it explicitly:

```bash
--glycine-bin-dir /path/to/Glycine/target/release
```

The directory must contain an executable named `glycine`. If full-length reads have already been generated, use `--skip-glycine` and provide `--full-length-fastq` instead.

## Reference preparation

`--ref-dir` expects the following files:

```text
reference/
├── genome.fa
├── genes.gtf
├── genes.bed
└── chrom_sizes.tsv
```

Generate the FASTA index, splice-junction BED, and chromosome sizes with:

```bash
samtools faidx genome.fa
paftools.js gff2bed -j genes.gtf > genes.bed
cut -f1,2 genome.fa.fai | sort -V > chrom_sizes.tsv
```

`paftools.js` is distributed with minimap2. If it is not directly on `PATH`, check the minimap2 installation directory.

## Run examples

### Raw FASTQ with Glycine

```bash
bash run_all.sh \
  --fastq /data/sample.fastq.gz \
  --barcode-list-10bp /data/BC_1536.txt \
  --ref-dir /data/GRCh38_strint \
  --out-dir ./sample_output \
  --sample-id sample \
  --glycine-bin-dir /software/Glycine/target/release \
  --threads 32 \
  --cluster-threads 8 \
  --top1-alpha 0.1 \
  --max-ed 2
```

### Existing full-length FASTQ without Glycine

```bash
bash run_all.sh \
  --skip-glycine \
  --full-length-fastq /data/sample.full-length-plus-rescued.fq.gz \
  --barcode-list-10bp /data/BC_1536.txt \
  --ref-dir /data/GRCh38_strint \
  --out-dir ./sample_output \
  --sample-id sample \
  --threads 32 \
  --cluster-threads 8 \
  --top1-alpha 0.1 \
  --max-ed 2
```

### Mixed-species workflow

Prepare a merged reference containing `genome.fa`, `genes.gtf`, `genes.bed`, and `chrom_sizes.tsv`, then run:

```bash
bash run_all_mixed_species.sh \
  --skip-glycine \
  --full-length-fastq /data/mixed.full-length-plus-rescued.fq.gz \
  --barcode-list-10bp /data/BC_1536.txt \
  --ref-dir /data/merged_hs_mm_ref \
  --out-dir ./mixed_output \
  --sample-id mixed_sample \
  --threads 32 \
  --cluster-threads 8 \
  --top1-alpha 0.1 \
  --max-ed 2
```

## SGE/qsub example

Create `strint_job.sh`:

```bash
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N strint_sample
#$ -o strint_sample.$JOB_ID.out
#$ -e strint_sample.$JOB_ID.err

set -euo pipefail
source ~/anaconda3/etc/profile.d/conda.sh
conda activate strint

bash /path/to/Strint-FL/run_all.sh \
  --fastq /data/sample.fastq.gz \
  --barcode-list-10bp /data/BC_1536.txt \
  --ref-dir /data/GRCh38_strint \
  --out-dir "$PWD/sample_output" \
  --sample-id sample \
  --glycine-bin-dir /software/Glycine/target/release \
  --threads 32 \
  --cluster-threads 8
```

Submit it with matching scheduler and workflow thread counts:

```bash
qsub -cwd -l vf=128G,p=32 -binding linear:32 -P PROJECT -q QUEUE strint_job.sh
```

Cluster resource names and policies vary; confirm `vf`, `p`, parallel-environment, queue, and binding syntax with the cluster administrator.

## Outputs

The main output directory contains `upstream/`, `alignment/`, `matrix/`, `qc/`, and `logs/`. Key outputs include final read-to-cell assignments, tagged BAM files, gene/isoform expression matrices, RNA QC tables, saturation results, and a self-contained single-cell HTML report.

RNA clustering produces `matrix/<sample-id>.rna_cluster.tsv`. It preserves every
final cell and records UMAP coordinates, Leiden cluster, raw UMI count, and
analysis status. The report renders the same coordinates as an RNA-cluster UMAP
and a raw-UMI UMAP in **Cells > RNA Cluster Analysis**. The first Scanpy run in a
new environment may take longer while Numba initializes its cache.

## Tests

```bash
python -m unittest discover -s tests -v
bash -n run_all.sh run_all_mixed_species.sh
cargo test --locked
```

## License

No open-source license is included yet. Add a license chosen by the repository owner before treating the project as open-source software.
