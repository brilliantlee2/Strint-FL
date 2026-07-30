[![English](https://img.shields.io/badge/Language-English-2563eb)](README.md)
[![中文](https://img.shields.io/badge/语言-中文-0f766e)](README_zh-CN.md)

# Strint Rust 0.0.3.5

Strint Rust 是一个面向全长单细胞 RNA 测序数据的分析流程。它整合了 Rust 实现的 barcode 矫正、read 标签添加、基因/转录本分配、UMI 聚类、表达矩阵生成，以及 Python 实现的流程辅助、QC、饱和度分析和离线 HTML 报告。

对外发布版本号为 `0.0.3.5`。由于 Cargo 遵循语义化版本规则，不接受四段纯数字版本号，因此 `Cargo.toml` 和 `Cargo.lock` 中对应写为合法的 `0.0.3+5`。

## GitHub 仓库内容

GitHub 仓库需要包含以下文件和目录：

```text
Cargo.toml                 Rust 项目及依赖定义
Cargo.lock                 固定 Rust 依赖版本
src/                       Rust 源代码和各个命令行程序
scripts/                   Python 回退脚本、QC、绘图和报告代码
vendor/rust-htslib/        Cargo.toml 指定的修订版 rust-htslib
run_all.sh                 单物种完整流程
run_all_mixed_species.sh   混合物种流程
main.py                    barcode 流程的 Python 回退入口
args_parser.py
utils.py
environment.yml            推荐的 Conda 环境
requirements.txt           Python 依赖列表
tests/                     发布结构和报告测试
```

不要上传 FASTQ/BAM、参考基因组、分析结果目录、`target/`、`report_new*/` 或 `vendor.zip`；这些内容已写入 `.gitignore`。

## 系统要求

- 推荐 Linux x86_64；流程和示例使用 Bash 及常见 HPC/Linux 工具。
- Conda、Miniforge、Mambaforge 或 Micromamba。
- 第一次创建 Conda 环境及第一次 Cargo 编译需要联网下载依赖。
- CPU 和内存需求取决于数据量。对于大型全长数据，建议先使用 16-32 线程，并为比对和 barcode/cell assignment 申请足够内存。
- Glycine 是独立的可选依赖；不使用 `--skip-glycine` 时必须安装。

## 安装教程

### 1. 下载源代码

使用 Git：

```bash
git clone https://github.com/brilliantlee2/Strint-FL.git
cd Strint-FL
```

也可以在 GitHub 下载 Source code ZIP，解压后进入项目目录。

确认本地修订依赖存在：

```bash
test -f vendor/rust-htslib/Cargo.toml
```

### 2. 创建推荐环境

Conda 中的 `rust` 包已经同时提供 `rustc` 和 `cargo`，不要再单独添加 `cargo` 包。

```bash
conda env create -f environment.yml
conda activate strint-rust
```

检查软件：

```bash
python --version
rustc --version
cargo --version
samtools --version
minimap2 --version
bedtools --version
```

如果受限集群中 Conda 无法解析 `bioframe` 或 `fast-edit-distance`，激活环境后执行：

```bash
python -m pip install bioframe fast-edit-distance
```

如果已经存在 Python 3.11 环境，也可以单独安装 Python 依赖：

```bash
python -m pip install -r requirements.txt
```

### 3. 编译 Rust 程序

在仓库根目录执行：

```bash
cargo build --release
```

如需严格使用 `Cargo.lock` 中的版本：

```bash
cargo build --release --locked
```

可执行文件生成在 `target/release/`，包括 `strint-rust`、`prepare_read_tags`、`assign_genes`、`cluster_umis_allbam`、`gene_expression`、`rna_qc_metrics` 等。

检查安装：

```bash
./target/release/strint-rust --help
bash run_all.sh -h
```

## Glycine 配置

需要运行 Glycine 时，请单独安装或编译 Glycine，并将其可执行文件目录加入 `PATH`，或者在命令中传入：

```bash
--glycine-bin-dir /path/to/Glycine/target/release
```

该目录中必须存在名为 `glycine` 的可执行文件。如果已经获得全长 FASTQ，可使用 `--skip-glycine` 并传入 `--full-length-fastq`。

## 参考基因组准备

`--ref-dir` 目录需要包含：

```text
reference/
├── genome.fa
├── genes.gtf
├── genes.bed
└── chrom_sizes.tsv
```

可使用以下命令生成 FASTA 索引、剪接位点 BED 和染色体长度文件：

```bash
samtools faidx genome.fa
paftools.js gff2bed -j genes.gtf > genes.bed
cut -f1,2 genome.fa.fai | sort -V > chrom_sizes.tsv
```

`paftools.js` 随 minimap2 一起提供。如果它不在 `PATH` 中，请检查 minimap2 的安装目录。

## 运行示例

### 从原始 FASTQ 开始并运行 Glycine

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

### 跳过 Glycine，直接使用全长 FASTQ

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

### 混合物种流程

先准备同时包含 `genome.fa`、`genes.gtf`、`genes.bed` 和 `chrom_sizes.tsv` 的合并参考目录，然后运行：

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

## SGE/qsub 示例

创建 `strint_job.sh`：

```bash
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N strint_sample
#$ -o strint_sample.$JOB_ID.out
#$ -e strint_sample.$JOB_ID.err

set -euo pipefail
source ~/anaconda3/etc/profile.d/conda.sh
conda activate strint-rust

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

投递命令：

```bash
qsub -cwd -l vf=128G,p=32 -binding linear:32 -P PROJECT -q QUEUE strint_job.sh
```

不同集群的资源名称和策略不同，请向管理员确认 `vf`、`p`、并行环境、队列和 CPU binding 的具体写法。

## 主要输出

输出目录包含 `upstream/`、`alignment/`、`matrix/`、`qc/` 和 `logs/`。主要结果包括最终 read-to-cell 分配、带标签 BAM、基因/转录本表达矩阵、RNA QC、饱和度分析和可独立打开的单细胞 HTML 报告。

## 测试

```bash
python -m unittest discover -s tests -v
bash -n run_all.sh run_all_mixed_species.sh
cargo test --locked
```

## 许可证

当前尚未添加开源许可证。在将项目作为开源软件发布前，请由仓库所有者选择并添加合适的许可证。
