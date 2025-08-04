# Single-Cell-Comparative-Genomics

A comprehensive pipeline for comparative genomics analysis of single-cell amplified genomes (SAGs) with statistical optimization and 2D sensitivity analysis.

## Features

- **Automated SNP calling** using Snippy
- **Statistical changepoint analysis** for optimal threshold determination
- **2D sensitivity analysis** for parameter optimization
- **Interactive visualization** with publication-ready plots
- **Consensus sequence building** with ANI analysis

## Quick Start

### 1. Environment Setup
```bash
conda env create -f environment.yml
conda activate single_cell_genomics_pipeline
```

### 2. Run Analysis
```bash
# Main analysis pipeline
python run_pipeline.py analyze \
  --reference reference.fasta \
  --reads_dir fastq_files/ \
  --output_dir results/

# Generate visualizations
python run_pipeline.py visualize \
  --analysis_dir results/

# Build consensus sequences (optional)
python run_pipeline.py consensus \
  --analysis_dir results/ \
  --reference reference.fasta \
  --cluster_dir cluster_definitions/ \
  --output_dir results/consensus/
```

## Pipeline Overview

1. **SNP Calling**: Identifies variants using Snippy
2. **Core Alignment**: Generates core genome alignment
3. **Matrix Conversion**: Creates SNP matrix for analysis
4. **Changepoint Analysis**: Statistically determines optimal thresholds
5. **2D Sensitivity Analysis**: Explores parameter space
6. **Visualization**: Generates comprehensive plots
7. **Consensus Building**: Creates representative sequences per cluster

## Requirements

- Python 3.9+
- R 4.3+
- Conda/Mamba package manager
- Bioinformatics tools (Snippy, snp-sites)

## Input Files

- **Reference genome**: FASTA format
- **Sequencing reads**: Paired-end FASTQ files
- **Cluster definitions**: TSV files (for consensus analysis)

## Output

- SNP matrices and frequency tables
- Statistical changepoint plots
- 2D sensitivity analysis results
- Publication-ready visualizations
- Consensus sequences with ANI values
