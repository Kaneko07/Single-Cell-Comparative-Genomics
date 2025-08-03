#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Enhanced Main Controller Script for the Single-Cell Comparative Genomics Pipeline
Now includes integrated 2D sensitivity analysis AND cluster consensus subcommand
Features: Real-time logging, efficient double filtering, 2D parameter exploration, cluster integration, and consensus building
"""
import subprocess
import argparse
import sys
import logging
from pathlib import Path

def setup_logging(log_file):
    """Sets up logging to file and console with real-time output."""
    log_file.parent.mkdir(parents=True, exist_ok=True)
    
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    
    file_handler = logging.FileHandler(log_file, mode='w')
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)
    root_logger.addHandler(file_handler)
    root_logger.addHandler(console_handler)
    
    logging.getLogger("matplotlib").setLevel(logging.WARNING)

def run_command(command, cwd=None, check=True, quiet_console=False):
    """Executes a command with optional console output suppression."""
    command_str = [str(c) for c in command]
    logging.info(f"Executing: {' '.join(command_str)}")
    if cwd:
        logging.debug(f"Working Directory: {cwd}")
    
    try:
        with subprocess.Popen(
            command_str, 
            cwd=cwd, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.STDOUT,
            text=True, 
            bufsize=1,
            universal_newlines=True,
            encoding='utf-8',
            errors='replace'
        ) as proc:
            
            for line in proc.stdout:
                line = line.rstrip()
                if line:
                    logging.info(f"[SUBPROCESS] {line}")
                    if not quiet_console:
                        print(line)
                        sys.stdout.flush()
            
            return_code = proc.wait()
            
            if check and return_code != 0:
                error_msg = f"Command failed with exit code {return_code}: {' '.join(command_str)}"
                logging.error(error_msg)
                print(f"ERROR: {error_msg}")
                raise subprocess.CalledProcessError(return_code, command_str)
            
            logging.info(f"Command completed successfully: {command_str[0]}")
        
    except FileNotFoundError:
        error_msg = f"Command not found: {command[0]}. Is it in your PATH?"
        logging.error(error_msg)
        print(f"ERROR: {error_msg}")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        error_msg = f"Command failed with exit code {e.returncode}"
        logging.error(error_msg)
        print(f"ERROR: {error_msg}")
        sys.exit(1)

def validate_snippy_directory(snippy_dir):
    """Validate that the provided snippy directory contains valid snippy outputs."""
    snippy_dir = Path(snippy_dir)
    if not snippy_dir.is_dir():
        logging.error(f"Snippy path is not a directory: {snippy_dir}")
        return False
    subdirs = [d for d in snippy_dir.iterdir() if d.is_dir()]
    if not subdirs:
        logging.error(f"No subdirectories found in snippy directory: {snippy_dir}")
        return False
    valid_outputs = 0
    required_files = ["snps.vcf", "snps.tab", "snps.log"]
    for subdir in subdirs:
        if all((subdir / fname).exists() for fname in required_files):
            valid_outputs += 1
    if valid_outputs == 0:
        logging.error(f"No valid snippy output directories found in {snippy_dir}")
        return False
    return True

def create_analyze_subparser(subparsers):
    """Create analyze subcommand parser (main pipeline)"""
    analyze_parser = subparsers.add_parser(
        'analyze',
        help='Run main comparative genomics analysis pipeline (generates data for visualization)',
        description="""Enhanced pipeline for single-cell comparative genomics with integrated 2D sensitivity analysis.""",
        formatter_class=argparse.RawTextHelpFormatter
    )
    analyze_parser.add_argument("--reference", required=True, type=Path, help="Path to the reference genome FASTA file.")
    analyze_parser.add_argument("--output_dir", required=True, type=Path, help="Main directory to store all output files.")
    input_group = analyze_parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--reads_dir", type=Path, help="Directory containing paired-end FASTQ files.")
    input_group.add_argument("--snippy_dir", type=Path, help="Directory containing existing snippy output subdirectories.")
    analyze_parser.add_argument("--threads", type=int, default=4, help="Number of parallel processes for snippy.")
    # MODIFIED: Restored the sensitivity analysis arguments
    analyze_parser.add_argument("--analysis_mode", choices=['single', '1d', '2d'], default='2d', help="Analysis mode: single (no sensitivity), 1d (SAG only), 2d (SNP+SAG) [default: 2d]")
    analyze_parser.add_argument("--snp_range", type=int, default=3, help="Range around SNP changepoint for 2D analysis (default: ±3)")
    analyze_parser.add_argument("--sag_range", type=int, default=3, help="Range around SAG changepoint for 2D analysis (default: ±3)")
    analyze_parser.set_defaults(func=run_main_analysis)

def create_consensus_subparser(subparsers):
    """Create consensus subcommand parser"""
    consensus_parser = subparsers.add_parser(
        'consensus',
        help='Build per-cluster consensus sequences and run ANI analysis',
        description="""Builds a representative sequence for each cluster and calculates its ANI value against the reference.
        
This subcommand should be run AFTER 'analyze' and cluster definition.

EXAMPLE:
  python run_pipeline.py consensus \\
    --analysis_dir results/ \\
    --reference ref.fasta \\
    --cluster_dir cluster_defs/ \\
    --output_dir results/consensus_output/
        """,
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    consensus_parser.add_argument("--analysis_dir", required=True, type=Path, help="The main output directory from a previous 'analyze' run.")
    consensus_parser.add_argument("--reference", required=True, type=Path, help="Path to the original reference FASTA file.")
    consensus_parser.add_argument("--cluster_dir", required=True, type=Path, help="Directory containing cluster definition files.")
    consensus_parser.add_argument("--output_dir", required=True, type=Path, help="Output directory for all consensus results.")
    consensus_parser.add_argument("--min_allele_freq", type=float, default=0.51, help="Minimum allele frequency for consensus.")
    consensus_parser.add_argument("--skip_ani", action='store_true', help="Skip OrthoANI calculation.")
    
    consensus_parser.set_defaults(func=run_consensus_analysis)

def create_visualize_subparser(subparsers):
    """Create visualize subcommand parser"""
    visualize_parser = subparsers.add_parser(
        'visualize',
        help='Generate visualizations from a completed sensitivity analysis run',
        description="""Generates plots from the results of the 'analyze' subcommand.
This should be run after you have performed phylogenetic analysis and filled in the
'cluster_analysis_template_true_2d.csv' file with the number of clusters for each MSA.

EXAMPLE:
  # After running 'analyze' and filling the template...
  python run_pipeline.py visualize --analysis_dir ./result07
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    visualize_parser.add_argument("--analysis_dir", required=True, type=Path, help="The main output directory from a previous 'analyze' run.")
    visualize_parser.add_argument("--cluster_file", type=Path, help="Optional: Path to the completed cluster template file if it's not in the default location.")
    visualize_parser.set_defaults(func=run_visualization)

def run_consensus_analysis(args):
    """Run consensus sequence analysis"""
    log_file = args.output_dir / "consensus_pipeline.log"
    setup_logging(log_file)
    
    print("="*80)
    print("🧬 PER-CLUSTER CONSENSUS & ANI ANALYSIS PIPELINE")
    print("="*80)
    
    snp_matrix_path = args.analysis_dir.resolve() / "03_snp_matrix" / "raw_snp_matrix.tsv"
    
    if not all([p.exists() for p in [args.reference, snp_matrix_path, args.cluster_dir]]):
        print("ERROR: One or more required input files/directories not found. Please check paths.")
        sys.exit(1)
    
    scripts_dir = Path(__file__).parent / "scripts"
    consensus_script = scripts_dir / "09_build_cluster_consensus.py"
    
    if not consensus_script.exists():
        print(f"ERROR: Consensus script not found: {consensus_script}")
        sys.exit(1)
    
    cmd = [
        "python", consensus_script,
        "--reference_fasta", args.reference.resolve(),
        "--snp_matrix", snp_matrix_path,
        "--cluster_dir", args.cluster_dir.resolve(),
        "--output_dir", args.output_dir.resolve(),
        "--min_allele_freq", str(args.min_allele_freq)
    ]
    
    if args.skip_ani:
        cmd.append("--skip_ani")
    
    run_command(cmd)
    
    print("\n" + "="*80)
    print("🎉 CONSENSUS WORKFLOW COMPLETED!")
    print("="*80)
    print(f"All results saved in: {args.output_dir}")

def run_main_analysis(args):
    """Run main comparative genomics analysis"""
    
    main_output_dir = args.output_dir.resolve()
    scripts_dir = Path(__file__).parent.resolve() / "scripts"
    setup_logging(main_output_dir / "pipeline.log")
    
    print("="*80)
    print("🧬 MAIN ANALYSIS PIPELINE (DATA GENERATION)")
    print("="*80)
    
    if args.reads_dir:
        snippy_outputs_dir = main_output_dir / "01_snippy_outputs"
    else:
        snippy_outputs_dir = args.snippy_dir.resolve()

    # Step 1: Snippy (Conditional)
    if args.reads_dir:
        print("\n--- Step 1: Snippy ---")
        cmd = ["python", scripts_dir / "01_run_snippy.py", "--reference", args.reference.resolve(), "--reads_dir", args.reads_dir.resolve(), "--output_dir", snippy_outputs_dir, "--threads", str(args.threads)]
        run_command(cmd)
    
    if not validate_snippy_directory(snippy_outputs_dir):
        sys.exit("Snippy outputs validation failed.")

    # Step 2: Core Alignment
    print("\n--- Step 2: Core Alignment ---")
    core_analysis_dir = main_output_dir / "02_core_analysis"
    snp_sites_vcf = core_analysis_dir / "core_snpsites.vcf"
    cmd = ["python", scripts_dir / "02_run_core_alignment.py", "--snippy_dir", snippy_outputs_dir, "--reference", args.reference.resolve(), "--output_dir", core_analysis_dir]
    run_command(cmd)

    # Step 3: VCF to Matrix
    print("\n--- Step 3: VCF to Matrix ---")
    matrix_dir = main_output_dir / "03_snp_matrix"
    raw_matrix_file = matrix_dir / "raw_snp_matrix.tsv"
    snp_freq_file = matrix_dir / "snp_frequency.tsv"
    cmd = ["python", scripts_dir / "03_vcf_to_matrix.py", "--input_vcf", snp_sites_vcf, "--output_matrix", raw_matrix_file, "--output_frequency", snp_freq_file]
    run_command(cmd)

    # --- Step 4 & 5: Changepoint Analysis ---
    print("\n--- Step 4 & 5: Changepoint Analysis ---")
    snp_analysis_dir = main_output_dir / "04_changepoint_snp"
    snp_threshold_file = snp_analysis_dir / "snp_threshold.txt"
    cmd = ["Rscript", scripts_dir / "04_run_changepoint.R", snp_freq_file, snp_threshold_file, snp_analysis_dir / "changepoint_snp_plot.png", "Number of Samples Sharing SNP", "Frequency of SNPs"]
    run_command(cmd)

    n_counts_dir = main_output_dir / "05_n_counts"
    sag_n_counts_freq_file = n_counts_dir / "sag_n_counts_frequency.tsv"
    cmd = ["python", scripts_dir / "05_calculate_n_counts.py", "--input_matrix", raw_matrix_file, "--snp_threshold_file", snp_threshold_file, "--output_frequency", sag_n_counts_freq_file]
    run_command(cmd)
    
    sag_analysis_dir = main_output_dir / "06_changepoint_sag"
    sag_threshold_file = sag_analysis_dir / "sag_threshold.txt"
    cmd = ["Rscript", scripts_dir / "04_run_changepoint.R", sag_n_counts_freq_file, sag_threshold_file, sag_analysis_dir / "changepoint_sag_plot.png", "Number of 'N's per SAG", "Frequency of SAGs"]
    run_command(cmd)

    # --- Step 6: Sensitivity Analysis or Final MSA ---
    print("\n--- Step 6: Final MSA Generation / Sensitivity Analysis ---")
    if args.analysis_mode == 'single':
        final_msa_dir = main_output_dir / "07_final_msa"
        cmd = ["python", scripts_dir / "06_generate_final_msa.py", "--input_matrix", raw_matrix_file, "--snp_threshold_file", snp_threshold_file, "--sag_threshold_file", sag_threshold_file, "--output_prefix", final_msa_dir / "final_filtered_msa"]
        run_command(cmd)
    elif args.analysis_mode == '2d':
        sensitivity_2d_dir = main_output_dir / "07_sensitivity_analysis_true_2d"
        cmd = ["python", scripts_dir / "07_true_2d_sensitivity.py", "--input_matrix", raw_matrix_file, "--snp_threshold_file", snp_threshold_file, "--output_dir", sensitivity_2d_dir, "--snp_range", str(args.snp_range), "--sag_range", str(args.sag_range)]
        run_command(cmd)
        
    print("\n" + "="*80)
    print("🎉 MAIN ANALYSIS (DATA GENERATION) FINISHED SUCCESSFULLY!")
    print("="*80)
    print("\nNEXT STEPS:")
    print("1. Perform phylogenetic analysis on the MSA files found in the '07_sensitivity_analysis_true_2d/phylogenetic_msa_true_2d' directory.")
    print("2. Fill in the 'cluster_count' column in '07_sensitivity_analysis_true_2d/cluster_analysis_template_true_2d.csv'.")
    print(f"3. Run the visualization command: python run_pipeline.py visualize --analysis_dir {args.output_dir}")

def run_visualization(args):
    """Generates plots from a completed analysis run."""
    main_output_dir = args.analysis_dir.resolve()
    scripts_dir = Path(__file__).parent.resolve() / "scripts"
    
    # 08_visualizationディレクトリを作成
    visualization_dir = main_output_dir / "08_visualization"
    visualization_dir.mkdir(parents=True, exist_ok=True)
    log_file = visualization_dir / "visualization.log"
    setup_logging(log_file)

    print("="*80)
    print("📊 VISUALIZATION PIPELINE")
    print("="*80)
    print(f"Analysis directory: {main_output_dir}")
    print(f"Visualization results will be saved in: {visualization_dir}")

    # Define paths based on the analysis directory
    sensitivity_2d_dir = main_output_dir / "07_sensitivity_analysis_true_2d"
    sensitivity_2d_csv = sensitivity_2d_dir / "true_2d_sensitivity_results.csv"
    snp_threshold_file = main_output_dir / "04_changepoint_snp" / "snp_threshold.txt"
    sag_threshold_file = main_output_dir / "06_changepoint_sag" / "sag_threshold.txt"
    
    # Use user-provided cluster file or find the default template
    cluster_file_to_use = args.cluster_file
    if not cluster_file_to_use:
        cluster_file_to_use = sensitivity_2d_dir / "cluster_analysis_template_true_2d.csv"

    # Check if all necessary files exist
    required_files = [sensitivity_2d_csv, snp_threshold_file, sag_threshold_file, cluster_file_to_use]
    if not all(f.exists() for f in required_files):
        print("ERROR: Not all required files for visualization were found in the analysis directory.")
        print("Please ensure the 'analyze' command completed successfully and the cluster file exists.")
        for f in required_files:
            if not f.exists():
                print(f"  Missing: {f}")
        sys.exit(1)

    print(f"\nInput files:")
    print(f"  Sensitivity data: {sensitivity_2d_csv}")
    print(f"  Cluster data: {cluster_file_to_use}")
    print(f"  SNP threshold: {snp_threshold_file}")
    print(f"  SAG threshold: {sag_threshold_file}")

    with open(snp_threshold_file, 'r') as f: snp_changepoint = f.read().strip()
    with open(sag_threshold_file, 'r') as f: sag_changepoint = f.read().strip()
    
    print(f"\nChangepoints:")
    print(f"  SNP changepoint: {snp_changepoint}")
    print(f"  SAG changepoint: {sag_changepoint}")
    
    cmd = [
        "Rscript", scripts_dir / "08_visualize_2d_sensitivity.R",
        sensitivity_2d_csv,
        visualization_dir,  # 08_visualizationディレクトリを指定
        snp_changepoint,
        sag_changepoint,
        cluster_file_to_use
    ]
    run_command(cmd)

    print("\n" + "="*80)
    print("🎉 VISUALIZATION FINISHED SUCCESSFULLY!")
    print("="*80)
    print(f"Plots saved in: {visualization_dir}")
    print("\nGenerated files:")
    if visualization_dir.exists():
        plot_files = []
        other_files = []
        for file_path in sorted(visualization_dir.glob("*")):
            if file_path.is_file():
                if file_path.suffix.lower() in ['.png', '.jpg', '.jpeg', '.pdf']:
                    plot_files.append(file_path.name)
                else:
                    other_files.append(file_path.name)
        
        if plot_files:
            print("  Plot files:")
            for file_name in plot_files:
                print(f"    - {file_name}")
        
        if other_files:
            print("  Other files:")
            for file_name in other_files:
                print(f"    - {file_name}")
    
    print(f"\nTo view the plots, navigate to: {visualization_dir}")

def main():
    parser = argparse.ArgumentParser(
        description="Enhanced Single-Cell Comparative Genomics Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    subparsers = parser.add_subparsers(dest='command', help='Available commands', required=True)
    
    create_analyze_subparser(subparsers)
    create_consensus_subparser(subparsers)
    create_visualize_subparser(subparsers)
    
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()