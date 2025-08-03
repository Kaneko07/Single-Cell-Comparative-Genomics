# scripts/01_run_snippy.py
# This script runs snippy for each sample individually.
import argparse
import subprocess
import sys
import logging
from pathlib import Path
from multiprocessing import Pool

# Setup logger for this script
logger = logging.getLogger(__name__)

def run_command_for_snippy(command):
    """A dedicated command runner for snippy to capture output correctly."""
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        return True, result.stdout, result.stderr
    except subprocess.CalledProcessError as e:
        return False, e.stdout, e.stderr

def run_snippy_for_sample(args_tuple):
    """Wrapper function to run snippy for a single sample."""
    sample_name, r1, r2, ref, out_dir = args_tuple
    
    # Check for final output to allow pipeline re-runs
    if (out_dir / "snps.vcf").exists():
        logger.info(f"Snippy output for '{sample_name}' already exists. Skipping.")
        return sample_name, True, "Exists"

    cmd = ["snippy", "--outdir", str(out_dir), "--ref", str(ref), "--R1", str(r1), "--R2", str(r2), "--force"]
    logger.info(f"Running Snippy for sample: {sample_name}")
    success, stdout, stderr = run_command_for_snippy(cmd)
    
    if not success:
        logger.error(f"Snippy failed for {sample_name}:\n{stderr}")
        return sample_name, False, stderr
    return sample_name, True, "Success"

def main():
    """Finds read pairs and runs snippy for each sample."""
    parser = argparse.ArgumentParser(description="Run Snippy for each sample.")
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--reads_dir", required=True, type=Path)
    parser.add_argument("--output_dir", required=True, type=Path, help="Directory to store individual snippy output folders.")
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    # Basic logger configuration for standalone script execution
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Find read pairs (supports .fq.gz and .fastq.gz)
    r1_files = sorted(list(args.reads_dir.glob("*_1.f*q.gz")))
    tasks = []
    for r1_path in r1_files:
        sample_name = r1_path.name.split('_1.f')[0]
        r2_path_str = f"{sample_name}_2.f{r1_path.name.split('_1.f')[1]}"
        r2_path = r1_path.parent / r2_path_str
        if r2_path.exists():
            tasks.append((sample_name, r1_path, r2_path, args.reference, args.output_dir / sample_name))
        else:
            logger.warning(f"Could not find R2 file for {r1_path.name}. Skipping this sample.")

    if not tasks:
        logger.error("No valid FASTQ pairs found in the reads directory.")
        sys.exit(1)

    logger.info(f"Found {len(tasks)} samples to process.")
    with Pool(processes=args.threads) as pool:
        results = pool.map(run_snippy_for_sample, tasks)
    
    failed_samples = [r[0] for r in results if not r[1]]
    if failed_samples:
        logger.error(f"Snippy failed for the following samples: {', '.join(failed_samples)}")
        sys.exit(1)
    
    logger.info("All snippy runs completed successfully.")

if __name__ == "__main__":
    main()

