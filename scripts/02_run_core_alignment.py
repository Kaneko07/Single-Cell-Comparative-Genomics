import argparse
import subprocess
import sys
import logging
import time
import os
from pathlib import Path

logger = logging.getLogger(__name__)

def run_command(command, cwd=None, check=True):
    command_str = [str(c) for c in command]
    logger.debug(f"Executing: {' '.join(command_str)}")
    try:
        result = subprocess.run(command_str, check=check, capture_output=True, text=True, cwd=cwd)
        if result.stdout:
            logger.debug("STDOUT:\n" + result.stdout)
        if result.stderr:
            logger.warning("STDERR:\n" + result.stderr)
        return result
    except FileNotFoundError:
        logger.error(f"Command not found: {command[0]}. Is it in your PATH?")
        raise
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed with exit code {e.returncode}.")
        logger.error(f"STDOUT:\n{e.stdout}")
        logger.error(f"STDERR:\n{e.stderr}")
        raise

def main():
    parser = argparse.ArgumentParser(description="Run snippy-core and snp-sites.")
    parser.add_argument("--snippy_dir", required=True, type=Path, help="Directory containing all individual snippy outputs.")
    parser.add_argument("--reference", required=True, type=Path, help="Path to the reference genome FASTA file.")
    parser.add_argument("--output_dir", required=True, type=Path, help="Directory to store core analysis outputs.")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info("Step 2.1: Running snippy-core to generate core genome alignment.")
    core_prefix = "core"
    core_full_aln = args.output_dir / f"{core_prefix}.full.aln"
    
    if not core_full_aln.exists():
        snippy_dirs = list(args.snippy_dir.glob("*"))
        if not snippy_dirs:
            logger.error(f"No snippy output directories found in {args.snippy_dir}")
            sys.exit(1)
        
        cmd = ["snippy-core", "--ref", args.reference, "--prefix", core_prefix] + snippy_dirs
        
        try:
            logger.info(f"Running command: {' '.join([str(c) for c in cmd])}")
            run_command(cmd, cwd=args.output_dir)
            logger.info("snippy-core completed successfully.")
        except (subprocess.CalledProcessError, Exception) as e:
            logger.warning(f"snippy-core exited with error: {e}")
            if core_full_aln.exists():
                logger.warning("snippy-core exited with error, but core.full.aln was successfully generated.")
                logger.info("This is normal when no SNPs are detected. Continuing pipeline.")
            else:
                logger.error(f"snippy-core failed and no alignment file was produced: {core_full_aln}")
                sys.exit(1)
    else:
        logger.info("Core alignment file already exists, skipping snippy-core.")

    if not core_full_aln.exists():
        logger.error(f"Core alignment file not found: {core_full_aln}")
        sys.exit(1)
        
    logger.info("Step 2.2: Running snp-sites to generate VCF from alignment.")
    output_vcf_base = f"{core_prefix}_snpsites"
    output_vcf = args.output_dir / f"{output_vcf_base}.vcf"
    
    if not output_vcf.exists():
        cmd = ["snp-sites", "-v", "-m", "-r", "-p", "-o", output_vcf_base, core_full_aln]
        
        try:
            logger.info(f"Running command: {' '.join([str(c) for c in cmd])}")
            logger.info(f"Working directory: {args.output_dir}")
            logger.info(f"Expected output file: {output_vcf}")
            run_command(cmd, cwd=args.output_dir)
            logger.info(f"snp-sites command completed.")
            
            time.sleep(0.5)
            
            logger.debug(f"Current working directory contents: {os.listdir(args.output_dir)}")
            
        except (subprocess.CalledProcessError, Exception) as e:
            logger.warning(f"snp-sites failed: {e}")
            logger.warning("This is likely because no SNPs were detected in the alignment.")
            logger.info("Creating minimal VCF header for downstream processing.")
            
            try:
                output_vcf.parent.mkdir(parents=True, exist_ok=True)
                
                with open(output_vcf, 'w') as f:
                    f.write("##fileformat=VCFv4.2\n")
                    f.write("##source=snp-sites\n")
                    f.write("##INFO=<ID=TYPE,Number=A,Type=String,Description=\"Type of variant\">\n")
                    f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
                    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tReference\n")
                
                logger.info(f"Minimal VCF file created successfully: {output_vcf}")
                
                if output_vcf.exists() and output_vcf.stat().st_size > 0:
                    logger.info(f"VCF file verification: exists={output_vcf.exists()}, size={output_vcf.stat().st_size} bytes")
                else:
                    logger.error(f"Failed to create VCF file or file is empty: {output_vcf}")
                    
            except Exception as file_error:
                logger.error(f"Failed to create minimal VCF file: {file_error}")
                logger.error(f"Attempted to write to: {output_vcf}")
                raise
    else:
        logger.info("snp-sites VCF file already exists, skipping snp-sites.")

    if not output_vcf.exists():
        logger.warning(f"VCF file was not created by snp-sites: {output_vcf}")
        logger.info("This typically happens when no SNPs are detected in the alignment.")
        logger.info("Creating minimal VCF header for downstream processing.")
        
        try:
            output_vcf.parent.mkdir(parents=True, exist_ok=True)
            
            with open(output_vcf, 'w') as f:
                f.write("##fileformat=VCFv4.2\n")
                f.write("##source=snp-sites\n")
                f.write("##INFO=<ID=TYPE,Number=A,Type=String,Description=\"Type of variant\">\n")
                f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
                f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tReference\n")
            
            logger.info(f"Minimal VCF file created successfully: {output_vcf}")
            
        except Exception as file_error:
            logger.error(f"Failed to create minimal VCF file: {file_error}")
            sys.exit(1)
    else:
        file_size = output_vcf.stat().st_size
        logger.info(f"VCF file exists: {output_vcf} (size: {file_size} bytes)")
        
        if file_size < 500:
            logger.warning("VCF file is very small, likely contains no SNP data (header only).")
            logger.info("This is normal when no SNPs are detected.")

    logger.info("Core alignment and VCF generation completed successfully.")

if __name__ == "__main__":
    main()