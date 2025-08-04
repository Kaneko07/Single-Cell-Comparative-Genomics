#!/usr/bin/env python3

import pandas as pd
import argparse
import sys
import logging
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from collections import Counter
import subprocess
import shutil
from io import StringIO

def setup_logging(log_file):
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file, mode='w'),
            logging.StreamHandler()
        ]
    )

def load_cluster_files(cluster_dir):
    clusters = {}
    cluster_files = sorted([p for p in Path(cluster_dir).glob("*") if p.suffix in {".tsv", ".txt"}])
    if not cluster_files:
        raise FileNotFoundError(f"No cluster files (*.tsv, *.txt) found in {cluster_dir}")
    logging.info(f"Found {len(cluster_files)} cluster files.")
    for f in cluster_files:
        try:
            base_name = f.stem
            cluster_name = base_name.split('_')[-1] if '_' in base_name else base_name
            df = pd.read_csv(f, sep='\t', header=None)
            samples = set(s for s in df.iloc[:, 0].dropna().astype(str).tolist() if s.strip())
            if samples:
                if cluster_name in clusters:
                    logging.warning(f"Duplicate cluster name '{cluster_name}' from file {f.name}. Samples will be merged.")
                    clusters[cluster_name].update(samples)
                else:
                    clusters[cluster_name] = samples
                logging.info(f"Loaded {len(samples)} samples for cluster '{cluster_name}' from file {f.name}.")
        except Exception as e:
            logging.error(f"Error reading {f}: {e}")
    
    for name in clusters:
        clusters[name] = list(clusters[name])
    return clusters

def build_per_cluster_consensus(clusters, reference_fasta, snp_matrix, min_allele_freq):
    logging.info("Building per-cluster consensus sequences...")
    ref_rec = SeqIO.read(reference_fasta, "fasta")
    ref_seq = MutableSeq(ref_rec.seq)
    snp_df = pd.read_csv(snp_matrix, sep='\t')
    
    rep_seqs, stats = [], []
    for name, samples in clusters.items():
        logging.info(f"  Processing cluster: {name}")
        valid_samples = [s for s in samples if s in snp_df.columns]
        if not valid_samples:
            logging.warning(f"  No samples from cluster '{name}' found in SNP matrix. Skipping.")
            continue
        
        cluster_seq = ref_seq[:]
        snps_applied = 0
        for _, row in snp_df.iterrows():
            counts = Counter(str(g) for g in row[valid_samples].values)
            n_ref, n_alt = counts.get('0', 0), counts.get('1', 0)
            n_valid = n_ref + n_alt
            if n_valid > 0 and (n_alt / n_valid) >= min_allele_freq:
                pos, r_base, a_base = int(row['POS']), str(row['REF']), str(row['ALT']).split(',')[0]
                if len(r_base) == 1 and len(a_base) == 1 and str(cluster_seq[pos - 1]).upper() == r_base.upper():
                    cluster_seq[pos - 1] = a_base
                    snps_applied += 1
        
        seq_rec = SeqIO.SeqRecord(Seq(cluster_seq), id=f"{name}_consensus", description=f"Consensus for cluster {name}")
        rep_seqs.append(seq_rec)
        stats.append({'cluster_name': name, 'found_samples': len(valid_samples), 'snps_applied': snps_applied})
    return rep_seqs, stats

def run_ani_analysis(reference_file, fasta_files, output_dir, orthoani_path):
    ani_results, all_ani_data = {}, []
    logging.info(f"Running OrthoANI for {len(fasta_files)} files...")
    for i, query_file in enumerate(fasta_files, 1):
        c_name = query_file.stem.replace('_consensus', '')
        logging.info(f"  ({i}/{len(fasta_files)}) ANI for cluster '{c_name}'...")
        try:
            cmd = [orthoani_path, "-r", str(reference_file), "-q", str(query_file)]
            res = subprocess.run(cmd, capture_output=True, text=True, check=True)
            if res.stdout.strip():
                ani_val = float(res.stdout.strip())
                ani_results[c_name] = ani_val
                all_ani_data.append({'reference': reference_file.name, 'query': query_file.name, 'ani': ani_val})
                logging.info(f"    -> ANI: {ani_val:.4f}")
        except Exception as e:
            logging.error(f"    -> OrthoANI failed for '{c_name}': {getattr(e, 'stderr', str(e)).strip()}")
    if all_ani_data:
        pd.DataFrame(all_ani_data).to_csv(output_dir / "ani_summary.tsv", sep='\t', index=False)
        logging.info(f"ANI summary table saved to: {output_dir / 'ani_summary.tsv'}")
    return ani_results

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reference_fasta", required=True, type=Path)
    parser.add_argument("--snp_matrix", required=True, type=Path)
    parser.add_argument("--cluster_dir", required=True, type=Path)
    parser.add_argument("--output_dir", required=True, type=Path)
    parser.add_argument("--min_allele_freq", type=float, default=0.51)
    parser.add_argument("--skip_ani", action='store_true')
    args = parser.parse_args()
    
    setup_logging(args.output_dir / "build_consensus.log")
    
    try:
        logging.info("Step 1: Loading cluster definitions...")
        clusters = load_cluster_files(args.cluster_dir)
        
        logging.info("\n--- Building Per-Cluster Consensus and Running ANI Analysis ---")
        per_cluster_seqs, per_cluster_stats = build_per_cluster_consensus(clusters, args.reference_fasta, args.snp_matrix, args.min_allele_freq)
        
        if not per_cluster_seqs:
            raise ValueError("Failed to generate any per-cluster consensus sequences.")

        fasta_paths = []
        for seq_rec in per_cluster_seqs:
            c_name = seq_rec.id.replace('_consensus', '')
            c_dir = args.output_dir / "per_cluster_consensus" / c_name
            c_dir.mkdir(parents=True, exist_ok=True)
            f_path = c_dir / f"{c_name}_consensus.fasta"
            SeqIO.write([seq_rec], f_path, "fasta")
            fasta_paths.append(f_path)
        logging.info(f"Saved {len(fasta_paths)} per-cluster consensus sequences.")
    
        orthoani_path = (shutil.which("pyorthoani") or shutil.which("orthoani")) if not args.skip_ani else None
        if orthoani_path:
            ani_results = run_ani_analysis(args.reference_fasta, fasta_paths, args.output_dir, orthoani_path)
            for s in per_cluster_stats:
                s['ani_to_reference'] = ani_results.get(s['cluster_name'])
        
        pd.DataFrame(per_cluster_stats).to_csv(args.output_dir / "per_cluster_stats.tsv", sep='\t', index=False)
        
        logging.info("\n" + "="*60)
        logging.info("🎉 CONSENSUS WORKFLOW COMPLETED")
        logging.info("="*60)
        
    except Exception as e:
        logging.error(f"An unexpected error occurred: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()