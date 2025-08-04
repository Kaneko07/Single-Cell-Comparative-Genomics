#!/usr/bin/env python3

import pandas as pd
import argparse
from pathlib import Path
import sys
import importlib.util

def import_msa_generator():
    script_dir = Path(__file__).parent
    msa_script = script_dir / "06_generate_final_msa.py"
    
    if not msa_script.exists():
        raise FileNotFoundError(f"MSA generation script not found: {msa_script}")
    
    spec = importlib.util.spec_from_file_location("msa_generator", msa_script)
    msa_generator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(msa_generator)
    
    return msa_generator.generate_msa_with_threshold

def create_cluster_template(output_dir, sensitivity_results, changepoint_threshold):
    successful_analyses = [r for r in sensitivity_results if r['success']]
    
    if not successful_analyses:
        print("WARNING: No successful analyses for cluster template creation")
        return None
    
    template_data = []
    
    for result in successful_analyses:
        threshold = result['threshold']
        is_changepoint = threshold == changepoint_threshold
        
        sag_count = result['sag_count']
        snp_count = result['final_snp_count']
        fasta_name = f"msa_threshold_{threshold:02d}_sags_{sag_count:03d}_snps_{snp_count:04d}.fasta"
        
        template_entry = {
            'threshold': threshold,
            'is_changepoint': is_changepoint,
            'fasta_file': fasta_name,
            'sag_count': sag_count,
            'snp_count': snp_count,
            'n_content_mean': round(result['n_content_mean'], 4),
            'cluster_count': None
        }
        template_data.append(template_entry)
    
    template_data.sort(key=lambda x: x['threshold'])
    
    template_csv = output_dir / "cluster_analysis_template.csv"
    template_df = pd.DataFrame(template_data)
    template_df.to_csv(template_csv, index=False)
    
    print(f"INFO: 1D cluster template created: {template_csv}")
    
    return template_csv

def main():
    parser = argparse.ArgumentParser(description="1D Sensitivity Analysis for SAG Pipeline")
    parser.add_argument("--input_matrix", required=True, type=Path, help="Input matrix file")
    parser.add_argument("--snp_threshold_file", required=True, type=Path, help="SNP threshold file")
    parser.add_argument("--sag_threshold_file", required=True, type=Path, help="SAG threshold file from changepoint analysis")
    parser.add_argument("--output_dir", required=True, type=Path, help="Output directory")
    parser.add_argument("--range_offset", type=int, default=10, help="Range offset around changepoint threshold (default: ±10)")
    args = parser.parse_args()
    
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    with open(args.snp_threshold_file, 'r') as f:
        snp_threshold = int(f.read().strip())
    with open(args.sag_threshold_file, 'r') as f:
        changepoint_threshold = int(f.read().strip())
    
    generate_msa_with_threshold = import_msa_generator()
    
    min_threshold = max(0, changepoint_threshold - args.range_offset)
    max_threshold = changepoint_threshold + args.range_offset
    thresholds = range(min_threshold, max_threshold + 1)
    
    print(f"INFO: 1D sensitivity analysis")
    print(f"INFO: Changepoint threshold: {changepoint_threshold}")
    print(f"INFO: Analysis range: {min_threshold} - {max_threshold}")
    print(f"INFO: SNP threshold: {snp_threshold}")
    
    msa_dir = args.output_dir / "phylogenetic_msa"
    msa_dir.mkdir(parents=True, exist_ok=True)
    
    sensitivity_results = []
    
    for threshold in thresholds:
        print(f"\n{'='*50}")
        print(f"INFO: Analyzing SAG threshold {threshold}...")
        print(f"{'='*50}")
        
        msa_prefix = msa_dir / f"msa_threshold_{threshold:02d}"
        stats_file = msa_dir / f"stats_threshold_{threshold:02d}.txt"
        
        try:
            stats = generate_msa_with_threshold(
                input_matrix_path=args.input_matrix,
                snp_threshold=snp_threshold,
                sag_threshold=threshold,
                output_prefix=msa_prefix,
                sag_n_counts_file=None,
                verbose=True,
                stats_output=stats_file
            )
            
            if stats['success']:
                print(f"✅ Threshold {threshold}: Success")
                print(f"   SAG count: {stats['sag_count']}, SNP count: {stats['final_snp_count']}")
                print(f"   Prototype effect: -{stats['prototype_effect']} SNPs")
                print(f"   Mean N content: {stats['n_content_mean']:.3f}")
                
                original_fasta = Path(stats['fasta_file'])
                original_phylip = Path(stats['phylip_file']) if stats['phylip_saved'] else None
                
                new_fasta_name = f"msa_threshold_{threshold:02d}_sags_{stats['sag_count']:03d}_snps_{stats['final_snp_count']:04d}.fasta"
                new_phylip_name = f"msa_threshold_{threshold:02d}_sags_{stats['sag_count']:03d}_snps_{stats['final_snp_count']:04d}.phylip"
                
                new_fasta_path = msa_dir / new_fasta_name
                new_phylip_path = msa_dir / new_phylip_name
                
                if original_fasta.exists():
                    original_fasta.rename(new_fasta_path)
                    stats['fasta_file'] = str(new_fasta_path)
                
                if original_phylip and original_phylip.exists():
                    original_phylip.rename(new_phylip_path)
                    stats['phylip_file'] = str(new_phylip_path)
                
            else:
                print(f"❌ Threshold {threshold}: Failed - {stats.get('error', 'Unknown error')}")
                
        except Exception as e:
            print(f"❌ Threshold {threshold}: Exception - {str(e)}")
            stats = {
                'success': False,
                'error': str(e),
                'threshold': threshold,
                'step05_snp_count': 0,
                'final_snp_count': 0,
                'prototype_effect': 0,
                'sag_count': 0,
                'special_count': 0,
                'total_sequences': 0,
                'sequence_length': 0,
                'n_content_mean': 0,
                'n_content_std': 0,
                'n_content_median': 0,
                'gc_content_mean': 0
            }
        
        result = {
            'threshold': threshold,
            'is_changepoint': threshold == changepoint_threshold,
            'success': stats['success'],
            'step05_snp_count': stats.get('step05_snp_count', 0),
            'final_snp_count': stats.get('final_snp_count', 0),
            'prototype_effect': stats.get('prototype_effect', 0),
            'sag_count': stats.get('sag_count', 0),
            'special_count': stats.get('special_count', 0),
            'total_sequences': stats.get('total_sequences', 0),
            'sequence_length': stats.get('sequence_length', 0),
            'n_content_mean': stats.get('n_content_mean', 0),
            'n_content_std': stats.get('n_content_std', 0),
            'n_content_median': stats.get('n_content_median', 0),
            'gc_content_mean': stats.get('gc_content_mean', 0)
        }
        
        sensitivity_results.append(result)
    
    results_df = pd.DataFrame(sensitivity_results)
    results_csv = args.output_dir / "sensitivity_analysis_results.csv"
    results_df.to_csv(results_csv, index=False)
    print(f"\nINFO: 1D sensitivity results saved: {results_csv}")
    
    create_cluster_template(args.output_dir, sensitivity_results, changepoint_threshold)
    
    successful = len([r for r in sensitivity_results if r['success']])
    total = len(sensitivity_results)
    
    print(f"\n=== 1D Sensitivity Analysis Summary ===")
    print(f"Total thresholds: {total}")
    print(f"Successful: {successful}")
    print(f"Failed: {total - successful}")
    print(f"Success rate: {successful/total*100:.1f}%")
    
    if successful > 0:
        successful_results = [r for r in sensitivity_results if r['success']]
        
        best_sag = max(successful_results, key=lambda x: x['sag_count'])
        best_snp = max(successful_results, key=lambda x: x['final_snp_count'])
        best_quality = min(successful_results, key=lambda x: x['n_content_mean'])
        
        print(f"\n📊 Best combinations:")
        print(f"Most SAGs: Threshold={best_sag['threshold']} ({best_sag['sag_count']} SAGs)")
        print(f"Most SNPs: Threshold={best_snp['threshold']} ({best_snp['final_snp_count']} SNPs)")
        print(f"Best quality: Threshold={best_quality['threshold']} (N={best_quality['n_content_mean']:.4f})")
    
    print(f"\n🎉 1D Analysis completed! Results saved in: {args.output_dir}")

if __name__ == "__main__":
    main()