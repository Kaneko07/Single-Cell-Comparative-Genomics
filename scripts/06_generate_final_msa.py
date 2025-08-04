import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
from pathlib import Path
import sys
import json

def calculate_msa_stats(sequences):
    import numpy as np
    
    if not sequences:
        return {
            'num_sequences': 0,
            'sequence_length': 0,
            'n_content_mean': 0,
            'n_content_std': 0,
            'n_content_median': 0,
            'gc_content_mean': 0
        }
    
    n_contents = []
    gc_contents = []
    
    for seq in sequences:
        seq_str = str(seq.seq)
        n_count = seq_str.count('N')
        n_content = n_count / len(seq_str) if len(seq_str) > 0 else 0
        n_contents.append(n_content)
        
        valid_bases = seq_str.replace('N', '')
        if len(valid_bases) > 0:
            gc_count = valid_bases.count('G') + valid_bases.count('C')
            gc_content = gc_count / len(valid_bases)
            gc_contents.append(gc_content)
    
    return {
        'num_sequences': len(sequences),
        'sequence_length': len(sequences[0].seq) if sequences else 0,
        'n_content_mean': np.mean(n_contents),
        'n_content_std': np.std(n_contents),
        'n_content_median': np.median(n_contents),
        'gc_content_mean': np.mean(gc_contents) if gc_contents else 0
    }

def generate_msa_with_threshold(input_matrix_path, snp_threshold, sag_threshold, output_prefix, 
                               sag_n_counts_file=None, verbose=True, stats_output=None):
    
    output_fasta = Path(f"{output_prefix}.fasta")
    output_phylip = Path(f"{output_prefix}.phylip")
    
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    if verbose:
        print(f"INFO: Generating MSA (SNP threshold: {snp_threshold}, SAG threshold: {sag_threshold})")
        sys.stdout.flush()

    if verbose:
        print(f"INFO: Loading data: {input_matrix_path}")
        sys.stdout.flush()
    
    df_raw = pd.read_csv(input_matrix_path, sep='\t')
    original_snp_count = len(df_raw)
    
    meta_cols = ['POS', 'Total', 'REF', 'ALT', 'QUAL']
    special_cols = ['Reference']
    
    existing_special_cols = [col for col in special_cols if col in df_raw.columns]
    
    all_non_sag_cols = meta_cols + existing_special_cols
    sag_cols = [col for col in df_raw.columns if col not in all_non_sag_cols]
    original_sag_count = len(sag_cols)

    if verbose:
        print(f"INFO: Original data - SNPs: {original_snp_count}, SAGs: {original_sag_count}")
        sys.stdout.flush()

    if verbose:
        print(f"INFO: First SNP filtering (Total > {snp_threshold})")
        sys.stdout.flush()
    
    df_snp_filtered = df_raw[df_raw['Total'] > snp_threshold].copy().reset_index(drop=True)
    first_snp_filtered_count = len(df_snp_filtered)
    
    if verbose:
        print(f"INFO: After first SNP filtering - SNPs: {first_snp_filtered_count}")
        sys.stdout.flush()

    if sag_n_counts_file is None:
        sag_n_counts_file = Path(input_matrix_path).parent.parent / "05_n_counts" / "sag_n_counts.tsv"
    
    if not sag_n_counts_file.exists():
        warning_msg = f"WARNING: Step 05 N-count file not found: {sag_n_counts_file}"
        if verbose:
            print(warning_msg)
        
        sags_to_keep = sag_cols.copy()
        if verbose:
            print(f"INFO: Using all SAGs: {len(sags_to_keep)} SAGs")
    else:
        sag_n_counts_df = pd.read_csv(sag_n_counts_file, sep='\t')
        
        valid_sag_entries = sag_n_counts_df[
            (sag_n_counts_df['N_count'] <= sag_threshold) & 
            (~sag_n_counts_df['ID'].isin(['QUAL', 'Reference']))
        ]
        candidate_sags = valid_sag_entries['ID'].tolist()
        
        sags_to_keep = [sag for sag in candidate_sags if sag in sag_cols]
        
        if verbose:
            print(f"INFO: SAG filtering (N ≤ {sag_threshold})")
            print(f"INFO: SAGs kept: {len(sags_to_keep)}")
            sys.stdout.flush()

    if len(sags_to_keep) == 0:
        if verbose:
            print(f"WARNING: No valid SAGs after filtering")
            print(f"INFO: Generating MSA with special columns only")
        
        if len(existing_special_cols) > 0:
            sags_to_keep = []
            msa_sample_cols = existing_special_cols
        else:
            if len(sag_cols) > 0:
                sags_to_keep = [sag_cols[0]]
                msa_sample_cols = sags_to_keep + existing_special_cols
                if verbose:
                    print(f"INFO: Emergency measure using first SAG: {sags_to_keep[0]}")
            else:
                return create_empty_msa_result(snp_threshold, sag_threshold, output_fasta, output_phylip, verbose, first_snp_filtered_count, original_sag_count, existing_special_cols, "No SAGs available")
    else:
        msa_sample_cols = sags_to_keep + existing_special_cols

    df_sag_filtered = df_snp_filtered[meta_cols + msa_sample_cols].copy()

    if verbose:
        print(f"INFO: Recalculating Total column")
        print(f"INFO: Recalculating SNP sharing numbers after SAG filtering...")
        sys.stdout.flush()
    
    if len(sags_to_keep) > 0:
        df_temp_calc = df_sag_filtered[sags_to_keep].replace('N', 0)
        
        for col in sags_to_keep:
            df_temp_calc[col] = pd.to_numeric(df_temp_calc[col], errors='coerce').fillna(0)
        
        recalculated_total = df_temp_calc.sum(axis=1)
        df_sag_filtered['Total'] = recalculated_total
        
        if verbose:
            print(f"INFO: Total column recalculation complete")
            print(f"INFO: Recalculated Total range: {df_sag_filtered['Total'].min()} - {df_sag_filtered['Total'].max()}")
            sys.stdout.flush()
    else:
        if verbose:
            print(f"INFO: No SAGs, maintaining Total column")

    if verbose:
        print(f"INFO: Second SNP filtering")
        print(f"INFO: Re-filtering with recalculated Total column (Total > {snp_threshold})")
        sys.stdout.flush()
    
    if len(sags_to_keep) > 0:
        df_double_filtered = df_sag_filtered[df_sag_filtered['Total'] > snp_threshold].copy().reset_index(drop=True)
    else:
        df_double_filtered = df_sag_filtered.copy()
    
    final_snp_count = len(df_double_filtered)
    
    prototype_effect = first_snp_filtered_count - final_snp_count
    
    if verbose:
        print(f"INFO: After second SNP filtering - final SNPs: {final_snp_count}")
        if prototype_effect > 0:
            print(f"INFO: Prototype effect - {prototype_effect} additional SNPs removed")
        else:
            print(f"INFO: Prototype effect - no additional SNP removal")
        sys.stdout.flush()

    if final_snp_count == 0:
        if verbose:
            print(f"WARNING: No SNPs after second filtering")
            print(f"INFO: Generating MSA with fixed sequences")
        
        return create_minimal_msa_result(msa_sample_cols, sags_to_keep, existing_special_cols, output_fasta, output_phylip, verbose, snp_threshold, sag_threshold, first_snp_filtered_count, final_snp_count, prototype_effect, original_sag_count)

    if verbose:
        print(f"INFO: Building final MSA...")
        sys.stdout.flush()
    
    sequences = []
    sequence_lengths = []
    
    for column in msa_sample_cols:            
        sequence_parts = []
        
        for i in range(len(df_double_filtered)):
            genotype = str(df_double_filtered.iloc[i][column])
            
            if genotype == '0':
                base = str(df_double_filtered.iloc[i]['REF'])
                sequence_parts.append(base)
            elif genotype == '1':
                alt_allele = str(df_double_filtered.iloc[i]['ALT'])
                base = alt_allele.split(',')[0] if ',' in alt_allele else alt_allele
                sequence_parts.append(base)
            else:
                sequence_parts.append('N')
        
        sequence = ''.join(sequence_parts)
        sequence_lengths.append(len(sequence))
        
        if column in existing_special_cols:
            description = f"special_sequence_{column}"
        else:
            description = f"SAG_{column}"
        
        seq_record = SeqRecord(Seq(sequence), id=column, description=description)
        sequences.append(seq_record)

    if verbose:
        print(f"INFO: Saving MSA files")
        sys.stdout.flush()
    
    unique_lengths = set(sequence_lengths)
    if verbose:
        print(f"INFO: Sequence length types: {unique_lengths}")
        sys.stdout.flush()
    
    if len(unique_lengths) > 1:
        if verbose:
            print(f"WARNING: Sequences have different lengths")
        from collections import Counter
        most_common_length = Counter(sequence_lengths).most_common(1)[0][0]
        if verbose:
            print(f"INFO: Most common sequence length: {most_common_length}")
        
        filtered_sequences = [seq for seq, length in zip(sequences, sequence_lengths) if length == most_common_length]
        if verbose:
            print(f"INFO: Sequences after filtering: {len(filtered_sequences)}")
        sequences = filtered_sequences
        sys.stdout.flush()

    if len(sequences) > 0:
        SeqIO.write(sequences, output_fasta, "fasta")
    else:
        with open(output_fasta, 'w') as f:
            f.write("")
    
    phylip_saved = False
    if len(unique_lengths) == 1 and len(sequences) > 0:
        try:
            SeqIO.write(sequences, output_phylip, "phylip")
            phylip_saved = True
            if verbose:
                print(f"INFO: MSA saved: {output_fasta} and {output_phylip}")
        except Exception as e:
            if verbose:
                print(f"WARNING: Failed to save PHYLIP format: {e}")
    else:
        if verbose:
            print(f"WARNING: Different sequence lengths or no sequences, skipping PHYLIP: {output_fasta} only")
        output_phylip.touch()
    
    msa_stats = calculate_msa_stats(sequences)
    
    stats = {
        'success': True,
        'threshold': sag_threshold,
        'original_snp_count': original_snp_count,
        'original_sag_count': original_sag_count,
        'step05_snp_count': first_snp_filtered_count,
        'final_snp_count': final_snp_count,
        'prototype_effect': prototype_effect,
        'sag_count': len(sags_to_keep),
        'special_count': len(existing_special_cols),
        'total_sequences': len(sequences),
        'sequence_length': msa_stats['sequence_length'],
        'n_content_mean': msa_stats['n_content_mean'],
        'n_content_std': msa_stats['n_content_std'],
        'n_content_median': msa_stats['n_content_median'],
        'gc_content_mean': msa_stats['gc_content_mean'],
        'fasta_file': str(output_fasta),
        'phylip_file': str(output_phylip) if phylip_saved else None,
        'phylip_saved': phylip_saved
    }
    
    if stats_output:
        stats_output = Path(stats_output)
        stats_output.parent.mkdir(parents=True, exist_ok=True)
        
        with open(stats_output, 'w') as f:
            f.write("=== Double Filtering Statistics ===\n")
            f.write(f"SAG threshold: {sag_threshold}\n")
            f.write(f"SNP threshold: {snp_threshold}\n\n")
            f.write(f"Original data: SNP={original_snp_count}, SAG={original_sag_count}\n")
            f.write(f"First SNP filtering (Total > {snp_threshold}): {first_snp_filtered_count} SNPs\n")
            f.write(f"SAG filtering (N ≤ {sag_threshold}): {len(sags_to_keep)} SAGs kept\n")
            f.write(f"Special columns (included in MSA): {len(existing_special_cols)} ({', '.join(existing_special_cols)})\n")
            f.write(f"Total column recalculation: Re-evaluated SNP sharing with SAGs only\n")
            f.write(f"Second SNP filtering (Total > {snp_threshold}): {final_snp_count} SNPs\n")
            f.write(f"Final MSA: SNP={final_snp_count}, sequences={len(sequences)} (SAG={len(sags_to_keep)} + special={len(existing_special_cols)})\n")
            f.write(f"Sequence length: {msa_stats['sequence_length']}\n")
            f.write(f"Prototype effect: {prototype_effect} additional SNPs removed\n")
            
            snp_retention_from_step05 = final_snp_count / first_snp_filtered_count if first_snp_filtered_count > 0 else 0
            sag_retention_rate = len(sags_to_keep) / original_sag_count if original_sag_count > 0 else 0
            
            f.write(f"Retention rate from first SNP filter: {snp_retention_from_step05:.3f}\n")
            f.write(f"SAG retention rate: {sag_retention_rate:.3f}\n")
            f.write(f"Average N content: {msa_stats['n_content_mean']:.4f}\n")
            f.write(f"Average GC content: {msa_stats['gc_content_mean']:.4f}\n")

    if verbose:
        print(f"\n=== Double Filtering Complete ===")
        print(f"INFO: SAG threshold: {sag_threshold}")
        print(f"INFO: Final result - SNPs: {final_snp_count}, SAGs: {len(sags_to_keep)}, sequences: {len(sequences)}")
        print(f"INFO: Prototype effect: +{prototype_effect} SNPs removed")
        snp_retention_from_step05 = final_snp_count / first_snp_filtered_count if first_snp_filtered_count > 0 else 0
        sag_retention_rate = len(sags_to_keep) / original_sag_count if original_sag_count > 0 else 0
        print(f"INFO: SAG retention rate: {sag_retention_rate:.1%}")
        print(f"INFO: Average N content: {msa_stats['n_content_mean']:.3f}")
        print(f"INFO: MSA files generated")
        sys.stdout.flush()
    
    return stats

def create_empty_msa_result(snp_threshold, sag_threshold, output_fasta, output_phylip, verbose, first_snp_filtered_count, original_sag_count, existing_special_cols, error_msg):
    if verbose:
        print(f"INFO: Creating empty MSA files: {error_msg}")
    
    with open(output_fasta, 'w') as f:
        f.write("")
    with open(output_phylip, 'w') as f:
        f.write("")
    
    return {
        'success': True,
        'threshold': sag_threshold,
        'original_snp_count': 0,
        'original_sag_count': original_sag_count,
        'step05_snp_count': first_snp_filtered_count,
        'final_snp_count': 0,
        'prototype_effect': first_snp_filtered_count,
        'sag_count': 0,
        'special_count': len(existing_special_cols),
        'total_sequences': 0,
        'sequence_length': 0,
        'n_content_mean': 0,
        'n_content_std': 0,
        'n_content_median': 0,
        'gc_content_mean': 0,
        'fasta_file': str(output_fasta),
        'phylip_file': str(output_phylip),
        'phylip_saved': False,
        'error_note': error_msg
    }

def create_minimal_msa_result(msa_sample_cols, sags_to_keep, existing_special_cols, output_fasta, output_phylip, verbose, snp_threshold, sag_threshold, first_snp_filtered_count, final_snp_count, prototype_effect, original_sag_count):
    if verbose:
        print(f"INFO: Creating minimal MSA files (fixed sequences)")
    
    sequences = []
    
    for column in msa_sample_cols:
        sequence = "A"
        
        if column in existing_special_cols:
            description = f"special_sequence_{column}_fixed"
        else:
            description = f"SAG_{column}_fixed"
        
        seq_record = SeqRecord(Seq(sequence), id=column, description=description)
        sequences.append(seq_record)
    
    if len(sequences) > 0:
        SeqIO.write(sequences, output_fasta, "fasta")
        try:
            SeqIO.write(sequences, output_phylip, "phylip")
            phylip_saved = True
        except:
            phylip_saved = False
            output_phylip.touch()
    else:
        with open(output_fasta, 'w') as f:
            f.write("")
        with open(output_phylip, 'w') as f:
            f.write("")
        phylip_saved = False
    
    msa_stats = calculate_msa_stats(sequences)
    
    return {
        'success': True,
        'threshold': sag_threshold,
        'original_snp_count': 0,
        'original_sag_count': original_sag_count,
        'step05_snp_count': first_snp_filtered_count,
        'final_snp_count': final_snp_count,
        'prototype_effect': prototype_effect,
        'sag_count': len(sags_to_keep),
        'special_count': len(existing_special_cols),
        'total_sequences': len(sequences),
        'sequence_length': msa_stats['sequence_length'],
        'n_content_mean': msa_stats['n_content_mean'],
        'n_content_std': msa_stats['n_content_std'],
        'n_content_median': msa_stats['n_content_median'],
        'gc_content_mean': msa_stats['gc_content_mean'],
        'fasta_file': str(output_fasta),
        'phylip_file': str(output_phylip) if phylip_saved else None,
        'phylip_saved': phylip_saved,
        'error_note': 'Minimal MSA with fixed sequences'
    }

def main():
    parser = argparse.ArgumentParser(description="Generate final MSA with parameterized thresholds.")
    parser.add_argument("--input_matrix", required=True, type=Path)
    parser.add_argument("--snp_threshold_file", required=True, type=Path)
    parser.add_argument("--sag_threshold_file", required=True, type=Path)
    parser.add_argument("--output_prefix", required=True, type=Path)
    parser.add_argument("--sag_n_counts_file", type=Path, help="Step 05 N-counts file (auto-detected if not provided)")
    parser.add_argument("--stats_output", type=Path, help="Statistics output file")
    args = parser.parse_args()

    with open(args.snp_threshold_file, 'r') as f:
        snp_threshold = int(f.read().strip())
    with open(args.sag_threshold_file, 'r') as f:
        sag_threshold = int(f.read().strip())

    if args.stats_output is None:
        args.stats_output = args.output_prefix.parent / "double_filtering_stats.txt"

    print("INFO: Force execution version - generating MSA files for all conditions")

    stats = generate_msa_with_threshold(
        input_matrix_path=args.input_matrix,
        snp_threshold=snp_threshold,
        sag_threshold=sag_threshold,
        output_prefix=args.output_prefix,
        sag_n_counts_file=args.sag_n_counts_file,
        verbose=True,
        stats_output=args.stats_output
    )

    print(f"INFO: Force execution complete - MSA files generated: {stats.get('fasta_file', 'N/A')}")

if __name__ == "__main__":
    main()