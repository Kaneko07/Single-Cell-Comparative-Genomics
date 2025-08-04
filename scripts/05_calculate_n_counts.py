import pandas as pd
import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Calculate N counts per SAG after applying SNP filter.")
    parser.add_argument("--input_matrix", required=True, type=Path)
    parser.add_argument("--snp_threshold_file", required=True, type=Path)
    parser.add_argument("--output_frequency", required=True, type=Path)
    args = parser.parse_args()

    args.output_frequency.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input_matrix, sep='\t')
    
    with open(args.snp_threshold_file, 'r') as f:
        snp_threshold = int(f.read().strip())

    print(f"INFO: Calculating N-counts")
    print(f"INFO: SNP threshold: {snp_threshold}")

    df_filtered = df[df['Total'] > snp_threshold]
    print(f"INFO: After SNP filtering: {len(df_filtered)} SNPs")

    meta_cols = ['POS', 'Total', 'REF', 'ALT', 'QUAL']
    special_cols = ['Reference']
    
    missing_meta_cols = [col for col in meta_cols if col not in df_filtered.columns]
    if missing_meta_cols:
        print(f"ERROR: Expected metadata columns not found: {missing_meta_cols}")
        print(f"ERROR: Actual columns: {list(df_filtered.columns)}")
        return
    
    existing_special_cols = [col for col in special_cols if col in df_filtered.columns]
    if len(existing_special_cols) != len(special_cols):
        missing = [col for col in special_cols if col not in existing_special_cols]
        print(f"WARNING: Special columns not found: {missing}")
    
    all_non_sag_cols = meta_cols + existing_special_cols
    sag_cols = [col for col in df_filtered.columns if col not in all_non_sag_cols]
    
    print(f"INFO: Metadata columns: {len(meta_cols)}")
    print(f"INFO: Special columns: {len(existing_special_cols)} ({existing_special_cols})")
    print(f"INFO: SAG columns: {len(sag_cols)}")
    print(f"INFO: SAG columns (first 5): {sag_cols[:5]}")

    n_counts = df_filtered[sag_cols].apply(lambda x: (x == 'N').sum(), axis=0)

    df_n_counts = pd.DataFrame({
        'ID': n_counts.index,
        'N_count': n_counts.values
    })

    intermediate_file = args.output_frequency.parent / "sag_n_counts.tsv"
    df_n_counts.to_csv(intermediate_file, sep='\t', index=False)
    print(f"INFO: Individual SAG N-counts saved to {intermediate_file}")

    n_count_freq = df_n_counts['N_count'].value_counts()

    df_n_count_freq = pd.DataFrame(n_count_freq).reset_index()
    df_n_count_freq.columns = ['N_Value', 'SAG_Count']

    df_n_count_freq = df_n_count_freq.sort_values('N_Value')

    df_n_count_freq.to_csv(args.output_frequency, sep='\t', index=False)
    
    print(f"INFO: SAG N-counts summary saved to {args.output_frequency}")
    print(f"INFO: Correctly excluded: QUAL, Reference (if present)")
    print(f"INFO: Pure SAG count: {len(sag_cols)}")

if __name__ == "__main__":
    main()