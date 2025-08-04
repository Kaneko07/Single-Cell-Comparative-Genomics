import pandas as pd
import numpy as np
import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Convert snp-sites VCF to SNP matrix.")
    parser.add_argument("--input_vcf", required=True, type=Path, help="Input VCF file from snp-sites")
    parser.add_argument("--output_matrix", required=True, type=Path, help="Output matrix file")
    parser.add_argument("--output_frequency", required=True, type=Path, help="Output SNP frequency file")
    args = parser.parse_args()

    args.output_matrix.parent.mkdir(parents=True, exist_ok=True)
    args.output_frequency.parent.mkdir(parents=True, exist_ok=True)

    print(f"INFO: Processing VCF file: {args.input_vcf}")

    with open(args.input_vcf, 'r') as file:
        data = file.read()

    data = data.replace('#CHROM', 'CHROM')

    temp_vcf = args.input_vcf.parent / "temp_original.vcv"
    with open(temp_vcf, 'w') as file:
        file.write(data)

    df = pd.read_csv(temp_vcf, sep='\t', comment='#')

    df_copy = df.copy()

    print(f"INFO: Data shape before processing: {df_copy.shape}")
    print(f"INFO: Columns before processing: {list(df_copy.columns)}")

    sample_columns = df_copy.columns[9:]
    for col in sample_columns:
        df_copy[col] = df_copy[col].astype(object)
    
    print(f"INFO: Sample columns converted to object type")
    
    for i, row in df_copy.iterrows():
        fourth_column_value = row.iloc[3]
        fifth_column_value = row.iloc[4]
        
        if str(fifth_column_value).startswith('*'):
            for col_idx in range(9, len(df_copy.columns)):
                if df_copy.iloc[i, col_idx] == 1:
                    df_copy.iloc[i, col_idx] = 'N'
        
        elif str(fifth_column_value).endswith('*'):
            for col_idx in range(9, len(df_copy.columns)):
                if df_copy.iloc[i, col_idx] == 2:
                    df_copy.iloc[i, col_idx] = 'N'

        if df_copy.iloc[i, -1] == 1 or df_copy.iloc[i, -1] == 2:
            df_copy.iloc[i, 3] = fifth_column_value
            df_copy.iloc[i, 4] = fourth_column_value

    for col in sample_columns:
        df_copy[col] = df_copy[col].replace(2, 1)

    reference_1_mask = df_copy.iloc[:, -1] == 1
    
    if reference_1_mask.any():
        for col in sample_columns:
            ref_1_data = df_copy.loc[reference_1_mask, col].copy()
            ref_1_data = ref_1_data.replace(0, 2)
            ref_1_data = ref_1_data.replace(1, 0)
            
            df_copy.loc[reference_1_mask, col] = ref_1_data

    df_copy.iloc[:, 3] = df_copy.iloc[:, 3].astype(str).str.replace(r'\*,', '', regex=True)
    df_copy.iloc[:, 3] = df_copy.iloc[:, 3].astype(str).str.replace(r',\*', '', regex=True)
    df_copy.iloc[:, 4] = df_copy.iloc[:, 4].astype(str).str.replace(r'\*,', '', regex=True)
    df_copy.iloc[:, 4] = df_copy.iloc[:, 4].astype(str).str.replace(r',\*', '', regex=True)

    columns_to_drop = [0, 2, 6, 7, 8]
    df_copy = df_copy.drop(df_copy.columns[columns_to_drop], axis=1)

    print(f"INFO: Data shape after column removal: {df_copy.shape}")
    print(f"INFO: Columns after removal: {list(df_copy.columns)}")

    df_temp = df_copy.replace('N', 0)
    
    calculation_cols = df_temp.columns[4:]
    for col in calculation_cols:
        df_temp[col] = pd.to_numeric(df_temp[col], errors='coerce').fillna(0)

    total = df_temp.iloc[:, 4:].sum(axis=1)

    df_copy.insert(1, 'Total', total)

    print(f"INFO: Final data shape: {df_copy.shape}")
    print(f"INFO: Final columns: {list(df_copy.columns)}")

    df_copy.to_csv(args.output_matrix, sep='\t', index=False)
    print(f"INFO: SNP matrix saved to {args.output_matrix}")

    total_counts = df_copy['Total'].value_counts().sort_index()
    freq_df = pd.DataFrame({
        'Num_Samples': total_counts.index, 
        'SNP_Count': total_counts.values
    })
    
    freq_df = freq_df[freq_df['Num_Samples'] > 0]

    if freq_df.empty:
        print("WARNING: No SNPs were shared by more than 0 samples. The SNP frequency table is empty.")
    
    freq_df.to_csv(args.output_frequency, sep='\t', index=False, header=True)
    print(f"INFO: SNP frequency table for changepoint analysis saved to {args.output_frequency}")

    temp_vcf.unlink()

    print(f"Processed SNPs: {len(df_copy)}")
    print(f"Number of samples: {len(df_copy.columns) - 5}")
    print(f"Total value range: {df_copy['Total'].min()} - {df_copy['Total'].max()}")

if __name__ == "__main__":
    main()