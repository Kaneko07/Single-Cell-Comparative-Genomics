#!/usr/bin/env python3

import argparse
import sys
import os
import json
import pandas as pd
import subprocess
from pathlib import Path
import time
from datetime import datetime
import logging
import importlib.util

class True2DSensitivityAnalysis:
    
    def __init__(self, args):
        self.args = args
        
        self.input_matrix = Path(args.input_matrix)
        self.output_dir = Path(args.output_dir)
        
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.logger = self.setup_logging()
        self.validate_inputs()
        
        self.load_initial_snp_threshold()
        
        self.snp_center = self.initial_snp_threshold
        self.snp_range = args.snp_range
        self.sag_range = args.sag_range
        
        self.sensitivity_results = []
        self.changepoint_cache = {}
        
        self.generate_msa_with_threshold = self.import_msa_generator()
        
    def setup_logging(self):
        level = logging.WARNING if self.args.quiet else logging.INFO
        
        log_dir = self.output_dir / "logs"
        log_dir.mkdir(parents=True, exist_ok=True)
        
        logging.basicConfig(
            level=level,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_dir / "true_2d_sensitivity.log"),
                logging.StreamHandler() if not self.args.quiet else logging.NullHandler()
            ]
        )
        return logging.getLogger(__name__)
    
    def validate_inputs(self):
        if not self.input_matrix.exists():
            raise FileNotFoundError(f"Input matrix not found: {self.input_matrix}")
            
        if not Path(self.args.snp_threshold_file).exists():
            raise FileNotFoundError(f"SNP threshold file not found: {self.args.snp_threshold_file}")
    
    def load_initial_snp_threshold(self):
        with open(self.args.snp_threshold_file, 'r') as f:
            self.initial_snp_threshold = int(f.read().strip())
        
        self.logger.info(f"Loaded initial SNP changepoint threshold: {self.initial_snp_threshold}")
    
    def import_msa_generator(self):
        script_dir = Path(__file__).parent
        msa_script = script_dir / "06_generate_final_msa.py"
        
        if not msa_script.exists():
            raise FileNotFoundError(f"MSA generation script not found: {msa_script}")
        
        spec = importlib.util.spec_from_file_location("msa_generator", msa_script)
        msa_generator = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(msa_generator)
        
        return msa_generator.generate_msa_with_threshold
    
    def generate_snp_values(self):
        snp_values = list(range(
            max(1, self.snp_center - self.snp_range),
            self.snp_center + self.snp_range + 1
        ))
        return snp_values
    
    def calculate_n_counts_for_snp_threshold(self, snp_threshold):
        self.logger.info(f"Calculating N-counts for SNP threshold: {snp_threshold}")
        
        temp_dir = self.output_dir / "temp_calculations" / f"snp_{snp_threshold}"
        temp_dir.mkdir(parents=True, exist_ok=True)
        
        try:
            df = pd.read_csv(self.input_matrix, sep='\t')
            
            df_filtered = df[df['Total'] > snp_threshold]
            
            if len(df_filtered) == 0:
                self.logger.warning(f"No SNPs remain after filtering with threshold {snp_threshold}")
                return None
            
            meta_cols = ['POS', 'Total', 'REF', 'ALT', 'QUAL']
            special_cols = ['Reference']
            existing_special_cols = [col for col in special_cols if col in df_filtered.columns]
            all_non_sag_cols = meta_cols + existing_special_cols
            sag_cols = [col for col in df_filtered.columns if col not in all_non_sag_cols]
            
            n_counts = df_filtered[sag_cols].apply(lambda x: (x == 'N').sum(), axis=0)
            
            n_count_freq = n_counts.value_counts().sort_index()
            
            freq_df = pd.DataFrame({
                'N_Value': n_count_freq.index,
                'SAG_Count': n_count_freq.values
            })
            
            freq_file = temp_dir / "sag_n_counts_frequency.tsv"
            freq_df.to_csv(freq_file, sep='\t', index=False)
            
            self.logger.info(f"N-counts calculated for SNP threshold {snp_threshold}: {len(sag_cols)} SAGs, {len(df_filtered)} SNPs")
            
            return freq_file
            
        except Exception as e:
            self.logger.error(f"Failed to calculate N-counts for SNP threshold {snp_threshold}: {str(e)}")
            return None
    
    def run_changepoint_analysis_for_sag(self, freq_file, snp_threshold):
        self.logger.info(f"Running changepoint analysis for SNP threshold: {snp_threshold}")
        
        temp_dir = freq_file.parent
        sag_threshold_file = temp_dir / "sag_threshold.txt"
        sag_plot_file = temp_dir / "changepoint_sag_plot.png"
        
        script_dir = Path(__file__).parent
        r_script = script_dir / "04_run_changepoint.R"
        
        if not r_script.exists():
            self.logger.error(f"Changepoint R script not found: {r_script}")
            return None
        
        try:
            cmd = [
                "Rscript", str(r_script),
                str(freq_file), str(sag_threshold_file), str(sag_plot_file),
                "Number of 'N's per SAG", "Frequency of SAGs"
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            
            if result.returncode == 0 and sag_threshold_file.exists():
                with open(sag_threshold_file, 'r') as f:
                    sag_threshold = int(f.read().strip())
                
                self.logger.info(f"SAG changepoint for SNP={snp_threshold}: {sag_threshold}")
                return sag_threshold
            else:
                self.logger.error(f"Changepoint analysis failed for SNP threshold {snp_threshold}: {result.stderr}")
                return None
                
        except Exception as e:
            self.logger.error(f"Exception in changepoint analysis for SNP threshold {snp_threshold}: {str(e)}")
            return None
    
    def get_sag_changepoint_for_snp(self, snp_threshold):
        
        if snp_threshold in self.changepoint_cache:
            return self.changepoint_cache[snp_threshold]
        
        freq_file = self.calculate_n_counts_for_snp_threshold(snp_threshold)
        if freq_file is None:
            return None
        
        sag_threshold = self.run_changepoint_analysis_for_sag(freq_file, snp_threshold)
        
        if sag_threshold is not None:
            self.changepoint_cache[snp_threshold] = sag_threshold
        
        return sag_threshold
    
    def generate_sag_values_for_snp(self, snp_threshold, sag_changepoint):
        if sag_changepoint is None:
            return []
        
        sag_values = list(range(
            max(0, sag_changepoint - self.sag_range),
            sag_changepoint + self.sag_range + 1
        ))
        return sag_values
    
    def run_single_analysis(self, snp_threshold, sag_threshold, is_snp_changepoint, is_sag_changepoint, combination_id, total_combinations):
        
        analysis_output = self.output_dir / "sensitivity_runs" / f"snp_{snp_threshold}_sag_{sag_threshold}"
        analysis_output.mkdir(parents=True, exist_ok=True)
        
        is_both_changepoint = (is_snp_changepoint and is_sag_changepoint)
        
        if is_both_changepoint:
            self.logger.info(f"🎯 [{combination_id}/{total_combinations}] TRUE CHANGEPOINT: SNP={snp_threshold}, SAG={sag_threshold}")
        else:
            self.logger.info(f"[{combination_id}/{total_combinations}] Starting: SNP={snp_threshold}, SAG={sag_threshold}")
        
        start_time = time.time()
        
        try:
            msa_prefix = analysis_output / f"msa_snp_{snp_threshold}_sag_{sag_threshold}"
            stats_file = analysis_output / "analysis_stats.txt"
            
            verbose_mode = False if not is_both_changepoint else True
            
            self.logger.info(f"Calling MSA generation with: matrix={self.input_matrix}, snp_thresh={snp_threshold}, sag_thresh={sag_threshold}")
            
            stats = self.generate_msa_with_threshold(
                input_matrix_path=self.input_matrix,
                snp_threshold=snp_threshold,
                sag_threshold=sag_threshold,
                output_prefix=msa_prefix,
                sag_n_counts_file=None,
                verbose=verbose_mode,
                stats_output=stats_file
            )
            
            elapsed_time = time.time() - start_time
            
            if stats['success']:
                results = self.enhance_analysis_results(
                    stats, snp_threshold, sag_threshold, is_snp_changepoint, is_sag_changepoint, elapsed_time
                )
                
                if is_both_changepoint:
                    self.logger.info(f"🎯✅ [{combination_id}/{total_combinations}] TRUE CHANGEPOINT SUCCESS: SNP={snp_threshold}, SAG={sag_threshold} → SAGs={stats['sag_count']}, SNPs={stats['final_snp_count']} ({elapsed_time:.1f}s)")
                else:
                    self.logger.info(f"✅ [{combination_id}/{total_combinations}] Success: SNP={snp_threshold}, SAG={sag_threshold} ({elapsed_time:.1f}s)")
                
                return results
            else:
                if is_both_changepoint:
                    self.logger.error(f"🎯❌ [{combination_id}/{total_combinations}] TRUE CHANGEPOINT FAILED: SNP={snp_threshold}, SAG={sag_threshold} → {stats.get('error', 'MSA generation failed')}")
                else:
                    self.logger.error(f"❌ [{combination_id}/{total_combinations}] Failed: SNP={snp_threshold}, SAG={sag_threshold}")
                
                return self.create_failure_result(snp_threshold, sag_threshold, is_snp_changepoint, is_sag_changepoint, elapsed_time, stats.get('error', 'MSA generation failed'))
                
        except Exception as e:
            elapsed_time = time.time() - start_time
            if is_both_changepoint:
                self.logger.error(f"🎯💥 [{combination_id}/{total_combinations}] TRUE CHANGEPOINT EXCEPTION: SNP={snp_threshold}, SAG={sag_threshold}: {str(e)}")
            else:
                self.logger.error(f"💥 [{combination_id}/{total_combinations}] Exception: SNP={snp_threshold}, SAG={sag_threshold}: {str(e)}")
            
            return self.create_failure_result(snp_threshold, sag_threshold, is_snp_changepoint, is_sag_changepoint, elapsed_time, str(e))
    
    def enhance_analysis_results(self, stats, snp_threshold, sag_threshold, is_snp_changepoint, is_sag_changepoint, elapsed_time):
        results = {
            'snp_threshold': snp_threshold,
            'sag_threshold': sag_threshold,
            'is_snp_changepoint': is_snp_changepoint,
            'is_sag_changepoint': is_sag_changepoint,
            'is_both_changepoint': (is_snp_changepoint and is_sag_changepoint),
            'success': True,
            'elapsed_time': elapsed_time,
            'timestamp': datetime.now().isoformat(),
            
            'sag_count': stats['sag_count'],
            'special_count': stats['special_count'],
            'total_sequences': stats['total_sequences'],
            'sequence_length': stats['sequence_length'],
            
            'original_snp_count': stats.get('original_snp_count', 0),
            'step05_snp_count': stats['step05_snp_count'],
            'final_snp_count': stats['final_snp_count'],
            'prototype_effect': stats['prototype_effect'],
            
            'n_content_mean': stats['n_content_mean'],
            'n_content_std': stats['n_content_std'],
            'n_content_median': stats['n_content_median'],
            'gc_content_mean': stats['gc_content_mean'],
            
            'fasta_file': stats.get('fasta_file'),
            'phylip_file': stats.get('phylip_file'),
            'phylip_saved': stats.get('phylip_saved', False),
            
            'calculated_sag_changepoint': self.changepoint_cache.get(snp_threshold),
            
            'cluster_count': None
        }
        
        return results
    
    def create_failure_result(self, snp_threshold, sag_threshold, is_snp_changepoint, is_sag_changepoint, elapsed_time, error_msg):
        return {
            'snp_threshold': snp_threshold,
            'sag_threshold': sag_threshold,
            'is_snp_changepoint': is_snp_changepoint,
            'is_sag_changepoint': is_sag_changepoint,
            'is_both_changepoint': False,
            'success': False,
            'elapsed_time': elapsed_time,
            'timestamp': datetime.now().isoformat(),
            'error': error_msg,
            'sag_count': 0,
            'special_count': 0,
            'total_sequences': 0,
            'sequence_length': 0,
            'original_snp_count': 0,
            'step05_snp_count': 0,
            'final_snp_count': 0,
            'prototype_effect': 0,
            'n_content_mean': 0,
            'n_content_std': 0,
            'n_content_median': 0,
            'gc_content_mean': 0,
            'fasta_file': None,
            'phylip_file': None,
            'phylip_saved': False,
            'calculated_sag_changepoint': self.changepoint_cache.get(snp_threshold),
            'cluster_count': None
        }
    
    def run_true_2d_sensitivity_analysis(self):
        
        snp_values = self.generate_snp_values()
        
        self.logger.info(f"Starting TRUE 2D sensitivity analysis")
        self.logger.info(f"Initial SNP changepoint: {self.initial_snp_threshold}")
        self.logger.info(f"SNP thresholds to analyze: {snp_values}")
        self.logger.info(f"SAG range around each changepoint: ±{self.sag_range}")
        
        start_time = time.time()
        results = []
        combination_id = 0
        
        self.logger.info("Phase 1: Calculating SAG changepoints for each SNP threshold...")
        for snp_threshold in snp_values:
            sag_changepoint = self.get_sag_changepoint_for_snp(snp_threshold)
            if sag_changepoint is not None:
                self.logger.info(f"SNP={snp_threshold} → SAG changepoint={sag_changepoint}")
                if snp_threshold == self.initial_snp_threshold:
                    self.logger.info(f"🎯 INITIAL SNP CHANGEPOINT: SNP={snp_threshold} → SAG={sag_changepoint}")
            else:
                self.logger.warning(f"Failed to calculate SAG changepoint for SNP={snp_threshold}")
        
        initial_sag_changepoint = self.changepoint_cache.get(self.initial_snp_threshold)
        if initial_sag_changepoint is not None:
            self.logger.info(f"🔍 TRUE CHANGEPOINT COMBINATION WILL BE: SNP={self.initial_snp_threshold}, SAG={initial_sag_changepoint}")
        else:
            self.logger.error(f"❌ CRITICAL: Could not determine SAG changepoint for initial SNP threshold {self.initial_snp_threshold}")
        
        self.logger.info("Phase 2: Running MSA analyses for all combinations...")
        
        total_combinations = 0
        for snp_threshold in snp_values:
            sag_changepoint = self.changepoint_cache.get(snp_threshold)
            if sag_changepoint is not None:
                sag_values = self.generate_sag_values_for_snp(snp_threshold, sag_changepoint)
                total_combinations += len(sag_values)
        
        self.logger.info(f"Total combinations to analyze: {total_combinations}")
        
        changepoint_combinations_found = []
        
        for snp_threshold in snp_values:
            is_snp_changepoint = (snp_threshold == self.initial_snp_threshold)
            sag_changepoint = self.changepoint_cache.get(snp_threshold)
            
            if sag_changepoint is None:
                self.logger.warning(f"Skipping SNP threshold {snp_threshold} (no SAG changepoint)")
                continue
            
            sag_values = self.generate_sag_values_for_snp(snp_threshold, sag_changepoint)
            
            self.logger.info(f"Analyzing SNP={snp_threshold} (is_changepoint={is_snp_changepoint}) with SAG values: {sag_values}")
            
            for sag_threshold in sag_values:
                combination_id += 1
                is_sag_changepoint = (sag_threshold == sag_changepoint)
                is_both_changepoint = (is_snp_changepoint and is_sag_changepoint)
                
                if is_both_changepoint:
                    self.logger.info(f"🎯 PROCESSING TRUE CHANGEPOINT COMBINATION: SNP={snp_threshold}, SAG={sag_threshold}")
                    changepoint_combinations_found.append((snp_threshold, sag_threshold))
                elif is_snp_changepoint:
                    self.logger.info(f"🔸 Processing SNP changepoint: SNP={snp_threshold}, SAG={sag_threshold}")
                elif is_sag_changepoint:
                    self.logger.info(f"🔹 Processing SAG changepoint: SNP={snp_threshold}, SAG={sag_threshold}")
                
                result = self.run_single_analysis(
                    snp_threshold, sag_threshold, 
                    is_snp_changepoint, is_sag_changepoint,
                    combination_id, total_combinations
                )
                results.append(result)
                
                if is_both_changepoint:
                    if result['success']:
                        self.logger.info(f"✅ TRUE CHANGEPOINT COMBINATION SUCCESS: SNP={snp_threshold}, SAG={sag_threshold} → SAGs={result['sag_count']}, SNPs={result['final_snp_count']}")
                    else:
                        self.logger.error(f"❌ TRUE CHANGEPOINT COMBINATION FAILED: SNP={snp_threshold}, SAG={sag_threshold} → {result.get('error', 'Unknown error')}")
                
                if combination_id % 5 == 0:
                    self.save_sensitivity_results(results)
        
        if changepoint_combinations_found:
            self.logger.info(f"✅ TRUE CHANGEPOINT COMBINATIONS PROCESSED: {changepoint_combinations_found}")
        else:
            self.logger.warning(f"⚠️ NO TRUE CHANGEPOINT COMBINATIONS FOUND!")
            self.logger.warning(f"Initial SNP threshold: {self.initial_snp_threshold}")
            self.logger.warning(f"Changepoint cache: {self.changepoint_cache}")
        
        total_elapsed = time.time() - start_time
        self.logger.info(f"TRUE 2D sensitivity analysis completed in {total_elapsed:.1f} seconds")
        
        self.sensitivity_results = results
        self.save_sensitivity_results(results)
        self.save_changepoint_summary()
        
        self.ensure_changepoint_combination_processed(results)
        
        self.create_cluster_template(results)
        self.create_phylogenetic_directory(results)
        
        return results
    
    def ensure_changepoint_combination_processed(self, results):
        
        expected_snp = self.initial_snp_threshold
        expected_sag = self.changepoint_cache.get(expected_snp)
        
        if expected_sag is None:
            self.logger.warning(f"⚠️ Cannot verify changepoint combination: SAG changepoint not calculated for SNP={expected_snp}")
            return
        
        changepoint_result = None
        for result in results:
            if (result['snp_threshold'] == expected_snp and 
                result['sag_threshold'] == expected_sag):
                changepoint_result = result
                break
        
        if changepoint_result is None:
            self.logger.error(f"🚨 CRITICAL ERROR: TRUE CHANGEPOINT COMBINATION NOT FOUND!")
            self.logger.error(f"Expected: SNP={expected_snp}, SAG={expected_sag}")
            self.logger.error(f"This combination should have been processed but is missing from results.")
            
            self.logger.info(f"🔧 FORCE EXECUTING TRUE CHANGEPOINT COMBINATION...")
            try:
                forced_result = self.run_single_analysis(
                    expected_snp, expected_sag, 
                    True, True,
                    len(results) + 1, len(results) + 1
                )
                
                results.append(forced_result)
                self.sensitivity_results = results
                self.save_sensitivity_results(results)
                
                if forced_result['success']:
                    self.logger.info(f"✅ FORCE EXECUTION SUCCESS: TRUE CHANGEPOINT COMBINATION COMPLETED")
                else:
                    self.logger.error(f"❌ FORCE EXECUTION FAILED: {forced_result.get('error', 'Unknown error')}")
                    
            except Exception as e:
                self.logger.error(f"💥 FORCE EXECUTION EXCEPTION: {str(e)}")
                
        else:
            if changepoint_result.get('is_both_changepoint', False):
                if changepoint_result['success']:
                    self.logger.info(f"✅ VERIFIED: TRUE CHANGEPOINT COMBINATION WAS PROCESSED SUCCESSFULLY")
                    self.logger.info(f"Result: SNP={expected_snp}, SAG={expected_sag} → SAGs={changepoint_result['sag_count']}, SNPs={changepoint_result['final_snp_count']}")
                else:
                    self.logger.warning(f"⚠️ TRUE CHANGEPOINT COMBINATION WAS PROCESSED BUT FAILED")
                    self.logger.warning(f"Error: {changepoint_result.get('error', 'Unknown error')}")
            else:
                self.logger.warning(f"⚠️ CHANGEPOINT COMBINATION FOUND BUT NOT MARKED AS BOTH_CHANGEPOINT")
                changepoint_result['is_both_changepoint'] = True
                self.logger.info(f"🔧 FIXED: Marked combination as both_changepoint")
    
    def save_changepoint_summary(self):
        summary_data = []
        for snp_threshold, sag_changepoint in self.changepoint_cache.items():
            summary_data.append({
                'snp_threshold': snp_threshold,
                'calculated_sag_changepoint': sag_changepoint,
                'is_initial_snp_changepoint': snp_threshold == self.initial_snp_threshold
            })
        
        summary_df = pd.DataFrame(summary_data)
        summary_file = self.output_dir / "changepoint_summary.csv"
        summary_df.to_csv(summary_file, index=False)
        
        self.logger.info(f"Changepoint summary saved: {summary_file}")
    
    def save_sensitivity_results(self, results):
        if not results:
            return
        
        df = pd.DataFrame(results)
        
        column_order = [
            'snp_threshold', 'sag_threshold', 'calculated_sag_changepoint',
            'is_snp_changepoint', 'is_sag_changepoint', 'is_both_changepoint',
            'success', 'elapsed_time', 'timestamp',
            'sag_count', 'special_count', 'total_sequences', 'sequence_length',
            'original_snp_count', 'step05_snp_count', 'final_snp_count', 'prototype_effect',
            'n_content_mean', 'n_content_std', 'n_content_median', 'gc_content_mean',
            'fasta_file', 'phylip_file', 'phylip_saved',
            'cluster_count', 'error'
        ]
        
        for col in column_order:
            if col not in df.columns:
                df[col] = None
        
        df = df[column_order]
        
        results_file = self.output_dir / "true_2d_sensitivity_results.csv"
        df.to_csv(results_file, index=False)
        
        if not self.args.quiet:
            self.logger.info(f"TRUE 2D sensitivity results saved: {results_file}")
    
    def create_cluster_template(self, results):
        successful = [r for r in results if r['success']]
        
        if not successful:
            return
        
        template_data = []
        for result in successful:
            snp_thresh = result['snp_threshold']
            sag_thresh = result['sag_threshold']
            
            sag_count = result['sag_count']
            snp_count = result['final_snp_count']
            fasta_name = f"msa_snp_{snp_thresh:02d}_sag_{sag_thresh:02d}_sags_{sag_count:03d}_snps_{snp_count:04d}.fasta"
            
            template_data.append({
                'snp_threshold': snp_thresh,
                'sag_threshold': sag_thresh,
                'calculated_sag_changepoint': result['calculated_sag_changepoint'],
                'parameter_pair': f"SNP{snp_thresh}_SAG{sag_thresh}",
                'is_snp_changepoint': result['is_snp_changepoint'],
                'is_sag_changepoint': result['is_sag_changepoint'],
                'is_both_changepoint': result['is_both_changepoint'],
                'fasta_file': fasta_name,
                'sag_count': sag_count,
                'final_snp_count': snp_count,
                'n_content_mean': round(result['n_content_mean'], 4),
                'cluster_count': None
            })
        
        template_data.sort(key=lambda x: (x['snp_threshold'], x['sag_threshold']))
        
        df = pd.DataFrame(template_data)
        template_file = self.output_dir / "cluster_analysis_template_true_2d.csv"
        df.to_csv(template_file, index=False)
        
        self.logger.info(f"TRUE 2D cluster template created: {template_file}")
    
    def create_phylogenetic_directory(self, results):
        successful = [r for r in results if r['success']]
        
        if not successful:
            self.logger.warning("No successful results for phylogenetic directory creation")
            return
        
        phylo_dir = self.output_dir / "phylogenetic_msa_true_2d"
        phylo_dir.mkdir(parents=True, exist_ok=True)
        
        copied_files = []
        missing_files = []
        
        self.logger.info(f"Creating phylogenetic directory with {len(successful)} successful results")
        
        for i, result in enumerate(successful, 1):
            snp_thresh = result['snp_threshold']
            sag_thresh = result['sag_threshold']
            
            self.logger.info(f"[{i}/{len(successful)}] Processing: SNP={snp_thresh}, SAG={sag_thresh}")
            
            if result['fasta_file'] and Path(result['fasta_file']).exists():
                source_file = Path(result['fasta_file'])
                
                sag_count = result['sag_count']
                snp_count = result['final_snp_count']
                
                priority = ""
                if result['is_both_changepoint']:
                    priority = "PRIORITY1_BOTH_CP_"
                elif result['is_snp_changepoint']:
                    priority = "PRIORITY2_SNP_CP_"
                elif result['is_sag_changepoint']:
                    priority = "PRIORITY2_SAG_CP_"
                elif (result['sag_count'] >= 50 and result['n_content_mean'] <= 0.15 and result['final_snp_count'] >= 100):
                    priority = "PRIORITY3_QUALITY_"
                
                new_filename = f"{priority}msa_snp_{snp_thresh:02d}_sag_{sag_thresh:02d}_sags_{sag_count:03d}_snps_{snp_count:04d}.fasta"
                dest_file = phylo_dir / new_filename
                
                try:
                    import shutil
                    shutil.copy2(source_file, dest_file)
                    copied_files.append(new_filename)
                    self.logger.info(f"  ✅ Copied: {new_filename}")
                except Exception as e:
                    self.logger.error(f"  ❌ Failed to copy {source_file}: {e}")
                    missing_files.append(f"SNP={snp_thresh}_SAG={sag_thresh} (copy failed)")
            else:
                missing_file_info = f"SNP={snp_thresh}_SAG={sag_thresh} (source file missing: {result.get('fasta_file', 'No file path')})"
                missing_files.append(missing_file_info)
                self.logger.warning(f"  ⚠️ Source file missing: {result.get('fasta_file', 'No file path')}")
        
        self.verify_combination_completeness(results, copied_files, missing_files)
        
        self.logger.info(f"TRUE 2D phylogenetic MSA directory created: {phylo_dir}")
        self.logger.info(f"Successfully copied: {len(copied_files)} MSA files")
        
        if missing_files:
            self.logger.warning(f"Missing files: {len(missing_files)}")
            for missing in missing_files:
                self.logger.warning(f"  - {missing}")
    
    def verify_combination_completeness(self, results, copied_files, missing_files):
        
        self.logger.info("\n=== COMBINATION COMPLETENESS VERIFICATION ===")
        
        snp_values = self.generate_snp_values()
        expected_combinations = []
        
        for snp_threshold in snp_values:
            sag_changepoint = self.changepoint_cache.get(snp_threshold)
            if sag_changepoint is not None:
                sag_values = self.generate_sag_values_for_snp(snp_threshold, sag_changepoint)
                for sag_threshold in sag_values:
                    expected_combinations.append((snp_threshold, sag_threshold))
        
        self.logger.info(f"Expected combinations: {len(expected_combinations)}")
        
        actual_combinations = [(r['snp_threshold'], r['sag_threshold']) for r in results]
        successful_combinations = [(r['snp_threshold'], r['sag_threshold']) for r in results if r['success']]
        
        self.logger.info(f"Actual processed: {len(actual_combinations)}")
        self.logger.info(f"Successful: {len(successful_combinations)}")
        self.logger.info(f"Failed: {len(actual_combinations) - len(successful_combinations)}")
        
        missing_combinations = []
        for expected in expected_combinations:
            if expected not in actual_combinations:
                missing_combinations.append(expected)
        
        if missing_combinations:
            self.logger.warning(f"⚠️ MISSING COMBINATIONS ({len(missing_combinations)}):")
            for snp, sag in missing_combinations:
                self.logger.warning(f"  - SNP={snp}, SAG={sag}")
                
            self.rescue_missing_combinations(missing_combinations)
        else:
            self.logger.info("✅ ALL EXPECTED COMBINATIONS WERE PROCESSED")
        
        failed_combinations = [(r['snp_threshold'], r['sag_threshold']) for r in results if not r['success']]
        if failed_combinations:
            self.logger.warning(f"⚠️ FAILED COMBINATIONS ({len(failed_combinations)}):")
            for snp, sag in failed_combinations:
                failed_result = next(r for r in results if r['snp_threshold'] == snp and r['sag_threshold'] == sag and not r['success'])
                error_msg = failed_result.get('error', 'Unknown error')
                self.logger.warning(f"  - SNP={snp}, SAG={sag}: {error_msg}")
    
    def rescue_missing_combinations(self, missing_combinations):
        
        if not missing_combinations:
            return
        
        self.logger.info(f"\n🔧 RESCUING MISSING COMBINATIONS ({len(missing_combinations)})...")
        
        rescued_count = 0
        
        for i, (snp_threshold, sag_threshold) in enumerate(missing_combinations, 1):
            self.logger.info(f"[RESCUE {i}/{len(missing_combinations)}] Attempting: SNP={snp_threshold}, SAG={sag_threshold}")
            
            try:
                is_snp_changepoint = (snp_threshold == self.initial_snp_threshold)
                sag_changepoint = self.changepoint_cache.get(snp_threshold)
                is_sag_changepoint = (sag_threshold == sag_changepoint) if sag_changepoint else False
                
                result = self.run_single_analysis(
                    snp_threshold, sag_threshold,
                    is_snp_changepoint, is_sag_changepoint,
                    f"RESCUE_{i}", len(missing_combinations)
                )
                
                self.sensitivity_results.append(result)
                
                if result['success']:
                    rescued_count += 1
                    self.logger.info(f"  ✅ RESCUED: SNP={snp_threshold}, SAG={sag_threshold}")
                    
                    if result['fasta_file'] and Path(result['fasta_file']).exists():
                        phylo_dir = self.output_dir / "phylogenetic_msa_true_2d"
                        source_file = Path(result['fasta_file'])
                        
                        priority = ""
                        if result['is_both_changepoint']:
                            priority = "PRIORITY1_BOTH_CP_"
                        elif result['is_snp_changepoint']:
                            priority = "PRIORITY2_SNP_CP_"
                        elif result['is_sag_changepoint']:
                            priority = "PRIORITY2_SAG_CP_"
                        
                        new_filename = f"{priority}msa_snp_{snp_threshold:02d}_sag_{sag_threshold:02d}_sags_{result['sag_count']:03d}_snps_{result['final_snp_count']:04d}.fasta"
                        dest_file = phylo_dir / new_filename
                        
                        try:
                            import shutil
                            shutil.copy2(source_file, dest_file)
                            self.logger.info(f"    📁 Copied rescued file: {new_filename}")
                        except Exception as e:
                            self.logger.warning(f"    ⚠️ Failed to copy rescued file: {e}")
                    
                else:
                    self.logger.error(f"  ❌ RESCUE FAILED: SNP={snp_threshold}, SAG={sag_threshold} → {result.get('error', 'Unknown error')}")
                    
            except Exception as e:
                self.logger.error(f"  💥 RESCUE EXCEPTION: SNP={snp_threshold}, SAG={sag_threshold} → {str(e)}")
        
        if rescued_count > 0:
            self.save_sensitivity_results(self.sensitivity_results)
            self.logger.info(f"✅ RESCUE COMPLETED: {rescued_count}/{len(missing_combinations)} combinations recovered")
        else:
            self.logger.warning(f"❌ RESCUE FAILED: No combinations could be recovered")
    
    def print_summary(self):
        
        if not self.sensitivity_results:
            return
        
        total = len(self.sensitivity_results)
        successful = len([r for r in self.sensitivity_results if r['success']])
        failed = total - successful
        
        self.logger.info(f"\n=== TRUE 2D Sensitivity Analysis Summary ===")
        self.logger.info(f"Total combinations: {total}")
        self.logger.info(f"Successful: {successful}")
        self.logger.info(f"Failed: {failed}")
        self.logger.info(f"Success rate: {successful/total*100:.1f}%")
        
        self.logger.info(f"\n📊 Changepoint Calculations:")
        for snp_threshold in sorted(self.changepoint_cache.keys()):
            sag_changepoint = self.changepoint_cache[snp_threshold]
            marker = " ★" if snp_threshold == self.initial_snp_threshold else ""
            self.logger.info(f"SNP={snp_threshold}{marker} → SAG changepoint={sag_changepoint}")
        
        if successful > 0:
            successful_results = [r for r in self.sensitivity_results if r['success']]
            
            both_changepoints = [r for r in successful_results if r.get('is_both_changepoint', False)]
            snp_changepoints = [r for r in successful_results if r.get('is_snp_changepoint', False)]
            sag_changepoints = [r for r in successful_results if r.get('is_sag_changepoint', False)]
            
            self.logger.info(f"\n🎯 CHANGEPOINT COMBINATION ANALYSIS:")
            self.logger.info(f"Both changepoints (TRUE combinations): {len(both_changepoints)}")
            self.logger.info(f"SNP changepoints only: {len(snp_changepoints)}")
            self.logger.info(f"SAG changepoints only: {len(sag_changepoints)}")
            
            if both_changepoints:
                self.logger.info(f"\n🏆 TRUE CHANGEPOINT COMBINATION RESULTS:")
                for result in both_changepoints:
                    self.logger.info(f"  ✅ SNP={result['snp_threshold']}, SAG={result['sag_threshold']} → "
                                   f"SAGs={result['sag_count']}, SNPs={result['final_snp_count']}, "
                                   f"N_content={result['n_content_mean']:.4f}")
            else:
                self.logger.warning(f"⚠️ NO TRUE CHANGEPOINT COMBINATIONS FOUND IN SUCCESSFUL RESULTS!")
                
                failed_results = [r for r in self.sensitivity_results if not r['success']]
                failed_changepoints = [r for r in failed_results 
                                     if (r['snp_threshold'] == self.initial_snp_threshold and 
                                         r['sag_threshold'] == self.changepoint_cache.get(self.initial_snp_threshold))]
                
                if failed_changepoints:
                    self.logger.error(f"❌ TRUE CHANGEPOINT COMBINATION FAILED:")
                    for result in failed_changepoints:
                        self.logger.error(f"  SNP={result['snp_threshold']}, SAG={result['sag_threshold']} → "
                                        f"Error: {result.get('error', 'Unknown error')}")
                else:
                    self.logger.error(f"🚨 TRUE CHANGEPOINT COMBINATION NOT FOUND AT ALL!")
            
            best_sag = max(successful_results, key=lambda x: x['sag_count'])
            best_snp = max(successful_results, key=lambda x: x['final_snp_count'])
            best_quality = min(successful_results, key=lambda x: x['n_content_mean'])
            
            self.logger.info(f"\n📈 Optimal combinations:")
            self.logger.info(f"Most SAGs: SNP={best_sag['snp_threshold']}, SAG={best_sag['sag_threshold']} ({best_sag['sag_count']} SAGs)")
            self.logger.info(f"Most SNPs: SNP={best_snp['snp_threshold']}, SAG={best_snp['sag_threshold']} ({best_snp['final_snp_count']} SNPs)")
            self.logger.info(f"Best quality: SNP={best_quality['snp_threshold']}, SAG={best_quality['sag_threshold']} (N={best_quality['n_content_mean']:.4f})")
        
        self.logger.info(f"\n🎉 TRUE 2D Analysis completed! Results saved in: {self.output_dir}")
        
        expected_snp = self.initial_snp_threshold
        expected_sag = self.changepoint_cache.get(expected_snp)
        if expected_sag is not None:
            self.logger.info(f"\n🔍 FINAL VERIFICATION:")
            self.logger.info(f"Expected TRUE changepoint combination: SNP={expected_snp}, SAG={expected_sag}")
            
            phylo_dir = self.output_dir / "phylogenetic_msa_true_2d"
            priority1_files = list(phylo_dir.glob("PRIORITY1_BOTH_CP_*.fasta")) if phylo_dir.exists() else []
            
            if priority1_files:
                self.logger.info(f"✅ PRIORITY1 MSA files found: {len(priority1_files)}")
                self.logger.info(f"First file: {priority1_files[0].name}")
            else:
                self.logger.warning(f"⚠️ No PRIORITY1 (both changepoint) MSA files found")
        else:
            self.logger.warning(f"⚠️ Cannot verify: SAG changepoint not calculated for SNP={expected_snp}")
    
    def execute(self):
        
        try:
            self.run_true_2d_sensitivity_analysis()
            
            self.print_summary()
            
        except KeyboardInterrupt:
            self.logger.info("TRUE 2D Analysis interrupted by user")
            sys.exit(1)
        except Exception as e:
            self.logger.error(f"TRUE 2D Analysis failed: {str(e)}")
            sys.exit(1)

def main():
    
    parser = argparse.ArgumentParser(description='TRUE 2D Sensitivity Analysis for SAG Pipeline')
    
    parser.add_argument('--input_matrix', required=True, type=Path, help='Input SNP matrix file')
    parser.add_argument('--snp_threshold_file', required=True, type=Path, help='Initial SNP threshold from changepoint analysis')
    parser.add_argument('--output_dir', required=True, type=Path, help='Output directory')
    
    parser.add_argument('--snp_range', type=int, default=3, 
                       help='Range around SNP changepoint (default: ±3)')
    parser.add_argument('--sag_range', type=int, default=3, 
                       help='Range around each calculated SAG changepoint (default: ±3)')
    
    parser.add_argument('--quiet', action='store_true', help='Suppress verbose output')
    
    args = parser.parse_args()
    
    print(f"DEBUG: Received arguments:")
    print(f"  input_matrix: {args.input_matrix}")
    print(f"  snp_threshold_file: {args.snp_threshold_file}")
    print(f"  output_dir: {args.output_dir}")
    print(f"  snp_range: {args.snp_range}")
    print(f"  sag_range: {args.sag_range}")
    print(f"  quiet: {args.quiet}")
    sys.stdout.flush()
    
    sensitivity_analysis = True2DSensitivityAnalysis(args)
    sensitivity_analysis.execute()

if __name__ == "__main__":
    main()