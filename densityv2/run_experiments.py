import os
import sys
import pandas as pd
import numpy as np
import time
import re
from pathlib import Path
from typing import Dict, List, Any, Tuple

# Fix import paths when running as main script
if __name__ == "__main__":
    # Add parent directory to Python path
    current_dir = Path(__file__).parent
    parent_dir = current_dir.parent
    sys.path.insert(0, str(parent_dir))
    
    # Now import with absolute paths
    from ilphaplo.get_data import get_data
    from ilphaplo.pipeline import run_pipeline
else:
    # Relative imports for when imported as module
    from .get_data import get_data
    from .pipeline import run_pipeline

def validate_density_constraint(matrix: np.ndarray, zone_rows: List[int], zone_cols: List[int], error_rate: float = 0.025) -> Tuple[bool, float, Dict[str, Any]]:
    """
    Validate density constraint for a detected zone.
    
    Constraint: Σ_{i,j} x_ij × (1-M[i,j]) ≤ error_rate × Σ_{i,j} x_ij
    
    Parameters
    ----------
    matrix : np.ndarray
        Binary matrix (0,1)
    zone_rows : List[int]
        Row indices of the zone
    zone_cols : List[int]
        Column indices of the zone
    error_rate : float
        Maximum allowed error rate
        
    Returns
    -------
    Tuple[bool, float, Dict[str, Any]]
        (constraint_satisfied, actual_error_rate, details)
    """
    if not zone_rows or not zone_cols:
        return True, 0.0, {"reason": "empty_zone"}
    
    # Extract zone submatrix
    zone_matrix = matrix[np.ix_(zone_rows, zone_cols)]
    
    # Calculate constraint terms
    total_selected_cells = zone_matrix.size  # All cells in zone are selected (x_ij = 1)
    zeros_in_zone = np.sum(zone_matrix == 0)  # Count of 0s in selected region
    
    # Calculate actual error rate
    actual_error_rate = zeros_in_zone / total_selected_cells if total_selected_cells > 0 else 0.0
    
    # Check constraint: zeros_in_zone ≤ error_rate × total_selected_cells
    constraint_satisfied = actual_error_rate <= error_rate
    
    details = {
        "zone_size": (len(zone_rows), len(zone_cols)),
        "total_cells": total_selected_cells,
        "zeros_count": zeros_in_zone,
        "ones_count": total_selected_cells - zeros_in_zone,
        "actual_error_rate": actual_error_rate,
        "allowed_error_rate": error_rate,
        "violation_amount": max(0, actual_error_rate - error_rate)
    }
    
    return constraint_satisfied, actual_error_rate, details

def load_target_files() -> List[str]:
    """Load list of files to process from filepaths_with_ilp_calls.txt"""
    target_files_path = Path(__file__).parent / "filepaths_with_ilp_calls.txt"
    
    if not target_files_path.exists():
        print(f"Warning: {target_files_path} not found")
        return []
    
    with open(target_files_path, 'r') as f:
        files = [line.strip() for line in f if line.strip()]
    
    print(f"Loaded {len(files)} target files from filepaths_with_ilp_calls.txt")
    return files

def find_matrice_directory(start_path: str = ".") -> Path:
    """Find the matrice directory by searching up the directory tree."""
    current = Path(start_path).resolve()
    
    # Check if a specific path was provided as command line argument
    if len(sys.argv) > 1:
        provided_path = Path(sys.argv[1])
        if provided_path.exists():
            print(f"Using provided matrix directory: {provided_path}")
            return provided_path
    
    # Look for matrices directory specifically
    possible_matrices_paths = [
        current / "matrices",           # ./matrices
        current.parent / "matrices",    # ../matrices
        Path("matrices"),               # relative matrices
    ]
    
    for matrices_path in possible_matrices_paths:
        if matrices_path.exists() and matrices_path.is_dir():
            print(f"Found matrices directory at: {matrices_path}")
            return matrices_path
    
    print("No matrices directory found.")
    return None

def extract_haplotype_count(filepath: str) -> int:
    """Extract haplotype count from file path or filename."""
    parts = Path(filepath).parts
    
    # First try to extract from directory structure (numbered directories)
    for part in parts:
        if part.isdigit():
            hap_count = int(part)
            if 2 <= hap_count <= 10:
                return hap_count
    
    # Try to extract from filename
    filename = Path(filepath).name
    import re
    
    # Look for various patterns in filename
    patterns = [
        r'(?:matrix|hap|hapl)_?(\d+)',
        r'(\d+)_?(?:hap|hapl)',
        r'(\d+)haplotype',
        r'h(\d+)',
        r'_(\d+)_',  # Pattern like matrix_2_1.csv
        r'^(\d+)_',  # Pattern like 2_matrix.csv
    ]
    
    for pattern in patterns:
        match = re.search(pattern, filename.lower())
        if match:
            hap_count = int(match.group(1))
            if 2 <= hap_count <= 10:
                return hap_count
    
    # Default to 0 if no haplotype count found (will be handled gracefully)
    return 0

def process_all_matrices(base_dir: str = None) -> List[Dict[str, Any]]:
    """
    Process only the CSV matrices listed in filepaths_with_ilp_calls.txt.
    """
    
    results = []
    violation_log = []
    
    # Load target files list
    target_files = load_target_files()
    if not target_files:
        print("No target files specified")
        return results
    
    # Find matrice directory
    if base_dir is None:
        base_path = find_matrice_directory()
        if base_path is None:
            print("Error: Could not find matrices directory")
            return results
    else:
        base_path = Path(base_dir)
        if not base_path.exists():
            print(f"Specified directory {base_path} not found")
            return results
    
    print(f"Using matrices directory: {base_path.resolve()}")
    
    # Process only target files
    processed_count = 0
    for i, target_file in enumerate(target_files):
        csv_file = base_path / target_file
        
        if not csv_file.exists():
            print(f"Warning: Target file not found: {csv_file}")
            continue
            
        try:
            print(f"Processing {i+1}/{len(target_files)}: {target_file}")
            
            # Extract haplotype count from directory structure
            haplotype_count = extract_haplotype_count(str(csv_file))
            
            # Load and process matrix
            start_time = time.time()
            rows_data, cols_data, edges_0, edges_1, row_names, col_names, df, comp_df = get_data(str(csv_file))
            load_time = time.time() - start_time
            
            # Convert to binary matrix (ensure 0/1 values)
            binary_matrix = df.values.astype(int)
            
            # Run pipeline and collect statistics
            pipeline_results = run_pipeline(
                binary_matrix,
                min_col_quality=3,
                min_row_quality=5,
                error_rate=0.025,
                filename=csv_file.name,
                haplotype_count=haplotype_count
            )
            
            # Validate density constraints for detected zones
            zones_found = pipeline_results.get('zones_detected', [])
            constraint_violations = []
            
            for zone_idx, zone in enumerate(zones_found):
                if len(zone) >= 3:  # zone format: (reads1, reads0, cols)
                    reads1, reads0, cols = zone
                    
                    # Validate density for both read groups
                    if reads1 and cols:
                        valid1, error_rate1, details1 = validate_density_constraint(
                            binary_matrix, reads1, cols, error_rate=0.025
                        )
                        if not valid1:
                            constraint_violations.append({
                                'zone_index': zone_idx,
                                'group': 'reads1',
                                'details': details1
                            })
                    
                    if reads0 and cols:
                        valid0, error_rate0, details0 = validate_density_constraint(
                            binary_matrix, reads0, cols, error_rate=0.025
                        )
                        if not valid0:
                            constraint_violations.append({
                                'zone_index': zone_idx,
                                'group': 'reads0',
                                'details': details0
                            })
            
            # Log violations if any
            if constraint_violations:
                violation_entry = {
                    'matrix_name': target_file,
                    'total_zones_found': len(zones_found),
                    'violations_count': len(constraint_violations),
                    'violating_zones': constraint_violations
                }
                violation_log.append(violation_entry)
                print(f"  -> CONSTRAINT VIOLATIONS: {len(constraint_violations)} violations found")
            
            # Compile comprehensive statistics
            stats = {
                # File information
                'filename': csv_file.name,
                'filepath': target_file,
                'haplotype_count': haplotype_count,
                'load_time': load_time,
                
                # Matrix characteristics
                'matrix_rows': binary_matrix.shape[0],
                'matrix_cols': binary_matrix.shape[1],
                'matrix_size': binary_matrix.size,
                'matrix_density': np.sum(binary_matrix) / binary_matrix.size,
                'ones_count': int(np.sum(binary_matrix)),
                'zeros_count': int(binary_matrix.size - np.sum(binary_matrix)),
                
                # Row statistics
                'row_min_sum': int(np.min(np.sum(binary_matrix, axis=1))),
                'row_max_sum': int(np.max(np.sum(binary_matrix, axis=1))),
                'row_mean_sum': float(np.mean(np.sum(binary_matrix, axis=1))),
                'row_std_sum': float(np.std(np.sum(binary_matrix, axis=1))),
                
                # Column statistics  
                'col_min_sum': int(np.min(np.sum(binary_matrix, axis=0))),
                'col_max_sum': int(np.max(np.sum(binary_matrix, axis=0))),
                'col_mean_sum': float(np.mean(np.sum(binary_matrix, axis=0))),
                'col_std_sum': float(np.std(np.sum(binary_matrix, axis=0))),
                
                # Pipeline performance
                'ilp_calls_total': pipeline_results.get('ilp_calls_total', 0),
                'preprocessing_time': pipeline_results.get('preprocessing_time', 0),
                'clustering_time': pipeline_results.get('clustering_time', 0),
                'total_time': pipeline_results.get('total_time', 0),
                'ilp_time_total': pipeline_results.get('ilp_time_total', 0),
                'ilp_time_ratio': pipeline_results.get('ilp_time_ratio', 0),
                'avg_ilp_time': pipeline_results.get('avg_ilp_time', 0),
                
                # Results
                'patterns_found': pipeline_results.get('patterns_found', 0),
                'regions_found': pipeline_results.get('regions_found', 0),
                'clustering_steps': pipeline_results.get('clustering_steps', 0),
                'regions_processed': pipeline_results.get('regions_processed', 0),
                'matrix_operations': pipeline_results.get('matrix_operations', 0),
                
                # Region characteristics
                'avg_region_size': pipeline_results.get('avg_region_size', 0),
                'largest_region': pipeline_results.get('largest_region', 0),
                'smallest_region': pipeline_results.get('smallest_region', 0),
                
                # Efficiency metrics
                'ilp_calls_per_pattern': pipeline_results.get('ilp_calls_total', 0) / max(1, pipeline_results.get('patterns_found', 1)),
                'time_per_ilp_call': pipeline_results.get('ilp_time_total', 0) / max(1, pipeline_results.get('ilp_calls_total', 1)),
                'patterns_per_second': pipeline_results.get('patterns_found', 0) / max(1, pipeline_results.get('total_time', 1)),
                'matrix_complexity': binary_matrix.shape[0] * binary_matrix.shape[1] * pipeline_results.get('matrix_density', 0),
                
                # Density validation results
                'constraint_violations': len(constraint_violations),
                'zones_validated': len(zones_found),
                'constraint_satisfied': len(constraint_violations) == 0,
            }
            
            # Add solver status counts if available
            solver_stats = pipeline_results.get('solver_status_counts', {})
            stats.update({
                'solver_success_count': solver_stats.get('success', 0),
                'solver_failed_count': solver_stats.get('failed', 0),
                'solver_error_count': solver_stats.get('error', 0),
            })
            
            results.append(stats)
            processed_count += 1
            print(f"  -> Completed: {stats['ilp_calls_total']} ILP calls, {stats['patterns_found']} patterns, {stats['total_time']:.3f}s")
            
        except Exception as e:
            print(f"  -> Error processing {target_file}: {str(e)}")
            # Add error entry
            error_stats = {
                'filename': csv_file.name,
                'filepath': str(csv_file.relative_to(base_path)),
                'haplotype_count': extract_haplotype_count(str(csv_file)),
                'error': str(e),
                'matrix_rows': 0,
                'matrix_cols': 0,
                'ilp_calls_total': 0,
                'total_time': 0,
                'patterns_found': 0
            }
            results.append(error_stats)
    
    # Save constraint violations log
    if violation_log:
        save_violations_log(violation_log)
    
    print(f"Processed {processed_count} target files successfully")
    return results

def save_results_to_csv(results: List[Dict[str, Any]], output_file: str = "experiment_results.csv"):
    """Save experiment results to CSV file."""
    
    if not results:
        print("No results to save")
        return
    
    # Convert to DataFrame
    df = pd.DataFrame(results)
    
    # Sort by haplotype count and matrix size
    df = df.sort_values(['haplotype_count', 'matrix_size'], ascending=[True, True])
    
    # Save to CSV
    df.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")
    
    # Print summary statistics with correct interpretation
    print("\n=== EXPERIMENT SUMMARY ===")
    print(f"Total matrices processed: {len(df)}")
    print(f"Efficient runs (no ILP needed): {len(df[df['ilp_calls_total'] == 0])}")
    print(f"Complex runs (ILP required): {len(df[df['ilp_calls_total'] > 0])}")
    print(f"Algorithm efficiency rate: {len(df[df['ilp_calls_total'] == 0]) / len(df) * 100:.1f}%")
    
    if len(df) > 0:
        all_runs_df = df
        print(f"\nHaplotype counts: {sorted(all_runs_df['haplotype_count'].unique())}")
        print(f"Matrix sizes range: {all_runs_df['matrix_size'].min()} - {all_runs_df['matrix_size'].max()}")
        print(f"Total ILP calls across all runs: {all_runs_df['ilp_calls_total'].sum()}")
        print(f"Total processing time: {all_runs_df['total_time'].sum():.2f}s")
        print(f"Average ILP calls per matrix: {all_runs_df['ilp_calls_total'].mean():.1f}")
        print(f"Average patterns found per matrix: {all_runs_df['patterns_found'].mean():.1f}")
        
        # Group by haplotype count
        print(f"\n=== BY HAPLOTYPE COUNT ===")
        for hap_count in sorted(all_runs_df['haplotype_count'].unique()):
            hap_df = all_runs_df[all_runs_df['haplotype_count'] == hap_count]
            efficient_count = len(hap_df[hap_df['ilp_calls_total'] == 0])
            efficiency_rate = efficient_count / len(hap_df) * 100
            print(f"Haplotypes {hap_count}: {len(hap_df)} matrices, "
                  f"avg ILP calls: {hap_df['ilp_calls_total'].mean():.1f}, "
                  f"avg time: {hap_df['total_time'].mean():.2f}s, "
                  f"efficiency: {efficiency_rate:.1f}%")

def save_violations_log(violations: List[Dict[str, Any]], output_file: str = None):
    """Save constraint violations to CSV log file."""
    if output_file is None:
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        output_file = f"constraint_violations_{timestamp}.csv"
    
    # Flatten violations for CSV format
    flattened_violations = []
    for violation in violations:
        matrix_name = violation['matrix_name']
        total_zones = violation['total_zones_found']
        
        for v in violation['violating_zones']:
            flattened_violations.append({
                'matrix_name': matrix_name,
                'total_zones_found': total_zones,
                'violating_zone_index': v['zone_index'],
                'violating_group': v['group'],
                'zone_size_rows': v['details']['zone_size'][0],
                'zone_size_cols': v['details']['zone_size'][1],
                'total_cells': v['details']['total_cells'],
                'zeros_count': v['details']['zeros_count'],
                'actual_error_rate': v['details']['actual_error_rate'],
                'allowed_error_rate': v['details']['allowed_error_rate'],
                'violation_amount': v['details']['violation_amount']
            })
    
    if flattened_violations:
        df = pd.DataFrame(flattened_violations)
        df.to_csv(output_file, index=False)
        print(f"Constraint violations logged to {output_file}")
        print(f"Total violations: {len(flattened_violations)} across {len(set(df['matrix_name']))} matrices")
    else:
        print("No constraint violations found")

def run_experiments():
    """Main function to run all experiments."""
    
    print("Starting matrix processing experiments...")
    print(f"Working directory: {os.getcwd()}")
    
    # Process all matrices - don't specify base_dir to trigger auto-discovery
    results = process_all_matrices()
    
    # Save results
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    output_file = f"experiment_results_{timestamp}.csv"
    save_results_to_csv(results, output_file)
    
    return results

if __name__ == "__main__":
    results = run_experiments()
