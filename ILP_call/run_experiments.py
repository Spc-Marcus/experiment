import os
import sys
import pandas as pd
import numpy as np
import time
import re
from pathlib import Path
from typing import Dict, List, Any
import argparse

# Fix import paths when running as main script
if __name__ == "__main__":
    # Add parent directories to Python path
    current_dir = Path(__file__).parent
    parent_dir = current_dir.parent
    experiment_dir = parent_dir
    sys.path.insert(0, str(experiment_dir))
    sys.path.insert(0, str(current_dir))
    
    # Import from current directory
    from get_data import get_data
    from pipeline import run_pipeline
else:
    # Relative imports for when imported as module
    from .get_data import get_data
    from .pipeline import run_pipeline

def find_matrice_directory(matrix_dir: str) -> Path:
    """Validate and return the provided matrix directory."""
    
    if not matrix_dir:
        print("Error: No matrix directory provided")
        return None
        
    provided_path = Path(matrix_dir)
    if not provided_path.exists():
        print(f"Error: Matrix directory {provided_path} does not exist")
        return None
        
    if not provided_path.is_dir():
        print(f"Error: {provided_path} is not a directory")
        return None
        
    # Check for haplotype subdirectories (numbered directories)
    numbered_subdirs = [d for d in provided_path.iterdir() 
                       if d.is_dir() and d.name.isdigit()]
    
    if not numbered_subdirs:
        print(f"Warning: No numbered subdirectories found in {provided_path}")
        print("Expected structure: matrix_dir/1/, matrix_dir/2/, etc.")
        
    # Count CSV files
    total_csv_files = sum(len(list(subdir.glob("*.csv"))) for subdir in numbered_subdirs)
    
    if total_csv_files == 0:
        print(f"Error: No CSV files found in numbered subdirectories of {provided_path}")
        return None
        
    print(f"Found matrix directory: {provided_path}")
    print(f"Haplotype subdirectories: {sorted([d.name for d in numbered_subdirs])}")
    print(f"Total CSV files: {total_csv_files}")
    
    return provided_path

def extract_haplotype_count(filepath: str) -> int:
    """Extract haplotype count from file path directory structure."""
    parts = Path(filepath).parts
    
    # Extract from directory structure (numbered directories)
    for part in parts:
        if part.isdigit():
            hap_count = int(part)
            if 1 <= hap_count <= 10:  # Support 1-10 haplotypes
                return hap_count
    
    # Default to 0 if no haplotype count found
    return 0

def process_all_matrices(matrix_dir: str, output_dir: str) -> List[Dict[str, Any]]:
    """
    Process all CSV matrices and collect statistics.
    
    Parameters
    ----------
    matrix_dir : str
        Directory containing haplotype subdirectories (1/, 2/, 3/, etc.)
    output_dir : str
        Directory where to save results
    """
    
    results = []
    
    # Validate matrix directory
    base_path = find_matrice_directory(matrix_dir)
    if base_path is None:
        return results
    
    print(f"Processing matrices from: {base_path.resolve()}")
    print(f"Results will be saved to: {Path(output_dir).resolve()}")
    
    # Find all CSV files in haplotype subdirectories
    csv_files = []
    for hap_dir in sorted(base_path.iterdir()):
        if hap_dir.is_dir() and hap_dir.name.isdigit():
            hap_csv_files = list(hap_dir.glob("*.csv"))
            csv_files.extend(hap_csv_files)
            print(f"Haplotype {hap_dir.name}: {len(hap_csv_files)} CSV files")
    
    print(f"\nTotal CSV files to process: {len(csv_files)}")
    
    if len(csv_files) == 0:
        print("No CSV files found in haplotype subdirectories")
        return results
    
    for i, csv_file in enumerate(csv_files):
        try:
            print(f"\nProcessing {i+1}/{len(csv_files)}: {csv_file.relative_to(base_path)}")
            
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
            
            # Compile comprehensive statistics
            stats = {
                # File information
                'filename': csv_file.name,
                'filepath': str(csv_file.relative_to(base_path)),
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
            }
            
            # Add solver status counts if available
            solver_stats = pipeline_results.get('solver_status_counts', {})
            stats.update({
                'solver_success_count': solver_stats.get('success', 0),
                'solver_failed_count': solver_stats.get('failed', 0),
                'solver_error_count': solver_stats.get('error', 0),
            })
            
            results.append(stats)
            print(f"  -> Completed: {stats['ilp_calls_total']} ILP calls, {stats['patterns_found']} patterns, {stats['total_time']:.3f}s")
            
        except Exception as e:
            print(f"  -> Error processing {csv_file}: {str(e)}")
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
    
    return results

def save_results_to_csv(results: List[Dict[str, Any]], output_dir: str, output_file: str = None):
    """Save experiment results to CSV file in specified directory."""
    
    if not results:
        print("No results to save")
        return None
    
    # Ensure output directory exists
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Generate filename if not provided
    if output_file is None:
        output_file = "ilp_call_stats.csv"
    
    # Full path for output file
    full_output_path = output_path / output_file
    
    # Convert to DataFrame
    df = pd.DataFrame(results)
    
    # Sort by haplotype count and matrix size
    df = df.sort_values(['haplotype_count', 'matrix_size'], ascending=[True, True])
    
    # Save to CSV
    df.to_csv(full_output_path, index=False)
    print(f"\nResults saved to: {full_output_path}")
    
    # Print summary statistics
    print("\n=== EXPERIMENT SUMMARY ===")
    print(f"Total matrices processed: {len(df)}")
    print(f"Efficient runs (no ILP needed): {len(df[df['ilp_calls_total'] == 0])}")
    print(f"Complex runs (ILP required): {len(df[df['ilp_calls_total'] > 0])}")
    print(f"Algorithm efficiency rate: {len(df[df['ilp_calls_total'] == 0]) / len(df) * 100:.1f}%")
    
    if len(df) > 0:
        print(f"\nHaplotype counts: {sorted(df['haplotype_count'].unique())}")
        print(f"Matrix sizes range: {df['matrix_size'].min()} - {df['matrix_size'].max()}")
        print(f"Total ILP calls across all runs: {df['ilp_calls_total'].sum()}")
        print(f"Total processing time: {df['total_time'].sum():.2f}s")
        print(f"Average ILP calls per matrix: {df['ilp_calls_total'].mean():.1f}")
        print(f"Average patterns found per matrix: {df['patterns_found'].mean():.1f}")
        
        # Group by haplotype count
        print(f"\n=== BY HAPLOTYPE COUNT ===")
        for hap_count in sorted(df['haplotype_count'].unique()):
            hap_df = df[df['haplotype_count'] == hap_count]
            efficient_count = len(hap_df[hap_df['ilp_calls_total'] == 0])
            efficiency_rate = efficient_count / len(hap_df) * 100
            print(f"Haplotypes {hap_count}: {len(hap_df)} matrices, "
                  f"avg ILP calls: {hap_df['ilp_calls_total'].mean():.1f}, "
                  f"avg time: {hap_df['total_time'].mean():.2f}s, "
                  f"efficiency: {efficiency_rate:.1f}%")
    
    return full_output_path

def run_experiments(matrix_dir: str, output_dir: str):
    """Main function to run all experiments with specified directories."""
    
    print("Starting ILP Haplotype Analysis Experiments")
    print("=" * 50)
    print(f"Matrix directory: {matrix_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Working directory: {os.getcwd()}")
    
    # Process all matrices
    results = process_all_matrices(matrix_dir, output_dir)
    
    # Save results in specified output directory
    results_file = save_results_to_csv(results, output_dir)
    
    if results_file:
        print(f"\n🎉 Experiment completed!")
        print(f"📊 Results file: {results_file}")
        print(f"📁 Output directory: {Path(output_dir).resolve()}")
        return results_file
    else:
        print("\n❌ Experiment failed - no results generated")
        return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Run ILP haplotype analysis experiments',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_experiments.py -m /path/to/matrices_bin -o /path/to/experiment
  python run_experiments.py --matrix-dir ./data/matrices_bin --output-dir ./results
        """
    )
    parser.add_argument('--matrix-dir', '-m', type=str, required=True,
                       help='Directory containing haplotype subdirectories (1/, 2/, 3/, etc.)')
    parser.add_argument('--output-dir', '-o', type=str, required=True,
                       help='Directory to save analysis results')
    
    args = parser.parse_args()
    
    results_file = run_experiments(args.matrix_dir, args.output_dir)
    
    if results_file:
        sys.exit(0)
    else:
        sys.exit(1)
    # Sort by haplotype count and matrix size
    df = df.sort_values(['haplotype_count', 'matrix_size'], ascending=[True, True])
    
    # Save to CSV
    df.to_csv(full_output_path, index=False)
    print(f"Results saved to {full_output_path}")
    
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

def run_experiments(matrix_dir: str = None, output_dir: str = "."):
    """Main function to run all experiments with specified directories."""
    
    print("Starting matrix processing experiments...")
    print(f"Working directory: {os.getcwd()}")
    print(f"Matrix directory: {matrix_dir if matrix_dir else 'auto-detect'}")
    print(f"Output directory: {output_dir}")
    
    # Process all matrices
    results = process_all_matrices(matrix_dir, output_dir)
    
    # Save results in specified output directory
    save_results_to_csv(results, output_dir)
    
    return results

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run ILP haplotype analysis experiments')
    parser.add_argument('--matrix-dir', '-m', type=str, 
                       help='Directory containing matrix subdirectories')
    parser.add_argument('--output-dir', '-o', type=str, default=".",
                       help='Directory to save analysis results (default: current directory)')
    
    args = parser.parse_args()
    
    results = run_experiments(args.matrix_dir, args.output_dir)
    results = run_experiments(args.matrix_dir, args.output_dir)
