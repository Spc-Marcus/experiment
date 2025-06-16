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

def find_matrice_directory(start_path: str = ".") -> Path:
    """Find the matrice directory by searching up the directory tree."""
    current = Path(start_path).resolve()
    
    # Check if a specific path was provided as command line argument
    if len(sys.argv) > 1:
        provided_path = Path(sys.argv[1])
        if provided_path.exists():
            print(f"Using provided matrix directory: {provided_path}")
            return provided_path
    
    # Look for KNN output directory with numbered subdirectories
    possible_matrices_paths = [
        current / "matrices",           # ./matrices
        current / "../matrices",        # ../matrices  
        current.parent / "matrices",    # parent/matrices
        Path("matrices"),               # relative matrices
        Path("../matrices"),            # relative ../matrices
    ]
    
    for matrices_path in possible_matrices_paths:
        if matrices_path.exists() and matrices_path.is_dir():
            # Check for numbered subdirectories (haplotype structure)
            numbered_subdirs = [d for d in matrices_path.iterdir() 
                              if d.is_dir() and d.name.isdigit()]
            
            # Check for CSV files in numbered subdirectories
            total_csv_files = 0
            for subdir in numbered_subdirs:
                csv_files_in_subdir = list(subdir.glob("*.csv"))
                total_csv_files += len(csv_files_in_subdir)
            
            # Also check for direct CSV files
            direct_csv_files = list(matrices_path.glob("*.csv"))
            total_csv_files += len(direct_csv_files)
            
            if total_csv_files > 0:
                print(f"Found KNN matrices directory at: {matrices_path}")
                if numbered_subdirs:
                    print(f"Contains numbered subdirectories: {sorted([d.name for d in numbered_subdirs])}")
                    print(f"Total CSV files in subdirectories: {total_csv_files - len(direct_csv_files)}")
                if direct_csv_files:
                    print(f"Direct CSV files: {len(direct_csv_files)}")
                print(f"Total CSV files: {total_csv_files}")
                return matrices_path
    
    # Look in current directory and parent directories for organized data
    for path in [current] + list(current.parents)[:2]:  # Don't go too far up
        # Check multiple possible directory names
        for dir_name in ["matrice", "matrices", "data", "datasets", "matrix_data", "haplotype_data"]:
            matrice_path = path / dir_name
            if matrice_path.exists() and matrice_path.is_dir():
                # Check for CSV files directly or in numbered subdirectories
                csv_files = list(matrice_path.glob("*.csv"))
                subdirs = [d for d in matrice_path.iterdir() if d.is_dir() and d.name.isdigit()]
                
                if csv_files:
                    print(f"Found matrix directory at: {matrice_path}")
                    print(f"Contains {len(csv_files)} CSV files directly")
                    return matrice_path
                elif subdirs:
                    # Check if subdirectories contain CSV files
                    subdir_csv_count = sum(len(list(d.glob("*.csv"))) for d in subdirs)
                    if subdir_csv_count > 0:
                        print(f"Found matrix directory at: {matrice_path}")
                        print(f"Contains haplotype subdirectories: {sorted([d.name for d in subdirs])}")
                        print(f"Total CSV files in subdirectories: {subdir_csv_count}")
                        return matrice_path
    
    # If not found, return None instead of creating
    print("No existing matrix directory found.")
    print("Searched in the following locations:")
    for path in possible_matrices_paths:
        print(f"  - {path}")
    print("")
    print("Please ensure your data is organized as:")
    print("matrices/ (KNN output with haplotype structure)")
    print("├── 2/")
    print("│   ├── matrix1.csv")
    print("│   └── matrix2.csv")
    print("├── 3/")
    print("│   └── matrix3.csv")
    print("└── ...")
    print("OR")
    print("matrices/ (direct CSV files)")
    print("├── matrix1.csv")
    print("├── matrix2.csv")
    print("└── ...")
    
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

def process_all_matrices(base_dir: str = None, solver_type: str = "gurobi") -> List[Dict[str, Any]]:  # NEW: Add solver parameter
    """
    Process all CSV matrices and collect statistics.
    
    Parameters
    ----------
    base_dir : str, optional
        Base directory containing matrices
    solver_type : str, optional
        Solver to use: "gurobi", "max_ones", or "max_ones_comp"
    """
    
    results = []
    
    # Find matrice directory
    if base_dir is None:
        base_path = find_matrice_directory()
        if base_path is None:
            print("Error: Could not find matrice directory with existing data")
            return results
    else:
        base_path = Path(base_dir)
        if not base_path.exists():
            print(f"Specified directory {base_path} not found")
            return results
    
    print(f"Using matrice directory: {base_path.resolve()}")
    print(f"Using solver: {solver_type}")
    
    # Find all CSV files in directory and subdirectories
    csv_files = list(base_path.rglob("*.csv"))
    
    print(f"Found {len(csv_files)} CSV files to process")
    
    if len(csv_files) == 0:
        print("No CSV files found in matrice directory")
        print("Please check that your data files are in the correct location:")
        print(f"  {base_path.resolve()}")
        return results
    
    for i, csv_file in enumerate(csv_files):
        try:
            print(f"Processing {i+1}/{len(csv_files)}: {csv_file} with {solver_type}")
            
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
                haplotype_count=haplotype_count,
                solver_type=solver_type  # NEW: Pass solver type
            )
            
            # Compile comprehensive statistics
            stats = {
                # File information
                'filename': csv_file.name,
                'filepath': str(csv_file.relative_to(base_path)),
                'haplotype_count': haplotype_count,
                'load_time': load_time,
                'solver_type': solver_type,  # NEW: Track solver used
                
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
                'solver_type': solver_type,  # NEW: Track solver even for errors
                'error': str(e),
                'matrix_rows': 0,
                'matrix_cols': 0,
                'ilp_calls_total': 0,
                'total_time': 0,
                'patterns_found': 0
            }
            results.append(error_stats)
    
    return results

def save_results_to_csv(results: List[Dict[str, Any]], output_file: str = "experiment_results.csv"):
    """Save experiment results to CSV file."""
    
    if not results:
        print("No results to save")
        return
    
    # Convert to DataFrame
    df = pd.DataFrame(results)
    
    # Sort by solver type, haplotype count and matrix size
    df = df.sort_values(['solver_type', 'haplotype_count', 'matrix_size'], ascending=[True, True, True])
    
    # Save to CSV
    df.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")
    
    # Print summary statistics with solver breakdown
    print("\n=== EXPERIMENT SUMMARY ===")
    print(f"Total matrices processed: {len(df)}")
    
    # Group by solver type
    solver_summary = df.groupby('solver_type').agg({
        'ilp_calls_total': ['count', lambda x: (x == 0).sum(), lambda x: (x > 0).sum()],
        'total_time': ['mean', 'sum'],
        'patterns_found': 'mean'
    }).round(3)
    
    print(f"\n=== SOLVER COMPARISON ===")
    for solver in df['solver_type'].unique():
        solver_df = df[df['solver_type'] == solver]
        efficient_count = len(solver_df[solver_df['ilp_calls_total'] == 0])
        efficiency_rate = efficient_count / len(solver_df) * 100
        avg_time = solver_df['total_time'].mean()
        avg_ilp_calls = solver_df['ilp_calls_total'].mean()
        
        print(f"{solver}: {len(solver_df)} runs, "
              f"efficiency: {efficiency_rate:.1f}%, "
              f"avg ILP calls: {avg_ilp_calls:.1f}, "
              f"avg time: {avg_time:.3f}s")
    
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

def run_experiments(solver_type: str = "gurobi"):  # NEW: Add solver parameter
    """Main function to run all experiments."""
    
    print(f"Starting matrix processing experiments with {solver_type} solver...")
    print(f"Working directory: {os.getcwd()}")
    
    # Process all matrices
    results = process_all_matrices(solver_type=solver_type)
    
    # Save results
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    output_file = f"experiment_results_{solver_type}_{timestamp}.csv"
    save_results_to_csv(results, output_file)
    
    return results

if __name__ == "__main__":
    # NEW: Add command line argument parsing
    parser = argparse.ArgumentParser(description='Run ILP haplotype experiments')
    parser.add_argument('--solver', choices=['gurobi', 'max_ones', 'max_ones_comp'], 
                       default='gurobi', help='Solver to use (default: gurobi)')
    parser.add_argument('base_dir', nargs='?', help='Base directory containing matrices (optional)')
    
    args = parser.parse_args()
    
    print(f"Using solver: {args.solver}")
    
    if args.base_dir:
        results = process_all_matrices(args.base_dir, args.solver)
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        output_file = f"experiment_results_{args.solver}_{timestamp}.csv"
        save_results_to_csv(results, output_file)
    else:
        results = run_experiments(args.solver)
