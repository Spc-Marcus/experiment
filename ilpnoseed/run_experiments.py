import os
import sys
import pandas as pd
import numpy as np
import time
import re
from pathlib import Path
from typing import Dict, List, Any

# Fix import paths when running as main script
if __name__ == "__main__":
    # Add parent directory to Python path
    current_dir = Path(__file__).parent
    parent_dir = current_dir.parent
    sys.path.insert(0, str(parent_dir))
    
    # Import from local ilpnoseed module instead of ilphaplo
    from ilpnoseed.get_data import get_data
    from ilpnoseed.pipeline import run_pipeline
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

def print_comparison_summary(comparison_summary: List[Dict[str, Any]]):
    """Print detailed comparison summary with new metrics."""
    
    if not comparison_summary:
        return
    
    print(f"\n🔍 COMPARISON SUMMARY:")
    print(f"Total files compared: {len(comparison_summary)}")
    
    # Basic success rates
    both_success = sum(1 for c in comparison_summary if c['results_identical'] is not None)
    identical_count = sum(1 for c in comparison_summary if c.get('results_identical', False))
    equivalent_count = sum(1 for c in comparison_summary if c.get('results_equivalent', False))
    
    print(f"Both versions succeeded: {both_success}/{len(comparison_summary)}")
    print(f"Identical results: {identical_count}/{both_success}")
    print(f"Equivalent results: {equivalent_count}/{both_success}")
    
    if both_success > 0:
        avg_noseed_time = np.mean([c['noseed_time'] for c in comparison_summary if c.get('noseed_time', 0) > 0])
        avg_seed_time = np.mean([c['seed_time'] for c in comparison_summary if c.get('seed_time', 0) > 0])
        print(f"Average time - NOSEED: {avg_noseed_time:.3f}s, SEED: {avg_seed_time:.3f}s")
    
    # NEW: Detailed cluster metrics summary
    detailed_comparisons = [c for c in comparison_summary if 'detailed_metrics' in c]
    if detailed_comparisons:
        print(f"\n=== DETAILED CLUSTER ANALYSIS ===")
        
        total_exchanged = sum(c['detailed_metrics'].get('reads_exchanged', 0) for c in detailed_comparisons)
        total_unique_noseed = sum(c['detailed_metrics'].get('reads_unique_noseed', 0) for c in detailed_comparisons)
        total_unique_seed = sum(c['detailed_metrics'].get('reads_unique_seed', 0) for c in detailed_comparisons)
        total_assignment_changes = sum(c['detailed_metrics'].get('read_assignment_changes', 0) for c in detailed_comparisons)
        
        avg_stability = np.mean([c['detailed_metrics'].get('cluster_stability_score', 1.0) for c in detailed_comparisons])
        
        print(f"Total reads exchanged between clusters: {total_exchanged}")
        print(f"Total reads unique to NOSEED: {total_unique_noseed}")
        print(f"Total reads unique to SEED: {total_unique_seed}")
        print(f"Total read assignment changes: {total_assignment_changes}")
        print(f"Average cluster stability score: {avg_stability:.3f}")
        
        # Perfect vs partial matches
        perfect_matches = sum(c['detailed_metrics'].get('perfect_cluster_matches', 0) for c in detailed_comparisons)
        partial_matches = sum(c['detailed_metrics'].get('partial_cluster_matches', 0) for c in detailed_comparisons)
        different_clusters = sum(c['detailed_metrics'].get('completely_different_clusters', 0) for c in detailed_comparisons)
        
        print(f"Perfect cluster matches across all files: {perfect_matches}")
        print(f"Partial cluster matches across all files: {partial_matches}")
        print(f"Completely different clusters across all files: {different_clusters}")
        
        # Files with significant differences
        significant_diffs = [c for c in detailed_comparisons 
                           if c['detailed_metrics'].get('reads_exchanged', 0) > 0 
                           or c['detailed_metrics'].get('cluster_stability_score', 1.0) < 0.9]
        
        if significant_diffs:
            print(f"\n=== FILES WITH SIGNIFICANT DIFFERENCES ({len(significant_diffs)}) ===")
            for comp in significant_diffs:
                metrics = comp['detailed_metrics']
                print(f"  {comp['filename']}:")
                print(f"    Reads exchanged: {metrics.get('reads_exchanged', 0)}")
                print(f"    Stability score: {metrics.get('cluster_stability_score', 1.0):.3f}")
                print(f"    Assignment changes: {metrics.get('read_assignment_changes', 0)}")
        
        # Time analysis
        time_differences = []
        for c in detailed_comparisons:
            noseed_time = c.get('noseed_time', 0)
            seed_time = c.get('seed_time', 0)
            if noseed_time > 0 and seed_time > 0:
                ratio = max(noseed_time, seed_time) / min(noseed_time, seed_time)
                time_differences.append({
                    'filename': c['filename'],
                    'noseed_time': noseed_time,
                    'seed_time': seed_time,
                    'ratio': ratio
                })
        
        if time_differences:
            avg_ratio = np.mean([td['ratio'] for td in time_differences])
            print(f"\n=== EXECUTION TIME ANALYSIS ===")
            print(f"Average time ratio (slower/faster): {avg_ratio:.2f}")
            
            significant_time_diffs = [td for td in time_differences if td['ratio'] > 2.0]
            if significant_time_diffs:
                print(f"Files with >2x time difference ({len(significant_time_diffs)}):")
                for td in significant_time_diffs:
                    faster = "NOSEED" if td['noseed_time'] < td['seed_time'] else "SEED"
                    print(f"  {td['filename']}: {td['ratio']:.1f}x ({faster} faster)")

def process_all_matrices(base_dir: str = None, enable_comparison: bool = False) -> List[Dict[str, Any]]:
    """
    Process all CSV matrices and collect statistics.
    
    Parameters
    ----------
    base_dir : str, optional
        Base directory containing matrices
    enable_comparison : bool, optional
        If True, compare noseed vs seed versions
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
    if enable_comparison:
        print("🔍 COMPARISON MODE ENABLED - Will compare noseed vs seed versions")
    
    # Find all CSV files in directory and subdirectories
    csv_files = list(base_path.rglob("*.csv"))
    
    print(f"Found {len(csv_files)} CSV files to process")
    
    if len(csv_files) == 0:
        print("No CSV files found in matrice directory")
        print("Please check that your data files are in the correct location:")
        print(f"  {base_path.resolve()}")
        return results
    
    comparison_summary = []
    
    for i, csv_file in enumerate(csv_files):
        try:
            print(f"Processing {i+1}/{len(csv_files)}: {csv_file}")
            
            # Extract haplotype count from directory structure
            haplotype_count = extract_haplotype_count(str(csv_file))
            
            # Load and process matrix
            start_time = time.time()
            rows_data, cols_data, edges_0, edges_1, row_names, col_names, df, comp_df = get_data(str(csv_file))
            load_time = time.time() - start_time
            
            # Convert to binary matrix (ensure 0/1 values)
            binary_matrix = df.values.astype(int)
            
            # Run pipeline with optional comparison
            try:
                pipeline_results = run_pipeline(
                    binary_matrix,
                    min_col_quality=3,
                    min_row_quality=5,
                    error_rate=0.025,
                    filename=csv_file.name,
                    haplotype_count=haplotype_count,
                    enable_comparison=enable_comparison
                )
            except Exception as pipeline_error:
                print(f"  -> Pipeline error: {str(pipeline_error)}")
                # Create minimal pipeline results for error case
                pipeline_results = {
                    'ilp_calls_total': 0,
                    'preprocessing_time': 0,
                    'clustering_time': 0,
                    'total_time': 0,
                    'ilp_time_total': 0,
                    'ilp_time_ratio': 0,
                    'avg_ilp_time': 0,
                    'patterns_found': 0,
                    'regions_found': 0,
                    'clustering_steps': 0,
                    'regions_processed': 0,
                    'matrix_operations': 0,
                    'avg_region_size': 0,
                    'largest_region': 0,
                    'smallest_region': 0,
                    'solver_status_counts': {}
                }
            
            # Process comparison results if available
            if enable_comparison and 'comparison' in pipeline_results:
                comparison = pipeline_results['comparison']
                comp_summary = {
                    'filename': csv_file.name,
                    'noseed_success': comparison.get('noseed_success', False),
                    'seed_success': comparison.get('seed_success', False),
                    'results_identical': comparison.get('clustering_identical', False),
                    'results_equivalent': comparison.get('clustering_equivalent', False),
                    'noseed_steps': comparison.get('noseed_steps', 0),
                    'seed_steps': comparison.get('seed_steps', 0),
                    'noseed_time': comparison.get('noseed_time', 0),
                    'seed_time': comparison.get('seed_time', 0)
                }
                
                # Add detailed metrics if available
                if 'detailed_metrics' in comparison:
                    comp_summary['detailed_metrics'] = comparison['detailed_metrics']
                
                comparison_summary.append(comp_summary)
            
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
            
            # Add comparison results if available
            if enable_comparison and 'comparison' in pipeline_results:
                comparison = pipeline_results['comparison']
                stats.update({
                    'comparison_enabled': True,
                    'noseed_success': comparison.get('noseed_success', False),
                    'seed_success': comparison.get('seed_success', False),
                    'results_identical': comparison.get('clustering_identical', False),
                    'results_equivalent': comparison.get('clustering_equivalent', False),
                    'noseed_steps': comparison.get('noseed_steps', 0),
                    'seed_steps': comparison.get('seed_steps', 0),
                    'step_differences_count': len(comparison.get('step_differences', []))
                })
                
                # Add detailed metrics
                if 'detailed_metrics' in comparison:
                    detailed = comparison['detailed_metrics']
                    stats.update({
                        'reads_exchanged': detailed.get('reads_exchanged', 0),
                        'reads_unique_noseed': detailed.get('reads_unique_noseed', 0),
                        'reads_unique_seed': detailed.get('reads_unique_seed', 0),
                        'read_assignment_changes': detailed.get('read_assignment_changes', 0),
                        'cluster_stability_score': detailed.get('cluster_stability_score', 1.0),
                        'perfect_cluster_matches': detailed.get('perfect_cluster_matches', 0),
                        'partial_cluster_matches': detailed.get('partial_cluster_matches', 0),
                        'completely_different_clusters': detailed.get('completely_different_clusters', 0)
                    })
            
            # Add solver status counts if available
            solver_stats = pipeline_results.get('solver_status_counts', {})
            stats.update({
                'solver_success_count': solver_stats.get('success', 0),
                'solver_failed_count': solver_stats.get('failed', 0),
                'solver_error_count': solver_stats.get('error', 0),
            })
            
            results.append(stats)
            print(f"  -> Completed: {stats['ilp_calls_total']} ILP calls, {stats['patterns_found']} patterns, {stats['total_time']:.3f}s")
            
            if enable_comparison and 'comparison' in pipeline_results:
                comp = pipeline_results['comparison']
                print(f"  -> Comparison: NOSEED({comp.get('noseed_steps', 0)} steps), SEED({comp.get('seed_steps', 0)} steps), "
                      f"Identical: {comp.get('clustering_identical', False)}, Equivalent: {comp.get('clustering_equivalent', False)}")
            
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
    
    # Print enhanced comparison summary
    if enable_comparison and comparison_summary:
        print_comparison_summary(comparison_summary)
    
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

def run_experiments(enable_comparison: bool = False):
    """Main function to run all experiments."""
    
    print("Starting matrix processing experiments...")
    if enable_comparison:
        print("🔍 Comparison mode enabled - will compare noseed vs seed versions")
    print(f"Working directory: {os.getcwd()}")
    
    # Process all matrices
    results = process_all_matrices(enable_comparison=enable_comparison)
    
    # Save results
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    output_file = f"experiment_results_{timestamp}.csv"
    if enable_comparison:
        output_file = f"comparison_experiment_results_{timestamp}.csv"
    
    save_results_to_csv(results, output_file)
    
    return results

if __name__ == "__main__":
    # Check for comparison flag in command line arguments
    enable_comparison = '--compare' in sys.argv or '-c' in sys.argv
    
    if enable_comparison:
        print("🔍 Comparison mode activated via command line")
    
    results = run_experiments(enable_comparison=enable_comparison)
