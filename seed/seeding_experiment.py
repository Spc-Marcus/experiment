import os
import sys
import pandas as pd
import numpy as np
import time
import itertools
from pathlib import Path
from typing import Dict, List, Any, Tuple

# Fix import paths
if __name__ == "__main__":
    current_dir = Path(__file__).parent
    parent_dir = current_dir.parent
    sys.path.insert(0, str(parent_dir))
    sys.path.insert(0, str(current_dir))  # Add current directory too
    
    try:
        from ilphaplo.get_data import get_data
        from ilphaplo.preprocess import pre_processing  
        from ilphaplo.clustering import clustering_full_matrix
        from ilphaplo.run_experiments import find_matrice_directory, extract_haplotype_count
    except ImportError as e:
        print(f"Import error: {e}")
        print("Trying direct imports...")
        from get_data import get_data
        from preprocess import pre_processing
        from clustering import clustering_full_matrix
        from run_experiments import find_matrice_directory, extract_haplotype_count
else:
    from .get_data import get_data
    from .preprocess import pre_processing
    from .clustering import clustering_full_matrix
    from .run_experiments import find_matrice_directory, extract_haplotype_count

def run_seeding_experiment(
    binary_matrix: np.ndarray,
    X_factor: int,
    step_n: int,
    filename: str = "unknown",
    haplotype_count: int = 0,
    min_col_quality: int = 3,
    min_row_quality: int = 5,
    error_rate: float = 0.025
) -> Dict[str, Any]:
    """
    Run pipeline with specific seeding parameters by calling functions directly.
    """
    start_time = time.time()
    
    try:
        # Phase 1: Preprocessing
        preprocessing_start = time.time()
        matrix, inhomogeneous_regions, steps = pre_processing(
            binary_matrix, 
            min_col_quality=min_col_quality
        )
        preprocessing_time = time.time() - preprocessing_start
        
        # Phase 2: Clustering with seeding parameters
        clustering_start = time.time()
        
        # Check if clustering_full_matrix accepts the new parameters
        import inspect
        sig = inspect.signature(clustering_full_matrix)
        
        if 'X_factor' in sig.parameters and 'step_n' in sig.parameters:
            # Function supports seeding parameters
            final_steps = clustering_full_matrix(
                matrix,
                regions=inhomogeneous_regions,
                steps=steps,
                min_row_quality=min_row_quality,
                min_col_quality=min_col_quality,
                error_rate=error_rate,
                filename=filename,
                haplotype_count=haplotype_count,
                X_factor=X_factor,
                step_n=step_n
            )
            seeding_used = True
        else:
            # Function doesn't support seeding parameters, use default
            print(f"Info: Using default clustering (no seeding) for X={X_factor}, step={step_n}")
            final_steps = clustering_full_matrix(
                matrix,
                regions=inhomogeneous_regions,
                steps=steps,
                min_row_quality=min_row_quality,
                min_col_quality=min_col_quality,
                error_rate=error_rate,
                filename=filename,
                haplotype_count=haplotype_count
            )
            seeding_used = False
        
        clustering_time = time.time() - clustering_start
        total_time = time.time() - start_time
        
        # Count ILP calls from final steps
        ilp_calls = len([step for step in final_steps if len(step[0]) > 0 and len(step[1]) > 0])
        
        return {
            'read_name': filename,
            'nb_haplotypes': haplotype_count,
            'matrix_m': binary_matrix.shape[0],
            'matrix_n': binary_matrix.shape[1],
            'total_time': total_time,
            'ilp_time': clustering_time * 0.8,  # Approximation that ILP is 80% of clustering time
            'X_factor': X_factor,
            'step_n': step_n,
            'nb_ilp_calculated': ilp_calls,
            'patterns_found': len(final_steps),
            'regions_found': len(inhomogeneous_regions),
            'preprocessing_time': preprocessing_time,
            'clustering_time': clustering_time,
            'matrix_density': np.sum(binary_matrix) / binary_matrix.size,
            'seeding_used': seeding_used
        }
        
    except Exception as e:
        print(f"Error in seeding experiment: {e}")
        import traceback
        traceback.print_exc()
        return {
            'read_name': filename,
            'nb_haplotypes': haplotype_count,
            'matrix_m': binary_matrix.shape[0],
            'matrix_n': binary_matrix.shape[1],
            'total_time': 0,
            'ilp_time': 0,
            'X_factor': X_factor,
            'step_n': step_n,
            'nb_ilp_calculated': 0,
            'patterns_found': 0,
            'error': str(e),
            'seeding_used': False
        }

def process_seeding_experiments(base_dir: str = None) -> List[Dict[str, Any]]:
    """
    Process all matrices with different seeding configurations.
    Only processes matrices with haplotype count >= 4.
    """
    results = []
    
    # Find matrice directory
    if base_dir is None:
        base_path = find_matrice_directory()
        if base_path is None:
            print("Error: Could not find matrice directory")
            return results
    else:
        base_path = Path(base_dir)
    
    # Find all CSV files
    csv_files = list(base_path.rglob("*.csv"))
    print(f"Found {len(csv_files)} CSV files for seeding experiments")
    
    # Seeding parameter combinations
    X_factors = [2, 3, 4, 5, 6]
    step_values = [2, 4, 6, 8, 10, 12, 14]
    
    # Filter CSV files by size AND haplotype count
    suitable_files = []
    for csv_file in csv_files:
        try:
            # Quick size check
            df_temp = pd.read_csv(csv_file, index_col=0)
            haplotype_count = extract_haplotype_count(str(csv_file))
            
            # Filter by haplotype count >= 4 AND matrix size
            if (haplotype_count >= 4 and 
                df_temp.shape[0] >= 50 and 
                df_temp.shape[1] >= 20):
                suitable_files.append(csv_file)
                print(f"Including: {csv_file.name} (haplotypes: {haplotype_count}, size: {df_temp.shape})")
            else:
                if haplotype_count < 4:
                    print(f"Skipping {csv_file.name}: haplotype count {haplotype_count} < 4")
                else:
                    print(f"Skipping {csv_file.name}: small matrix {df_temp.shape}")
        except Exception as e:
            print(f"Error checking {csv_file}: {e}")
    
    print(f"Found {len(suitable_files)} suitable files (haplotypes >= 4) for seeding experiments")
    
    # Check if seeding is actually supported
    try:
        from clustering import clustering_full_matrix
        import inspect
        sig = inspect.signature(clustering_full_matrix)
        seeding_supported = 'X_factor' in sig.parameters and 'step_n' in sig.parameters
        
        if not seeding_supported:
            print("WARNING: Seeding parameters not supported in clustering_full_matrix!")
            print("All experiments will use default clustering parameters.")
            print("To enable seeding, ensure clustering.py has X_factor and step_n parameters.")
    except:
        seeding_supported = False
        print("Could not determine seeding support status.")
    
    total_experiments = len(suitable_files) * len(X_factors) * len(step_values)
    experiment_count = 0
    
    for csv_file in suitable_files:
        try:
            print(f"Processing file: {csv_file.name}")
            
            # Load matrix data
            rows_data, cols_data, edges_0, edges_1, row_names, col_names, df, comp_df = get_data(str(csv_file))
            binary_matrix = df.values.astype(int)
            haplotype_count = extract_haplotype_count(str(csv_file))
            
            # Double-check haplotype count (safety check)
            if haplotype_count < 4:
                print(f"  Skipping: haplotype count {haplotype_count} < 4")
                continue
            
            for X_factor in X_factors:
                for step_n in step_values:
                    experiment_count += 1
                    print(f"  Experiment {experiment_count}/{total_experiments}: X={X_factor}, step={step_n}")
                    
                    try:
                        result = run_seeding_experiment(
                            binary_matrix,
                            X_factor=X_factor,
                            step_n=step_n,
                            filename=csv_file.name,
                            haplotype_count=haplotype_count
                        )
                        
                        # Only save if ILP was actually used (nb_ilp_calculated > 1)
                        if result['nb_ilp_calculated'] > 1:
                            results.append(result)
                            seeding_status = "with seeding" if result.get('seeding_used', False) else "default clustering"
                            print(f"    -> Saved: {result['nb_ilp_calculated']} ILP calls, {result['total_time']:.3f}s ({seeding_status})")
                        else:
                            print(f"    -> Skipped: only {result['nb_ilp_calculated']} ILP calls")
                            
                    except Exception as e:
                        print(f"    -> Error: {str(e)}")
                        
        except Exception as e:
            print(f"Error processing {csv_file}: {str(e)}")
    
    return results

def analyze_seeding_results(results: List[Dict[str, Any]]) -> None:
    """
    Analyze seeding experiment results.
    """
    if not results:
        print("No results to analyze")
        return
    
    df = pd.DataFrame(results)
    
    print(f"\n=== SEEDING EXPERIMENT ANALYSIS ===")
    print(f"Total experiments with ILP > 1: {len(df)}")
    print(f"X_factor values tested: {sorted(df['X_factor'].unique())}")
    print(f"step_n values tested: {sorted(df['step_n'].unique())}")
    
    # Check if seeding was actually used
    if 'seeding_used' in df.columns:
        seeding_count = df['seeding_used'].sum()
        print(f"Experiments using seeding: {seeding_count} out of {len(df)}")
        if seeding_count == 0:
            print("WARNING: No experiments used seeding parameters!")
            print("This suggests the clustering function doesn't support seeding yet.")
            print("Results represent baseline performance only.")
    
    # Performance by X_factor
    print(f"\n=== PERFORMANCE BY X_FACTOR ===")
    x_factor_stats = df.groupby('X_factor').agg({
        'total_time': ['count', 'mean', 'std'],
        'nb_ilp_calculated': ['mean', 'std'],
        'ilp_time': ['mean', 'std']
    }).round(3)
    print(x_factor_stats)
    
    # Performance by step_n
    print(f"\n=== PERFORMANCE BY STEP_N ===")
    step_n_stats = df.groupby('step_n').agg({
        'total_time': ['count', 'mean', 'std'],
        'nb_ilp_calculated': ['mean', 'std'],
        'ilp_time': ['mean', 'std']
    }).round(3)
    print(step_n_stats)
    
    # If seeding was not used, note this in the analysis
    if 'seeding_used' in df.columns and df['seeding_used'].sum() == 0:
        print(f"\n=== IMPORTANT NOTE ===")
        print(f"All experiments used default clustering parameters.")
        print(f"To enable seeding experiments, ensure clustering.py supports X_factor and step_n parameters.")
        print(f"Current results represent baseline performance across different matrix files.")
    else:
        # Normal seeding analysis
        # Best combinations
        print(f"\n=== OPTIMAL COMBINATIONS ===")
        
        # Best by total time
        best_time = df.groupby(['X_factor', 'step_n'])['total_time'].mean().sort_values()
        print("Top 5 combinations by average total time:")
        for (x, s), time_val in best_time.head().items():
            count = len(df[(df['X_factor'] == x) & (df['step_n'] == s)])
            print(f"  X={x}, step_n={s}: {time_val:.3f}s (n={count})")
        
        # Best by ILP efficiency
        best_ilp = df.groupby(['X_factor', 'step_n'])['nb_ilp_calculated'].mean().sort_values()
        print("\nTop 5 combinations by fewest ILP calls:")
        for (x, s), ilp_val in best_ilp.head().items():
            count = len(df[(df['X_factor'] == x) & (df['step_n'] == s)])
            print(f"  X={x}, step_n={s}: {ilp_val:.1f} ILP calls (n={count})")
        
        # Matrix size impact
        print(f"\n=== IMPACT BY MATRIX SIZE ===")
        df['matrix_size'] = df['matrix_m'] * df['matrix_n']
        df['size_category'] = pd.cut(df['matrix_size'], 
                                   bins=3, 
                                   labels=['Small', 'Medium', 'Large'])
        
        size_impact = df.groupby(['size_category', 'X_factor'])['total_time'].mean().unstack()
        print("Average total time by matrix size and X_factor:")
        print(size_impact.round(3))

def save_seeding_results(results: List[Dict[str, Any]], output_file: str = None) -> None:
    """
    Save seeding experiment results to CSV.
    """
    if not results:
        print("No results to save")
        return
    
    if output_file is None:
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        output_file = f"seeding_experiment_results_{timestamp}.csv"
    
    df = pd.DataFrame(results)
    df.to_csv(output_file, index=False)
    print(f"Seeding experiment results saved to {output_file}")

def run_seeding_experiments():
    """
    Main function to run seeding experiments (haplotypes >= 4 only).
    """
    print("Starting ILP seeding parameter experiments...")
    print("Testing X_factors: [2, 3, 4, 5, 6]")
    print("Testing step_n values: [2, 4, 6, 8, 10, 12, 14]")
    print("Condition: Save only if nb_ilp_calculated > 1")
    print("Filter: Only processing matrices with haplotype count >= 4")
    
    # Run experiments
    results = process_seeding_experiments()
    
    # Save and analyze results
    save_seeding_results(results)
    analyze_seeding_results(results)
    
    return results

if __name__ == "__main__":
    results = run_seeding_experiments()
