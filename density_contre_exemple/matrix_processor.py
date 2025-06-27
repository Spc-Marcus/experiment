import numpy as np
import pandas as pd
import time
import argparse
from pathlib import Path
from typing import List, Tuple
from clustering import clustering_full_matrix, save_results_to_file

def load_csv_matrix(csv_file_path: str) -> np.ndarray:
    """
    Load CSV matrix with only 0s and 1s.
    
    Parameters
    ----------
    csv_file_path : str
        Path to CSV file containing binary matrix
        
    Returns
    -------
    np.ndarray
        Binary matrix with values 0 and 1
    """
    # Load CSV file
    df = pd.read_csv(csv_file_path, index_col=0)
    
    # Ensure binary values (0, 1)
    df = df.fillna(0).astype(int)
    df = df.clip(0, 1)  # Ensure only 0 and 1 values
    
    print(f"Loaded matrix from {csv_file_path}")
    print(f"Shape: {df.shape}")
    print(f"Density: {df.sum().sum() / df.size:.3f}")
    
    return df.values

def process_matrix_with_logging(
    matrix: np.ndarray,
    output_file: str = "clustering_results.json",
    regions: List[List[int]] = None,
    steps: List[Tuple[List[int], List[int], List[int]]] = None,
    error_rate: float = 0.025,
    min_row_quality: int = 5,
    min_col_quality: int = 3
) -> List[Tuple[List[int], List[int], List[int]]]:
    """
    Process matrix with full logging of quasi-biclique and clustering steps.
    
    Parameters
    ----------
    matrix : np.ndarray
        Input binary matrix (0s and 1s only)
    output_file : str
        Output file for logging results
    regions : List[List[int]], optional
        Column regions to process
    steps : List[Tuple[List[int], List[int], List[int]]], optional
        Initial clustering steps
    error_rate : float
        Error rate for clustering
    min_row_quality : int
        Minimum row quality
    min_col_quality : int
        Minimum column quality
        
    Returns
    -------
    List[Tuple[List[int], List[int], List[int]]]
        Final clustering steps
    """
    
    start_time = time.time()
    
    # Default parameters
    if regions is None:
        regions = [list(range(matrix.shape[1]))]  # Use all columns as one region
    if steps is None:
        steps = []
    
    print(f"Processing matrix of shape {matrix.shape}")
    print(f"Using {len(regions)} regions")
    print(f"Error rate: {error_rate}")
    print(f"Min row quality: {min_row_quality}")
    print(f"Min column quality: {min_col_quality}")
    
    # Process with clustering and logging
    final_steps = clustering_full_matrix(
        matrix,
        regions,
        steps,
        min_row_quality=min_row_quality,
        min_col_quality=min_col_quality,
        error_rate=error_rate,
        log_results=True,
        filename=Path(output_file).stem
    )
    
    total_time = time.time() - start_time
    
    # Save results to file
    save_results_to_file(output_file)
    
    print(f"\nProcessing completed in {total_time:.2f} seconds")
    print(f"Found {len(final_steps)} valid clustering steps")
    print(f"Results saved to {output_file}")
    
    return final_steps

def main():
    """Main program to process CSV matrix with clustering and logging."""
    parser = argparse.ArgumentParser(
        description="Process CSV matrix (0s and 1s) with clustering and save all intermediate results"
    )
    
    parser.add_argument("csv_file", type=str, help="Path to CSV file containing binary matrix")
    parser.add_argument("--output", type=str, default="clustering_results_base.json", 
                       help="Output file for results (default: clustering_results_base.json)")
    parser.add_argument("--error-rate", type=float, default=0.025, 
                       help="Error rate for clustering (default: 0.025)")
    parser.add_argument("--min-row-quality", type=int, default=5, 
                       help="Minimum row quality (default: 5)")
    parser.add_argument("--min-col-quality", type=int, default=3, 
                       help="Minimum column quality (default: 3)")
    
    args = parser.parse_args()
    
    # Verify CSV file exists
    csv_path = Path(args.csv_file)
    if not csv_path.exists():
        print(f"Error: CSV file {args.csv_file} not found")
        return
    
    try:
        # Load matrix from CSV
        matrix = load_csv_matrix(args.csv_file)
        
        # Process matrix with logging
        results = process_matrix_with_logging(
            matrix,
            output_file=args.output,
            error_rate=args.error_rate,
            min_row_quality=args.min_row_quality,
            min_col_quality=args.min_col_quality
        )
        
        print(f"\nFinal results:")
        for i, (reads1, reads0, cols) in enumerate(results):
            print(f"  Step {i+1}: {len(reads1)} vs {len(reads0)} reads, {len(cols)} columns")
            
    except Exception as e:
        print(f"Error processing matrix: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
