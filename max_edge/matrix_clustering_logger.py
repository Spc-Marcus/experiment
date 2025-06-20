import numpy as np
import json
import time
from pathlib import Path
from typing import List, Tuple, Dict, Any
from clustering import clustering_full_matrix, find_quasi_biclique, clustering_step_with_stats

class ClusteringLogger:
    """Logger for clustering operations that saves all intermediate results."""
    
    def __init__(self, output_file: str = "clustering_results.json"):
        self.output_file = output_file
        self.results = {
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "matrix_info": {},
            "quasi_biclique_results": [],
            "clustering_steps": [],
            "final_results": {}
        }
    
    def save_matrix_to_txt(self, matrix: np.ndarray, filename: str = None):
        """Save matrix to text file for visualization."""
        if filename is None:
            filename = self.output_file.replace('.json', '_matrix.txt')
        
        try:
            with open(filename, 'w') as f:
                f.write(f"Matrix shape: {matrix.shape[0]} rows x {matrix.shape[1]} columns\n")
                f.write(f"Total elements: {matrix.size}\n")
                f.write(f"Ones count: {np.sum(matrix == 1)}\n")
                f.write(f"Zeros count: {np.sum(matrix == 0)}\n")
                f.write(f"Minus ones count: {np.sum(matrix == -1)}\n")
                f.write(f"Density: {np.sum(matrix == 1) / matrix.size:.4f}\n")
                f.write(f"Timestamp: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write("\nMatrix content:\n")
                f.write("=" * 50 + "\n")
                
                # Column headers
                f.write("    ")
                for j in range(matrix.shape[1]):
                    f.write(f"{j:3d}")
                f.write("\n")
                
                # Matrix rows with row indices
                for i in range(matrix.shape[0]):
                    f.write(f"{i:3d} ")
                    for j in range(matrix.shape[1]):
                        if matrix[i, j] == -1:
                            f.write("  -")
                        else:
                            f.write(f"{matrix[i, j]:3d}")
                    f.write("\n")
            
            print(f"Matrix saved to {filename}")
        except Exception as e:
            print(f"Error saving matrix to file: {e}")
    
    def log_matrix_info(self, matrix: np.ndarray):
        """Log basic matrix information."""
        self.results["matrix_info"] = {
            "shape": list(matrix.shape),
            "total_elements": int(matrix.size),
            "ones_count": int(np.sum(matrix == 1)),
            "zeros_count": int(np.sum(matrix == 0)),
            "minus_ones_count": int(np.sum(matrix == -1)),
            "density": float(np.sum(matrix == 1) / matrix.size)
        }
    
    def log_quasi_biclique(self, matrix: np.ndarray, rows: List[int], cols: List[int], success: bool, phase: str = ""):
        """Log quasi-biclique detection results."""
        result = {
            "phase": phase,
            "input_shape": list(matrix.shape),
            "selected_rows": rows,
            "selected_cols": cols,
            "success": success,
            "selected_area": len(rows) * len(cols) if rows and cols else 0,
            "timestamp": time.strftime("%H:%M:%S")
        }
        
        if rows and cols:
            # Calculate density of selected region
            selected_matrix = matrix[np.ix_(rows, cols)]
            result["selected_density"] = float(np.sum(selected_matrix == 1) / selected_matrix.size)
            result["selected_ones"] = int(np.sum(selected_matrix == 1))
            result["selected_zeros"] = int(np.sum(selected_matrix == 0))
        
        self.results["quasi_biclique_results"].append(result)
    
    def log_clustering_step(self, reads1: List[int], reads0: List[int], cols: List[int], ilp_calls: int, step_num: int):
        """Log clustering step results."""
        step_result = {
            "step_number": step_num,
            "reads_group1": reads1,
            "reads_group0": reads0,
            "columns": cols,
            "ilp_calls": ilp_calls,
            "group1_size": len(reads1),
            "group0_size": len(reads0),
            "columns_count": len(cols),
            "timestamp": time.strftime("%H:%M:%S")
        }
        self.results["clustering_steps"].append(step_result)
    
    def log_final_results(self, final_steps: List[Tuple[List[int], List[int], List[int]]], total_time: float):
        """Log final clustering results."""
        self.results["final_results"] = {
            "total_steps": len(final_steps),
            "total_time_seconds": total_time,
            "valid_steps": []
        }
        
        for i, (reads1, reads0, cols) in enumerate(final_steps):
            step_info = {
                "step_id": i,
                "reads_group1": reads1,
                "reads_group0": reads0,
                "columns": cols,
                "group1_size": len(reads1),
                "group0_size": len(reads0),
                "columns_count": len(cols)
            }
            self.results["final_results"]["valid_steps"].append(step_info)
    
    def save_results(self):
        """Save all logged results to file."""
        with open(self.output_file, 'w') as f:
            json.dump(self.results, f, indent=2)
        print(f"Results saved to {self.output_file}")

def process_matrix_with_logging(matrix: np.ndarray, 
                               regions: List[List[int]] = None,
                               steps: List[Tuple[List[int], List[int], List[int]]] = None,
                               output_file: str = "clustering_results.json",
                               **kwargs) -> List[Tuple[List[int], List[int], List[int]]]:
    """
    Process matrix with full logging of all intermediate results.
    
    Parameters:
    -----------
    matrix : np.ndarray
        Input binary matrix
    regions : List[List[int]], optional
        Column regions to process
    steps : List[Tuple[List[int], List[int], List[int]]], optional
        Initial clustering steps
    output_file : str
        Output file for results
    **kwargs : additional arguments for clustering_full_matrix
    
    Returns:
    --------
    List[Tuple[List[int], List[int], List[int]]]
        Final clustering steps
    """
    
    logger = ClusteringLogger(output_file)
    start_time = time.time()
    
    # Save matrix to text file first
    logger.save_matrix_to_txt(matrix)
    
    # Log matrix information
    logger.log_matrix_info(matrix)
    print(f"Processing matrix of shape {matrix.shape}")
    
    # Default parameters
    if regions is None:
        regions = [list(range(matrix.shape[1]))]  # Use all columns as one region
    if steps is None:
        steps = []
    
    # Process with full matrix clustering
    final_steps = clustering_full_matrix(
        matrix, 
        regions, 
        steps, 
        **kwargs
    )
    
    total_time = time.time() - start_time
    
    # Log final results
    logger.log_final_results(final_steps, total_time)
    
    # Save all results
    logger.save_results()
    
    print(f"Clustering completed in {total_time:.2f} seconds")
    print(f"Found {len(final_steps)} valid clustering steps")
    
    return final_steps

def main():
    """Example usage of the matrix clustering logger."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Process matrix with clustering and save all intermediate results")
    parser.add_argument("--matrix-file", type=str, help="Path to numpy matrix file (.npy)")
    parser.add_argument("--output", type=str, default="clustering_results.json", help="Output file for results")
    parser.add_argument("--error-rate", type=float, default=0.025, help="Error rate for clustering")
    parser.add_argument("--min-row-quality", type=int, default=5, help="Minimum row quality")
    parser.add_argument("--min-col-quality", type=int, default=3, help="Minimum column quality")
    
    args = parser.parse_args()
    
    if args.matrix_file:
        # Load matrix from file
        matrix = np.load(args.matrix_file)
        print(f"Loaded matrix from {args.matrix_file}")
    else:
        # Create example matrix
        np.random.seed(42)
        matrix = np.random.choice([0, 1, -1], size=(20, 15), p=[0.3, 0.6, 0.1])
        print("Using example random matrix")
    
    # Process matrix with logging
    results = process_matrix_with_logging(
        matrix,
        output_file=args.output,
        error_rate=args.error_rate,
        min_row_quality=args.min_row_quality,
        min_col_quality=args.min_col_quality
    )
    
    print(f"\nProcessing complete. Check {args.output} for detailed results.")

if __name__ == "__main__":
    main()

