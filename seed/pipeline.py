import numpy as np
import sys
from pathlib import Path
from typing import Dict, Any

# Fix import paths
try:
    from .preprocess import pre_processing
    from .clustering import clustering_full_matrix
except ImportError:
    # Add parent directory to path and try absolute imports
    current_dir = Path(__file__).parent
    parent_dir = current_dir.parent
    sys.path.insert(0, str(parent_dir))
    
    from ilphaplo.preprocess import pre_processing
    from ilphaplo.clustering import clustering_full_matrix

# Import decorators with fallback
try:
    from ..decorators.ilp_tracker import PipelineTracker, get_performance_summary
except ImportError:
    try:
        sys.path.insert(0, str(Path(__file__).parent.parent))
        from decorators.ilp_tracker import PipelineTracker, get_performance_summary
    except ImportError:
        # Fallback implementation
        class PipelineTracker:
            def __init__(self):
                self.start_time = None
                
            def __enter__(self):
                import time
                self.start_time = time.time()
                return self
                
            def __exit__(self, *args):
                pass
                
            def get_results(self):
                import time
                total_time = time.time() - self.start_time if self.start_time else 0
                return {
                    'ilp_calls_total': 0, 
                    'total_time': total_time, 
                    'patterns_found': 0,
                    'preprocessing_time': 0,
                    'clustering_time': 0,
                    'ilp_time_total': 0,
                    'regions_processed': 0,
                    'clustering_steps': 0,
                    'matrix_operations': 0,
                    'solver_status_counts': {},
                    'ilp_time_ratio': 0,
                    'avg_ilp_time': 0
                }
        
        def get_performance_summary():
            return {}

def run_pipeline(binary_matrix: np.ndarray, 
                min_col_quality: int = 3,
                min_row_quality: int = 5,
                error_rate: float = 0.025,
                filename: str = "unknown",
                haplotype_count: int = 0,
                X_factor: int = 3,
                step_n: int = None) -> Dict[str, Any]:
    """
    Execute complete preprocessing and clustering pipeline with seeding parameters.
    
    Parameters
    ----------
    binary_matrix : np.ndarray
        Input binary matrix (0,1) to process
    min_col_quality : int, optional
        Minimum columns for valid regions (default: 3)
    min_row_quality : int, optional
        Minimum rows for valid clusters (default: 5)
    error_rate : float, optional
        ILP optimization error tolerance (default: 0.025)
    filename : str, optional
        Original filename for logging (default: "unknown")
    haplotype_count : int, optional
        Number of haplotypes for organization (default: 0)
    X_factor : int, optional
        Seeding division factor (default: 3)
    step_n : int, optional
        Seeding step size (default: auto-calculated)
        
    Returns
    -------
    Dict[str, Any]
        Performance metrics automatically captured by decorators
    """
    
    with PipelineTracker() as tracker:
        # Phase 1: Preprocessing (tracked automatically)
        matrix, inhomogeneous_regions, steps = pre_processing(
            binary_matrix, 
            min_col_quality=min_col_quality
        )
        
        # Phase 2: Clustering with seeding parameters
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
        
        # Get comprehensive results from tracker
        results = tracker.get_results()
        
        # Add additional computed metrics
        results.update({
            'matrix_shape': binary_matrix.shape,
            'regions_found': len(inhomogeneous_regions),
            'inhomogeneous_regions': len(inhomogeneous_regions),
            'matrix_density': np.sum(binary_matrix) / binary_matrix.size,
            'preprocessing_regions': len(steps),
            'avg_region_size': np.mean([len(r) for r in inhomogeneous_regions]) if inhomogeneous_regions else 0,
            'largest_region': max([len(r) for r in inhomogeneous_regions]) if inhomogeneous_regions else 0,
            'smallest_region': min([len(r) for r in inhomogeneous_regions]) if inhomogeneous_regions else 0
        })
        
    return results


def benchmark_pipeline(binary_matrix: np.ndarray, 
                      iterations: int = 1,
                      **kwargs) -> Dict[str, Any]:
    """
    Run pipeline multiple times and return aggregated statistics.
    
    Parameters
    ----------
    binary_matrix : np.ndarray
        Input binary matrix to process
    iterations : int, optional
        Number of times to run pipeline (default: 1)
    **kwargs
        Additional parameters passed to run_pipeline
        
    Returns
    -------
    Dict[str, Any]
        Aggregated performance metrics over all iterations
    """
    
    results = []
    
    for i in range(iterations):
        result = run_pipeline(binary_matrix, **kwargs)
        results.append(result)
    
    # Aggregate results
    aggregated = {
        'iterations': iterations,
        'matrix_shape': binary_matrix.shape,
        'avg_ilp_calls': np.mean([r['ilp_calls_total'] for r in results]),
        'total_ilp_calls': np.sum([r['ilp_calls_total'] for r in results]),
        'avg_preprocessing_time': np.mean([r['preprocessing_time'] for r in results]),
        'avg_clustering_time': np.mean([r['clustering_time'] for r in results]),
        'avg_total_time': np.mean([r['total_time'] for r in results]),
        'total_execution_time': np.sum([r['total_time'] for r in results]),
        'avg_patterns_found': np.mean([r['patterns_found'] for r in results]),
        'avg_regions_found': np.mean([r['regions_found'] for r in results]),
        'matrix_density': results[0]['matrix_density'],  # Same for all iterations
        'std_ilp_calls': np.std([r['ilp_calls_total'] for r in results]),
        'std_total_time': np.std([r['total_time'] for r in results]),
        'min_ilp_calls': np.min([r['ilp_calls_total'] for r in results]),
        'max_ilp_calls': np.max([r['ilp_calls_total'] for r in results])
    }
    
    return aggregated


# Example usage
if __name__ == "__main__":
    # Generate sample binary matrix
    np.random.seed(42)
    test_matrix = np.random.choice([0, 1], size=(100, 50), p=[0.7, 0.3])
    
    # Single run
    results = run_pipeline(test_matrix)
    print(f"ILP calls: {results['ilp_calls_total']}")
    print(f"Total time: {results['total_time']:.3f}s")
    print(f"Patterns found: {results['patterns_found']}")
    
    # Benchmark multiple runs
    benchmark_results = benchmark_pipeline(test_matrix, iterations=5)
    print(f"Average ILP calls: {benchmark_results['avg_ilp_calls']:.1f}")
    print(f"Average time: {benchmark_results['avg_total_time']:.3f}s")
