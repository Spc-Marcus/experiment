import time
import functools
from typing import Dict, Any, Callable
from threading import local
import numpy as np

# Thread-local storage for tracking metrics
_tracker = local()

def get_tracker() -> Dict[str, Any]:
    """Get the current tracker state."""
    if not hasattr(_tracker, 'data'):
        reset_tracker()
    return _tracker.data

def reset_tracker():
    """Reset all tracking counters."""
    _tracker.data = {
        'ilp_calls_total': 0,
        'ilp_time_total': 0.0,
        'preprocessing_time': 0.0,
        'clustering_time': 0.0,
        'patterns_found': 0,
        'regions_processed': 0,
        'clustering_steps': 0,
        'optimization_phases': 0,
        'solver_status_counts': {},
        'matrix_operations': 0,
        'last_matrix_shape': None
    }

def track_ilp_call(phase: str = 'unknown'):
    """Decorator to track ILP solver calls and timing."""
    def decorator(func: Callable) -> Callable:
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            tracker = get_tracker()
            start_time = time.time()
            
            try:
                result = func(*args, **kwargs)
                execution_time = time.time() - start_time
                
                # Track the call
                tracker['ilp_calls_total'] += 1
                tracker['ilp_time_total'] += execution_time
                
                # Track solver status if available in result
                if isinstance(result, tuple) and len(result) >= 3:
                    # Assume third element indicates success
                    status = 'success' if result[2] else 'failed'
                    tracker['solver_status_counts'][status] = tracker['solver_status_counts'].get(status, 0) + 1
                
                # Track patterns found
                if isinstance(result, tuple) and len(result) >= 2:
                    if len(result[0]) > 0 and len(result[1]) > 0:
                        tracker['patterns_found'] += 1
                
                return result
            except Exception as e:
                execution_time = time.time() - start_time
                tracker['ilp_calls_total'] += 1
                tracker['ilp_time_total'] += execution_time
                tracker['solver_status_counts']['error'] = tracker['solver_status_counts'].get('error', 0) + 1
                raise e
                
        return wrapper
    return decorator

def track_preprocessing():
    """Decorator to track preprocessing operations."""
    def decorator(func: Callable) -> Callable:
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            tracker = get_tracker()
            start_time = time.time()
            
            result = func(*args, **kwargs)
            execution_time = time.time() - start_time
            
            tracker['preprocessing_time'] += execution_time
            
            # Track matrix shape if available
            if args and hasattr(args[0], 'shape'):
                tracker['last_matrix_shape'] = args[0].shape
            
            # Track regions found
            if isinstance(result, tuple) and len(result) >= 2:
                if isinstance(result[1], list):
                    tracker['regions_processed'] = len(result[1])
                    
            return result
        return wrapper
    return decorator

def track_clustering():
    """Decorator to track clustering operations."""
    def decorator(func: Callable) -> Callable:
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            tracker = get_tracker()
            start_time = time.time()
            
            result = func(*args, **kwargs)
            execution_time = time.time() - start_time
            
            tracker['clustering_time'] += execution_time
            
            # Track clustering steps
            if isinstance(result, list):
                tracker['clustering_steps'] = len(result)
                
            return result
        return wrapper
    return decorator

def get_performance_summary() -> Dict[str, Any]:
    """Get a comprehensive performance summary."""
    tracker = get_tracker()
    
    total_time = tracker['preprocessing_time'] + tracker['clustering_time']
    
    summary = {
        'ilp_calls_total': tracker['ilp_calls_total'],
        'ilp_time_total': tracker['ilp_time_total'],
        'preprocessing_time': tracker['preprocessing_time'],
        'clustering_time': tracker['clustering_time'],
        'total_time': total_time,
        'patterns_found': tracker['patterns_found'],
        'regions_processed': tracker['regions_processed'],
        'clustering_steps': tracker['clustering_steps'],
        'optimization_phases': tracker['optimization_phases'],
        'matrix_operations': tracker['matrix_operations'],
        'solver_status_counts': tracker['solver_status_counts'].copy(),
        'matrix_shape': tracker['last_matrix_shape'],
        'avg_ilp_time': tracker['ilp_time_total'] / max(1, tracker['ilp_calls_total']) if tracker['ilp_calls_total'] > 0 else 0,
        'ilp_time_ratio': tracker['ilp_time_total'] / max(total_time, 0.001),
    }
    
    # Add density info if matrix shape is available
    if tracker['last_matrix_shape']:
        summary['matrix_size'] = tracker['last_matrix_shape'][0] * tracker['last_matrix_shape'][1]
        
    return summary

# Context manager for tracking a complete pipeline run
class PipelineTracker:
    """Context manager for tracking complete pipeline execution."""
    
    def __enter__(self):
        reset_tracker()
        self.start_time = time.time()
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.total_time = time.time() - self.start_time
        
    def get_results(self) -> Dict[str, Any]:
        """Get final results including total execution time."""
        summary = get_performance_summary()
        summary['pipeline_total_time'] = getattr(self, 'total_time', 0)
        return summary

# Example usage functions
if __name__ == "__main__":
    # Example of how to use the decorators
    
    @track_ilp_call('test_phase')
    def mock_ilp_solver(matrix: np.ndarray) -> tuple:
        time.sleep(0.1)  # Simulate solving time
        rows = [0, 1, 2]
        cols = [0, 1]
        success = True
        return rows, cols, success
    
    @track_preprocessing()
    def mock_preprocessing(matrix: np.ndarray) -> tuple:
        time.sleep(0.05)
        regions = [[0, 1, 2], [3, 4]]
        steps = []
        return matrix, regions, steps
    
    @track_clustering()
    def mock_clustering(matrix: np.ndarray, regions, steps) -> list:
        time.sleep(0.2)
        # Mock some ILP calls
        for _ in range(3):
            mock_ilp_solver(matrix)
        return [([0, 1], [2, 3], [0, 1])]
    
    # Test the tracking
    with PipelineTracker() as tracker:
        test_matrix = np.random.choice([0, 1], size=(50, 20))
        
        matrix, regions, steps = mock_preprocessing(test_matrix)
        final_steps = mock_clustering(matrix, regions, steps)
        
        results = tracker.get_results()
        
    print("Performance Summary:")
    for key, value in results.items():
        print(f"{key}: {value}")
