from typing import List, Tuple
import numpy as np
import pulp as plp
import logging
import os
import sys
import time
import contextlib
from io import StringIO
from pathlib import Path
import json

# Fix decorator imports with multiple fallback attempts
try:
    from ..decorators.ilp_tracker import track_clustering, track_ilp_call
except ImportError:
    try:
        sys.path.insert(0, str(Path(__file__).parent.parent))
        from decorators.ilp_tracker import track_clustering, track_ilp_call
    except ImportError:
        def track_clustering():
            def decorator(func):
                return func
            return decorator
        
        def track_ilp_call(phase):
            def decorator(func):
                return func
            return decorator

# Fix decorateur import path
try:
    from ..decorateur.perf import print_decorator
except ImportError:
    try:
        sys.path.insert(0, str(Path(__file__).parent.parent))
        from decorateur.perf import print_decorator
    except ImportError:
        # Fallback if decorateur module is not available
        def print_decorator(name):
            def decorator(func):
                return func
            return decorator

# Remove undefined decorators and functions
def get_stats():
    return None

def timed_matrix_operation(name):
    def decorator(func):
        return func
    return decorator

def timed_ilp_call(name):
    def decorator(func):
        return func
    return decorator

logger = logging.getLogger(__name__)

# Suppress ALL Gurobi output
# os.environ['GRB_LICENSE_FILE'] = ''  # Remove this line
# gurobi_logger = logging.getLogger('gurobipy')  # Remove this line
# gurobi_logger.setLevel(logging.CRITICAL)  # Remove this line
# gurobi_logger.propagate = False  # Remove this line

# Remove Gurobi options
# options = {
# 	"WLSACCESSID":"af4b8280-70cd-47bc-aeef-69ecf14ecd10",
# 	"WLSSECRET":"04da6102-8eb3-4e38-ba06-660ea8f87bf2",
# 	"LICENSEID":2669217
# }

@contextlib.contextmanager
def suppress_pulp_output():
    """Context manager to completely suppress PuLP output"""
    old_stdout = sys.stdout
    old_stderr = sys.stderr
    try:
        sys.stdout = StringIO()
        sys.stderr = StringIO()
        yield
    finally:
        sys.stdout = old_stdout
        sys.stderr = old_stderr

# Add global logger for results
results_logger = {
    "quasi_biclique_results": [],
    "clustering_steps": [],
    "matrix_info": {}
}

def log_matrix_info(matrix: np.ndarray):
    """Log basic matrix information."""
    global results_logger
    results_logger["matrix_info"] = {
        "shape": list(matrix.shape),
        "total_elements": int(matrix.size),
        "ones_count": int(np.sum(matrix == 1)),
        "zeros_count": int(np.sum(matrix == 0)),
        "minus_ones_count": int(np.sum(matrix == -1)),
        "density": float(np.sum(matrix == 1) / matrix.size)
    }

def log_quasi_biclique_result(matrix: np.ndarray, rows: List[int], cols: List[int], success: bool, phase: str = ""):
    """Log quasi-biclique detection results."""
    global results_logger
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
        selected_matrix = matrix[np.ix_(rows, cols)]
        result["selected_density"] = float(np.sum(selected_matrix == 1) / selected_matrix.size)
        result["selected_ones"] = int(np.sum(selected_matrix == 1))
        result["selected_zeros"] = int(np.sum(selected_matrix == 0))
        
        # Add positions of zeros in the selected quasi-biclique
        zeros_positions = []
        for i, row_idx in enumerate(rows):
            for j, col_idx in enumerate(cols):
                if selected_matrix[i, j] == 0:
                    zeros_positions.append([row_idx, col_idx])
        
        result["zeros_positions"] = zeros_positions
        result["zeros_positions_count"] = len(zeros_positions)
        

        # Add positions of minus ones if present
        if np.any(selected_matrix == -1):
            minus_ones_positions = []
            for i, row_idx in enumerate(rows):
                for j, col_idx in enumerate(cols):
                    if selected_matrix[i, j] == -1:
                        minus_ones_positions.append([row_idx, col_idx])
            result["minus_ones_positions"] = minus_ones_positions
            result["minus_ones_positions_count"] = len(minus_ones_positions)
    
    results_logger["quasi_biclique_results"].append(result)

def log_clustering_step_result(reads1: List[int], reads0: List[int], cols: List[int], ilp_calls: int, step_num: int, input_matrix: np.ndarray = None):
    """Log clustering step results with density calculations."""
    global results_logger
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
    
    # Calculate densities if matrix is provided
    if input_matrix is not None and len(reads1) > 0 and len(reads0) > 0 and len(cols) > 0:
        # Extract submatrices for each group
        group1_matrix = input_matrix[np.ix_(reads1, cols)]
        group0_matrix = input_matrix[np.ix_(reads0, cols)]
        
        # Calculate densities for group 1
        group1_total = group1_matrix.size
        group1_ones = np.sum(group1_matrix == 1)
        group1_zeros = np.sum(group1_matrix == 0)
        group1_minus_ones = np.sum(group1_matrix == -1)
        
        step_result["group1_density"] = {
            "ones_count": int(group1_ones),
            "zeros_count": int(group1_zeros),
            "minus_ones_count": int(group1_minus_ones),
            "ones_density": float(group1_ones / group1_total) if group1_total > 0 else 0.0,
            "zeros_density": float(group1_zeros / group1_total) if group1_total > 0 else 0.0,
            "minus_ones_density": float(group1_minus_ones / group1_total) if group1_total > 0 else 0.0
        }
        
        # Calculate densities for group 0
        group0_total = group0_matrix.size
        group0_ones = np.sum(group0_matrix == 1)
        group0_zeros = np.sum(group0_matrix == 0)
        group0_minus_ones = np.sum(group0_matrix == -1)
        
        step_result["group0_density"] = {
            "ones_count": int(group0_ones),
            "zeros_count": int(group0_zeros),
            "minus_ones_count": int(group0_minus_ones),
            "ones_density": float(group0_ones / group0_total) if group0_total > 0 else 0.0,
            "zeros_density": float(group0_zeros / group0_total) if group0_total > 0 else 0.0,
            "minus_ones_density": float(group0_minus_ones / group0_total) if group0_total > 0 else 0.0
        }
        
        # Calculate contrast between groups
        step_result["density_contrast"] = {
            "ones_difference": step_result["group1_density"]["ones_density"] - step_result["group0_density"]["ones_density"],
            "zeros_difference": step_result["group1_density"]["zeros_density"] - step_result["group0_density"]["zeros_density"]
        }
    
    results_logger["clustering_steps"].append(step_result)

def save_results_to_file(filename: str):
    """Save all logged results to file."""
    global results_logger
    results_logger["timestamp"] = time.strftime("%Y-%m-%d %H:%M:%S")
    with open(filename, 'w') as f:
        json.dump(results_logger, f, indent=2)
    print(f"Results saved to {filename}")

def save_matrix_to_txt(matrix: np.ndarray, filename: str = "matrix_display.txt"):
    """Save matrix to text file for visualization."""
    try:
        with open(filename, 'w') as f:
            f.write(f"Matrix shape: {matrix.shape[0]} rows x {matrix.shape[1]} columns\n")
            f.write(f"Total elements: {matrix.size}\n")
            f.write(f"Ones count: {np.sum(matrix == 1)}\n")
            f.write(f"Zeros count: {np.sum(matrix == 0)}\n")
            f.write(f"Minus ones count: {np.sum(matrix == -1)}\n")
            f.write(f"Density: {np.sum(matrix == 1) / matrix.size:.4f}\n")
            f.write("\nMatrix content:\n")
            f.write("=" * 50 + "\n")
            
            # Write matrix with row and column indices
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

@print_decorator('clustering')
@track_clustering()
@timed_matrix_operation('full_clustering')
def clustering_full_matrix(
        input_matrix:np.ndarray, 
        regions :List[List[int]],
        steps : List[Tuple[List[int], List[int], List[int]]], 
        min_row_quality:int=5,
        min_col_quality:int = 3,
        error_rate : float = 0.025,
        filename: str = "unknown",
        haplotype_count: int = 0,
        log_results: bool = False,
        save_matrix_txt: bool = True
        ) -> List[Tuple[List[int], List[int], List[int]]]:
    """
    Perform exhaustive iterative biclustering on a binary matrix to extract all significant patterns.
    
    This function systematically processes predefined regions of a binary matrix to identify 
    all possible separations between rows based on their column patterns. It applies binary 
    clustering iteratively until no more significant patterns can be found.
    
    Parameters
    ----------
    input_matrix : np.ndarray
        Input binary matrix with values (0, 1) where rows and columns represent 
        data points and features respectively. Values indicate feature presence (1), 
        absence (0).
    regions : List[List[int]]
        List of column index groups to process separately. Each region contains 
        column indices that should be analyzed together as a coherent unit.
    steps : List[Tuple[List[int], List[int], List[int]]]
        Pre-existing clustering results to preserve. Each tuple contains
        (row_indices_group1, row_indices_group2, column_indices).
    min_row_quality : int, optional
        Minimum number of rows required for a cluster to be considered valid.
        Default is 5.
    min_col_quality : int, optional
        Minimum number of columns required for a region to be processed.
        Regions with fewer columns are skipped. Default is 3.
    error_rate : float, optional
        Tolerance level for pattern detection, allowing for noise and imperfections
        in the binary patterns. Default is 0.025 (2.5%).
    filename : str, optional
        Original filename for logging purposes (unused in processing).
    haplotype_count : int, optional
        Number of haplotypes for organization (unused in processing).
    
    Returns
    -------
    List[Tuple[List[int], List[int], List[int]]]
        Complete list of all valid clustering steps found. Each tuple contains:
        - [0] : List[int] - Row indices in first group (pattern match)
        - [1] : List[int] - Row indices in second group (pattern opposite)  
        - [2] : List[int] - Column indices where this separation is significant
        
        Only returns steps where both groups are non-empty and column count
        meets minimum quality threshold.
    
    Algorithm
    ---------
    1. **Initialization**: Start with existing clustering steps
    
    2. **Region Iteration**: Process each column region independently:
       - Skip regions with insufficient columns (< min_col_quality)
       - Initialize remaining columns for processing
       
    3. **Exhaustive Pattern Extraction**: For each region:
       - Apply binary clustering to find one significant row separation
       - Convert local column indices to global matrix coordinates
       - Save valid separations to results
       - Remove processed columns from remaining set
       - Continue until no significant patterns remain
    
    4. **Result Filtering**: Return only clustering steps that satisfy:
       - Both row groups contain at least one element
       - Column set meets minimum quality requirements
    
    See Also
    --------
    clustering_step : Single iteration of binary pattern detection
    find_quasi_biclique : Core optimization algorithm for dense pattern detection  
    pre_processing : Preprocessing step that identifies initial regions and steps
    post_processing : Final step that converts clustering steps to read groups
    
    """
    # Save matrix to text file at the beginning
    if save_matrix_txt:
        matrix_filename = f"matrix_{filename}_{time.strftime('%H%M%S')}.txt"
        save_matrix_to_txt(input_matrix, matrix_filename)
    
    # Log matrix info if requested
    if log_results:
        log_matrix_info(input_matrix)
    
    # Initialize result list with existing steps
    steps_result = steps.copy() if steps else []
    total_clustering_steps = 0
    total_ilp_calls = 0
    start_time = time.time()
    
    # Process each region if any regions are provided
    if len(regions) > 0:
        for idx, region in enumerate(regions):   
            # Initialize remaining columns for this region
            remain_cols = region
            status = True
            region_start_time = time.time()
            region_steps = 0
            region_ilp_calls = 0
            
            # Only process regions that meet minimum quality threshold
            if len(remain_cols) >= min_col_quality:
                logger.debug(f"Processing region {idx} with {len(remain_cols)} columns")  # GARDER CE LOG
                
                # Iteratively extract patterns until no more significant ones found
                while len(remain_cols) >= min_col_quality and status:
                    # Apply clustering to current remaining columns
                    sub_matrix = input_matrix[:, remain_cols]
                    reads1, reads0, cols, ilp_calls = clustering_step_with_stats(
                        sub_matrix, 
                        error_rate=error_rate,
                        min_row_quality=min_row_quality, 
                        min_col_quality=min_col_quality,
                        log_results=log_results
                    )
                    
                    region_ilp_calls += ilp_calls
                    
                    # Convert local column indices back to global matrix coordinates
                    cols = [remain_cols[c] for c in cols]
                    
                    # Check if valid pattern was found
                    if len(cols) == 0:
                        status = False  # No more patterns, stop processing this region
                    else:
                        # Save valid clustering step
                        steps_result.append((reads1, reads0, cols))
                        region_steps += 1
                        logger.debug(f"Found clustering step: {len(reads1)} vs {len(reads0)} reads, {len(cols)} columns")  # GARDER CE LOG
                        
                        # Log step if requested
                        if log_results:
                            log_clustering_step_result(reads1, reads0, cols, ilp_calls, region_steps, input_matrix)
                        
                        # Remove processed columns from remaining set
                        remain_cols = [c for c in remain_cols if c not in cols]
            else:
                logger.debug(f"Skipping region {idx}: insufficient columns ({len(remain_cols)} < {min_col_quality})")  # GARDER CE LOG
            
            total_clustering_steps += region_steps
            total_ilp_calls += region_ilp_calls
    
    # Log clustering results for debugging (GARDER CES LOGS)
    logger.debug(f"Clustering completed: {len(steps_result)} total steps found")
    
    # Filter and return only valid steps with non-empty groups and sufficient columns
    valid_steps = [step for step in steps_result if len(step[0]) > 0 and len(step[1]) > 0 and len(step[2]) >= min_col_quality]
    logger.info(f"Found {len(valid_steps)} valid clustering steps")  # GARDER CE LOG
    
    return valid_steps

def clustering_step_with_stats(input_matrix: np.ndarray,
                             error_rate: float = 0.025,
                             min_row_quality: int = 5,
                             min_col_quality: int = 3,
                             log_results: bool = False
                             ) -> Tuple[List[int], List[int], List[int], int]:
    """Version de clustering_step qui retourne aussi le nombre d'appels ILP."""
    ilp_calls = 0
    
    # Create binary matrices for pattern detection
    matrix1 = input_matrix.copy()
    matrix1[matrix1 == -1] = 0  # Convert missing values to 0 for positive patterns
    
    matrix0 = input_matrix.copy()
    matrix0[matrix0 == -1] = 1
    matrix0 = (matrix0 - 1) * -1  # Invert matrix for negative pattern detection
    
    logger.debug(f"Starting clustering step on {matrix1.shape[0]} rows, {matrix1.shape[1]} columns")  # GARDER CE LOG
    
    # Initialize tracking variables for iterative clustering
    remain_rows = range(matrix1.shape[0])  # All rows initially available
    current_cols = range(matrix1.shape[1])  # All columns initially available
    clustering_1 = True  # Alternate between positive (True) and negative (False) patterns
    status = True  # Continue while valid patterns are found
    rw1, rw0 = [], []  # Accumulate rows for positive and negative groups
    
    # Iteratively extract patterns until insufficient data remains
    iteration = 0
    start_time = time.time()
    while len(remain_rows) >= min_row_quality and len(current_cols) >= min_col_quality and status:
        iteration += 1
        logger.debug(f"Clustering iteration {iteration}: {len(remain_rows)} rows, {len(current_cols)} columns, pattern={'positive' if clustering_1 else 'negative'}")
        
        # Apply quasi-biclique detection on appropriate matrix
        if clustering_1:
            # Search for positive patterns (dense regions of 1s)
            rw, cl, status = find_quasi_biclique(matrix1[remain_rows][:, current_cols], error_rate, log_results, f"iteration_{iteration}_positive")
        else:
            # Search for negative patterns (dense regions of 0s in original)
            rw, cl, status = find_quasi_biclique(matrix0[remain_rows][:, current_cols], error_rate, log_results, f"iteration_{iteration}_negative")
        
        if status:
            ilp_calls += 1
             
        # Convert local indices back to global matrix coordinates
        rw = [remain_rows[r] for r in rw]  # Map row indices to original matrix
        cl = [current_cols[c] for c in cl]  # Map column indices to original matrix
        
        current_cols = cl  # Update working column set to detected significant columns
        
        # Accumulate rows into appropriate pattern group if valid pattern found
        if status and len(cl) > 0:
            if clustering_1:
                rw1.extend(rw)  # Add to positive pattern group
                logger.debug(f"Added {len(rw)} rows to positive group")
            else:
                rw0.extend(rw)  # Add to negative pattern group
                logger.debug(f"Added {len(rw)} rows to negative group")
        else:
            logger.debug(f"No valid pattern found in iteration {iteration}")
                
        # Remove processed rows from remaining set for next iteration
        remain_rows = [r for r in remain_rows if r not in rw]
        # Alternate pattern detection type for next iteration
        clustering_1 = not clustering_1
    
    # Log final clustering statistics (GARDER CE LOG)
    logger.debug(f"Clustering step completed: {len(rw1)} positive reads, {len(rw0)} negative reads, {len(current_cols)} columns")
    
    # Log densities if requested
    if log_results and len(rw1) > 0 and len(rw0) > 0 and len(current_cols) > 0:
        group1_matrix = input_matrix[np.ix_(rw1, current_cols)]
        group0_matrix = input_matrix[np.ix_(rw0, current_cols)]
        
        group1_ones_density = np.sum(group1_matrix == 1) / group1_matrix.size
        group1_zeros_density = np.sum(group1_matrix == 0) / group1_matrix.size
        group0_ones_density = np.sum(group0_matrix == 1) / group0_matrix.size
        group0_zeros_density = np.sum(group0_matrix == 0) / group0_matrix.size
        
        logger.debug(f"Group 1 densities: ones={group1_ones_density:.3f}, zeros={group1_zeros_density:.3f}")
        logger.debug(f"Group 0 densities: ones={group0_ones_density:.3f}, zeros={group0_zeros_density:.3f}")
    
    return rw1, rw0, current_cols, ilp_calls

@print_decorator('clustering')
@track_ilp_call('quasi_biclique_detection')
def find_quasi_biclique(
    input_matrix: np.ndarray,
    error_rate: float = 0.025,
    log_results: bool = False,
    phase: str = ""
) -> Tuple[List[int], List[int], bool]:
    """
    Find a quasi-biclique in a binary matrix using the max_Ones_comp function with a 
    seed-and-extend strategy.
    
    Parameters
    ----------
    input_matrix : np.ndarray
        Input binary matrix with values (0, 1)
    error_rate : float, optional
        Maximum fraction of 0s allowed in the quasi-biclique. Default is 0.025.
    log_results : bool, optional
        Whether to log detailed results. Default is False.
    phase : str, optional
        Name of the detection phase for logging. Default is "".
    
    Returns
    -------
    Tuple[List[int], List[int], bool]
        Triple containing:
        - [0] : List[int] - Row indices included in the quasi-biclique
        - [1] : List[int] - Column indices included in the quasi-biclique  
        - [2] : bool - Success status (True if valid solution found, False otherwise)
    """
    from density_test_max_one.max_one import max_Ones_comp
    
    # Copy input matrix to avoid modifying original
    X_problem = input_matrix.copy()
    
    # Get matrix dimensions
    n_rows, n_cols = X_problem.shape
    
    # Handle edge case: empty matrix
    if n_rows == 0 or n_cols == 0:
        logger.debug("Empty matrix provided to quasi-biclique detection")
        return [], [], False
    
    logger.debug(f"Starting quasi-biclique detection on {n_rows}x{n_cols} matrix")
    
    # Sort rows and columns by decreasing number of 1s (highest density first)
    cols_sorted = np.argsort(X_problem.sum(axis=0))[::-1]
    rows_sorted = np.argsort(X_problem.sum(axis=1))[::-1]
    
    # Phase 1: Seed region selection - find dense area to start optimization
    seed_rows = n_rows // 3
    seed_cols = n_cols // 3
    
    # Adjust step size based on matrix width
    step_n = 10 if n_cols > 50 else 2
    
    # Search for the largest sub-region with >99% density of 1s
    for x in range(n_rows // 3, n_rows, 10):
        for y in range(seed_cols, n_cols, step_n):
            nb_of_ones = 0
            for row in rows_sorted[:x]:
                for col in cols_sorted[:y]:
                    nb_of_ones += X_problem[row, col]
            ratio_ones = nb_of_ones / (x * y)
            if ratio_ones > 0.99:
                seed_rows = x
                seed_cols = y
    
    logger.debug(f"Using seed region: {seed_rows} rows x {seed_cols} columns")
    
    # Calculate degree for each row and column in the seed region
    seed_row_indices = rows_sorted[:seed_rows]
    seed_col_indices = cols_sorted[:seed_cols]
    
    # Calculate degree for seed rows/columns (number of 1s in each row/column)
    seed_matrix = X_problem[np.ix_(seed_row_indices, seed_col_indices)]
    row_degrees = np.sum(seed_matrix == 1, axis=1)
    col_degrees = np.sum(seed_matrix == 1, axis=0)
    
    # Create rows_data and cols_data for seed regions
    rows_data = [(i, int(row_degrees[idx])) for idx, i in enumerate(seed_row_indices)]
    cols_data = [(j, int(col_degrees[idx])) for idx, j in enumerate(seed_col_indices)]
    
    # Create edges: list of (row, col) pairs where the matrix has a 1
    edges = []
    for i_idx, i in enumerate(seed_row_indices):
        for j_idx, j in enumerate(seed_col_indices):
            if seed_matrix[i_idx, j_idx] == 1:
                edges.append((i, j))
    
    # Solve initial optimization problem with seed rows and columns
    try:
        with suppress_pulp_output():
            # Create the model
            model = max_Ones_comp(rows_data, cols_data, edges, error_rate)
            # Solve the model
            model.solve(plp.PULP_CBC_CMD(msg=0))
        
        # Extract seed solution results
        rw = []
        cl = []
        for var in model.variables():
            if var.name.startswith('row_') and var.varValue > 0.5:
                rw.append(int(var.name.split('_')[1]))
            elif var.name.startswith('col_') and var.varValue > 0.5:
                cl.append(int(var.name.split('_')[1]))
        
        logger.debug(f"Initial seed solution: {len(rw)} rows, {len(cl)} columns")
        
        # Phase 2: Row extension - add compatible rows
        rem_rows = [r for r in range(n_rows) if r not in seed_row_indices]
        if len(cl) > 0:
            # Find rows with >50% compatibility with selected columns
            rem_rows_sum = X_problem[rem_rows][:, cl].sum(axis=1)
            potential_rows = [r for idx, r in enumerate(rem_rows) 
                             if rem_rows_sum[idx] > 0.5 * len(cl)]
        else:
            potential_rows = []
        
        # If there are potential rows to add, create a new optimization
        if potential_rows and len(rw) > 0 and len(cl) > 0:
            logger.debug(f"Extending with {len(potential_rows)} compatible rows")
            
            # All rows to consider (selected rows + potential rows)
            extended_rows = rw + potential_rows
            
            # Calculate degrees for extended rows
            extended_row_degrees = np.sum(X_problem[extended_rows][:, cl] == 1, axis=1)
            
            # Create new data for extended optimization
            extended_rows_data = [(r, int(extended_row_degrees[idx])) for idx, r in enumerate(extended_rows)]
            extended_cols_data = [(c, int(np.sum(X_problem[extended_rows, c] == 1))) for c in cl]
            
            # Create edges for extended optimization
            extended_edges = []
            for i in extended_rows:
                for j in cl:
                    if X_problem[i, j] == 1:
                        extended_edges.append((i, j))
            
            # Run extended optimization with rows
            with suppress_pulp_output():
                extended_model = max_Ones_comp(extended_rows_data, extended_cols_data, extended_edges, error_rate)
                extended_model.solve(plp.PULP_CBC_CMD(msg=0))
            
            # Extract results after row extension
            rw = []
            cl = []
            for var in extended_model.variables():
                if var.name.startswith('row_') and var.varValue > 0.5:
                    rw.append(int(var.name.split('_')[1]))
                elif var.name.startswith('col_') and var.varValue > 0.5:
                    cl.append(int(var.name.split('_')[1]))
        
        # Phase 3: Column extension - add compatible columns
        rem_cols = [c for c in range(n_cols) if c not in cl]
        if len(rw) > 0:
            # Find columns with >90% compatibility with selected rows
            rem_cols_sum = X_problem[rw][:, rem_cols].sum(axis=0)
            potential_cols = [c for idx, c in enumerate(rem_cols) 
                             if rem_cols_sum[idx] > 0.9 * len(rw)]
        else:
            potential_cols = []
        
        # If there are potential columns to add, create a final optimization
        if potential_cols and len(rw) > 0:
            logger.debug(f"Extending with {len(potential_cols)} compatible columns")
            
            # All columns to consider (selected columns + potential columns)
            final_cols = cl + potential_cols
            
            # Calculate degrees for final optimization
            final_row_degrees = np.sum(X_problem[rw][:, final_cols] == 1, axis=1)
            final_col_degrees = np.sum(X_problem[rw][:, final_cols] == 1, axis=0)
            
            # Create data for final optimization
            final_rows_data = [(r, int(final_row_degrees[idx])) for idx, r in enumerate(rw)]
            final_cols_data = [(c, int(np.sum(X_problem[rw, c] == 1))) for c in final_cols]
            
            # Create edges for final optimization
            final_edges = []
            for i in rw:
                for j in final_cols:
                    if X_problem[i, j] == 1:
                        final_edges.append((i, j))
            
            # Run final optimization
            with suppress_pulp_output():
                final_model = max_Ones_comp(final_rows_data, final_cols_data, final_edges, error_rate)
                final_model.solve(plp.PULP_CBC_CMD(msg=0))
            
            # Extract final results
            rw = []
            cl = []
            for var in final_model.variables():
                if var.name.startswith('row_') and var.varValue > 0.5:
                    rw.append(int(var.name.split('_')[1]))
                elif var.name.startswith('col_') and var.varValue > 0.5:
                    cl.append(int(var.name.split('_')[1]))
        
        logger.debug(f"Final quasi-biclique: {len(rw)} rows, {len(cl)} columns")
        
        # Check if we found a valid solution
        if len(rw) > 0 and len(cl) > 0:
            # Calculate density of final selection
            selected = X_problem[np.ix_(rw, cl)]
            density = np.sum(selected == 1) / selected.size
            logger.debug(f"Final density: {density:.4f}")
            
            # Log results if requested
            if log_results:
                log_quasi_biclique_result(input_matrix, rw, cl, True, phase)
            
            return rw, cl, True
        else:
            # Log failed attempt if requested
            if log_results:
                log_quasi_biclique_result(input_matrix, [], [], False, phase)
            return [], [], False
        
    except Exception as e:
        logger.error(f"Error solving optimization problem: {e}")
        # Log failed attempt if requested
        if log_results:
            log_quasi_biclique_result(input_matrix, [], [], False, phase)
        return [], [], False