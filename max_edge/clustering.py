from typing import List, Tuple
import numpy as np
import gurobipy as grb
import logging
import os
import sys
import time
import contextlib
from io import StringIO
from pathlib import Path
import json
# Import max_e_r functions
from max_e_r import max_e_r, max_e_wr
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
os.environ['GRB_LICENSE_FILE'] = ''  # Prevent license file messages
gurobi_logger = logging.getLogger('gurobipy')
gurobi_logger.setLevel(logging.CRITICAL)
gurobi_logger.propagate = False

options = {
	"WLSACCESSID":"af4b8280-70cd-47bc-aeef-69ecf14ecd10",
	"WLSSECRET":"04da6102-8eb3-4e38-ba06-660ea8f87bf2",
	"LICENSEID":2669217
}

@contextlib.contextmanager
def suppress_gurobi_output():
    """Context manager to completely suppress Gurobi output"""
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
        "selected_rows": [int(r) for r in rows],  # Convert to int
        "selected_cols": [int(c) for c in cols],  # Convert to int
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
                    zeros_positions.append([int(row_idx), int(col_idx)])  # Convert to int
        
        result["zeros_positions"] = zeros_positions
        result["zeros_positions_count"] = len(zeros_positions)
        

        # Add positions of minus ones if present
        if np.any(selected_matrix == -1):
            minus_ones_positions = []
            for i, row_idx in enumerate(rows):
                for j, col_idx in enumerate(cols):
                    if selected_matrix[i, j] == -1:
                        minus_ones_positions.append([int(row_idx), int(col_idx)])  # Convert to int
            result["minus_ones_positions"] = minus_ones_positions
            result["minus_ones_positions_count"] = len(minus_ones_positions)
    
    results_logger["quasi_biclique_results"].append(result)

def log_clustering_step_result(reads1: List[int], reads0: List[int], cols: List[int], ilp_calls: int, step_num: int, input_matrix: np.ndarray = None):
    """Log clustering step results with density calculations."""
    global results_logger
    step_result = {
        "step_number": step_num,
        "reads_group1": [int(r) for r in reads1],  # Convert to int
        "reads_group0": [int(r) for r in reads0],  # Convert to int
        "columns": [int(c) for c in cols],  # Convert to int
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
            result = find_quasi_biclique(matrix1[remain_rows][:, current_cols], error_rate, log_results, f"iteration_{iteration}_positive")
        else:
            result = find_quasi_biclique(matrix0[remain_rows][:, current_cols], error_rate, log_results, f"iteration_{iteration}_negative")
        
        # Defensive unpacking to avoid "cannot unpack non-iterable int object"
        if not (isinstance(result, tuple) and len(result) == 3):
            logger.warning(f"find_quasi_biclique returned unexpected result: {result}")
            rw, cl, status = [], [], False
        else:
            rw, cl, status = result

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
    Find a quasi-biclique in a binary matrix using integer linear programming optimization.
    
    This function identifies the largest dense sub-matrix where most elements are 1s,
    with tolerance for noise defined by the error rate. It uses a three-phase approach:
    seeding with a high-density region, then iteratively extending by rows and columns
    to maximize the objective while maintaining density constraints.
    
    Parameters
    ----------
    input_matrix : np.ndarray
        Input binary matrix with values (0, 1) where rows and columns represent 
        data points and features respectively. Values indicate feature presence (1) 
        or absence (0).
    error_rate : float, optional
        Maximum fraction of 0s allowed in the quasi-biclique, defining noise tolerance.
        A value of 0.025 means up to 2.5% of cells can be 0s. Default is 0.025.
    
    Returns
    -------
    Tuple[List[int], List[int], bool]
        Triple containing the quasi-biclique results:
        - [0] : List[int] - Row indices included in the quasi-biclique
        - [1] : List[int] - Column indices included in the quasi-biclique  
        - [2] : bool - Success status (True if valid solution found, False otherwise)
        
        Empty lists are returned when no significant quasi-biclique is found or
        optimization fails.
    
    Algorithm
    ---------
    1. **Preprocessing and Sorting**:
       - Sort rows and columns by decreasing number of 1s
       - Handle edge cases (empty matrix)
       
    2. **Seed Region Selection**:
       - Start with top 1/3 of rows and columns by density
       - Use max_e_r with delta=0 for initial seeding
       
    3. **Phase 1 - Row Extension**:
       - Identify remaining rows with >50% compatibility with selected columns
       - Use max_e_wr to extend with compatible rows
       
    4. **Phase 2 - Column Extension**:
       - Identify remaining columns with >90% compatibility with selected rows
       - Use max_e_wr for final optimization
       
    5. **Result Extraction**:
       - Extract selected rows and columns from optimal solution
       - Return success status based on optimization outcome
    
    """
    # Copy input matrix to avoid modifying original
    X_problem = input_matrix.copy()

    # Sort rows and columns by decreasing number of 1s (highest density first)
    cols_sorted = np.argsort(X_problem.sum(axis=0))[::-1]
    rows_sorted = np.argsort(X_problem.sum(axis=1))[::-1]

    m = len(rows_sorted)
    n = len(cols_sorted)

    # Handle edge case: empty matrix
    if m == 0 or n == 0:
        logger.debug("Empty matrix provided to quasi-biclique detection")
        return [], [], False

    logger.debug(f"Starting quasi-biclique detection on {m}x{n} matrix")
    
    # Helper function to convert matrix to edge format
    def matrix_to_edges(matrix, rows, cols):
        """Convert matrix subset to edge format required by max_e_r"""
        edges = []
        rows_data = []
        cols_data = []
        
        # Create rows_data with degrees
        for i, row in enumerate(rows):
            degree = np.sum(matrix[row, cols])
            rows_data.append((i, degree))  # Use local indices
        
        # Create cols_data with degrees  
        for j, col in enumerate(cols):
            degree = np.sum(matrix[rows, col])
            cols_data.append((j, degree))  # Use local indices
        
        # Create edges list
        for i, row in enumerate(rows):
            for j, col in enumerate(cols):
                if matrix[row, col] == 1:
                    edges.append((i, j))  # Use local indices
        
        return rows_data, cols_data, edges

    # Phase 1: Find seed using max_e_r with delta=0
    rows_data, cols_data, edges = matrix_to_edges(X_problem, rows_sorted, cols_sorted)

    if not edges:
        logger.debug("No edges found in matrix")
        return [], [], False

    # Initialize final result variables
    final_rw = []
    final_cl = []

    try:
        # Create and solve max_e_r model for seeding with delta=0 (perfect density)
        with suppress_gurobi_output():
            model = max_e_r(rows_data, cols_data, edges, delta=0.0, debug=0)
            if model is None:
                logger.debug("Failed to create max_e_r model")
                return [], [], False
            
            model.Params.OutputFlag = 0
            model.Params.LogToConsole = 0
            model.Params.MIPGAP = 0.05
            model.Params.TimeLimit = 20
            model.optimize()

        if model.Status not in [grb.GRB.OPTIMAL, grb.GRB.TIME_LIMIT]:
            logger.debug(f"Seed optimization failed with status {model.Status}")
            return [], [], False

        # Extract seed solution
        seed_rw = []
        seed_cl = []
        for var in model.getVars():
            if var.VarName.startswith('row_') and var.X > 0.5:
                local_idx = int(var.VarName.split('_')[1])
                seed_rw.append(local_idx)
            elif var.VarName.startswith('col_') and var.X > 0.5:
                local_idx = int(var.VarName.split('_')[1])
                seed_cl.append(local_idx)

        prev_obj = int(model.ObjVal) if model.Status == grb.GRB.OPTIMAL else 0
        logger.debug(f"Initial seed solution: {len(seed_rw)} rows, {len(seed_cl)} columns, obj={prev_obj}")

        if not seed_rw or not seed_cl:
            logger.debug("Empty seed solution found")
            return [], [], False

    except Exception as e:
        logger.warning(f"Failed to solve seed optimization: {e}")
        return [], [], False

    # Phase 2: Extend with all rows and columns using max_e_wr
    try:
        with suppress_gurobi_output():
            model = max_e_wr(rows_data, cols_data, edges, seed_rw, seed_cl, prev_obj, delta=error_rate, debug=0)
            if model is None:
                logger.debug("Failed to create max_e_wr model for extension")
                # Fall back to seed solution
                final_rw = [rows_sorted[i] for i in seed_rw]
                final_cl = [cols_sorted[i] for i in seed_cl]
            else:
                model.Params.OutputFlag = 0
                model.Params.LogToConsole = 0
                model.Params.MIPGAP = 0.05
                model.Params.TimeLimit = 30
                model.optimize()

                if model.Status in [grb.GRB.OPTIMAL, grb.GRB.TIME_LIMIT]:
                    # Extract final solution
                    for var in model.getVars():
                        if var.VarName.startswith('row_') and var.X > 0.5:
                            local_idx = int(var.VarName.split('_')[1])
                            final_rw.append(rows_sorted[local_idx])
                        elif var.VarName.startswith('col_') and var.X > 0.5:
                            local_idx = int(var.VarName.split('_')[1])
                            final_cl.append(cols_sorted[local_idx])
                    
                    logger.debug(f"Extension completed: {len(final_rw)} rows, {len(final_cl)} columns")
                else:
                    logger.debug(f"Extension optimization failed with status {model.Status}")
                    # Fall back to seed solution
                    final_rw = [rows_sorted[i] for i in seed_rw]
                    final_cl = [cols_sorted[i] for i in seed_cl]

    except Exception as e:
        logger.warning(f"Extension failed: {e}")
        # Fall back to seed solution
        final_rw = [rows_sorted[i] for i in seed_rw]
        final_cl = [cols_sorted[i] for i in seed_cl]

    logger.debug(f"Quasi-biclique found: {len(final_rw)} rows, {len(final_cl)} columns")
    
    # Log results if requested
    if log_results:
        log_quasi_biclique_result(input_matrix, final_rw, final_cl, True, phase)
           
    return final_rw, final_cl, True
