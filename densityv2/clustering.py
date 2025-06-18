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
import itertools

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
        haplotype_count: int = 0
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
                        min_col_quality=min_col_quality
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
            rw, cl, status = find_quasi_biclique(matrix1[remain_rows][:, current_cols], error_rate)
        else:
            # Search for negative patterns (dense regions of 0s in original)
            rw, cl, status = find_quasi_biclique(matrix0[remain_rows][:, current_cols], error_rate)
        
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
    return rw1, rw0, current_cols, ilp_calls

debug = 0

@print_decorator('clustering')
@track_ilp_call('quasi_biclique_detection')
def find_quasi_biclique(
    input_matrix: np.ndarray,
    error_rate: float = 0.025
) -> Tuple[List[int], List[int], bool]:
    """
    Find a quasi-biclique in a binary matrix using the max_e_wr integer linear programming approach.
    
    This function identifies the largest dense sub-matrix where most elements are 1s,
    with tolerance for noise defined by the error rate. It uses the improved model from
    Chang et al. with iterative warm restart optimization using Gurobi.
    
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
    """
    # Handle edge case: empty matrix
    m, n = input_matrix.shape
    if m == 0 or n == 0:
        logger.debug("Empty matrix provided to quasi-biclique detection")
        return [], [], False

    logger.debug(f"Starting quasi-biclique detection on {m}x{n} matrix using max_e_wr")
    
    # Prepare data for the max_e_wr model
    rows_data = [(i, int(input_matrix[i, :].sum())) for i in range(m)]
    cols_data = [(j, int(input_matrix[:, j].sum())) for j in range(n)]
    edges = [(i, j) for i in range(m) for j in range(n) if input_matrix[i, j] == 1]
    
    if not edges:
        logger.debug("No edges (1s) found in matrix")
        return [], [], False
    
    delta = error_rate
    
    try:
        # Phase 1: Get initial solution with max_e_r
        initial_model = max_e_r(rows_data, cols_data, edges, delta)
        
        with suppress_gurobi_output():
            initial_model.optimize()
        
        if initial_model.Status != grb.GRB.OPTIMAL:
            logger.warning(f"Initial optimization failed with status: {initial_model.Status}")
            return [], [], False
        
        # Extract initial solution
        current_rows = []
        current_cols = []
        current_obj = initial_model.ObjVal
        
        for var in initial_model.getVars():
            if var.VarName.startswith('row_') and var.X > 0.5:
                row_idx = int(var.VarName.split('_')[1])
                current_rows.append(row_idx)
            elif var.VarName.startswith('col_') and var.X > 0.5:
                col_idx = int(var.VarName.split('_')[1])
                current_cols.append(col_idx)
        
        logger.debug(f"Initial solution: {len(current_rows)} rows, {len(current_cols)} columns, obj={current_obj}")
        
        # Phase 2: Iterative improvement with max_e_wr
        max_iterations = 5
        iteration = 0
        improved = True
        
        while improved and iteration < max_iterations:
            iteration += 1
            logger.debug(f"Improvement iteration {iteration}")
            
            # Create warm restart model
            wr_model = max_e_wr(rows_data, cols_data, edges, current_rows, current_cols, current_obj, delta)
            
            with suppress_gurobi_output():
                wr_model.optimize()
            
            if wr_model.Status != grb.GRB.OPTIMAL:
                logger.debug(f"Warm restart iteration {iteration} failed, keeping previous solution")
                improved = False
                continue
            
            # Extract improved solution
            new_rows = []
            new_cols = []
            new_obj = wr_model.ObjVal
            
            for var in wr_model.getVars():
                if var.VarName.startswith('row_') and var.X > 0.5:
                    row_idx = int(var.VarName.split('_')[1])
                    new_rows.append(row_idx)
                elif var.VarName.startswith('col_') and var.X > 0.5:
                    col_idx = int(var.VarName.split('_')[1])
                    new_cols.append(col_idx)
            
            # Check if improvement was achieved
            if new_obj > current_obj:
                current_rows = new_rows
                current_cols = new_cols
                current_obj = new_obj
                logger.debug(f"Improved solution: {len(current_rows)} rows, {len(current_cols)} columns, obj={current_obj}")
            else:
                improved = False
                logger.debug(f"No improvement in iteration {iteration}, stopping")
        
        logger.debug(f"Final max_e_wr solution: {len(current_rows)} rows, {len(current_cols)} columns")
        return current_rows, current_cols, True
        
    except Exception as e:
        logger.error(f"Error in max_e_wr optimization: {e}")
        return [], [], False

def max_e_r(rows_data, cols_data, edges, delta):
    """
    Arguments:
    ----------
    rows_data: list of tuples (row, degree) of rows in the matrix.
    cols_data: list of tuples (col, degree) of columns in the matrix.
    edges: list of tuples (row, col) corresponding to the ones of the matrix.
    delta: float, error tolerance for density constraints.

    Returns:
    --------
    grb.Model:
    The ILP model.
    
    This is an improved model from Chang et al. 
    """
    # Initialize Gurobi environment
    with suppress_gurobi_output():
        env = grb.Env(params=options)
        env.setParam('OutputFlag', 0)
        env.setParam('LogToConsole', 0)
        env.setParam('LogFile', "")
        model = grb.Model('max_e_r', env=env)
        model.Params.OutputFlag = 0
        model.Params.LogToConsole = 0
        model.Params.LogFile = ""
        model.Params.MIPGAP = 0.05
        model.Params.TimeLimit = 20

    # Variables for rows and columns
    lpRows = {}
    for row, degree in rows_data:
        lpRows[row] = (model.addVar(lb=0, ub=1, vtype=grb.GRB.INTEGER, name=f'row_{row}'), degree)
    
    lpCols = {}
    for col, degree in cols_data:
        lpCols[col] = (model.addVar(lb=0, ub=1, vtype=grb.GRB.INTEGER, name=f'col_{col}'), degree)
    
    lpCells = {}
    for row, col in edges:
        lpCells[(row, col)] = model.addVar(lb=0, ub=1, vtype=grb.GRB.CONTINUOUS, name=f'cell_{row}_{col}')

    # Objective: maximize sum of selected cells
    model.setObjective(grb.quicksum(lpCells.values()), grb.GRB.MAXIMIZE)
    
    # Constraints for row and column thresholds
    row_threshold = 1
    col_threshold = 1
    
    if debug >= 3:
        print()
        print('-' * 40) 
        print(f"******** Solving model ******** max_e_r with delta = {delta}")
        print(f' # rows_data = {len(rows_data)}, # cols_data = {len(cols_data)}, # edges = {len(edges)}') 
        print(f'row_threshold = {row_threshold}')
        print(f'col_threshold = {col_threshold}')
        print()
        print('-' * 40)
        
    model.addConstr(grb.quicksum(lpvar for lpvar, _ in lpRows.values()) >= row_threshold, "row_threshold")
    model.addConstr(grb.quicksum(lpvar for lpvar, _ in lpCols.values()) >= col_threshold, "col_threshold")

    # Constraints for matrix structure
    for row, col in edges:
        if row in lpRows and col in lpCols and (row, col) in lpCells:
            model.addConstr(lpRows[row][0] >= lpCells[(row, col)], f'cell_{row}_{col}_1')
            model.addConstr(lpCols[col][0] >= lpCells[(row, col)], f'cell_{row}_{col}_2')
        else:
            if debug >= 1:
                print(f"Warning: row={row}, col={col} not found in lpRows/lpCols!")

    # Add density constraints
    __row_density_gurobi(rows_data, cols_data, edges, model, lpRows, lpCols, delta)
    __col_density_gurobi(rows_data, cols_data, edges, model, lpRows, lpCols, delta)

    return model 

def max_e_wr(rows_data, cols_data, edges, rows_res, cols_res, prev_obj, delta):
    """
    Arguments:
    ----------
    rows_data: list of tuples (row, degree) of rows in the matrix.
    cols_data: list of tuples (col, degree) of columns in the matrix.
    edges: list of tuples (row, col) corresponding to the ones of the matrix.
    rows_res : list of indices of rows (input solution)
    cols_res : list of indices of cols (input solution)
    prev_obj: previous objective value
    delta: float, error tolerance for density constraints.

    Returns:
    --------
    grb.Model:
    The ILP model.
    """
    # Initialize Gurobi environment
    with suppress_gurobi_output():
        env = grb.Env(params=options)
        env.setParam('OutputFlag', 0)
        env.setParam('LogToConsole', 0)
        env.setParam('LogFile', "")
        model = grb.Model('max_e_wr', env=env)
        model.Params.OutputFlag = 0
        model.Params.LogToConsole = 0
        model.Params.LogFile = ""
        model.Params.MIPGAP = 0.05
        model.Params.TimeLimit = 20

    # Variables for rows and columns
    lpRows = {}
    for row, degree in rows_data:
        lpRows[row] = (model.addVar(lb=0, ub=1, vtype=grb.GRB.INTEGER, name=f'row_{row}'), degree)
    
    lpCols = {}
    for col, degree in cols_data:
        lpCols[col] = (model.addVar(lb=0, ub=1, vtype=grb.GRB.INTEGER, name=f'col_{col}'), degree)
    
    lpCells = {}
    for row, col in edges:
        lpCells[(row, col)] = model.addVar(lb=0, ub=1, vtype=grb.GRB.CONTINUOUS, name=f'cell_{row}_{col}')

    cols_res_set = set(map(int, cols_res)) 
    rows_res_set = set(map(int, rows_res))

    if debug >= 2:
        print(f" !!!!!!!!!!!!!!!!!! I got a lower bound {prev_obj}")

    # Warm start initialization
    for row, (var, _) in lpRows.items():
        if row in rows_res_set:
            var.Start = 1
        else:
            var.Start = 0

    for col, (var, _) in lpCols.items():
        if col in cols_res_set:
            var.Start = 1
        else:
            var.Start = 0

    for (row, col), var in lpCells.items():
        if row in rows_res_set and col in cols_res_set:
            var.Start = 1
        else:
            var.Start = 0
            
    if debug >= 3:
        print("\n Initial point before solving:")
        print("ROWS:", {r: (1 if r in rows_res_set else 0) for r in [r for r, _ in rows_data]})
        print("COLS:", {c: (1 if c in cols_res_set else 0) for c in [c for c, _ in cols_data]})

    # Objective 
    model.setObjective(grb.quicksum(lpCells.values()), grb.GRB.MAXIMIZE)

    row_threshold = 2
    col_threshold = 2

    if debug >= 3:
        print()
        print('-' * 40) 
        print(f"******** Solving model ******** max_e_wr with delta = {delta}")
        print(f"# rows_data = {len(rows_data)}, # cols_data = {len(cols_data)}, # edges = {len(edges)}") 
        print(f"row_threshold = {row_threshold}")
        print(f"col_threshold = {col_threshold}")
        print()
        print('-' * 40)

    # Add constraints for row and column thresholds
    model.addConstr(grb.quicksum(lpvar for lpvar, _ in lpRows.values()) >= row_threshold, "row_threshold")
    model.addConstr(grb.quicksum(lpvar for lpvar, _ in lpCols.values()) >= col_threshold, "col_threshold")

    # Constraint for objective improvement
    model.addConstr(grb.quicksum(lpCells.values()) >= prev_obj + 1, "improvement")

    # Constraints for matrix structure
    for row, col in edges:
        if (row, col) in lpCells:
            model.addConstr(lpRows[row][0] >= lpCells[(row, col)], f'cell_{row}_{col}_1')
            model.addConstr(lpCols[col][0] >= lpCells[(row, col)], f'cell_{row}_{col}_2')
        else:
            if debug >= 1:
                print(f"Warning: ({row}, {col}) not found in lpCells!")

    __row_density_gurobi(rows_data, cols_data, edges, model, lpRows, lpCols, delta)
    __col_density_gurobi(rows_data, cols_data, edges, model, lpRows, lpCols, delta)

    return model 

def __col_density_gurobi(rows_data, cols_data, edges, model, lpRows, lpCols, delta):
    """
    Adds col density constraints to the Gurobi model.

    Arguments:
    ----------
    rows_data: list of tuples (row, degree) of rows in the matrix.
    cols_data: list of tuples (col, degree) of columns in the matrix.
    edges: list of tuples (row, col) corresponding to the ones of the matrix.
    model: Gurobi Model to add constraints to.
    lpRows: dict of row variables and their degrees.
    lpCols: dict of column variables and their degrees.
    delta: float, error tolerance for density constraints. 
    """
    Big_R = len(rows_data) + 1
    Big_C = len(cols_data) + 1
    Big_M = Big_R + Big_C 
    
    for col, _ in cols_data:
        col_edges = [u for u, v in edges if v == col]
        if col_edges:
            model.addConstr(
                grb.quicksum(lpRows[row][0] for row in col_edges) - (1 - delta) * grb.quicksum(lpRows[row][0] for row, _ in rows_data) >= 
                (lpCols[col][0] - 1) * Big_M,
                f"col_err_rate_1_{col}"
            )

def __row_density_gurobi(rows_data, cols_data, edges, model, lpRows, lpCols, delta):
    """
    Adds row density constraints to the Gurobi model.

    Arguments:
    ----------
    rows_data: list of tuples (row, degree) of rows in the matrix.
    cols_data: list of tuples (col, degree) of columns in the matrix.
    edges: list of tuples (row, col) corresponding to the ones of the matrix.
    model: Gurobi Model to add constraints to.
    lpRows: dict of row variables and their degrees.
    lpCols: dict of column variables and their degrees.
    delta: float, error tolerance for density constraints.
    """
    Big_R = len(rows_data) + 1
    Big_C = len(cols_data) + 1
    Big_M = Big_R + Big_C

    for row, _ in rows_data:
        row_edges = [v for u, v in edges if u == row]
        if row_edges:
            model.addConstr(
                grb.quicksum(lpCols[col][0] for col in row_edges) - (1 - delta) * grb.quicksum(lpCols[col][0] for col, _ in cols_data) >= 
                (lpRows[row][0] - 1) * Big_M,
                f"row_err_rate_1_{row}"
            )
        if row_edges:
            model.addConstr(
                grb.quicksum(lpCols[col][0] for col in row_edges) - (1 - delta) * grb.quicksum(lpCols[col][0] for col, _ in cols_data) >= 
                (lpRows[row][0] - 1) * Big_M,
                f"row_err_rate_1_{row}"
            )
