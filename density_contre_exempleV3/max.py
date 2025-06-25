import pulp as pl
from pulp import LpProblem, LpMaximize, LpVariable, lpSum, LpStatusOptimal
import itertools

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
    LpProblem:
    The ILP model.
    
    This is an improved model from Chnag et al. 
    """
    model = LpProblem(name="max_e_r", sense=LpMaximize)

    # Variables for rows and columns
    lpRows = {row: (LpVariable(f'row_{row}', cat='Integer',
                              lowBound=0, upBound=1), degree) for row, degree in rows_data}
    lpCols = {col: (LpVariable(f'col_{col}', cat='Integer',
                              lowBound=0, upBound=1), degree) for col, degree in cols_data}
    lpCells = {(row, col): LpVariable(
        #f'cell_{row}_{col}', cat="Binary")
        f'cell_{row}_{col}', cat='Continuous', lowBound=0, upBound=1) 
        for (row, col) in edges
    }

    #
    model += lpSum(lpCells), 'maximize_sum'
    # Constraints for row and column thresholds
    row_threshold = 1
    col_threshold = 1
   
    model += lpSum(lpvar for lpvar, _ in lpRows.values()) >= row_threshold, "row_threshold"
    model += lpSum(lpvar for lpvar, _ in lpCols.values()) >= col_threshold, "col_threshold"
    
    # Constraints for matrix structure
    # Check for missing columns in cols_data
    missing_cols = {col for _, col in edges} - {col for col, _ in cols_data}
    if missing_cols:
        print(f"Warning: These columns are missing from cols_data but appear in edges: {missing_cols}")

    for row, col in edges:
        if row in lpRows and col in lpCols and (row, col) in lpCells:
            model += (lpRows[row][0] >= lpCells[(row, col)]), f'cell_{row}_{col}_1'
            model += (lpCols[col][0] >= lpCells[(row, col)]), f'cell_{row}_{col}_2'
        else:
            print(f"Warning: row={row}, col={col} not found in lpRows/lpCols!")

    # Add row density constraints
    __row_density(rows_data, cols_data, edges, model, lpRows, lpCols, delta)
    __col_density(rows_data, cols_data, edges, model, lpRows, lpCols, delta)

    return model 


def max_e_wr(rows_data, cols_data, edges, rows_res, cols_res, prev_obj, delta):
    """
    Arguments:
    ----------
    rows_data: list of tuples (row, degree) of rows in the matrix.
    cols_data: list of tuples (col, degree) of columns in the matrix.
    edges: list of tuples (row, col) corresponding to the ones of the matrix.
    rows_res : list of indexies of rows (input solution)
    cols_res : list of indexies of cols (input solution)
    delta: float, error tolerance for density constraints.

    Returns:
    --------
    LpProblem:
    The ILP model.
    
    """
    model = LpProblem(name="max_e_wr", sense=LpMaximize)

    # Variables for rows and columns
    lpRows = {row: (LpVariable(f'row_{row}', cat='Integer',
                              lowBound=0, upBound=1), degree) for row, degree in rows_data}
    lpCols = {col: (LpVariable(f'col_{col}', cat='Integer',
                              lowBound=0, upBound=1), degree) for col, degree in cols_data}
    lpCells = {(row, col): LpVariable(
        #f'cell_{row}_{col}', cat="Binary")
        f'cell_{row}_{col}', cat='Continuous', lowBound=0, upBound=1) 
        for (row, col) in edges
    }

    cols_res_set = set(map(int, cols_res)) 
    rows_res_set = set(map(int, rows_res))

    

    # Warm start (PuLP way: use setInitialValue if needed)
    # Set and store initial values
    cells_indices = []
    row_initial_values = {}
    for row, (var, _) in lpRows.items():
        if row in rows_res_set:
            val = 1 
            var.setInitialValue(val)
            row_initial_values[row] = val # Store for later use
        else:
            val = 0
            row_initial_values[row] = val # Store for later use 

    col_initial_values = {}
    for col, (var, _) in lpCols.items():
        if col in cols_res_set:
            val = 1 
            var.setInitialValue(val)
            col_initial_values[col] = val
        else: 
            val = 0
            col_initial_values[col] = val

    cell_initial_values = {}
    for (row, col), var in lpCells.items():
        if (row, col) in itertools.product(row_initial_values, col_initial_values):
            var.setInitialValue(1)
            val = 1 
            cell_initial_values[(row, col)] = val
        else:
            var.setInitialValue(0)
            val = 0
            col_initial_values[col] = val 
    
    
    # Objective 
    model += lpSum(lpCells), 'current value'

    row_threshold = 2
    col_threshold = 2

    

    # Add constraints for row and column thresholds
    model += lpSum(lpvar for lpvar, _ in lpRows.values()) >= row_threshold, "row_threshold"
    model += lpSum(lpvar for lpvar, _ in lpCols.values()) >= col_threshold, "col_threshold"

    # Constraint for objective improvement
    model += lpSum(lpCells.values()) >= prev_obj + 1, "improvement"

    # Constraints for matrix structure
    for row, col in edges:
        if (row, col) in lpCells:
            model += (lpRows[row][0] >= lpCells[(row, col)]), f'cell_{row}_{col}_1'
            model += (lpCols[col][0] >= lpCells[(row, col)]), f'cell_{row}_{col}_2'
        else:
            print(f"Warning: ({row}, {col}) not found in lpCells!")

    #__row_density(rows_data, cols_data, edges, model, lpRows, lpCols, delta)
    __col_density(rows_data, cols_data, edges, model, lpRows, lpCols, delta)

    return model 


def __col_density(rows_data, cols_data, edges, model, lpRows, lpCols, delta):
    """
    Adds col density constraints to the model.

    Arguments:
    ----------
    rows_data: list of tuples (row, degree) of rows in the matrix.
    cols_data: list of tuples (col, degree) of columns in the matrix.
    edges: list of tuples (row, col) corresponding to the zeros of the matrix.
    model: LpProblem to add constraints to.
    lpRows: dict of row variables and their degrees.
    lpCols: dict of column variables and their degrees.
    delta: float, error tolerance for density constraints. 
    Attention : in the text it is written as ALPHA TO DO: (ALPHA : = 1-DELTA) !!!!
    """
    mu = 0.0001
    Big_R = len(rows_data) + 1
    Big_C = len(cols_data) + 1
    Big_M= Big_R+Big_C 
    for col, _ in cols_data:
        col_edges = [u for u, v in edges if v == col]
        model += (
            lpSum(lpRows[row][0] for row in col_edges) - (1 - delta) * lpSum(lpRows[row][0] for row, _ in rows_data) >= 
            (lpCols[col][0]-1) * Big_C
        ), f"col_err_rate_1_{col}"


def __row_density(rows_data, cols_data, edges, model, lpRows, lpCols, delta):
    """
    Adds row density constraints to the model.

    Arguments:
    ----------
    rows_data: list of tuples (row, degree) of rows in the matrix.
    cols_data: list of tuples (col, degree) of columns in the matrix.
    edges: list of tuples (row, col) corresponding to the zeros of the matrix.
    model: LpProblem to add constraints to.
    lpRows: dict of row variables and their degrees.
    lpCols: dict of column variables and their degrees.
    delta: float, error tolerance for density constraints.
    Attention : in the text it is written as ALPHA TO DO: (ALPHA : = 1-DELTA) !!!!
    """
    mu = 0.0001
    Big_R = len(rows_data) + 1
    Big_C = len(cols_data) + 1
    Big_M = Big_R+Big_C

    for row, _ in rows_data:
        row_edges = [v for u, v in edges if u == row]
        model += (
            lpSum(lpCols[col][0] for col in row_edges) - (1 - delta) * lpSum(lpCols[col][0] for col, _ in cols_data) >= 
            (lpRows[row][0]-1) * Big_R
        ), f"row_err_rate_1_{row}"
