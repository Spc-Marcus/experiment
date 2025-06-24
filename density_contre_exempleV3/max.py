import gurobipy as gp
from gurobipy import GRB
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
    gp.Model:
    The ILP model.
    
    This is an improved model from Chang et al. 
    """
    model = gp.Model("max_e_r")
    model.setAttr("ModelSense", GRB.MAXIMIZE)

    # Variables for rows and columns
    lpRows = {row: (model.addVar(vtype=GRB.BINARY, name=f'row_{row}'), degree) 
              for row, degree in rows_data}
    lpCols = {col: (model.addVar(vtype=GRB.BINARY, name=f'col_{col}'), degree) 
              for col, degree in cols_data}
    lpCells = {(row, col): model.addVar(vtype=GRB.CONTINUOUS, lb=0, ub=1, 
                                        name=f'cell_{row}_{col}') 
               for (row, col) in edges}

    # Objective
    model.setObjective(gp.quicksum(lpCells.values()), GRB.MAXIMIZE)

    # Constraints for row and column thresholds
    row_threshold = 1
    col_threshold = 1
    
    
    
    model.addConstr(gp.quicksum(lpvar for lpvar, _ in lpRows.values()) >= row_threshold, 
                    "row_threshold")
    model.addConstr(gp.quicksum(lpvar for lpvar, _ in lpCols.values()) >= col_threshold, 
                    "col_threshold")

    # Constraints for matrix structure
    # Check for missing columns in cols_data
    missing_cols = {col for _, col in edges} - {col for col, _ in cols_data}
    if missing_cols:
        print(f"Warning: These columns are missing from cols_data but appear in edges: {missing_cols}")

    for row, col in edges:
        if row in lpRows and col in lpCols and (row, col) in lpCells:
            model.addConstr(lpRows[row][0] >= lpCells[(row, col)], f'cell_{row}_{col}_1')
            model.addConstr(lpCols[col][0] >= lpCells[(row, col)], f'cell_{row}_{col}_2')
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
    rows_res : list of indices of rows (input solution)
    cols_res : list of indices of cols (input solution)
    delta: float, error tolerance for density constraints.

    Returns:
    --------
    gp.Model:
    The ILP model.
    """
    model = gp.Model("max_e_wr")
    model.setAttr("ModelSense", GRB.MAXIMIZE)

    # Variables for rows and columns
    lpRows = {row: (model.addVar(vtype=GRB.BINARY, name=f'row_{row}'), degree) 
              for row, degree in rows_data}
    lpCols = {col: (model.addVar(vtype=GRB.BINARY, name=f'col_{col}'), degree) 
              for col, degree in cols_data}
    lpCells = {(row, col): model.addVar(vtype=GRB.CONTINUOUS, lb=0, ub=1, 
                                        name=f'cell_{row}_{col}') 
               for (row, col) in edges}

    cols_res_set = set(map(int, cols_res)) 
    rows_res_set = set(map(int, rows_res))

    
    # Warm start (Gurobi way: use Start attribute)
    cells_indices = []
    row_initial_values = {}
    for row, (var, _) in lpRows.items():
        if row in rows_res_set:
            val = 1 
            var.Start = val
            row_initial_values[row] = val
        else:
            val = 0
            row_initial_values[row] = val

    col_initial_values = {}
    for col, (var, _) in lpCols.items():
        if col in cols_res_set:
            val = 1 
            var.Start = val
            col_initial_values[col] = val
        else: 
            val = 0
            col_initial_values[col] = val

    cell_initial_values = {}
    for (row, col), var in lpCells.items():
        if (row, col) in itertools.product(row_initial_values, col_initial_values):
            var.Start = 1
            val = 1 
            cell_initial_values[(row, col)] = val
        else:
            var.Start = 0
            val = 0
            cell_initial_values[(row, col)] = val
    
   

    # Objective 
    model.setObjective(gp.quicksum(lpCells.values()), GRB.MAXIMIZE)

    row_threshold = 2
    col_threshold = 2

    

    # Add constraints for row and column thresholds
    model.addConstr(gp.quicksum(lpvar for lpvar, _ in lpRows.values()) >= row_threshold, 
                    "row_threshold")
    model.addConstr(gp.quicksum(lpvar for lpvar, _ in lpCols.values()) >= col_threshold, 
                    "col_threshold")

    # Constraint for objective improvement
    model.addConstr(gp.quicksum(lpCells.values()) >= prev_obj + 1, "improvement")

    # Constraints for matrix structure
    for row, col in edges:
        if (row, col) in lpCells:
            model.addConstr(lpRows[row][0] >= lpCells[(row, col)], f'cell_{row}_{col}_1')
            model.addConstr(lpCols[col][0] >= lpCells[(row, col)], f'cell_{row}_{col}_2')
        else:
            print(f"Warning: ({row}, {col}) not found in lpCells!")

    __row_density(rows_data, cols_data, edges, model, lpRows, lpCols, delta)
    __col_density(rows_data, cols_data, edges, model, lpRows, lpCols, delta)

    return model 

def __col_density(rows_data, cols_data, edges, model, lpRows, lpCols, delta):
    """
    Adds col density constraints to the model.

    Arguments:
    ----------
    rows_data: list of tuples (row, degree) of rows in the matrix.
    cols_data: list of tuples (col, degree) of columns in the matrix.
    edges: list of tuples (row, col) corresponding to the one of the matrix.
    model: gp.Model to add constraints to.
    lpRows: dict of row variables and their degrees.
    lpCols: dict of column variables and their degrees.
    delta: float, error tolerance for density constraints. 
    Attention : in the text it is written as ALPHA TO DO: (ALPHA : = 1-DELTA) !!!!
    """
    mu = 0.0001
    Big_R = len(rows_data)
    Big_C = len(cols_data)
    Big_M = Big_R + Big_C 
    
    for col, _ in cols_data:
        col_edges = [u for u, v in edges if v == col]
        model.addConstr(
            gp.quicksum(lpRows[row][0] for row in col_edges) - 
            (1 - delta) * gp.quicksum(lpRows[row][0] for row, _ in rows_data) >= 
            (lpCols[col][0] - 1) * Big_C,
            f"col_err_rate_1_{col}"
        )

def __row_density(rows_data, cols_data, edges, model, lpRows, lpCols, delta):
    """
    Adds row density constraints to the model.

    Arguments:
    ----------
    rows_data: list of tuples (row, degree) of rows in the matrix.
    cols_data: list of tuples (col, degree) of columns in the matrix.
    edges: list of tuples (row, col) corresponding to the one of the matrix.
    model: gp.Model to add constraints to.
    lpRows: dict of row variables and their degrees.
    lpCols: dict of column variables and their degrees.
    delta: float, error tolerance for density constraints.
    Attention : in the text it is written as ALPHA TO DO: (ALPHA : = 1-DELTA) !!!!
    """
    mu = 0.0001
    Big_R = len(rows_data)
    Big_C = len(cols_data)
    Big_M = Big_R + Big_C

    for row, _ in rows_data:
        row_edges = [v for u, v in edges if u == row]
        model.addConstr(
            gp.quicksum(lpCols[col][0] for col in row_edges) - 
            (1 - delta) * gp.quicksum(lpCols[col][0] for col, _ in cols_data) >= 
            (lpRows[row][0] - 1) * Big_R,
            f"row_err_rate_1_{row}"
        )
