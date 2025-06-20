import gurobipy as gp
from gurobipy import GRB
import itertools

def max_e_r(rows_data, cols_data, edges, delta, debug=0):
    """
    Create an ILP model for maximum edge-row optimization.
    
    Arguments:
    ----------
    rows_data: list of tuples (row, degree) of rows in the matrix.
    cols_data: list of tuples (col, degree) of columns in the matrix.
    edges: list of tuples (row, col) corresponding to the ones of the matrix.
    delta: float, error tolerance for density constraints.
    debug: int, debug level for output verbosity.

    Returns:
    --------
    gp.Model: The ILP model.
    
    This is an improved model from Chang et al.
    """
    try:
        model = gp.Model("max_e_r")
        model.setAttr("ModelSense", GRB.MAXIMIZE)

        # Create variables
        lp_rows, lp_cols, lp_cells = _create_variables(model, rows_data, cols_data, edges)

        # Set objective
        model.setObjective(gp.quicksum(lp_cells.values()), GRB.MAXIMIZE)

        # Add constraints
        _add_threshold_constraints(model, lp_rows, lp_cols, row_threshold=1, col_threshold=1)
        _add_matrix_structure_constraints(model, lp_rows, lp_cols, lp_cells, edges)
        _add_density_constraints(model, rows_data, cols_data, edges, lp_rows, lp_cols, delta)

        if debug >= 3:
            _print_model_info(model, rows_data, cols_data, edges, delta, 1, 1)

        return model
    except Exception as e:
        print(f"Error creating max_e_r model: {e}")
        return None


def max_e_wr(rows_data, cols_data, edges, rows_res, cols_res, prev_obj, delta, debug=0):
    """
    Create an ILP model for maximum edge-row optimization with warm restart.
    
    Arguments:
    ----------
    rows_data: list of tuples (row, degree) of rows in the matrix.
    cols_data: list of tuples (col, degree) of columns in the matrix.
    edges: list of tuples (row, col) corresponding to the ones of the matrix.
    rows_res: list of indices of rows (input solution).
    cols_res: list of indices of cols (input solution).
    prev_obj: previous objective value for improvement constraint.
    delta: float, error tolerance for density constraints.
    debug: int, debug level for output verbosity.

    Returns:
    --------
    gp.Model: The ILP model.
    """
    try:
        model = gp.Model("max_e_wr")
        model.setAttr("ModelSense", GRB.MAXIMIZE)

        # Create variables
        lp_rows, lp_cols, lp_cells = _create_variables(model, rows_data, cols_data, edges)

        # Set warm start values
        _set_warm_start_values(lp_rows, lp_cols, lp_cells, rows_res, cols_res, debug)

        # Set objective
        model.setObjective(gp.quicksum(lp_cells.values()), GRB.MAXIMIZE)

        # Add constraints
        row_threshold, col_threshold = 2, 2
        _add_threshold_constraints(model, lp_rows, lp_cols, row_threshold, col_threshold)
        _add_improvement_constraint(model, lp_cells, prev_obj)
        _add_matrix_structure_constraints(model, lp_rows, lp_cols, lp_cells, edges)
        _add_density_constraints(model, rows_data, cols_data, edges, lp_rows, lp_cols, delta)

        if debug >= 2:
            print(f"Got a lower bound: {prev_obj}")
        
        if debug >= 3:
            _print_model_info(model, rows_data, cols_data, edges, delta, row_threshold, col_threshold)

        return model
    except Exception as e:
        print(f"Error creating max_e_wr model: {e}")
        return None


def _create_variables(model, rows_data, cols_data, edges):
    """Create Gurobi variables for rows, columns, and cells."""
    lp_rows = {
        row: (model.addVar(vtype=GRB.BINARY, name=f'row_{row}'), degree)
        for row, degree in rows_data
    }
    
    lp_cols = {
        col: (model.addVar(vtype=GRB.BINARY, name=f'col_{col}'), degree)
        for col, degree in cols_data
    }
    
    lp_cells = {
        (row, col): model.addVar(vtype=GRB.CONTINUOUS, lb=0, ub=1, name=f'cell_{row}_{col}')
        for (row, col) in edges
    }
    
    return lp_rows, lp_cols, lp_cells


def _set_warm_start_values(lp_rows, lp_cols, lp_cells, rows_res, cols_res, debug):
    """Set initial values for warm start."""
    rows_res_set = set(map(int, rows_res))
    cols_res_set = set(map(int, cols_res))
    
    # Set row initial values
    row_initial_values = {}
    for row, (var, _) in lp_rows.items():
        val = 1 if row in rows_res_set else 0
        var.setAttr("Start", val)
        row_initial_values[row] = val

    # Set column initial values
    col_initial_values = {}
    for col, (var, _) in lp_cols.items():
        val = 1 if col in cols_res_set else 0
        var.setAttr("Start", val)
        col_initial_values[col] = val

    # Set cell initial values
    cell_initial_values = {}
    for (row, col), var in lp_cells.items():
        val = 1 if row in rows_res_set and col in cols_res_set else 0
        var.setAttr("Start", val)
        cell_initial_values[(row, col)] = val

    if debug >= 3:
        print("\nInitial point before solving:")
        print(f"ROWS: {row_initial_values}")
        print(f"COLS: {col_initial_values}")
        
    if debug >= 4:
        print(f"CELLS: {cell_initial_values}")


def _add_threshold_constraints(model, lp_rows, lp_cols, row_threshold, col_threshold):
    """Add row and column threshold constraints."""
    model.addConstr(gp.quicksum(lpvar for lpvar, _ in lp_rows.values()) >= row_threshold, "row_threshold")
    model.addConstr(gp.quicksum(lpvar for lpvar, _ in lp_cols.values()) >= col_threshold, "col_threshold")


def _add_improvement_constraint(model, lp_cells, prev_obj):
    """Add constraint for objective improvement."""
    model.addConstr(gp.quicksum(lp_cells.values()) >= prev_obj + 1, "improvement")


def _add_matrix_structure_constraints(model, lp_rows, lp_cols, lp_cells, edges):
    """Add constraints for matrix structure."""
    try:
        # Check for missing columns
        edge_cols = {col for _, col in edges}
        available_cols = set(lp_cols.keys())
        missing_cols = edge_cols - available_cols
        
        if missing_cols:
            print(f"Warning: Missing columns in cols_data: {missing_cols}")

        for row, col in edges:
            if row in lp_rows and col in lp_cols and (row, col) in lp_cells:
                model.addConstr(lp_rows[row][0] >= lp_cells[(row, col)], f'cell_{row}_{col}_row')
                model.addConstr(lp_cols[col][0] >= lp_cells[(row, col)], f'cell_{row}_{col}_col')
            else:
                print(f"Warning: row={row}, col={col} not found in variables!")
    except Exception as e:
        print(f"Error adding matrix structure constraints: {e}")


def _add_density_constraints(model, rows_data, cols_data, edges, lp_rows, lp_cols, delta):
    """Add row and column density constraints."""
    _add_row_density_constraints(rows_data, cols_data, edges, model, lp_rows, lp_cols, delta)
    _add_col_density_constraints(rows_data, cols_data, edges, model, lp_rows, lp_cols, delta)


def _add_row_density_constraints(rows_data, cols_data, edges, model, lp_rows, lp_cols, delta):
    """Add row density constraints to the model."""
    big_m = len(rows_data) + len(cols_data)
    
    for row, _ in rows_data:
        row_edges = [v for u, v in edges if u == row]
        model.addConstr(
            gp.quicksum(lp_cols[col][0] for col in row_edges) - 
            (1 - delta) * gp.quicksum(lp_cols[col][0] for col, _ in cols_data) >= 
            (lp_rows[row][0] - 1) * big_m,
            f"row_density_{row}"
        )


def _add_col_density_constraints(rows_data, cols_data, edges, model, lp_rows, lp_cols, delta):
    """Add column density constraints to the model."""
    big_m = len(rows_data) + len(cols_data)
    
    for col, _ in cols_data:
        col_edges = [u for u, v in edges if v == col]
        model.addConstr(
            gp.quicksum(lp_rows[row][0] for row in col_edges) - 
            (1 - delta) * gp.quicksum(lp_rows[row][0] for row, _ in rows_data) >= 
            (lp_cols[col][0] - 1) * big_m,
            f"col_density_{col}"
        )


def _print_model_info(model, rows_data, cols_data, edges, delta, row_threshold, col_threshold):
    """Print model information for debugging."""
    print()
    print('-' * 40)
    print(f"Solving model: {model.getAttr('ModelName')} with delta = {delta}")
    print(f"# rows_data = {len(rows_data)}, # cols_data = {len(cols_data)}, # edges = {len(edges)}")
    print(f"row_threshold = {row_threshold}")
    print(f"col_threshold = {col_threshold}")
    print('-' * 40)
