from pulp import LpProblem, LpVariable, LpMaximize, lpSum
import gurobipy as gp
from gurobipy import GRB

def max_Ones(rows_data, cols_data, edges, rho):
    """
    ARGUMENTS:
    ----------
    * rows_data: list of the tuples (row, degree) of rows the matrix.
    * cols_data: list of the tuples (col, degree) of columns the matrix.
    * edges: list of tuple (row,col) corresponding to the zeros of the matrix.
    * rho: percentage of accepted zeros 
    """

    # ------------------------------------------------------------------------ #
    # Model with maximization
    # ------------------------------------------------------------------------ #
    model = LpProblem(name='maximize_ones', sense=LpMaximize)

    # ------------------------------------------------------------------------ #
    # Variables
    # ------------------------------------------------------------------------ #
    lpRows = {row: (LpVariable(f'row_{row}', cat='Integer',
                              lowBound=0, upBound=1), degree) for row, degree in rows_data}
    lpCols = {col: (LpVariable(f'col_{col}', cat='Integer',
                              lowBound=0, upBound=1), degree) for col, degree in cols_data}
    lpCells = {}
    for row, _ in rows_data:
        for col, _ in cols_data:
            if (row, col) in edges:
                lpCells[(row, col)] = (LpVariable(
                    f'cell_{row}_{col}', cat='Integer', lowBound=0, upBound=1), 1)
            else:
                lpCells[(row, col)] = (LpVariable( 
                    f'cell_{row}_{col}', cat='Integer', lowBound=0, upBound=1), 0)

    # ------------------------------------------------------------------------ #
    # Objective function
    # ------------------------------------------------------------------------ #
    model += lpSum([cellValue*lpvar for lpvar,
                   cellValue in lpCells.values()]), 'maximize_ones_in_matrix'

    # ------------------------------------------------------------------------ #
    # Constraints
    # ------------------------------------------------------------------------ #

    for row, col in lpCells:
        model += (lpRows[row][0] >= lpCells[(row, col)][0]), f'cell_{row}_{col}_1'
        model += (lpCols[col][0] >= lpCells[(row, col)][0]), f'cell_{row}_{col}_2'
        model += (lpRows[row][0]+lpCols[col][0] -1 <= lpCells[(row, col)][0]), f'cell_{row}_{col}_3'

    model += (lpSum([(1-cellValue)*lpvar for lpvar, cellValue in lpCells.values()]) <= rho *
              lpSum([lpvar for lpvar, _ in lpCells.values()])), f'err_rate'

    return model


def max_Ones_comp(rows_data, cols_data, edges, rho):
    """
    ARGUMENTS:
    ----------
    * rows_data: list of the tuples (row, degree) of rows the matrix.
    * cols_data: list of the tuples (col, degree) of columns the matrix.
    * edges: list of tuple (row,col) corresponding to the zeros of the matrix.
    * rho: percentage of accepted zeros 
    """

    # ------------------------------------------------------------------------ #
    # Model with maximization
    # ------------------------------------------------------------------------ #
    model = LpProblem(name='maximize_ones_compacted', sense=LpMaximize)

    # ------------------------------------------------------------------------ #
    # Variables
    # ------------------------------------------------------------------------ #
    lpRows = {row: (LpVariable(f'row_{row}', cat='Integer',
                              lowBound=0, upBound=1), degree) for row, degree in rows_data}
    lpCols = {col: (LpVariable(f'col_{col}', cat='Integer',
                              lowBound=0, upBound=1), degree) for col, degree in cols_data}
    lpCells = {}
    for row, _ in rows_data:
        for col, _ in cols_data:
            if (row, col) in edges:
                lpCells[(row, col)] = (LpVariable(
                    f'cell_{row}_{col}', cat='Integer', lowBound=0, upBound=1), 1)
            else:
                lpCells[(row, col)] = (LpVariable( 
                    f'cell_{row}_{col}', cat='Integer', lowBound=0, upBound=1), 0)

    # ------------------------------------------------------------------------ #
    # Objective function
    # ------------------------------------------------------------------------ #
    model += lpSum([cellValue*lpvar for lpvar,
                   cellValue in lpCells.values()]), 'maximize_ones_compacted'

    # ------------------------------------------------------------------------ #
    # Constraints
    # ------------------------------------------------------------------------ #
    row_threshold = 2
    col_threshold = 2
    model += (lpSum(lpvar for lpvar, _ in lpRows.values()) >= row_threshold), "row_threshold"
    model += (lpSum(lpvar for lpvar, _ in lpCols.values()) >= col_threshold), "col_threshold"
    
    for row, col in lpCells:
        model += (lpRows[row][0]+lpCols[col][0] -1 <= lpCells[(row, col)][0]), f'cell_{row}_{col}_3'
    
    # compacting with degree 
    #########################################
    for col, _ in cols_data:
        col_edges = [u for u, v in edges if v == col] 
        model += (
            lpSum(lpCells[(row, col)][0] for row in col_edges) <= lpCols[col][1]*lpCols[col][0]
        ), f"col_degre_{col}"
    
    for row, _ in rows_data:
        row_edges = [v for u, v in edges if u == row] 
        model += (
            lpSum(lpCells[(row, col)][0] for col in row_edges) <= lpRows[row][1]*lpRows[row][0]
        ), f"row_degre_{row}"
    #########################################

    model += (lpSum([(1-cellValue)*lpvar for lpvar, cellValue in lpCells.values()]) <= rho *
              lpSum([lpvar for lpvar, _ in lpCells.values()])), f'err_rate'

    return model


def max_Ones_gurobi(rows_data, cols_data, edges, rho):
    """
    ARGUMENTS:
    ----------
    * rows_data: list of the tuples (row, degree) of rows the matrix.
    * cols_data: list of the tuples (col, degree) of columns the matrix.
    * edges: list of tuple (row,col) corresponding to the zeros of the matrix.
    * rho: percentage of accepted zeros 
    """

    # ------------------------------------------------------------------------ #
    # Model with maximization
    # ------------------------------------------------------------------------ #
    model = gp.Model(name='maximize_ones')
    model.setAttr(GRB.Attr.ModelSense, GRB.MAXIMIZE)

    # ------------------------------------------------------------------------ #
    # Variables
    # ------------------------------------------------------------------------ #
    lpRows = {}
    for row, degree in rows_data:
        var = model.addVar(name=f'row_{row}', vtype=GRB.BINARY, lb=0, ub=1)
        lpRows[row] = (var, degree)
    
    lpCols = {}
    for col, degree in cols_data:
        var = model.addVar(name=f'col_{col}', vtype=GRB.BINARY, lb=0, ub=1)
        lpCols[col] = (var, degree)
    
    lpCells = {}
    for row, _ in rows_data:
        for col, _ in cols_data:
            var = model.addVar(name=f'cell_{row}_{col}', vtype=GRB.BINARY, lb=0, ub=1)
            if (row, col) in edges:
                lpCells[(row, col)] = (var, 1)
            else:
                lpCells[(row, col)] = (var, 0)

    # ------------------------------------------------------------------------ #
    # Objective function
    # ------------------------------------------------------------------------ #
    obj_expr = gp.quicksum(cellValue * lpvar for lpvar, cellValue in lpCells.values())
    model.setObjective(obj_expr, GRB.MAXIMIZE)

    # ------------------------------------------------------------------------ #
    # Constraints
    # ------------------------------------------------------------------------ #
    for row, col in lpCells:
        model.addConstr(lpRows[row][0] >= lpCells[(row, col)][0], name=f'cell_{row}_{col}_1')
        model.addConstr(lpCols[col][0] >= lpCells[(row, col)][0], name=f'cell_{row}_{col}_2')
        model.addConstr(lpRows[row][0] + lpCols[col][0] - 1 <= lpCells[(row, col)][0], name=f'cell_{row}_{col}_3')

    err_expr = gp.quicksum((1-cellValue) * lpvar for lpvar, cellValue in lpCells.values())
    total_expr = gp.quicksum(lpvar for lpvar, _ in lpCells.values())
    model.addConstr(err_expr <= rho * total_expr, name='err_rate')

    return model


def max_Ones_comp_gurobi(rows_data, cols_data, edges, rho):
    """
    ARGUMENTS:
    ----------
    * rows_data: list of the tuples (row, degree) of rows the matrix.
    * cols_data: list of the tuples (col, degree) of columns the matrix.
    * edges: list of tuple (row,col) corresponding to the zeros of the matrix.
    * rho: percentage of accepted zeros 
    """

    # ------------------------------------------------------------------------ #
    # Model with maximization
    # ------------------------------------------------------------------------ #
    model = gp.Model(name='maximize_ones_compacted')
    model.setAttr(GRB.Attr.ModelSense, GRB.MAXIMIZE)

    # ------------------------------------------------------------------------ #
    # Variables
    # ------------------------------------------------------------------------ #
    lpRows = {}
    for row, degree in rows_data:
        var = model.addVar(name=f'row_{row}', vtype=GRB.BINARY, lb=0, ub=1)
        lpRows[row] = (var, degree)
    
    lpCols = {}
    for col, degree in cols_data:
        var = model.addVar(name=f'col_{col}', vtype=GRB.BINARY, lb=0, ub=1)
        lpCols[col] = (var, degree)
    
    lpCells = {}
    for row, _ in rows_data:
        for col, _ in cols_data:
            var = model.addVar(name=f'cell_{row}_{col}', vtype=GRB.BINARY, lb=0, ub=1)
            if (row, col) in edges:
                lpCells[(row, col)] = (var, 1)
            else:
                lpCells[(row, col)] = (var, 0)

    # ------------------------------------------------------------------------ #
    # Objective function
    # ------------------------------------------------------------------------ #
    obj_expr = gp.quicksum(cellValue * lpvar for lpvar, cellValue in lpCells.values())
    model.setObjective(obj_expr, GRB.MAXIMIZE)

    # ------------------------------------------------------------------------ #
    # Constraints
    # ------------------------------------------------------------------------ #
    row_threshold = 2
    col_threshold = 2
    
    row_sum = gp.quicksum(lpvar for lpvar, _ in lpRows.values())
    model.addConstr(row_sum >= row_threshold, name="row_threshold")
    
    col_sum = gp.quicksum(lpvar for lpvar, _ in lpCols.values())
    model.addConstr(col_sum >= col_threshold, name="col_threshold")
    
    for row, col in lpCells:
        model.addConstr(lpRows[row][0] + lpCols[col][0] - 1 <= lpCells[(row, col)][0], name=f'cell_{row}_{col}_3')
    
    # compacting with degree 
    #########################################
    for col, _ in cols_data:
        col_edges = [u for u, v in edges if v == col] 
        cell_sum = gp.quicksum(lpCells[(row, col)][0] for row in col_edges)
        model.addConstr(cell_sum <= lpCols[col][1] * lpCols[col][0], name=f"col_degre_{col}")
    
    for row, _ in rows_data:
        row_edges = [v for u, v in edges if u == row] 
        cell_sum = gp.quicksum(lpCells[(row, col)][0] for col in row_edges)
        model.addConstr(cell_sum <= lpRows[row][1] * lpRows[row][0], name=f"row_degre_{row}")
    #########################################

    err_expr = gp.quicksum((1-cellValue) * lpvar for lpvar, cellValue in lpCells.values())
    total_expr = gp.quicksum(lpvar for lpvar, _ in lpCells.values())
    model.addConstr(err_expr <= rho * total_expr, name='err_rate')

    return model