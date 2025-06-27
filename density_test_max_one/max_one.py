import pulp as plp

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
    model = plp.LpProblem(name='maximize_ones', sense=plp.LpMaximize)

    # ------------------------------------------------------------------------ #
    # Variables
    # ------------------------------------------------------------------------ #
    lpRows = {}
    for row, degree in rows_data:
        var = plp.LpVariable(name=f'row_{row}', cat='Binary', lowBound=0, upBound=1)
        lpRows[row] = (var, degree)
    
    lpCols = {}
    for col, degree in cols_data:
        var = plp.LpVariable(name=f'col_{col}', cat='Binary', lowBound=0, upBound=1)
        lpCols[col] = (var, degree)
    
    lpCells = {}
    for row, _ in rows_data:
        for col, _ in cols_data:
            var = plp.LpVariable(name=f'cell_{row}_{col}', cat='Binary', lowBound=0, upBound=1)
            if (row, col) in edges:
                lpCells[(row, col)] = (var, 1)
            else:
                lpCells[(row, col)] = (var, 0)

    # ------------------------------------------------------------------------ #
    # Objective function
    # ------------------------------------------------------------------------ #
    obj_expr = plp.lpSum([cellValue * lpvar for lpvar, cellValue in lpCells.values()])
    model += obj_expr

    # ------------------------------------------------------------------------ #
    # Constraints
    # ------------------------------------------------------------------------ #
    for row, col in lpCells:
        model += lpRows[row][0] >= lpCells[(row, col)][0], f'cell_{row}_{col}_1'
        model += lpCols[col][0] >= lpCells[(row, col)][0], f'cell_{row}_{col}_2'
        model += lpRows[row][0] + lpCols[col][0] - 1 <= lpCells[(row, col)][0], f'cell_{row}_{col}_3'

    err_expr = plp.lpSum([(1-cellValue) * lpvar for lpvar, cellValue in lpCells.values()])
    total_expr = plp.lpSum([lpvar for lpvar, _ in lpCells.values()])
    model += err_expr <= rho * total_expr, 'err_rate'

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
    model = plp.LpProblem(name='maximize_ones_compacted', sense=plp.LpMaximize)

    # ------------------------------------------------------------------------ #
    # Variables
    # ------------------------------------------------------------------------ #
    lpRows = {}
    for row, degree in rows_data:
        var = plp.LpVariable(name=f'row_{row}', cat='Binary', lowBound=0, upBound=1)
        lpRows[row] = (var, degree)
    
    lpCols = {}
    for col, degree in cols_data:
        var = plp.LpVariable(name=f'col_{col}', cat='Binary', lowBound=0, upBound=1)
        lpCols[col] = (var, degree)
    
    lpCells = {}
    for row, _ in rows_data:
        for col, _ in cols_data:
            var = plp.LpVariable(name=f'cell_{row}_{col}', cat='Binary', lowBound=0, upBound=1)
            if (row, col) in edges:
                lpCells[(row, col)] = (var, 1)
            else:
                lpCells[(row, col)] = (var, 0)

    # ------------------------------------------------------------------------ #
    # Objective function
    # ------------------------------------------------------------------------ #
    obj_expr = plp.lpSum([cellValue * lpvar for lpvar, cellValue in lpCells.values()])
    model += obj_expr

    # ------------------------------------------------------------------------ #
    # Constraints
    # ------------------------------------------------------------------------ #
    row_threshold = 1
    col_threshold = 1
    
    row_sum = plp.lpSum([lpvar for lpvar, _ in lpRows.values()])
    model += row_sum >= row_threshold, "row_threshold"
    
    col_sum = plp.lpSum([lpvar for lpvar, _ in lpCols.values()])
    model += col_sum >= col_threshold, "col_threshold"
    
    for row, col in lpCells:
        model += lpRows[row][0] + lpCols[col][0] - 1 <= lpCells[(row, col)][0], f'cell_{row}_{col}_3'
    
    # compacting with degree 
    #########################################
    for col, _ in cols_data:
        col_edges = [u for u, v in edges if v == col] 
        cell_sum = plp.lpSum([lpCells[(row, col)][0] for row in col_edges])
        model += cell_sum <= lpCols[col][1] * lpCols[col][0], f"col_degre_{col}"
    
    for row, _ in rows_data:
        row_edges = [v for u, v in edges if u == row] 
        cell_sum = plp.lpSum([lpCells[(row, col)][0] for col in row_edges])
        model += cell_sum <= lpRows[row][1] * lpRows[row][0], f"row_degre_{row}"
    #########################################

    err_expr = plp.lpSum([(1-cellValue) * lpvar for lpvar, cellValue in lpCells.values()])
    total_expr = plp.lpSum([lpvar for lpvar, _ in lpCells.values()])
    model += err_expr <= rho * total_expr, 'err_rate'

    return model
    model.addConstr(err_expr <= rho * total_expr, name='err_rate')

    return model
