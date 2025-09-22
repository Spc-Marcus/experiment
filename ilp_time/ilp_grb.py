import os
import numpy as np
from typing import List, Tuple
from model.max_one_grb import max_Ones_gurobi
from model.max_one_grb_v2 import MaxOneModel
import contextlib
import sys
import gurobipy as grb

@contextlib.contextmanager
def suppress_gurobi_output():
    old_stdout = sys.stdout
    old_stderr = sys.stderr
    try:
        sys.stdout = open(os.devnull, 'w')
        sys.stderr = open(os.devnull, 'w')
        yield
    finally:
        sys.stdout.close()
        sys.stderr.close()
        sys.stdout = old_stdout
        sys.stderr = old_stderr

def find_quasi_dens_matrix_max_ones(
    input_matrix: np.ndarray,
    error_rate: float = 0.025
) -> Tuple[List[int], List[int], bool]:
    """
    Find a quasi-biclique in a binary matrix using integer linear programming optimization (Gurobi).
    """
    X_problem = input_matrix.copy()
    cols_sorted = np.argsort(X_problem.sum(axis=0))[::-1]
    rows_sorted = np.argsort(X_problem.sum(axis=1))[::-1]
    m = len(rows_sorted)
    n = len(cols_sorted)
    if m == 0 or n == 0:
        return [], [], False
    
    seed_cols = max(n // 3, 2)
    if n > 50:
        step_n = 10
    else:
        step_n = 2
    
    for x in range(m // 3, m, 10):
        for y in range(seed_cols, n, step_n):
            nb_of_ones = 0
            for row in rows_sorted[:x]:
                for col in cols_sorted[:y]:
                    nb_of_ones += X_problem[row, col]
            ratio_ones = nb_of_ones / (x * y) if (x * y) > 0 else 0
            if ratio_ones > 0.99:
                seed_cols = y
    
    try:
        with suppress_gurobi_output():
            # --- PHASE 1: SEED ---
            seed_row_indices = rows_sorted
            seed_col_indices = cols_sorted[:seed_cols]
            row_degrees = np.sum(X_problem[seed_row_indices, :][:, seed_col_indices] == 1, axis=1)
            col_degrees = np.sum(X_problem[seed_row_indices, :][:, seed_col_indices] == 1, axis=0)
            rows_data = [(int(r), int(row_degrees[i])) for i, r in enumerate(seed_row_indices)]
            cols_data = [(int(c), int(col_degrees[i])) for i, c in enumerate(seed_col_indices)]
            edges = []
            for i, r in enumerate(seed_row_indices):
                for j, c in enumerate(seed_col_indices):
                    if X_problem[r, c] == 1:
                        edges.append((int(r), int(c)))
            
            model = max_Ones_gurobi(rows_data, cols_data, edges, 0)
            model.setParam('OutputFlag', 0)
            model.setParam('MIPGap', 0.05)
            model.setParam('TimeLimit', 20)
            model.optimize()
            
            status = model.Status
            if status in (grb.GRB.INF_OR_UNBD, grb.GRB.INFEASIBLE, grb.GRB.UNBOUNDED):
                return [], [], False
            elif status == grb.GRB.TIME_LIMIT or status == grb.GRB.OPTIMAL:
                rw = []
                cl = []
                for v in model.getVars():
                    if v.VarName.startswith('row_') and v.X > 0.5:
                        rw.append(int(v.VarName.split('_')[1]))
                    elif v.VarName.startswith('col_') and v.X > 0.5:
                        cl.append(int(v.VarName.split('_')[1]))
            else:
                return [], [], False
            
            # --- PHASE 2: EXTENSION COLONNES ---
            rem_cols = [c for c in cols_sorted if c not in cl]
            if len(rw) > 0:
                rem_cols_sum = X_problem[rw][:, rem_cols].sum(axis=0)
                potential_cols = [c for idx, c in enumerate(rem_cols) if rem_cols_sum[idx] > 0.9 * len(rw)]
            else:
                potential_cols = []
            
            if potential_cols:
                all_col_indices = cl + potential_cols
                row_degrees = np.sum(X_problem[rw, :][:, all_col_indices] == 1, axis=1)
                rows_data = [(int(r), int(row_degrees[i])) for i, r in enumerate(rw)]
                col_degrees = np.sum(X_problem[rw, :][:, all_col_indices] == 1, axis=0)
                cols_data = [(int(c), int(col_degrees[i])) for i, c in enumerate(all_col_indices)]
                edges = []
                for i, r in enumerate(rw):
                    for j, c in enumerate(all_col_indices):
                        if X_problem[r, c] == 1:
                            edges.append((int(r), int(c)))
                
                model = max_Ones_gurobi(rows_data, cols_data, edges, error_rate)
                model.setParam('OutputFlag', 0)
                model.setParam('MIPGap', 0.05)
                model.setParam('TimeLimit', 180)
                model.optimize()
                
                status = model.Status
                if status in (grb.GRB.INF_OR_UNBD, grb.GRB.INFEASIBLE, grb.GRB.UNBOUNDED):
                    return [], [], False
                elif status == grb.GRB.TIME_LIMIT or status == grb.GRB.OPTIMAL:
                    rw = []
                    cl = []
                    for v in model.getVars():
                        if v.VarName.startswith('row_') and v.X > 0.5:
                            rw.append(int(v.VarName.split('_')[1]))
                        elif v.VarName.startswith('col_') and v.X > 0.5:
                            cl.append(int(v.VarName.split('_')[1]))
                else:
                    return [], [], False
            
            # --- PHASE 3: EXTENSION LIGNES ---
            rem_rows = [r for r in rows_sorted if r not in rw]
            if len(cl) > 0:
                rem_rows_sum = X_problem[rem_rows][:, cl].sum(axis=1)
                potential_rows = [r for idx, r in enumerate(rem_rows) if rem_rows_sum[idx] > 0.5 * len(cl)]
            else:
                potential_rows = []
            
            if potential_rows:
                all_row_indices = rw + potential_rows
                row_degrees = np.sum(X_problem[all_row_indices, :][:, cl] == 1, axis=1)
                rows_data = [(int(r), int(row_degrees[i])) for i, r in enumerate(all_row_indices)]
                col_degrees = np.sum(X_problem[all_row_indices, :][:, cl] == 1, axis=0)
                cols_data = [(int(c), int(col_degrees[i])) for i, c in enumerate(cl)]
                edges = []
                for i, r in enumerate(all_row_indices):
                    for j, c in enumerate(cl):
                        if X_problem[r, c] == 1:
                            edges.append((int(r), int(c)))
                
                model = max_Ones_gurobi(rows_data, cols_data, edges, error_rate)
                model.setParam('OutputFlag', 0)
                model.setParam('MIPGap', 0.05)
                model.setParam('TimeLimit', 20)
                model.optimize()
                
                status = model.Status
                if status in (grb.GRB.INF_OR_UNBD, grb.GRB.INFEASIBLE, grb.GRB.UNBOUNDED):
                    return [], [], False
                elif status == grb.GRB.TIME_LIMIT or status == grb.GRB.OPTIMAL:
                    rw = []
                    cl = []
                    for v in model.getVars():
                        if v.VarName.startswith('row_') and v.X > 0.5:
                            rw.append(int(v.VarName.split('_')[1]))
                        elif v.VarName.startswith('col_') and v.X > 0.5:
                            cl.append(int(v.VarName.split('_')[1]))
                else:
                    return [], [], False
            
            status = model.Status
            if status in (grb.GRB.INF_OR_UNBD, grb.GRB.INFEASIBLE, grb.GRB.UNBOUNDED):
                return [], [], False
            elif status == grb.GRB.TIME_LIMIT:
                return rw, cl, True
            elif status != grb.GRB.OPTIMAL:
                return [], [], False
            return rw, cl, True
    except Exception:
        return [], [], False

def find_quasi_biclique_max_one_V2(
    input_matrix: np.ndarray,
    error_rate: float = 0.025,
) -> Tuple[List[int], List[int], bool]:
    """
    Find a quasi-biclique using MaxOneModel.
    """
    X_problem = input_matrix.copy()
    cols_sorted = np.argsort(X_problem.sum(axis=0))[::-1]
    rows_sorted = np.argsort(X_problem.sum(axis=1))[::-1]
    m = len(rows_sorted)
    n = len(cols_sorted)
    if m == 0 or n == 0:
        return [], [], False
    
    seed_cols = max(n // 3, 2)
    if n > 50:
        step_n = 10
    else:
        step_n = 2
    
    for x in range(m // 3, m, 10):
        for y in range(seed_cols, n, step_n):
            nb_of_ones = 0
            for row in rows_sorted[:x]:
                for col in cols_sorted[:y]:
                    nb_of_ones += X_problem[row, col]
            ratio_ones = nb_of_ones / (x * y) if (x * y) > 0 else 0
            if ratio_ones > 0.99:
                seed_cols = y
    
    try:
        with suppress_gurobi_output():
            # --- PHASE 1: SEED ---
            row_degrees = np.sum(X_problem == 1, axis=1)
            rows_data = [(r, int(row_degrees[r])) for r in range(m)]
            col_degrees = np.sum(X_problem == 1, axis=0)
            cols_data = [(int(c), int(col_degrees[c])) for c in range(n)]
            edges = []
            cols_sums = X_problem.sum(axis=0)
            cols_sorted = np.argsort(cols_sums)[::-1]
            seed_cols = min(max(n // 3, min(n, 5)), 50)
            no_use_cols_seed = cols_sorted[seed_cols:]

            for r in range(m):
                for c in range(n):
                    if X_problem[r, c] == 1:
                        edges.append((int(r), int(c)))
            
            model = MaxOneModel(rows_data, cols_data, edges, error_rate)
            model.model.setParam('OutputFlag', 0)
            model.build(error_rate)
            model.model.setParam('MIPGap', 0.05)
            model.model.setParam('TimeLimit', 20)
            model.add_forced_cols_zero(no_use_cols_seed)
            model.add_density_constraint(0.0)
            model.optimize()
           
            if model.status == 2:
                rw = model.get_selected_rows()
                cl = model.get_selected_cols()
            else:
                return [], [], False
            
            # --- PHASE 2: EXTENSION COLONNES ---
            no_use_rows_seed = [r for r in range(m) if r not in rw]
            potential_cols = [c for idx,c in enumerate(no_use_cols_seed) if np.sum(X_problem[rw, :][:, c]) > 0.9 * len(rw)]
            no_use_cols_seed = [c for c in range(n) if c not in cl and c not in potential_cols]
            model.remove_forced_cols_zero(no_use_cols_seed)
            model.add_forced_cols_zero(no_use_cols_seed)
            model.add_forced_rows_zero(no_use_rows_seed)
            model.add_improvement_constraint(model.objVal)
            model.update_density_constraints(error_rate)
            model.optimize()

            if model.status == 2:
                rw = model.get_selected_rows()
                cl = model.get_selected_cols()
            else: 
                return rw, cl, True

            # --- PHASE 3: EXTENSION LIGNES ---

            no_use_cols_seed = [c for c in range(n) if c not in cl]
            potential_rows = [r for idx,r in enumerate(no_use_rows_seed) if np.sum(X_problem[:, cl][:, cl][r, :]) > 0.5 * len(cl)]
            no_use_rows_seed = [r for r in range(m) if r not in rw and r not in potential_rows]
            model.remove_forced_rows_zero(no_use_rows_seed)
            model.add_forced_rows_zero(no_use_rows_seed)
            model.add_forced_cols_zero(no_use_cols_seed)
            model.add_improvement_constraint(model.objVal)
            model.update_density_constraints(error_rate)
            model.optimize()   

            if model.status == 2:
                rw = model.get_selected_rows()
                cl = model.get_selected_cols()
            else: 
                return rw, cl, True

            if rw and cl:
                return rw, cl, True
            else:
                return [], [], False
            
    except Exception:
        return [], [], False
    