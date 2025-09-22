import gurobipy as gp
from gurobipy import GRB
import math


class MaxOneModel:
    def __init__(self, rows_data, cols_data, edges, error_rate):
        self.rows_data = rows_data
        self.cols_data = cols_data
        self.edges = set(edges)
        self.error_rate = float(error_rate)
        self.model = gp.Model("max_one_grb_v2")
        self.model.setAttr('ModelSense', GRB.MAXIMIZE)

        self.lp_rows = self._create_row_variables()
        self.lp_cols = self._create_col_variables()
        self.lp_cells = self._create_cell_variables()

        self._is_built = False
        self._forced_row_constrs = {}
        self._forced_col_constrs = {}
        self._improvement_constr = None
        self._density_constrs = []
        self._threshold_constrs = []

        self._err_expr = None
        self._tot_expr = None
        self._ones_expr = None

    def _create_row_variables(self):
        return {row: (self.model.addVar(vtype=GRB.BINARY, lb=0, ub=1, name=f'row_{row}'), degree)
                for row, degree in self.rows_data}

    def _create_col_variables(self):
        return {col: (self.model.addVar(vtype=GRB.BINARY, lb=0, ub=1, name=f'col_{col}'), degree)
                for col, degree in self.cols_data}

    def _create_cell_variables(self):
        lpCells = {}
        for row, _ in self.rows_data:
            for col, _ in self.cols_data:
                var = self.model.addVar(vtype=GRB.BINARY, lb=0, ub=1, name=f'cell_{row}_{col}')
                cell_value = 1 if (row, col) in self.edges else 0
                lpCells[(row, col)] = (var, cell_value)
        return lpCells

    def add_forced_rows_zero(self, rows):
        for row in rows:
            if row in self.lp_rows and row not in self._forced_row_constrs:
                self._forced_row_constrs[row] = self.model.addConstr(self.lp_rows[row][0] == 0,
                                                                     f"force_row_{row}_zero")

    def remove_forced_rows_zero(self, rows):
        for row in rows:
            if row in self._forced_row_constrs:
                self.model.remove(self._forced_row_constrs[row])
                del self._forced_row_constrs[row]

    def add_forced_cols_zero(self, cols):
        for col in cols:
            if col in self.lp_cols and col not in self._forced_col_constrs:
                self._forced_col_constrs[col] = self.model.addConstr(self.lp_cols[col][0] == 0,
                                                                     f"force_col_{col}_zero")

    def remove_forced_cols_zero(self, cols):
        for col in cols:
            if col in self._forced_col_constrs:
                self.model.remove(self._forced_col_constrs[col])
                del self._forced_col_constrs[col]

    def add_improvement_constraint(self, prev_obj):
        if self._improvement_constr is not None:
            self.model.remove(self._improvement_constr)
        self._improvement_constr = self.model.addConstr(self._ones_expr >= prev_obj, "improvement")

    def remove_improvement_constraint(self):
        if self._improvement_constr is not None:
            self.model.remove(self._improvement_constr)
            self._improvement_constr = None

    def add_threshold_constraint(self, threshold):
        # Interprète threshold comme une fraction minimale de lignes/colonnes sélectionnées
        # et impose aussi un minimum de 2 pour éviter les solutions triviales
        self.remove_threshold_constraints()
        m = len(self.rows_data)
        n = len(self.cols_data)
        min_rows = max(2, math.ceil(float(threshold) * m)) if threshold is not None else 2
        min_cols = max(2, math.ceil(float(threshold) * n)) if threshold is not None else 2
        c1 = self.model.addConstr(gp.quicksum(v for v, _ in self.lp_rows.values()) >= min_rows, "row_threshold")
        c2 = self.model.addConstr(gp.quicksum(v for v, _ in self.lp_cols.values()) >= min_cols, "col_threshold")
        self._threshold_constrs = [c1, c2]

    def remove_threshold_constraints(self):
        for c in self._threshold_constrs:
            self.model.remove(c)
        self._threshold_constrs = []

    def _add_matrix_structure_constraints(self):
        for (row, col), (cell_var, cell_value) in self.lp_cells.items():
            self.model.addConstr(self.lp_rows[row][0] >= cell_var, f'cell_{row}_{col}_1')
            self.model.addConstr(self.lp_cols[col][0] >= cell_var, f'cell_{row}_{col}_2')
            self.model.addConstr(self.lp_rows[row][0] + self.lp_cols[col][0] - 1 <= cell_var, f'cell_{row}_{col}_3')

    def add_density_constraint(self, density):
        # Ajoute contrainte: erreurs <= density * total
        # Erreurs = sélection de cellules dont la vraie valeur est 0
        if self._err_expr is None or self._tot_expr is None:
            self._prepare_core_expressions()
        constr = self.model.addConstr(self._err_expr <= float(density) * self._tot_expr, "density_constraint")
        self._density_constrs.append(constr)

    def remove_density_constraints(self):
        for c in self._density_constrs:
            self.model.remove(c)
        self._density_constrs = []

    def update_density_constraints(self, density):
        self.remove_density_constraints()
        self.add_density_constraint(density)

    def _prepare_core_expressions(self):
        ones_terms = []
        err_terms = []
        tot_terms = []
        for (row, col), (cell_var, cell_value) in self.lp_cells.items():
            tot_terms.append(cell_var)
            if cell_value == 1:
                ones_terms.append(cell_var)
            else:
                err_terms.append(cell_var)
        self._ones_expr = gp.quicksum(ones_terms) if ones_terms else gp.LinExpr(0)
        self._err_expr = gp.quicksum(err_terms) if err_terms else gp.LinExpr(0)
        self._tot_expr = gp.quicksum(tot_terms) if tot_terms else gp.LinExpr(0)

    def build(self, delta=0.0):
        if self._is_built:
            raise RuntimeError("Model already built")
        self._prepare_core_expressions()
        self.model.setObjective(self._ones_expr, GRB.MAXIMIZE)
        self._add_matrix_structure_constraints()
        if self._tot_expr is None or self._err_expr is None:
            self._prepare_core_expressions()
        self.model.addConstr(self._err_expr <= float(self.error_rate) * self._tot_expr, name='err_rate')
        if delta is not None and float(delta) > 0:
            self.add_density_constraint(delta)
        self._is_built = True

    def optimize(self):
        self.model.optimize()

    def getVars(self):
        return self.model.getVars()

    def setParam(self, param, value):
        self.model.setParam(param, value)

    @property
    def objVal(self):
        return self.model.ObjVal

    @property
    def status(self):
        return self.model.Status

    def get_selected_rows(self):
        return [int(v.VarName.split('_')[1]) for v in self.getVars() if v.VarName.startswith('row_') and v.X > 0.5]

    def get_selected_cols(self):
        return [int(v.VarName.split('_')[1]) for v in self.getVars() if v.VarName.startswith('col_') and v.X > 0.5]

    def reset_model(self):
        self.model.remove(self.model.getConstrs())
        self.model.remove(self.model.getVars())
        self.lp_rows = self._create_row_variables()
        self.lp_cols = self._create_col_variables()
        self.lp_cells = self._create_cell_variables()
        self._forced_row_constrs = {}
        self._forced_col_constrs = {}
        self._improvement_constr = None
        self._density_constrs = []
        self._threshold_constrs = []
        self._is_built = False


# Backward-compat alias if someone imports old name
max_one_grb_v2 = MaxOneModel