import gurobipy as gp
from gurobipy import GRB

class max_one_grb_v2:
    def __init__(self,rows_data, cols_data, edges, error_rate):
        self.rows_data = rows_data
        self.cols_data = cols_data
        self.edges = edges
        self.error_rate = error_rate
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
    
    def _create_row_variables(self):
        return {
            row: (self.model.addVar(vtype=GRB.BINARY, lb=0, ub=1, name=f'row_{row}'), degree) 
            for row, degree in self.rows_data
            }
    def _create_col_variables(self):
        return {
            col: (self.model.addVar(vtype=GRB.BINARY, lb=0, ub=1, name=f'col_{col}'), degree) 
            for col, degree in self.cols_data
            }
    def _create_cell_variables(self):
        lpCells = {}
        for row, _ in self.rows_data:
            for col, _ in self.cols_data:
                var = self.model.addVar(name=f'cell_{row}_{col}', vtype=GRB.BINARY, lb=0, ub=1)
                if (row, col) in self.edges:
                    lpCells[(row, col)] = (var, 1)
                else:
                    lpCells[(row, col)] = (var, 0)
        return lpCells
    
    def add_forced_rows_zero(self, rows):
        pass
    def remove_forced_rows_zero(self, rows):
        pass
    def add_forced_cols_zero(self, cols):
        pass
    def remove_forced_cols_zero(self, cols):
        pass
    def add_improvement_constraint(self, prev_obj):
        pass
    def remove_improvement_constraint(self):
        pass
    def add_threshold_constraint(self, threshold):
        pass
    def _add_matrix_structure_constraints(self):
        pass
    def add_density_constraint(self, density):
        pass
    def remove_density_constraints(self):
        pass
    def update_density_constraints(self, density):
        pass
    def build(self):
        pass
    def optimize(self):
        pass
    def getVars(self):
        pass
    def setParam(self, param, value):
        self.model.setParam(param, value)
    @property
    def objVal(self):
        return self.model.ObjVal
    @property
    def status(self):
        return self.model.Status
    def get_selected_rows(self):
        pass
    def get_selected_cols(self):
        pass
    def reset_model(self):
        pass