import time
import os
import csv
from typing import Tuple
import numpy as np
from create_matrix import create_simple_matrix, extend_matrix, add_noise_to_matrix, mix_matrices
from model.max_one_grb import max_Ones_gurobi, max_Ones_comp_gurobi

def matrix(rows, cols, size_cols=None, size_rows=None, error_rate=None, randomize=False):
    """
    Crée une matrice de 0 et 1 avec des dimensions spécifiées.
    Si size_cols et size_rows sont fournis, la matrice est étendue en conséquence.
    Si error_rate est fourni, du bruit est ajouté à la matrice.
    """
    # Créer la matrice de base
    matrix = create_simple_matrix(rows, cols)
    
    # Étendre la matrice si des tailles de colonnes et de lignes sont spécifiées
    if size_cols and size_rows:
        matrix = extend_matrix(matrix, size_cols, size_rows)
    
    # Ajouter du bruit si un taux d'erreur est spécifié
    if error_rate is not None:
        matrix = add_noise_to_matrix(matrix, error_rate)
    
    # Mélanger la matrice si randomize est True
    if randomize:
        matrix = mix_matrices(matrix)
    
    # S'assurer que la matrice est un numpy array
    return np.array(matrix)


def _build_max_one_inputs(X: np.ndarray) -> Tuple[list[tuple[int, int]], list[tuple[int, int]], list[tuple[int, int]]]:
    """Prépare rows_data, cols_data, edges (positions à 1) pour les modèles max_one."""
    m, n = X.shape
    row_degrees = np.sum(X == 1, axis=1)
    col_degrees = np.sum(X == 1, axis=0)
    rows_data = [(int(r), int(row_degrees[r])) for r in range(m)]
    cols_data = [(int(c), int(col_degrees[c])) for c in range(n)]
    edges = []
    # edges = positions de 1 (cf. ilp_grb: edges ajoutés quand X[r,c]==1)
    ones = np.argwhere(X == 1)
    for r, c in ones:
        edges.append((int(r), int(c)))
    return rows_data, cols_data, edges


def run_max_one(X: np.ndarray, error_rate: float, mip_gap: float = 0.05, time_limit: int = 60):
    """Exécute max_Ones_gurobi sur toute la matrice et retourne temps et sélection."""
    rows_data, cols_data, edges = _build_max_one_inputs(X)
    model = max_Ones_gurobi(rows_data, cols_data, edges, error_rate)
    model.setParam('OutputFlag', 0)
    model.setParam('MIPGap', mip_gap)
    model.setParam('TimeLimit', time_limit)
    t0 = time.time()
    model.optimize()
    t1 = time.time()
    status = model.Status
    sel_rows, sel_cols = [], []
    if status in (2, 9):  # OPTIMAL=2, TIME_LIMIT=9
        for v in model.getVars():
            if v.VarName.startswith('row_') and v.X > 0.5:
                sel_rows.append(int(v.VarName.split('_')[1]))
            elif v.VarName.startswith('col_') and v.X > 0.5:
                sel_cols.append(int(v.VarName.split('_')[1]))
    # Calcul métriques
    ones = int(X[np.ix_(sel_rows, sel_cols)].sum()) if (sel_rows and sel_cols) else 0
    return {
        'time': t1 - t0,
        'status': status,
        'rows': sel_rows,
        'cols': sel_cols,
        'ones': ones,
        'obj': getattr(model, 'ObjVal', None),
    }


def run_max_one_v2(X: np.ndarray, error_rate: float, mip_gap: float = 0.05, time_limit: int = 60):
    """Exécute la version compactée (v2) via max_Ones_comp_gurobi et retourne temps et sélection."""
    rows_data, cols_data, edges = _build_max_one_inputs(X)
    model = max_Ones_comp_gurobi(rows_data, cols_data, edges, error_rate)
    model.setParam('OutputFlag', 0)
    model.setParam('MIPGap', mip_gap)
    model.setParam('TimeLimit', time_limit)
    t0 = time.time()
    model.optimize()
    t1 = time.time()
    status = model.Status
    sel_rows, sel_cols = [], []
    if status in (2, 9):  # OPTIMAL=2, TIME_LIMIT=9
        for v in model.getVars():
            if v.VarName.startswith('row_') and v.X > 0.5:
                sel_rows.append(int(v.VarName.split('_')[1]))
            elif v.VarName.startswith('col_') and v.X > 0.5:
                sel_cols.append(int(v.VarName.split('_')[1]))
    ones = int(X[np.ix_(sel_rows, sel_cols)].sum()) if (sel_rows and sel_cols) else 0
    return {
        'time': t1 - t0,
        'status': status,
        'rows': sel_rows,
        'cols': sel_cols,
        'ones': ones,
        'obj': getattr(model, 'ObjVal', None),
    }


def load_best_by_error_rate(csv_path: str) -> dict:
    """Charge le mapping Error-Rate -> (Threshold, Distance) depuis un CSV.
    Retourne un dict: {float(error_rate): {'threshold': float, 'distance': float}}
    """
    mapping = {}
    if not os.path.exists(csv_path):
        return mapping
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                er = float(row.get('Error-Rate'))
                thr = float(row.get('Threshold'))
                dist = float(row.get('Distance'))
            except Exception:
                continue
            mapping[er] = {'threshold': thr, 'distance': dist}
    return mapping

def test_all(thresholds, error_rates, strips, haplotypes, nb_matrix_permutations,
             min_rows, min_cols, max_rows, max_cols, csv_file=None, distance_thresh=None):
    """
    Compare max_one (classique) vs max_one_v2 (compactée) sans préprocessing.
    Note: paramètres `thresholds` et `distance_thresh` sont ignorés (compatibilité interface).
    """
    csv_file = csv_file or "results.csv"
    print(f"Début des tests (sans préprocessing) avec {len(error_rates)} error_rates, {len(strips)} strips, {len(haplotypes)} haplotypes")

    # En-tête CSV
    header = (
        "Error-Rate,Best-Threshold,Best-Distance,Strip,Haplotype,Matrix-Cols,Matrix-Rows,"
        "Time-MaxOne,Rows-MaxOne,Cols-MaxOne,Ones-MaxOne,Status-MaxOne,Obj-MaxOne,"
        "Time-MaxOneV2,Rows-MaxOneV2,Cols-MaxOneV2,Ones-MaxOneV2,Status-MaxOneV2,Obj-MaxOneV2\n"
    )
    open(csv_file, 'w').close()
    with open(csv_file, 'a') as f:
        f.write(header)
        f.flush()

        # Charger meilleurs paramètres par error-rate si disponible
        best_map = load_best_by_error_rate(os.path.join(os.path.dirname(__file__), 'best_by_error_rate_ilp.csv'))

        total_iterations = 0
        for error_rate in error_rates:
            for strip in strips:
                for haplotype in haplotypes:
                    print(f"Testing: error_rate={error_rate}, strip={strip}, haplotype={haplotype}")
                    for iteration in range(nb_matrix_permutations):
                        try:
                            # Générer matrice
                            X = matrix(haplotype, strip,
                                       size_cols=[min_cols, max_cols],
                                       size_rows=[min_rows, max_rows],
                                       error_rate=error_rate,
                                       randomize=True)
                            m, n = X.shape

                            # max_one (classique)
                            res_v1 = run_max_one(X, error_rate)

                            # max_one_v2 (compactée)
                            res_v2 = run_max_one_v2(X, error_rate)

                            best_thr = best_map.get(error_rate, {}).get('threshold', '')
                            best_dist = best_map.get(error_rate, {}).get('distance', '')

                            line = (
                                f"{error_rate},{best_thr},{best_dist},{strip},{haplotype},{n},{m},"
                                f"{res_v1['time']:.6f},{len(res_v1['rows'])},{len(res_v1['cols'])},{res_v1['ones']},{res_v1['status']},{res_v1['obj'] if res_v1['obj'] is not None else ''},"
                                f"{res_v2['time']:.6f},{len(res_v2['rows'])},{len(res_v2['cols'])},{res_v2['ones']},{res_v2['status']},{res_v2['obj'] if res_v2['obj'] is not None else ''}\n"
                            )
                            f.write(line)
                            f.flush()

                            total_iterations += 1
                            print(f"  Iteration {iteration+1}/{nb_matrix_permutations} completed")
                        except ValueError as e:
                            if "Limite fondamentale dépassée" in str(e):
                                print(f"  Skipping combination: {e}")
                                break
                            else:
                                print(f"  Error in iteration {iteration+1}: {e}")
                        except Exception as e:
                            print(f"  Error in iteration {iteration+1}: {e}")

        print(f"Total iterations completed: {total_iterations}")


if __name__ == "__main__":
    # Nb de matrices à tester (par combinaison strip/haplotype/error_rate)
    nb_matrix_permutations = 10
    # Taille des matrices (paramètres d'extension)
    min_rows = 3
    min_cols = 3
    max_rows = 20
    max_cols = 12
    # Paramètres d'erreur
    thresholds = []  # ignoré
    error_rates = [0.0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.05]
    distance_thresh = None  # ignoré
    # Dimensions de base (avant extension)
    strips = [4, 5, 6, 7, 8]
    haplotypes = [5, 6, 7, 8, 9, 10]

    # Lancer les tests
    test_all(thresholds, error_rates, strips, haplotypes, nb_matrix_permutations,
             min_rows, min_cols, max_rows, max_cols,
             csv_file="results.csv", distance_thresh=distance_thresh)
    print("Tests terminés. Les résultats sont enregistrés dans 'results.csv'.")