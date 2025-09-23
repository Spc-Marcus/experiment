import time
import os
import csv
from typing import Tuple
import numpy as np
from create_matrix import create_simple_matrix, extend_matrix, add_noise_to_matrix, mix_matrices
from model.max_one_grb import max_Ones_gurobi, max_Ones_comp_gurobi
from post_processing import post_processing
try:
    # Prétraitement réutilisé depuis efficience_ALL (utilisé en interne, sans l'exposer dans les colonnes)
    from efficience_ALL.pre_processing import pre_processing
except ImportError:
    # Fallback si l'import en espace de noms échoue (ex., exécution différente)
    from pre_processing import pre_processing  # type: ignore

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
    print(f"Début des tests avec {len(error_rates)} error_rates, {len(strips)} strips, {len(haplotypes)} haplotypes")

    # En-tête CSV
    header = (
        "Error-Rate,Best-Threshold,Best-Distance,Strip,Haplotype,Matrix-Cols,Matrix-Rows,Cols-ILP-Input,"
        "Time-MaxOne,Rows-MaxOne,Cols-MaxOne,Ones-MaxOne,Status-MaxOne,Obj-MaxOne,Final-Clusters-MaxOne,"
        "Time-MaxOneV2,Rows-MaxOneV2,Cols-MaxOneV2,Ones-MaxOneV2,Status-MaxOneV2,Obj-MaxOneV2,Final-Clusters-MaxOneV2,"
        "Equivalent-Final\n"
    )
    # Écrire l'en-tête une fois
    with open(csv_file, 'w') as f:
        f.write(header)

    # Charger meilleurs paramètres par error-rate depuis Param_best si disponible
    here = os.path.dirname(__file__)
    local_best = os.path.join(here, 'best_by_error_rate_ilp.csv')
    param_best = os.path.join(os.path.dirname(here), 'Param_best', 'best_by_error_rate_ilp.csv')
    best_csv = local_best if os.path.exists(local_best) else param_best
    best_map = load_best_by_error_rate(best_csv)

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

                        # Paramètres optimaux par error-rate
                        best_thr = best_map.get(error_rate, {}).get('threshold', None)
                        best_dist = best_map.get(error_rate, {}).get('distance', None)
                        if best_thr is None:
                            best_thr = 0.1
                        if best_dist is None:
                            best_dist = 0.1

                        # Prétraitement (interne)
                        inhomogeneous_regions, steps_pre = pre_processing(
                            X, min_col_quality=min_cols, certitude=float(best_thr), error_rate=error_rate
                        )
                        # Colonnes utilisées pour l'ILP
                        cols_ilp = inhomogeneous_regions if inhomogeneous_regions else list(range(n))
                        X_ilp = X[:, cols_ilp]
                        cols_ilp_count = len(cols_ilp)

                        # max_one (classique) sur les colonnes retenues
                        res_v1 = run_max_one(X_ilp, error_rate)

                        # max_one_v2 (compactée) sur les colonnes retenues
                        res_v2 = run_max_one_v2(X_ilp, error_rate)

                        # Post-processing: nombre d'haplotypes finaux (clusters) par modèle
                        def compute_final_clusters_count(sel_rows, sel_cols_local, dist):
                            try:
                                # Noms simples pour chaque read
                                read_names = [f"R{i}" for i in range(m)]
                                # Construire un step unique à partir de la bicluster ILP (indices colonnes globaux)
                                sel_rows_list = list(sel_rows) if isinstance(sel_rows, (list, tuple, np.ndarray)) else []
                                sel_cols_local_list = list(sel_cols_local) if isinstance(sel_cols_local, (list, tuple, np.ndarray)) else []
                                sel_cols_global = [cols_ilp[c] for c in sel_cols_local_list]
                                other_rows = [r for r in range(m) if r not in set(sel_rows_list)]
                                steps = list(steps_pre) + [(sel_rows_list, other_rows, sel_cols_global)]
                                # Choisir un seuil de distance (par défaut 0.1 si inconnu)
                                distance = dist if isinstance(dist, (int, float)) and not isinstance(dist, bool) else 0.1
                                clusters, _, orphan_reads_names, _ = post_processing(
                                    X, steps, read_names, distance_thresh=float(distance), min_reads_per_cluster=min_rows
                                )
                                # Ne compter que les clusters non vides
                                return len([c for c in clusters if len(c) > 0])
                            except Exception:
                                # Fallback simple si le post-traitement échoue
                                sel_rows_list = list(sel_rows) if isinstance(sel_rows, (list, tuple, np.ndarray)) else []
                                other_rows = [r for r in range(m) if r not in set(sel_rows_list)]
                                clusters = []
                                if len(sel_rows_list) >= (min_rows or 1):
                                    clusters.append(sel_rows_list)
                                if len(other_rows) >= (min_rows or 1):
                                    clusters.append(other_rows)
                                return len(clusters)

                        final_clusters_v1 = compute_final_clusters_count(res_v1['rows'], res_v1['cols'], best_dist)
                        final_clusters_v2 = compute_final_clusters_count(res_v2['rows'], res_v2['cols'], best_dist)

                        # Équivalence des résultats finaux (nb d'haplotypes)
                        equivalent_final = 1 if final_clusters_v1 == final_clusters_v2 else 0

                        line = (
                            f"{error_rate},{best_thr},{best_dist},{strip},{haplotype},{n},{m},{cols_ilp_count},"
                            f"{res_v1['time']:.6f},{len(res_v1['rows'])},{len(res_v1['cols'])},{res_v1['ones']},{res_v1['status']},{res_v1['obj'] if res_v1['obj'] is not None else ''},{final_clusters_v1},"
                            f"{res_v2['time']:.6f},{len(res_v2['rows'])},{len(res_v2['cols'])},{res_v2['ones']},{res_v2['status']},{res_v2['obj'] if res_v2['obj'] is not None else ''},{final_clusters_v2},{equivalent_final}\n"
                        )
                        with open(csv_file, 'a') as f:
                            f.write(line)

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
    error_rates = [0.03, 0.035, 0.04, 0.05]
    distance_thresh = None  # ignoré
    # Dimensions de base (avant extension)
    strips = [4, 5, 6, 7, 8]
    haplotypes = [5, 6, 7, 8, 9, 10]

    # Lancer les tests
    test_all(thresholds, error_rates, strips, haplotypes, nb_matrix_permutations,
             min_rows, min_cols, max_rows, max_cols,
             csv_file="results.csv", distance_thresh=distance_thresh)
    print("Tests terminés. Les résultats sont enregistrés dans 'results.csv'.")