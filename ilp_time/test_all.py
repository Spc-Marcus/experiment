import time
import os
import csv
from typing import Tuple
import numpy as np
import gurobipy as gp
from create_matrix import create_simple_matrix, extend_matrix, add_noise_to_matrix, mix_matrices
from post_processing import post_processing
from ilp import clustering_full_matrix
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


def run_max_one(X: np.ndarray, inhomogeneous_regions: list[ int], error_rate: float,min_rows:int=3, min_cols:int=3):
    """Exécute le pipline avec max one v2"""
    t0 = time.time()
    steps,data = clustering_full_matrix(X, regions=inhomogeneous_regions, error_rate=error_rate, version=1, min_row_quality=min_rows, min_col_quality=min_cols)
    gp.disposeDefaultEnv()
    t1 = time.time()
    return {
        'time': t1 - t0,
        'steps': steps,
        'data': data
    }


def run_max_one_v2(X: np.ndarray, inhomogeneous_regions: list[ int], error_rate: float,min_rows:int=3, min_cols:int=3):
    """Exécute le pipline avec max one v2"""
    t0 = time.time()
    steps,data = clustering_full_matrix(X, regions=inhomogeneous_regions, error_rate=error_rate, version=2, min_row_quality=min_rows, min_col_quality=min_cols)
    gp.disposeDefaultEnv()
    t1 = time.time()
    return {
        'time': t1 - t0,
        'steps': steps,
        'data': data
    }
def run_max_e_r_v2(X: np.ndarray, inhomogeneous_regions: list[ int], error_rate: float,min_rows:int=3, min_cols:int=3): 
    t0 = time.time()
    steps,data = clustering_full_matrix(X, regions=inhomogeneous_regions, error_rate=error_rate, version=3, min_row_quality=min_rows, min_col_quality=min_cols)
    gp.disposeDefaultEnv()
    t1 = time.time()
    return {
        'time': t1 - t0,
        'steps': steps,
        'data': data
    }
def run_max_one_v1_2(X: np.ndarray, inhomogeneous_regions: list[ int], error_rate: float,min_rows:int=3, min_cols:int=3):
    """Exécute le pipline avec max one v1.2"""
    t0 = time.time()
    steps,data = clustering_full_matrix(X, regions=inhomogeneous_regions, error_rate=error_rate, version=4, min_row_quality=min_rows, min_col_quality=min_cols)
    gp.disposeDefaultEnv()
    t1 = time.time()
    return {
        'time': t1 - t0,
        'steps': steps,
        'data': data
    }
def clusters_equivalent(clusters_a: list[np.ndarray], clusters_b: list[np.ndarray]) -> bool:
    """Compare two clusterings by membership sets (ignoring cluster order)."""
    if clusters_a is None or clusters_b is None:
        return False
    a_sets = {tuple(sorted(c.tolist())) for c in clusters_a}
    b_sets = {tuple(sorted(c.tolist())) for c in clusters_b}
    return a_sets == b_sets


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
             min_rows, min_cols, max_rows, max_cols, csv_file=None, distance_thresh=None, use_preprocessing=False):
    """
    Compare max_one (classique) vs max_one_v2 (compactée) avec ou sans préprocessing.
    Note: paramètres `thresholds` et `distance_thresh` sont ignorés (compatibilité interface).
    use_preprocessing: bool, si True, applique le préprocessing pour filtrer les colonnes inhomogènes.
    """
    csv_file = csv_file or "results.csv"
    print(f"Début des tests avec {len(error_rates)} error_rates, {len(strips)} strips, {len(haplotypes)} haplotypes, preprocessing={'activé' if use_preprocessing else 'désactivé'}")

    # En-tête CSV
    header_cols = [
        "Error-Rate","Best-Threshold","Best-Distance","Strip","Haplotype",
        "Matrix-Rows","Matrix-Cols","Cols-ILP-Input",
        "Time-ILP-V1","ILP-Steps-V1","Strips-Found-V1","Final-Clusters-V1","Orphans-V1","Unused-Cols-V1",
        "Time-ILP-V2","ILP-Steps-V2","Strips-Found-V2","Final-Clusters-V2","Orphans-V2","Unused-Cols-V2",
        "Clusters-Equivalent"
    ]
    # Écrire l'en-tête une fois
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header_cols)

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

                        # Prétraitement (optionnel)
                        if use_preprocessing:
                            inhomogeneous_regions, steps_pre = pre_processing(
                                X, min_col_quality=min_cols, certitude=float(best_thr), error_rate=error_rate
                            )
                        else:
                            inhomogeneous_regions = list(range(n))
                            steps_pre = []
                        # Colonnes utilisées pour l'ILP
                        cols_ilp = inhomogeneous_regions if inhomogeneous_regions else list(range(n))
                        cols_ilp_count = len(cols_ilp)
                        res_v1 = run_max_one(X, inhomogeneous_regions, error_rate, min_rows=min_rows, min_cols=min_cols)
                        # max_one (classique) sur les colonnes retenues
                        res_v3 = run_max_e_r_v2(X,inhomogeneous_regions, error_rate, min_rows=min_rows, min_cols=min_cols)
                        # max_one_v2 (compactée) sur les colonnes retenues
                        #res_v2 = run_max_one_v2(X,inhomogeneous_regions, error_rate, min_rows=min_rows, min_cols=min_cols)
                        res_v2 = run_max_one_v1_2(X, inhomogeneous_regions, error_rate, min_rows=min_rows, min_cols=min_cols) 
                        # Post-processing pour chaque version
                        read_names = [f"r{i}" for i in range(m)]
                        dist_used = distance_thresh if distance_thresh is not None else float(best_dist)
                        steps1 = steps_pre+res_v1['steps']
                        steps2 = steps_pre+res_v2['steps']
                        clusters_v1, reduced_v1, orphans_v1, unused_cols_v1 = post_processing(
                            X, steps1, read_names, distance_thresh=dist_used, min_reads_per_cluster=min_rows
                        )
                        clusters_v2, reduced_v2, orphans_v2, unused_cols_v2 = post_processing(
                            X, steps2, read_names, distance_thresh=dist_used, min_reads_per_cluster=min_rows
                        )

                        # Collecter métriques pour CSV
                        row = [
                            error_rate, best_thr, best_dist, strip, haplotype,
                            m, n, cols_ilp_count,
                            res_v3['time'], res_v3['data'].get('nb_ilp_steps', -1), res_v3['data'].get('nb_strips_from_ilp', -1),
                            len(clusters_v1), len(orphans_v1), len(unused_cols_v1),
                            res_v2['time'], res_v2['data'].get('nb_ilp_steps', -1), res_v2['data'].get('nb_strips_from_ilp', -1),
                            len(clusters_v2), len(orphans_v2), len(unused_cols_v2),
                            clusters_equivalent(clusters_v1, clusters_v2)
                        ]
                        with open(csv_file, 'a', newline='') as f:
                            writer = csv.writer(f)
                            writer.writerow(row)

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
    nb_matrix_permutations = 100
    # Taille des matrices (paramètres d'extension)
    min_rows = 3
    min_cols = 3
    max_rows = 45
    max_cols = 35
    # Paramètres d'erreur
    thresholds = []  # ignoré
    error_rates = [0.01, 0.02, 0.03]
    distance_thresh = None  # ignoré
    # Dimensions de base (avant extension)
    strips = [4, 5, 6, 7]
    haplotypes = [5, 7, 9]
    # Lancer les tests
    test_all(thresholds, error_rates, strips, haplotypes, nb_matrix_permutations,
             min_rows, min_cols, max_rows, max_cols,
             csv_file="results.csv", distance_thresh=distance_thresh, use_preprocessing=False)
    print("Tests terminés. Les résultats sont enregistrés dans 'results.csv'.")