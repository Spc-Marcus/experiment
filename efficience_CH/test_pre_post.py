from create_matrix import create_simple_matrix, extend_matrix, add_noise_to_matrix, mix_matrices
from pre_processing import pre_processing
from post_processing import post_processing
import time

def matrix(rows, cols,size_cols=None, size_rows=None,error_rate=None, randomize=False):
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
    
    return matrix

def test_pre_post(thresholds, error_rates, strips, haplotypes, nb_matrix_permutations,min_rows, min_cols, max_rows, max_cols,csv_file=None):
    csv_file = csv_file or "results.csv"
    print(f"Début des tests avec {len(thresholds)} thresholds, {len(error_rates)} error_rates, {len(strips)} strips, {len(haplotypes)} haplotypes")
    
    open(csv_file, 'w').close()  # Clear the file before writing
    with open(csv_file, 'a') as f:
        f.write("Threshold,Error-Rate,Strip,Haplotype,Matrix-Size,Time-Pre-Processing,Time-Post-Processing,Steps-Count,Unused-Cols,Final-Clusters,Orphan-Reads\n")
        f.flush()
        
        total_iterations = 0
        for threshold in thresholds:
            for error_rate in error_rates:
                for strip in strips:
                    for haplotype in haplotypes:
                        print(f"Testing: threshold={threshold}, error_rate={error_rate}, strip={strip}, haplotype={haplotype}")
                        for iteration in range(nb_matrix_permutations):
                            try:
                                # Créer une matrice avec les paramètres spécifiés
                                matrix_data = matrix(haplotype, strip, size_cols=[min_cols, max_cols], size_rows=[min_rows, max_rows], error_rate=error_rate, randomize=True)
                                matrix_size = matrix_data.shape
                                read_names = [f"read_{i}" for i in range(matrix_size[0])]
                                
                                # Prétraitement de la matrice
                                pre_start_time = time.time()
                                inhomogeneous_regions, steps = pre_processing(matrix_data, min_col_quality=min_cols, certitude=threshold, error_rate=error_rate)
                                pre_end_time = time.time()
                                
                                # Post-traitement des strips
                                post_start_time = time.time()
                                clusters, reduced_matrix, orphan_reads, unused_columns = post_processing(
                                    matrix_data, steps, read_names, distance_thresh=0.1, min_reads_per_cluster=min_rows
                                )
                                post_end_time = time.time()
                                
                                # Enregistrer les résultats dans le fichier CSV
                                line = f"{threshold},{error_rate},{strip},{haplotype},{matrix_size},{pre_end_time - pre_start_time},{post_end_time - post_start_time},{len(steps)},{len(unused_columns)},{len(clusters)},{len(orphan_reads)}\n"
                                f.write(line)
                                f.flush()
                                total_iterations += 1
                                print(f"  Iteration {iteration+1}/{nb_matrix_permutations} completed")
                            except ValueError as e:
                                if "Limite fondamentale dépassée" in str(e):
                                    print(f"  Skipping combination: {e}")
                                    break  # Skip all iterations for this combination
                                else:
                                    print(f"  Error in iteration {iteration+1}: {e}")
                            except Exception as e:
                                print(f"  Error in iteration {iteration+1}: {e}")
        
        print(f"Total iterations completed: {total_iterations}")


if __name__ == "__main__":
    # Nb de matrice a tester par position
    nb_matrix_permutations = 25
    # Taille des matrices
    min_rows = 5
    min_cols = 4
    max_rows = 10
    max_cols = 9
    # Valeurs de seuil et d'erreur pour les tests
    thresholds = [0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.125,0.15,0.175,0.2]
    error_rates = [0.0,0.005,0.01,0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05,0.06,0.07,0.08,0.09,0.1]
    # Valeur de strip et d'haplotype pour les tests
    strips = [3,4,5,6,7,8,9]
    haplotypes = [4,5,6,7,8,9,10,11,12]

    # Lancer les tests
    test_pre_post(thresholds, error_rates, strips, haplotypes, nb_matrix_permutations, min_rows, min_cols, max_rows, max_cols, csv_file="results.csv")
    print("Tests terminés. Les résultats sont enregistrés dans 'results.csv'.")