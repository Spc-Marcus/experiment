import numpy as np
import random
import csv
from itertools import permutations

def create_simple_matrix(rows, cols) -> list[list[int]]:
    '''
    Crée une matrice binaire de taille rows x cols
    telle que toutes les lignes ET toutes les colonnes sont uniques.
    Méthode rapide : tire aléatoirement des lignes uniques, puis vérifie les colonnes.
    '''
    if rows > 2 ** cols or cols > 2 ** rows:
        raise ValueError("Impossible de garantir unicité des lignes et colonnes avec ces dimensions.")

    all_possible_rows = [list(map(int, format(i, f'0{cols}b'))) for i in range(2 ** cols)]
    tries = 10000
    for _ in range(tries):
        candidate_rows = random.sample(all_possible_rows, rows)
        columns = list(zip(*candidate_rows))
        if len({tuple(col) for col in columns}) == cols:
            return candidate_rows
    raise ValueError("Impossible de générer une matrice avec lignes et colonnes toutes uniques pour ces dimensions après plusieurs essais.")

def extend_matrix(matrix: list[list[int]], size_cols: list[int], size_rows: list[int]) -> np.ndarray:
    '''
    Etend une matrice binaire tel que chaque colonne est multiplier de [min_cols, max_cols]
    et chaque ligne est multiplier de [min_rows, max_rows].
    '''
    min_cols, max_cols = size_cols
    min_rows, max_rows = size_rows
    
    # Étendre les lignes
    extended_rows = []
    for row in matrix:
        nb_copies = random.randint(min_rows, max_rows)
        for _ in range(nb_copies):
            extended_rows.append(row[:])
    
    # Étendre les colonnes
    if not extended_rows:
        return np.array([])
    
    num_cols = len(extended_rows[0])
    col_copies = [random.randint(min_cols, max_cols) for _ in range(num_cols)]
    
    final_matrix = []
    for row in extended_rows:
        new_row = []
        for j, col_val in enumerate(row):
            new_row.extend([col_val] * col_copies[j])
        final_matrix.append(new_row)
    
    return np.array(final_matrix)

def add_noise_to_matrix(matrix: np.ndarray, error_rate: float) -> np.ndarray:
    """Ajoute du bruit à une matrice selon un taux d'erreur"""
    noisy_matrix = matrix.copy()
    n_errors = int(error_rate * matrix.size)
    
    # Choisir des positions aléatoirement et inverser leur valeur
    flat_indices = np.random.choice(matrix.size, n_errors, replace=False)
    row_indices, col_indices = np.unravel_index(flat_indices, matrix.shape)
    
    for r, c in zip(row_indices, col_indices):
        noisy_matrix[r, c] = 1 - noisy_matrix[r, c]
    
    return noisy_matrix

def mix_matrices(matrix1: np.ndarray) -> np.ndarray:
    """
    Mélange les valeurs de la matrice en permutant aléatoirement les lignes et les colonnes.
    """
    shuffled_matrix = matrix1.copy()
    
    # Mélanger les lignes
    np.random.shuffle(shuffled_matrix)
    
    # Mélanger les colonnes
    shuffled_matrix = shuffled_matrix[:, np.random.permutation(shuffled_matrix.shape[1])]
    
    return shuffled_matrix

def save_matrix_with_labels(matrix, filepath: str):
    """Sauvegarde une matrice (list[list[int]] ou np.ndarray) dans un CSV
    avec noms de lignes rX et colonnes cY.
    Format:
        ,c0,c1,c2,...
        r0,0,1,0,...
        r1,1,0,1,...
    """
    if not isinstance(matrix, np.ndarray):
        matrix = np.array(matrix)
    if matrix.size == 0:
        raise ValueError("Matrice vide, rien à sauvegarder.")
    n_rows, n_cols = matrix.shape
    header = [''] + [f'c{j}' for j in range(n_cols)]
    with open(filepath, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for i, row in enumerate(matrix):
            writer.writerow([f'r{i}'] + list(row))
    return filepath

if __name__ == "__main__":
    # Exemple d'utilisation
    rows = 12  # Changé pour respecter la contrainte
    cols = 10
    matrix = create_simple_matrix(rows, cols)
    for row in matrix:
        print(row)
    
    size_cols = [5, 30]  # Exemple de taille de colonnes
    size_rows = [5, 30]  # Exemple de taille de lignes
    extended_matrix = extend_matrix(matrix, size_cols, size_rows)
    print("\nMatrice étendue :")
    print(extended_matrix)
    print("\nTaille de la matrice étendue :", extended_matrix.shape)
    # Ajouter du bruit à la matrice
    error_rate = 0.025  # 1% de bruit
    noisy_matrix = add_noise_to_matrix(extended_matrix, error_rate)
    print("\nMatrice avec bruit :")
    print(noisy_matrix)
    # Mélanger la matrice
    mixed_matrix = mix_matrices(noisy_matrix)
    print("\nMatrice mélangée :")
    print(mixed_matrix)

    # Sauvegardes
    save_matrix_with_labels(matrix, 'matrix_initial.csv')
    save_matrix_with_labels(extended_matrix, 'matrix_etendue.csv')
    save_matrix_with_labels(noisy_matrix, 'matrix_bruitee.csv')
    save_matrix_with_labels(mixed_matrix, 'matrix_melangee.csv')
    print("\nMatrices sauvegardées dans: matrix_initial.csv, matrix_etendue.csv, matrix_bruitee.csv, matrix_melangee.csv")