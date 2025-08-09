import numpy as np
import random

def create_simple_matrix(rows, cols) -> list[list[int]]:
    '''
    Cree une matrice simple binaire de taille rows x cols
    tel que les lignes sont uniques et les colonnes sont uniques.
    '''
    max_unique_rows = 2 ** cols
    if rows > max_unique_rows:
        raise ValueError(f"Impossible: rows ({rows}) > 2^cols (2^{cols} = {max_unique_rows}). Limite fondamentale dépassée.")
    
    matrix = [[0 for _ in range(cols)] for _ in range(rows)]
    for i in range(rows):
        # Convertir i en binaire pour créer une ligne unique
        binary_repr = format(i, f'0{cols}b')
        for j in range(cols):
            matrix[i][j] = int(binary_repr[j])
    return matrix

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

if __name__ == "__main__":
    # Exemple d'utilisation
    rows = 3  # Changé pour respecter la contrainte
    cols = 3
    matrix = create_simple_matrix(rows, cols)
    for row in matrix:
        print(row)
    
    size_cols = [2, 5]  # Exemple de taille de colonnes
    size_rows = [3, 5]  # Exemple de taille de lignes
    extended_matrix = extend_matrix(matrix, size_cols, size_rows)
    print("\nMatrice étendue :")
    print(extended_matrix)
    # Ajouter du bruit à la matrice
    error_rate = 0.01  # 1% de bruit
    noisy_matrix = add_noise_to_matrix(extended_matrix, error_rate)
    print("\nMatrice avec bruit :")
    print(noisy_matrix)
    # Mélanger la matrice
    mixed_matrix = mix_matrices(noisy_matrix)
    print("\nMatrice mélangée :")
    print(mixed_matrix)