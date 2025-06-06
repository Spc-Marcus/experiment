import pandas as pd
import itertools
from typing import List, Tuple, Any

def get_data(path: str):
    """
    Lit un fichier CSV et retourne les données de lignes, colonnes et arêtes
    """
    rows_data = []
    cols_data = []
    edges_0 = []
    edges_1 = []

    df = pd.read_csv(path, header=0, index_col=0)
    df[df == -1] = 0  # Ensure all values are 0 or 1

    rows = df.sum(axis=1)
    row_names = rows.index
    rows_data = list(zip(range(len(row_names)), rows))

    cols = df.sum(axis=0)
    col_names = cols.index
    cols_data = list(zip(range(len(col_names)), cols))

    df = df.reset_index(drop=True)
    df = df.T.reset_index(drop=True).T
    edges_1 = list(df[df == 1].stack().index)
    edges_0 = list(df[df == 0].stack().index)

    # Compute the complement matrix
    comp_df = 1 - df  # Flipping 0s to 1s and 1s to 0s

    return rows_data, cols_data, edges_0, edges_1, row_names, col_names, df, comp_df


def get_complement_edges(num_row: int, num_col: int, edges: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """
    Calcule les arêtes complémentaires
    """
    # Create a set of all possible edges using Cartesian product
    all_edges = set(itertools.product(range(num_row), range(num_col)))

    # Convert the original edges list to a set for efficient subtraction
    original_edges = set(edges)

    # Complement edges = all_edges - original_edges
    complement_edges = list(all_edges - original_edges)

    return complement_edges


def get_complement_rowcols(rows: List[Tuple[int, int]], 
                          cols: List[Tuple[int, int]], 
                          edges: List[Tuple[int, int]]) -> Tuple[List[Tuple[int, int]], 
                                                                List[Tuple[int, int]], 
                                                                List[Tuple[int, int]]]:
    """
    Calcule les degrés complémentaires des lignes et colonnes
    """
    cols_compl_map = {col: len(rows) - degree for col, degree in cols}
    rows_compl_map = {row: len(cols) - degree for row, degree in rows}
    cols_compl = [(col, degree) for col, degree in cols_compl_map.items()]
    rows_compl = [(row, degree) for row, degree in rows_compl_map.items()]
    edges_compl = [(r, c) for r, _ in rows_compl for c, _ in cols_compl if (r, c) not in edges]
    
    return rows_compl, cols_compl, edges_compl


def get_data_txt_file(path: str) -> Tuple[List[Tuple[int, int]], 
                                         List[Tuple[int, int]], 
                                         List[Tuple[int, int]], 
                                         range, 
                                         range]:
    """
    Lit un fichier texte avec un format spécifique et retourne les données
    """
    with open(path, 'r') as file:
        content = file.readlines()
    
    name = content[0][2:-1].strip()
    num_row = int(content[1].split(":")[1].strip())
    num_col = int(content[2].split(":")[1].strip())
    num_edge = int(content[3].split(":")[1].strip())

    deg_row = [0] * num_row
    deg_col = [0] * num_col
    edges = []
    df = pd.DataFrame(0, index=range(num_row), columns=range(num_col))

    for line in content[4:]:
        splitted_line = line.strip().split()
        u, v = int(splitted_line[0]), int(splitted_line[1])
        
        edges.append((u, v))
        deg_row[u] += 1
        deg_col[v] += 1
        # df.iloc[u, v] = 1 

    rows_data = list(zip(range(num_row), deg_row))
    cols_data = list(zip(range(num_col), deg_col))

    # Compute the complement matrix
    # comp_df = 1 - df # Flip 0s to 1s and 1s to 0s

    return rows_data, cols_data, edges, range(num_row), range(num_col)  # , df, comp_df


# Exemple d'utilisation
if __name__ == "__main__":
    # Test avec un fichier CSV
    try:
        rows, cols, edges_0, edges_1, row_names, col_names, df, comp_df = get_data("data.csv")
        print(f"CSV: {len(rows)} rows, {len(cols)} cols, {len(edges_1)} edges with value 1")
    except FileNotFoundError:
        print("Fichier CSV non trouvé")
    
    # Test avec un fichier texte
    try:
        rows, cols, edges, row_range, col_range = get_data_txt_file("data.txt")
        print(f"TXT: {len(rows)} rows, {len(cols)} cols, {len(edges)} edges")
        
        # Test des fonctions complémentaires
        complement_edges = get_complement_edges(len(rows), len(cols), edges)
        rows_compl, cols_compl, edges_compl = get_complement_rowcols(rows, cols, edges)
        
        print(f"Complement: {len(complement_edges)} edges, {len(edges_compl)} edges_compl")
    except FileNotFoundError:
        print("Fichier TXT non trouvé")