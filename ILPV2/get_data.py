import pandas as pd
import numpy as np
from pathlib import Path
from typing import Tuple, List, Any

def get_data(csv_file_path: str) -> Tuple[List, List, List, List, List, List, pd.DataFrame, pd.DataFrame]:
    """
    Load and process CSV matrix data.
    
    Parameters
    ----------
    csv_file_path : str
        Path to the CSV file containing the binary matrix
        
    Returns
    -------
    tuple
        (rows_data, cols_data, edges_0, edges_1, row_names, col_names, df, comp_df)
    """
    
    # Load CSV file
    df = pd.read_csv(csv_file_path, index_col=0)
    
    # Ensure binary values (0, 1)
    df = df.fillna(0).astype(int)
    df = df.clip(0, 1)  # Ensure only 0 and 1 values
    
    # Extract basic information
    row_names = list(df.index)
    col_names = list(df.columns)
    
    # Create simple data structures for compatibility
    rows_data = list(range(len(row_names)))
    cols_data = list(range(len(col_names)))
    
    # Create edge lists (for compatibility with existing code)
    edges_0 = []  # Positions with value 0
    edges_1 = []  # Positions with value 1
    
    for i, row_name in enumerate(row_names):
        for j, col_name in enumerate(col_names):
            if df.iloc[i, j] == 1:
                edges_1.append((i, j))
            else:
                edges_0.append((i, j))
    
    # Create complement dataframe (inverted values)
    comp_df = 1 - df
    
    return rows_data, cols_data, edges_0, edges_1, row_names, col_names, df, comp_df
