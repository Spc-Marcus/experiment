import pandas as pd
import numpy as np
import time
from sklearn.cluster import FeatureAgglomeration
from sklearn.cluster import AgglomerativeClustering
from sklearn.impute import KNNImputer
import json
from pathlib import Path



def save_matrix_as_csv(matrix: np.ndarray, reads: list, positions: list, 
                       output_path: Path, contig_name: str, start_pos: int) -> None:
    """
    Save matrix as CSV file with read names as index and positions as columns.
    
    Parameters
    ----------
    matrix : np.ndarray
        The variant matrix to save
    reads : list
        List of read names corresponding to matrix rows
    positions : list
        List of genomic positions corresponding to matrix columns
    output_path : Path
        Directory where to save the CSV file
    contig_name : str
        Name of the contig being processed
    start_pos : int
        Starting position of the genomic window
    """
    if matrix.size == 0 or len(reads) == 0:
        return
        
    # Create DataFrame with proper labels
    df = pd.DataFrame(matrix, index=reads, columns=positions)
    
    # Generate output filename
    csv_filename = f"matrix_{contig_name}_{start_pos}.csv"
    csv_path = output_path / csv_filename
    
    # Save to CSV
    df.to_csv(csv_path, index=True)

def create_matrix(dict_of_sus_pos : dict, min_coverage_threshold :int =0.6, 
                  contig_size: int = 0, contig_name: str = "", start_pos: int = 0,
                  csv_output_dir: Path = None, cover_threshold: float = 0.6,
                  min_rows: int = 20, min_cols: int = 20) -> tuple:
    """
    Preprocess a dictionary of suspicious positions to create a filtered matrix.
    
    This function converts a dictionary of genomic positions and their read values
    into a filtered pandas DataFrame and numpy matrix. It applies filtering to ensure
    reads span the entire window and positions have sufficient coverage.
    
    Parameters:
    -----------
    dict_of_sus_pos : dict
        Dictionary containing suspicious positions, where keys are genomic positions
        and values are dictionaries mapping {read_name: binary_value}
        Example: {1000: {'read1': 1, 'read2': 0}, 1001: {'read1': 0, 'read2': 1}}
    min_coverage_threshold : float, optional
        Minimum coverage threshold for columns (0.0-1.0). Columns with coverage
        below this threshold will be removed. Default is 0.6 (60% coverage).
    cover_threshold : float, optional
        Coverage threshold for read spanning (0.0-1.0). Reads must have data in both
        the beginning and end thirds of positions. Default is 0.6 (60% coverage).
    min_rows : int, optional
        Minimum number of rows (reads) required for saving matrix. Default is 20.
    min_cols : int, optional
        Minimum number of columns (positions) required for saving matrix. Default is 20.
    
    Additional Parameters:
    ---------------------
    contig_size : int, optional
        Size of the contig being processed. Used to determine if CSV saving is needed.
    contig_name : str, optional
        Name of the contig for CSV filename generation.
    start_pos : int, optional
        Starting position for CSV filename generation.
    csv_output_dir : Path, optional
        Directory to save CSV files when conditions are met.
    
    Returns:
    --------
    tuple
        A tuple containing:
        - matrix (numpy.ndarray): Filtered matrix of genomic variants with shape
          (n_reads, n_positions). Empty array if no data passes filtering.
        - reads (list): List of read names corresponding to matrix rows after
          filtering. Empty list if no reads pass filtering.
    
    Notes:
    ------
    The function applies filtering steps:
    1. Removes reads that don't have sufficient coverage in the first third of positions
    2. Removes reads that don't have sufficient coverage in the last third of positions  
    3. Removes positions (columns) with insufficient overall coverage
    4. Only saves matrices that meet the minimum size requirements (min_rows x min_cols)
    
    This ensures that retained reads span the entire genomic window and that
    positions have enough data for reliable variant calling.
    """
    if not dict_of_sus_pos:
        return np.array([]), []
    
    # Create DataFrame from dictionary
    df = pd.DataFrame(dict_of_sus_pos)
    
    if df.empty:
        return np.array([]), []
    
    # Store original positions for CSV export
    original_positions = list(df.columns)
    
    variant_matrix = df.copy()
    
    for col in df.columns:
        position_data = df[col].dropna()  # Ignorer les NaN
        
        if len(position_data) == 0:
            continue
            
        allele_counts = position_data.value_counts()
        
        if len(allele_counts) == 0:
            continue
        elif len(allele_counts) == 1:
            # Position monomorphe - tous les reads ont la même valeur
            majority_allele = allele_counts.index[0]  # Store the single allele
            variant_matrix[col] = 1  # Tous identiques à la référence = 1
        else:
            majority_allele = allele_counts.index[0]  # Plus fréquent (référence)
            
            # CORRECTION: 1 = identique à la référence (majority), 0 = variant
            variant_matrix[col] = (df[col] == majority_allele).astype(int)
            
            # Conserver les NaN pour les positions manquantes
            variant_matrix.loc[df[col].isna(), col] = np.nan
            
    # Filter reads that span the window using the new cover_threshold
    if len(variant_matrix.columns) >= 3:
        # Check coverage in first third of positions
        first_third = variant_matrix.iloc[:, :len(variant_matrix.columns)//3]
        first_third_coverage = first_third.notna().sum(axis=1) / len(first_third.columns)
        reads_with_good_start = first_third_coverage >= cover_threshold
        
        # Check coverage in last third of positions  
        last_third = variant_matrix.iloc[:, 2*len(variant_matrix.columns)//3:]
        last_third_coverage = last_third.notna().sum(axis=1) / len(last_third.columns)
        reads_with_good_end = last_third_coverage >= cover_threshold
        
        # Keep only reads that meet both criteria
        reads_spanning_window = reads_with_good_start & reads_with_good_end
        variant_matrix = variant_matrix.loc[reads_spanning_window, :]
    
    # Filter columns with insufficient coverage
    if not variant_matrix.empty:
        variant_matrix = variant_matrix.dropna(axis=1, thresh=min_coverage_threshold * len(variant_matrix.index))
    
    # Extract filtered reads and positions
    reads = list(variant_matrix.index)
    filtered_positions = list(variant_matrix.columns)
    
    # Convert to numpy matrix
    matrix = variant_matrix.to_numpy() if not variant_matrix.empty else np.array([])
    
    # Save as CSV if conditions are met
    if (csv_output_dir is not None and 
        contig_size > 90000 and 
        matrix.shape[0] >= min_rows and 
        matrix.shape[1] >= min_cols):
        
        save_matrix_as_csv(matrix, reads, filtered_positions, 
                          csv_output_dir, contig_name, start_pos)
    
    return matrix, reads


