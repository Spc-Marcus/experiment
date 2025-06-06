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
                  csv_output_dir: Path = None) -> tuple:
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
    The function applies three main filtering steps:
    1. Removes reads that have no data in the first third of positions
    2. Removes reads that have no data in the last third of positions  
    3. Removes positions (columns) with insufficient coverage
    
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
            
    # Filter reads that span the window (existing logic)
    if len(variant_matrix.columns) >= 3:
        tmp_idx = variant_matrix.iloc[:, :len(variant_matrix.columns)//3].dropna(axis=0, how='all')
        variant_matrix = variant_matrix.loc[tmp_idx.index, :]
        
        tmp_idx = variant_matrix.iloc[:, 2*len(variant_matrix.columns)//3:].dropna(axis=0, how='all')
        variant_matrix = variant_matrix.loc[tmp_idx.index, :]
    
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
        matrix.size > 0 and
        matrix.shape[0] >= 20 and 
        matrix.shape[1] >= 20):
        
        save_matrix_as_csv(matrix, reads, filtered_positions, 
                          csv_output_dir, contig_name, start_pos)
    
    return matrix, reads


def impute_missing_values(matrix: np.ndarray, n_neighbors: int = 10) -> np.ndarray:
    """
    Impute missing values in a genomic variant matrix using K-Nearest Neighbors.
    
    This function fills missing or uncertain values in a variant matrix by leveraging
    patterns from similar reads. It uses the K-Nearest Neighbors algorithm to identify
    reads with similar variant patterns and imputes missing values based on the
    consensus of these neighbors.

    Parameters
    ----------
    matrix : np.ndarray
        Input matrix with shape (n_reads, n_positions) containing variant calls.
        Values can include:
        - 1: Reference allele or major variant
        - 0: Alternative allele or minor variant  
        - np.nan: Missing or uncertain values to be imputed
        - Other numerical values: Intermediate probabilities or quality scores
        
        Missing values (np.nan, None) will be replaced with estimated values
        based on K nearest neighbors.
        
    n_neighbors : int, optional
        Number of nearest neighbors to use for imputation. Must be positive
        and less than the number of reads with complete data.
        Default is 10.
        
    Returns
    -------
    np.ndarray
        Matrix with same shape as input but with missing values imputed.
        All np.nan values will be replaced with estimated values based on
        the K-nearest neighbors algorithm. Original non-missing values are
        preserved unchanged.
        
        **Value range**: Imputed values will be in the range of existing
        non-missing values in the dataset (typically 0-1 for binary variants).

    Raises
    ------
    ValueError
        If n_neighbors is not positive or exceeds the number of complete reads
    MemoryError
        If matrix is too large for available system memory
    RuntimeError
        If imputation fails due to insufficient non-missing data
        
    See Also
    --------
    binarize_matrix : Often applied after imputation to ensure binary values
    pre_processing : Comprehensive preprocessing pipeline including imputation
    sklearn.impute.KNNImputer : Underlying scikit-learn implementation
    """
    # Input validation
    if matrix.size == 0:
        return matrix.copy()
        
    if n_neighbors <= 0:
        raise ValueError(f"n_neighbors must be positive, got {n_neighbors}")
    
    # Check if there are any missing values to impute
    if not np.isnan(matrix).any():
        return matrix.copy()
    
    # Check for sufficient non-missing data
    n_reads = matrix.shape[0]
    if n_neighbors >= n_reads:
        n_neighbors = max(1, n_reads - 1)
    
    try:
        # Initialize KNN imputer with specified parameters
        imputer = KNNImputer(
            n_neighbors=n_neighbors)
        
        # Perform imputation
        imputed_matrix = imputer.fit_transform(matrix)
        
        return imputed_matrix
        
    except Exception as e:
        raise RuntimeError(f"Could not impute missing values: {e}") from e


def binarize_matrix(matrix: np.ndarray, certitude: float = 0.3, default: int = 0) -> np.ndarray:
    """
    Binarize a continuous variant matrix using certainty-based thresholds.
    
    This function converts continuous or probabilistic variant values into binary
    classifications (0/1) based on confidence thresholds. Values with high certainty
    are assigned definitive binary values, while uncertain values receive a default
    assignment. This is essential for clustering algorithms that require binary input.
    
    The binarization strategy uses symmetric thresholds around 0.5 to classify
    values as definitely 0, definitely 1, or uncertain. This approach preserves
    high-confidence variant calls while handling ambiguous cases consistently.
    
    Parameters
    ----------
    matrix : np.ndarray
        Input matrix with shape (n_reads, n_positions) containing continuous
        variant values. Typically contains:
        - Values near 0: Strong evidence for reference allele
        - Values near 1: Strong evidence for alternative allele  
        - Values near 0.5: Uncertain/ambiguous evidence
        - Values between 0-1: Probabilistic variant calls
        
        **Expected input range**: Most algorithms assume values in [0, 1] range,
        though the function handles arbitrary numeric values.
        
    certitude : float, optional
        Certainty threshold for binary classification (0.0 ≤ certitude < 0.5).
        Default is 0.3 (30% certainty requirement).
        
    default : int, optional
        Default value assigned to uncertain entries that fall between thresholds
        Default is 0 (conservative approach).
        
    Returns
    -------
    np.ndarray
        Binarized matrix with same shape as input. All values converted to:
        - 0: Reference allele or confident non-variant
        - 1: Alternative allele or confident variant
        - default: Uncertain values (as specified by default parameter)
        
        **Data type**: Returns integer array for efficient processing in
        clustering algorithms.
        
    Raises
    ------
    ValueError
        If certitude is not in valid range [0, 0.5)
    TypeError
        If matrix contains non-numeric values
        
    See Also
    --------
    impute_missing_values : Often applied before binarization
    pre_processing : Complete preprocessing pipeline including binarization
    """
    # Input validation
    if not 0.0 <= certitude < 0.5:
        raise ValueError(f"Certitude must be in range [0, 0.5), got {certitude}")
    
    if matrix.size == 0:
        return matrix.astype(int)
    
    # Check for non-numeric values
    if not np.issubdtype(matrix.dtype, np.number):
        raise TypeError("Matrix must contain numeric values")
    
    # Calculate thresholds
    lower_threshold = certitude
    upper_threshold = 1.0 - certitude
    
    # Apply three-way classification using numpy vectorized operations
    # This is equivalent to nested np.where but more readable
    result = np.full_like(matrix, default, dtype=int)
    result[matrix <= lower_threshold] = 0
    result[matrix >= upper_threshold] = 1
    
    return result