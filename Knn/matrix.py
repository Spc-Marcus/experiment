import pandas as pd
import numpy as np
import time
from sklearn.impute import KNNImputer
import logging

logger = logging.getLogger(__name__)

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
    """
    # Input validation
    if matrix.size == 0:
        return matrix.copy()
        
    if n_neighbors <= 0:
        raise ValueError(f"n_neighbors must be positive, got {n_neighbors}")
    
    # Check if there are any missing values to impute
    if not np.isnan(matrix).any():
        logger.debug("No missing values found, returning original matrix")
        return matrix.copy()
    
    # Check for sufficient non-missing data
    n_reads = matrix.shape[0]
    if n_neighbors >= n_reads:
        logger.warning(f"n_neighbors ({n_neighbors}) >= n_reads ({n_reads}), using {n_reads-1}")
        n_neighbors = max(1, n_reads - 1)
    
    try:
        # Initialize KNN imputer with specified parameters
        imputer = KNNImputer(n_neighbors=n_neighbors)
        
        # Perform imputation
        logger.debug(f"Imputing missing values using {n_neighbors} neighbors")
        imputed_matrix = imputer.fit_transform(matrix)
        
        # Log imputation statistics
        original_missing = np.isnan(matrix).sum()
        remaining_missing = np.isnan(imputed_matrix).sum()
        success_rate = (original_missing - remaining_missing) / max(1, original_missing) * 100
        
        logger.debug(f"Imputation completed: {original_missing} → {remaining_missing} missing values "
                   f"({success_rate:.1f}% success rate)")
        
        return imputed_matrix
        
    except Exception as e:
        logger.error(f"Imputation failed: {e}")
        raise RuntimeError(f"Could not impute missing values: {e}") from e
