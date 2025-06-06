import numpy as np
import time
import sys
from pathlib import Path
from sklearn.cluster import FeatureAgglomeration, AgglomerativeClustering

try:
    from ..decorators.ilp_tracker import track_preprocessing
except ImportError:
    try:
        sys.path.insert(0, str(Path(__file__).parent.parent))
        from decorators.ilp_tracker import track_preprocessing
    except ImportError:
        def track_preprocessing():
            def decorator(func):
                return func
            return decorator

@track_preprocessing()
def pre_processing(input_matrix: np.ndarray, min_col_quality: int = 3) -> tuple:
    """
    Pre-processes the input matrix by identifying inhomogeneous regions.
    
    Steps performed:
    1. Splits the matrix into regions (groups of columns) using hierarchical clustering (FeatureAgglomeration) based on Hamming distance.
    2. Further splits large or noisy regions and identifies inhomogeneous regions based on the distribution of values.
    3. Optionally, splits homogeneous regions into two clusters using AgglomerativeClustering.
    
    Parameters
    ----------
    input_matrix : numpy.ndarray
        The input binary matrix (0,1) to be pre-processed, with shape (m, n).
    min_col_quality : int, optional
        Minimum number of columns required for a region to be considered significant (default is 3).
    
    Returns
    -------
    matrix : numpy.ndarray
        The input matrix (unchanged since already binary).
    inhomogenious_regions : list of list of int
        List of regions (each a list of column indices) identified as inhomogeneous.
    steps : list of tuple
        List of tuples describing the steps taken to split homogeneous regions. Each tuple contains:
            (cluster1_indices, cluster0_indices, region_columns)
    """
    start_time = time.time()
    m, n = input_matrix.shape
    matrix = input_matrix  # No processing needed, already binary
    
    # Init regions and steps
    regions = [list(range(n))]
    inhomogenious_regions = [list(range(n))]
    steps = []
    homogeneous_regions = []

    if m > 5 and n > 15:
        # Initial clustering of columns using FeatureAgglomeration
        agglo = FeatureAgglomeration(n_clusters=None, metric='hamming', linkage='complete', distance_threshold=0.35)
        agglo.fit(matrix)
        labels = agglo.labels_
        splitted_cols = {}
    
        # Group columns by their cluster labels
        for idx, label in enumerate(labels):
            if label in splitted_cols:
                splitted_cols[label].append(idx)
            else:
                splitted_cols[label] = [idx]
        
        regions = []        
        for cols in splitted_cols.values():
            # If the region is too large, remove the noisiest reads by splitting them with a strict distance
            # threshold, then take only the significant clusters
            if len(cols) > 15:
                matrix_reg = matrix[:, cols].copy()
                agglo = FeatureAgglomeration(n_clusters=None, metric='hamming', linkage='complete', distance_threshold=0.025)
                agglo.fit(matrix_reg)
                labels = agglo.labels_
                groups = {}

                for idx, label in enumerate(labels):
                    if label in groups:
                        groups[label].append(cols[idx])
                    else:
                        groups[label] = [cols[idx]]

                for c in groups.values():
                    if len(c) >= min_col_quality:
                        regions.append(c)
            elif len(cols) >= min_col_quality:
                regions.append(cols)
                
        # Identify inhomogeneous regions        
        inhomogenious_regions = []        
        for region in regions:
            thres = min_col_quality / len(region)
            matrix_reg = matrix[:, region].copy()
            x_matrix = matrix_reg.sum(axis=1) / len(region)  # Proportion of variants per read
            
            # Check if region is heterogeneous (many reads with intermediate proportions)
            if len(x_matrix[(x_matrix >= thres) * (x_matrix <= 1 - thres)]) > 10:
                inhomogenious_regions.append(region)
            else:
                # Cut into 2 regions of 1 and 0 for homogeneous regions
                homogeneous_regions.append(region)
                agglo = AgglomerativeClustering(n_clusters=2, metric='hamming', linkage='complete')
                agglo.fit(matrix_reg)  
                labels = agglo.labels_  
                cluster1, cluster0 = [], []
                
                for idx, label in enumerate(labels):
                    if label == 0:
                        cluster0.append(idx)
                    else:
                        cluster1.append(idx)
                
                # Calculate cluster averages
                mat_cl1 = matrix_reg[cluster1, :].sum() / (len(cluster1) * len(region))
                mat_cl0 = matrix_reg[cluster0, :].sum() / (len(cluster0) * len(region))

                # Only keep well-separated clusters (>90% or <10%)
                if (mat_cl0 > 0.9 or mat_cl0 < 0.1) and (mat_cl1 > 0.9 or mat_cl1 < 0.1):
                    steps.append((cluster1, cluster0, region))
    
    processing_time = time.time() - start_time
    
    # Record stats if enabled
    try:
        from .stats import get_stats
        stats = get_stats()
        if stats and stats.enabled:
            stats.record_preprocessing(
                input_shape=input_matrix.shape,
                output_shape=matrix.shape,
                execution_time=processing_time,
                regions_found=len(regions),
                inhomogeneous_regions=len(inhomogenious_regions),
                preprocessing_steps=len(steps)
            )
    except ImportError:
        pass
    
    return matrix, inhomogenious_regions, steps