import numpy as np
import logging
import time
from typing import List, Tuple
from sklearn.cluster import AgglomerativeClustering

logger = logging.getLogger(__name__)

def post_processing(matrix: np.ndarray, steps: List[Tuple[List[int], List[int], List[int]]], 
                   read_names: List[str], distance_thresh: float = 0.1) -> List[np.ndarray]:
    """
    Post-process clustering results and create final strain clusters.
    
    This function takes the output of the biclustering algorithm and creates final
    read clusters representing different microbial strains. It processes hierarchical
    clustering steps, filters small clusters, assigns orphaned reads, and merges
    similar clusters to produce biologically meaningful strain groups.
    
    Parameters
    ----------
    matrix : np.ndarray
        Binary variant matrix with shape (n_reads, n_positions) containing
        processed variant calls (0, 1, or -1 for uncertain)
    steps : List[Tuple[List[int], List[int], List[int]]]
        List of clustering steps from biclustering algorithm. Each step contains:
        - reads1: List of read indices assigned to cluster 1
        - reads0: List of read indices assigned to cluster 0  
        - cols: List of column indices used for this clustering decision
    read_names : List[str]
        List of read names corresponding to matrix rows
    distance_thresh : float, optional
        Distance threshold for merging similar clusters using Hamming distance.
        Default is 0.1 (10% difference threshold)
        
    Returns
    -------
    List[np.ndarray]
        List of strain clusters, where each cluster is a numpy array of read names
        belonging to that strain. Empty list if no valid clusters found.
        
    Raises
    ------
    ValueError
        If matrix dimensions don't match read_names length
    RuntimeError
        If clustering post-processing fails
        
    See Also
    --------
    clustering_full_matrix : Previous step in the pipeline
    """
    start_time = time.time()
    
    try:
        if hasattr(steps, '__iter__') and not isinstance(steps, (list, tuple)):
            steps = list(steps)
            
        # Input validation
        if len(read_names) != matrix.shape[0]:
            raise ValueError(f"Matrix rows ({matrix.shape[0]}) don't match read_names length ({len(read_names)})")
        
        if matrix.size == 0 or len(read_names) == 0:
            logger.warning("Empty matrix or read names, returning empty clusters")
            return []
        
        logger.debug(f"Starting post-processing: {matrix.shape[0]} reads, {len(steps)} clustering steps")
        
        # Begin with all reads in the same group
        clusters = [list(range(len(read_names)))]
        
        # Go through each step and separate reads based on clustering decisions
        for step_idx, step in enumerate(steps):
            reads1, reads0, cols = step
            if hasattr(reads1, '__iter__') and not isinstance(reads1, (list, tuple, np.ndarray)):
                reads1 = list(reads1)
            if hasattr(reads0, '__iter__') and not isinstance(reads0, (list, tuple, np.ndarray)):
                reads0 = list(reads0)
            
            if len(reads1) == 0 or len(reads0) == 0:
                logger.debug(f"Skipping step {step_idx}: empty group (reads1={len(reads1)}, reads0={len(reads0)})")
                continue
                
            logger.debug(f"Processing step {step_idx}: {len(reads1)} vs {len(reads0)} reads, {len(cols)} columns")
            new_clusters = []
            
            # For each existing cluster, split it based on current step
            for cluster_idx, cluster in enumerate(clusters):
                clust1 = [c for c in cluster if c in reads1]
                clust0 = [c for c in cluster if c in reads0]
                
                # Only keep non-empty clusters
                if len(clust1) > 0:
                    new_clusters.append(clust1)
                if len(clust0) > 0:
                    new_clusters.append(clust0)
            
            clusters = new_clusters
        
        logger.debug(f"After hierarchical splitting: {len(clusters)} clusters")
        
        # Remove small clusters and collect orphaned reads
        min_cluster_size = 5
        orphaned_reads = []
        large_clusters = []
        
        for cluster_idx, cluster in enumerate(clusters):
            if len(cluster) <= min_cluster_size:
                orphaned_reads.extend(cluster)
                logger.debug(f"Orphaning small cluster {cluster_idx} with {len(cluster)} reads")
            else:
                large_clusters.append(cluster)
        
        clusters = large_clusters
        logger.info(f"Kept {len(clusters)} large clusters, orphaned {len(orphaned_reads)} reads")
        
        if len(clusters) == 0:
            logger.warning("No large clusters found after filtering")
            return []
        
        logger.debug("Calculating cluster means...")
        
        # Calculate mean vectors for each cluster
        cluster_means = []
        for cluster_idx, cluster in enumerate(clusters):
            logger.debug(f"Calculating mean for cluster {cluster_idx} with {len(cluster)} reads")
            
            # Create mean vector ignoring uncertain values (-1)
            cluster_matrix = matrix[cluster]
            cluster_mean = np.zeros(cluster_matrix.shape[1])
            
            for col in range(cluster_matrix.shape[1]):
                col_data = cluster_matrix[:, col]
                valid_data = col_data[col_data != -1]  # Exclude uncertain values
                if len(valid_data) > 0:
                    cluster_mean[col] = np.round(np.mean(valid_data))
                else:
                    cluster_mean[col] = 0  # Default for columns with only uncertain values
            
            cluster_means.append(cluster_mean)
        
        logger.debug("Assigning orphaned reads...")
        
        # Assign orphaned reads to closest clusters
        reassigned_count = 0
        for read_idx in orphaned_reads:
            if read_idx >= matrix.shape[0]:
                continue
                
            read_vector = matrix[read_idx]
            
            # Calculate distances to all cluster means
            distances = []
            for cluster_mean in cluster_means:
                # Only consider positions where read has definite calls (not -1)
                valid_positions = read_vector != -1
                if np.any(valid_positions):
                    read_valid = read_vector[valid_positions]
                    mean_valid = cluster_mean[valid_positions]
                    if len(read_valid) > 0:
                        dist = np.mean(read_valid != mean_valid)  # Hamming distance
                    else:
                        dist = 1.0  # Maximum distance if no valid positions
                else:
                    dist = 1.0  # Maximum distance if all uncertain
                distances.append(dist)
            
            # Assign to closest cluster if distance is reasonable
            if distances:
                min_dist_idx = np.argmin(distances)
                if distances[min_dist_idx] < 0.3:  # 30% difference threshold
                    clusters[min_dist_idx].append(read_idx)
                    reassigned_count += 1
        
        logger.info(f"Reassigned {reassigned_count} orphaned reads to existing clusters")
        
        logger.debug("Merging similar clusters...")
        
        # Merge similar clusters based on their mean vectors
        if len(clusters) > 1:
            # Recalculate means after reassignment
            final_cluster_means = []
            for cluster in clusters:
                cluster_matrix = matrix[cluster]
                cluster_mean = np.zeros(cluster_matrix.shape[1])
                
                for col in range(cluster_matrix.shape[1]):
                    col_data = cluster_matrix[:, col]
                    valid_data = col_data[col_data != -1]
                    if len(valid_data) > 0:
                        cluster_mean[col] = np.round(np.mean(valid_data))
                    else:
                        cluster_mean[col] = 0
                
                final_cluster_means.append(cluster_mean)
            
            # Hierarchical clustering of cluster means
            if len(final_cluster_means) > 1:
                agglo_cl = AgglomerativeClustering(
                    n_clusters=None,
                    metric='hamming',
                    linkage='complete',
                    distance_threshold=distance_thresh
                )
                cluster_labels = agglo_cl.fit_predict(final_cluster_means)
                
                # Merge clusters with same labels
                merged_clusters = {}
                for idx, label in enumerate(cluster_labels):
                    if label in merged_clusters:
                        merged_clusters[label].extend(clusters[idx])
                    else:
                        merged_clusters[label] = list(clusters[idx])
                
                clusters = list(merged_clusters.values())
                logger.info(f"Merged similar clusters to {len(clusters)} final strain groups")
        
        logger.debug("Converting to read names...")
        
        # Convert read indices to read names and sort
        result_clusters = []
        for cluster_idx, cluster in enumerate(clusters):
            logger.debug(f"Converting cluster {cluster_idx} with {len(cluster)} reads")
            
            if len(cluster) > 0:  # Only include non-empty clusters
                # CORRECTION: Ensure cluster is a list
                cluster_list = list(cluster) if not isinstance(cluster, list) else cluster
                
                # Filter valid indices and convert to read names
                valid_indices = [r for r in cluster_list if isinstance(r, (int, np.integer)) and 0 <= r < len(read_names)]
                
                if len(valid_indices) > 0:
                    cluster_names = [read_names[r] for r in valid_indices]
                    cluster_names_sorted = np.array(sorted(cluster_names))
                    result_clusters.append(cluster_names_sorted)
                    logger.debug(f"Final cluster {cluster_idx}: {len(cluster_names)} reads")
        
        logger.info(f"Post-processing completed: {len(result_clusters)} distinct strain clusters identified")
        
        return result_clusters
        
    except Exception as e:
        logger.critical(f"Post-processing failed: {e}")
        raise RuntimeError(f"Could not complete post-processing: {e}") from e