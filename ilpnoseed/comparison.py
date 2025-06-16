import numpy as np
import sys
import time
import logging
from pathlib import Path
from typing import List, Tuple, Dict, Any

# Import noseed version (current)
from .clustering import clustering_full_matrix as clustering_noseed
from .postprocess import post_processing

# Import seed version with careful path handling
try:
    # Try to import from the seed version
    sys.path.insert(0, str(Path(__file__).parent.parent / "ilphaplo"))
    from clustering_seed import clustering_full_matrix as clustering_seed
    from postprocess import post_processing as post_processing_seed
except ImportError:
    try:
        # Alternative path
        sys.path.insert(0, str(Path(__file__).parent.parent))
        from ilphaplo.clustering_seed import clustering_full_matrix as clustering_seed
        from ilphaplo.postprocess import post_processing as post_processing_seed
    except ImportError:
        try:
            # Try from current directory - use the seed version we created
            from .clustering_seed import clustering_full_matrix as clustering_seed
            # Use same postprocess for both if seed version not found
            post_processing_seed = post_processing
        except ImportError:
            # Fallback - create a dummy function
            def clustering_seed(*args, **kwargs):
                print("WARNING: clustering_seed not available")
                return []
            post_processing_seed = post_processing

logger = logging.getLogger(__name__)

def compare_final_clusters(clusters1: List[np.ndarray], clusters2: List[np.ndarray], strict: bool = True) -> bool:
    """
    Compare two sets of final clusters.
    
    Parameters
    ----------
    clusters1, clusters2 : List[np.ndarray]
        Lists of clusters (each cluster is array of read names)
    strict : bool
        If True, require exact match. If False, allow for minor differences.
        
    Returns
    -------
    bool
        True if clusters are identical/equivalent, False otherwise
    """
    
    if len(clusters1) != len(clusters2):
        return False
    
    if len(clusters1) == 0:
        return True  # Both empty
    
    # Convert to sets of frozensets for easier comparison
    set1 = {frozenset(cluster.tolist()) for cluster in clusters1}
    set2 = {frozenset(cluster.tolist()) for cluster in clusters2}
    
    if strict:
        return set1 == set2
    else:
        # Allow for small differences (e.g., 1-2 reads difference per cluster)
        if len(set1) != len(set2):
            return False
        
        # For each cluster in set1, find best match in set2
        threshold = 0.1  # 10% difference allowed
        
        for cluster1 in set1:
            best_match_ratio = 0
            for cluster2 in set2:
                if len(cluster1) == 0 and len(cluster2) == 0:
                    ratio = 1.0
                else:
                    common = len(cluster1 & cluster2)
                    total = len(cluster1 | cluster2)
                    ratio = common / total if total > 0 else 0
                
                best_match_ratio = max(best_match_ratio, ratio)
            
            if best_match_ratio < (1 - threshold):
                return False
        
        return True

def analyze_cluster_differences(clusters1: List[np.ndarray], clusters2: List[np.ndarray]) -> List[str]:
    """
    Analyze differences between two sets of clusters with detailed metrics.
    
    Returns
    -------
    List[str]
        List of difference descriptions
    """
    
    differences = []
    
    # Basic counts
    differences.append(f"Cluster count: NOSEED={len(clusters1)}, SEED={len(clusters2)}")
    
    if len(clusters1) > 0:
        noseed_sizes = [len(cluster) for cluster in clusters1]
        differences.append(f"NOSEED cluster sizes: {sorted(noseed_sizes)}")
        differences.append(f"NOSEED total reads: {sum(noseed_sizes)}")
    
    if len(clusters2) > 0:
        seed_sizes = [len(cluster) for cluster in clusters2]
        differences.append(f"SEED cluster sizes: {sorted(seed_sizes)}")
        differences.append(f"SEED total reads: {sum(seed_sizes)}")
    
    # Detailed cluster analysis if both have clusters
    if len(clusters1) > 0 and len(clusters2) > 0:
        # Find overlaps
        for i, cluster1 in enumerate(clusters1):
            set1 = set(cluster1.tolist())
            best_overlap = 0
            best_j = -1
            
            for j, cluster2 in enumerate(clusters2):
                set2 = set(cluster2.tolist())
                overlap = len(set1 & set2)
                if overlap > best_overlap:
                    best_overlap = overlap
                    best_j = j
            
            if best_j >= 0:
                set2 = set(clusters2[best_j].tolist())
                common = len(set1 & set2)
                only_noseed = len(set1 - set2)
                only_seed = len(set2 - set1)
                
                differences.append(f"Cluster {i+1} vs Cluster {best_j+1}: {common} common, {only_noseed} only NOSEED, {only_seed} only SEED")
    
    return differences

def calculate_detailed_cluster_metrics(clusters1: List[np.ndarray], clusters2: List[np.ndarray]) -> Dict[str, Any]:
    """
    Calculate detailed metrics about cluster differences.
    
    Returns
    -------
    Dict[str, Any]
        Detailed metrics about cluster differences
    """
    
    metrics = {
        'reads_exchanged': 0,
        'reads_unique_noseed': 0,
        'reads_unique_seed': 0,
        'reads_lost_noseed': 0,
        'reads_lost_seed': 0,
        'cluster_stability_score': 0.0,
        'read_assignment_changes': 0,
        'perfect_cluster_matches': 0,
        'partial_cluster_matches': 0,
        'completely_different_clusters': 0
    }
    
    if len(clusters1) == 0 and len(clusters2) == 0:
        metrics['cluster_stability_score'] = 1.0
        return metrics
    
    if len(clusters1) == 0 or len(clusters2) == 0:
        # One version has no clusters
        if len(clusters1) > 0:
            metrics['reads_lost_seed'] = sum(len(c) for c in clusters1)
        if len(clusters2) > 0:
            metrics['reads_lost_noseed'] = sum(len(c) for c in clusters2)
        return metrics
    
    # Get all reads from both versions
    all_reads_noseed = set()
    all_reads_seed = set()
    
    for cluster in clusters1:
        all_reads_noseed.update(cluster.tolist())
    
    for cluster in clusters2:
        all_reads_seed.update(cluster.tolist())
    
    # Calculate unique reads
    metrics['reads_unique_noseed'] = len(all_reads_noseed - all_reads_seed)
    metrics['reads_unique_seed'] = len(all_reads_seed - all_reads_noseed)
    
    # Calculate lost reads (reads that were clustered in one version but not the other)
    common_reads = all_reads_noseed & all_reads_seed
    metrics['reads_lost_noseed'] = len(all_reads_noseed - common_reads)
    metrics['reads_lost_seed'] = len(all_reads_seed - common_reads)
    
    # Analyze cluster assignments for common reads
    noseed_assignments = {}  # read -> cluster_index
    seed_assignments = {}    # read -> cluster_index
    
    for i, cluster in enumerate(clusters1):
        for read in cluster.tolist():
            if read in common_reads:
                noseed_assignments[read] = i
    
    for i, cluster in enumerate(clusters2):
        for read in cluster.tolist():
            if read in common_reads:
                seed_assignments[read] = i
    
    # Count reads that changed cluster assignment
    assignment_changes = 0
    for read in common_reads:
        if read in noseed_assignments and read in seed_assignments:
            # We can't directly compare cluster indices, so we'll use a different approach
            # Find which other reads this read is clustered with in each version
            noseed_cluster_reads = set(clusters1[noseed_assignments[read]].tolist())
            seed_cluster_reads = set(clusters2[seed_assignments[read]].tolist())
            
            # Check if the read's cluster context has changed significantly
            noseed_context = noseed_cluster_reads & common_reads
            seed_context = seed_cluster_reads & common_reads
            
            # Calculate Jaccard similarity of cluster contexts
            if len(noseed_context | seed_context) > 0:
                similarity = len(noseed_context & seed_context) / len(noseed_context | seed_context)
                if similarity < 0.8:  # 80% similarity threshold
                    assignment_changes += 1
    
    metrics['read_assignment_changes'] = assignment_changes
    
    # Calculate reads exchanged (approximation based on assignment changes)
    metrics['reads_exchanged'] = assignment_changes
    
    # Analyze cluster-level similarities
    perfect_matches = 0
    partial_matches = 0
    no_matches = 0
    
    used_seed_clusters = set()
    
    for i, cluster1 in enumerate(clusters1):
        set1 = set(cluster1.tolist()) & common_reads
        if len(set1) == 0:
            continue
            
        best_similarity = 0
        best_j = -1
        
        for j, cluster2 in enumerate(clusters2):
            if j in used_seed_clusters:
                continue
                
            set2 = set(cluster2.tolist()) & common_reads
            if len(set2) == 0:
                continue
            
            # Calculate Jaccard similarity
            similarity = len(set1 & set2) / len(set1 | set2) if len(set1 | set2) > 0 else 0
            
            if similarity > best_similarity:
                best_similarity = similarity
                best_j = j
        
        if best_similarity > 0.95:
            perfect_matches += 1
            if best_j >= 0:
                used_seed_clusters.add(best_j)
        elif best_similarity > 0.7:
            partial_matches += 1
            if best_j >= 0:
                used_seed_clusters.add(best_j)
        else:
            no_matches += 1
    
    metrics['perfect_cluster_matches'] = perfect_matches
    metrics['partial_cluster_matches'] = partial_matches
    metrics['completely_different_clusters'] = no_matches
    
    # Calculate overall cluster stability score
    total_clusters = max(len(clusters1), len(clusters2))
    if total_clusters > 0:
        stability = (perfect_matches + 0.5 * partial_matches) / total_clusters
        metrics['cluster_stability_score'] = min(1.0, stability)
    
    return metrics

def should_show_detailed_comparison(comparison: Dict[str, Any]) -> bool:
    """
    Determine if a detailed comparison should be shown based on the significance of differences.
    
    Parameters
    ----------
    comparison : Dict[str, Any]
        Comparison results
        
    Returns
    -------
    bool
        True if detailed comparison should be shown
    """
    
    # Always show if there were errors
    if not comparison.get('noseed_success', True) or not comparison.get('seed_success', True):
        return True
    
    # Show if clustering is not equivalent
    if not comparison.get('clustering_equivalent', True):
        return True
    
    # Show if final clusters are not equivalent  
    if not comparison.get('clusters_equivalent', True):
        return True
    
    # Show if there are significant differences in execution
    noseed_time = comparison.get('noseed_time', 0)
    seed_time = comparison.get('seed_time', 0)
    if max(noseed_time, seed_time) > 0:
        time_ratio = max(noseed_time, seed_time) / max(min(noseed_time, seed_time), 0.001)
        if time_ratio > 2.0:  # One version is more than 2x slower
            return True
    
    # Show if there are step differences
    if len(comparison.get('step_differences', [])) > 0:
        return True
    
    # Show if there are significant cluster differences
    if len(comparison.get('cluster_differences', [])) > 3:  # More than basic stats
        return True
    
    return False

def compare_clustering_results(
    input_matrix: np.ndarray,
    regions: List[List[int]],
    steps: List[Tuple[List[int], List[int], List[int]]],
    min_row_quality: int = 5,
    min_col_quality: int = 3,
    error_rate: float = 0.025,
    filename: str = "unknown",
    haplotype_count: int = 0,
    verbose: bool = True
) -> Dict[str, Any]:
    """
    Compare clustering results between noseed and seed versions.
    Includes post-processing to compare final strain clusters.
    
    Parameters
    ----------
    Same as clustering_full_matrix parameters
    verbose : bool
        If True, print detailed comparison results
        
    Returns
    -------
    Dict containing comparison results and both outputs
    """
    
    if verbose:
        print(f"\n🔍 COMPARISON MODE ENABLED for {filename}")
    
    # Generate synthetic read names for post-processing
    read_names = [f"read_{i}" for i in range(input_matrix.shape[0])]
    
    # Run noseed version (clustering + post-processing)
    start_time = time.time()
    
    try:
        noseed_clustering_result = clustering_noseed(
            input_matrix=input_matrix,
            regions=regions,
            steps=steps,
            min_row_quality=min_row_quality,
            min_col_quality=min_col_quality,
            error_rate=error_rate,
            filename=filename,
            haplotype_count=haplotype_count
        )
        
        # Post-processing for noseed
        noseed_clusters = post_processing(
            matrix=input_matrix,
            steps=noseed_clustering_result,
            read_names=read_names,
            distance_thresh=0.1
        )
        
        noseed_time = time.time() - start_time
        noseed_success = True
        noseed_error = None
        
    except Exception as e:
        noseed_clustering_result = []
        noseed_clusters = []
        noseed_time = time.time() - start_time
        noseed_success = False
        noseed_error = str(e)
    
    # Run seed version (clustering + post-processing)
    start_time = time.time()
    
    try:
        seed_clustering_result = clustering_seed(
            input_matrix=input_matrix,
            regions=regions,
            steps=steps,
            min_row_quality=min_row_quality,
            min_col_quality=min_col_quality,
            error_rate=error_rate,
            filename=filename,
            haplotype_count=haplotype_count
        )
        
        # Post-processing for seed
        seed_clusters = post_processing_seed(
            matrix=input_matrix,
            steps=seed_clustering_result,
            read_names=read_names,
            distance_thresh=0.1
        )
        
        seed_time = time.time() - start_time
        seed_success = True
        seed_error = None
        
    except Exception as e:
        seed_clustering_result = []
        seed_clusters = []
        seed_time = time.time() - start_time
        seed_success = False
        seed_error = str(e)
    
    # Compare results
    comparison = {
        'filename': filename,
        'matrix_shape': input_matrix.shape,
        'noseed_success': noseed_success,
        'seed_success': seed_success,
        'noseed_error': noseed_error,
        'seed_error': seed_error,
        'noseed_steps': len(noseed_clustering_result) if noseed_success else 0,
        'seed_steps': len(seed_clustering_result) if seed_success else 0,
        'noseed_clusters': len(noseed_clusters) if noseed_success else 0,
        'seed_clusters': len(seed_clusters) if seed_success else 0,
        'noseed_time': noseed_time,
        'seed_time': seed_time,
        'noseed_clustering_result': noseed_clustering_result,
        'seed_clustering_result': seed_clustering_result,
        'noseed_final_clusters': noseed_clusters,
        'seed_final_clusters': seed_clusters,
        'clustering_identical': False,
        'clustering_equivalent': False,
        'clusters_identical': False,
        'clusters_equivalent': False,
        'step_differences': [],
        'cluster_differences': [],
        'detailed_metrics': {}
    }
    
    # Detailed comparison if both succeeded
    if noseed_success and seed_success:
        # 1. Compare clustering steps (existing logic)
        if len(noseed_clustering_result) == len(seed_clustering_result):
            clustering_identical = True
            for i, (noseed_step, seed_step) in enumerate(zip(noseed_clustering_result, seed_clustering_result)):
                noseed_reads1, noseed_reads0, noseed_cols = noseed_step
                seed_reads1, seed_reads0, seed_cols = seed_step
                
                # Sort for comparison
                noseed_reads1_sorted = sorted(noseed_reads1)
                noseed_reads0_sorted = sorted(noseed_reads0)
                noseed_cols_sorted = sorted(noseed_cols)
                seed_reads1_sorted = sorted(seed_reads1)
                seed_reads0_sorted = sorted(seed_reads0)
                seed_cols_sorted = sorted(seed_cols)
                
                if not (noseed_reads1_sorted == seed_reads1_sorted and 
                       noseed_reads0_sorted == seed_reads0_sorted and 
                       noseed_cols_sorted == seed_cols_sorted):
                    clustering_identical = False
                    comparison['step_differences'].append({
                        'step': i,
                        'noseed_reads1': noseed_reads1,
                        'noseed_reads0': noseed_reads0,
                        'noseed_cols': noseed_cols,
                        'seed_reads1': seed_reads1,
                        'seed_reads0': seed_reads0,
                        'seed_cols': seed_cols
                    })
            
            comparison['clustering_identical'] = clustering_identical
        
        # Check clustering equivalence
        if len(noseed_clustering_result) == len(seed_clustering_result):
            clustering_equivalent = True
            total_diff_threshold = 0.1
            
            for i, (noseed_step, seed_step) in enumerate(zip(noseed_clustering_result, seed_clustering_result)):
                noseed_total = len(noseed_step[0]) + len(noseed_step[1])
                seed_total = len(seed_step[0]) + len(seed_step[1])
                
                if abs(noseed_total - seed_total) / max(noseed_total, seed_total, 1) > total_diff_threshold:
                    clustering_equivalent = False
                    break
            
            comparison['clustering_equivalent'] = clustering_equivalent
        
        # 2. Compare final clusters
        clusters_identical = compare_final_clusters(noseed_clusters, seed_clusters, strict=True)
        clusters_equivalent = compare_final_clusters(noseed_clusters, seed_clusters, strict=False)
        
        comparison['clusters_identical'] = clusters_identical
        comparison['clusters_equivalent'] = clusters_equivalent
        
        # 3. Calculate detailed metrics
        detailed_metrics = calculate_detailed_cluster_metrics(noseed_clusters, seed_clusters)
        comparison['detailed_metrics'] = detailed_metrics
        
        # Analyze cluster differences
        if not clusters_identical:
            comparison['cluster_differences'] = analyze_cluster_differences(noseed_clusters, seed_clusters)
        
        # Display results based on significance
        if verbose:
            # Always show basic comparison summary
            print(f"\n🔍 COMPARISON: {filename}")
            print(f"  NOSEED: {len(noseed_clustering_result)} steps → {len(noseed_clusters)} clusters, {noseed_time:.3f}s")
            print(f"  SEED:   {len(seed_clustering_result)} steps → {len(seed_clusters)} clusters, {seed_time:.3f}s")
            
            # Only show detailed output if there are significant differences
            if should_show_detailed_comparison(comparison):
                print(f"\n{'='*80}")
                print(f"⚠️  DETAILED COMPARISON NEEDED: {filename}")
                print(f"Matrix shape: {input_matrix.shape}")
                print(f"Regions: {len(regions)}, Initial steps: {len(steps)}")
                print(f"{'='*80}")
                
                # Show NOSEED results
                print(f"\n--- NOSEED VERSION RESULTS ---")
                if noseed_success:
                    print(f"NOSEED CLUSTERING: {len(noseed_clustering_result)} steps found")
                    print(f"NOSEED CLUSTERS: {len(noseed_clusters)} final strain clusters")
                    print(f"NOSEED TIME: {noseed_time:.3f}s")
                else:
                    print(f"NOSEED ERROR: {noseed_error}")
                
                # Show SEED results
                print(f"--- SEED VERSION RESULTS ---")
                if seed_success:
                    print(f"SEED CLUSTERING: {len(seed_clustering_result)} steps found")
                    print(f"SEED CLUSTERS: {len(seed_clusters)} final strain clusters")
                    print(f"SEED TIME: {seed_time:.3f}s")
                else:
                    print(f"SEED ERROR: {seed_error}")
                
                # Show comprehensive comparison
                print(f"--- COMPREHENSIVE COMPARISON ---")
                print(f"Clustering steps - NOSEED: {len(noseed_clustering_result)}, SEED: {len(seed_clustering_result)}")
                print(f"Final clusters - NOSEED: {len(noseed_clusters)}, SEED: {len(seed_clusters)}")
                print(f"Execution time - NOSEED: {noseed_time:.3f}s, SEED: {seed_time:.3f}s")
                print(f"Clustering steps identical: {comparison['clustering_identical']}")
                print(f"Clustering steps equivalent: {comparison['clustering_equivalent']}")
                print(f"Final clusters identical: {comparison['clusters_identical']}")
                print(f"Final clusters equivalent: {comparison['clusters_equivalent']}")
                
                # Show detailed metrics
                if detailed_metrics:
                    print(f"\n--- DETAILED CLUSTER METRICS ---")
                    print(f"Reads exchanged between clusters: {detailed_metrics['reads_exchanged']}")
                    print(f"Reads unique to NOSEED: {detailed_metrics['reads_unique_noseed']}")
                    print(f"Reads unique to SEED: {detailed_metrics['reads_unique_seed']}")
                    print(f"Read assignment changes: {detailed_metrics['read_assignment_changes']}")
                    print(f"Cluster stability score: {detailed_metrics['cluster_stability_score']:.3f}")
                    print(f"Perfect cluster matches: {detailed_metrics['perfect_cluster_matches']}")
                    print(f"Partial cluster matches: {detailed_metrics['partial_cluster_matches']}")
                    print(f"Completely different clusters: {detailed_metrics['completely_different_clusters']}")
                
                # Show clustering step differences (if not too many)
                if comparison['step_differences'] and len(comparison['step_differences']) <= 3:
                    print(f"\n🔍 CLUSTERING STEP DIFFERENCES:")
                    for diff in comparison['step_differences']:
                        step_num = diff['step']
                        print(f"  Step {step_num + 1}: Groups differ in composition")
                
                # Show cluster differences (if not too many)
                if comparison['cluster_differences'] and len(comparison['cluster_differences']) <= 5:
                    print(f"\n🔍 FINAL CLUSTER DIFFERENCES:")
                    for diff in comparison['cluster_differences']:
                        print(f"  {diff}")
                
                print(f"{'='*80}")
            else:
                # Results are equivalent - just show brief status
                clustering_status = "✅" if comparison['clustering_identical'] else "≈"
                cluster_status = "✅" if comparison['clusters_identical'] else "≈"
                
                clustering_equiv = "identical" if comparison['clustering_identical'] else "equivalent"
                cluster_equiv = "identical" if comparison['clusters_identical'] else "equivalent"
                
                print(f"  {clustering_status} Clustering {clustering_equiv}, {cluster_status} Clusters {cluster_equiv}")
    
    elif verbose:
        # At least one version failed - show details
        print(f"\n{'='*60}")
        print(f"❌ EXECUTION ERRORS: {filename}")
        print(f"{'='*60}")
        
        if not noseed_success:
            print(f"NOSEED failed: {noseed_error}")
        if not seed_success:
            print(f"SEED failed: {seed_error}")
        
        print(f"{'='*60}")
    
    return comparison

def run_comparison_test(
    input_matrix: np.ndarray,
    regions: List[List[int]],
    steps: List[Tuple[List[int], List[int], List[int]]],
    **kwargs
) -> Dict[str, Any]:
    """
    Run a single comparison test and return results.
    
    This is a wrapper around compare_clustering_results for easier testing.
    """
    
    return compare_clustering_results(
        input_matrix=input_matrix,
        regions=regions,
        steps=steps,
        verbose=True,
        **kwargs
    )
