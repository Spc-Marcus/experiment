"""
Simple stats module for compatibility.
"""

class StatsTracker:
    def __init__(self):
        self.enabled = False
    
    def record_matrix_processing(self, **kwargs):
        pass
    
    def record_clustering_step(self, **kwargs):
        pass
    
    def record_preprocessing(self, **kwargs):
        pass

# Global stats instance
_stats = StatsTracker()

def get_stats():
    """Return global stats tracker."""
    return _stats
