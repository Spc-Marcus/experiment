"""
Quick test to verify seeding parameter support in clustering functions.
"""

import sys
import inspect
from pathlib import Path

# Add current directory to path
current_dir = Path(__file__).parent
sys.path.insert(0, str(current_dir))

def test_seeding_support():
    """Test if seeding parameters are supported."""
    
    print("Testing seeding parameter support...")
    
    try:
        from clustering import clustering_full_matrix, find_quasi_biclique
        
        # Check clustering_full_matrix
        sig_clustering = inspect.signature(clustering_full_matrix)
        clustering_params = list(sig_clustering.parameters.keys())
        
        print(f"\nclustering_full_matrix parameters:")
        for param in clustering_params:
            print(f"  - {param}")
        
        has_X_factor = 'X_factor' in clustering_params
        has_step_n = 'step_n' in clustering_params
        
        print(f"\nSeeding parameter support in clustering_full_matrix:")
        print(f"  ✅ X_factor: {has_X_factor}")
        print(f"  ✅ step_n: {has_step_n}")
        
        # Check find_quasi_biclique
        sig_quasi = inspect.signature(find_quasi_biclique)
        quasi_params = list(sig_quasi.parameters.keys())
        
        print(f"\nfind_quasi_biclique parameters:")
        for param in quasi_params:
            print(f"  - {param}")
        
        has_X_factor_quasi = 'X_factor' in quasi_params
        has_step_n_quasi = 'step_n' in quasi_params
        
        print(f"\nSeeding parameter support in find_quasi_biclique:")
        print(f"  ✅ X_factor: {has_X_factor_quasi}")
        print(f"  ✅ step_n: {has_step_n_quasi}")
        
        # Overall status
        full_support = (has_X_factor and has_step_n and 
                       has_X_factor_quasi and has_step_n_quasi)
        
        print(f"\n{'='*50}")
        if full_support:
            print("✅ SEEDING FULLY SUPPORTED!")
            print("You can run seeding experiments with:")
            print("python seeding_experiment.py")
        else:
            print("❌ SEEDING NOT FULLY SUPPORTED")
            print("Missing parameters in clustering functions.")
            print("Seeding experiments will use default clustering.")
            
            if not (has_X_factor and has_step_n):
                print("\nTo fix: Add X_factor and step_n parameters to clustering_full_matrix()")
            if not (has_X_factor_quasi and has_step_n_quasi):
                print("To fix: Add X_factor and step_n parameters to find_quasi_biclique()")
        
        return full_support
        
    except ImportError as e:
        print(f"❌ Import error: {e}")
        return False
    except Exception as e:
        print(f"❌ Error: {e}")
        return False

if __name__ == "__main__":
    test_seeding_support()
