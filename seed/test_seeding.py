"""
Test script to verify seeding parameters are working correctly.
"""

import sys
import numpy as np
from pathlib import Path

# Add current directory to path
current_dir = Path(__file__).parent
sys.path.insert(0, str(current_dir))

def test_seeding_parameters():
    """Test that seeding parameters are properly passed through the pipeline."""
    
    print("Testing seeding parameter integration...")
    
    try:
        # Import functions
        from clustering import clustering_full_matrix, find_quasi_biclique
        from preprocess import pre_processing
        import inspect
        
        # Check function signatures
        print("\n=== CHECKING FUNCTION SIGNATURES ===")
        
        # Check clustering_full_matrix
        sig_clustering = inspect.signature(clustering_full_matrix)
        print(f"clustering_full_matrix parameters: {list(sig_clustering.parameters.keys())}")
        
        has_X_factor = 'X_factor' in sig_clustering.parameters
        has_step_n = 'step_n' in sig_clustering.parameters
        
        print(f"  - Has X_factor parameter: {has_X_factor}")
        print(f"  - Has step_n parameter: {has_step_n}")
        
        # Check find_quasi_biclique
        sig_quasi = inspect.signature(find_quasi_biclique)
        print(f"find_quasi_biclique parameters: {list(sig_quasi.parameters.keys())}")
        
        has_X_factor_quasi = 'X_factor' in sig_quasi.parameters
        has_step_n_quasi = 'step_n' in sig_quasi.parameters
        
        print(f"  - Has X_factor parameter: {has_X_factor_quasi}")
        print(f"  - Has step_n parameter: {has_step_n_quasi}")
        
        if not (has_X_factor and has_step_n and has_X_factor_quasi and has_step_n_quasi):
            print("\n❌ ERROR: Missing seeding parameters in function signatures!")
            return False
        
        print("\n✅ All function signatures have seeding parameters")
        
        # Test with small matrix
        print("\n=== TESTING WITH SMALL MATRIX ===")
        
        # Create test matrix
        np.random.seed(42)
        test_matrix = np.random.choice([0, 1], size=(20, 15), p=[0.3, 0.7])
        
        print(f"Test matrix shape: {test_matrix.shape}")
        print(f"Test matrix density: {np.sum(test_matrix) / test_matrix.size:.3f}")
        
        # Run preprocessing
        matrix, regions, steps = pre_processing(test_matrix, min_col_quality=3)
        print(f"Preprocessing: {len(regions)} regions, {len(steps)} steps")
        
        if len(regions) == 0:
            print("No regions found, creating artificial region for testing")
            regions = [list(range(min(10, test_matrix.shape[1])))]
        
        # Test clustering with different seeding parameters
        print("\n=== TESTING DIFFERENT SEEDING PARAMETERS ===")
        
        test_params = [
            (2, 2),
            (3, 4), 
            (4, 6),
            (5, 8)
        ]
        
        for X_factor, step_n in test_params:
            try:
                print(f"\nTesting X_factor={X_factor}, step_n={step_n}")
                
                result_steps = clustering_full_matrix(
                    matrix,
                    regions=regions,
                    steps=steps,
                    min_row_quality=3,
                    min_col_quality=3,
                    error_rate=0.025,
                    filename="test_matrix",
                    haplotype_count=2,
                    X_factor=X_factor,
                    step_n=step_n
                )
                
                print(f"  ✅ Success: Found {len(result_steps)} clustering steps")
                
            except Exception as e:
                print(f"  ❌ Error with X_factor={X_factor}, step_n={step_n}: {e}")
                return False
        
        print("\n✅ All seeding parameter tests passed!")
        return True
        
    except Exception as e:
        print(f"\n❌ Test failed with error: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_seeding_parameters()
    if success:
        print("\n🎉 Seeding parameter integration is working correctly!")
        print("You can now run the seeding experiments with:")
        print("python seeding_experiment.py")
    else:
        print("\n💥 Seeding parameter integration has issues that need to be fixed.")
        sys.exit(1)
