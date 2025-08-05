import numpy as np
import matplotlib.pyplot as plt
from src.MatrixStriper.pre_processing import pre_processing
import time

def create_perfect_strip_matrix(n_reads, n_positions, n_strips=3, min_cols_per_strip=5) -> np.ndarray:
    """Crée une matrice avec des strips parfaits (sans erreur)"""
    matrix = np.zeros((n_reads, n_positions), dtype=int)
    
    # Vérifier que nous avons assez de colonnes pour respecter le minimum
    if n_positions < n_strips * min_cols_per_strip:
        raise ValueError(f"Pas assez de colonnes ({n_positions}) pour {n_strips} strips avec minimum {min_cols_per_strip} colonnes chacun")
    
    # Diviser les reads en groupes pour chaque strip
    reads_per_strip = n_reads // n_strips
    positions_per_strip = max(min_cols_per_strip, n_positions // n_strips)
    
    for i in range(n_strips):
        start_read = i * reads_per_strip
        end_read = start_read + reads_per_strip
        start_pos = i * positions_per_strip
        end_pos = min(start_pos + positions_per_strip, n_positions)
        
        # Créer un strip parfait : première moitié des reads à 0, seconde moitié à 1
        half_reads = reads_per_strip // 2
        matrix[start_read + half_reads:end_read, start_pos:end_pos] = 1
    
    return matrix

def add_noise_to_matrix(matrix, error_rate):
    """Ajoute du bruit à une matrice selon un taux d'erreur"""
    noisy_matrix = matrix.copy()
    n_errors = int(error_rate * matrix.size)
    
    # Choisir des positions aléatoirement et inverser leur valeur
    flat_indices = np.random.choice(matrix.size, n_errors, replace=False)
    row_indices, col_indices = np.unravel_index(flat_indices, matrix.shape)
    
    for r, c in zip(row_indices, col_indices):
        noisy_matrix[r, c] = 1 - noisy_matrix[r, c]
    
    return noisy_matrix

def test_preprocessing_performance():
    """Test la performance du preprocessing avec différents niveaux d'erreur"""
    print("=== Test de Performance du Preprocessing ===\n")
    
    # Paramètres de test
    matrix_sizes = [(200, 150)]
    error_rates = [0.0,0.005,0.01,0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05,0.06,0.07,0.08,0.09,0.1]
    nb_strips = [3, 4, 5, 6, 7, 8]
    #nb_strips = [6]
    certitudes = [0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.125,0.15,0.175,0.2]
    n_tests_per_config = 20
    
    results = []
    
    for n_reads, n_positions in matrix_sizes:
        for n_strips in nb_strips:
            print(f"Test avec matrice {n_reads}x{n_positions}, {n_strips} strips")
            print("=" * 70)
            
            for error_rate in error_rates:
                print(f"\nError rate {error_rate:.3f}:")
                print("-" * 40)
                
                for certitude in certitudes:
                    success_rate = 0
                    total_time = 0
                    strips_found_sum = 0
                    
                    for test_i in range(n_tests_per_config):
                        # Créer matrice parfaite puis ajouter du bruit
                        perfect_matrix = create_perfect_strip_matrix(n_reads, n_positions, n_strips=n_strips)
                        if error_rate > 0:
                            test_matrix = add_noise_to_matrix(perfect_matrix, error_rate)
                        else:
                            test_matrix = perfect_matrix
                        
                        # Tester le preprocessing
                        start_time = time.time()
                        try:
                            inhomogeneous, strips = pre_processing(
                                test_matrix, 
                                min_col_quality=3, 
                                error_rate=error_rate,
                                certitude=certitude
                            )
                            processing_time = time.time() - start_time
                            
                            # Vérifier si on a trouvé le nombre de strips attendu
                            if len(strips) == n_strips:
                                success_rate += 1
                            
                            strips_found_sum += len(strips)
                            total_time += processing_time
                            
                        except Exception as e:
                            print(f"Erreur lors du test avec error_rate={error_rate}, certitude={certitude}: {e}")
                            continue
                    
                    avg_success_rate = success_rate / n_tests_per_config
                    avg_strips_found = strips_found_sum / n_tests_per_config
                    avg_time = total_time / n_tests_per_config
                    
                    results.append({
                        'matrix_size': f"{n_reads}x{n_positions}",
                        'n_strips': n_strips,
                        'certitude': certitude,
                        'error_rate': error_rate,
                        'success_rate': avg_success_rate,
                        'avg_strips_found': avg_strips_found,
                        'avg_time': avg_time
                    })
                    
                    print(f"  Certitude {certitude:5.3f}: "
                          f"Success {avg_success_rate*100:5.1f}%, "
                          f"Strips {avg_strips_found:4.1f}/{n_strips}, "
                          f"Temps {avg_time*1000:6.1f}ms")
            
            print("\n")
    
    return results

def plot_results(results):
    """Visualise les résultats"""
    import matplotlib.pyplot as plt
    
    # Grouper par taille de matrice, nombre de strips et certitude
    matrix_sizes = list(set([r['matrix_size'] for r in results]))
    nb_strips = list(set([r['n_strips'] for r in results]))
    certitudes = sorted(list(set([r['certitude'] for r in results])))
    
    # Graphique 1: Performance vs Error Rate pour différentes certitudes (3 strips, matrice 200x150)
    fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6), constrained_layout=True)
    fig1.suptitle('Performance vs Taux d\'Erreur (3 strips, 200x150)')
    
    colors = ['blue', 'red', 'green', 'orange', 'purple', 'brown', 'pink', 'gray']
    target_matrix = '200x150'
    target_strips = 6
    
    for j, certitude in enumerate(certitudes):
        data = [r for r in results if r['matrix_size'] == target_matrix 
               and r['n_strips'] == target_strips and r['certitude'] == certitude]
        if not data:
            continue
            
        error_rates = [d['error_rate'] for d in data]
        success_rates = [d['success_rate'] for d in data]
        avg_strips = [d['avg_strips_found'] for d in data]
        
        color = colors[j % len(colors)]
        label = f'cert={certitude:.3f}'
        
        ax1.plot(error_rates, success_rates, 'o-', 
               color=color, linewidth=2, markersize=4, label=label)
        ax2.plot(error_rates, avg_strips, 'o-', 
               color=color, linewidth=2, markersize=4, label=label)
    
    ax1.set_title('Taux de Succès vs Taux d\'Erreur')
    ax1.set_xlabel('Taux d\'Erreur')
    ax1.set_ylabel('Taux de Succès')
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(-0.1, 1.1)
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    ax2.set_title('Strips Trouvés vs Taux d\'Erreur')
    ax2.set_xlabel('Taux d\'Erreur')
    ax2.set_ylabel('Nombre de Strips')
    ax2.grid(True, alpha=0.3)
    ax2.axhline(y=target_strips, color='black', linestyle='--', alpha=0.3, label='Attendu')
    ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.savefig('/home/mafoin/stage/MatrixStriper/preprocessing_performance_certitude.png', dpi=300, bbox_inches='tight')
    
    # Graphique 2: Heatmap certitude vs error_rate pour le taux de succès
    # Calculer le layout optimal (carré ou ligne)
    if len(nb_strips) <= 3:
        nrows, ncols = 1, len(nb_strips)
        figsize = (4*len(nb_strips), 6)
    else:
        # Pour plus de 3 strips, utiliser un layout en grille
        ncols = int(np.ceil(np.sqrt(len(nb_strips))))
        nrows = int(np.ceil(len(nb_strips) / ncols))
        figsize = (4*ncols, 6*nrows)
    
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, constrained_layout=True)
    # S'assurer que axes est toujours un tableau 1D
    if len(nb_strips) == 1:
        axes = [axes]
    elif len(nb_strips) > 1 and nrows > 1:
        axes = axes.flatten()
    fig.suptitle('Heatmap: Taux de Succès (Certitude vs Error Rate) - 200x150')
    if len(nb_strips) == 1:
        axes = [axes]
    
    for i, n_strips in enumerate(nb_strips):
        # Créer une matrice pour le heatmap
        error_rates_subset = sorted([r['error_rate'] for r in results if r['error_rate'] <= 0.1])
        error_rates_subset = sorted(list(set(error_rates_subset)))
        
        heatmap_data = np.zeros((len(certitudes), len(error_rates_subset)))
        
        for j, certitude in enumerate(certitudes):
            for k, error_rate in enumerate(error_rates_subset):
                data = [r for r in results if r['matrix_size'] == target_matrix 
                       and r['n_strips'] == n_strips and r['certitude'] == certitude 
                       and abs(r['error_rate'] - error_rate) < 0.001]
                if data:
                    heatmap_data[j, k] = data[0]['success_rate']
        
        im = axes[i].imshow(heatmap_data, cmap='RdYlGn', aspect='auto', vmin=0, vmax=1)
        axes[i].set_title(f'{n_strips} strips')
        axes[i].set_xlabel('Error Rate')
        axes[i].set_ylabel('Certitude')
        axes[i].set_xticks(range(len(error_rates_subset)))
        axes[i].set_xticklabels([f'{er:.4f}' for er in error_rates_subset], rotation=45)
        axes[i].set_yticks(range(len(certitudes)))
        axes[i].set_yticklabels([f'{c:.4f}' for c in certitudes])
        
        # Ajouter les valeurs dans les cellules
        for y in range(len(certitudes)):
            for x in range(len(error_rates_subset)):
                if heatmap_data[y, x] > 0:
                    text = axes[i].text(x, y, f'{heatmap_data[y, x]:.1f}',
                                       ha="center", va="center", color="black", fontsize=8)
    

    plt.savefig('/home/mafoin/stage/MatrixStriper/heatmap_certitude_error.png', dpi=300, bbox_inches='tight')
    
    plt.close('all')  # Fermer les figures pour libérer la mémoire

if __name__ == "__main__":
    np.random.seed(42)  # Pour la reproductibilité
    
    
    # Test de performance
    results = test_preprocessing_performance()
    
    # Visualisation
    plot_results(results)

