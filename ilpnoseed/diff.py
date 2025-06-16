import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def load_and_process_data(file1, file2):
    """Charge et traite les données des deux fichiers CSV"""
    
    # Charger les données
    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)
    
    print(f"Colonnes disponibles dans {file1}: {list(df1.columns)}")
    print(f"Colonnes disponibles dans {file2}: {list(df2.columns)}")
    
    # Filtrer uniquement les données avec des appels ILP (ilp_calls_total > 0)
    df1_ilp = df1[df1['ilp_calls_total'] > 0].copy()
    df2_ilp = df2[df2['ilp_calls_total'] > 0].copy()
    
    print(f"Données avec ILP dans {file1}: {len(df1_ilp)} / {len(df1)} instances")
    print(f"Données avec ILP dans {file2}: {len(df2_ilp)} / {len(df2)} instances")
    
    # Trouver les instances communes (basé sur filename)
    common_files = set(df1_ilp['filename']) & set(df2_ilp['filename'])
    print(f"Instances communes avec ILP: {len(common_files)}")
    
    # Filtrer pour ne garder que les instances communes
    df1_common = df1_ilp[df1_ilp['filename'].isin(common_files)].copy()
    df2_common = df2_ilp[df2_ilp['filename'].isin(common_files)].copy()
    
    print(f"Instances finales - Sans seed: {len(df1_common)}")
    print(f"Instances finales - Avec seed: {len(df2_common)}")
    
    # Déterminer quelle colonne de temps utiliser
    time_col1 = 'execution_time' if 'execution_time' in df1_common.columns else 'total_time'
    time_col2 = 'execution_time' if 'execution_time' in df2_common.columns else 'total_time'
    
    if time_col1 not in df1_common.columns:
        raise ValueError(f"Colonne '{time_col1}' manquante dans {file1}")
    if time_col2 not in df2_common.columns:
        raise ValueError(f"Colonne '{time_col2}' manquante dans {file2}")
    
    print(f"Utilisation de '{time_col1}' comme temps d'exécution pour {file1}")
    print(f"Utilisation de '{time_col2}' comme temps d'exécution pour {file2}")
    
    return df1_common, df2_common, time_col1, time_col2

def compare_execution_times(df1, df2, time_col1, time_col2, label1="Dataset 1", label2="Dataset 2"):
    """Compare les temps d'exécution entre deux datasets"""
    
    # Créer une figure avec plusieurs sous-graphiques
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('Comparaison des temps de calcul ILP (instances communes)', fontsize=16)
    
    # 1. Comparaison des temps moyens par nombre d'haplotypes
    ax1 = axes[0, 0]
    if 'haplotype_count' in df1.columns and 'haplotype_count' in df2.columns:
        time_by_haplo1 = df1.groupby('haplotype_count')[time_col1].mean()
        time_by_haplo2 = df2.groupby('haplotype_count')[time_col2].mean()
        
        haplos = sorted(set(df1['haplotype_count'].unique()) | set(df2['haplotype_count'].unique()))
        times1 = [time_by_haplo1.get(h, 0) for h in haplos]
        times2 = [time_by_haplo2.get(h, 0) for h in haplos]
        
        x = np.arange(len(haplos))
        width = 0.35
        
        ax1.bar(x - width/2, times1, width, label=label1, alpha=0.8)
        ax1.bar(x + width/2, times2, width, label=label2, alpha=0.8)
        ax1.set_xlabel('Nombre d\'haplotypes')
        ax1.set_ylabel('Temps moyen (s)')
        ax1.set_title('Temps moyen par nombre d\'haplotypes')
        ax1.set_xticks(x)
        ax1.set_xticklabels(haplos)
        ax1.legend()
    else:
        ax1.text(0.5, 0.5, 'Colonne haplotype_count\nnon disponible', 
                ha='center', va='center', transform=ax1.transAxes)
    
    # 2. Distribution des temps d'exécution
    ax2 = axes[0, 1]
    ax2.hist(df1[time_col1], bins=20, alpha=0.7, label=label1, density=True)
    ax2.hist(df2[time_col2], bins=20, alpha=0.7, label=label2, density=True)
    ax2.set_xlabel('Temps d\'exécution (s)')
    ax2.set_ylabel('Densité')
    ax2.set_title('Distribution des temps d\'exécution')
    ax2.legend()
    
    # 3. Boxplot comparatif
    ax3 = axes[1, 0]
    data_to_plot = [df1[time_col1], df2[time_col2]]
    ax3.boxplot(data_to_plot, labels=[label1, label2])
    ax3.set_ylabel('Temps d\'exécution (s)')
    ax3.set_title('Comparaison des distributions (boxplot)')
    
    # 4. Scatter plot pour comparaison directe des instances communes
    ax4 = axes[1, 1]
    # Merger sur filename pour comparer les mêmes instances
    merged = pd.merge(df1[['filename', time_col1]], 
                     df2[['filename', time_col2]], 
                     on='filename', 
                     suffixes=('_1', '_2'))
    
    if not merged.empty:
        ax4.scatter(merged[f'{time_col1}_1'], merged[f'{time_col2}_2'], alpha=0.6)
        # Ligne y=x pour référence
        max_time = max(merged[f'{time_col1}_1'].max(), merged[f'{time_col2}_2'].max())
        ax4.plot([0, max_time], [0, max_time], 'r--', alpha=0.8)
        ax4.set_xlabel(f'Temps {label1} (s)')
        ax4.set_ylabel(f'Temps {label2} (s)')
        ax4.set_title('Comparaison directe (instances communes)')
        
        # Ajouter corrélation
        correlation = merged[f'{time_col1}_1'].corr(merged[f'{time_col2}_2'])
        ax4.text(0.05, 0.95, f'Corrélation: {correlation:.3f}', 
                transform=ax4.transAxes, bbox=dict(boxstyle="round", facecolor='wheat'))
    else:
        ax4.text(0.5, 0.5, 'Pas d\'instances\ncommunes trouvées', 
                ha='center', va='center', transform=ax4.transAxes)
    
    plt.tight_layout()
    return fig

def generate_statistics(df1, df2, time_col1, time_col2, label1="Dataset 1", label2="Dataset 2"):
    """Génère des statistiques comparatives"""
    
    stats = {
        label1: {
            'count': len(df1),
            'mean_time': df1[time_col1].mean(),
            'median_time': df1[time_col1].median(),
            'std_time': df1[time_col1].std(),
            'min_time': df1[time_col1].min(),
            'max_time': df1[time_col1].max(),
            'q25': df1[time_col1].quantile(0.25),
            'q75': df1[time_col1].quantile(0.75)
        },
        label2: {
            'count': len(df2),
            'mean_time': df2[time_col2].mean(),
            'median_time': df2[time_col2].median(),
            'std_time': df2[time_col2].std(),
            'min_time': df2[time_col2].min(),
            'max_time': df2[time_col2].max(),
            'q25': df2[time_col2].quantile(0.25),
            'q75': df2[time_col2].quantile(0.75)
        }
    }
    
    return pd.DataFrame(stats).round(4)

def main():
    """Fonction principale"""
    
    # Chemins vers vos fichiers CSV
    file1 = "exp_no_seed.csv"  # Remplacez par le chemin de votre premier fichier
    file2 = "experiment_results_20250613_190853.csv"  # Remplacez par le chemin de votre second fichier
    
    try:
        # Charger les données
        print("Chargement des données...")
        df1_ilp, df2_ilp, time_col1, time_col2 = load_and_process_data(file1, file2)
        
        print(f"Sans seed: {len(df1_ilp)} instances")
        print(f"Avec seed: {len(df2_ilp)} instances")
        
        # Générer les statistiques
        print("\n=== STATISTIQUES COMPARATIVES ===")
        stats_df = generate_statistics(df1_ilp, df2_ilp, time_col1, time_col2, "Sans seed", "Avec seed")
        print(stats_df)
        
        # Créer les graphiques de comparaison
        print("\nGénération des graphiques...")
        fig = compare_execution_times(df1_ilp, df2_ilp, time_col1, time_col2, "Sans seed", "Avec seed")
        
        # Sauvegarder les résultats
        fig.savefig('comparison_ilp_times.png', dpi=300, bbox_inches='tight')
        stats_df.to_csv('comparison_statistics.csv')
        
        print("Comparaison terminée!")
        print("- Graphiques sauvegardés: comparison_ilp_times.png")
        print("- Statistiques sauvegardées: comparison_statistics.csv")
        
        # Afficher les graphiques
        plt.show()
        
        # Analyse rapide
        print("\n=== ANALYSE RAPIDE ===")
        ratio_mean = stats_df.loc['mean_time', 'Avec seed'] / stats_df.loc['mean_time', 'Sans seed']
        print(f"Ratio temps moyen (Avec seed / Sans seed): {ratio_mean:.2f}")
        
        if ratio_mean > 1:
            print(f"ILP avec seed est {ratio_mean:.2f}x plus lent en moyenne")
        else:
            print(f"ILP avec seed est {1/ratio_mean:.2f}x plus rapide en moyenne")
        
    except FileNotFoundError as e:
        print(f"Erreur: Fichier non trouvé - {e}")
        print("Assurez-vous que les chemins vers vos fichiers CSV sont corrects")
    except Exception as e:
        print(f"Erreur: {e}")

if __name__ == "__main__":
    main()