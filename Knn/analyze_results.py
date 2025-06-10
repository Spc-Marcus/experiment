import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import argparse
from pathlib import Path

def load_results(csv_path):
    """Load experiment results from CSV file."""
    df = pd.read_csv(csv_path)
    # Convertir sous_dossier en numérique pour représenter le nombre d'haplotypes
    df['nb_haplotypes'] = pd.to_numeric(df['sous_dossier'])
    
    # Extraire la taille totale de la matrice
    df[['rows', 'cols']] = df['taille_matrice'].str.extract(r'(\d+)×(\d+)').astype(int)
    df['taille_totale'] = df['rows'] * df['cols']
    
    # Calculer des métriques dérivées
    df['taux_reduction'] = ((df['nan_avant'] - df['nan_après']) / df['nan_avant'] * 100).round(2)
    df['densite_nan_avant'] = (df['nan_avant'] / df['taille_totale'] * 100).round(2)
    df['densite_nan_apres'] = (df['nan_après'] / df['taille_totale'] * 100).round(2)
    
    return df

def analyze_basic_stats(df):
    """Analyze basic statistics of the results."""
    print("=== ANALYSE DES RÉSULTATS ===\n")
    
    # Statistics par nombre d'haplotypes
    print("1. STATISTIQUES PAR NOMBRE D'HAPLOTYPES:")
    by_haplotypes = df.groupby('nb_haplotypes').agg({
        'fichier': 'nunique',
        'nan_avant': 'mean',
        'nan_après': 'mean',
        'taux_reduction': 'mean'
    }).round(2)
    by_haplotypes.columns = ['nb_matrices', 'nan_avant_moyen', 'nan_après_moyen', 'taux_reduction_moyen']
    print(by_haplotypes)
    print()
    
    # Statistics par valeur k
    print("2. PERFORMANCE PAR VALEUR K:")
    by_k = df.groupby('k_value').agg({
        'nan_après': ['mean', 'std', 'min', 'max'],
        'taux_reduction': 'mean'
    }).round(2)
    by_k.columns = ['moyenne', 'std', 'min', 'max', 'taux_reduction_moyen']
    print(by_k)
    print()
    
    # Meilleure valeur k globalement
    best_k = df.groupby('k_value')['nan_après'].mean().idxmin()
    best_mean = df.groupby('k_value')['nan_après'].mean().min()
    print(f"3. MEILLEURE VALEUR K GLOBALE: k={best_k} (moyenne: {best_mean:.2f} valeurs incertaines)")
    print()
    
    # Analyser les matrices avec 0 valeurs incertaines (cas parfait)
    perfect_cases = df[df['nan_après'] == 0]
    if len(perfect_cases) > 0:
        print(f"4. MATRICES PARFAITEMENT IMPUTÉES: {len(perfect_cases)} cas ({len(perfect_cases)/len(df)*100:.1f}%)")
        perfect_by_k = perfect_cases.groupby('k_value').size()
        print("   Distribution par k:")
        for k, count in perfect_by_k.items():
            print(f"   k={k}: {count} matrices parfaites")
        print()
        
    # Top 10 des matrices les plus problématiques (excluant les cas parfaits)
    problematic = df[(df['k_value'] == 11) & (df['nan_après'] > 0)]
    if len(problematic) > 0:
        print("5. TOP 10 MATRICES AVEC PLUS DE VALEURS INCERTAINES (k=11, excluant les parfaites):")
        top_problematic = problematic.nlargest(10, 'nan_après')
        for _, row in top_problematic.iterrows():
            print(f"   {row['nb_haplotypes']} haplotypes/{row['fichier']}: {row['nan_après']} incertaines ({row['taille_matrice']})")
    print()

def analyze_correlations(df):
    """Analyze detailed correlations between variables."""
    print("5. ANALYSES DE CORRÉLATION DÉTAILLÉES:")
    
    # Variables numériques pour l'analyse
    numeric_vars = ['nb_haplotypes', 'k_value', 'taille_totale', 'rows', 'cols', 
                   'nan_avant', 'nan_après', 'taux_reduction', 'densite_nan_avant', 'densite_nan_apres']
    
    # Matrice de corrélation complète
    correlation_matrix = df[numeric_vars].corr()
    print("Matrice de corrélation complète:")
    print(correlation_matrix.round(3))
    print()
    
    # Corrélations spécifiques d'intérêt
    print("CORRÉLATIONS SPÉCIFIQUES:")
    
    print("a) Taille de matrice vs NaN:")
    print(f"   Taille totale ↔ NaN avant: {df['taille_totale'].corr(df['nan_avant']):.3f}")
    print(f"   Taille totale ↔ NaN après: {df['taille_totale'].corr(df['nan_après']):.3f}")
    print(f"   Lignes ↔ NaN avant: {df['rows'].corr(df['nan_avant']):.3f}")
    print(f"   Colonnes ↔ NaN avant: {df['cols'].corr(df['nan_avant']):.3f}")
    print()
    
    print("b) K-value vs Performance:")
    print(f"   K-value ↔ NaN après: {df['k_value'].corr(df['nan_après']):.3f}")
    print(f"   K-value ↔ Taux réduction: {df['k_value'].corr(df['taux_reduction']):.3f}")
    print()
    
    print("c) Haplotypes vs Performance:")
    print(f"   Haplotypes ↔ NaN avant: {df['nb_haplotypes'].corr(df['nan_avant']):.3f}")
    print(f"   Haplotypes ↔ NaN après: {df['nb_haplotypes'].corr(df['nan_après']):.3f}")
    print(f"   Haplotypes ↔ Taux réduction: {df['nb_haplotypes'].corr(df['taux_reduction']):.3f}")
    print()
    
    print("d) Densité NaN vs Complexité:")
    print(f"   Haplotypes ↔ Densité NaN avant: {df['nb_haplotypes'].corr(df['densite_nan_avant']):.3f}")
    print(f"   Taille ↔ Densité NaN avant: {df['taille_totale'].corr(df['densite_nan_avant']):.3f}")
    print()

def analyze_efficiency_by_groups(df):
    """Analyze efficiency by different groupings."""
    print("6. EFFICACITÉ PAR GROUPES:")
    
    # Efficacité par combinaison haplotypes-k
    print("a) Efficacité par Haplotypes et K-value:")
    pivot_efficiency = df.pivot_table(values='taux_reduction', 
                                    index='nb_haplotypes', 
                                    columns='k_value', 
                                    aggfunc='mean').round(1)
    print(pivot_efficiency)
    print()
    
    # Analyse par catégories de taille
    df['taille_categorie'] = pd.cut(df['taille_totale'], 
                                   bins=[0, 5000, 10000, 20000, np.inf], 
                                   labels=['Petite', 'Moyenne', 'Grande', 'Très grande'])
    
    print("b) Performance par Catégorie de Taille:")
    by_size = df.groupby('taille_categorie').agg({
        'nan_après': 'mean',
        'taux_reduction': 'mean',
        'fichier': 'count'
    }).round(2)
    by_size.columns = ['nan_après_moyen', 'taux_reduction_moyen', 'nb_matrices']
    print(by_size)
    print()

def create_plots(df):
    """Create visualization plots."""
    plt.style.use('seaborn-v0_8')
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('Analyse des Résultats d\'Imputation KNN par Nombre d\'Haplotypes', fontsize=16, fontweight='bold')
    
    # Vérifier les valeurs k disponibles et choisir la meilleure
    available_k = sorted(df['k_value'].unique())
    best_k = df.groupby('k_value')['nan_après'].mean().idxmin()
    print(f"Utilisation de k={best_k} pour les graphiques détaillés")
    
    # Plot 1: Performance par valeur k
    ax1 = axes[0, 0]
    k_performance = df.groupby('k_value')['nan_après'].agg(['mean', 'std'])
    ax1.errorbar(k_performance.index, k_performance['mean'], 
                yerr=k_performance['std'], marker='o', capsize=5, capthick=2)
    ax1.set_xlabel('Valeur K')
    ax1.set_ylabel('Valeurs Incertaines Moyennes')
    ax1.set_title('Performance Moyenne par Valeur K')
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Performance par nombre d'haplotypes (remplace le graphique vide)
    ax2 = axes[0, 1]
    haplotype_perf = df[df['k_value'] == best_k].groupby('nb_haplotypes')['nan_après'].agg(['mean', 'std', 'count'])
    
    # Filtrer les haplotypes avec assez de données
    haplotype_perf = haplotype_perf[haplotype_perf['count'] >= 3]
    
    if not haplotype_perf.empty:
        bars = ax2.bar(haplotype_perf.index, haplotype_perf['mean'], 
                      yerr=haplotype_perf['std'], capsize=5, alpha=0.7)
        ax2.set_xlabel('Nombre d\'Haplotypes')
        ax2.set_ylabel('Valeurs Incertaines Moyennes')
        ax2.set_title(f'Performance par Nombre d\'Haplotypes (k={best_k})')
        ax2.grid(True, alpha=0.3, axis='y')
        
        # Ajouter le nombre de matrices sur chaque barre
        for bar, count in zip(bars, haplotype_perf['count']):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height + bar.get_height()*0.01,
                    f'n={int(count)}', ha='center', va='bottom', fontsize=8)
    else:
        ax2.text(0.5, 0.5, 'Données insuffisantes\npour l\'analyse', 
                ha='center', va='center', transform=ax2.transAxes, fontsize=12)
        ax2.set_title('Performance par Haplotypes')
    
    # Plot 3: Heatmap de performance
    ax3 = axes[1, 0]
    pivot_data = df.pivot_table(values='nan_après', index='nb_haplotypes', 
                               columns='k_value', aggfunc='mean')
    if not pivot_data.empty and pivot_data.max().max() > 0:
        sns.heatmap(pivot_data, annot=True, fmt='.1f', cmap='YlOrRd', ax=ax3)
        ax3.set_title('Heatmap: Valeurs Incertaines Moyennes')
        ax3.set_xlabel('Valeur K')
        ax3.set_ylabel('Nombre d\'Haplotypes')
    else:
        ax3.text(0.5, 0.5, 'Résultats parfaits\npour tous les cas', 
                ha='center', va='center', transform=ax3.transAxes, fontsize=12)
        ax3.set_title('Imputation Parfaite')
    
    # Plot 4: Évolution de l'efficacité par K (remplace le graphique vide)
    ax4 = axes[1, 1]
    
    # Calculer l'efficacité moyenne par k
    efficiency_by_k = df.groupby('k_value').agg({
        'taux_reduction': 'mean',
        'nan_après': 'mean'
    })
    
    # Double axe pour montrer taux de réduction et valeurs incertaines
    ax4_twin = ax4.twinx()
    
    line1 = ax4.plot(efficiency_by_k.index, efficiency_by_k['taux_reduction'], 
                    'b-o', label='Taux Réduction (%)', linewidth=2)
    line2 = ax4_twin.plot(efficiency_by_k.index, efficiency_by_k['nan_après'], 
                         'r-s', label='Valeurs Incertaines', linewidth=2)
    
    ax4.set_xlabel('Valeur K')
    ax4.set_ylabel('Taux de Réduction Moyen (%)', color='b')
    ax4_twin.set_ylabel('Valeurs Incertaines Moyennes', color='r')
    ax4.set_title('Évolution de l\'Efficacité par Valeur K')
    ax4.grid(True, alpha=0.3)
    
    # Combiner les légendes
    lines = line1 + line2
    labels = [l.get_label() for l in lines]
    ax4.legend(lines, labels, loc='center right')
    
    plt.tight_layout()
    plt.savefig('knn_analysis_plots.png', dpi=300, bbox_inches='tight')
    print("Graphiques sauvegardés dans 'knn_analysis_plots.png'")

def create_detailed_plots(df):
    """Create additional detailed plots."""
    available_k = sorted(df['k_value'].unique())
    best_k = df.groupby('k_value')['nan_après'].mean().idxmin()
    
    # Plot de l'évolution par k pour chaque nombre d'haplotypes
    plt.figure(figsize=(12, 8))
    colors = plt.cm.Set1(np.linspace(0, 1, len(df['nb_haplotypes'].unique())))
    
    for i, haplotypes in enumerate(sorted(df['nb_haplotypes'].unique())):
        haplotype_data = df[df['nb_haplotypes'] == haplotypes]
        k_means = haplotype_data.groupby('k_value')['nan_après'].mean()
        if len(k_means) > 1:  # Seulement si on a plusieurs points
            plt.plot(k_means.index, k_means.values, marker='o', 
                    label=f'{haplotypes} haplotypes', color=colors[i])
    
    plt.xlabel('Valeur K')
    plt.ylabel('Valeurs Incertaines Moyennes')
    plt.title('Évolution de la Performance par Valeur K et Nombre d\'Haplotypes')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig('knn_evolution_by_haplotypes.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Distribution des taux de réduction (seulement pour les cas non-parfaits)
    df['taux_reduction'] = ((df['nan_avant'] - df['nan_après']) / df['nan_avant'] * 100)
    non_perfect = df[(df['k_value'] == best_k) & (df['nan_après'] > 0)]
    
    if len(non_perfect) > 0:
        plt.figure(figsize=(10, 6))
        plt.hist(non_perfect['taux_reduction'], bins=30, alpha=0.7, edgecolor='black')
        plt.xlabel('Taux de Réduction d\'Incertitude (%)')
        plt.ylabel('Nombre de Matrices')
        plt.title(f'Distribution des Taux de Réduction d\'Incertitude (k={best_k}, cas non-parfaits)')
        plt.grid(True, alpha=0.3)
        plt.savefig('knn_reduction_distribution.png', dpi=300, bbox_inches='tight')
        plt.show()
    else:
        print("Toutes les matrices ont été parfaitement imputées - pas de distribution à afficher")

def create_correlation_plots(df):
    """Create correlation and scatter plots."""
    plt.style.use('default')
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Analyses de Corrélation - Imputation KNN par Haplotypes', fontsize=16, fontweight='bold')
    
    # Choisir automatiquement la meilleure valeur de k
    best_k = df.groupby('k_value')['nan_après'].mean().idxmin()
    
    # Plot 1: Matrice de corrélation (heatmap)
    ax1 = axes[0, 0]
    corr_vars = ['nb_haplotypes', 'k_value', 'taille_totale', 'nan_avant', 'nan_après', 'taux_reduction']
    corr_matrix = df[corr_vars].corr()
    
    # Masquer la partie supérieure pour éviter la redondance
    mask = np.triu(np.ones_like(corr_matrix, dtype=bool))
    sns.heatmap(corr_matrix, mask=mask, annot=True, fmt='.2f', cmap='RdBu_r', 
                center=0, ax=ax1, square=True)
    ax1.set_title('Matrice de Corrélation')
    
    # Plot 2: Taille vs NaN avant
    ax2 = axes[0, 1]
    scatter = ax2.scatter(df['taille_totale'], df['nan_avant'], 
                         c=df['nb_haplotypes'], cmap='viridis', alpha=0.6)
    ax2.set_xlabel('Taille Totale de la Matrice')
    ax2.set_ylabel('NaN Avant Imputation')
    ax2.set_title('Relation Taille vs NaN Avant (couleur = haplotypes)')
    ax2.set_xscale('log')
    plt.colorbar(scatter, ax=ax2, label='Nb Haplotypes')
    
    # Plot 3: K-value vs Efficacité par haplotypes
    ax3 = axes[1, 0]
    for haplotype in sorted(df['nb_haplotypes'].unique()):
        data = df[df['nb_haplotypes'] == haplotype]
        k_means = data.groupby('k_value')['taux_reduction'].mean()
        if len(k_means) > 1:  # Seulement tracer si on a plusieurs points
            ax3.plot(k_means.index, k_means.values, 'o-', label=f'{haplotype} haplotypes', alpha=0.7)
    
    ax3.set_xlabel('Valeur K')
    ax3.set_ylabel('Taux de Réduction Moyen (%)')
    ax3.set_title('Efficacité par K-value et Haplotypes')
    ax3.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Densité NaN vs Performance
    ax4 = axes[1, 1]
    best_k_data = df[df['k_value'] == best_k]  # Utiliser la meilleure valeur k
    
    if len(best_k_data) > 0:
        scatter = ax4.scatter(best_k_data['densite_nan_avant'], best_k_data['densite_nan_apres'], 
                             c=best_k_data['nb_haplotypes'], cmap='plasma', alpha=0.6, s=50)
        
        # Ligne de référence y=x (pas d'amélioration)
        max_val = max(best_k_data['densite_nan_avant'].max(), best_k_data['densite_nan_apres'].max())
        ax4.plot([0, max_val], [0, max_val], 'r--', alpha=0.5, label='Pas d\'amélioration')
        
        ax4.set_xlabel('Densité NaN Avant (%)')
        ax4.set_ylabel('Densité NaN Après (%)')
        ax4.set_title(f'Amélioration de la Densité NaN (k={best_k})')
        ax4.legend()
        plt.colorbar(scatter, ax=ax4, label='Nb Haplotypes')
    else:
        ax4.text(0.5, 0.5, 'Pas de données\ndisponibles', 
                ha='center', va='center', transform=ax4.transAxes, fontsize=12)
        ax4.set_title('Pas de Données')
    
    plt.tight_layout()
    plt.savefig('correlation_analysis.png', dpi=300, bbox_inches='tight')
    print("Graphiques de corrélation sauvegardés dans 'correlation_analysis.png'")

def create_performance_plots(df):
    """Create performance-focused plots."""
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Analyse de Performance - Imputation KNN', fontsize=16, fontweight='bold')
    
    # Choisir automatiquement la meilleure valeur de k
    best_k = df.groupby('k_value')['nan_après'].mean().idxmin()
    
    # Plot 1: Boxplot performance par haplotypes
    ax1 = axes[0, 0]
    best_k_data = df[df['k_value'] == best_k]
    
    if len(best_k_data) > 0:
        haplotypes_with_data = []
        box_data = []
        
        for h in sorted(best_k_data['nb_haplotypes'].unique()):
            data = best_k_data[best_k_data['nb_haplotypes'] == h]['taux_reduction']
            if len(data) >= 3:  # Au moins 3 points pour un boxplot significatif
                haplotypes_with_data.append(h)
                box_data.append(data)
        
        if box_data:
            ax1.boxplot(box_data, labels=haplotypes_with_data)
            ax1.set_xlabel('Nombre d\'Haplotypes')
            ax1.set_ylabel('Taux de Réduction (%)')
            ax1.set_title(f'Distribution du Taux de Réduction par Haplotypes (k={best_k})')
            ax1.grid(True, alpha=0.3)
        else:
            ax1.text(0.5, 0.5, 'Données insuffisantes\npour les boxplots', 
                    ha='center', va='center', transform=ax1.transAxes, fontsize=12)
            ax1.set_title('Données Insuffisantes')
    
    # Plot 2: Performance moyenne par K
    ax2 = axes[0, 1]
    k_performance = df.groupby('k_value').agg({
        'nan_après': 'mean',
        'taux_reduction': 'mean'
    })
    
    ax2_twin = ax2.twinx()
    line1 = ax2.plot(k_performance.index, k_performance['nan_après'], 'b-o', label='NaN Après')
    line2 = ax2_twin.plot(k_performance.index, k_performance['taux_reduction'], 'r-s', label='Taux Réduction (%)')
    
    ax2.set_xlabel('Valeur K')
    ax2.set_ylabel('NaN Après Imputation', color='b')
    ax2_twin.set_ylabel('Taux de Réduction (%)', color='r')
    ax2.set_title('Performance Globale par Valeur K')
    
    # Combiner les légendes
    lines = line1 + line2
    labels = [l.get_label() for l in lines]
    ax2.legend(lines, labels, loc='center right')
    
    # Plot 3: Heatmap performance par haplotype et k
    ax3 = axes[1, 0]
    pivot_perf = df.pivot_table(values='nan_après', index='nb_haplotypes', 
                               columns='k_value', aggfunc='mean')
    if not pivot_perf.empty:
        sns.heatmap(pivot_perf, annot=True, fmt='.0f', cmap='YlOrRd', ax=ax3)
        ax3.set_xlabel('Valeur K')
        ax3.set_ylabel('Nombre d\'Haplotypes')
        ax3.set_title('NaN Après Imputation (Heatmap)')
    
    # Plot 4: Scatter taille vs efficacité avec régression
    ax4 = axes[1, 1]
    k_best = df.groupby(['nb_haplotypes', 'fichier'])['nan_après'].idxmin()
    best_results = df.loc[k_best]
    
    if len(best_results) > 5:  # Au moins 5 points pour une régression
        scatter = ax4.scatter(best_results['taille_totale'], best_results['taux_reduction'], 
                             c=best_results['nb_haplotypes'], cmap='viridis', alpha=0.7)
        
        # Ligne de régression
        valid_data = best_results.dropna(subset=['taille_totale', 'taux_reduction'])
        if len(valid_data) > 2:
            z = np.polyfit(np.log10(valid_data['taille_totale']), valid_data['taux_reduction'], 1)
            p = np.poly1d(z)
            x_reg = np.logspace(np.log10(valid_data['taille_totale'].min()), 
                               np.log10(valid_data['taille_totale'].max()), 100)
            ax4.plot(x_reg, p(np.log10(x_reg)), 'r--', alpha=0.8, 
                    label=f'Régression: y={z[0]:.2f}*log(x)+{z[1]:.1f}')
        
        ax4.set_xlabel('Taille Totale (log scale)')
        ax4.set_ylabel('Taux de Réduction (%)')
        ax4.set_title('Efficacité vs Taille (meilleur K par matrice)')
        ax4.set_xscale('log')
        ax4.legend()
        plt.colorbar(scatter, ax=ax4, label='Nb Haplotypes')
    else:
        ax4.text(0.5, 0.5, 'Données insuffisantes\npour la régression', 
                ha='center', va='center', transform=ax4.transAxes, fontsize=12)
        ax4.set_title('Données Insuffisantes')
    
    plt.tight_layout()
    plt.savefig('performance_analysis.png', dpi=300, bbox_inches='tight')
    print("Graphiques de performance sauvegardés dans 'performance_analysis.png'")

def find_optimal_k(df):
    """Find optimal k value for each matrix type."""
    print("7. VALEUR K OPTIMALE PAR MATRICE:")
    
    # Pour chaque matrice, trouver le k optimal
    optimal_k = df.loc[df.groupby(['nb_haplotypes', 'fichier'])['nan_après'].idxmin()]
    
    # Distribution des k optimaux
    k_distribution = optimal_k['k_value'].value_counts().sort_index()
    print("Distribution des k optimaux:")
    for k, count in k_distribution.items():
        print(f"   k={k}: {count} matrices ({count/len(optimal_k)*100:.1f}%)")
    print()
    
    # K optimal moyen par nombre d'haplotypes
    optimal_by_haplotypes = optimal_k.groupby('nb_haplotypes')['k_value'].agg(['mean', 'std', 'count']).round(2)
    print("K optimal par nombre d'haplotypes:")
    print(optimal_by_haplotypes)
    print()

def main():
    parser = argparse.ArgumentParser(description="Analyser les résultats d'expérimentation KNN par nombre d'haplotypes")
    parser.add_argument("results_file", nargs='?', default="experiment_results.csv",
                       help="Fichier CSV des résultats (défaut: experiment_results.csv)")
    parser.add_argument("--plots", action="store_true", 
                       help="Générer des graphiques de visualisation")
    parser.add_argument("--correlations", action="store_true",
                       help="Afficher les analyses de corrélation détaillées")
    
    args = parser.parse_args()
    
    # Vérifier que le fichier existe
    if not Path(args.results_file).exists():
        print(f"Erreur: Le fichier {args.results_file} n'existe pas.")
        return
    
    # Charger les données
    df = load_results(args.results_file)
    print(f"Données chargées: {len(df)} entrées, {df['fichier'].nunique()} matrices uniques")
    print(f"Nombres d'haplotypes analysés: {sorted(df['nb_haplotypes'].unique())}")
    print(f"Valeurs K testées: {sorted(df['k_value'].unique())}")
    print()
    
    # Analyses de base
    analyze_basic_stats(df)
    
    # Analyses de corrélation
    if args.correlations:
        analyze_correlations(df)
    
    analyze_efficiency_by_groups(df)
    find_optimal_k(df)
    
    # Génération des graphiques
    if args.plots:
        print("Génération des graphiques...")
        create_correlation_plots(df)
        create_performance_plots(df)
        print("Graphiques générés avec succès!")
    
    # Résumé final
    print("=== RÉSUMÉ EXÉCUTIF ===")
    best_k_global = df.groupby('k_value')['nan_après'].mean().idxmin()
    best_perf = df.groupby('k_value')['nan_après'].mean().min()
    perfect_rate = (df['nan_après'] == 0).mean() * 100
    
    print(f"• Meilleure valeur K globale: {best_k_global}")
    print(f"• Performance moyenne optimale: {best_perf:.1f} valeurs incertaines")
    print(f"• Taux de succès parfait: {perfect_rate:.1f}%")
    
    # Corrélation la plus forte
    corr_matrix = df[['nb_haplotypes', 'k_value', 'taille_totale', 'nan_avant', 'nan_après']].corr()
    strongest_corr = corr_matrix.abs().unstack().sort_values(ascending=False)
    # Exclure les corrélations parfaites (variable avec elle-même)
    strongest_corr = strongest_corr[strongest_corr < 1.0].head(1)
    print(f"• Corrélation la plus forte: {strongest_corr.index[0]} = {strongest_corr.iloc[0]:.3f}")

if __name__ == "__main__":
    main()
