import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import argparse
from pathlib import Path

def load_results(csv_path):
    """Load experiment results from CSV file."""
    df = pd.read_csv(csv_path)
    # Convert sous_dossier to numeric to represent number of haplotypes
    df['nb_haplotypes'] = pd.to_numeric(df['sous_dossier'])
    
    # Extract total matrix size
    df[['rows', 'cols']] = df['taille_matrice'].str.extract(r'(\d+)×(\d+)').astype(int)
    df['taille_totale'] = df['rows'] * df['cols']
    
    # Calculate derived metrics
    df['taux_reduction'] = ((df['nan_avant'] - df['nan_après']) / df['nan_avant'] * 100).round(2)
    df['densite_nan_avant'] = (df['nan_avant'] / df['taille_totale'] * 100).round(2)
    df['densite_nan_apres'] = (df['nan_après'] / df['taille_totale'] * 100).round(2)
    
    return df

def analyze_basic_stats(df):
    """Analyze basic statistics of the results."""
    print("=== RESULTS ANALYSIS ===\n")
    
    # Statistics by number of haplotypes
    print("1. STATISTICS BY NUMBER OF HAPLOTYPES:")
    by_haplotypes = df.groupby('nb_haplotypes').agg({
        'fichier': 'nunique',
        'nan_avant': 'mean',
        'nan_après': 'mean',
        'taux_reduction': 'mean'
    }).round(2)
    by_haplotypes.columns = ['nb_matrices', 'mean_nan_before', 'mean_nan_after', 'mean_reduction_rate']
    print(by_haplotypes)
    print()
    
    # Statistics by k value
    print("2. PERFORMANCE BY K VALUE:")
    by_k = df.groupby('k_value').agg({
        'nan_après': ['mean', 'std', 'min', 'max'],
        'taux_reduction': 'mean'
    }).round(2)
    by_k.columns = ['mean', 'std', 'min', 'max', 'mean_reduction_rate']
    print(by_k)
    print()
    
    # Best k value globally
    best_k = df.groupby('k_value')['nan_après'].mean().idxmin()
    best_mean = df.groupby('k_value')['nan_après'].mean().min()
    print(f"3. BEST GLOBAL K VALUE: k={best_k} (mean: {best_mean:.2f} uncertain values)")
    print()
    
    # Analyze matrices with 0 uncertain values (perfect case)
    perfect_cases = df[df['nan_après'] == 0]
    if len(perfect_cases) > 0:
        print(f"4. PERFECTLY IMPUTED MATRICES: {len(perfect_cases)} cases ({len(perfect_cases)/len(df)*100:.1f}%)")
        perfect_by_k = perfect_cases.groupby('k_value').size()
        print("   Distribution by k:")
        for k, count in perfect_by_k.items():
            print(f"   k={k}: {count} perfect matrices")
        print()
        
    # Top 10 most problematic matrices (excluding perfect cases)
    problematic = df[(df['k_value'] == 11) & (df['nan_après'] > 0)]
    if len(problematic) > 0:
        print("5. TOP 10 MATRICES WITH MOST UNCERTAIN VALUES (k=11, excluding perfect ones):")
        top_problematic = problematic.nlargest(10, 'nan_après')
        for _, row in top_problematic.iterrows():
            print(f"   {row['nb_haplotypes']} haplotypes/{row['fichier']}: {row['nan_après']} uncertain ({row['taille_matrice']})")
    print()

def analyze_correlations(df):
    """Analyze detailed correlations between variables."""
    print("5. DETAILED CORRELATION ANALYSES:")
    
    # Numeric variables for analysis
    numeric_vars = ['nb_haplotypes', 'k_value', 'taille_totale', 'rows', 'cols', 
                   'nan_avant', 'nan_après', 'taux_reduction', 'densite_nan_avant', 'densite_nan_apres']
    
    # Complete correlation matrix
    correlation_matrix = df[numeric_vars].corr()
    print("Complete correlation matrix:")
    print(correlation_matrix.round(3))
    print()
    
    # Specific correlations of interest
    print("SPECIFIC CORRELATIONS:")
    
    print("a) Matrix size vs NaN:")
    print(f"   Total size ↔ NaN before: {df['taille_totale'].corr(df['nan_avant']):.3f}")
    print(f"   Total size ↔ NaN after: {df['taille_totale'].corr(df['nan_après']):.3f}")
    print(f"   Rows ↔ NaN before: {df['rows'].corr(df['nan_avant']):.3f}")
    print(f"   Columns ↔ NaN before: {df['cols'].corr(df['nan_avant']):.3f}")
    print()
    
    print("b) K-value vs Performance:")
    print(f"   K-value ↔ NaN after: {df['k_value'].corr(df['nan_après']):.3f}")
    print(f"   K-value ↔ Reduction rate: {df['k_value'].corr(df['taux_reduction']):.3f}")
    print()
    
    print("c) Haplotypes vs Performance:")
    print(f"   Haplotypes ↔ NaN before: {df['nb_haplotypes'].corr(df['nan_avant']):.3f}")
    print(f"   Haplotypes ↔ NaN after: {df['nb_haplotypes'].corr(df['nan_après']):.3f}")
    print(f"   Haplotypes ↔ Reduction rate: {df['nb_haplotypes'].corr(df['taux_reduction']):.3f}")
    print()
    
    print("d) NaN density vs Complexity:")
    print(f"   Haplotypes ↔ NaN density before: {df['nb_haplotypes'].corr(df['densite_nan_avant']):.3f}")
    print(f"   Size ↔ NaN density before: {df['taille_totale'].corr(df['densite_nan_avant']):.3f}")
    print()

def analyze_efficiency_by_groups(df):
    """Analyze efficiency by different groupings."""
    print("6. EFFICIENCY BY GROUPS:")
    
    # Efficiency by haplotype-k combination
    print("a) Efficiency by Haplotypes and K-value:")
    pivot_efficiency = df.pivot_table(values='taux_reduction', 
                                    index='nb_haplotypes', 
                                    columns='k_value', 
                                    aggfunc='mean').round(1)
    print(pivot_efficiency)
    print()
    
    # Analysis by size categories
    df['taille_categorie'] = pd.cut(df['taille_totale'], 
                                   bins=[0, 5000, 10000, 20000, np.inf], 
                                   labels=['Small', 'Medium', 'Large', 'Very large'])
    
    print("b) Performance by Size Category:")
    by_size = df.groupby('taille_categorie').agg({
        'nan_après': 'mean',
        'taux_reduction': 'mean',
        'fichier': 'count'
    }).round(2)
    by_size.columns = ['mean_nan_after', 'mean_reduction_rate', 'nb_matrices']
    print(by_size)
    print()

def create_plots(df):
    """Create visualization plots."""
    plt.style.use('seaborn-v0_8')
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('KNN Imputation Results Analysis by Number of Haplotypes', fontsize=16, fontweight='bold')
    
    # Check available k values and choose the best
    available_k = sorted(df['k_value'].unique())
    best_k = df.groupby('k_value')['nan_après'].mean().idxmin()
    print(f"Using k={best_k} for detailed graphs")
    
    # Plot 1: Performance by k value
    ax1 = axes[0, 0]
    k_performance = df.groupby('k_value')['nan_après'].agg(['mean', 'std'])
    ax1.errorbar(k_performance.index, k_performance['mean'], 
                yerr=k_performance['std'], marker='o', capsize=5, capthick=2)
    ax1.set_xlabel('K Value')
    ax1.set_ylabel('Average Uncertain Values')
    ax1.set_title('Average Performance by K Value')
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Performance by number of haplotypes
    ax2 = axes[0, 1]
    haplotype_perf = df[df['k_value'] == best_k].groupby('nb_haplotypes')['nan_après'].agg(['mean', 'std', 'count'])
    
    # Filter haplotypes with enough data
    haplotype_perf = haplotype_perf[haplotype_perf['count'] >= 3]
    
    if not haplotype_perf.empty:
        bars = ax2.bar(haplotype_perf.index, haplotype_perf['mean'], 
                      yerr=haplotype_perf['std'], capsize=5, alpha=0.7)
        ax2.set_xlabel('Number of Haplotypes')
        ax2.set_ylabel('Average Uncertain Values')
        ax2.set_title(f'Performance by Number of Haplotypes (k={best_k})')
        ax2.grid(True, alpha=0.3, axis='y')
        
        # Add number of matrices on each bar
        for bar, count in zip(bars, haplotype_perf['count']):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height + bar.get_height()*0.01,
                    f'n={int(count)}', ha='center', va='bottom', fontsize=8)
    else:
        ax2.text(0.5, 0.5, 'Insufficient data\nfor analysis', 
                ha='center', va='center', transform=ax2.transAxes, fontsize=12)
        ax2.set_title('Performance by Haplotypes')
    
    # Plot 3: Performance heatmap
    ax3 = axes[1, 0]
    pivot_data = df.pivot_table(values='nan_après', index='nb_haplotypes', 
                               columns='k_value', aggfunc='mean')
    if not pivot_data.empty and pivot_data.max().max() > 0:
        sns.heatmap(pivot_data, annot=True, fmt='.1f', cmap='YlOrRd', ax=ax3)
        ax3.set_title('Heatmap: Average Uncertain Values')
        ax3.set_xlabel('K Value')
        ax3.set_ylabel('Number of Haplotypes')
    else:
        ax3.text(0.5, 0.5, 'Perfect results\nfor all cases', 
                ha='center', va='center', transform=ax3.transAxes, fontsize=12)
        ax3.set_title('Perfect Imputation')
    
    # Plot 4: Efficiency evolution by K
    ax4 = axes[1, 1]
    
    # Calculate average efficiency by k
    efficiency_by_k = df.groupby('k_value').agg({
        'taux_reduction': 'mean',
        'nan_après': 'mean'
    })
    
    # Double axis to show reduction rate and uncertain values
    ax4_twin = ax4.twinx()
    
    line1 = ax4.plot(efficiency_by_k.index, efficiency_by_k['taux_reduction'], 
                    'b-o', label='Reduction Rate (%)', linewidth=2)
    line2 = ax4_twin.plot(efficiency_by_k.index, efficiency_by_k['nan_après'], 
                         'r-s', label='Uncertain Values', linewidth=2)
    
    ax4.set_xlabel('K Value')
    ax4.set_ylabel('Average Reduction Rate (%)', color='b')
    ax4_twin.set_ylabel('Average Uncertain Values', color='r')
    ax4.set_title('Efficiency Evolution by K Value')
    ax4.grid(True, alpha=0.3)
    
    # Combine legends
    lines = line1 + line2
    labels = [l.get_label() for l in lines]
    ax4.legend(lines, labels, loc='center right')
    
    plt.tight_layout()
    plt.savefig('knn_analysis_plots.png', dpi=300, bbox_inches='tight')
    print("Graphs saved in 'knn_analysis_plots.png'")

def create_detailed_plots(df):
    """Create additional detailed plots."""
    available_k = sorted(df['k_value'].unique())
    best_k = df.groupby('k_value')['nan_après'].mean().idxmin()
    
    # Plot evolution by k for each number of haplotypes
    plt.figure(figsize=(12, 8))
    colors = plt.cm.Set1(np.linspace(0, 1, len(df['nb_haplotypes'].unique())))
    
    for i, haplotypes in enumerate(sorted(df['nb_haplotypes'].unique())):
        haplotype_data = df[df['nb_haplotypes'] == haplotypes]
        k_means = haplotype_data.groupby('k_value')['nan_après'].mean()
        if len(k_means) > 1:  # Only if we have multiple points
            plt.plot(k_means.index, k_means.values, marker='o', 
                    label=f'{haplotypes} haplotypes', color=colors[i])
    
    plt.xlabel('K Value')
    plt.ylabel('Average Uncertain Values')
    plt.title('Performance Evolution by K Value and Number of Haplotypes')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig('knn_evolution_by_haplotypes.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Distribution of reduction rates (only for non-perfect cases)
    df['taux_reduction'] = ((df['nan_avant'] - df['nan_après']) / df['nan_avant'] * 100)
    non_perfect = df[(df['k_value'] == best_k) & (df['nan_après'] > 0)]
    
    if len(non_perfect) > 0:
        plt.figure(figsize=(10, 6))
        plt.hist(non_perfect['taux_reduction'], bins=30, alpha=0.7, edgecolor='black')
        plt.xlabel('Uncertainty Reduction Rate (%)')
        plt.ylabel('Number of Matrices')
        plt.title(f'Distribution of Uncertainty Reduction Rates (k={best_k}, non-perfect cases)')
        plt.grid(True, alpha=0.3)
        plt.savefig('knn_reduction_distribution.png', dpi=300, bbox_inches='tight')
        plt.show()
    else:
        print("All matrices were perfectly imputed - no distribution to display")

def create_correlation_plots(df):
    """Create correlation and scatter plots."""
    plt.style.use('default')
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Correlation Analyses - KNN Imputation by Haplotypes', fontsize=16, fontweight='bold')
    
    # Automatically choose the best k value
    best_k = df.groupby('k_value')['nan_après'].mean().idxmin()
    
    # Plot 1: Correlation matrix (heatmap)
    ax1 = axes[0, 0]
    corr_vars = ['nb_haplotypes', 'k_value', 'taille_totale', 'nan_avant', 'nan_après', 'taux_reduction']
    corr_matrix = df[corr_vars].corr()
    
    # Mask upper part to avoid redundancy
    mask = np.triu(np.ones_like(corr_matrix, dtype=bool))
    sns.heatmap(corr_matrix, mask=mask, annot=True, fmt='.2f', cmap='RdBu_r', 
                center=0, ax=ax1, square=True)
    ax1.set_title('Correlation Matrix')
    
    # Plot 2: Size vs NaN before
    ax2 = axes[0, 1]
    scatter = ax2.scatter(df['taille_totale'], df['nan_avant'], 
                         c=df['nb_haplotypes'], cmap='viridis', alpha=0.6)
    ax2.set_xlabel('Total Matrix Size')
    ax2.set_ylabel('NaN Before Imputation')
    ax2.set_title('Size vs NaN Before Relationship (color = haplotypes)')
    ax2.set_xscale('log')
    plt.colorbar(scatter, ax=ax2, label='Nb Haplotypes')
    
    # Plot 3: K-value vs Efficiency by haplotypes
    ax3 = axes[1, 0]
    for haplotype in sorted(df['nb_haplotypes'].unique()):
        data = df[df['nb_haplotypes'] == haplotype]
        k_means = data.groupby('k_value')['taux_reduction'].mean()
        if len(k_means) > 1:  # Only plot if we have multiple points
            ax3.plot(k_means.index, k_means.values, 'o-', label=f'{haplotype} haplotypes', alpha=0.7)
    
    ax3.set_xlabel('K Value')
    ax3.set_ylabel('Average Reduction Rate (%)')
    ax3.set_title('Efficiency by K-value and Haplotypes')
    ax3.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: NaN density vs Performance
    ax4 = axes[1, 1]
    best_k_data = df[df['k_value'] == best_k]  # Use best k value
    
    if len(best_k_data) > 0:
        scatter = ax4.scatter(best_k_data['densite_nan_avant'], best_k_data['densite_nan_apres'], 
                             c=best_k_data['nb_haplotypes'], cmap='plasma', alpha=0.6, s=50)
        
        # Reference line y=x (no improvement)
        max_val = max(best_k_data['densite_nan_avant'].max(), best_k_data['densite_nan_apres'].max())
        ax4.plot([0, max_val], [0, max_val], 'r--', alpha=0.5, label='No improvement')
        
        ax4.set_xlabel('NaN Density Before (%)')
        ax4.set_ylabel('NaN Density After (%)')
        ax4.set_title(f'NaN Density Improvement (k={best_k})')
        ax4.legend()
        plt.colorbar(scatter, ax=ax4, label='Nb Haplotypes')
    else:
        ax4.text(0.5, 0.5, 'No data\navailable', 
                ha='center', va='center', transform=ax4.transAxes, fontsize=12)
        ax4.set_title('No Data')
    
    plt.tight_layout()
    plt.savefig('correlation_analysis.png', dpi=300, bbox_inches='tight')
    print("Correlation graphs saved in 'correlation_analysis.png'")

def create_performance_plots(df):
    """Create performance-focused plots."""
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Performance Analysis - KNN Imputation', fontsize=16, fontweight='bold')
    
    # Automatically choose the best k value
    best_k = df.groupby('k_value')['nan_après'].mean().idxmin()
    
    # Plot 1: Boxplot performance by haplotypes
    ax1 = axes[0, 0]
    best_k_data = df[df['k_value'] == best_k]
    
    if len(best_k_data) > 0:
        haplotypes_with_data = []
        box_data = []
        
        for h in sorted(best_k_data['nb_haplotypes'].unique()):
            data = best_k_data[best_k_data['nb_haplotypes'] == h]['taux_reduction']
            if len(data) >= 3:  # At least 3 points for meaningful boxplot
                haplotypes_with_data.append(h)
                box_data.append(data)
        
        if box_data:
            ax1.boxplot(box_data, labels=haplotypes_with_data)
            ax1.set_xlabel('Number of Haplotypes')
            ax1.set_ylabel('Reduction Rate (%)')
            ax1.set_title(f'Reduction Rate Distribution by Haplotypes (k={best_k})')
            ax1.grid(True, alpha=0.3)
        else:
            ax1.text(0.5, 0.5, 'Insufficient data\nfor boxplots', 
                    ha='center', va='center', transform=ax1.transAxes, fontsize=12)
            ax1.set_title('Insufficient Data')
    
    # Plot 2: Average performance by K
    ax2 = axes[0, 1]
    k_performance = df.groupby('k_value').agg({
        'nan_après': 'mean',
        'taux_reduction': 'mean'
    })
    
    ax2_twin = ax2.twinx()
    line1 = ax2.plot(k_performance.index, k_performance['nan_après'], 'b-o', label='NaN After')
    line2 = ax2_twin.plot(k_performance.index, k_performance['taux_reduction'], 'r-s', label='Reduction Rate (%)')
    
    ax2.set_xlabel('K Value')
    ax2.set_ylabel('NaN After Imputation', color='b')
    ax2_twin.set_ylabel('Reduction Rate (%)', color='r')
    ax2.set_title('Global Performance by K Value')
    
    # Combine legends
    lines = line1 + line2
    labels = [l.get_label() for l in lines]
    ax2.legend(lines, labels, loc='center right')
    
    # Plot 3: Performance heatmap by haplotype and k
    ax3 = axes[1, 0]
    pivot_perf = df.pivot_table(values='nan_après', index='nb_haplotypes', 
                               columns='k_value', aggfunc='mean')
    if not pivot_perf.empty:
        sns.heatmap(pivot_perf, annot=True, fmt='.0f', cmap='YlOrRd', ax=ax3)
        ax3.set_xlabel('K Value')
        ax3.set_ylabel('Number of Haplotypes')
        ax3.set_title('NaN After Imputation (Heatmap)')
    
    # Plot 4: Size vs efficiency scatter with regression
    ax4 = axes[1, 1]
    k_best = df.groupby(['nb_haplotypes', 'fichier'])['nan_après'].idxmin()
    best_results = df.loc[k_best]
    
    if len(best_results) > 5:  # At least 5 points for regression
        scatter = ax4.scatter(best_results['taille_totale'], best_results['taux_reduction'], 
                             c=best_results['nb_haplotypes'], cmap='viridis', alpha=0.7)
        
        # Regression line
        valid_data = best_results.dropna(subset=['taille_totale', 'taux_reduction'])
        if len(valid_data) > 2:
            z = np.polyfit(np.log10(valid_data['taille_totale']), valid_data['taux_reduction'], 1)
            p = np.poly1d(z)
            x_reg = np.logspace(np.log10(valid_data['taille_totale'].min()), 
                               np.log10(valid_data['taille_totale'].max()), 100)
            ax4.plot(x_reg, p(np.log10(x_reg)), 'r--', alpha=0.8, 
                    label=f'Regression: y={z[0]:.2f}*log(x)+{z[1]:.1f}')
        
        ax4.set_xlabel('Total Size (log scale)')
        ax4.set_ylabel('Reduction Rate (%)')
        ax4.set_title('Efficiency vs Size (best K per matrix)')
        ax4.set_xscale('log')
        ax4.legend()
        plt.colorbar(scatter, ax=ax4, label='Nb Haplotypes')
    else:
        ax4.text(0.5, 0.5, 'Insufficient data\nfor regression', 
                ha='center', va='center', transform=ax4.transAxes, fontsize=12)
        ax4.set_title('Insufficient Data')
    
    plt.tight_layout()
    plt.savefig('performance_analysis.png', dpi=300, bbox_inches='tight')
    print("Performance graphs saved in 'performance_analysis.png'")

def find_optimal_k(df):
    """Find optimal k value for each matrix type."""
    print("7. OPTIMAL K VALUE PER MATRIX:")
    
    # For each matrix, find optimal k
    optimal_k = df.loc[df.groupby(['nb_haplotypes', 'fichier'])['nan_après'].idxmin()]
    
    # Distribution of optimal k values
    k_distribution = optimal_k['k_value'].value_counts().sort_index()
    print("Distribution of optimal k values:")
    for k, count in k_distribution.items():
        print(f"   k={k}: {count} matrices ({count/len(optimal_k)*100:.1f}%)")
    print()
    
    # Average optimal K by number of haplotypes
    optimal_by_haplotypes = optimal_k.groupby('nb_haplotypes')['k_value'].agg(['mean', 'std', 'count']).round(2)
    print("Optimal K by number of haplotypes:")
    print(optimal_by_haplotypes)
    print()

def main():
    parser = argparse.ArgumentParser(description="Analyze KNN experimentation results by number of haplotypes")
    parser.add_argument("results_file", nargs='?', default="experiment_results.csv",
                       help="CSV results file (default: experiment_results.csv)")
    parser.add_argument("--plots", action="store_true", 
                       help="Generate visualization plots")
    parser.add_argument("--correlations", action="store_true",
                       help="Display detailed correlation analyses")
    
    args = parser.parse_args()
    
    # Check if file exists
    if not Path(args.results_file).exists():
        print(f"Error: File {args.results_file} does not exist.")
        return
    
    # Load data
    df = load_results(args.results_file)
    print(f"Data loaded: {len(df)} entries, {df['fichier'].nunique()} unique matrices")
    print(f"Haplotype numbers analyzed: {sorted(df['nb_haplotypes'].unique())}")
    print(f"K values tested: {sorted(df['k_value'].unique())}")
    print()
    
    # Basic analyses
    analyze_basic_stats(df)
    
    # Correlation analyses
    if args.correlations:
        analyze_correlations(df)
    
    analyze_efficiency_by_groups(df)
    find_optimal_k(df)
    
    # Generate plots
    if args.plots:
        print("Generating graphs...")
        create_correlation_plots(df)
        create_performance_plots(df)
        print("Graphs generated successfully!")
    
    # Final summary
    print("=== EXECUTIVE SUMMARY ===")
    best_k_global = df.groupby('k_value')['nan_après'].mean().idxmin()
    best_perf = df.groupby('k_value')['nan_après'].mean().min()
    perfect_rate = (df['nan_après'] == 0).mean() * 100
    
    print(f"• Best global K value: {best_k_global}")
    print(f"• Optimal average performance: {best_perf:.1f} uncertain values")
    print(f"• Perfect success rate: {perfect_rate:.1f}%")
    
    # Strongest correlation
    corr_matrix = df[['nb_haplotypes', 'k_value', 'taille_totale', 'nan_avant', 'nan_après']].corr()
    strongest_corr = corr_matrix.abs().unstack().sort_values(ascending=False)
    # Exclude perfect correlations (variable with itself)
    strongest_corr = strongest_corr[strongest_corr < 1.0].head(1)
    print(f"• Strongest correlation: {strongest_corr.index[0]} = {strongest_corr.iloc[0]:.3f}")

if __name__ == "__main__":
    main()
