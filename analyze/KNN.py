import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import argparse
from pathlib import Path
import sys
from io import StringIO

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

def analyze_basic_stats(df, output_file=None):
    """Analyze basic statistics of the results."""
    output = StringIO() if output_file else sys.stdout
    
    print("=== KNN IMPUTATION RESULTS ===\n", file=output)
    
    # Best k value globally
    best_k = df.groupby('k_value')['nan_après'].mean().idxmin()
    best_mean = df.groupby('k_value')['nan_après'].mean().min()
    print(f"BEST GLOBAL K VALUE: k={best_k} (average: {best_mean:.1f} remaining NaN)", file=output)
    print(file=output)
    
    # Performance summary by k value
    print("PERFORMANCE BY K VALUE:", file=output)
    by_k = df.groupby('k_value').agg({
        'nan_après': 'mean',
        'taux_reduction': 'mean'
    }).round(1)
    by_k.columns = ['avg_remaining_nan', 'avg_reduction_rate_%']
    print(by_k, file=output)
    print(file=output)
    
    # Perfect cases
    perfect_cases = df[df['nan_après'] == 0]
    perfect_rate = len(perfect_cases) / len(df) * 100
    print(f"PERFECT IMPUTATION RATE: {perfect_rate:.1f}% ({len(perfect_cases)}/{len(df)} cases)", file=output)
    print(file=output)
    
    if output_file:
        return output.getvalue()

def analyze_correlations(df, output_file=None):
    """Analyze key correlations."""
    output = StringIO() if output_file else sys.stdout
    
    print("KEY CORRELATIONS:", file=output)
    print(f"• Matrix size ↔ Initial NaN: {df['taille_totale'].corr(df['nan_avant']):.3f}", file=output)
    print(f"• Haplotypes ↔ Reduction rate: {df['nb_haplotypes'].corr(df['taux_reduction']):.3f}", file=output)
    print(f"• K-value ↔ Remaining NaN: {df['k_value'].corr(df['nan_après']):.3f}", file=output)
    print(file=output)
    
    if output_file:
        return output.getvalue()

def create_correlation_plots(df, output_dir):
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
    output_path = output_dir / 'correlation_analysis.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Correlation graphs saved in '{output_path}'")

def create_performance_plots(df, output_dir):
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
            ax1.boxplot(box_data, tick_labels=haplotypes_with_data)
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
    output_path = output_dir / 'performance_analysis.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Performance graphs saved in '{output_path}'")

def analyze_efficiency_by_groups(df, output_file=None):
    """Analyze efficiency by haplotype groups."""
    output = StringIO() if output_file else sys.stdout
    
    print("PERFORMANCE BY HAPLOTYPE COUNT:", file=output)
    by_haplotypes = df.groupby('nb_haplotypes').agg({
        'fichier': 'nunique',
        'taux_reduction': 'mean',
        'nan_après': 'mean'
    }).round(1)
    by_haplotypes.columns = ['matrices_count', 'avg_reduction_%', 'avg_remaining_nan']
    print(by_haplotypes, file=output)
    print(file=output)
    
    if output_file:
        return output.getvalue()

def find_optimal_k(df, output_file=None):
    """Find optimal k distribution."""
    output = StringIO() if output_file else sys.stdout
    
    print("OPTIMAL K DISTRIBUTION:", file=output)
    
    # For each matrix, find optimal k
    optimal_k = df.loc[df.groupby(['nb_haplotypes', 'fichier'])['nan_après'].idxmin()]
    k_distribution = optimal_k['k_value'].value_counts().sort_index()
    
    for k, count in k_distribution.items():
        percentage = count/len(optimal_k)*100
        print(f"k={k}: {count} matrices ({percentage:.1f}%)", file=output)
    print(file=output)
    
    if output_file:
        return output.getvalue()

def generate_summary(df, output_dir):
    """Generate a concise summary text file."""
    summary_content = []
    
    # Header
    summary_content.append("=" * 60)
    summary_content.append("KNN IMPUTATION - ANALYSIS SUMMARY")
    summary_content.append("=" * 60)
    summary_content.append(f"Dataset: {len(df)} experiments on {df['fichier'].nunique()} matrices")
    summary_content.append(f"Haplotypes tested: {sorted(df['nb_haplotypes'].unique())}")
    summary_content.append(f"K values tested: {sorted(df['k_value'].unique())}")
    summary_content.append("")
    
    # Add analysis sections
    summary_content.append(analyze_basic_stats(df, output_file=True))
    summary_content.append(analyze_correlations(df, output_file=True))
    summary_content.append(analyze_efficiency_by_groups(df, output_file=True))
    summary_content.append(find_optimal_k(df, output_file=True))
    
    # Executive summary
    best_k_global = df.groupby('k_value')['nan_après'].mean().idxmin()
    best_perf = df.groupby('k_value')['nan_après'].mean().min()
    perfect_rate = (df['nan_après'] == 0).mean() * 100
    avg_reduction = df['taux_reduction'].mean()
    
    summary_content.append("=" * 60)
    summary_content.append("EXECUTIVE SUMMARY")
    summary_content.append("=" * 60)
    summary_content.append(f"• Best K value globally: {best_k_global}")
    summary_content.append(f"• Average remaining NaN with best K: {best_perf:.1f}")
    summary_content.append(f"• Perfect imputation success rate: {perfect_rate:.1f}%")
    summary_content.append(f"• Overall average reduction rate: {avg_reduction:.1f}%")
    summary_content.append("=" * 60)
    
    # Save summary
    summary_path = output_dir / 'analysis_summary.txt'
    with open(summary_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(summary_content))
    
    print(f"Analysis summary saved in '{summary_path}'")
    return summary_path

def main():
    parser = argparse.ArgumentParser(description="Analyze KNN experimentation results by number of haplotypes")
    parser.add_argument("results_file", nargs='?', default="experiment_results.csv",
                       help="CSV results file (default: experiment_results.csv)")
    parser.add_argument("--plots", action="store_true", 
                       help="Generate visualization plots")
    parser.add_argument("--output", "-o", type=str, default="./results",
                       help="Output directory for plots and summary (default: ./results)")
    
    args = parser.parse_args()
    
    # Check if file exists
    if not Path(args.results_file).exists():
        print(f"Error: File {args.results_file} does not exist.")
        return
    
    # Create output directory structure
    output_base = Path(args.output)
    knn_output_dir = output_base / "KNN"
    knn_output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Output directory: {knn_output_dir}")
    
    # Load data
    df = load_results(args.results_file)
    print(f"Data loaded: {len(df)} entries, {df['fichier'].nunique()} unique matrices")
    print(f"Analysis results will be saved to summary file only")
    
    # Generate plots (only if requested)
    if args.plots:
        print("Generating graphs...")
        create_correlation_plots(df, knn_output_dir)
        create_performance_plots(df, knn_output_dir)
        print("Graphs generated successfully!")
    
    # Generate complete summary file (all analysis goes here, nothing to console)
    print("Generating analysis summary...")
    generate_summary(df, knn_output_dir)
    
    print(f"All results saved in: {knn_output_dir}")

if __name__ == "__main__":
    main()
