import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse
import sys
from io import StringIO
import contextlib

def load_results(csv_file: str) -> pd.DataFrame:
    """Load experiment results from CSV file."""
    csv_path = Path(csv_file)
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV file not found: {csv_path}")
    if not csv_path.suffix.lower() == '.csv':
        raise ValueError(f"File must be a CSV file: {csv_path}")
    
    df = pd.read_csv(csv_file)
    
    # Validate required columns
    required_columns = ['haplotype_count', 'ilp_calls_total', 'matrix_size', 'total_time']
    missing_columns = [col for col in required_columns if col not in df.columns]
    
    if missing_columns:
        print(f"Warning: Missing required columns: {missing_columns}")
        print(f"Available columns: {list(df.columns)}")
        raise ValueError(f"CSV file missing required columns: {missing_columns}")
    
    # Handle missing optional columns by creating them with default values
    if 'patterns_found' not in df.columns:
        df['patterns_found'] = 0
    if 'matrix_density' not in df.columns:
        df['matrix_density'] = df.get('ones_count', 0) / df.get('matrix_size', 1)
    if 'regions_found' not in df.columns:
        df['regions_found'] = 0
    if 'ilp_calls_per_pattern' not in df.columns:
        df['ilp_calls_per_pattern'] = df['ilp_calls_total'] / df['patterns_found'].replace(0, 1)
    if 'patterns_per_second' not in df.columns:
        df['patterns_per_second'] = df['patterns_found'] / df['total_time'].replace(0, 1)
    
    return df

def analyze_performance_by_haplotype(df: pd.DataFrame):
    """Analyze performance metrics by haplotype count."""
    
    print("\n=== PERFORMANCE BY HAPLOTYPE COUNT ===")
    
    # Analyze ALL data with correct interpretation
    print(f"Total runs: {len(df)}")
    print(f"Runs requiring ILP optimization: {len(df[df['ilp_calls_total'] > 0])}")
    print(f"Runs solved without ILP (efficient): {len(df[df['ilp_calls_total'] == 0])}")
    print(f"ILP utilization rate: {len(df[df['ilp_calls_total'] > 0]) / len(df) * 100:.1f}%")
    print(f"Haplotype counts found: {sorted(df['haplotype_count'].unique())}")
    
    # ILP utilization rate by haplotype count (higher = more ILP usage)
    utilization_by_hap = df.groupby('haplotype_count').agg({
        'ilp_calls_total': ['count', lambda x: (x > 0).sum(), lambda x: (x == 0).sum()]
    })
    utilization_by_hap.columns = ['total_runs', 'runs_with_ilp', 'runs_without_ilp']
    utilization_by_hap['ilp_utilization_rate_%'] = utilization_by_hap['runs_with_ilp'] / utilization_by_hap['total_runs'] * 100
    
    print(f"\n=== ILP UTILIZATION RATE BY HAPLOTYPE COUNT ===")
    print("(Higher utilization = more problems requiring ILP optimization)")
    print(utilization_by_hap.round(2))
    
    # Detailed statistics for ALL runs
    all_stats = df.groupby('haplotype_count').agg({
        'matrix_rows': ['count', 'mean', 'std', 'min', 'max'],
        'matrix_cols': ['mean', 'std', 'min', 'max'],
        'matrix_size': ['mean', 'std', 'min', 'max'],
        'matrix_density': ['mean', 'std', 'min', 'max'],
        'total_time': ['mean', 'std', 'min', 'max'],
        'ilp_calls_total': ['mean', 'std', 'min', 'max'],
        'patterns_found': ['mean', 'std', 'min', 'max']
    }).round(3)
    
    print(f"\n=== ALL RUNS STATISTICS BY HAPLOTYPE COUNT ===")
    print(all_stats)
    
    # Statistics for ILP-requiring runs only
    ilp_runs_df = df[df['ilp_calls_total'] > 0].copy()
    if len(ilp_runs_df) > 0:
        ilp_stats = ilp_runs_df.groupby('haplotype_count').agg({
            'matrix_rows': ['count', 'mean', 'std'],
            'matrix_cols': ['mean', 'std'],
            'matrix_density': ['mean', 'std'],
            'ilp_calls_total': ['mean', 'std', 'min', 'max'],
            'total_time': ['mean', 'std', 'min', 'max'],
            'patterns_found': ['mean', 'std'],
            'regions_found': ['mean', 'std']
        }).round(3)
        
        print(f"\n=== ILP-REQUIRING RUNS STATISTICS ===")
        print("(Only runs that needed ILP optimization)")
        print(ilp_stats)
    
    # Statistics for efficient runs (no ILP needed)
    efficient_df = df[df['ilp_calls_total'] == 0].copy()
    if len(efficient_df) > 0:
        print(f"\n=== EFFICIENT RUNS CHARACTERISTICS (No ILP Needed) ===")
        efficient_stats = efficient_df.groupby('haplotype_count').agg({
            'matrix_size': ['count', 'mean', 'std', 'min', 'max'],
            'matrix_density': ['mean', 'std'],
            'total_time': ['mean', 'std'],
            'patterns_found': ['mean', 'std']
        }).round(3)
        print(efficient_stats)
        
        # Show examples of efficient runs
        print(f"\nExamples of efficient runs (solved without ILP):")
        for hap_count in sorted(efficient_df['haplotype_count'].unique()):
            hap_efficient = efficient_df[efficient_df['haplotype_count'] == hap_count].head(3)
            print(f"\nHaplotype {hap_count}:")
            for _, row in hap_efficient.iterrows():
                print(f"  {row['filename']}: {row['matrix_rows']}x{row['matrix_cols']}, density={row['matrix_density']:.3f}, "
                      f"patterns={row['patterns_found']}, time={row['total_time']:.3f}s")
    
    return all_stats

def analyze_scalability(df: pd.DataFrame):
    """Analyze scalability with matrix size - including all data points."""
    
    print("\n=== SCALABILITY ANALYSIS (ALL DATA) ===")
    
    # Add matrix complexity categories for ALL data
    df['size_category'] = pd.cut(df['matrix_size'], 
                               bins=5, 
                               labels=['Very Small', 'Small', 'Medium', 'Large', 'Very Large'])
    
    # Analysis by size category for all runs
    size_stats_all = df.groupby('size_category').agg({
        'matrix_size': ['count', 'mean'],
        'haplotype_count': ['mean', 'std'],
        'matrix_density': ['mean', 'std'],
        'total_time': ['mean', 'std'],
        'ilp_calls_total': ['mean', 'std', lambda x: (x == 0).sum(), lambda x: (x > 0).sum()],
        'patterns_found': ['mean', 'std']
    }).round(3)
    
    print("Statistics by matrix size (all runs):")
    print(size_stats_all)
    
    # Analysis for ILP-requiring runs only
    ilp_runs_df = df[df['ilp_calls_total'] > 0].copy()
    if len(ilp_runs_df) > 0:
        ilp_runs_df['size_category'] = pd.cut(ilp_runs_df['matrix_size'], 
                                           bins=5, 
                                           labels=['Very Small', 'Small', 'Medium', 'Large', 'Very Large'])
        
        # Only include columns that exist in the dataframe
        agg_dict = {
            'matrix_size': ['count', 'mean'],
            'ilp_calls_total': ['mean', 'std'],
            'total_time': ['mean', 'std']
        }
        
        if 'ilp_calls_per_pattern' in ilp_runs_df.columns:
            agg_dict['ilp_calls_per_pattern'] = ['mean', 'std']
        if 'patterns_per_second' in ilp_runs_df.columns:
            agg_dict['patterns_per_second'] = ['mean', 'std']
            
        size_stats_ilp = ilp_runs_df.groupby('size_category').agg(agg_dict).round(3)
        
        print(f"\nStatistics by matrix size (ILP-requiring runs only):")
        print(size_stats_ilp)
        
        # Correlation analysis for ILP-requiring runs
        print(f"\n=== CORRELATIONS (ILP-REQUIRING RUNS) ===")
        corr_columns = ['matrix_size', 'matrix_density', 'haplotype_count', 
                       'ilp_calls_total', 'total_time', 'patterns_found']
        # Only use columns that exist
        corr_columns = [col for col in corr_columns if col in ilp_runs_df.columns]
        
        if len(corr_columns) >= 2:
            correlations = ilp_runs_df[corr_columns].corr()
            print(correlations.round(3))
        else:
            print("Insufficient columns for correlation analysis")
    
    # Overall correlations including all runs
    print(f"\n=== CORRELATIONS (ALL RUNS) ===")
    corr_columns = ['matrix_size', 'matrix_density', 'haplotype_count', 
                   'ilp_calls_total', 'total_time', 'patterns_found']
    # Only use columns that exist
    corr_columns = [col for col in corr_columns if col in df.columns]
    
    if len(corr_columns) >= 2:
        correlations_all = df[corr_columns].corr()
        print(correlations_all.round(3))
    else:
        print("Insufficient columns for correlation analysis")

@contextlib.contextmanager
def capture_output():
    """Context manager to capture print output."""
    old_stdout = sys.stdout
    sys.stdout = captured_output = StringIO()
    try:
        yield captured_output
    finally:
        sys.stdout = old_stdout

def save_summary_to_file(summary_text: str, output_dir: str):
    """Save analysis summary to text file."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    summary_file = output_path / "analysis_summary.txt"
    
    with open(summary_file, 'w', encoding='utf-8') as f:
        f.write("ILP HAPLOTYPE ANALYSIS - EXPERIMENT SUMMARY\n")
        f.write("=" * 60 + "\n\n")
        f.write(summary_text)
    
    print(f"Analysis summary saved to: {summary_file}")

def create_visualizations(df: pd.DataFrame, output_dir: str = "plots", subdir: str = "ilp_call"):
    """Create performance visualization plots for ALL data with correct interpretation."""
    
    # Create custom subdirectory
    ilp_output_dir = Path(output_dir) / subdir
    ilp_output_dir.mkdir(parents=True, exist_ok=True)
    
    # Set up plotting style
    plt.style.use('default')
    sns.set_palette("husl")
    
    print(f"Creating visualizations for {len(df)} total data points...")
    print(f"Saving plots to: {ilp_output_dir.resolve()}")
    
    # 1. Matrix size distribution by haplotype count (ALL data)
    plt.figure(figsize=(12, 8))
    for hap_count in sorted(df['haplotype_count'].unique()):
        hap_data = df[df['haplotype_count'] == hap_count]
        plt.scatter(hap_data['matrix_size'], hap_data['total_time'], 
                   label=f'{hap_count} haplotypes', alpha=0.7, s=60)
    
    plt.xlabel('Matrix Size (rows × cols)')
    plt.ylabel('Total Execution Time (seconds)')
    plt.title('Execution Time vs Matrix Size (All Runs)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(f"{ilp_output_dir}/time_vs_size_all_data.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Separate ILP-requiring and efficient runs
    ilp_runs_df = df[df['ilp_calls_total'] > 0].copy()
    efficient_df = df[df['ilp_calls_total'] == 0].copy()
    
    if len(ilp_runs_df) > 0 and len(efficient_df) > 0:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
        
        # ILP-requiring runs
        scatter1 = ax1.scatter(ilp_runs_df['matrix_size'], ilp_runs_df['total_time'], 
                             c=ilp_runs_df['haplotype_count'], cmap='viridis', alpha=0.7, s=60)
        ax1.set_xlabel('Matrix Size')
        ax1.set_ylabel('Total Execution Time (seconds)')
        ax1.set_title(f'Complex Runs (ILP Required, n={len(ilp_runs_df)})')
        ax1.grid(True, alpha=0.3)
        plt.colorbar(scatter1, ax=ax1, label='Haplotype Count')
        
        # Efficient runs
        scatter2 = ax2.scatter(efficient_df['matrix_size'], efficient_df['total_time'], 
                             c=efficient_df['haplotype_count'], cmap='viridis', alpha=0.7, s=60)
        ax2.set_xlabel('Matrix Size')
        ax2.set_ylabel('Total Execution Time (seconds)')
        ax2.set_title(f'Efficient Runs (No ILP Needed, n={len(efficient_df)})')
        ax2.grid(True, alpha=0.3)
        plt.colorbar(scatter2, ax=ax2, label='Haplotype Count')
        
        plt.tight_layout()
        plt.savefig(f"{ilp_output_dir}/complex_vs_efficient_comparison.png", dpi=300, bbox_inches='tight')
        plt.close()
    
    # 3. ILP calls distribution
    plt.figure(figsize=(12, 8))
    
    # Create bins for ILP calls (including 0)
    ilp_bins = [0, 1, 5, 10, 25, 50, 100, float('inf')]
    ilp_labels = ['0 (Efficient)', '1-4', '5-9', '10-24', '25-49', '50-99', '100+']
    
    df['ilp_category'] = pd.cut(df['ilp_calls_total'], bins=ilp_bins, labels=ilp_labels, include_lowest=True)
    
    ilp_dist = df['ilp_category'].value_counts().sort_index()
    
    colors = ['green'] + ['orange'] * (len(ilp_labels) - 1)  # Green for efficient (0), orange for others
    plt.bar(range(len(ilp_dist)), ilp_dist.values, color=colors, alpha=0.7)
    plt.xlabel('ILP Calls Required')
    plt.ylabel('Number of Runs')
    plt.title('Distribution of ILP Optimization Requirements')
    plt.xticks(range(len(ilp_dist)), ilp_dist.index, rotation=45)
    plt.grid(True, alpha=0.3)
    
    # Add count labels on bars
    for i, v in enumerate(ilp_dist.values):
        plt.text(i, v + 0.5, str(v), ha='center')
    
    plt.tight_layout()
    plt.savefig(f"{ilp_output_dir}/ilp_calls_distribution.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # 4. ILP Usage vs Matrix Characteristics
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # ILP Utilization vs Matrix Size
    df['requires_ilp'] = (df['ilp_calls_total'] > 0).astype(int)
    
    size_utilization = df.groupby(pd.cut(df['matrix_size'], bins=10))['requires_ilp'].mean()
    
    size_bins = []
    size_utilization_values = []
    
    for interval, utilization in size_utilization.items():
        if not pd.isna(utilization):  # Skip NaN values
            bin_center = (interval.left + interval.right) / 2
            size_bins.append(bin_center)
            size_utilization_values.append(utilization)
    
    axes[0,0].plot(size_bins, np.array(size_utilization_values) * 100, 'ro-', linewidth=2, markersize=8)
    axes[0,0].set_xlabel('Matrix Size')
    axes[0,0].set_ylabel('ILP Utilization Rate (%)')
    axes[0,0].set_title('ILP Utilization vs Matrix Size')
    axes[0,0].grid(True, alpha=0.3)
    axes[0,0].set_ylim(bottom=0)  # Start from 0%
    
    # ILP Utilization vs Matrix Density
    density_utilization = df.groupby(pd.cut(df['matrix_density'], bins=10))['requires_ilp'].mean()
    
    # Get bin centers for plotting
    density_bins = []
    utilization_values = []
    
    for interval, utilization in density_utilization.items():
        if not pd.isna(utilization):  # Skip NaN values
            bin_center = (interval.left + interval.right) / 2
            density_bins.append(bin_center)
            utilization_values.append(utilization)
    
    axes[0,1].plot(density_bins, np.array(utilization_values) * 100, 'bo-', linewidth=2, markersize=8)
    axes[0,1].set_xlabel('Matrix Density')
    axes[0,1].set_ylabel('ILP Utilization Rate (%)')
    axes[0,1].set_title('ILP Utilization vs Matrix Density')
    axes[0,1].grid(True, alpha=0.3)
    axes[0,1].set_ylim(bottom=0)  # Start from 0%
    
    # ILP Utilization vs Haplotype Count
    hap_utilization = df.groupby('haplotype_count')['requires_ilp'].mean()
    axes[1,0].bar(hap_utilization.index, hap_utilization.values * 100, color='purple', alpha=0.7)
    axes[1,0].set_xlabel('Haplotype Count')
    axes[1,0].set_ylabel('ILP Utilization Rate (%)')
    axes[1,0].set_title('ILP Utilization vs Haplotype Count')
    axes[1,0].grid(True, alpha=0.3)
    axes[1,0].set_ylim(bottom=0)  # Start from 0%
    
    # ILP Calls vs Matrix Complexity (for non-zero cases)
    if len(ilp_runs_df) > 0:
        ilp_runs_df['complexity'] = ilp_runs_df['matrix_size'] * ilp_runs_df['matrix_density']
        axes[1,1].scatter(ilp_runs_df['complexity'], ilp_runs_df['ilp_calls_total'], alpha=0.6)
        axes[1,1].set_xlabel('Matrix Complexity (Size × Density)')
        axes[1,1].set_ylabel('ILP Calls Required')
        axes[1,1].set_title('ILP Calls vs Matrix Complexity')
        axes[1,1].grid(True, alpha=0.3)
        axes[1,1].set_ylim(bottom=0)  # Start from 0
    
    plt.tight_layout()
    plt.savefig(f"{ilp_output_dir}/ilp_utilization_analysis.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Visualization plots saved to {ilp_output_dir}/")
    print(f"- time_vs_size_all_data.png: All data points execution times")
    print(f"- complex_vs_efficient_comparison.png: ILP-requiring vs efficient runs")
    print(f"- ilp_calls_distribution.png: Distribution of ILP optimization requirements")
    print(f"- ilp_utilization_analysis.png: ILP utilization patterns vs matrix characteristics")

def analyze_algorithm_efficiency(df: pd.DataFrame):
    """Analyze algorithm efficiency patterns."""
    
    print("\n=== ALGORITHM EFFICIENCY ANALYSIS ===")
    
    efficient_df = df[df['ilp_calls_total'] == 0].copy()
    ilp_runs_df = df[df['ilp_calls_total'] > 0].copy()
    
    print(f"Total runs using ILP: {len(ilp_runs_df)} out of {len(df)} runs ({len(ilp_runs_df)/len(df)*100:.1f}%)")
    
    if len(ilp_runs_df) == 0:
        print("No ILP-requiring runs to analyze!")
        return
    
    # 1. Characteristics that lead to ILP usage
    print(f"\n=== CHARACTERISTICS OF ILP-REQUIRING RUNS ===")
    
    # Matrix size analysis
    size_utilization_analysis = df.groupby(pd.cut(df['matrix_size'], 
                                            bins=[0, 10000, 20000, 40000, 80000, float('inf')],
                                            labels=['Very Small (<10K)', 'Small (10-20K)', 'Medium (20-40K)', 'Large (40-80K)', 'Very Large (>80K)'])).agg({
        'ilp_calls_total': ['count', lambda x: (x > 0).sum(), lambda x: (x > 0).sum() / len(x) * 100]
    }).round(2)
    size_utilization_analysis.columns = ['total_runs', 'ilp_requiring_runs', 'ilp_utilization_rate_%']
    
    print("ILP utilization rate by matrix size:")
    print(size_utilization_analysis)
    
    # Matrix density analysis
    density_utilization_analysis = df.groupby(pd.cut(df['matrix_density'], 
                                               bins=[0, 0.7, 0.8, 0.9, 0.95, 1.0],
                                               labels=['Low (<0.7)', 'Medium (0.7-0.8)', 'High (0.8-0.9)', 'Very High (0.9-0.95)', 'Extreme (>0.95)'])).agg({
        'ilp_calls_total': ['count', lambda x: (x > 0).sum(), lambda x: (x > 0).sum() / len(x) * 100]
    }).round(2)
    density_utilization_analysis.columns = ['total_runs', 'ilp_requiring_runs', 'ilp_utilization_rate_%']
    
    print(f"\nILP utilization rate by matrix density:")
    print(density_utilization_analysis)
    
    # 2. Compare efficient vs ILP-requiring runs characteristics
    if len(ilp_runs_df) > 0:
        print(f"\n=== COMPARISON: EFFICIENT vs ILP-REQUIRING RUNS ===")
        
        comparison_stats = pd.DataFrame({
            'Efficient_Runs_Mean': efficient_df[['matrix_size', 'matrix_density', 'matrix_rows', 'matrix_cols', 'total_time', 'patterns_found']].mean(),
            'Efficient_Runs_Std': efficient_df[['matrix_size', 'matrix_density', 'matrix_rows', 'matrix_cols', 'total_time', 'patterns_found']].std(),
            'ILP_Runs_Mean': ilp_runs_df[['matrix_size', 'matrix_density', 'matrix_rows', 'matrix_cols', 'total_time', 'patterns_found']].mean(),
            'ILP_Runs_Std': ilp_runs_df[['matrix_size', 'matrix_density', 'matrix_rows', 'matrix_cols', 'total_time', 'patterns_found']].std()
        }).round(3)
        
        print(comparison_stats)
        
        # Statistical significance test
        from scipy import stats
        
        print(f"\n=== STATISTICAL TESTS (Efficient vs ILP-requiring) ===")
        
        for metric in ['matrix_size', 'matrix_density', 'total_time']:
            if metric in efficient_df.columns and metric in ilp_runs_df.columns:
                eff_values = efficient_df[metric].dropna()
                ilp_values = ilp_runs_df[metric].dropna()
                
                if len(eff_values) > 1 and len(ilp_values) > 1:
                    statistic, p_value = stats.mannwhitneyu(eff_values, ilp_values, alternative='two-sided')
                    print(f"{metric}: p-value = {p_value:.6f} {'(significant)' if p_value < 0.05 else '(not significant)'}")
    
    # 3. Identify optimal matrix characteristics for efficiency
    print(f"\n=== OPTIMAL CHARACTERISTICS FOR EFFICIENCY ===")
    
    if len(efficient_df) > 5:  # Need enough data points
        print("Matrix characteristics that promote efficiency (based on efficient runs):")
        print(f"  Typical matrix size: {efficient_df['matrix_size'].mean():.0f} ± {efficient_df['matrix_size'].std():.0f}")
        print(f"  Typical matrix density: {efficient_df['matrix_density'].mean():.3f} ± {efficient_df['matrix_density'].std():.3f}")
        print(f"  Typical dimensions: {efficient_df['matrix_rows'].mean():.0f} × {efficient_df['matrix_cols'].mean():.0f}")
        print(f"  Average execution time: {efficient_df['total_time'].mean():.3f}s ± {efficient_df['total_time'].std():.3f}s")
        print(f"  Average patterns found: {efficient_df['patterns_found'].mean():.1f} ± {efficient_df['patterns_found'].std():.1f}")
        
        # Find the "sweet spot" ranges
        print(f"\n=== SWEET SPOT RANGES (95% of efficient runs) ===")
        for metric, name in [('matrix_size', 'Matrix Size'), ('matrix_density', 'Matrix Density'), 
                           ('matrix_rows', 'Matrix Rows'), ('matrix_cols', 'Matrix Columns')]:
            q2_5 = efficient_df[metric].quantile(0.025)
            q97_5 = efficient_df[metric].quantile(0.975)
            print(f"  {name}: {q2_5:.0f} to {q97_5:.0f}" if metric != 'matrix_density' 
                  else f"  {name}: {q2_5:.3f} to {q97_5:.3f}")

def main():
    parser = argparse.ArgumentParser(
        description='Analyze ILP haplotype analysis experiment results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python ILP_call.py ../data/experiments/ilp_call_stats.csv
  python ILP_call.py ./results/experiment_results.csv --plots
  python ILP_call.py /path/to/results.csv --plots --output ./analysis_plots
  python ILP_call.py data.csv --plots --subdir my_ilp_analysis
        """
    )
    parser.add_argument('csv_file', type=str,
                       help='Path to CSV file with experiment results')
    parser.add_argument('--plots', action='store_true', 
                       help='Generate visualization plots')
    parser.add_argument('--output', '-o', type=str, default="./results",
                       help='Output directory (default: ./results)')
    parser.add_argument('--subdir', type=str, default="ilp_call",
                       help='Subdirectory name for ILP analysis results (default: ilp_call)')
    
    args = parser.parse_args()
    
    # Validate CSV file path
    try:
        csv_path = Path(args.csv_file)
        if not csv_path.exists():
            print(f"Error: CSV file not found: {csv_path}")
            return 1
        if not csv_path.suffix.lower() == '.csv':
            print(f"Error: File must be a CSV file: {csv_path}")
            return 1
            
        print(f"Loading experiment results from: {csv_path.resolve()}")
        
        # Load and analyze results
        df = load_results(args.csv_file)
        
    except Exception as e:
        print(f"Error loading CSV file: {e}")
        return 1
    
    print(f"Loaded {len(df)} experiment results from {args.csv_file}")
    
    # Create output directory for ILP analysis with custom subdirectory name
    ilp_output_dir = Path(args.output) / args.subdir
    ilp_output_dir.mkdir(parents=True, exist_ok=True)
    
    # Capture all analysis output to save to file
    with capture_output() as captured:
        print(f"Total runs: {len(df)}")
        print(f"Runs without ILP: {len(df[df['ilp_calls_total'] == 0])}")
        print(f"Runs requiring ILP: {len(df[df['ilp_calls_total'] > 0])}")
        print(f"ILP utilization rate: {len(df[df['ilp_calls_total'] > 0]) / len(df) * 100:.1f}%")
        
        # Detailed analysis
        analyze_performance_by_haplotype(df)
        analyze_scalability(df)
        analyze_algorithm_efficiency(df)
    
    # Save summary to text file
    summary_text = captured.getvalue()
    save_summary_to_file(summary_text, ilp_output_dir)
    
    # Generate plots if requested
    if args.plots:
        create_visualizations(df, args.output, args.subdir)
    
    print(f"Analysis completed. Results saved to: {ilp_output_dir.resolve()}")
    
    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())
