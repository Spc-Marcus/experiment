import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def load_and_process_data(file1, file2):
    """Load and process data from two CSV files"""
    
    # Load data
    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)
    
    print(f"Available columns in {file1}: {list(df1.columns)}")
    print(f"Available columns in {file2}: {list(df2.columns)}")
    
    # Filter only data with ILP calls (ilp_calls_total > 0)
    df1_ilp = df1[df1['ilp_calls_total'] > 0].copy()
    df2_ilp = df2[df2['ilp_calls_total'] > 0].copy()
    
    print(f"Data with ILP in {file1}: {len(df1_ilp)} / {len(df1)} instances")
    print(f"Data with ILP in {file2}: {len(df2_ilp)} / {len(df2)} instances")
    
    # Find common instances (based on filename)
    common_files = set(df1_ilp['filename']) & set(df2_ilp['filename'])
    print(f"Common instances with ILP: {len(common_files)}")
    
    # Filter to keep only common instances
    df1_common = df1_ilp[df1_ilp['filename'].isin(common_files)].copy()
    df2_common = df2_ilp[df2_ilp['filename'].isin(common_files)].copy()
    
    print(f"Final instances - max_one_comp: {len(df1_common)}")
    print(f"Final instances - MaxOne: {len(df2_common)}")
    
    # Determine which time column to use
    time_col1 = 'execution_time' if 'execution_time' in df1_common.columns else 'ilp_time_total'
    time_col2 = 'execution_time' if 'execution_time' in df2_common.columns else 'ilp_time_total'
    
    if time_col1 not in df1_common.columns:
        raise ValueError(f"Column '{time_col1}' missing in {file1}")
    if time_col2 not in df2_common.columns:
        raise ValueError(f"Column '{time_col2}' missing in {file2}")
    
    print(f"Using '{time_col1}' as execution time for {file1}")
    print(f"Using '{time_col2}' as execution time for {file2}")
    
    return df1_common, df2_common, time_col1, time_col2

def compare_execution_times(df1, df2, time_col1, time_col2, label1="Dataset 1", label2="Dataset 2"):
    """Compare execution times between two datasets"""
    
    # Create figure with multiple subplots
    fig, axes = plt.subplots(2, 3, figsize=(20, 12))
    fig.suptitle(f'ILP Computation Time Comparison - {len(df1)} Common Instances', fontsize=16)
    
    # 1. Average time comparison by haplotype count
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
        ax1.set_xlabel('Number of Haplotypes')
        ax1.set_ylabel('Average Time (s)')
        ax1.set_title('Average Time by Haplotype Count')
        ax1.set_xticks(x)
        ax1.set_xticklabels(haplos)
        ax1.legend()
    else:
        ax1.text(0.5, 0.5, 'Haplotype_count column\nnot available', 
                ha='center', va='center', transform=ax1.transAxes)
    
    # 2. Execution time vs Matrix Density (line plot with bins)
    ax2 = axes[0, 1]
    
    # Create density bins and calculate average execution time
    def calculate_time_by_density_bins(df, time_col, bin_size=0.005):
        df_temp = df.copy()
        # Create smaller density bins (e.g., 0.005 width instead of automatic)
        min_density = df_temp['matrix_density'].min()
        max_density = df_temp['matrix_density'].max()
        bins = np.arange(min_density, max_density + bin_size, bin_size)
        df_temp['density_bins'] = pd.cut(df_temp['matrix_density'], bins=bins)
        grouped = df_temp.groupby('density_bins', observed=False).agg({
            time_col: ['mean', 'count'],
            'matrix_density': 'mean'
        }).reset_index()
        # Flatten column names
        grouped.columns = ['density_bins', f'{time_col}_mean', f'{time_col}_count', 'matrix_density_mean']
        # Only keep bins with at least 2 data points
        return grouped[(grouped[f'{time_col}_mean'].notna()) & (grouped[f'{time_col}_count'] >= 2)]
    
    time_by_density1 = calculate_time_by_density_bins(df1, time_col1)
    time_by_density2 = calculate_time_by_density_bins(df2, time_col2)
    
    if not time_by_density1.empty:
        ax2.plot(time_by_density1['matrix_density_mean'], time_by_density1[f'{time_col1}_mean'], 
                'o-', color='red', label=f'{label1} ({len(time_by_density1)} bins)', linewidth=2, markersize=6)
    
    if not time_by_density2.empty:
        ax2.plot(time_by_density2['matrix_density_mean'], time_by_density2[f'{time_col2}_mean'], 
                'o-', color='blue', label=f'{label2} ({len(time_by_density2)} bins)', linewidth=2, markersize=6)
    
    ax2.set_xlabel('Matrix Density')
    ax2.set_ylabel('Average Execution Time (s)')
    ax2.set_title('Average Execution Time vs Matrix Density')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 3. ILP calls comparison
    ax3 = axes[0, 2]
    merged_ilp = pd.merge(df1[['filename', 'ilp_calls_total']], 
                         df2[['filename', 'ilp_calls_total']], 
                         on='filename', 
                         suffixes=('_1', '_2'))
    
    if not merged_ilp.empty:
        ax3.scatter(merged_ilp['ilp_calls_total_1'], merged_ilp['ilp_calls_total_2'], alpha=0.6)
        # y=x reference line
        max_calls = max(merged_ilp['ilp_calls_total_1'].max(), merged_ilp['ilp_calls_total_2'].max())
        ax3.plot([0, max_calls], [0, max_calls], 'r--', alpha=0.8)
        ax3.set_xlabel(f'ILP Calls {label1}')
        ax3.set_ylabel(f'ILP Calls {label2}')
        ax3.set_title('ILP Calls Comparison')
        
        # Add correlation
        correlation_ilp = merged_ilp['ilp_calls_total_1'].corr(merged_ilp['ilp_calls_total_2'])
        ax3.text(0.05, 0.95, f'Correlation: {correlation_ilp:.3f}', 
                transform=ax3.transAxes, bbox=dict(boxstyle="round", facecolor='lightblue'))
    else:
        ax3.text(0.5, 0.5, 'No common\ninstances', 
                ha='center', va='center', transform=ax3.transAxes)
    
    # 4. Time boxplot comparison
    ax4 = axes[1, 0]
    data_to_plot = [df1[time_col1], df2[time_col2]]
    ax4.boxplot(data_to_plot, labels=[label1, label2])
    ax4.set_ylabel('Execution Time (s)')
    ax4.set_title('Time Distribution Comparison (Boxplot)')
    
    # 5. Execution time vs Matrix Size (replacing ILP calls boxplot)
    ax5 = axes[1, 1]
    
    # Create size bins and calculate average execution time
    def calculate_time_by_size_bins(df, time_col, n_bins=100):
        df_temp = df.copy()
        # Create size bins
        df_temp['size_bins'] = pd.cut(df_temp['matrix_size'], bins=n_bins)
        grouped = df_temp.groupby('size_bins', observed=False).agg({
            time_col: ['mean', 'count'],
            'matrix_size': 'mean'
        }).reset_index()
        # Flatten column names
        grouped.columns = ['size_bins', f'{time_col}_mean', f'{time_col}_count', 'matrix_size_mean']
        # Only keep bins with at least 1 data points
        return grouped[(grouped[f'{time_col}_mean'].notna()) & (grouped[f'{time_col}_count'] >= 1)]
    
    time_by_size1 = calculate_time_by_size_bins(df1, time_col1)
    time_by_size2 = calculate_time_by_size_bins(df2, time_col2)
    
    if not time_by_size1.empty:
        ax5.plot(time_by_size1['matrix_size_mean'], time_by_size1[f'{time_col1}_mean'], 
                'o-', color='red', label=f'{label1} ({len(time_by_size1)} bins)', linewidth=1.5, markersize=4)
    
    if not time_by_size2.empty:
        ax5.plot(time_by_size2['matrix_size_mean'], time_by_size2[f'{time_col2}_mean'], 
                'o-', color='blue', label=f'{label2} ({len(time_by_size2)} bins)', linewidth=1.5, markersize=4)
    
    ax5.set_xlabel('Matrix Size (rows × cols)')
    ax5.set_ylabel('Average Execution Time (s)')
    ax5.set_title('Average Execution Time vs Matrix Size')
    ax5.legend()
    ax5.grid(True, alpha=0.8)
    
    # 6. Direct time comparison scatter plot
    ax6 = axes[1, 2]
    # Merge on filename to compare same instances
    merged = pd.merge(df1[['filename', time_col1]], 
                     df2[['filename', time_col2]], 
                     on='filename', 
                     suffixes=('_1', '_2'))
    
    if not merged.empty:
        ax6.scatter(merged[f'{time_col1}_1'], merged[f'{time_col2}_2'], alpha=0.6)
        # y=x reference line
        max_time = max(merged[f'{time_col1}_1'].max(), merged[f'{time_col2}_2'].max())
        ax6.plot([0, max_time], [0, max_time], 'r--', alpha=0.8)
        ax6.set_xlabel(f'Time {label1} (s)')
        ax6.set_ylabel(f'Time {label2} (s)')
        ax6.set_title('Direct Time Comparison')
        
        # Add correlation
        correlation = merged[f'{time_col1}_1'].corr(merged[f'{time_col2}_2'])
        ax6.text(0.05, 0.95, f'Correlation: {correlation:.3f}', 
                transform=ax6.transAxes, bbox=dict(boxstyle="round", facecolor='wheat'))
    else:
        ax6.text(0.5, 0.5, 'No common\ninstances found', 
                ha='center', va='center', transform=ax6.transAxes)
    
    plt.tight_layout()
    return fig

def generate_statistics(df1, df2, time_col1, time_col2, label1="Dataset 1", label2="Dataset 2"):
    """Generate comparative statistics"""
    
    stats = {
        label1: {
            'count': len(df1),
            'mean_time': df1[time_col1].mean(),
            'median_time': df1[time_col1].median(),
            'std_time': df1[time_col1].std(),
            'min_time': df1[time_col1].min(),
            'max_time': df1[time_col1].max(),
            'q25': df1[time_col1].quantile(0.25),
            'q75': df1[time_col1].quantile(0.75),
            'mean_ilp_calls': df1['ilp_calls_total'].mean(),
            'median_ilp_calls': df1['ilp_calls_total'].median(),
            'max_ilp_calls': df1['ilp_calls_total'].max()
        },
        label2: {
            'count': len(df2),
            'mean_time': df2[time_col2].mean(),
            'median_time': df2[time_col2].median(),
            'std_time': df2[time_col2].std(),
            'min_time': df2[time_col2].min(),
            'max_time': df2[time_col2].max(),
            'q25': df2[time_col2].quantile(0.25),
            'q75': df2[time_col2].quantile(0.75),
            'mean_ilp_calls': df2['ilp_calls_total'].mean(),
            'median_ilp_calls': df2['ilp_calls_total'].median(),
            'max_ilp_calls': df2['ilp_calls_total'].max()
        }
    }
    
    return pd.DataFrame(stats).round(4)

def main():
    """Main function"""
    
    # Paths to your CSV files
    file1 = "experiment_results_max_one_comp.csv"  # Replace with path to your first file
    file2 = "experiment_results_max_one.csv"  # Replace with path to your second file
    
    try:
        # Load data
        print("Loading data...")
        df1_ilp, df2_ilp, time_col1, time_col2 = load_and_process_data(file1, file2)
        
        print(f"max_one_comp: {len(df1_ilp)} instances")
        print(f"MaxOne: {len(df2_ilp)} instances")
        
        # Generate statistics
        print("\n=== COMPARATIVE STATISTICS ===")
        stats_df = generate_statistics(df1_ilp, df2_ilp, time_col1, time_col2, "max_one_comp", "MaxOne")
        print(stats_df)
        
        # Create comparison plots
        print("\nGenerating plots...")
        fig = compare_execution_times(df1_ilp, df2_ilp, time_col1, time_col2, "max_one_comp", "MaxOne")
        
        # Save results
        fig.savefig('comparison_ilp_times.png', dpi=300, bbox_inches='tight')
        stats_df.to_csv('comparison_statistics.csv')
        
        print("Comparison completed!")
        print("- Plots saved: comparison_ilp_times.png")
        print("- Statistics saved: comparison_statistics.csv")
        
        # Quick analysis
        print("\n=== QUICK ANALYSIS ===")
        ratio_mean = stats_df.loc['mean_time', 'MaxOne'] / stats_df.loc['mean_time', 'max_one_comp']
        print(f"Mean time ratio (MaxOne / max_one_comp): {ratio_mean:.2f}")
        
        if ratio_mean > 1:
            print(f"MaxOne is {ratio_mean:.2f}x slower on average")
        else:
            print(f"MaxOne is {1/ratio_mean:.2f}x faster on average")
        
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
        print("Make sure the paths to your CSV files are correct")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()