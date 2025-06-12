import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional
import warnings
warnings.filterwarnings('ignore')

# Graphics configuration
plt.style.use('default')
sns.set_palette("Set2")
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 11

class SeedingAnalyzer:
    """Analyzer for seeding experimentation results."""
    
    def __init__(self, csv_file: str, output_dir: str = "plots"):
        """
        Initialize the analyzer.
        
        Parameters
        ----------
        csv_file : str
            Path to the CSV results file
        output_dir : str
            Output directory for graphics
        """
        self.csv_file = csv_file
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Load data
        self.df = pd.read_csv(csv_file)
        print(f"Data loaded: {len(self.df)} experiments")
        
        # Clean and prepare data
        self._prepare_data()
        
    def _prepare_data(self):
        """Prepare and clean the data."""
        # Calculate derived metrics
        self.df['matrix_size'] = self.df['matrix_m'] * self.df['matrix_n']
        self.df['ilp_efficiency'] = self.df['patterns_found'] / np.maximum(self.df['nb_ilp_calculated'], 1)
        
        print(f"Data prepared: {len(self.df)} valid experiments")
        print(f"X_factor range: {self.df['X_factor'].min()} - {self.df['X_factor'].max()}")
        print(f"step_n range: {self.df['step_n'].min()} - {self.df['step_n'].max()}")
        
    def generate_summary_stats(self):
        """Generate descriptive statistics."""
        print("\n" + "="*60)
        print("DESCRIPTIVE STATISTICS")
        print("="*60)
        
        # General statistics
        print(f"\nTotal number of experiments: {len(self.df)}")
        print(f"Number of unique files: {self.df['read_name'].nunique()}")
        print(f"Number of haplotypes: {sorted(self.df['nb_haplotypes'].unique())}")
        
        # Seeding parameters tested
        print(f"\nX_factor parameters tested: {sorted(self.df['X_factor'].unique())}")
        print(f"step_n parameters tested: {sorted(self.df['step_n'].unique())}")
        
        # Key performance metrics
        print(f"\nPerformance metrics:")
        print(f"Average total time: {self.df['total_time'].mean():.3f}s (±{self.df['total_time'].std():.3f})")
        print(f"Average ILP calls: {self.df['nb_ilp_calculated'].mean():.1f} (±{self.df['nb_ilp_calculated'].std():.1f})")
        print(f"Average patterns found: {self.df['patterns_found'].mean():.1f} (±{self.df['patterns_found'].std():.1f})")
        print(f"Average matrix density: {self.df['matrix_density'].mean():.3f}")
        
        # Seeding usage
        if 'seeding_used' in self.df.columns:
            seeding_used = self.df['seeding_used'].sum()
            print(f"\nExperiments using seeding: {seeding_used}/{len(self.df)} ({seeding_used/len(self.df)*100:.1f}%)")
            
        # Top configurations
        print(f"\n" + "="*40)
        print("TOP 3 CONFIGURATIONS - TOTAL TIME")
        print("="*40)
        best_time = self.df.groupby(['X_factor', 'step_n'])['total_time'].mean().sort_values().head(3)
        for i, ((x, s), time_val) in enumerate(best_time.items()):
            count = len(self.df[(self.df['X_factor'] == x) & (self.df['step_n'] == s)])
            print(f"{i+1}. X_factor={x}, step_n={s}: {time_val:.3f}s (n={count})")
            
        print(f"\n" + "="*40)
        print("TOP 3 CONFIGURATIONS - ILP EFFICIENCY")
        print("="*40)
        best_eff = self.df.groupby(['X_factor', 'step_n'])['ilp_efficiency'].mean().sort_values(ascending=False).head(3)
        for i, ((x, s), eff_val) in enumerate(best_eff.items()):
            count = len(self.df[(self.df['X_factor'] == x) & (self.df['step_n'] == s)])
            print(f"{i+1}. X_factor={x}, step_n={s}: {eff_val:.3f} patterns/ILP (n={count})")
    
    def plot_optimization_landscape(self):
        """Graphics 1: Optimization landscape with time heatmap and efficiency contours."""
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))
        
        # Prepare pivot data
        pivot_time = self.df.groupby(['X_factor', 'step_n'])['total_time'].mean().unstack()
        pivot_eff = self.df.groupby(['X_factor', 'step_n'])['ilp_efficiency'].mean().unstack()
        
        # Time total heatmap
        im = ax.imshow(pivot_time.values, cmap='RdYlBu_r', aspect='auto', alpha=0.8)
        
        # Add time values in heatmap
        for i in range(len(pivot_time.index)):
            for j in range(len(pivot_time.columns)):
                if not np.isnan(pivot_time.iloc[i, j]):
                    text = ax.text(j, i, f'{pivot_time.iloc[i, j]:.2f}s',
                                  ha="center", va="center", color="black", fontsize=10, fontweight='bold')
        
        # Add efficiency contours
        if not pivot_eff.empty and not pivot_eff.isna().all().all():
            X, Y = np.meshgrid(range(len(pivot_eff.columns)), range(len(pivot_eff.index)))
            contours = ax.contour(X, Y, pivot_eff.values, levels=5, colors='white', alpha=0.7, linewidths=2)
            ax.clabel(contours, inline=True, fontsize=9, fmt='%.2f eff')
        
        # Axes configuration
        ax.set_title('Seeding Parameters Optimization Landscape\n(Color=Time, Contours=ILP Efficiency)', 
                    fontsize=14, fontweight='bold', pad=20)
        ax.set_xlabel('step_n', fontsize=12)
        ax.set_ylabel('X_factor', fontsize=12)
        ax.set_xticks(range(len(pivot_time.columns)))
        ax.set_xticklabels(pivot_time.columns)
        ax.set_yticks(range(len(pivot_time.index)))
        ax.set_yticklabels(pivot_time.index)
        
        # Colorbar for time
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Average Total Time (seconds)', fontsize=12)
        
        # Identify and mark optimal points
        best_time_idx = np.unravel_index(np.nanargmin(pivot_time.values), pivot_time.shape)
        best_eff_idx = np.unravel_index(np.nanargmax(pivot_eff.values), pivot_eff.shape)
        
        # Mark best time
        ax.scatter(best_time_idx[1], best_time_idx[0], marker='*', s=200, 
                  c='red', edgecolors='white', linewidth=2, label='Optimal Time')
        
        # Mark best efficiency
        ax.scatter(best_eff_idx[1], best_eff_idx[0], marker='D', s=150, 
                  c='lime', edgecolors='black', linewidth=2, label='Optimal Efficiency')
        
        ax.legend(loc='upper right', bbox_to_anchor=(1, 1))
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'optimization_landscape.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"✅ Optimization landscape saved: {self.output_dir / 'optimization_landscape.png'}")
    
    def plot_3d_seeding_analysis(self):
        """Graphics 2: 3D analysis of step_n, X_factor, and ILP execution time."""
        from mpl_toolkits.mplot3d import Axes3D
        
        fig = plt.figure(figsize=(16, 12))
        
        # Create 2x2 subplot grid for different 3D views
        
        # 1. 3D Surface: step_n, X_factor, ILP time
        ax1 = fig.add_subplot(221, projection='3d')
        pivot_ilp_time = self.df.groupby(['X_factor', 'step_n'])['ilp_time'].mean().unstack()
        if not pivot_ilp_time.empty:
            X, Y = np.meshgrid(pivot_ilp_time.columns, pivot_ilp_time.index)
            Z_ilp = pivot_ilp_time.values
            
            surf1 = ax1.plot_surface(X, Y, Z_ilp, cmap='viridis', alpha=0.8)
            ax1.set_xlabel('step_n')
            ax1.set_ylabel('X_factor')
            ax1.set_zlabel('ILP Time (s)')
            ax1.set_title('ILP Execution Time Surface', fontweight='bold')
            fig.colorbar(surf1, ax=ax1, shrink=0.5)
        
        # 2. 3D Scatter: step_n, X_factor, Total time (colored by efficiency)
        ax2 = fig.add_subplot(222, projection='3d')
        param_data = self.df.groupby(['X_factor', 'step_n']).agg({
            'total_time': 'mean',
            'ilp_efficiency': 'mean',
            'nb_ilp_calculated': 'mean'
        }).reset_index()
        
        scatter = ax2.scatter(param_data['step_n'], param_data['X_factor'], param_data['total_time'],
                             c=param_data['ilp_efficiency'], cmap='plasma', s=100, alpha=0.8)
        ax2.set_xlabel('step_n')
        ax2.set_ylabel('X_factor')
        ax2.set_zlabel('Total Time (s)')
        ax2.set_title('Total Time vs Parameters\n(Color=Efficiency)', fontweight='bold')
        fig.colorbar(scatter, ax=ax2, shrink=0.5, label='ILP Efficiency')
        
        # 3. Analysis of haplotypes impact on ILP performance
        ax3 = fig.add_subplot(223, projection='3d')
        
        # Group by haplotypes and show ILP efficiency vs seeding parameters
        hap_data = []
        for hap_count in sorted(self.df['nb_haplotypes'].unique()):
            hap_df = self.df[self.df['nb_haplotypes'] == hap_count]
            hap_summary = hap_df.groupby(['X_factor', 'step_n']).agg({
                'ilp_efficiency': 'mean',
                'total_time': 'mean'
            }).reset_index()
            hap_summary['haplotypes'] = hap_count
            hap_data.append(hap_summary)
        
        if hap_data:
            combined_hap_data = pd.concat(hap_data, ignore_index=True)
            scatter3 = ax3.scatter(combined_hap_data['X_factor'], 
                                  combined_hap_data['step_n'], 
                                  combined_hap_data['ilp_efficiency'],
                                  c=combined_hap_data['haplotypes'], 
                                  cmap='tab10', s=80, alpha=0.8)
            ax3.set_xlabel('X_factor')
            ax3.set_ylabel('step_n')
            ax3.set_zlabel('ILP Efficiency')
            ax3.set_title('ILP Efficiency by Haplotypes\n(Color=Haplotype Count)', fontweight='bold')
            fig.colorbar(scatter3, ax=ax3, shrink=0.5, label='Haplotype Count')
        
        # 4. Seeding effectiveness: ILP calls reduction analysis
        ax4 = fig.add_subplot(224, projection='3d')
        
        # Calculate seeding benefit (fewer ILP calls = better seeding)
        param_ilp_data = self.df.groupby(['X_factor', 'step_n']).agg({
            'nb_ilp_calculated': 'mean',
            'matrix_size': 'mean',
            'patterns_found': 'mean'
        }).reset_index()
        
        # Normalize ILP calls by matrix complexity
        param_ilp_data['ilp_density'] = param_ilp_data['nb_ilp_calculated'] / param_ilp_data['matrix_size']
        
        # Create 3D surface for ILP density
        if len(param_ilp_data) > 0:
            pivot_ilp_density = param_ilp_data.pivot(index='X_factor', columns='step_n', values='ilp_density')
            if not pivot_ilp_density.empty:
                X, Y = np.meshgrid(pivot_ilp_density.columns, pivot_ilp_density.index)
                Z_density = pivot_ilp_density.values
                
                surf4 = ax4.plot_surface(X, Y, Z_density, cmap='RdYlBu_r', alpha=0.8)
                ax4.set_xlabel('step_n')
                ax4.set_ylabel('X_factor')
                ax4.set_zlabel('ILP Calls per Matrix Element')
                ax4.set_title('Seeding Efficiency Surface\n(Lower = Better)', fontweight='bold')
                fig.colorbar(surf4, ax=ax4, shrink=0.5)
        
        plt.tight_layout()
        plt.savefig(self.output_dir / '3d_seeding_analysis.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"✅ 3D seeding analysis saved: {self.output_dir / '3d_seeding_analysis.png'}")
    
    def generate_performance_report(self):
        """Generate a detailed performance report."""
        report_file = self.output_dir / 'performance_report.txt'
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("="*80 + "\n")
            f.write("SEEDING EXPERIMENTS ANALYSIS REPORT\n")
            f.write("="*80 + "\n\n")
            
            # General information
            f.write("GENERAL INFORMATION\n")
            f.write("-"*30 + "\n")
            f.write(f"Source file: {self.csv_file}\n")
            f.write(f"Number of experiments: {len(self.df)}\n")
            f.write(f"Unique files: {self.df['read_name'].nunique()}\n")
            f.write(f"Haplotypes: {sorted(self.df['nb_haplotypes'].unique())}\n")
            f.write(f"X_factor: {sorted(self.df['X_factor'].unique())}\n")
            f.write(f"step_n: {sorted(self.df['step_n'].unique())}\n\n")
            
            # Global metrics
            f.write("GLOBAL METRICS\n")
            f.write("-"*20 + "\n")
            f.write(f"Average total time: {self.df['total_time'].mean():.3f}s (±{self.df['total_time'].std():.3f})\n")
            f.write(f"Average ILP calls: {self.df['nb_ilp_calculated'].mean():.1f} (±{self.df['nb_ilp_calculated'].std():.1f})\n")
            f.write(f"Average ILP efficiency: {self.df['ilp_efficiency'].mean():.3f} (±{self.df['ilp_efficiency'].std():.3f})\n")
            f.write(f"Average patterns: {self.df['patterns_found'].mean():.1f} (±{self.df['patterns_found'].std():.1f})\n\n")
            
            # Optimal configurations
            f.write("OPTIMAL CONFIGURATIONS\n")
            f.write("-"*25 + "\n")
            
            # Calculate optimal configurations safely
            best_overall = self.df.groupby(['X_factor', 'step_n']).agg({
                'total_time': 'mean',
                'ilp_efficiency': 'mean'
            }).round(3)
            
            # Handle edge cases where there might be only one configuration
            if len(best_overall) == 0:
                f.write("No valid configurations found.\n")
                return
            
            # Combined score calculation with safe handling
            time_range = best_overall['total_time'].max() - best_overall['total_time'].min()
            eff_range = best_overall['ilp_efficiency'].max() - best_overall['ilp_efficiency'].min()
            
            if time_range > 0:
                time_normalized = 1 - (best_overall['total_time'] - best_overall['total_time'].min()) / time_range
            else:
                time_normalized = pd.Series([1.0] * len(best_overall), index=best_overall.index)
                
            if eff_range > 0:
                eff_normalized = (best_overall['ilp_efficiency'] - best_overall['ilp_efficiency'].min()) / eff_range
            else:
                eff_normalized = pd.Series([1.0] * len(best_overall), index=best_overall.index)
            
            combined_score = (time_normalized + eff_normalized) / 2
            
            # Get optimal configurations safely
            optimal_time = self.df.groupby(['X_factor', 'step_n'])['total_time'].mean().idxmin()
            optimal_eff = self.df.groupby(['X_factor', 'step_n'])['ilp_efficiency'].mean().idxmax()
            
            # Handle combined score safely
            if len(combined_score) > 0:
                optimal_combined = combined_score.idxmax()
                f.write(f"⭐ Optimal TIME: X_factor={optimal_time[0]}, step_n={optimal_time[1]}\n")
                f.write(f"⭐ Optimal EFFICIENCY: X_factor={optimal_eff[0]}, step_n={optimal_eff[1]}\n")
                f.write(f"⭐ Optimal COMBINED (time+efficiency): X_factor={optimal_combined[0]}, step_n={optimal_combined[1]}\n\n")
            else:
                f.write(f"⭐ Optimal TIME: X_factor={optimal_time[0]}, step_n={optimal_time[1]}\n")
                f.write(f"⭐ Optimal EFFICIENCY: X_factor={optimal_eff[0]}, step_n={optimal_eff[1]}\n")
                f.write(f"⭐ Optimal COMBINED: Same as time optimal\n\n")
            
            f.write("COMBINED SCORE CALCULATION:\n")
            f.write("Combined Score = (Normalized_Inverted_Time + Normalized_Efficiency) / 2\n")
            f.write("Where:\n")
            f.write("  - Normalized_Inverted_Time = 1 - (time - min_time) / (max_time - min_time)\n")
            f.write("  - Normalized_Efficiency = (efficiency - min_efficiency) / (max_efficiency - min_efficiency)\n\n")
            
            # Top 3 configurations
            f.write("TOP 3 CONFIGURATIONS\n")
            f.write("-"*20 + "\n")
            
            if len(combined_score) > 0:
                top3 = combined_score.sort_values(ascending=False).head(3)
                for i, ((x, s), score) in enumerate(top3.items()):
                    time_val = best_overall.loc[(x, s), 'total_time']
                    eff_val = best_overall.loc[(x, s), 'ilp_efficiency']
                    count = len(self.df[(self.df['X_factor'] == x) & (self.df['step_n'] == s)])
                    f.write(f"{i+1}. X_factor={x}, step_n={s}: {time_val:.3f}s, {eff_val:.3f}eff (score={score:.3f}, n={count})\n")
            else:
                f.write("Insufficient data for ranking configurations.\n")
            
        print(f"✅ Performance report saved: {report_file}")
    
    def run_full_analysis(self):
        """Execute complete analysis with only 2 essential graphics."""
        print("🚀 Starting complete seeding results analysis...")
        
        # Descriptive statistics
        self.generate_summary_stats()
        
        # Generate 2 essential graphics
        print(f"\n📊 Generating visualizations in {self.output_dir}/...")
        
        self.plot_optimization_landscape()
        self.plot_3d_seeding_analysis()
        
        # Final report
        self.generate_performance_report()
        
        print(f"\n✅ Complete analysis finished!")
        print(f"📁 Files saved in: {self.output_dir}/")
        print(f"📋 Detailed report: {self.output_dir}/performance_report.txt")
        
        # Main conclusions summary
        print(f"\n🎯 MAIN CONCLUSIONS:")
        optimal_time = self.df.groupby(['X_factor', 'step_n'])['total_time'].mean().idxmin()
        optimal_eff = self.df.groupby(['X_factor', 'step_n'])['ilp_efficiency'].mean().idxmax()
        print(f"  • Best time config: X_factor={optimal_time[0]}, step_n={optimal_time[1]}")
        print(f"  • Best efficiency config: X_factor={optimal_eff[0]}, step_n={optimal_eff[1]}")
        print(f"  • 2 graphics generated: optimization landscape + 3D seeding analysis")

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='Analyze seeding experimentation results')
    parser.add_argument('csv_file', help='CSV results file')
    parser.add_argument('--output-dir', '-o', default='seeding_analysis_plots', 
                       help='Output directory for graphics')
    parser.add_argument('--summary-only', '-s', action='store_true', 
                       help='Display only statistics, without graphics')
    
    args = parser.parse_args()
    
    # Check file exists
    if not Path(args.csv_file).exists():
        print(f"❌ Error: File {args.csv_file} does not exist")
        sys.exit(1)
    
    # Create analyzer
    analyzer = SeedingAnalyzer(args.csv_file, args.output_dir)
    
    if args.summary_only:
        analyzer.generate_summary_stats()
    else:
        analyzer.run_full_analysis()

if __name__ == "__main__":
    main()
