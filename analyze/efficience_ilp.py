import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

class ILPEfficiencyAnalyzer:
    def __init__(self, csv_file_path=None):
        """Initialize the ILP analyzer with CSV data if provided"""
        self.csv_file_path = csv_file_path
        self.df = None
        if csv_file_path and Path(csv_file_path).exists():
            self.load_data()
    
    def load_data(self):
        """Load data from CSV file with ILP-specific columns"""
        print(f"Loading ILP data from {self.csv_file_path}")
        
        try:
            # Define expected data types for ILP analysis
            dtype_dict = {
                'Threshold': 'float64',
                'Error-Rate': 'float64',
                'Strip': 'int64',
                'Haplotype': 'int64',
                'Matrix-Size-cols': 'int64',
                'Matrix-Size-rows': 'int64',
                'Time-Pre-Processing': 'float64',
                'Steps-Count-Pre': 'int64',
                'Unused-Cols-Pre': 'int64',
                'Final-Clusters-Pre': 'int64',
                'Orphan-Reads-Pre': 'int64',
                'Steps-Count-ILP': 'int64',
                'Final-Clusters-ILP': 'int64',
                'Orphan-Reads-ILP': 'int64',
                'Distance': 'float64'
            }
            
            self.df = pd.read_csv(self.csv_file_path, dtype=dtype_dict)
            print(f"Successfully loaded {len(self.df)} records")
            
            # Calculate derived metrics
            self.calculate_efficiency_metrics()
            
            print("Data overview:")
            print(self.df.dtypes)
            print("\nSample data:")
            print(self.df.head())
            
            print(f"\nData ranges:")
            print(f"Haplotype: {self.df['Haplotype'].min()} to {self.df['Haplotype'].max()}")
            print(f"Error-Rate: {self.df['Error-Rate'].min()} to {self.df['Error-Rate'].max()}")
            print(f"Threshold: {self.df['Threshold'].min()} to {self.df['Threshold'].max()}")
            print(f"Distance: {self.df['Distance'].min()} to {self.df['Distance'].max()}")
            
        except Exception as e:
            print(f"Error loading data: {e}")
            self.df = None
    
    def calculate_efficiency_metrics(self):
        """Calculate efficiency metrics comparing pre and post ILP"""
        if self.df is None:
            return
        
        # Success rates (clusters match expected haplotypes)
        self.df['Success_Rate_Pre'] = (self.df['Final-Clusters-Pre'] == self.df['Haplotype']).astype(int)
        self.df['Success_Rate_ILP'] = (self.df['Final-Clusters-ILP'] == self.df['Haplotype']).astype(int)
        
        # Improvement metrics
        self.df['Clusters_Improvement'] = self.df['Final-Clusters-ILP'] - self.df['Final-Clusters-Pre']
        self.df['Orphan_Reduction'] = self.df['Orphan-Reads-Pre'] - self.df['Orphan-Reads-ILP']
        self.df['Steps_Reduction'] = self.df['Steps-Count-Pre'] - self.df['Steps-Count-ILP']
        
        # Relative improvements
        self.df['Orphan_Reduction_Rate'] = np.where(
            self.df['Orphan-Reads-Pre'] > 0,
            self.df['Orphan_Reduction'] / self.df['Orphan-Reads-Pre'],
            0
        )
        
        print("\nEfficiency metrics calculated:")
        print(f"Pre-ILP success rate: {self.df['Success_Rate_Pre'].mean():.3f}")
        print(f"Post-ILP success rate: {self.df['Success_Rate_ILP'].mean():.3f}")
        print(f"Average cluster improvement: {self.df['Clusters_Improvement'].mean():.3f}")
        print(f"Average orphan reduction: {self.df['Orphan_Reduction'].mean():.3f}")
    
    def plot_success_rate_comparison(self, save_path=None, threshold=None, error_rate=None, distance=None, show=True):
        """Compare success rates before and after ILP by haplotype count.

        Étapes:
        1) (Optionnel) Filtrer sur Threshold, Error-Rate, Distance pour comparer à paramètres constants.
        2) S'assurer que les métriques de succès existent (Pre/ILP).
        3) Agréger la moyenne par Haplotype (Pre vs Post ILP).
        4) Tracer des barres groupées, avec un titre indiquant les filtres (ou la moyenne sur plusieurs valeurs).
        """
        if self.df is None or len(self.df) == 0:
            print("No data available for success rate comparison")
            return

        # 1) Appliquer les filtres optionnels (comparaison à paramètres constants)
        df_plot = self.df.copy()
        filters_desc = []
        if threshold is not None:
            df_plot = df_plot[df_plot['Threshold'] == threshold]
            filters_desc.append(f"Threshold={threshold}")
        if error_rate is not None:
            df_plot = df_plot[df_plot['Error-Rate'] == error_rate]
            filters_desc.append(f"Error-Rate={error_rate}")
        if distance is not None:
            df_plot = df_plot[df_plot['Distance'] == distance]
            filters_desc.append(f"Distance={distance}")

        if df_plot.empty:
            print("Aucune donnée après filtrage (vérifier Threshold/Error-Rate/Distance).")
            return

        # Information si on moyenne encore sur des valeurs multiples non fixées
        notes = []
        if threshold is None and df_plot['Threshold'].nunique() > 1:
            notes.append(f"moy. sur {df_plot['Threshold'].nunique()} Thresholds")
        if error_rate is None and df_plot['Error-Rate'].nunique() > 1:
            notes.append(f"moy. sur {df_plot['Error-Rate'].nunique()} Error-Rates")
        if distance is None and df_plot['Distance'].nunique() > 1:
            notes.append(f"moy. sur {df_plot['Distance'].nunique()} Distances")

        # 2) S'assurer que les colonnes de succès existent
        if 'Success_Rate_Pre' not in df_plot.columns or 'Success_Rate_ILP' not in df_plot.columns:
            self.calculate_efficiency_metrics()
            df_plot = self.df.copy()  # recalc appliqué au DataFrame complet, réappliquer filtres
            if threshold is not None:
                df_plot = df_plot[df_plot['Threshold'] == threshold]
            if error_rate is not None:
                df_plot = df_plot[df_plot['Error-Rate'] == error_rate]
            if distance is not None:
                df_plot = df_plot[df_plot['Distance'] == distance]
            if df_plot.empty:
                print("Aucune donnée après recalcul et filtrage.")
                return

        # 3) Agrégation par haplotype
        haplotype_counts = sorted(df_plot['Haplotype'].unique())
        if len(haplotype_counts) == 0:
            print("No haplotype values found")
            return

        success_data = []
        for hap_count in haplotype_counts:
            hap_data = df_plot[df_plot['Haplotype'] == hap_count]
            success_data.append({
                'Haplotype': hap_count,
                'Pre-ILP': hap_data['Success_Rate_Pre'].mean(),
                'Post-ILP': hap_data['Success_Rate_ILP'].mean(),
                'Count': len(hap_data)
            })
        success_df = pd.DataFrame(success_data)

        # 4) Tracé
        plt.figure(figsize=(12, 6))
        x = np.arange(len(haplotype_counts))
        width = 0.35

        plt.bar(x - width/2, success_df['Pre-ILP'], width, label='Pre-ILP', alpha=0.8)
        plt.bar(x + width/2, success_df['Post-ILP'], width, label='Post-ILP', alpha=0.8)

        base_title = "Comparaison du Taux de Succès: Avant vs Après ILP"
        if filters_desc:
            base_title += " (" + ", ".join(filters_desc) + ")"
        plt.title(base_title)
        if notes:
            plt.suptitle(" ; ".join(notes), y=0.98, fontsize=9)

        plt.xlabel("Nombre d'Haplotypes")
        plt.ylabel('Taux de Succès')
        plt.xticks(x, haplotype_counts)
        plt.ylim(0, 1)
        plt.legend()
        plt.grid(axis='y', alpha=0.3)

        # Labels
        for i, (pre, post) in enumerate(zip(success_df['Pre-ILP'], success_df['Post-ILP'])):
            plt.text(i - width/2, min(pre + 0.01, 0.99), f'{pre:.2f}', ha='center', va='bottom')
            plt.text(i + width/2, min(post + 0.01, 0.99), f'{post:.2f}', ha='center', va='bottom')

        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        if show:
            plt.show()
        else:
            plt.close()

    def export_all_success_rate_by_params(self, output_dir="/home/mafoin/stage/experiment/analyze/plots_ilp/success_rate_by_param"):
        """Exporte des graphiques de 'success rate' pour:
        - chaque valeur de Threshold (moyenne sur Error-Rate/Distance),
        - chaque valeur de Error-Rate (moyenne sur Threshold/Distance),
        - chaque valeur de Distance (moyenne sur Threshold/Error-Rate).
        Les graphiques sont enregistrés, non affichés.
        """
        if self.df is None or len(self.df) == 0:
            print("No valid data to export")
            return

        out_base = Path(output_dir)
        (out_base / "threshold").mkdir(parents=True, exist_ok=True)
        (out_base / "error_rate").mkdir(parents=True, exist_ok=True)
        (out_base / "distance").mkdir(parents=True, exist_ok=True)

        # Threshold sweep
        if 'Threshold' in self.df.columns:
            for val in sorted(self.df['Threshold'].dropna().unique()):
                fname = f"threshold_{str(val).replace('.', '_')}.png"
                self.plot_success_rate_comparison(
                    save_path=out_base / "threshold" / fname,
                    threshold=val, error_rate=None, distance=None,
                    show=False
                )

        # Error-Rate sweep
        if 'Error-Rate' in self.df.columns:
            for val in sorted(self.df['Error-Rate'].dropna().unique()):
                fname = f"error_rate_{str(val).replace('.', '_')}.png"
                self.plot_success_rate_comparison(
                    save_path=out_base / "error_rate" / fname,
                    threshold=None, error_rate=val, distance=None,
                    show=False
                )

        # Distance sweep
        if 'Distance' in self.df.columns:
            for val in sorted(self.df['Distance'].dropna().unique()):
                fname = f"distance_{str(val).replace('.', '_')}.png"
                self.plot_success_rate_comparison(
                    save_path=out_base / "distance" / fname,
                    threshold=None, error_rate=None, distance=val,
                    show=False
                )

        print(f"Export terminé: {out_base}")

    def print_summary_statistics(self):
        """Print comprehensive summary statistics"""
        if self.df is None or len(self.df) == 0:
            print("No data available for summary statistics")
            return
        
        print("\n" + "="*60)
        print("ANALYSE COMPARATIVE ILP - STATISTIQUES RÉSUMÉES")
        print("="*60)
        
        print(f"\nNombre total d'expériences: {len(self.df)}")
        print(f"Paramètres uniques:")
        print(f"  - Haplotypes: {sorted(self.df['Haplotype'].unique())}")
        print(f"  - Error Rates: {sorted(self.df['Error-Rate'].unique())}")
        print(f"  - Thresholds: {sorted(self.df['Threshold'].unique())}")
        print(f"  - Distances: {sorted(self.df['Distance'].unique())}")
        
        print(f"\nPERFORMANCE GLOBALE:")
        print(f"  Taux de succès Pre-ILP:  {self.df['Success_Rate_Pre'].mean():.3f}")
        print(f"  Taux de succès Post-ILP: {self.df['Success_Rate_ILP'].mean():.3f}")
        improvement = self.df['Success_Rate_ILP'].mean() - self.df['Success_Rate_Pre'].mean()
        print(f"  Amélioration absolue:    {improvement:.3f}")
        
        print(f"\nRÉDUCTION DES ORPHELINS:")
        print(f"  Moyenne: {self.df['Orphan_Reduction'].mean():.2f}")
        print(f"  Médiane: {self.df['Orphan_Reduction'].median():.2f}")
        print(f"  Max: {self.df['Orphan_Reduction'].max():.0f}")
        
        print(f"\nRÉDUCTION DES ÉTAPES:")
        print(f"  Moyenne: {self.df['Steps_Reduction'].mean():.2f}")
        print(f"  Médiane: {self.df['Steps_Reduction'].median():.2f}")
        print(f"  Max: {self.df['Steps_Reduction'].max():.0f}")
        
        # Performance by haplotype count
        print(f"\nPERFORMANCE PAR NOMBRE D'HAPLOTYPES:")
        for hap_count in sorted(self.df['Haplotype'].unique()):
            hap_data = self.df[self.df['Haplotype'] == hap_count]
            pre_success = hap_data['Success_Rate_Pre'].mean()
            post_success = hap_data['Success_Rate_ILP'].mean()
            print(f"  {hap_count} haplotypes: {pre_success:.3f} → {post_success:.3f} (+{post_success-pre_success:.3f})")
    
    def generate_all_plots(self, output_dir="/home/mafoin/stage/experiment/analyze/plots_ilp"):
        """Generate only the success rate comparison plot"""
        if self.df is None or len(self.df) == 0:
            print("No valid data to analyze")
            return

        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)

        print("Génération du graphique de succès (ILP vs Pre-ILP)...")

        # (Optionnel) vous pouvez passer des filtres ici: threshold=..., error_rate=..., distance=...
        self.plot_success_rate_comparison(output_path / "success_rate_comparison.png")
        self.export_all_success_rate_by_params(output_dir=output_dir)
        print(f"Graphique sauvegardé dans: {output_path}")

def main():
    """Main function to run ILP analysis"""
    # Default CSV file path - update this to your ILP results file
    csv_file = "/home/mafoin/stage/experiment/efficience_ALL/results.csv"
    
    # Initialize analyzer
    analyzer = ILPEfficiencyAnalyzer(csv_file)
    
    if analyzer.df is not None:
        # Generate all plots and analyses
        analyzer.generate_all_plots()
    else:
        print(f"Impossible de charger les données depuis {csv_file}")
        print("Vérifiez que le fichier existe et a le bon format.")

if __name__ == "__main__":
    main()
