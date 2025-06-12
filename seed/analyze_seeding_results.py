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

# Configuration des graphiques
plt.style.use('default')
sns.set_palette("Set2")
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 11

class SeedingAnalyzer:
    """Analyseur pour les résultats d'expérimentation de seeding."""
    
    def __init__(self, csv_file: str, output_dir: str = "plots"):
        """
        Initialise l'analyseur.
        
        Parameters
        ----------
        csv_file : str
            Chemin vers le fichier CSV des résultats
        output_dir : str
            Répertoire de sortie pour les graphiques
        """
        self.csv_file = csv_file
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Charger les données
        self.df = pd.read_csv(csv_file)
        print(f"Données chargées: {len(self.df)} expériences")
        
        # Nettoyer et préparer les données
        self._prepare_data()
        
    def _prepare_data(self):
        """Prépare et nettoie les données."""
        # Calculer des métriques dérivées
        self.df['matrix_size'] = self.df['matrix_m'] * self.df['matrix_n']
        self.df['ilp_efficiency'] = self.df['patterns_found'] / np.maximum(self.df['nb_ilp_calculated'], 1)
        
        print(f"Données préparées: {len(self.df)} expériences valides")
        print(f"Plage X_factor: {self.df['X_factor'].min()} - {self.df['X_factor'].max()}")
        print(f"Plage step_n: {self.df['step_n'].min()} - {self.df['step_n'].max()}")
        
    def generate_summary_stats(self):
        """Génère des statistiques descriptives."""
        print("\n" + "="*60)
        print("STATISTIQUES DESCRIPTIVES")
        print("="*60)
        
        # Statistiques générales
        print(f"\nNombre total d'expériences: {len(self.df)}")
        print(f"Nombre de fichiers uniques: {self.df['read_name'].nunique()}")
        print(f"Nombre d'haplotypes: {sorted(self.df['nb_haplotypes'].unique())}")
        
        # Paramètres de seeding testés
        print(f"\nParamètres X_factor testés: {sorted(self.df['X_factor'].unique())}")
        print(f"Paramètres step_n testés: {sorted(self.df['step_n'].unique())}")
        
        # Métriques de performance clés
        print(f"\nMétriques de performance:")
        print(f"Temps total moyen: {self.df['total_time'].mean():.3f}s (±{self.df['total_time'].std():.3f})")
        print(f"Appels ILP moyens: {self.df['nb_ilp_calculated'].mean():.1f} (±{self.df['nb_ilp_calculated'].std():.1f})")
        print(f"Patterns trouvés moyens: {self.df['patterns_found'].mean():.1f} (±{self.df['patterns_found'].std():.1f})")
        print(f"Densité de matrice moyenne: {self.df['matrix_density'].mean():.3f}")
        
        # Utilisation du seeding
        if 'seeding_used' in self.df.columns:
            seeding_used = self.df['seeding_used'].sum()
            print(f"\nExpériences utilisant le seeding: {seeding_used}/{len(self.df)} ({seeding_used/len(self.df)*100:.1f}%)")
            
        # Top configurations
        print(f"\n" + "="*40)
        print("TOP 3 CONFIGURATIONS - TEMPS TOTAL")
        print("="*40)
        best_time = self.df.groupby(['X_factor', 'step_n'])['total_time'].mean().sort_values().head(3)
        for i, ((x, s), time_val) in enumerate(best_time.items()):
            count = len(self.df[(self.df['X_factor'] == x) & (self.df['step_n'] == s)])
            print(f"{i+1}. X_factor={x}, step_n={s}: {time_val:.3f}s (n={count})")
            
        print(f"\n" + "="*40)
        print("TOP 3 CONFIGURATIONS - EFFICACITÉ ILP")
        print("="*40)
        best_eff = self.df.groupby(['X_factor', 'step_n'])['ilp_efficiency'].mean().sort_values(ascending=False).head(3)
        for i, ((x, s), eff_val) in enumerate(best_eff.items()):
            count = len(self.df[(self.df['X_factor'] == x) & (self.df['step_n'] == s)])
            print(f"{i+1}. X_factor={x}, step_n={s}: {eff_val:.3f} patterns/ILP (n={count})")
    
    def plot_optimization_landscape(self):
        """Graphique 1: Paysage d'optimisation avec heatmap du temps et contours d'efficacité."""
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))
        
        # Préparer les données pivot
        pivot_time = self.df.groupby(['X_factor', 'step_n'])['total_time'].mean().unstack()
        pivot_eff = self.df.groupby(['X_factor', 'step_n'])['ilp_efficiency'].mean().unstack()
        
        # Heatmap du temps total
        im = ax.imshow(pivot_time.values, cmap='RdYlBu_r', aspect='auto', alpha=0.8)
        
        # Ajouter les valeurs de temps dans la heatmap
        for i in range(len(pivot_time.index)):
            for j in range(len(pivot_time.columns)):
                if not np.isnan(pivot_time.iloc[i, j]):
                    text = ax.text(j, i, f'{pivot_time.iloc[i, j]:.2f}s',
                                  ha="center", va="center", color="black", fontsize=10, fontweight='bold')
        
        # Ajouter les contours d'efficacité
        if not pivot_eff.empty and not pivot_eff.isna().all().all():
            X, Y = np.meshgrid(range(len(pivot_eff.columns)), range(len(pivot_eff.index)))
            contours = ax.contour(X, Y, pivot_eff.values, levels=5, colors='white', alpha=0.7, linewidths=2)
            ax.clabel(contours, inline=True, fontsize=9, fmt='%.2f eff')
        
        # Configuration des axes
        ax.set_title('Paysage d\'Optimisation des Paramètres de Seeding\n(Couleur=Temps, Contours=Efficacité ILP)', 
                    fontsize=14, fontweight='bold', pad=20)
        ax.set_xlabel('step_n', fontsize=12)
        ax.set_ylabel('X_factor', fontsize=12)
        ax.set_xticks(range(len(pivot_time.columns)))
        ax.set_xticklabels(pivot_time.columns)
        ax.set_yticks(range(len(pivot_time.index)))
        ax.set_yticklabels(pivot_time.index)
        
        # Colorbar pour le temps
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Temps Total Moyen (secondes)', fontsize=12)
        
        # Identifier et marquer les points optimaux
        best_time_idx = np.unravel_index(np.nanargmin(pivot_time.values), pivot_time.shape)
        best_eff_idx = np.unravel_index(np.nanargmax(pivot_eff.values), pivot_eff.shape)
        
        # Marquer le meilleur temps
        ax.scatter(best_time_idx[1], best_time_idx[0], marker='*', s=200, 
                  c='red', edgecolors='white', linewidth=2, label='Optimal Temps')
        
        # Marquer la meilleure efficacité
        ax.scatter(best_eff_idx[1], best_eff_idx[0], marker='D', s=150, 
                  c='lime', edgecolors='black', linewidth=2, label='Optimal Efficacité')
        
        ax.legend(loc='upper right', bbox_to_anchor=(1, 1))
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'optimization_landscape.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"✅ Paysage d'optimisation sauvegardé: {self.output_dir / 'optimization_landscape.png'}")
    
    def plot_performance_vs_complexity(self):
        """Graphique 2: Performance vs Complexité avec analyse multi-dimensionnelle."""
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('Analyse Performance vs Complexité des Matrices', fontsize=16, fontweight='bold')
        
        # 1. Temps vs Taille de matrice (coloré par nombre d'haplotypes, taille par efficacité)
        scatter1 = ax1.scatter(self.df['matrix_size'], self.df['total_time'], 
                              c=self.df['nb_haplotypes'], s=self.df['ilp_efficiency']*50, 
                              cmap='viridis', alpha=0.7, edgecolors='black', linewidth=0.5)
        ax1.set_xlabel('Taille de la matrice (m×n)')
        ax1.set_ylabel('Temps total (s)')
        ax1.set_title('Temps vs Taille\n(Couleur=Haplotypes, Taille=Efficacité)', fontweight='bold')
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.grid(True, alpha=0.3)
        
        # Ligne de tendance
        z1 = np.polyfit(np.log(self.df['matrix_size']), np.log(self.df['total_time']), 1)
        x_trend = np.logspace(np.log10(self.df['matrix_size'].min()), 
                             np.log10(self.df['matrix_size'].max()), 100)
        y_trend = np.exp(z1[1]) * x_trend ** z1[0]
        ax1.plot(x_trend, y_trend, 'r--', alpha=0.8, linewidth=2, 
                label=f'Tendance (pente={z1[0]:.2f})')
        ax1.legend()
        
        cbar1 = plt.colorbar(scatter1, ax=ax1)
        cbar1.set_label('Nb Haplotypes')
        
        # 2. ILP vs Densité (coloré par X_factor, taille par step_n)
        scatter2 = ax2.scatter(self.df['matrix_density'], self.df['nb_ilp_calculated'], 
                              c=self.df['X_factor'], s=self.df['step_n']*3, 
                              cmap='plasma', alpha=0.7, edgecolors='black', linewidth=0.5)
        ax2.set_xlabel('Densité de la matrice')
        ax2.set_ylabel('Nombre d\'appels ILP')
        ax2.set_title('ILP vs Densité\n(Couleur=X_factor, Taille=step_n)', fontweight='bold')
        ax2.grid(True, alpha=0.3)
        
        cbar2 = plt.colorbar(scatter2, ax=ax2)
        cbar2.set_label('X_factor')
        
        # 3. Efficacité par configuration (violin plot)
        # Créer une variable combinée pour les configurations
        self.df['config'] = self.df['X_factor'].astype(str) + '_' + self.df['step_n'].astype(str)
        top_configs = self.df.groupby('config')['total_time'].mean().sort_values().head(8).index
        df_top = self.df[self.df['config'].isin(top_configs)]
        
        # Violin plot de l'efficacité
        parts = ax3.violinplot([df_top[df_top['config'] == config]['ilp_efficiency'].values 
                               for config in top_configs], 
                              positions=range(len(top_configs)), widths=0.6, showmeans=True)
        
        for pc in parts['bodies']:
            pc.set_facecolor('lightblue')
            pc.set_alpha(0.7)
        
        ax3.set_xlabel('Configuration (X_factor_step_n)')
        ax3.set_ylabel('Efficacité ILP (patterns/ILP)')
        ax3.set_title('Distribution de l\'Efficacité\n(Top 8 configurations)', fontweight='bold')
        ax3.set_xticks(range(len(top_configs)))
        ax3.set_xticklabels(top_configs, rotation=45)
        ax3.grid(True, alpha=0.3)
        
        # 4. Analyse comparative des paramètres
        param_analysis = []
        for x in sorted(self.df['X_factor'].unique()):
            for s in sorted(self.df['step_n'].unique()):
                subset = self.df[(self.df['X_factor'] == x) & (self.df['step_n'] == s)]
                if len(subset) > 0:
                    param_analysis.append({
                        'X_factor': x,
                        'step_n': s,
                        'avg_time': subset['total_time'].mean(),
                        'avg_ilp': subset['nb_ilp_calculated'].mean(),
                        'avg_eff': subset['ilp_efficiency'].mean(),
                        'count': len(subset)
                    })
        
        param_df = pd.DataFrame(param_analysis)
        
        # Bubble chart: X_factor vs step_n, couleur=temps, taille=efficacité
        bubble = ax4.scatter(param_df['X_factor'], param_df['step_n'], 
                           c=param_df['avg_time'], s=param_df['avg_eff']*100, 
                           cmap='RdYlBu_r', alpha=0.8, edgecolors='black', linewidth=1)
        
        # Ajouter les valeurs de temps sur chaque bulle
        for _, row in param_df.iterrows():
            ax4.text(row['X_factor'], row['step_n'], f'{row["avg_time"]:.2f}s',
                    ha='center', va='center', fontsize=8, fontweight='bold')
        
        ax4.set_xlabel('X_factor')
        ax4.set_ylabel('step_n')
        ax4.set_title('Vue d\'Ensemble des Paramètres\n(Couleur=Temps, Taille=Efficacité)', fontweight='bold')
        ax4.grid(True, alpha=0.3)
        
        cbar4 = plt.colorbar(bubble, ax=ax4)
        cbar4.set_label('Temps Moyen (s)')
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'performance_vs_complexity.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"✅ Analyse performance vs complexité sauvegardée: {self.output_dir / 'performance_vs_complexity.png'}")
    
    def generate_performance_report(self):
        """Génère un rapport de performance concis mais complet."""
        report_file = self.output_dir / 'performance_report.txt'
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("="*80 + "\n")
            f.write("RAPPORT D'ANALYSE DES EXPÉRIENCES DE SEEDING\n")
            f.write("="*80 + "\n\n")
            
            # Informations générales
            f.write("INFORMATIONS GÉNÉRALES\n")
            f.write("-"*30 + "\n")
            f.write(f"Fichier source: {self.csv_file}\n")
            f.write(f"Nombre d'expériences: {len(self.df)}\n")
            f.write(f"Fichiers uniques: {self.df['read_name'].nunique()}\n")
            f.write(f"Haplotypes: {sorted(self.df['nb_haplotypes'].unique())}\n")
            f.write(f"X_factor: {sorted(self.df['X_factor'].unique())}\n")
            f.write(f"step_n: {sorted(self.df['step_n'].unique())}\n\n")
            
            # Métriques globales
            f.write("MÉTRIQUES GLOBALES\n")
            f.write("-"*20 + "\n")
            f.write(f"Temps total moyen: {self.df['total_time'].mean():.3f}s (±{self.df['total_time'].std():.3f})\n")
            f.write(f"Appels ILP moyens: {self.df['nb_ilp_calculated'].mean():.1f} (±{self.df['nb_ilp_calculated'].std():.1f})\n")
            f.write(f"Efficacité ILP moyenne: {self.df['ilp_efficiency'].mean():.3f} (±{self.df['ilp_efficiency'].std():.3f})\n")
            f.write(f"Patterns moyens: {self.df['patterns_found'].mean():.1f} (±{self.df['patterns_found'].std():.1f})\n\n")
            
            # Configurations optimales
            f.write("CONFIGURATIONS OPTIMALES\n")
            f.write("-"*25 + "\n")
            
            optimal_time = self.df.groupby(['X_factor', 'step_n'])['total_time'].mean().idxmin()
            optimal_eff = self.df.groupby(['X_factor', 'step_n'])['ilp_efficiency'].mean().idxmax()
            
            f.write(f"⭐ Optimal TEMPS: X_factor={optimal_time[0]}, step_n={optimal_time[1]}\n")
            f.write(f"⭐ Optimal EFFICACITÉ: X_factor={optimal_eff[0]}, step_n={optimal_eff[1]}\n\n")
            
            # Top 3 configurations
            f.write("TOP 3 CONFIGURATIONS\n")
            f.write("-"*20 + "\n")
            best_overall = self.df.groupby(['X_factor', 'step_n']).agg({
                'total_time': 'mean',
                'ilp_efficiency': 'mean'
            }).round(3)
            
            # Score combiné (temps normalisé inversé + efficacité normalisée)
            time_normalized = 1 - (best_overall['total_time'] - best_overall['total_time'].min()) / (best_overall['total_time'].max() - best_overall['total_time'].min())
            eff_normalized = (best_overall['ilp_efficiency'] - best_overall['ilp_efficiency'].min()) / (best_overall['ilp_efficiency'].max() - best_overall['ilp_efficiency'].min())
            combined_score = (time_normalized + eff_normalized) / 2
            
            top3 = combined_score.sort_values(ascending=False).head(3)
            for i, ((x, s), score) in enumerate(top3.items()):
                time_val = best_overall.loc[(x, s), 'total_time']
                eff_val = best_overall.loc[(x, s), 'ilp_efficiency']
                count = len(self.df[(self.df['X_factor'] == x) & (self.df['step_n'] == s)])
                f.write(f"{i+1}. X_factor={x}, step_n={s}: {time_val:.3f}s, {eff_val:.3f}eff (score={score:.3f}, n={count})\n")
            
        print(f"✅ Rapport de performance sauvegardé: {report_file}")
    
    def run_full_analysis(self):
        """Exécute l'analyse complète avec seulement 2 graphiques essentiels."""
        print("🚀 Démarrage de l'analyse des résultats de seeding...")
        
        # Statistiques descriptives
        self.generate_summary_stats()
        
        # Génération des 2 graphiques essentiels
        print(f"\n📊 Génération des visualisations dans {self.output_dir}/...")
        
        self.plot_optimization_landscape()
        self.plot_performance_vs_complexity()
        
        # Rapport final
        self.generate_performance_report()
        
        print(f"\n✅ Analyse complète terminée!")
        print(f"📁 Fichiers sauvegardés dans: {self.output_dir}/")
        print(f"📋 Rapport détaillé: {self.output_dir}/performance_report.txt")
        
        # Résumé des principales conclusions
        print(f"\n🎯 CONCLUSIONS PRINCIPALES:")
        optimal_time = self.df.groupby(['X_factor', 'step_n'])['total_time'].mean().idxmin()
        optimal_eff = self.df.groupby(['X_factor', 'step_n'])['ilp_efficiency'].mean().idxmax()
        print(f"  • Meilleure config temps: X_factor={optimal_time[0]}, step_n={optimal_time[1]}")
        print(f"  • Meilleure config efficacité: X_factor={optimal_eff[0]}, step_n={optimal_eff[1]}")
        print(f"  • 2 graphiques générés: paysage d'optimisation + analyse multi-dimensionnelle")

def main():
    """Fonction principale."""
    parser = argparse.ArgumentParser(description='Analyse des résultats d\'expérimentation de seeding')
    parser.add_argument('csv_file', help='Fichier CSV des résultats')
    parser.add_argument('--output-dir', '-o', default='seeding_analysis_plots', 
                       help='Répertoire de sortie pour les graphiques')
    parser.add_argument('--summary-only', '-s', action='store_true', 
                       help='Affiche seulement les statistiques, sans graphiques')
    
    args = parser.parse_args()
    
    # Vérifier que le fichier existe
    if not Path(args.csv_file).exists():
        print(f"❌ Erreur: Le fichier {args.csv_file} n'existe pas")
        sys.exit(1)
    
    # Créer l'analyseur
    analyzer = SeedingAnalyzer(args.csv_file, args.output_dir)
    
    if args.summary_only:
        analyzer.generate_summary_stats()
    else:
        analyzer.run_full_analysis()

if __name__ == "__main__":
    main()
