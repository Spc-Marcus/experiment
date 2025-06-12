"""Script rapide pour analyser les résultats de seeding sans interface graphique."""

import pandas as pd
import numpy as np
import sys
from pathlib import Path

def quick_analysis(csv_file: str):
    """Analyse rapide des données de seeding."""
    
    if not Path(csv_file).exists():
        print(f"❌ Fichier {csv_file} introuvable")
        return
    
    # Charger les données
    df = pd.read_csv(csv_file)
    print(f"📊 Données chargées: {len(df)} expériences")
    
    # Statistiques de base
    print(f"\n{'='*50}")
    print("RÉSUMÉ RAPIDE")
    print(f"{'='*50}")
    
    print(f"Fichiers traités: {df['read_name'].nunique()}")
    print(f"Haplotypes: {sorted(df['nb_haplotypes'].unique())}")
    print(f"X_factor testé: {sorted(df['X_factor'].unique())}")
    print(f"step_n testé: {sorted(df['step_n'].unique())}")
    
    # Performance moyenne
    print(f"\n📈 PERFORMANCE MOYENNE:")
    print(f"Temps total: {df['total_time'].mean():.3f}s (±{df['total_time'].std():.3f})")
    print(f"Appels ILP: {df['nb_ilp_calculated'].mean():.1f} (±{df['nb_ilp_calculated'].std():.1f})")
    print(f"Patterns trouvés: {df['patterns_found'].mean():.1f} (±{df['patterns_found'].std():.1f})")
    
    # Efficacité calculée
    df['efficiency'] = df['patterns_found'] / np.maximum(df['nb_ilp_calculated'], 1)
    print(f"Efficacité ILP: {df['efficiency'].mean():.3f} patterns/ILP")
    
    # Meilleures configurations
    print(f"\n🏆 TOP 3 CONFIGURATIONS (par temps):")
    best_time = df.groupby(['X_factor', 'step_n'])['total_time'].mean().sort_values()
    for i, ((x, s), time_val) in enumerate(best_time.head(3).items()):
        count = len(df[(df['X_factor'] == x) & (df['step_n'] == s)])
        print(f"  {i+1}. X_factor={x}, step_n={s}: {time_val:.3f}s (n={count})")
    
    print(f"\n🎯 TOP 3 CONFIGURATIONS (par efficacité ILP):")
    best_eff = df.groupby(['X_factor', 'step_n'])['efficiency'].mean().sort_values(ascending=False)
    for i, ((x, s), eff_val) in enumerate(best_eff.head(3).items()):
        count = len(df[(df['X_factor'] == x) & (df['step_n'] == s)])
        print(f"  {i+1}. X_factor={x}, step_n={s}: {eff_val:.3f} patterns/ILP (n={count})")
    
    # Vérification de l'utilisation du seeding
    if 'seeding_used' in df.columns:
        seeding_used = df['seeding_used'].sum()
        print(f"\n⚙️  UTILISATION DU SEEDING:")
        print(f"Expériences avec seeding: {seeding_used}/{len(df)} ({seeding_used/len(df)*100:.1f}%)")
        if seeding_used == 0:
            print("⚠️  ATTENTION: Aucune expérience n'a utilisé les paramètres de seeding!")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python quick_analysis.py <fichier_csv>")
        sys.exit(1)
    
    quick_analysis(sys.argv[1])
