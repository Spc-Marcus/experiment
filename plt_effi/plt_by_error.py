import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

csv_files = glob.glob(os.path.join(os.path.dirname(__file__), '../data/best_by_error_rate_*.csv'))

plt.figure(figsize=(10, 6))

for csv_path in csv_files:
    base = os.path.basename(csv_path)
    label = base.replace('best_by_error_rate_', '').replace('.csv', '')

    # Essayer de lire avec différents séparateurs
    try:
        df = pd.read_csv(csv_path, sep='\t')
    except Exception:
        df = pd.read_csv(csv_path, sep=',')

    # Afficher les colonnes pour debug
    print(f"{csv_path} colonnes: {df.columns.tolist()}")

    # Nettoyer les noms de colonnes
    df.columns = [col.strip() for col in df.columns]

    # Trouver les colonnes correspondantes (insensible à la casse et aux tirets)
    col_error = next((c for c in df.columns if c.lower().replace('-', '').replace('_', '') == 'errorrate'), None)
    col_success = next((c for c in df.columns if c.lower().replace('-', '').replace('_', '') == 'meansuccess'), None)

    if col_error and col_success:
        plt.plot(df[col_error], df[col_success], marker='o', label=label)
    else:
        print(f"Colonnes non trouvées dans {csv_path}")

plt.xlabel('Error Rate')
plt.ylabel('Mean Success')
plt.title('Succès moyen en fonction du taux d\'erreur')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()