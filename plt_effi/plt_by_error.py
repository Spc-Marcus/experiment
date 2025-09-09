import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

# Chercher tous les fichiers CSV correspondants à la racine
csv_files = glob.glob(os.path.join(os.path.dirname(__file__), '../data/best_by_error_rate_*.csv'))

plt.figure(figsize=(10, 6))

for csv_path in csv_files:
    # Extraire le nom XXXX pour la légende
    base = os.path.basename(csv_path)
    label = base.replace('best_by_error_rate_', '').replace('.csv', '')

    # Charger le CSV
    df = pd.read_csv(csv_path, sep='\t')
    # Tracer la courbe
    plt.plot(df['Error-Rate'], df['Mean_Success'], marker='o', label=label)

plt.xlabel('Error Rate')
plt.ylabel('Mean Success')
plt.title('Succès moyen en fonction du taux d\'erreur')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()