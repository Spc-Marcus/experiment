import pandas as pd

def extract_filepaths_with_ilp_calls(csv_file, output_file):
    """
    Lit un CSV et extrait les filepath où ilp_calls_total > 0
    
    Args:
        csv_file (str): Chemin vers le fichier CSV d'entrée
        output_file (str): Chemin vers le fichier texte de sortie
    """
    try:
        # Lire le fichier CSV
        df = pd.read_csv(csv_file)
        
        # Filtrer les lignes où ilp_calls_total > 0
        filtered_df = df[df['ilp_calls_total'] > 0]
        
        # Extraire les filepath
        filepaths = filtered_df['filepath'].tolist()
        
        # Écrire dans le fichier texte
        with open(output_file, 'w') as f:
            for filepath in filepaths:
                f.write(filepath + '\n')
        
        print(f"✓ {len(filepaths)} filepath(s) trouvé(s) avec ilp_calls_total > 0")
        print(f"✓ Résultats sauvegardés dans: {output_file}")
        
    except FileNotFoundError:
        print(f"Erreur: Le fichier {csv_file} n'existe pas")
    except KeyError as e:
        print(f"Erreur: Colonne manquante dans le CSV: {e}")
    except Exception as e:
        print(f"Erreur inattendue: {e}")

# Utilisation
if __name__ == "__main__":
    csv_input = "ILPV2/experiment_results_gurobi.csv"  # Remplacez par le nom de votre fichier CSV
    txt_output = "filepaths_with_ilp_calls.txt"
    
    extract_filepaths_with_ilp_calls(csv_input, txt_output)