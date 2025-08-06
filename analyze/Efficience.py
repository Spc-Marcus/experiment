import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
import ast
warnings.filterwarnings('ignore')

class EfficiencyAnalyzer:
    def __init__(self, csv_file_path=None):
        """Initialize the analyzer with CSV data if provided"""
        self.csv_file_path = csv_file_path
        self.df = None
        if csv_file_path and Path(csv_file_path).exists():
            self.load_data()
    
    def load_data(self):
        """Load data from CSV file with proper data type specification"""
        print(f"Loading data from {self.csv_file_path}")
        
        # First, let's examine the actual CSV structure
        with open(self.csv_file_path, 'r') as f:
            lines = f.readlines()[:10]  # Read first 10 lines
            print("First few lines of CSV:")
            for i, line in enumerate(lines):
                print(f"Line {i}: {line.strip()}")
        
        # The issue is that Matrix-Size is not properly quoted, causing column shift
        # Let's try a different approach: read as text and manually parse
        try:
            # Read the CSV line by line and manually parse
            data_rows = []
            with open(self.csv_file_path, 'r') as f:
                header = f.readline().strip()
                print(f"Header: {header}")
                
                for line_num, line in enumerate(f):
                    if line_num >= 10:  # Limit for debugging
                        break
                    parts = line.strip().split(',')
                    print(f"Line {line_num}: {len(parts)} parts: {parts[:8]}...")
            
            # It seems the Matrix-Size column is causing issues
            # Let's try reading with a different approach
            
            # Read the CSV and handle the Matrix-Size column specially
            import io
            
            # Read the file and fix the Matrix-Size formatting
            with open(self.csv_file_path, 'r') as f:
                content = f.read()
            
            # Fix the Matrix-Size column by ensuring it's properly quoted
            # Pattern: replace ,(digits, digits), with ,"(digits, digits)",
            import re
            fixed_content = re.sub(r',\((\d+), (\d+)\),', r',"(\1, \2)",', content)
            
            # Create a StringIO object to read the fixed content
            fixed_csv = io.StringIO(fixed_content)
            
            # Now try to read with proper dtypes
            dtype_dict = {
                'Threshold': 'float64',
                'Error-Rate': 'float64',
                'Strip': 'int64',
                'Haplotype': 'int64',
                'Steps-Count': 'int64',
                'Unused-Cols': 'int64',
                'Final-Clusters': 'int64',
                'Orphan-Reads': 'int64',
                'Time-Pre-Processing': 'float64',
                'Time-Post-Processing': 'float64'
            }
            
            self.df = pd.read_csv(fixed_csv, dtype=dtype_dict)
            print(f"Successfully loaded {len(self.df)} records after fixing CSV format")
            
            # Handle Matrix-Size column (convert string representation of tuple to actual tuple)
            if 'Matrix-Size' in self.df.columns:
                self.df['Matrix-Size'] = self.df['Matrix-Size'].apply(ast.literal_eval)
            
            print("Data types after correction:")
            print(self.df.dtypes)
            print("\nSample data:")
            print(self.df.head())
            
            print(f"Haplotype range: {self.df['Haplotype'].min()} to {self.df['Haplotype'].max()}")
            print(f"Unique haplotype counts: {sorted(self.df['Haplotype'].unique())}")
            print(f"Error-Rate range: {self.df['Error-Rate'].min()} to {self.df['Error-Rate'].max()}")
            print(f"Threshold range: {self.df['Threshold'].min()} to {self.df['Threshold'].max()}")
            
        except Exception as e:
            print(f"Error with regex fix approach: {e}")
            print("Trying manual column-by-column parsing...")
            
            # Manual parsing approach
            try:
                data_rows = []
                with open(self.csv_file_path, 'r') as f:
                    header_line = f.readline().strip()
                    expected_columns = header_line.split(',')
                    print(f"Expected columns: {expected_columns}")
                    
                    for line_num, line in enumerate(f):
                        line = line.strip()
                        if not line:
                            continue
                            
                        # Split by comma, but handle the Matrix-Size tuple
                        parts = []
                        in_tuple = False
                        current_part = ""
                        
                        for char in line:
                            if char == '(' and not in_tuple:
                                in_tuple = True
                                current_part += char
                            elif char == ')' and in_tuple:
                                in_tuple = False
                                current_part += char
                            elif char == ',' and not in_tuple:
                                parts.append(current_part)
                                current_part = ""
                            else:
                                current_part += char
                        
                        # Add the last part
                        if current_part:
                            parts.append(current_part)
                        
                        if len(parts) == len(expected_columns):
                            data_rows.append(parts)
                        else:
                            print(f"Skipping malformed line {line_num}: expected {len(expected_columns)} columns, got {len(parts)}")
                            if line_num < 5:  # Show first few problematic lines
                                print(f"  Line: {line}")
                                print(f"  Parts: {parts}")
                
                # Create DataFrame from parsed data
                if data_rows:
                    self.df = pd.DataFrame(data_rows, columns=expected_columns)
                    
                    # Convert data types
                    self.df['Threshold'] = pd.to_numeric(self.df['Threshold'], errors='coerce')
                    self.df['Error-Rate'] = pd.to_numeric(self.df['Error-Rate'], errors='coerce')
                    self.df['Strip'] = pd.to_numeric(self.df['Strip'], errors='coerce')
                    self.df['Haplotype'] = pd.to_numeric(self.df['Haplotype'], errors='coerce')
                    self.df['Steps-Count'] = pd.to_numeric(self.df['Steps-Count'], errors='coerce')
                    self.df['Unused-Cols'] = pd.to_numeric(self.df['Unused-Cols'], errors='coerce')
                    self.df['Final-Clusters'] = pd.to_numeric(self.df['Final-Clusters'], errors='coerce')
                    self.df['Orphan-Reads'] = pd.to_numeric(self.df['Orphan-Reads'], errors='coerce')
                    self.df['Time-Pre-Processing'] = pd.to_numeric(self.df['Time-Pre-Processing'], errors='coerce')
                    self.df['Time-Post-Processing'] = pd.to_numeric(self.df['Time-Post-Processing'], errors='coerce')
                    
                    # Remove rows with any NaN values
                    before_count = len(self.df)
                    self.df = self.df.dropna()
                    after_count = len(self.df)
                    
                    if before_count != after_count:
                        print(f"Removed {before_count - after_count} rows with invalid values")
                    
                    # Convert to proper types
                    self.df['Haplotype'] = self.df['Haplotype'].astype(int)
                    self.df['Strip'] = self.df['Strip'].astype(int)
                    self.df['Steps-Count'] = self.df['Steps-Count'].astype(int)
                    self.df['Unused-Cols'] = self.df['Unused-Cols'].astype(int)
                    self.df['Final-Clusters'] = self.df['Final-Clusters'].astype(int)
                    self.df['Orphan-Reads'] = self.df['Orphan-Reads'].astype(int)
                    
                    print(f"Successfully parsed {len(self.df)} valid records")
                    print("Final data types:")
                    print(self.df.dtypes)
                    print("\nSample data:")
                    print(self.df.head())
                    
                    if len(self.df) > 0:
                        print(f"Haplotype range: {self.df['Haplotype'].min()} to {self.df['Haplotype'].max()}")
                        print(f"Unique haplotype counts: {sorted(self.df['Haplotype'].unique())}")
                        print(f"Error-Rate range: {self.df['Error-Rate'].min()} to {self.df['Error-Rate'].max()}")
                        print(f"Threshold range: {self.df['Threshold'].min()} to {self.df['Threshold'].max()}")
                else:
                    print("No valid data rows found")
                    
            except Exception as manual_error:
                print(f"Manual parsing also failed: {manual_error}")
                self.df = None

    def plot_certainty_error_heatmap(self, save_path=None):
        """Generate heatmap of success rate (clusters = haplotypes) for certainty vs error rate"""
        if self.df is None or len(self.df) == 0:
            print("No data loaded or no valid data")
            return
        
        # Calculate success rate (Final-Clusters == Haplotype)
        self.df['Success_Rate'] = (self.df['Final-Clusters'] == self.df['Haplotype']).astype(int)
        
        # Get unique haplotype counts and limit to reasonable range
        haplotype_counts = sorted(self.df['Haplotype'].unique())
        
        if len(haplotype_counts) == 0:
            print("No valid haplotype data found")
            return
        
        # Limit to first 7 haplotype counts for visualization
        if len(haplotype_counts) > 7:
            haplotype_counts = haplotype_counts[:7]
            print(f"Limiting to first 7 haplotype counts: {haplotype_counts}")
        
        print(f"Processing haplotype counts: {haplotype_counts}")

        # Calculate layout for subplots
        n_haplotypes = len(haplotype_counts)
        if n_haplotypes <= 3:
            nrows, ncols = 1, n_haplotypes
            figsize = (6*n_haplotypes, 5)
        elif n_haplotypes <= 6:
            nrows, ncols = 2, 3
            figsize = (18, 10)
        else:
            ncols = int(np.ceil(np.sqrt(n_haplotypes)))
            nrows = int(np.ceil(n_haplotypes / ncols))
            figsize = (6*ncols, 5*nrows)
        
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize, constrained_layout=True)
        
        # Ensure axes is always a list
        if n_haplotypes == 1:
            axes = [axes]
        elif nrows > 1:
            axes = axes.flatten()
        
        fig.suptitle('Heatmap: Taux de Succès (Certitude vs Error Rate) par Haplotype', fontsize=16)
        
        for i, haplotype_count in enumerate(haplotype_counts):
            if i >= len(axes):
                break
            
            print(f"Processing {haplotype_count} haplotypes...")
                
            # Filter data for current haplotype count
            haplotype_data = self.df[self.df['Haplotype'] == haplotype_count]
            
            if len(haplotype_data) == 0:
                axes[i].text(0.5, 0.5, f'No data for\n{haplotype_count} haplotypes', 
                           ha='center', va='center', transform=axes[i].transAxes)
                axes[i].set_title(f'{haplotype_count} Haplotypes')
                continue
            
            print(f"  Found {len(haplotype_data)} records for {haplotype_count} haplotypes")
            
            # Create pivot table for heatmap
            try:
                pivot_data = haplotype_data.groupby(['Threshold', 'Error-Rate'])['Success_Rate'].mean().unstack(fill_value=np.nan)
                
                if pivot_data.empty:
                    axes[i].text(0.5, 0.5, f'No valid data for\n{haplotype_count} haplotypes', 
                               ha='center', va='center', transform=axes[i].transAxes)
                    axes[i].set_title(f'{haplotype_count} Haplotypes')
                    continue
                
                # Create heatmap
                im = sns.heatmap(pivot_data, 
                               ax=axes[i], 
                               cmap='RdYlGn', 
                               vmin=0, 
                               vmax=1,
                               annot=True, 
                               fmt='.2f',
                               cbar_kws={'label': 'Taux de Succès'})
                
                axes[i].set_title(f'{haplotype_count} Haplotypes')
                axes[i].set_xlabel('Error Rate')
                axes[i].set_ylabel('Certitude (Threshold)')
                
                # Rotate x-axis labels for better readability
                axes[i].tick_params(axis='x', rotation=45)
                axes[i].tick_params(axis='y', rotation=0)
                
            except Exception as e:
                print(f"Error creating heatmap for {haplotype_count} haplotypes: {e}")
                axes[i].text(0.5, 0.5, f'Error creating heatmap\nfor {haplotype_count} haplotypes', 
                           ha='center', va='center', transform=axes[i].transAxes)
                axes[i].set_title(f'{haplotype_count} Haplotypes')
        
        # Hide unused subplots
        if n_haplotypes < len(axes):
            for j in range(n_haplotypes, len(axes)):
                axes[j].set_visible(False)
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.show()

    
    def generate_all_plots(self, output_dir="/home/mafoin/stage/experiment/analyze/plots"):
        """Generate all analysis plots"""
        if self.df is None or len(self.df) == 0:
            print("No valid data to analyze")
            return
            
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        print("Génération de tous les graphiques d'analyse...")
        
        # Generate each plot
        self.plot_certainty_error_heatmap(output_path / "heatmap_certitude_error_by_haplotype.png")

        
        print(f"Tous les graphiques sauvegardés dans: {output_path}")

def main():
    """Main function to run all analyses"""
    # Default CSV file path
    csv_file = "/home/mafoin/stage/experiment/efficience_CH/results.csv"
    
    # Initialize analyzer
    analyzer = EfficiencyAnalyzer(csv_file)
    
    if analyzer.df is not None:
        # Generate all plots and analyses
        analyzer.generate_all_plots()
    else:
        print(f"Impossible de charger les données depuis {csv_file}")
        print("Vérifiez que le fichier existe et a le bon format.")

if __name__ == "__main__":
    main()
