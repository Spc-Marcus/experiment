import json
import pandas as pd
import numpy as np
import argparse
import re
from pathlib import Path

class MatrixParser:
    def __init__(self, results_path):
        """Initialize parser with clustering results file (JSON or text format)."""
        self.results_path = results_path
        
        # Detect file format and load data
        if self._is_json_format():
            with open(results_path, 'r') as f:
                self.clustering_data = json.load(f)
            self.format_type = 'json'
        else:
            self.lists_data = self._parse_text_format()
            self.format_type = 'text'
    
    def _is_json_format(self):
        """Check if the file is in JSON format."""
        try:
            with open(self.results_path, 'r') as f:
                json.load(f)
            return True
        except (json.JSONDecodeError, UnicodeDecodeError):
            return False
    
    def _parse_text_format(self):
        """Parse the simple text format with lists."""
        with open(self.results_path, 'r') as f:
            content = f.read()
        
        # Extract lists using regex
        list_pattern = r'\[([0-9, ]+)\]'
        matches = re.findall(list_pattern, content)
        
        lists = []
        for match in matches:
            # Convert string of numbers to list of integers
            numbers = [int(x.strip()) for x in match.split(',') if x.strip()]
            lists.append(numbers)
        
        return lists
    
    def extract_rows_cols_from_step(self, step_number=1):
        """Extract rows and columns from a specific clustering step."""
        if self.format_type != 'json':
            raise ValueError("Clustering steps are only available in JSON format")
            
        clustering_steps = self.clustering_data.get('clustering_steps', [])
        
        for step in clustering_steps:
            if step.get('step_number') == step_number:
                rows = step.get('reads_group1', []) + step.get('reads_group0', [])
                cols = step.get('columns', [])
                return sorted(rows), sorted(cols)
        
        raise ValueError(f"Step {step_number} not found in clustering results")
    
    def extract_rows_cols_from_biclique(self, phase_index=0):
        """Extract rows and columns from quasi-biclique results."""
        if self.format_type != 'json':
            raise ValueError("Biclique results are only available in JSON format")
            
        biclique_results = self.clustering_data.get('quasi_biclique_results', [])
        
        if phase_index >= len(biclique_results):
            raise ValueError(f"Phase index {phase_index} out of range")
        
        result = biclique_results[phase_index]
        rows = result.get('selected_rows', [])
        cols = result.get('selected_cols', [])
        
        return sorted(rows), sorted(cols)
    
    def extract_rows_cols_from_lists(self, list1_index=0, list2_index=1):
        """Extract rows and columns from text format lists."""
        if self.format_type != 'text':
            raise ValueError("List extraction is only available for text format")
        
        if list1_index >= len(self.lists_data) or list2_index >= len(self.lists_data):
            raise ValueError(f"List indices out of range. Available lists: {len(self.lists_data)}")
        
        rows = self.lists_data[list1_index]
        cols = self.lists_data[list2_index]
        
        return sorted(rows), sorted(cols)
    
    def get_available_options(self):
        """Get available extraction options based on format."""
        if self.format_type == 'json':
            options = []
            if hasattr(self, 'clustering_data'):
                steps = self.clustering_data.get('clustering_steps', [])
                if steps:
                    options.append(f"Clustering steps: {[s.get('step_number') for s in steps]}")
                
                bicliques = self.clustering_data.get('quasi_biclique_results', [])
                if bicliques:
                    options.append(f"Biclique results: {len(bicliques)} available (indices 0-{len(bicliques)-1})")
            return options
        else:
            return [f"Text format: {len(self.lists_data)} lists available (indices 0-{len(self.lists_data)-1})"]

    def filter_matrix(self, csv_path, rows, cols, output_path):
        """Filter CSV matrix by specified rows and columns."""
        # Read CSV matrix
        matrix = pd.read_csv(csv_path, header=None)
        
        # Validate indices
        max_row = matrix.shape[0] - 1
        max_col = matrix.shape[1] - 1
        
        valid_rows = [r for r in rows if 0 <= r <= max_row]
        valid_cols = [c for c in cols if 0 <= c <= max_col]
        
        if len(valid_rows) != len(rows):
            print(f"Warning: Some row indices out of range. Using {len(valid_rows)}/{len(rows)} rows.")
        if len(valid_cols) != len(cols):
            print(f"Warning: Some column indices out of range. Using {len(valid_cols)}/{len(cols)} columns.")
        
        # Extract submatrix
        filtered_matrix = matrix.iloc[valid_rows, valid_cols]
        
        # Save to CSV
        filtered_matrix.to_csv(output_path, header=False, index=False)
        
        print(f"Filtered matrix saved to {output_path}")
        print(f"Original shape: {matrix.shape}")
        print(f"Filtered shape: {filtered_matrix.shape}")
        print(f"Selected rows: {valid_rows}")
        print(f"Selected columns: {valid_cols}")
        
        return filtered_matrix

def main():
    parser = argparse.ArgumentParser(description='Extract matrix subsets based on clustering results')
    parser.add_argument('results_file', help='Path to results file (JSON or text format)')
    parser.add_argument('input_csv', help='Path to input CSV matrix')
    parser.add_argument('output_csv', help='Path to output filtered CSV matrix')
    parser.add_argument('--step', type=int, help='Clustering step number to use (JSON format only)')
    parser.add_argument('--biclique', type=int, help='Use biclique result at specified index (JSON format only)')
    parser.add_argument('--list1', type=int, default=0, help='Index of first list for rows (text format, default: 0)')
    parser.add_argument('--list2', type=int, default=1, help='Index of second list for cols (text format, default: 1)')
    parser.add_argument('--info', action='store_true', help='Show available options and exit')
    
    args = parser.parse_args()
    
    # Initialize parser
    matrix_parser = MatrixParser(args.results_file)
    
    # Show info if requested
    if args.info:
        print(f"File format: {matrix_parser.format_type}")
        print("Available options:")
        for option in matrix_parser.get_available_options():
            print(f"  {option}")
        return
    
    try:
        # Extract rows and columns based on format
        if matrix_parser.format_type == 'json':
            if args.biclique is not None:
                rows, cols = matrix_parser.extract_rows_cols_from_biclique(args.biclique)
                print(f"Using biclique result {args.biclique}")
            elif args.step is not None:
                rows, cols = matrix_parser.extract_rows_cols_from_step(args.step)
                print(f"Using clustering step {args.step}")
            else:
                print("For JSON format, specify --step or --biclique")
                return
        else:  # text format
            rows, cols = matrix_parser.extract_rows_cols_from_lists(args.list1, args.list2)
            print(f"Using list {args.list1} as rows and list {args.list2} as columns")
        
        # Convert to 1-based indexing by adding 1 to each element
        rows = [r + 1 for r in rows]
        cols = [c + 1 for c in cols]
        rows.append(0)  # Add 0 to rows for 1-based indexing
        cols.append(0)  # Add 0 to cols for 1-based indexing
        cols.sort()  # Ensure columns are sorted
        rows.sort()  # Ensure rows are sorted
        # Filter and save matrix
        matrix_parser.filter_matrix(args.input_csv, rows, cols, args.output_csv)
        
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()