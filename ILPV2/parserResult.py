import re
from dataclasses import dataclass
from typing import List, Dict, Optional

@dataclass
class Step:
    step_number: int
    columns_used: List[int]
    group_1_reads: List[str]
    group_1_count: int
    group_0_reads: List[str]
    group_0_count: int

@dataclass
class ClusteringResult:
    total_steps: int
    total_reads: int
    steps: List[Step]

class ClusteringParser:
    def __init__(self):
        self.step_pattern = re.compile(r'^Step (\d+):$')
        self.columns_pattern = re.compile(r'Columns used: \[([0-9, ]+)\]')
        self.group_pattern = re.compile(r'Group ([01]) reads \((\d+)\): (.+)')
    
    def parse_file(self, filepath: str) -> ClusteringResult:
        """Parse un fichier de résultats de clustering"""
        with open(filepath, 'r') as file:
            content = file.read()
        
        return self.parse_content(content)
    
    def parse_content(self, content: str) -> ClusteringResult:
        """Parse le contenu d'un fichier de clustering"""
        lines = content.strip().split('\n')
        
        # Extraire les informations générales
        total_steps = 0
        total_reads = 0
        
        for line in lines:
            if line.startswith('Total steps:'):
                total_steps = int(line.split(':')[1].strip())
            elif line.startswith('Total reads:'):
                total_reads = int(line.split(':')[1].strip())
        
        # Parser les étapes
        steps = []
        current_step = None
        
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            
            # Détecter une nouvelle étape
            step_match = self.step_pattern.match(line)
            if step_match:
                step_number = int(step_match.group(1))
                
                # Lire les colonnes utilisées
                i += 1
                columns_line = lines[i].strip()
                columns_match = self.columns_pattern.search(columns_line)
                columns_used = []
                if columns_match:
                    columns_str = columns_match.group(1)
                    columns_used = [int(x.strip()) for x in columns_str.split(',')]
                
                # Lire les groupes
                group_1_reads = []
                group_1_count = 0
                group_0_reads = []
                group_0_count = 0
                
                # Group 1
                i += 1
                group1_line = lines[i].strip()
                group1_match = self.group_pattern.match(group1_line)
                if group1_match:
                    group_1_count = int(group1_match.group(2))
                    reads_str = group1_match.group(3)
                    group_1_reads = self._parse_reads_list(reads_str)
                
                # Group 0
                i += 1
                group0_line = lines[i].strip()
                group0_match = self.group_pattern.match(group0_line)
                if group0_match:
                    group_0_count = int(group0_match.group(2))
                    reads_str = group0_match.group(3)
                    group_0_reads = self._parse_reads_list(reads_str)
                
                step = Step(
                    step_number=step_number,
                    columns_used=columns_used,
                    group_1_reads=group_1_reads,
                    group_1_count=group_1_count,
                    group_0_reads=group_0_reads,
                    group_0_count=group_0_count
                )
                steps.append(step)
            
            i += 1
        
        return ClusteringResult(
            total_steps=total_steps,
            total_reads=total_reads,
            steps=steps
        )
    
    def _parse_reads_list(self, reads_str: str) -> List[str]:
        """Parse une liste de reads depuis une chaîne"""
        # Format: "read_548, read_566, read_623 ... and 26 more"
        reads = []
        
        # Diviser par "... and"
        if "... and" in reads_str:
            main_part = reads_str.split("... and")[0]
        else:
            main_part = reads_str
        
        # Extraire les noms de reads
        read_names = [name.strip() for name in main_part.split(',')]
        reads = [name for name in read_names if name.startswith('read_')]
        
        return reads

# Exemple d'utilisation
if __name__ == "__main__":
    parser = ClusteringParser()
    
    # Parser un fichier
    result = parser.parse_file("resultats/steps_ctg0_75000_80000.txt")
    
    print(f"Total steps: {result.total_steps}")
    print(f"Total reads: {result.total_reads}")
    print(f"Parsed {len(result.steps)} steps")
    
    for step in result.steps:
        print(f"\nStep {step.step_number}:")
        print(f"  Columns: {step.columns_used}")
        print(f"  Group 1: {step.group_1_count} reads")
        print(f"  Group 0: {step.group_0_count} reads")
        print(f"  First few Group 1 reads: {step.group_1_reads[:5]}")