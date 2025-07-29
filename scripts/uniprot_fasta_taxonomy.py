import requests
import pandas as pd
from typing import List, Dict
import time
import os
from collections import defaultdict
import re

def read_uniprot_ids(file_path: str) -> List[str]:
    """
    Read UniProt IDs from a file
    """
    try:
        with open(file_path, 'r') as f:
            # Read lines and remove any whitespace
            ids = [line.strip() for line in f if line.strip()]
        return ids
    except FileNotFoundError:
        print(f"Error: File {file_path} not found")
        return []
    except Exception as e:
        print(f"Error reading file: {e}")
        return []

def get_uniprot_data(uniprot_ids: List[str]) -> Dict[str, Dict]:
    """
    Fetch FASTA and taxonomy information for a list of UniProt IDs
    """
    base_url = "https://rest.uniprot.org/uniprotkb/"
    results = {}
    
    for uniprot_id in uniprot_ids:
        # Get FASTA
        fasta_url = f"{base_url}{uniprot_id}.fasta"
        fasta_response = requests.get(fasta_url)
        
        # Get taxonomy information
        entry_url = f"{base_url}{uniprot_id}"
        entry_response = requests.get(entry_url)
        
        if fasta_response.status_code == 200 and entry_response.status_code == 200:
            fasta_data = fasta_response.text
            entry_data = entry_response.json()
            
            # Extract taxonomy information
            try:
                # Get the full taxonomy lineage
                taxonomy = entry_data.get('organism', {}).get('taxonId', 0)
                organism = entry_data.get('organism', {}).get('scientificName', 'Unknown')
                lineage = []
                
                # Get the full lineage if available
                if 'lineage' in entry_data.get('organism', {}):
                    lineage = [taxon.get('scientificName', '') for taxon in entry_data['organism']['lineage']]
                
                results[uniprot_id] = {
                    'fasta': fasta_data,
                    'taxonomy_id': taxonomy,
                    'organism': organism,
                    'lineage': lineage
                }
                print(f"Successfully processed {uniprot_id}")
            except Exception as e:
                print(f"Error processing taxonomy for {uniprot_id}: {e}")
                results[uniprot_id] = {
                    'fasta': fasta_data,
                    'taxonomy_id': 0,
                    'organism': 'Unknown',
                    'lineage': []
                }
        else:
            print(f"Failed to process {uniprot_id}")
        
        # Be nice to the API
        time.sleep(0.1)
    
    return results

def sort_by_taxonomy(uniprot_data: Dict[str, Dict]) -> List[Dict]:
    """
    Sort the UniProt data by taxonomy ID
    """
    sorted_data = sorted(
        uniprot_data.items(),
        key=lambda x: x[1]['taxonomy_id']
    )
    return [{'id': k, **v} for k, v in sorted_data]

def save_results(sorted_data: List[Dict], output_file: str = 'sorted_uniprot.fasta', species_only: bool = False):
    """
    Save the sorted FASTA sequences to a file
    Args:
        sorted_data: List of dictionaries containing UniProt data
        output_file: Path to output file
        species_only: If True, only write species name instead of full taxonomy information
    """
    with open(output_file, 'w') as f:
        for entry in sorted_data:
            if species_only:
                # Write only the species name (organism)
                f.write(f"# Species: {entry['organism']}\n")
            else:
                # Write full taxonomy information
                f.write(f"# Taxonomy ID: {entry['taxonomy_id']}\n")
                f.write(f"# Organism: {entry['organism']}\n")
                if entry['lineage']:
                    f.write(f"# Lineage: {'; '.join(entry['lineage'])}\n")
            f.write(entry['fasta'])
            f.write('\n')

def parse_taxonomy_from_fasta(fasta_file):
    taxonomy_dict = {}
    current_id = None
    
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Parse UniProt ID and taxonomy
                match = re.search(r'\|([A-Z0-9]+)\|.*?([A-Za-z]+/.*?)(?:\s|$)', line)
                if match:
                    uniprot_id = match.group(1)
                    taxonomy = match.group(2).strip('/ ').split('/')
                    taxonomy_dict[uniprot_id] = taxonomy
    
    return taxonomy_dict

def create_tree_structure(taxonomy_dict):
    # Create nested dictionary structure
    tree = defaultdict(lambda: defaultdict(dict))
    
    for uniprot_id, taxonomy in taxonomy_dict.items():
        current = tree
        for rank in taxonomy:
            current = current[rank]
        # Add the leaf (protein ID)
        current[uniprot_id] = {}
    
    return tree

def to_newick(tree, level=0):
    """Convert tree structure to Newick format"""
    if not tree:
        return ""
    
    elements = []
    for key, subtree in tree.items():
        if subtree:
            # Internal node
            subnet = to_newick(subtree, level + 1)
            elements.append(f"{subnet}{key}")
        else:
            # Leaf node (protein ID)
            elements.append(key)
    
    if level == 0:
        return f"({','.join(elements)});"
    else:
        return f"({','.join(elements)})"

def main():
    # Read UniProt IDs from file
    input_file = "NPM/ALphaFoldDB_fetched/member-accessions-A0A0A1TVP1.txt"
    uniprot_ids = read_uniprot_ids(input_file)
    
    if not uniprot_ids:
        print("No UniProt IDs found in the input file")
        return
    
    print(f"Found {len(uniprot_ids)} UniProt IDs to process")
    print("Fetching UniProt data...")
    uniprot_data = get_uniprot_data(uniprot_ids)
    
    print("Sorting by taxonomy...")
    sorted_data = sort_by_taxonomy(uniprot_data)
    
    print("Saving results...")
    save_results(sorted_data, species_only=False)
    
    print("Done! Results saved to sorted_uniprot.fasta")

    # Parse FASTA file
    taxonomy_dict = parse_taxonomy_from_fasta("NPM/mobidb/sorted_uniprot.fasta")
    
    # Create tree structure
    tree = create_tree_structure(taxonomy_dict)
    
    # Convert to Newick format
    newick = to_newick(tree)
    
    # Write to file
    with open("taxonomy_tree.nwk", "w") as f:
        f.write(newick)
    
    # Print some statistics
    print(f"Processed {len(taxonomy_dict)} sequences")
    print("Tree has been written to taxonomy_tree.nwk")

if __name__ == "__main__":
    main() 