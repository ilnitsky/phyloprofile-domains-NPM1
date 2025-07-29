#!/usr/bin/env python3

import csv
import os

def extract_domains(input_file):
    # Create output file names
    base_name = os.path.splitext(input_file)[0]
    output_file = f"{base_name}_extracted_domains.tsv"
    n_terminal_fasta = f"{base_name}_N_terminal.fasta"
    idr_fasta = f"{base_name}_IDR.fasta"
    c_terminal_fasta = f"{base_name}_C_terminal.fasta"
    
    results = []
    
    with open(input_file, 'r') as f:
        # Read as tab-separated file
        reader = csv.reader(f, delimiter='\t')
        # Get header
        header = next(reader)
        
        # Process each row
        for row in reader:
            if len(row) < 10:  # Skip rows without enough columns
                continue
                
            # Skip rows with "-" for protein or "NOT_FOUND" in N-terminal
            # But we'll keep rows where only C-terminal is "NOT_FOUND"
            if row[4] == "-" or row[7] == "NOT_FOUND":
                continue
                
            species = row[2]
            common_name = row[3]
            protein = row[4]
            uniprot = row[5]
            
            # Get domain coordinates
            try:
                n_terminal_coords = row[7]
                c_terminal_coords = row[8]
                full_seq = row[9]
                
                # Extract N-terminal domain
                if "-" in n_terminal_coords:
                    n_start, n_end = map(int, n_terminal_coords.split("-"))
                    # Adjust for 0-based indexing in Python
                    n_terminal_seq = full_seq[n_start-1:n_end]
                else:
                    n_terminal_seq = ""
                    continue  # Skip if no valid N-terminal
                
                # Extract IDR (region between domains or after N-terminal if C-terminal not found)
                idr_seq = ""
                idr_coords = ""
                
                # Extract C-terminal domain if it's available
                if c_terminal_coords != "NOT_FOUND" and "-" in c_terminal_coords:
                    c_start, c_end = map(int, c_terminal_coords.split("-"))
                    # Adjust for 0-based indexing in Python
                    c_terminal_seq = full_seq[c_start-1:c_end]
                    
                    # IDR is region between N and C terminals
                    if n_end < c_start - 1:  # Check if there is a gap
                        idr_seq = full_seq[n_end:c_start-1]
                        idr_coords = f"{n_end+1}-{c_start-1}"
                else:
                    c_terminal_seq = ""
                    
                    # IDR is everything after N-terminal
                    if n_end < len(full_seq):
                        idr_seq = full_seq[n_end:]
                        idr_coords = f"{n_end+1}-{len(full_seq)}"
                
                # Add to results if at least N-terminal is available
                if n_terminal_seq:
                    results.append({
                        'species': species,
                        'common_name': common_name,
                        'protein': protein,
                        'uniprot': uniprot,
                        'n_terminal_coords': n_terminal_coords,
                        'n_terminal_seq': n_terminal_seq,
                        'idr_coords': idr_coords,
                        'idr_seq': idr_seq,
                        'c_terminal_coords': c_terminal_coords,
                        'c_terminal_seq': c_terminal_seq
                    })
            except (ValueError, IndexError) as e:
                print(f"Error processing row for {species} {protein}: {e}")
    
    # Write results to TSV file
    with open(output_file, 'w') as out_f:
        writer = csv.writer(out_f, delimiter='\t')
        # Write header
        writer.writerow(['Species', 'Common Name', 'Protein', 'Uniprot ID', 
                         'N-terminal Coordinates', 'N-terminal Sequence',
                         'IDR Coordinates', 'IDR Sequence',
                         'C-terminal Coordinates', 'C-terminal Sequence'])
        
        # Write data
        for item in results:
            writer.writerow([
                item['species'],
                item['common_name'],
                item['protein'],
                item['uniprot'],
                item['n_terminal_coords'],
                item['n_terminal_seq'],
                item['idr_coords'],
                item['idr_seq'],
                item['c_terminal_coords'],
                item['c_terminal_seq']
            ])
    
    # Write FASTA files
    # N-terminal domains
    with open(n_terminal_fasta, 'w') as n_fasta:
        for item in results:
            if item['n_terminal_seq']:
                header = f">{item['species']}|{item['protein']}|{item['uniprot']}|N_terminal|{item['n_terminal_coords']}"
                n_fasta.write(f"{header}\n{item['n_terminal_seq']}\n")
    
    # IDR regions
    with open(idr_fasta, 'w') as idr_f:
        for item in results:
            if item['idr_seq']:
                header = f">{item['species']}|{item['protein']}|{item['uniprot']}|IDR|{item['idr_coords']}"
                idr_f.write(f"{header}\n{item['idr_seq']}\n")
    
    # C-terminal domains
    with open(c_terminal_fasta, 'w') as c_fasta:
        for item in results:
            if item['c_terminal_seq']:
                header = f">{item['species']}|{item['protein']}|{item['uniprot']}|C_terminal|{item['c_terminal_coords']}"
                c_fasta.write(f"{header}\n{item['c_terminal_seq']}\n")
    
    print(f"Extracted domains written to {output_file}")
    print(f"N-terminal domains written to {n_terminal_fasta}")
    print(f"IDR regions written to {idr_fasta}")
    print(f"C-terminal domains written to {c_terminal_fasta}")
    
    return results

if __name__ == "__main__":
    input_file = "45_Species/Zoopark_Domains_NPM.tsv"
    extract_domains(input_file) 