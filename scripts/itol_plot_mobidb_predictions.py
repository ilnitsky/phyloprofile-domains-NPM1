from collections import defaultdict
import re

# Define domain type to shape/color mapping
DOMAIN_STYLES = {
    'Low complexity': ('RE', '#FFB6C1'),  # Light pink
    'Polyampholyte': ('HH', '#87CEEB'),   # Sky blue
    'Polar': ('EL', '#98FB98'),           # Pale green
    'Negative Polyelectrolyte': ('DI', '#FFA07A'),  # Light salmon
    'Positive Polyelectrolyte': ('TR', '#DDA0DD'),  # Plum
    'Proline-rich': ('OC', '#F0E68C'),    # Khaki
    'Glycine-rich': ('HV', '#E6E6FA'),    # Lavender
    'Beta-strand': ('RE', '#4682B4'),     # Steel blue
    'Helix': ('RE', '#FF8C00'),           # Dark orange
}

# Default style for unknown domain types
DEFAULT_STYLE = ('RE', '#CCCCCC')  # Grey rectangle

def parse_mobidb(file_path):
    proteins = defaultdict(list)
    current_protein = None
    
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                protein_id = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                domain_type = parts[3]
                
                # Only store entries with actual domain types (not just '-')
                if domain_type != '-':
                    proteins[protein_id].append((start, end, domain_type))
    
    return proteins

def format_itol_domains(proteins):
    output = []
    
    # Add header
    output.append("DATASET_DOMAINS")
    output.append("SEPARATOR COMMA")
    output.append("DATASET_LABEL,Protein domains")
    output.append("COLOR,#000000")
    output.append("LABEL_AUTO_COLOR,1")
    output.append("BACKBONE_COLOR,#aaaaaa")
    output.append("BACKBONE_HEIGHT,10")
    output.append("BORDER_WIDTH,0")
    output.append("GRADIENT_FILL,0")
    output.append("DATA")
    
    # Process each protein
    for protein_id, domains in proteins.items():
        # Find the maximum end position to determine protein length
        protein_length = max(end for _, end, _ in domains)
        
        # Format domains
        domain_str = []
        for start, end, domain_type in domains:
            shape, color = DOMAIN_STYLES.get(domain_type, DEFAULT_STYLE)
            domain_str.append(f"{shape}|{start}|{end}|{color}|{domain_type}")
        
        # Combine into final string
        protein_line = f"{protein_id},{protein_length},{','.join(domain_str)}"
        output.append(protein_line)
    
    return "\n".join(output)

# Main execution
proteins = parse_mobidb("mobidb_sorted_uniprot.fasta")
itol_output = format_itol_domains(proteins)

# Write output to file
with open("protein_domains.txt", "w") as f:
    f.write(itol_output)