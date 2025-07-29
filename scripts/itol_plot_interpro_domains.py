from collections import defaultdict

# Define domain type to shape/color mapping based on Source_DB
DOMAIN_STYLES = {
    'pfam': ('RE', '#FF6347'),        # Tomato (highest priority)
    'cathgene3d': ('RE', '#4682B4'),  # Steel Blue
    'panther': ('RE', '#32CD32'),     # Lime Green
    'ssf': ('RE', '#FFD700'),         # Gold
    'profile': ('RE', '#DDA0DD'),     # Plum
    'pirsf': ('RE', '#FFA07A'),       # Light Salmon (lowest priority)
}

# Priority order for overlapping domains
PRIORITY = {
    'pfam': 5,
    'cathgene3d': 4,
    'panther': 3,
    'ssf': 2,
    'profile': 1,
    'pirsf': 0
}

# Default style for unknown source databases
DEFAULT_STYLE = ('RE', '#CCCCCC')  # Grey rectangle

def parse_uniprot_domains(file_path):
    proteins = defaultdict(list)
    
    with open(file_path, 'r') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 7:
                protein_id = parts[1]  # UniProt_ID
                source_db = parts[2]   # Source_DB
                name = parts[4]        # Name (instead of Accession)
                start = int(parts[5])  # Start
                end = int(parts[6])    # End
                
                # Store domain entry with source_db and name
                proteins[protein_id].append((start, end, source_db, name))
    
    return proteins

def resolve_overlaps(domains):
    """Resolve overlapping domains based on priority."""
    if not domains:
        return []
    
    # Sort by start position, then by priority (descending)
    domains.sort(key=lambda x: (x[0], -PRIORITY.get(x[2], 0)))
    
    resolved = []
    current_end = -1
    
    for start, end, source_db, name in domains:
        # If no overlap with previous, add directly
        if start > current_end:
            resolved.append((start, end, source_db, name))
            current_end = end
        # If overlap, only add if higher priority (already sorted by priority)
        elif PRIORITY.get(source_db, 0) > PRIORITY.get(resolved[-1][2], 0):
            resolved[-1] = (start, end, source_db, name)
            current_end = end
    
    return resolved

def format_itol_domains(proteins):
    output = []
    
    # Add iTOL header
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
        # Resolve overlapping domains
        resolved_domains = resolve_overlaps(domains)
        
        # Find the maximum end position to determine protein length
        protein_length = max(end for _, end, _, _ in resolved_domains) if resolved_domains else 0
        
        # Format domains
        domain_str = []
        for start, end, source_db, name in resolved_domains:
            shape, color = DOMAIN_STYLES.get(source_db, DEFAULT_STYLE)
            domain_str.append(f"{shape}|{start}|{end}|{color}|{name}")
        
        # Combine into final string
        if domain_str:  # Only include proteins with domains
            protein_line = f"{protein_id},{protein_length},{','.join(domain_str)}"
            output.append(protein_line)
    
    return "\n".join(output)

# Main execution
file_path = "sorted_uniprot_domains.txt"
proteins = parse_uniprot_domains(file_path)
itol_output = format_itol_domains(proteins)

# Write output to file
with open("uniprot_domains_itol.txt", "w") as f:
    f.write(itol_output)

print("iTOL domain file written to 'uniprot_domains_itol.txt'")