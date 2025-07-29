import os
import requests
import subprocess
import pandas as pd
from typing import List, Dict, Tuple, Any
import re
import time
import xml.etree.ElementTree as ET
import json
from collections import defaultdict
import argparse
import sys
from pathlib import Path

# Define domain type to shape/color mapping based on Source_DB (InterPro)
INTERPRO_DOMAIN_STYLES = {
    'pfam': ('RE', '#FF6347'),        # Tomato (highest priority)
    'cathgene3d': ('RE', '#4682B4'),  # Steel Blue
    'panther': ('RE', '#32CD32'),     # Lime Green
    'ssf': ('RE', '#FFD700'),         # Gold
    'profile': ('RE', '#DDA0DD'),     # Plum
    'pirsf': ('RE', '#FFA07A'),       # Light Salmon (lowest priority)
}

# Priority order for overlapping domains (InterPro)
PRIORITY = {
    'pfam': 5,
    'cathgene3d': 4,
    'panther': 3,
    'ssf': 2,
    'profile': 1,
    'pirsf': 0
}

# Default style for unknown source databases (InterPro)
INTERPRO_DEFAULT_STYLE = ('RE', '#CCCCCC')  # Grey rectangle

# Define domain type to shape/color mapping for MobiDB
MOBIDB_DOMAIN_STYLES = {
    'Low complexity': ('RE', '#FFB6C1'),          # Light pink
    'Polyampholyte': ('HH', '#87CEEB'),           # Sky blue
    'Polar': ('EL', '#98FB98'),                   # Pale green
    'Negative Polyelectrolyte': ('DI', '#FFA07A'), # Light salmon
    'Positive Polyelectrolyte': ('TR', '#DDA0DD'), # Plum
    'Proline-rich': ('OC', '#F0E68C'),            # Khaki
    'Glycine-rich': ('HV', '#E6E6FA'),            # Lavender
    'Beta-strand': ('RE', '#1E90FF'),             # DodgerBlue 
    'Helix': ('RE', '#87CEFA'),                   # LightSkyBlue 
}

# Default style for unknown domain types (MobiDB)
MOBIDB_DEFAULT_STYLE = ('RE', '#CCCCCC')  # Grey rectangle

def read_uniprot_ids(file_path: str) -> List[str]:
    """Read UniProt IDs from a file"""
    try:
        with open(file_path, 'r') as f:
            ids = [line.strip() for line in f if line.strip()]
        return ids
    except FileNotFoundError:
        print(f"Error: File {file_path} not found")
        return []
    except Exception as e:
        print(f"Error reading file: {e}")
        return []

def parse_fasta_header(header: str) -> Dict:
    """Parse UniProt FASTA header to extract taxonomy information"""
    organism_match = re.search(r'OS=(.*?)(?=\sOX=|\sPE=|\sGN=|\s\w+=|$)', header)
    taxonomy_match = re.search(r'OX=(\d+)', header)
    
    organism = organism_match.group(1) if organism_match else "Unknown"
    taxonomy_id = int(taxonomy_match.group(1)) if taxonomy_match else 0
    
    return {
        'organism': organism,
        'taxonomy_id': taxonomy_id
    }

def get_uniprot_data(uniprot_ids: List[str]) -> Dict[str, Dict]:
    """Fetch FASTA sequences and extract taxonomy information from headers"""
    base_url = "https://rest.uniprot.org/uniprotkb/"
    results = {}
    
    for uniprot_id in uniprot_ids:
        try:
            fasta_url = f"{base_url}{uniprot_id}.fasta"
            fasta_response = requests.get(fasta_url)
            
            if fasta_response.status_code == 200:
                fasta_text = fasta_response.text
                header = fasta_text.split('\n')[0]
                tax_info = parse_fasta_header(header)
                
                results[uniprot_id] = {
                    'fasta': fasta_text,
                    'taxonomy_id': tax_info['taxonomy_id'],
                    'organism': tax_info['organism']
                }
                print(f"Successfully processed {uniprot_id}")
            else:
                print(f"Failed to get FASTA for {uniprot_id}")
            
            time.sleep(0.1)
            
        except Exception as e:
            print(f"Error processing {uniprot_id}: {e}")
            continue
    
    return results

def sort_by_taxonomy(uniprot_data: Dict[str, Dict]) -> List[Dict]:
    """Sort the UniProt data by taxonomy ID and organism name"""
    sorted_data = sorted(
        uniprot_data.items(),
        key=lambda x: (x[1]['taxonomy_id'], x[1]['organism'])
    )
    return [{'id': k, **v} for k, v in sorted_data]

def get_taxonomy_lineage(taxid: int) -> Dict:
    """Get full taxonomic lineage from NCBI Taxonomy using E-utilities"""
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    efetch_url = f"{base_url}efetch.fcgi?db=taxonomy&id={taxid}&format=xml"
    
    try:
        response = requests.get(efetch_url)
        if response.status_code == 200:
            root = ET.fromstring(response.text)
            lineage_elem = root.find('.//Lineage')
            if lineage_elem is not None:
                lineage = lineage_elem.text.split('; ')
                if lineage[0] == "cellular organisms":
                    lineage = lineage[1:]
                formatted_lineage = '/'.join(lineage)
                return {
                    'taxid': taxid,
                    'lineage': formatted_lineage
                }
    except Exception as e:
        print(f"Error fetching taxonomy for {taxid}: {e}")
    
    return {'taxid': taxid, 'lineage': ''}

def get_secondary_structure_regions(uniprot_id: str, output_dir: str = "pdbs") -> List[Tuple[int, int, str]]:
    """runs alphafold_disorder.py, and returns beta strand and helix regions.
    Returns empty list if PDB is not available."""
    os.makedirs(output_dir, exist_ok=True)
    
    pdb_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    pdb_path = os.path.join(output_dir, f"AF-{uniprot_id}-F1-model_v4.pdb")
    
    print(f"Downloading PDB for {uniprot_id}...")
    response = requests.get(pdb_url)
    if response.status_code != 200:
        print(f"No PDB available for {uniprot_id} (HTTP {response.status_code}), skipping secondary structure analysis.")
        return []
    
    with open(pdb_path, "w") as f:
        f.write(response.text)
    
    disorder_script = "/home/ilnitsky/tools/AlphaFold-disorder/alphafold_disorder.py"
    output_tsv = os.path.join(output_dir, f"{uniprot_id}_af_disorder.tsv")
    output_tsv_data = os.path.join(output_dir, f"{uniprot_id}_af_disorder_data.tsv")
    
    print(f"Running alphafold_disorder.py for {uniprot_id}...")
    result = subprocess.run(
        [disorder_script, "-i", pdb_path, "-o", output_tsv],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"alphafold_disorder.py failed for {uniprot_id}: {result.stderr}, skipping.")
        return []
    
    print(f"Parsing {output_tsv}...")
    if not os.path.exists(output_tsv_data):
        print(f"Expected output file {output_tsv_data} not found for {uniprot_id}, skipping.")
        return []
    
    df = pd.read_csv(output_tsv_data, sep="\t")
    
    regions = []
    current_start = None
    current_type = None
    
    for index, row in df.iterrows():
        pos = row["pos"]
        ss = row["ss"]
        
        if ss in ["E", "H"]:
            if current_type is None:
                current_start = pos
                current_type = "Beta-strand" if ss == "E" else "Helix"
            elif ss != ("E" if current_type == "Beta-strand" else "H"):
                regions.append((current_start, pos - 1, current_type))
                current_start = pos
                current_type = "Beta-strand" if ss == "E" else "Helix"
        elif current_type is not None:
            regions.append((current_start, pos - 1, current_type))
            current_start = None
            current_type = None
    
    if current_type is not None:
        regions.append((current_start, df["pos"].iloc[-1], current_type))
    
    return regions

def get_interpro_domains(uniprot_id: str) -> List[Dict[str, Any]]:
    """Fetch domain information from InterPro for a given UniProt ID, including protein length."""
    databases = ["pfam", "panther", "cathgene3d", "ssf", "pirsf", "profile"]
    base_url = "https://www.ebi.ac.uk/interpro/api/entry"
    domains = []
    
    for db in databases:
        url = f"{base_url}/{db}/protein/uniprot/{uniprot_id}?extra_fields=description,go_terms"
        try:
            response = requests.get(url)
            if response.status_code == 200:
                data = response.json()
                for result in data.get('results', []):
                    metadata = result.get('metadata', {})
                    proteins = result.get('proteins', [])
                    
                    if proteins and len(proteins) > 0:
                        protein = proteins[0]
                        protein_length = protein.get('protein_length', 0)
                        locations = protein.get('entry_protein_locations', [])
                        
                        if locations and len(locations) > 0:
                            location = locations[0]
                            fragments = location.get('fragments', [])
                            
                            if fragments and len(fragments) > 0:
                                fragment = fragments[0]
                                domain_info = {
                                    'source_database': metadata.get('source_database', ''),
                                    'accession': metadata.get('accession', ''),
                                    'name': metadata.get('name', ''),
                                    'start': fragment.get('start', 0),
                                    'end': fragment.get('end', 0),
                                    'score': location.get('score', None),
                                    'protein_length': protein_length
                                }
                                domains.append(domain_info)
            
            time.sleep(0.1)
            
        except Exception as e:
            print(f"Error fetching {db} domains for {uniprot_id}: {e}")
    
    return domains

def resolve_overlaps(domains: List[Tuple[int, int, str, str]]) -> List[Tuple[int, int, str, str]]:
    """Resolve overlapping domains based on priority."""
    if not domains:
        return []
    
    domains.sort(key=lambda x: (x[0], -PRIORITY.get(x[2], 0)))
    
    resolved = []
    current_end = -1
    
    for start, end, source_db, name in domains:
        if start > current_end:
            resolved.append((start, end, source_db, name))
            current_end = end
        elif PRIORITY.get(source_db, 0) > PRIORITY.get(resolved[-1][2], 0):
            resolved[-1] = (start, end, source_db, name)
            current_end = end
    
    return resolved

def format_interpro_itol_domains(proteins: Dict[str, Tuple[List[Tuple[int, int, str, str]], int]]) -> str:
    """Format InterPro domains into iTOL DATASET_DOMAINS format using actual protein length."""
    output = []
    
    output.append("DATASET_DOMAINS")
    output.append("SEPARATOR COMMA")
    output.append("DATASET_LABEL,InterPro Domains")
    output.append("COLOR,#000000")
    output.append("LABEL_AUTO_COLOR,1")
    output.append("BACKBONE_COLOR,#aaaaaa")
    output.append("BACKBONE_HEIGHT,10")
    output.append("BORDER_WIDTH,0")
    output.append("GRADIENT_FILL,0")
    output.append("DATA")
    
    for protein_id, (domains, protein_length) in proteins.items():
        resolved_domains = resolve_overlaps(domains)
        
        domain_str = []
        for start, end, source_db, name in resolved_domains:
            shape, color = INTERPRO_DOMAIN_STYLES.get(source_db, INTERPRO_DEFAULT_STYLE)
            domain_str.append(f"{shape}|{start}|{end}|{color}|{name}")
        
        if domain_str:
            protein_line = f"{protein_id},{protein_length},{','.join(domain_str)}"
            output.append(protein_line)
    
    return "\n".join(output)

def parse_mobidb(file_path: str) -> Dict[str, List[Tuple[int, int, str]]]:
    """Parse the full_annotation.txt file for MobiDB and secondary structure regions."""
    proteins = defaultdict(list)
    
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

def format_mobidb_itol_domains(proteins: Dict[str, List[Tuple[int, int, str]]]) -> str:
    """Format MobiDB and secondary structure domains into iTOL DATASET_DOMAINS format with short labels."""
    # Mapping of full domain types to short labels
    SHORT_LABELS = {
        'Low complexity': 'Low comp',
        'Polyampholyte': 'Polyamph',
        'Polar': 'Polar',
        'Negative Polyelectrolyte': 'NPE',
        'Positive Polyelectrolyte': 'PPE',
        'Proline-rich': 'Pro-rich',
        'Glycine-rich': 'Gly-rich',  
        'Beta-strand': 'B',
        'Helix': 'H'
    }
    
    output = []
    
    output.append("DATASET_DOMAINS")
    output.append("SEPARATOR COMMA")
    output.append("DATASET_LABEL,MobiDB and Secondary Structure")
    output.append("COLOR,#000000")
    output.append("LABEL_AUTO_COLOR,1")
    output.append("BACKBONE_COLOR,#aaaaaa")
    output.append("BACKBONE_HEIGHT,10")
    output.append("BORDER_WIDTH,0")
    output.append("GRADIENT_FILL,0")
    output.append("DATA")
    
    for protein_id, domains in proteins.items():
        # Find the maximum end position to determine protein length
        protein_length = max(end for _, end, _ in domains)
        
        domain_str = []
        for start, end, domain_type in domains:
            shape, color = MOBIDB_DOMAIN_STYLES.get(domain_type, MOBIDB_DEFAULT_STYLE)
            # Use the short label if available, otherwise use the full domain_type
            short_label = SHORT_LABELS.get(domain_type, domain_type)
            domain_str.append(f"{shape}|{start}|{end}|{color}|{short_label}")
        
        if domain_str:
            protein_line = f"{protein_id},{protein_length},{','.join(domain_str)}"
            output.append(protein_line)
    
    return "\n".join(output)

def save_results(sorted_data: List[Dict], output_dir: str = 'output', output_prefix: str = 'sorted_uniprot', species_only: bool = False):
    """Save the sorted FASTA sequences with lineage information, secondary structure regions,
    domain information, and iTOL domain datasets, skipping Unknown entries."""
    os.makedirs(output_dir, exist_ok=True)
    
    fasta_dir = os.path.join(output_dir, "fasta")
    structure_dir = os.path.join(output_dir, "structure")
    domains_dir = os.path.join(output_dir, "domains")
    tree_dir = os.path.join(output_dir, "tree")
    
    os.makedirs(fasta_dir, exist_ok=True)
    os.makedirs(structure_dir, exist_ok=True)
    os.makedirs(domains_dir, exist_ok=True)
    os.makedirs(tree_dir, exist_ok=True)
    
    output_file = os.path.join(fasta_dir, f"{output_prefix}.fasta")
    ss_output_file = os.path.join(structure_dir, f"{output_prefix}_secondary_structure.txt")
    domain_output_file = os.path.join(domains_dir, f"{output_prefix}_domains.txt")
    interpro_itol_file = os.path.join(domains_dir, f"{output_prefix}_interpro_domains_itol.txt")
    tree_file = os.path.join(tree_dir, f"{output_prefix}.nwk")
    
    taxonomy_dict = {}
    proteins_domains = defaultdict(lambda: ([], 0))  # (domains_list, protein_length)
    
    with open(output_file, 'w') as f_fasta, \
         open(ss_output_file, 'w') as f_ss, \
         open(domain_output_file, 'w') as f_domain:
        
        f_domain.write("Organism\tUniProt_ID\tSource_DB\tAccession\tName\tStart\tEnd\tScore\n")
        
        for i, entry in enumerate(sorted_data):
            if entry['organism'].startswith('Unknown'):
                print(f"Skipping {entry['id']} as it is marked as Unknown")
                continue
            
            tax_info = get_taxonomy_lineage(entry['taxonomy_id'])
            taxonomy_dict[entry['id']] = {
                'taxid': tax_info['taxid'],
                'lineage': tax_info['lineage'],
                'organism': entry['organism']
            }
            
            fasta_lines = entry['fasta'].split('\n')
            header = fasta_lines[0]
            sequence = '\n'.join(line for line in fasta_lines[1:] if line.strip())
            
            organism_name = entry['organism'].replace(' ', '_')
            uniprot_id = entry['id']
            combined_id = f"{organism_name}-{uniprot_id}"
            
            try:
                if header.startswith('>'):
                    header = header[1:]
                
                id_part = header.split(' ')[0]
                
                if species_only:
                    new_header = f">{combined_id}"
                else:
                    desc_part = ' '.join(header.split(' ')[1:])
                    new_header = f">{id_part} {desc_part.split(' OS=')[0]} {tax_info['lineage']}/ {entry['organism']} OX={entry['taxonomy_id']}"
                    if 'GN=' in desc_part:
                        new_header += f" {' '.join(part for part in desc_part.split(' ') if part.startswith(('GN=', 'PE=', 'SV=')))}"
                
            except Exception as e:
                print(f"Error formatting header for entry {entry.get('id', 'unknown')}: {e}")
                new_header = f">{combined_id}" if species_only else f"{header} {tax_info['lineage']}/"
            
            f_fasta.write(f"{new_header}\n")
            f_fasta.write(f"{sequence}\n")

            try:
                regions = get_secondary_structure_regions(uniprot_id, output_dir=os.path.join(structure_dir, "pdbs"))
                if regions:  # Only write if there are regions
                    for start, end, struct_type in regions:
                        f_ss.write(f"{combined_id}\t{start}\t{end}\t{struct_type}\n")
                else:
                    print(f"No secondary structure regions found for {uniprot_id}, continuing with other data.")
            except Exception as e:
                print(f"Error processing secondary structure for {uniprot_id}: {e}, continuing with other data.")
                
            try:
                print(f"Fetching domain information for {uniprot_id}...")
                domains = get_interpro_domains(uniprot_id)
                for domain in domains:
                    f_domain.write(f"{organism_name}\t{uniprot_id}\t{domain['source_database']}\t{domain['accession']}\t{domain['name']}\t{domain['start']}\t{domain['end']}\t{domain['score'] if domain['score'] is not None else 'null'}\n")
                    domain_list, _ = proteins_domains[combined_id]
                    domain_list.append((domain['start'], domain['end'], domain['source_database'], domain['name']))
                    proteins_domains[combined_id] = (domain_list, domain['protein_length'])
            except Exception as e:
                print(f"Error processing domain information for {uniprot_id}: {e}")

    with open(interpro_itol_file, 'w') as f_itol:
        itol_output = format_interpro_itol_domains(proteins_domains)
        f_itol.write(itol_output)
        print(f"InterPro iTOL domain dataset written to {interpro_itol_file}")

    with open(tree_file, 'w') as f:
        tree = {}
        organisms = defaultdict(list)
        for entry_id, tax_data in taxonomy_dict.items():
            organism = tax_data['organism'].replace(' ', '_')
            organisms[organism].append(entry_id)
        
        for organism, uniprot_ids in organisms.items():
            lineage = taxonomy_dict[uniprot_ids[0]]['lineage'].split('/')
            current_node = tree
            
            for level in lineage:
                if level:
                    level = level.strip().replace(' ', '_')
                    if level not in current_node:
                        current_node[level] = {}
                    current_node = current_node[level]
            
            if len(uniprot_ids) > 1:
                paralogue_str = ','.join(f"{organism}-{uid}" for uid in uniprot_ids)
                current_node[organism] = {f"({paralogue_str}){organism}": {}}
            else:
                current_node[organism] = {f"{organism}-{uniprot_ids[0]}": {}}

        def dict_to_newick(node):
            if not node:
                return ""
            children = []
            for key, value in node.items():
                child_str = dict_to_newick(value)
                if child_str:
                    children.append(f"({child_str}){key}")
                else:
                    children.append(key)
            return ','.join(children)

        newick_str = dict_to_newick(tree)
        f.write(f"({newick_str});")
    
    return {
        'fasta': output_file,
        'ss': ss_output_file,
        'domains': domain_output_file,
        'interpro_itol': interpro_itol_file,
        'tree': tree_file
    }

def run_mobidb_and_create_annotation(
    mobidb_path="/home/ilnitsky/tools/MobiDB-lite/src/mobidb_lite/__main__.py",
    fasta_file="sorted_uniprot.fasta",
    output_dir="output",
    ss_file="sorted_uniprot_secondary_structure.txt",
    output_prefix="sorted_uniprot"
):
    """Run MobiDB-lite, create full annotation, and generate MobiDB iTOL file."""
    structure_dir = os.path.join(output_dir, "structure")
    domains_dir = os.path.join(output_dir, "domains")
    os.makedirs(structure_dir, exist_ok=True)
    os.makedirs(domains_dir, exist_ok=True)
    
    mobidb_output = os.path.join(structure_dir, f"{output_prefix}_mobidb.txt")
    full_annotation_file = os.path.join(structure_dir, f"{output_prefix}_full_annotation.txt")
    mobidb_itol_file = os.path.join(domains_dir, f"{output_prefix}_mobidb_domains_itol.txt")
    
    print(f"Running MobiDB-lite on {fasta_file}...")
    try:
        subprocess.run(
            ["python3", mobidb_path, fasta_file, mobidb_output],
            check=True,
            capture_output=True,
            text=True
        )
        print(f"MobiDB-lite analysis completed and saved to {mobidb_output}")
    except subprocess.CalledProcessError as e:
        print(f"Error running MobiDB-lite: {e}")
        print(f"Error output: {e.stderr}")
        return None
    
    print(f"Creating combined annotation file {full_annotation_file}...")
    try:
        result = subprocess.run(
            f"cat {mobidb_output} {ss_file} | sort",
            check=True,
            capture_output=True,
            text=True,
            shell=True
        )
        
        with open(full_annotation_file, "w") as f:
            f.write(result.stdout)
        
        print(f"Successfully created {full_annotation_file}")
        
        # Parse MobiDB and secondary structure data and create iTOL file
        print(f"Parsing {full_annotation_file} for iTOL...")
        mobidb_proteins = parse_mobidb(full_annotation_file)
        mobidb_itol_output = format_mobidb_itol_domains(mobidb_proteins)
        
        with open(mobidb_itol_file, "w") as f_itol:
            f_itol.write(mobidb_itol_output)
        print(f"MobiDB iTOL domain dataset written to {mobidb_itol_file}")
        
        return full_annotation_file, mobidb_itol_file
    except subprocess.CalledProcessError as e:
        print(f"Error creating combined annotation: {e}")
        print(f"Error output: {e.stderr}")
        return None, None

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Process UniProt IDs to generate structured protein data")
    
    parser.add_argument("-i", "--input", required=True,
                        help="Input file containing UniProt IDs (one per line)")
    
    parser.add_argument("-o", "--output-dir", default="output",
                        help="Output directory for all generated files (default: 'output')")
    
    parser.add_argument("-p", "--prefix", default="sorted_uniprot",
                        help="Prefix for all output files (default: 'sorted_uniprot')")
    
    parser.add_argument("-m", "--mobidb-path", 
                        default="/home/ilnitsky/tools/MobiDB-lite/src/mobidb_lite/__main__.py",
                        help="Path to MobiDB-lite __main__.py script")
    
    parser.add_argument("-s", "--species-only", action="store_true",
                        help="Use only species name in FASTA headers")
    
    return parser.parse_args()

def main():
    """Main function with command-line argument support"""
    args = parse_arguments()
    
    if not os.path.exists(args.input):
        print(f"Error: Input file {args.input} not found")
        sys.exit(1)
    
    uniprot_ids = read_uniprot_ids(args.input)
    
    if not uniprot_ids:
        print("No UniProt IDs found in the input file")
        sys.exit(1)
    
    print(f"Found {len(uniprot_ids)} UniProt IDs to process")
    print("Fetching UniProt data...")
    uniprot_data = get_uniprot_data(uniprot_ids)
    
    if not uniprot_data:
        print("No data was retrieved from UniProt")
        sys.exit(1)
        
    print("Sorting by taxonomy...")
    sorted_data = sort_by_taxonomy(uniprot_data)
    
    print(f"Saving results to {args.output_dir}...")
    output_files = save_results(
        sorted_data, 
        output_dir=args.output_dir,
        output_prefix=args.prefix,
        species_only=args.species_only
    )
    
    print("Running MobiDB-lite and creating full annotation...")
    annotation_file, mobidb_itol_file = run_mobidb_and_create_annotation(
        mobidb_path=args.mobidb_path,
        fasta_file=output_files['fasta'],
        output_dir=args.output_dir,
        ss_file=output_files['ss'],
        output_prefix=args.prefix
    )
    
    print("\nProcessing complete. Files generated:")
    print(f"FASTA: {output_files['fasta']}")
    print(f"Secondary Structure: {output_files['ss']}")
    print(f"InterPro Domains: {output_files['domains']}")
    print(f"InterPro iTOL Domains: {output_files['interpro_itol']}")
    print(f"Phylogenetic Tree: {output_files['tree']}")
    if annotation_file:
        print(f"Full Annotation: {annotation_file}")
    if mobidb_itol_file:
        print(f"MobiDB iTOL Domains: {mobidb_itol_file}")

if __name__ == "__main__":
    main()