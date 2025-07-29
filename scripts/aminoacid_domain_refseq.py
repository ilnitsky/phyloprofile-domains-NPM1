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
import glob
from Bio import SeqIO

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

def get_refseq_info(accession: str) -> Dict:
    """Get taxonomy information for RefSeq accession using NCBI E-utilities."""
    taxonomy_id = 0
    organism = "Unknown"
    
    try:
        # Get the protein GI number
        esearch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term={accession}"
        response = requests.get(esearch_url)
        if response.status_code != 200:
            print(f"Failed to query NCBI for {accession}: HTTP {response.status_code}")
            return {'organism': organism, 'taxonomy_id': taxonomy_id}
            
        protein_id_match = re.search(r'<Id>(\d+)</Id>', response.text)
        if not protein_id_match:
            print(f"No protein ID found for {accession}")
            return {'organism': organism, 'taxonomy_id': taxonomy_id}
            
        protein_id = protein_id_match.group(1)
        
        # Get the protein record in GenPept XML format
        efetch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={protein_id}&rettype=gp&retmode=xml"
        response = requests.get(efetch_url)
        if response.status_code != 200:
            print(f"Failed to fetch protein record for {accession}: HTTP {response.status_code}")
            return {'organism': organism, 'taxonomy_id': taxonomy_id}
            
        # Parse the XML to extract taxonomy ID and organism name
        root = ET.fromstring(response.text)
        
        # Find the taxonomy ID using a simpler approach
        for qualifier in root.findall(".//GBQualifier"):
            name_elem = qualifier.find("GBQualifier_name")
            value_elem = qualifier.find("GBQualifier_value")
            
            if name_elem is not None and value_elem is not None:
                if name_elem.text == "db_xref" and "taxon:" in value_elem.text:
                    taxid_match = re.search(r'taxon:(\d+)', value_elem.text)
                    if taxid_match:
                        taxonomy_id = int(taxid_match.group(1))
                        break
        
        # Find the organism name
        organism_elem = root.find(".//GBSeq_organism")
        if organism_elem is not None:
            organism = organism_elem.text
        
        time.sleep(0.1)  # Be nice to NCBI servers
    
    except Exception as e:
        print(f"Error retrieving information for {accession}: {e}")
    
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

def create_accession_to_folder_mapping(iprscan_dir: str) -> Dict[str, str]:
    """Create a mapping from RefSeq accessions to their corresponding folder names
    by parsing the InterProScan TSV file names and contents."""
    mapping = {}
    
    # Get all TSV files in the InterProScan directory
    iprscan_files = glob.glob(os.path.join(iprscan_dir, "*.tsv"))
    
    for iprscan_file in iprscan_files:
        # Extract the folder name from the file name (remove .fasta.tsv)
        folder_name = os.path.basename(iprscan_file).replace('.fasta.tsv', '')
        
        # Read the first column of the TSV file to get all accessions
        try:
            with open(iprscan_file, 'r') as f:
                for line in f:
                    if line.strip():
                        accession = line.split('\t')[0]
                        mapping[accession] = folder_name
        except Exception as e:
            print(f"Error reading {iprscan_file}: {e}")
    
    return mapping

def get_secondary_structure_from_fold_res(organism: str, protein_id: str, fold_res_dir: str = "fold_res", accession_mapping: Dict[str, str] = None) -> List[Tuple[int, int, str]]:
    """Parse PDB files from fold_res directory structure and extract secondary structure.
    Directory structure: fold_res/Genus_species/Genus_species/ranked_3.pdb
    Or for multiple proteins: fold_res/Genus_species_1/Genus_species_1/ranked_3.pdb
    Returns beta strand and helix regions as a list of tuples."""
    
    # Use the mapping if available, otherwise fall back to organism name
    if accession_mapping and protein_id in accession_mapping:
        organism_dir = accession_mapping[protein_id]
    else:
        # Format organism name to match directory structure (replace spaces with underscores)
        organism_dir = organism.replace(' ', '_')
        
        # Check if this is a multi-protein organism (Genus_species_1, Genus_species_2, etc.)
        if '_' in protein_id and protein_id.split('_')[-1].isdigit():
            # Extract the number suffix
            suffix = protein_id.split('_')[-1]
            organism_dir = f"{organism_dir}_{suffix}"
    
    # Path to the PDB file
    pdb_path = os.path.join(fold_res_dir, organism_dir, organism_dir, "ranked_3.pdb")
    
    if not os.path.exists(pdb_path):
        print(f"No PDB file found at {pdb_path}, skipping secondary structure analysis.")
        return []
    
    print(f"Found PDB file for {organism_dir} at {pdb_path}")
    
    # Create output directory for temporary files
    output_dir = os.path.join("structure", "pdbs")
    os.makedirs(output_dir, exist_ok=True)
    
    # Run alphafold_disorder.py to analyze secondary structure
    disorder_script = "/home/ilnitsky/tools/AlphaFold-disorder/alphafold_disorder.py"
    output_tsv = os.path.join(output_dir, f"{protein_id}_af_disorder.tsv")
    output_tsv_data = os.path.join(output_dir, f"{protein_id}_af_disorder_data.tsv")
    
    print(f"Running alphafold_disorder.py for {protein_id}...")
    result = subprocess.run(
        [disorder_script, "-i", pdb_path, "-o", output_tsv],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"alphafold_disorder.py failed for {protein_id}: {result.stderr}, skipping.")
        return []
    
    print(f"Parsing {output_tsv}...")
    if not os.path.exists(output_tsv_data):
        print(f"Expected output file {output_tsv_data} not found for {protein_id}, skipping.")
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

def save_results(sorted_data: List[Dict], output_dir: str = 'output', output_prefix: str = 'sorted_uniprot', species_only: bool = False, fold_res_dir: str = 'fold_res', accession_mapping: Dict[str, str] = None):
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
                regions = get_secondary_structure_from_fold_res(entry['organism'], uniprot_id, fold_res_dir=fold_res_dir, accession_mapping=accession_mapping)
                if regions:  # Only write if there are regions
                    for start, end, struct_type in regions:
                        f_ss.write(f"{combined_id}\t{start}\t{end}\t{struct_type}\n")
                else:
                    print(f"No secondary structure regions found for {uniprot_id}, continuing with other data.")
            except Exception as e:
                print(f"Error processing secondary structure for {uniprot_id}: {e}, continuing with other data.")
                
            # Use pre-parsed domain information if available
            try:
                if 'domains' in entry and entry['domains']:
                    print(f"Using pre-parsed domain information for {uniprot_id}...")
                    domains = entry['domains']
                    for domain in domains:
                        f_domain.write(f"{organism_name}\t{uniprot_id}\t{domain['source_database']}\t{domain['accession']}\t{domain['name']}\t{domain['start']}\t{domain['end']}\t{domain['score'] if domain['score'] is not None else 'null'}\n")
                        domain_list, _ = proteins_domains[combined_id]
                        domain_list.append((domain['start'], domain['end'], domain['source_database'], domain['name']))
                        proteins_domains[combined_id] = (domain_list, domain['protein_length'])
                else:
                    # Fallback to fetching from InterPro API (for UniProt IDs)
                    print(f"Fetching domain information from InterPro API for {uniprot_id}...")
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
    parser = argparse.ArgumentParser(description="Process RefSeq accessions or FASTA files to generate structured protein data")
    
    parser.add_argument("-i", "--input", required=True,
                        help="Input file containing RefSeq accessions (one per line) or directory containing FASTA files")
    
    parser.add_argument("-o", "--output-dir", default="output",
                        help="Output directory for all generated files (default: 'output')")
    
    parser.add_argument("-p", "--prefix", default="sorted_refseq",
                        help="Prefix for all output files (default: 'sorted_refseq')")
    
    parser.add_argument("-m", "--mobidb-path", 
                        default="/home/ilnitsky/tools/MobiDB-lite/src/mobidb_lite/__main__.py",
                        help="Path to MobiDB-lite __main__.py script")
    
    parser.add_argument("-s", "--species-only", action="store_true",
                        help="Use only species name in FASTA headers")
    
    parser.add_argument("--iprscan-dir", default="iprscan_results",
                        help="Directory containing InterProScan TSV result files")
    
    parser.add_argument("--fold-res-dir", default="fold_res",
                        help="Directory containing folding results in the structure Genus_species/Genus_species/ranked_3.pdb")
    
    return parser.parse_args()

def parse_interproscan_tsv(file_path: str) -> Dict[str, List[Dict[str, Any]]]:
    """Parse InterProScan TSV output file for domain information.
    
    Format: 
    0: Protein accession
    1: MD5 hash
    2: Sequence length
    3: Analysis (source database)
    4: Signature accession
    5: Signature description
    6: Start location
    7: Stop location
    8: Score (E-value or other)
    9: Status
    10: Date
    11: InterPro accession
    12: InterPro description
    13: GO terms
    14: Pathway annotations
    
    Returns a dictionary mapping protein IDs to lists of domain dictionaries.
    """
    # Dictionary to store domains grouped by protein ID
    protein_domains = defaultdict(list)
    
    # Dictionary to track protein lengths
    protein_lengths = {}
    
    try:
        with open(file_path, 'r') as f:
            for line in f:
                # Skip empty lines
                if not line.strip():
                    continue
                    
                parts = line.strip().split('\t')
                
                if len(parts) < 11:
                    print(f"Warning: Skipping line with insufficient columns: {line}")
                    continue
                
                # Extract protein ID and length
                protein_id = parts[0]
                
                # Store protein length if available
                if parts[2].isdigit():
                    protein_length = int(parts[2])
                    protein_lengths[protein_id] = protein_length
                
                # Skip MobiDBLite entries as they're handled separately
                if parts[3] == "MobiDBLite":
                    continue
                
                # Map source database names to our supported types
                source_db = parts[3].lower()
                if source_db == "gene3d":
                    source_db = "cathgene3d"
                elif source_db == "superfamily":
                    source_db = "ssf"
                
                # Only process domains from supported databases
                if source_db in PRIORITY:
                    domain_info = {
                        'source_database': source_db,
                        'accession': parts[4],
                        'name': parts[5],
                        'start': int(parts[6]),
                        'end': int(parts[7]),
                        'score': float(parts[8]) if parts[8] != '-' else None,
                        'protein_length': protein_lengths.get(protein_id, 0)
                    }
                    protein_domains[protein_id].append(domain_info)
    
    except Exception as e:
        print(f"Error parsing InterProScan results from {file_path}: {e}")
        print(f"Line being processed: {line if 'line' in locals() else 'Unknown'}")
    
    return protein_domains

def get_protein_length(fasta_file: str, protein_id: str) -> int:
    """Get protein length from FASTA file."""
    with open(fasta_file, 'r') as f:
        current_id = None
        sequence = []
        
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id == protein_id:
                    break
                current_id = line[1:].split()[0]
                sequence = []
            elif current_id == protein_id:
                sequence.append(line)
        
        if sequence:
            return len(''.join(sequence))
    
    return 0

def main():
    """Main function with command-line argument support"""
    args = parse_arguments()
    
    # Check if input is a directory with FASTA files or a text file with accessions
    if os.path.isdir(args.input):
        # Original behavior: process FASTA files in directory
        fasta_files = glob.glob(os.path.join(args.input, "*.fasta"))
        all_results = {}
        
        for fasta_file in fasta_files:
            species_name = os.path.basename(fasta_file).replace(".fasta", "")
            print(f"Processing {species_name}...")
            
            # Look for corresponding InterProScan results for this species
            iprscan_file = os.path.join(args.iprscan_dir, f"{species_name}.fasta.tsv")
            
            # Parse the InterProScan file for all proteins in this species
            protein_domains = {}
            if os.path.exists(iprscan_file):
                print(f"Found InterProScan results at {iprscan_file}")
                protein_domains = parse_interproscan_tsv(iprscan_file)
            else:
                print(f"Warning: No InterProScan results found at {iprscan_file}")
            
            # Read FASTA file
            with open(fasta_file, 'r') as f:
                for record in SeqIO.parse(f, "fasta"):
                    protein_id = record.id
                    # Get taxonomy information directly from NCBI
                    tax_info = get_refseq_info(protein_id)
                    
                    # Get domains for this protein ID
                    domains = protein_domains.get(protein_id, [])
                    
                    # Add protein length to domain info if not already set
                    protein_length = len(record.seq)
                    for domain in domains:
                        if domain.get('protein_length', 0) == 0:
                            domain['protein_length'] = protein_length
                    
                    all_results[protein_id] = {
                        'fasta': str(record.format("fasta")),
                        'taxonomy_id': tax_info['taxonomy_id'],
                        'organism': tax_info['organism'],
                        'domains': domains
                    }
    
    else:
        # New behavior: process RefSeq accessions from text file
        print(f"Processing RefSeq accessions from {args.input}...")
        refseq_ids = read_uniprot_ids(args.input)  # Reuse this function for reading IDs
        all_results = {}
        
        # Look for a single InterProScan file for all accessions
        iprscan_file = os.path.join(args.iprscan_dir, "all_proteins.tsv")
        protein_domains = {}
        if os.path.exists(iprscan_file):
            print(f"Found InterProScan results at {iprscan_file}")
            protein_domains = parse_interproscan_tsv(iprscan_file)
        else:
            print(f"Warning: No InterProScan results found at {iprscan_file}")
        
        # Collect all InterProScan files and parse them
        protein_domains = {}
        iprscan_files = glob.glob(os.path.join(args.iprscan_dir, "*.tsv"))
        
        for iprscan_file in iprscan_files:
            print(f"Found InterProScan results at {iprscan_file}")
            file_domains = parse_interproscan_tsv(iprscan_file)
            protein_domains.update(file_domains)
        
        if not protein_domains:
            print(f"Warning: No InterProScan results found in {args.iprscan_dir}")
        
        for refseq_id in refseq_ids:
            print(f"Processing {refseq_id}...")
            
            # Get taxonomy information directly from NCBI
            tax_info = get_refseq_info(refseq_id)
            
            # Get FASTA sequence from NCBI
            try:
                esearch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term={refseq_id}"
                response = requests.get(esearch_url)
                if response.status_code == 200:
                    protein_id_match = re.search(r'<Id>(\d+)</Id>', response.text)
                    if protein_id_match:
                        protein_id = protein_id_match.group(1)
                        
                        # Get FASTA sequence
                        efetch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={protein_id}&rettype=fasta&retmode=text"
                        fasta_response = requests.get(efetch_url)
                        if fasta_response.status_code == 200:
                            fasta_text = fasta_response.text
                            
                            # Get domains for this protein ID
                            domains = protein_domains.get(refseq_id, [])
                            
                            # Calculate protein length from FASTA
                            sequence_lines = fasta_text.split('\n')[1:]
                            sequence = ''.join(line.strip() for line in sequence_lines if line.strip())
                            protein_length = len(sequence)
                            
                            # Add protein length to domain info if not already set
                            for domain in domains:
                                if domain.get('protein_length', 0) == 0:
                                    domain['protein_length'] = protein_length
                            
                            all_results[refseq_id] = {
                                'fasta': fasta_text,
                                'taxonomy_id': tax_info['taxonomy_id'],
                                'organism': tax_info['organism'],
                                'domains': domains
                            }
                            
                            time.sleep(0.2)  # Be nice to NCBI servers
                        else:
                            print(f"Failed to get FASTA for {refseq_id}")
                    else:
                        print(f"No protein ID found for {refseq_id}")
                else:
                    print(f"Failed to search for {refseq_id}")
                    
            except Exception as e:
                print(f"Error processing {refseq_id}: {e}")
                continue
    
    # Continue with existing processing...
    sorted_data = sort_by_taxonomy(all_results)
    
    print(f"Saving results to {args.output_dir}...")
    accession_mapping = create_accession_to_folder_mapping(args.iprscan_dir)
    output_files = save_results(
        sorted_data, 
        output_dir=args.output_dir,
        output_prefix=args.prefix,
        species_only=args.species_only,
        fold_res_dir=args.fold_res_dir,
        accession_mapping=accession_mapping
    )
    
    # Run MobiDB-lite and create annotation files
    mobidb_output = run_mobidb_and_create_annotation(
        mobidb_path=args.mobidb_path,
        fasta_file=output_files['fasta'],
        output_dir=args.output_dir,
        ss_file=output_files['ss'],
        output_prefix=args.prefix
    )
    
    print("Processing completed successfully!")
    print(f"Output files:")
    for file_type, file_path in output_files.items():
        print(f"  - {file_type}: {file_path}")
    if mobidb_output and all(mobidb_output):
        full_annotation, mobidb_itol = mobidb_output
        print(f"  - MobiDB annotation: {full_annotation}")
        print(f"  - MobiDB iTOL: {mobidb_itol}")

if __name__ == "__main__":
    main()