import pandas as pd
import subprocess
import os
df = pd.read_csv('/home/ilnitsky/NPM/Ecdysozoa.txt', sep="\t")
df = df.drop_duplicates(subset='Organism Name', keep='first')

def ensure_directory_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
        print(f"Directory created: {directory}")
    else:
        print(f"Directory already exists: {directory}")

base_dir = "/home/ilnitsky/NPM/genomes/ecdysozoa"
results_dir = "/home/ilnitsky/NPM/results_ecdysozoa"

query_proteins = "/home/ilnitsky/NPM/all_proteins/N_term/N_term_Arthropoda_subsample.fasta"
c_term_query_proteins = "/home/ilnitsky/NPM/all_proteins/C_term/C_term_Echinodermata.fasta"


subdirs = ['regions', 'extracted', 'nterm_blast_results', 'cterm_blast_results', 'cterm_psiblast']
for subdir in subdirs:
    ensure_directory_exists(os.path.join(results_dir, subdir))

ensure_directory_exists(base_dir)
ensure_directory_exists(results_dir)

for index, row in df.iterrows():
    try:
        assembly_accesion = row['Assembly Accession']
        assembly_name = row['Assembly Name']
        organism_name = row['Organism Name'].replace(' ', '_')
        genome_fna = f"{base_dir}/{organism_name}/ncbi_dataset/data/{assembly_accesion}/{assembly_accesion}_{assembly_name}_genomic.fna"
        
        print(f"\nProcessing organism: {organism_name}")
        print(f"Assembly name: {assembly_name}")
        print(f"Genome file path: {genome_fna}")

        # Check if organism directory exists
        organism_dir = f"{base_dir}/{organism_name}"
        if not os.path.exists(organism_dir):
            print(f"Creating directory for organism: {organism_dir}")
            os.makedirs(organism_dir)
            
            curl_command = f"""
                curl -o {base_dir}/{organism_name}.zip 'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{assembly_accesion}/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED'
                unzip {base_dir}/{organism_name}.zip -d {base_dir}/{organism_name}
            """
            print(f"Downloading and extracting genome for: {organism_name}")
            subprocess.run(curl_command, shell=True, check=True)

        # Check if genome file exists
        if not os.path.exists(genome_fna):
            raise FileNotFoundError(f"Genome file not found: {genome_fna}")

        make_indexes_command = f"""
            set -x 
            makeblastdb -in {genome_fna} -dbtype nucl -parse_seqids || echo "makeblastdb failed with status $?"
            samtools faidx {genome_fna} || echo "samtools faidx failed with status $?"
        """
        print("Making indexes...")
        try:
            result = subprocess.run(make_indexes_command, shell=True, check=True, capture_output=True, text=True)
            print("STDOUT:", result.stdout)
            print("STDERR:", result.stderr)
        except subprocess.CalledProcessError as e:
            print(f"Command failed with exit status {e.returncode}")
            print("STDOUT:", e.stdout)
            print("STDERR:", e.stderr)


        blast_command = f"""
            tblastn -query {query_proteins} \
                -db {genome_fna}  \
                -outfmt 6 \
                -evalue 0.0001 \
            | tee -a {results_dir}/nterm_blast_results/{organism_name}_blast_results.txt \
            | awk 'BEGIN{{OFS="\\t"}} {{ if ($9 < $10) {{start = $9 - 10000; end = $10 + 50000; strand="+"}} \
                else {{start = $10 - 10000; end = $9 + 50000; strand="-"}} \
                if (start < 0) start = 0;  print $2, start, end, strand;}}' \
                > {results_dir}/regions/{organism_name}.bed
            """
        
        subprocess.run(blast_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("Running command:", blast_command)

        get_strand = f"awk 'NR==1 {{print $4}}' {results_dir}/regions/{organism_name}.bed"

        result = subprocess.run(get_strand, shell=True, capture_output=True, text=True)
        print("Running command:", get_strand)

        strand = "R" if result.stdout.strip() == "-" else "F"

        extract_fasta_command = f"""
            bedtools getfasta -fi {genome_fna} -bed {results_dir}/regions/{organism_name}.bed -fo {results_dir}/extracted/{organism_name}.extracted.fasta
            """

        subprocess.run(extract_fasta_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("Running command:", extract_fasta_command)
        

        translate_in_3_frames_command = f"""
            transeq {results_dir}/extracted/{organism_name}.extracted.fasta {results_dir}/extracted/{organism_name}.aminoacid.extracted.fasta  -frame={strand}
        """

        subprocess.run(translate_in_3_frames_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("Running command:", translate_in_3_frames_command)

        makeblastdb_extracted_command = f"""
            makeblastdb -in {results_dir}/extracted/{organism_name}.aminoacid.extracted.fasta -dbtype prot
            """

        subprocess.run(makeblastdb_extracted_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("Running command:", makeblastdb_extracted_command)

        c_terminal_blast = f"""
            blastp -query {c_term_query_proteins} \
                -db {results_dir}/extracted/{organism_name}.aminoacid.extracted.fasta  \
                -outfmt 6 \
                -evalue 0.001 \
                > {results_dir}/cterm_blast_results/{organism_name}.txt
        """
        subprocess.run(c_terminal_blast, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("Running command:", c_terminal_blast)

        #PSI-BLAST search
        psi_blast_command = f"""
            psiblast -query {c_term_query_proteins} \
                -db {results_dir}/extracted/{organism_name}.aminoacid.extracted.fasta \
                -num_iterations 3 \
                -outfmt 6 \
                -evalue 0.001 \
                -out {results_dir}/cterm_psiblast/{organism_name}.txt
        """
        subprocess.run(psi_blast_command, shell=True, check=True)
        print("Running command:", psi_blast_command)

    except Exception as e:
        print(f"Error processing {organism_name}: {str(e)}")
        continue
    
    # augustus_CDS_command = f"""
    #     augustus --species=amphimedon {results_dir}/extracted/{organism_name}.extracted.fasta > {results_dir}/augustus/{organism_name}.augustus --AUGUSTUS_CONFIG_PATH=/home/ilnitsky/anaconda3/envs/npm/config/
    #     """

    # hmmer_command = f"""
    #     hmmsearch --domtblout {results_dir}/hmmer/{organism_name}.C_term.txt ~/gpfs/NPM/ncbi_blast_search/N_term_Cnidaria.hmm {results_dir}/extracted/{organism_name}.aminoacid.extracted.fasta > {results_dir}/hmmer/NPM_ALL_hmm_search.{organism_name}.txt
    # """
    # rm_command = f"rm {genome_dir}/{organism_name}.*"





