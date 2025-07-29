#!/usr/bin/env python3
"""
Script to create a proper multiple sequence alignment by combining
full sequences with structural alignment information.
"""

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import re

def read_fasta_file(filepath):
    """Read a FASTA file and return the sequence and header."""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    header = lines[0].strip()
    # Join sequence lines
    sequence = ''.join(line.strip() for line in lines[1:])
    return header, sequence

def get_uniprot_id_from_filename(filename):
    """Extract UniProt ID from filename."""
    match = re.search(r'([A-Z0-9]+)\.fasta', filename)
    return match.group(1) if match else None

def create_proper_msa():
    """Create a proper MSA with full sequences and structural alignment gaps."""
    
    # Directory containing the FASTA files
    fasta_dir = "aligned_structures"
    
    # Read all FASTA files
    fasta_files = [f for f in os.listdir(fasta_dir) if f.endswith('.fasta') and not f.startswith('aligned_')]
    
    # Dictionary to store full sequences
    full_sequences = {}
    
    for fasta_file in fasta_files:
        uniprot_id = get_uniprot_id_from_filename(fasta_file)
        if uniprot_id:
            filepath = os.path.join(fasta_dir, fasta_file)
            header, sequence = read_fasta_file(filepath)
            full_sequences[uniprot_id] = {"header": header, "sequence": sequence}
            print(f"{uniprot_id}: {len(sequence)} residues")
    
    # Read the structural alignment to understand the common region
    structural_msa_path = os.path.join(fasta_dir, "aligned_sequences_structural.fasta")
    structural_msa = AlignIO.read(structural_msa_path, "fasta")
    
    print(f"\nStructural alignment length: {structural_msa.get_alignment_length()}")
    
    # Create a mapping from structural alignment to full sequences
    # We need to identify where the structural alignment fits in each full sequence
    
    def find_best_match(structural_seq, full_seq):
        """Find the best match of structural sequence within full sequence."""
        structural_seq_str = str(structural_seq.seq).replace('-', '')  # Remove gaps
        full_seq_str = str(full_seq)
        
        # Simple sliding window approach
        best_score = 0
        best_start = 0
        
        for i in range(len(full_seq_str) - len(structural_seq_str) + 1):
            window = full_seq_str[i:i+len(structural_seq_str)]
            score = sum(1 for a, b in zip(window, structural_seq_str) if a == b and a != '-')
            if score > best_score:
                best_score = score
                best_start = i
        
        return best_start, best_score
    
    # Create proper MSA
    msa_records = []
    
    for record in structural_msa:
        # Extract UniProt ID from the record ID (structure ID in format "E9FX34 Structure 1")
        uniprot_id = record.id.split()[0]  # Get the first part before the space
        
        if uniprot_id in full_sequences:
            full_seq = full_sequences[uniprot_id]["sequence"]
            
            # Find where the structural alignment fits in the full sequence
            # We need to compare without gaps
            struct_seq_nogaps = str(record.seq).replace('-', '')
            
            # Find the start position in the full sequence
            start_pos, score = find_best_match(record, full_seq)
            
            print(f"{uniprot_id}: structural alignment starts at position {start_pos} (score: {score})")
            
            # Create aligned sequence with gaps
            aligned_seq = []
            
            # Add the portion of the full sequence before the structural alignment
            if start_pos > 0:
                aligned_seq.extend(full_seq[:start_pos])
            
            # Add the structural alignment segment
            # We need to track position in the full sequence as we process the structural alignment
            full_seq_pos = start_pos
            for aa in record.seq:
                if aa != '-':  # If not a gap in structural alignment
                    aligned_seq.append(full_seq[full_seq_pos])
                    full_seq_pos += 1
                else:  # If it's a gap in structural alignment
                    aligned_seq.append('-')
            
            # Add the remainder of the full sequence
            if full_seq_pos < len(full_seq):
                aligned_seq.extend(full_seq[full_seq_pos:])
            
            # Create SeqRecord
            seq_record = SeqRecord(
                Seq(''.join(aligned_seq)),
                id=uniprot_id,
                description=f"Full sequence aligned - {len(full_seq)} residues"
            )
            msa_records.append(seq_record)
    
    # Create MSA
    msa = MultipleSeqAlignment(msa_records)
    
    # Save the proper MSA
    output_path = os.path.join(fasta_dir, "proper_aligned_sequences.fasta")
    AlignIO.write(msa, output_path, "fasta")
    
    print(f"\nProper MSA with gaps saved to {output_path}")
    print(f"Alignment length: {msa.get_alignment_length()}")
    
    # Print the alignment
    print("\nProper Multiple Sequence Alignment:")
    for record in msa:
        print(f">{record.id} {record.description}")
        # Print sequence in chunks of 60 characters
        seq_str = str(record.seq)
        for i in range(0, len(seq_str), 60):
            print(seq_str[i:i+60])
        print()

if __name__ == "__main__":
    create_proper_msa() 