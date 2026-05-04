"""
Concatenate output alignments from Seq-Gen and generate a single FASTA.

This script processes the output alignments produced by Seq-Gen, 
where the outputs are organized as gene1, gene2, ..., geneM.nex files. 
The script concatenates all gene files within each rep folder 
and produces a single concatenated FASTA file for each rep directory. 
RepID is provided as one argument in the final output. 

The script handles missing gene.nexus files as follows:
- If hidden paralogy is simulated, some geneM.nex files may be missing. 
    The script will ignore such cases. 
- If gene loss happens and one accession/species/taxa got lost, 
    the script wll append "-" x length of the sequence to substitute 
    the loss with empty sequence. 

### Example Usage:
A test kit is saved in `$ROOT/example/test_concatenate_seq`.
Run the following shell commands in `$ROOT`:
#!/bin/bash
cd ~ `$ROOT/example/` 
for i in 1 2; do
    python ../scripts/concatenate_seq.py \
        test_concatenate_seq/rep$i \
        test_concatenate_seq/ \
        $i
done
### 
    This above code will concatenate nexus files within folder rep1 and rep2. 
    It creates two output files concatenated_alignment_rep1.fasta 
    and concatenated_alignment_rep2.fasta in the specific output directory. 
    Each fasta is the concatenated sequences from each rep folder.  
"""

from Bio import SeqIO
import os
import sys
import re
from collections import defaultdict

def concatenate_nexus_to_fasta(nexus_dir, fasta_dir, rep_id):
    """
    Concatenate sequences from multiple Nexus in a directory into a FASTA file.
    Paras:
    nexus_dir (str): Path to input Nexus files.
    fasta_dir (str): Path to output FASTA files.
    rep_id (char/int): Rep ID in the nexus file name. 
    Return: The output file is called "concated_alignment_rep$ID.fasta"
    """
    # If a key doesn't exist, it initializes the key 
    concate_sequences = defaultdict(str) 
    # Keep track of all taxa -> help with handle missing taxa in a gene tree file
    all_taxa = set() 

    # Get a set for all taxa 
    for file_name in sorted(os.listdir(nexus_dir)):
        if file_name.endswith(".nex"):  # loop through .nex files in the directory 
            nexus_path = os.path.join(nexus_dir, file_name) # Input nexus files 
            with open(nexus_path, "r") as nexus_file: 
                for seq_record in SeqIO.parse(nexus_file, "nexus"):
                    all_taxa.add(seq_record.id) # record all taxa across all .nexus files 

    # Add seq from a particular nexus file into a dictionary.
    # Missing taxa in this .nexus will be processed with the "-"
    for file_name in sorted(os.listdir(nexus_dir)): 
        if file_name.endswith(".nex"):  
            nexus_path = os.path.join(nexus_dir, file_name) 
            current_taxa = set() # set for the current taxa for a particular nexus file
            seq_length = None 
            
            with open(nexus_path, "r") as nexus_file: 
                for seq_record in SeqIO.parse(nexus_file, "nexus"):
                    seq_length = len(seq_record.seq)
                    current_taxa.add(seq_record.id) # add current taxa into the set 
                    # Concatenate all taxa present in the file: 
                    concate_sequences[seq_record.id] += str(seq_record.seq) 
                
                # Add "-" * seq_length string to missing taxa in this file: 
                # Check if there is missing between current taxa and all taxa 
                missing_taxa = all_taxa - current_taxa 
                for taxon in missing_taxa: 
                    # Set the length to be the same as the last sequence: 
                    concate_sequences[taxon] += "-" * seq_length 

    # Reorder sequences alphabetically by letter first, then by number
    reordered_sequences, reorder_message = reorder_sequences_alphabetically(concate_sequences)
    
    # Write the reordered dic to a fasta files: 
    # Output are in the output folder named as concatenated_alignment_rep$id.fasta
    fasta_path = os.path.join(fasta_dir, f"concate_alignment_rep{rep_id}.fasta") 
    with open(fasta_path, "w") as fasta_file: 
        for taxon, sequence in reordered_sequences.items():
            fasta_file.write(f">{taxon}\n{sequence}\n")
    
    # Return the message for the calling script
    return reorder_message

def reorder_sequences_alphabetically(sequences_dict):
    """
    Reorder the sequences dictionary alphabetically by letter first, then by number.
    For example: A_0, A_1, B_0, B_1, B_2, C_0, etc.
    Args:
        sequences_dict (dict): Dictionary with taxon names as keys and sequences as values
    Returns:
        tuple: (reordered_dict, message) - Reordered dictionary and status message
    """
    message = ""  # Initialize message string
    
    def sort_key(taxon_name):
        """
        Create a sort key that separates letter and number parts.
        Returns a tuple (letter_part, number_part) for proper sorting.
        """
        # Match pattern like A_0, B_1, etc.
        match = re.match(r'^([A-Za-z]+)_?(\d+)?', taxon_name)
        if match:
            letter_part = match.group(1)
            number_part = int(match.group(2)) if match.group(2) else 0
            return (letter_part, number_part)
        else:
            # If pattern doesn't match, sort by the original name
            return (taxon_name, 0)
    
    # Sort all taxa alphabetically by letter first, then by number
    sorted_taxa = sorted(sequences_dict.keys(), key=sort_key)
    
    reordered_dict = {}
    for taxon in sorted_taxa:
        reordered_dict[taxon] = sequences_dict[taxon]
    
    # Create a message if there is no "A_" in the keys 
    if not any(taxon.startswith("A_") for taxon in sorted_taxa):
        message = "No 'A_' taxon found. No outgroup for following-up analyses. " 
    else: 
        pass # If trees are starting with "A_", no message is needed 

    return reordered_dict, message

if __name__ == "__main__":
    # First argument is the input nexus folder: 
    nexus_dir = sys.argv[1] 
    # Second argument is the folder to store the output fasta: 
    fasta_dir = sys.argv[2] 
    rep_id = sys.argv[3] # replicate ID
    
    # Call the function and get the message
    message = concatenate_nexus_to_fasta(nexus_dir, fasta_dir, rep_id)
    
    # Print only the message (this will be captured by Julia)
    print(message)

