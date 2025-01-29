"""
Concatenate output alignments from Seq-Gen and generate a single FASTA file.

This script processes the output alignments produced by Seq-Gen, where the outputs 
are organized as gene1, gene2, ..., geneM.nex files within a directory. The script concatenates all gene files within each rep folder and produces a single concatenated FASTA file for each rep directory. RepID is provided as one argument in the final output. 

The script handles missing gene.nexus files as follows:
- If hidden paralogy is simulated, some geneM.nex files may be missing. The script will ignore such cases. 
- If gene loss happens and one accession/species/taxa got lost, the script wll append "-" x length of the sequence to substitute the loss with empty sequence. 

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
This above code will concatenate nexus files within folder rep1 and rep2. It creates two output files concatenated_alignment_rep1.fasta and concatenated_alignment_rep2.fasta in the specific output directory. Each fasta is the concatenated sequences from each rep folder.  
"""

from Bio import SeqIO
import os
import sys
from collections import defaultdict

def concatenate_nexus_to_fasta(nexus_dir, fasta_dir, rep_id):
    """
    Concatenate sequences from multiple Nexus files in a directory into a single FASTA file.
    Paras:
    nexus_dir (str): Path to input Nexus files.
    fasta_dir (str): Path to output FASTA files.
    rep_id (char/int): Rep ID in the nexus file name. The nexus files are saved in rep1, rep2...
    Return: The output file is called "concated_alignment_rep$ID.fasta"
    """
    concate_sequences = defaultdict(str) # If a key doesn't exist, it initializes the key
    all_taxa = set() # Keep track of all taxa -- This will help with handle missing taxa in one gene tree files

    # Get a set for all taxa 
    for file_name in sorted(os.listdir(nexus_dir)):
        if file_name.endswith(".nex"):  # loop through .nex files in the directory 
            nexus_path = os.path.join(nexus_dir, file_name) # Input nexus files 
            with open(nexus_path, "r") as nexus_file: 
                for seq_record in SeqIO.parse(nexus_file, "nexus"):
                    all_taxa.add(seq_record.id) # record all taxa across all .nexus files 

    # Add seq from a particular nexus file into a dictionary. Missing taxa in this .nexus will be processed with the "-"
    for file_name in sorted(os.listdir(nexus_dir)): 
        if file_name.endswith(".nex"):  
            nexus_path = os.path.join(nexus_dir, file_name) 
            current_taxa = set() # set for the current taxa for a particular nexus file
            seq_length = None 
            
            with open(nexus_path, "r") as nexus_file: 
                for seq_record in SeqIO.parse(nexus_file, "nexus"):
                    seq_length = len(seq_record.seq)
                    current_taxa.add(seq_record.id) # add current taxa into the set 
                    concate_sequences[seq_record.id] += str(seq_record.seq) # Concatenate all taxa present in the file 
                
                # Add "-" * seq_length string to missing taxa in this file: 
                missing_taxa = all_taxa - current_taxa # Check if there is missing betweeb current taxa and all taxa 
                for taxon in missing_taxa: 
                    print(f"Notice: missing taxa {taxon} within the file {file_name}. Let's add an empty string!")
                    concate_sequences[taxon] += "-" * seq_length # Set the length to be the same as the last sequence 

     # Write the dic to a fasta files: 
    fasta_path = os.path.join(fasta_dir, f"concate_alignment_rep{rep_id}.fasta") # Output fasta files: The output file is within the assigned output folder named as concatenated_alignment_rep$id.fasta
    with open(fasta_path, "w") as fasta_file: 
        for taxon, sequence in concate_sequences.items():
            fasta_file.write(f">{taxon}\n{sequence}\n")

if __name__ == "__main__":
    nexus_dir = sys.argv[1] # First argument is the input nexus folder
    fasta_dir = sys.argv[2] # Second argument is the folder to store the output fasta 
    rep_id = sys.argv[3] # replicate ID
    concatenate_nexus_to_fasta(nexus_dir, fasta_dir, rep_id)
    print("Aww! Concatenated nexus files into fasta!")
