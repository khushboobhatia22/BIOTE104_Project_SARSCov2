from Bio import SeqIO
import random

fasta_files = ['B.1.1.529_cds.fna', 'B.1.351_cds.fna',	'B.1.525_cds.fna', 'B.1.617.1_cds.fna',	'C.37_cds.fna',	'P.2_cds.fna', 'B.1.1.7_cds.fna', 'B.1.427_cds.fna', 'B.1.526_cds.fna',	'B.1.617.2_cds.fna', 'P.1_cds.fna',	'P.3_cds.fna', 'JN.1.11.1_cds.fna', 'KP.3.1.1_cds.fna', 'XBB.1.5_cds.fna', 'XEC_cds.fna']

path = '' #Path to Fasta files

def select_random_fasta_records(input_fasta, output_fasta, count=2000):
    """
    Select random records from a FASTA file and save them to a new file.

    Args:
    - input_fasta (str): Path to the input FASTA file.
    - output_fasta (str): Path to save the selected random records.
    - count (int): Number of random records to select.

    Returns:
    - None
    """
    # Read all records from the input FASTA file
    records = list(SeqIO.parse(input_fasta, "fasta"))
    
    # Check if the input file has enough records
    if len(records) > count:
        random_records = random.sample(records, count)
    else:
        random_records = records
    
    # Write the selected records to the output FASTA file
    SeqIO.write(random_records, output_fasta, "fasta")
    print(f"{count} random records written to {output_fasta}.")


for fasta_file in fasta_files:
    input_fasta = path+fasta_file
    output_fasta = path+fasta_file+'_trimmed'
    try:
        select_random_fasta_records(input_fasta, output_fasta, count=2000)
    except Exception as e:
        print(f"Error: {e}")
