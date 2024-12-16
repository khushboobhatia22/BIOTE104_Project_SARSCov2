from Bio import SeqIO
import math

fasta_file = 'cds_s_consensus_aligned.fna'
path = '' #Path to fasta_file

def find_entropy_by_position(input_fasta):

    # Read all records from the input FASTA file
    records = list(SeqIO.parse(input_fasta, "fasta"))
    length = len(records[0].seq)
    entropy = [0]*length
    for i in range(length):
        n_count = dict()
        for record in records:
            base =  record.seq[i].upper()
            if base not in n_count:
                n_count[base] = 0
            n_count[base] += 1
        entropy_base = 0
        for base in n_count:
            base_probability = n_count[base] / len(records)
            entropy_base += (-1)*base_probability * math.log(base_probability, 2)
        entropy[i] = entropy_base

    for i in range(len(entropy)):
        pos_set = set()
        for record in records:
            pos_set.add(record.seq[i])
        print(i+1, entropy[i], pos_set)   
    
try:
    input_fasta = path+fasta_file
    find_entropy_by_position(input_fasta)
except Exception as e:
    print(f"Error: {e}")


