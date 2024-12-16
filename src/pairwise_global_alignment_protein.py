from Bio.Align import PairwiseAligner
from Bio import SeqIO
from Bio.Align.substitution_matrices import load
import json

base_path = '' #Path for fasta_files. Same will be used for output files
fasta_files = ['B.1.1.529_protein.faa', 'B.1.351_protein.faa', 'B.1.525_protein.faa', 'B.1.617.1_protein.faa', 'C.37_protein.faa', 'P.2_protein.faa', \
    'B.1.1.7_protein.faa', 'B.1.427_protein.faa', 'B.1.526_protein.faa', 'B.1.617.2_protein.faa', 'P.1_protein.faa', 'P.3_protein.faa', \
    'JN.1.11.1_protein.faa', 'KP.3.1.1_protein.faa', 'XEC_protein.faa', 'XBB.1.5_protein.faa']

ref_path = '' #Path to reference sequence files
ref_file = ref_path+'protein.faa'

protien_name = 'surface glycoprotein'

def pairwise_global_alignment(seq1, seq2, mutation_dict):

    # Perform global alignment (Needleman-Wunsch algorithm)
    aligner = PairwiseAligner()
    blosum62 = load("BLOSUM62")
    aligner.substitution_matrix = blosum62
    aligner.mode = 'global'
    # Set gap penalties for the query (source) sequence
    aligner.query_open_gap_score = -2  # Query gap opening penalty
    aligner.query_extend_gap_score = -0.5  # Query gap extension penalty

    # Set gap penalties for the reference (target) sequence
    aligner.target_open_gap_score = -10    # Reference gap opening penalty
    aligner.target_extend_gap_score = -10   # Reference gap extension penalty
    try:
        alignments = aligner.align(seq1, seq2)
    except:
        print("Error")
        return
    res1 = None
    res2 = None

    for alignment in alignments:
        res1 = alignment[0]
        res2 = alignment[1]
        break
    for i in range(len(res1)):
        if res1[i] != res2[i]:           
            mut = res1[i]+str(i+1)+res2[i]
            if mut not in mutation_dict:
                mutation_dict[mut] = 0
            mutation_dict[mut] += 1
            
                

def get_ref_seq(fasta_file, protein):
    result = None
    for record in SeqIO.parse(fasta_file, "fasta"):
        if protein in record.description :
            result = record.seq
            break

    return result

def get_seq(fasta_file, protein):
    result = []
    for record in SeqIO.parse(fasta_file, "fasta"):
            result.append(record.seq)
    return result

ref_seq = get_ref_seq(ref_file, protien_name)
for fasta_file in fasta_files:
    mutation_dict = dict()
    seqs = get_seq(base_path+fasta_file, protien_name)
    n_seqs = len(seqs)
    print(n_seqs)
    for seq in seqs:
        pairwise_global_alignment(ref_seq,seq, mutation_dict)

    mutation_dict = {key: (value / n_seqs) * 100 for key, value in mutation_dict.items()}
    mutation_dict = dict(sorted(mutation_dict.items(), key=lambda item: item[1], reverse=True))

    with open(base_path+fasta_file+'_mutations', 'w') as outhandle:
        json.dump(mutation_dict, outhandle, indent=4)

    del mutation_dict

