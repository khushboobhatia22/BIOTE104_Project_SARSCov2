from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
from Bio.Align import AlignInfo
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

fasta_files = ['B.1.1.7_cds.fna', 'B.1.1.529_cds.fna', 'B.1.351_cds.fna',	'B.1.525_cds.fna', 'B.1.617.1_cds.fna',	'C.37_cds.fna',	'P.2_cds.fna', 'B.1.427_cds.fna', 'B.1.526_cds.fna',	'B.1.617.2_cds.fna', 'P.1_cds.fna',	'P.3_cds.fna', 'JN.1.11.1_cds.fna', 'KP.3.1.1_cds.fna', 'XBB.1.5_cds.fna', 'XEC_cds.fna']
path = '' #Path to fasta_files. Same path will be used for output file
outfile = 'cds_s_consensus.fna'

mafft_path = '' #Path to MAFFT executable

def multiple_alignment_with_consensus(input_fasta, output_aln, mafft_path):
    """
    Perform multiple sequence alignment using MAFFT and calculate the consensus sequence.

    Args:
    - input_fasta (str): Path to the input FASTA file containing genome sequences.
    - output_aln (str): Path to save the aligned output in FASTA format.
    - mafft_path (str): Path to the MAFFT executable (default: assumes `mafft` is in PATH).

    Returns:
    - alignment: A MultipleSeqAlignment object containing the aligned sequences.
    - consensus: The consensus sequence as a string.
    """
    # Ensure the input file exists
    if not os.path.exists(input_fasta):
        raise FileNotFoundError(f"Input file {input_fasta} does not exist.")
    
    # Run MAFFT using Biopython's wrapper
    mafft_cline = MafftCommandline(mafft_path, input=input_fasta)
    stdout, stderr = mafft_cline()
    
    # Save the alignment to a file
    with open(output_aln, "w") as output_file:
        output_file.write(stdout)
    
    # Parse the alignment
    alignment = AlignIO.read(output_aln, "fasta")
    
    # Calculate the consensus sequence
    summary_info = AlignInfo.SummaryInfo(alignment)
    consensus = summary_info.dumb_consensus(threshold=0.3)
    
    return alignment, str(consensus)

try:
    out_sequences = []
    if mafft_path == '' :
        print('ERROR: MAFFT Path is empty. Please provide path to MAFFT executable')
    
    for fasta_file in fasta_files:
        input_fasta = path+fasta_file+'_trimmed'
        output_aln = path+fasta_file+'_trimmed_aligned'
        alignment, consensus = multiple_alignment_with_consensus(input_fasta, output_aln, mafft_path)
        print(f"Alignment completed. Saved to {output_aln}.")
        print("Consensus Sequence:")
        print(consensus)
        variant = fasta_file.split("_")[0]
        seq = SeqRecord(Seq(str(consensus)), id=variant, description="")
        out_sequences.append(seq)

    with open(path+outfile, "w") as output_handle:
        SeqIO.write(out_sequences, output_handle, "fasta")
except Exception as e:
    print(f"Error: {e}")
