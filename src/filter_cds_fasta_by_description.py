from Bio import SeqIO

# Specify the path to your FASTA file
fasta_files = ['B.1.1.529_cds.fna', 'B.1.351_cds.fna',	'B.1.525_cds.fna', 'B.1.617.1_cds.fna',	'C.37_cds.fna',	'P.2_cds.fna', 'B.1.1.7_cds.fna', 'B.1.427_cds.fna', 'B.1.526_cds.fna',	'B.1.617.2_cds.fna', 'P.1_cds.fna',	'P.3_cds.fna', 'JN.1.11.1_cds.fna', 'KP.3.1.1_cds.fna', 'XBB.1.5_cds.fna', 'XEC_cds.fna']
path = '' #Path to above files

result = []
for fasta_file in fasta_files:
    # Parse the FASTA file and iterate through the sequences
    for record in SeqIO.parse(path+fasta_file, "fasta"):
        if 'surface glycoprotein' in record.description :
            result.append(record)

    outfile = path+fasta_file+'_filtered'
    with open(outfile, 'w') as outhandle:
        SeqIO.write(result, outhandle, "fasta")

