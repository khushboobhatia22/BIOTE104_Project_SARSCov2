from pyfaidx import Fasta
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import json

va_file = 'variant_accessionjson'
np_file = 'nucleotide_protein_map.json'
fasta_file = 'protein.faa'  # Input FASTA file
path = '' #Path to above files. Same will be used for output file

# Function to load, trim headers, and save to a new FASTA file
def filter_fasta(fasta_file, va_file, np_file, outpath):
    variant_data = None
    np_data = None
    acc_dict = dict()

    np_dict = dict()
    with open(va_file, 'r') as file:
        variant_data = json.load(file)

    for variant in variant_data:
        for acc in variant_data[variant]:
            acc_dict[acc] = variant

    with open(np_file, 'r') as file:
        np_data = json.load(file)

    for n in np_data:
        p = np_data[n]['prot_acc']
        np_dict[p] = n

    # Load the FASTA file using pyfaidx
    fasta_data = Fasta(fasta_file)
    
    trimmed_records = dict()
    
    # Iterate through each sequence in the FASTA file
    for header, sequence in fasta_data.items():
        # Trim the header to keep the part before the colon
        trimmed_header = header.split(":")[0]
        trimmed_header = trimmed_header.strip()
        #print(trimmed_header)
        #print(np_dict)

        if trimmed_header in np_dict:
            nucl = np_dict[trimmed_header]
            acc_variant = acc_dict[nucl]

            if acc_variant not in trimmed_records:
                trimmed_records[acc_variant] = []
        
            # Create a SeqRecord with the trimmed header and sequence
            seq_record = SeqRecord(Seq(str(sequence)), id=header, description="")

            trimmed_records[acc_variant].append(seq_record)
            #print(acc_variant, header, sequence)
            #print(seq_record)

    for variant in trimmed_records:
        #print(variant, trimmed_records[variant])
        filename = outpath+variant+'_protein.faa'
        # Write the trimmed sequences to the output FASTA file using Biopython
        with open(filename, "w") as output_handle:
            SeqIO.write(trimmed_records[variant], output_handle, "fasta")


filter_fasta(path+fasta_file, path+va_file, path+np_file, path)
