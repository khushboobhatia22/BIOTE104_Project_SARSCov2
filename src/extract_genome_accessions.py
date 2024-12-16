import json

# Path to the file containing JSON objects on each line
#file_path = '/Users/khushboobhatia/bioinformatics/project/ncbi/protein_annotation/ncbi_dataset/data/data_report.jsonl'
file_path = 'variant_accession.json'
path = ''  #Path to above file

'''
B.1.1.7 (Alpha), B.1.351 (Beta), B.1.525 (Eta), B.1.427/B.1.429 (Epsilon), B.1.526 (Iota), B.1.617.1 (Kappa), B.1.617.2 (Delta), C.37 (Lambda), P.1 (Gamma), P.2 (Zeta), P.3 (Theta), B.1.1.529 (Omicron)
'''

variants_interest = ['B.1.1.7', 'B.1.351', 'B.1.525', 'B.1.427', 'B.1.526', 'B.1.617.1', 'B.1.617.2', 'C.37', 'P.1', 'P.2', 'P.3', 'B.1.1.529', 'XBB.1.5', 'JN.1.11.1', 'KP.3.1.1', 'XEC']

# Open the file and parse each line as a separate JSON object
with open(file_path, 'r') as file:
    data = json.load(file)
    for variant in data.keys():
        with open(path+variant+'_n_accession', 'w') as outfile:
            for accession in data[variant]:
                outfile.write(accession)
                outfile.write('\n')
            