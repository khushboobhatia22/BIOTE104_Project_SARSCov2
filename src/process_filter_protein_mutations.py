import json

protein_mut_files = ['B.1.1.529_protein.faa_mutations',\
'B.1.1.7_protein.faa_mutations',\
'B.1.351_protein.faa_mutations',\
'B.1.427_protein.faa_mutations',\
'B.1.525_protein.faa_mutations',\
'B.1.526_protein.faa_mutations',\
'B.1.617.1_protein.faa_mutations',\
'B.1.617.2_protein.faa_mutations',\
'C.37_protein.faa_mutations',\
'JN.1.11.1_protein.faa_mutations',\
'KP.3.1.1_protein.faa_mutations',\
'P.1_protein.faa_mutations',\
'P.2_protein.faa_mutations',\
'P.3_protein.faa_mutations',\
'XBB.1.5_protein.faa_mutations',\
'XEC_protein.faa_mutations'\
]

path = '' #Path to above files. Same will be used for output file

outfile = 'mutations.csv'

mutation_dict = dict()

i=0
header = 'pos,AA change'
for mut_file in protein_mut_files:
    with open(path+mut_file, 'r') as file:
        data = json.load(file)
    variant = mut_file.split("_")[0]
    header += ','+variant
    
    for mut in data.keys():
        if mut not in mutation_dict:
            mutation_dict[mut] = [0]*len(protein_mut_files)
        if data[mut] > 50:
            mutation_dict[mut][i] = 1
    i += 1
header += '\n'

with open(path+outfile, 'w') as outhandle:
    outhandle.write(header)
    for mut in mutation_dict:
        if 1 in mutation_dict[mut]:
            pos = mut[1:-1]
            line = pos+','+mut
            for item in mutation_dict[mut]:
                line += ','+str(item)
            line += '\n'
            outhandle.write(line)

        


