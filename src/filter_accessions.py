import json

# Path to the file containing JSON objects on each line
path = '' #Path to below files
data_record_file = path+'data_report.jsonl'
annotation_file = path+'annotation_report.jsonl'
data_record_outfile = path+'variant_accession.json'
n_p_map = path+'nucleotide_protein_map.json'

variants_interest = ['B.1.1.7', 'B.1.351', 'B.1.525', 'B.1.427', 'B.1.526', 'B.1.617.1', 'B.1.617.2', 'C.37', 'P.1', 'P.2', 'P.3', 'B.1.1.529', 'XBB.1.5', 'JN.1.11.1', 'KP.3.1.1', 'XEC']

variants = dict()
n_accessions = set()
# Open the file and parse each line as a separate JSON object
with open(data_record_file, 'r') as file:
    for line in file:
        line_1 = line.strip()
        try:
            data = json.loads(line_1)
            if 'completeness' in data and data['completeness'] == 'COMPLETE' and 'host' in data and 'organismName' in data['host'] and data['host']['organismName'] == 'Homo sapiens' and 'virus' in data and 'pangolinClassification' in data['virus']:
                variant_name = data['virus']['pangolinClassification']
                if variant_name in variants_interest:
                    if variant_name not in variants:
                        variants[variant_name] = []
                    variants[variant_name].append(data['accession'])
                    n_accessions.add(data['accession'])
        except json.JSONDecodeError as e:
            print(f"Error decoding data_record JSON on line: {line}. Error: {e}")

annotation_records = dict()

# Open the file and parse each line as a separate JSON object
with open(annotation_file, 'r') as file:
    for line in file:
        line = line.strip()
        try:
            json_data = json.loads(line)
            if json_data['accession'] in n_accessions and 'genes' in json_data:
                for gene in json_data['genes']:
                    if gene['name'] == 'S' and 'cds' in gene:
                        for cd in gene['cds']:
                            if 'name' in cd and cd['name'] == 'surface glycoprotein':
                                annotation_records[cd['nucleotide']['accessionVersion']] = { \
                                    'nucl_range_begin' : cd['nucleotide']['range'][0]['begin'] , \
                                    'nucl_range_end' : cd['nucleotide']['range'][0]['end'] , \
                                    'prot_acc' : cd['protein']['accessionVersion'] , \
                                    'prot_range_begin' : cd['protein']['range'][0]['begin'] , \
                                    'prot_range_end' : cd['protein']['range'][0]['end'] \
                                    }
        except json.JSONDecodeError as e:
            print(f"Error decoding annotation JSON on line: {line}. Error: {e}")

with open(data_record_outfile, 'w') as file:
    json.dump(variants, file, indent=4)

with open(n_p_map, 'w') as file:
    json.dump(annotation_records, file, indent=4)
