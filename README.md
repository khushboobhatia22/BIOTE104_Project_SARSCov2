Project: 'Evolutionary Dynamics of the SARS-CoV-2 Spike Protein: Sequence Analysis, Phylogenetic Relationships, and Vaccine-Induced Selective Pressure'

This repository contains python source files and data files for the project. Below is a brief description of important files and directories:

src/ #Source code used for data preparation and sequence analysis

src/filter_accessions.py #Processes metadata from NCBI dataset and creates a filetered list of accession numbers for each variant.

src/filter_protein_fasta.py #Reads a protein FASTA file and creates filtered FASTA files by variant

src/pairwise_global_alignment_protein.py #Reads a protein FASTA file and performs Pairwsise Alignment of all sequences with reference sequence.

src/multiple_alignment_consensus_nucl.py #Reads a FASTA file and performs multiple alignment of all the sequences in the file. Also computes a consensus sequence.

src/entropy_nucl.py #Reads a FASTA file containing multiple aligned sequences and calculates entropy for each position.


filtered_dataset/ #Filtered protein and coding sequences and json files containing mapping of protein and genome accession numbers per variant.


aligned_sequences/ #FASTA files for each variant containing aligned sequences of all the samples used for analysis.


consensus_sequences/ #FASTA file containing consensus coding sequence for each variant


mutations/ #Mutations observed in each variant



Link to the paper will be updated soon.
