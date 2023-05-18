####################
## IMPORT MODULES ##
####################

import os, argparse, glob
from tqdm import tqdm
from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline


##################################
## Parse command line arguments ##
##################################

# Create parser
parser = argparse.ArgumentParser(description='Merging and subsequent MSA of multiple multifasta files')

# Add arguments
parser.add_argument('-b', '--base_dir', required = True, help = 'The base directory path for the analysis.')
parser.add_argument('-o', '--output_dir', required = True, help = 'The output directory path.')

# Parse the arguments
args = parser.parse_args()


##################################
## COMBINE ALL MULTIFASTA FILES ##
##################################

# Define the search directory
search_directory = args.base_dir

# Specify the file pattern to search for
file_pattern = os.path.join(search_directory, '**/COI_ASVS.fasta')

# Find all matching files in the specified directory and subdirectories
matching_files = glob.glob(file_pattern, recursive=True)

if len(matching_files) > 0:
    # Initialize an empty list to store the sequences
    sequences = []

    print('Reading sequences')
    
    for file_path in tqdm(matching_files, total=len(matching_files),   ):

        # Read sequences from each COI_ASVS.fasta file
        with open(file_path, "r") as fasta_file:
            
            records = SeqIO.parse(fasta_file, "fasta")
            
            for record in records:
            
                sequences.append(record)
    
    # Write the combined sequences to a new multifasta file
    output_file = f'{args.output_dir.rstrip("/")}/CombinedCOIASVS.fasta'
    
    with open(output_file, "w") as outfile:
        
        SeqIO.write(sequences, outfile, "fasta")

    print(f"Combined sequences written to {output_file}.")

    #########################################
    ## PERFORM MULTIPLE SEQUENCE ALIGNMENT ##
    #########################################

    # Define the input file name
    input_file = output_file

    # Define the output file name
    output_file = f'{args.output_dir.rstrip("/")}/AlignedCOIASVS.fasta'

    print('\nPerforming MSA')

    # Run Clustal Omega multiple sequence alignment command
    clustalomega_cline = ClustalOmegaCommandline(infile=input_file, outfile=output_file, verbose=True, auto=True, force = True)
    stdout, stderr = clustalomega_cline()

    print(f"Aligned sequences written to {output_file}.")

else:
    exit('Error: No multifasta files called COI_ASVS.fasta were found.')
