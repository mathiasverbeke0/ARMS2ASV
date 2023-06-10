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
parser.add_argument('-m', '--multifasta', required = True, help = 'The multifasta file used for MSA.')
parser.add_argument('-o', '--output_dir', required = True, help = 'The output directory path.')
parser.add_argument('-e', '--environment', required = True, choices = ['linux', 'windows'], help = 'The operating system you are currently using.')
parser.add_argument('-c', '--clustal_omega', required = False, help = 'The absolute path to the clustalo executable')

# Parse the arguments
args = parser.parse_args()


#########################################
## PERFORM MULTIPLE SEQUENCE ALIGNMENT ##
#########################################

# Define the input file name
input_file = args.multifasta

# Define the output file name
output_file = f'{args.output_dir.rstrip("/")}/Alignment.fasta'

if args.environment == 'linux' :
    # Run Clustal Omega multiple sequence alignment command
    clustalomega_cline = ClustalOmegaCommandline(infile=input_file, outfile=output_file, verbose=True, auto=True, force = True)
    stdout, stderr = clustalomega_cline()

elif args.environment == 'windows':
    clustalomega_cline = ClustalOmegaCommandline(args.clustal_omega, infile=input_file, outfile=output_file, verbose=True, auto=True, force = True)
    stdout, stderr = clustalomega_cline()

print(f"Aligned sequences written to {output_file}")
