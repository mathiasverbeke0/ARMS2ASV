import random
from Bio import SeqIO

def subset_multifasta(input_file, output_file_name, num_files):
    # Read the input multifasta file
    records = list(SeqIO.parse(input_file, 'fasta'))

    # Shuffle the records randomly
    random.shuffle(records)

    # Calculate the number of records per file
    records_per_file = len(records) // num_files

    # Create and write to the output files
    for i in range(num_files):
        output_file = f'{output_file_name}{i+1}.fasta'
        subset_records = records[i * records_per_file: (i+1) * records_per_file]

        # Write the subset of fasta records to the current file
        SeqIO.write(subset_records, output_file, 'fasta')

    print(f'{num_files} output files created.')

# Usage example
input_file = '/home/guest/MetabarcodeReflib/COI/MIDORI2_COI_LONGEST_MARINE_SPECIES.fas'
output_file_name = '/home/guest/MetabarcodeReflib/COI/MIDORI2_COI_LONGEST_MARINE_SPECIES_SUBSET'
num_files = 5
subset_multifasta(input_file, output_file_name, num_files)