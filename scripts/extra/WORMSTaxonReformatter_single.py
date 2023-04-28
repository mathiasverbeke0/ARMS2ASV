import requests
from Bio import SeqIO
from tqdm import tqdm
import argparse

# Create parser
parser = argparse.ArgumentParser(description='Reformat multifasta taxonomy description lines')

# Add arguments
parser.add_argument('-i', '--input', required = True, help = 'Multifasta input file')
parser.add_argument('-o', '--output', required = True, help = 'Multifasta output file')

# Parse the arguments
args = parser.parse_args()

# Variables for end message
l = 0
n = 0

# Count the number of description lines in the input file
with open(args.input, 'r') as in_file:
    tot = 0
    for line in in_file.readlines():
        if line.startswith('>'):
            tot += 1

# Open the input and output file
with open(args.input, "r") as in_file, open(args.output, "w") as out_file:

    # Loop over all description lines and corresponding sequences
    for record in tqdm(SeqIO.parse(in_file, "fasta"), total=tot):
        # Create list of all levels in the description line
        level_list = record.description.split(';')
        
        # Fetch the sequence
        sequence = record.seq

        flag = False

        # Iterate over every level (i.e. rank) in the level_list
        for level in level_list:
            
            scientific_name = level

            # Construct the link that will be used to query WORMS
            url = f'http://www.marinespecies.org/rest/AphiaRecordsByMatchNames?scientificnames[]={scientific_name}'
            
            # Send request to WORMS, wait a maximum of 100 seconds
            response = requests.get(url, timeout = 100)

            try:
                # Decode the response from WORMS
                data = response.json()[0][0]

            # If the response does not contain any info, go to the next level (i.e. rank)
            except Exception as e:
                continue

            # If the rank in the acquired data equals species, extract the lineage from the response        
            if data['rank'] == 'Species':
                
                # Set the flag to indicate that a species hit was found
                flag = True

                Kingdom = data['kingdom'] 
                Phylum = data['phylum']
                Class = data['class']
                Order = data['order']
                Family = data['family']
                Genus = data['genus']

                break
        
        if flag == True:
            # Write the new description line and corresponding sequence to output file
            out_file.write(f'>{";".join([Kingdom, Phylum, Class, Order, Family, Genus, scientific_name])}\n{sequence}\n')
        
        else: 
            # Extract strings with 2 or more words from the level_list
            # Note: Most species names consist out of 2 or more words
            two_word_string = [s for s in level_list if len(s.split()) >= 2]
            l += 1

            if len(two_word_string) == 0:
                two_word_string = ['Failed auto-species detection'] 
                n+= 1

            # Write a backup description line and corresponding sequence to output file
            out_file.write(f'>{";".join(["WORMS", "failed", "to", "find", "lineage", two_word_string[0]])}\n{sequence}\n')

# Print statistics after reformatting of the description lines has been performed
print(f'\nSequences: {tot}\nLineages not found: {l} ({round((l/tot)*100,2)}%)\nSpecies not found: {n}  ({round((n/tot)*100,2)}%)')