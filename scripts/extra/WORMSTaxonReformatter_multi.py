import requests
from Bio import SeqIO
from tqdm import tqdm
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import sys

# Create parser
parser = argparse.ArgumentParser(description='Reformat multifasta taxonomy description lines')

# Add arguments
parser.add_argument('-i', '--input', required = True, help = 'Multifasta input file')
parser.add_argument('-o', '--output', required = True, help = 'Multifasta output file')
parser.add_argument('-d', '--delimiter', required = True, choices = ['colon', 'semicolon', 'pipe', 'comma'], help = 'Description line delimiter of input file')

# Parse the arguments
args = parser.parse_args()

# Description line delimiter

if args.delimiter == 'colon':
     delimiter = ':'
elif args.delimiter == 'semicolon':
     delimiter = ';'
elif args.delimiter == 'pipe':
     delimiter = '|'
elif args.delimiter == 'comma':
    delimiter = ','

# Variables for end message
display = True
l = 0
n = 0

# WORMFetch Function
def WORMFetch(record):
        # Create list of all levels in a single description line
        level_list = record.description.split(delimiter)
        
        # Fetch the sequence
        sequence = record.seq

        # Remove empty strings from the list
        level_list = [level for level in level_list if level != '']

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

                Kingdom = str(data['kingdom'])
                Phylum = str(data['phylum'])
                Class = str(data['class'])
                Order = str(data['order'])
                Family = str(data['family'])
                Genus = str(data['genus'])

                break
        
        if flag == True:
            # Write the new description line and corresponding sequence to output file
            return f'>{";".join([Kingdom, Phylum, Class, Order, Family, Genus, scientific_name])};\n{sequence.upper()}\n'
        
        else: 
            # Extract strings with 2 or more words from the level_list
            # Note: Most species names consist out of 2 or more words
            two_word_string = [s for s in level_list if len(s.split()) >= 2]
            global l 
            l += 1

            if len(two_word_string) == 0:
                two_word_string = ['Failed auto-species detection'] 
                global n
                n += 1

            # Write a backup description line and corresponding sequence to output file
            return f'>{";".join(["WORMS", "failed", "to", "find", "lineage", two_word_string[0]])};\n{sequence.upper()}\n'

# Count the number of description lines in the input file
with open(args.input, 'r') as in_file:
    tot = 0
    for line in in_file.readlines():
        if line.startswith('>'):
            tot += 1

# Open the input and output file
with open(args.input, "r") as in_file, open(args.output, "w") as out_file:
    
    # Use multithreading to send multiple requests to WORMS at the same time 
    with ThreadPoolExecutor(max_workers=8) as executor:
        # Submit tasks to the executor and store the future objects in a list
            futures = [executor.submit(WORMFetch, record) for record in SeqIO.parse(in_file, "fasta")]

            # Wait for all the tasks to complete and print the results
            for future in tqdm(as_completed(futures), total = tot):
                try:
                    result = future.result()
                    out_file.write(result)
                
                except Exception as e:
                    
                    print(f"An error occurred: {e}\nHalting all subsequent executions. This might take some time.")

                    # Cancel all remaining futures
                    for remaining_future in futures:
                        if not remaining_future.done():
                            remaining_future.cancel()
                    
                    sys.exit()

# Print statistics after reformatting of the description lines has been performed
print(f'\nSequences: {tot}\nLineages not found: {l} ({round((l/tot)*100,2)}%)\nSpecies not found: {n}  ({round((n/tot)*100,2)}%)')