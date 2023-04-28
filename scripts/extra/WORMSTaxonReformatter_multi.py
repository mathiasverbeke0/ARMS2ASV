import requests
import json
from Bio import SeqIO
from tqdm import tqdm
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed

# Create parser
parser = argparse.ArgumentParser(description='Reformat multifasta taxonomy description lines')

# Add arguments
parser.add_argument('-i', '--input', required = True, help = 'Multifasta input file')
parser.add_argument('-o', '--output', required = True, help = 'Multifasta output file')

# Parse the arguments
args = parser.parse_args()

# Variables for end message
s = 0
l = 0
n = 0

# WORMFetch Function
def WORMFetch(record):
        # Create list of all levels in the description line
        level_list = record.description.split(';')
        sequence = record.seq
        
        global s 
        s += 1

        flag = False

        for level in level_list:
            
            scientific_name = level

            url = f'http://www.marinespecies.org/rest/AphiaRecordsByMatchNames?scientificnames[]={scientific_name}'

            response = requests.get(url)

            try:
                data = response.json()[0][0]

            except Exception as e:
                continue

            if data['rank'] == 'Species':
                
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
            return f'>{";".join([Kingdom, Phylum, Class, Order, Family, Genus, scientific_name])}\n{sequence}\n'
        
        else: 
            # Extract strings with 2 or more words from the level_list
            # Note: species names consist out of 2 or more words
            two_word_string = [s for s in level_list if len(s.split()) >= 2]
            global l 
            l += 1

            if len(two_word_string) == 0:
                two_word_string = ['Failed auto-species detection'] 
                global n
                n += 1

            # Write a backup description line and corresponding sequence to output file
            return f'>{";".join(["WORMS", "failed", "to", "find", "lineage", two_word_string[0]])}\n{sequence}\n'

with open(args.input, 'r') as in_file:
    tot = 0
    for line in in_file.readlines():
        if line.startswith('>'):
            tot += 1

with open(args.input, "r") as in_file, open(args.output, "w") as out_file:
    with ThreadPoolExecutor(max_workers=3) as executor:
        # Submit tasks to the executor and store the future objects in a list
            futures = [executor.submit(WORMFetch, record) for record in SeqIO.parse(in_file, "fasta")]

            # Wait for all the tasks to complete and print the results
            for future in tqdm(as_completed(futures), total = tot):
                result = future.result()
                out_file.write(result)

print(f'\nSequences: {s}\nLineages not found: {l} ({round((l/s)*100,2)}%)\nSpecies not found: {n}  ({round((n/s)*100,2)}%)')