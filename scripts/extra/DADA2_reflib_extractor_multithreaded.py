from Bio import SeqIO
import re, sys
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed

# Function to check if a line contains more than one word
def has_multiple_words(line):
    words = line.split()
    return len(words) > 1

# Read lines from a file and store lines with more than one word in a list
def read_lines_with_multiple_words(filename):
    lines_with_multiple_words = []
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if has_multiple_words(line):
                lines_with_multiple_words.append(line)
    return lines_with_multiple_words

# Checking the description line for desired organisms
def check_description_line(record, organisms):
    for organism in organisms:

        description = record.description
        description = description.replace('_', ' ') ## THIS IS ONLY FOR THE 18S DATABASE

        if organism in description:
            #description = re.sub(r'_\d+', '', record.description) ## THIS IS ONLY FOR THE MIDORI COI DATABASES
            #description = description.replace('order_', '').replace('class_', '') ## THIS IS ONLY FOR THE MIDORI COI DATABASES
               
            
            return description, record.seq, organism
            
    return None, None, None

# Specify the filename
filename = '/home/guest/Traineeship/ARMSProject/data/marineSpecies.txt'  # Replace with your file name

# Read lines with multiple words from the file
organisms = read_lines_with_multiple_words(filename)

with open('/home/guest/Traineeship/ARMSProject/data/pr2_version_5.0.0_SSU_dada2.fasta', 'r') as infile:
    records = SeqIO.parse(infile, "fasta")

    print('Counting sequences...')

    count = 0
    for record in records:
        count += 1


with open('/home/guest/Traineeship/ARMSProject/data/pr2_version_5.0.0_SSU_dada2.fasta', 'r') as infile, open('/home/guest/Traineeship/ARMSProject/data/PR2_18S_rRNA.fas', 'w') as outfile:
    records = SeqIO.parse(infile, "fasta")

    print('Preparing sequence extraction...')

    count2 = 0
    
    with ThreadPoolExecutor(max_workers=8) as executor:
    # Submit tasks to the executor and store the future objects in a list
        futures = [executor.submit(check_description_line, record, organisms) for record in tqdm(records, total = count)]
        
        print('Extracting sequences...')
        
        # Wait for all the tasks to complete and print the results
        for future in tqdm(as_completed(futures), total = count):
            try:
                description, sequence, organism = future.result()

                if description == None:
                        continue
                
                else:
                        count2 +=1
                        outfile.write(f'>{description}\n{sequence}\n')
            
            except Exception as e:
                    print(f'An error occurred: {e}\Halting all subsequent executions. This might take some time.')
                    
                    # Cancel all remaining futures
                    for remaining_future in futures:
                        if not remaining_future.done():
                            remaining_future.cancel()

                    sys.exit('')

print(f'Extracted {count2} out of {count} sequences ({round((count2/count)*100, 2)}%)')