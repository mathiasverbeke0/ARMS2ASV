from Bio import SeqIO
import re, sys
from tqdm import tqdm

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

# Specify the filename
filename = '/home/guest/Traineeship/ARMSProject/data/species.txt'  # Replace with your file name

# Read lines with multiple words from the file
organisms = read_lines_with_multiple_words(filename)

with open('/home/guest/Traineeship/ARMSProject/data/MIDORI2_UNIQ_NUC_GB255_CO1_DADA2.fasta', 'r') as infile:
    records = SeqIO.parse(infile, "fasta")

    print('Counting sequences...')

    count = 0
    for record in records:
        count += 1


with open('/home/guest/Traineeship/ARMSProject/data/MIDORI2_LONGEST_NUC_GB255_CO1_DADA2.fasta', 'r') as infile, open('/home/guest/Traineeship/ARMSProject/data/MIDORI2_COI_LONGEST.fas', 'w') as outfile:
    records = SeqIO.parse(infile, "fasta")

    print('Extracting sequences...')

    count2 = 0
    for record in tqdm(records, total = count):
        for organism in organisms:
            if organism in record.description:
                count2 +=1

                description = re.sub(r'\d+', '', record.description)
                description = description.replace('_', '')
                outfile.write(f'>{description}\n{record.seq}\n')

print(f'Extracted {count2} sequences')