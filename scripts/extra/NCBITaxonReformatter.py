from Bio import SeqIO
from ete3 import NCBITaxa
import argparse
from tqdm import tqdm

# Create parser
parser = argparse.ArgumentParser(description='Reformat multifasta taxonomy description lines')

# Add arguments
parser.add_argument('-i', '--input', required = True, help = 'Multifasta input file')
parser.add_argument('-o', '--output', required = True, help = 'Multifasta output file')

# Parse the arguments
args = parser.parse_args()

# Load the NCBI taxonomy database
ncbi = NCBITaxa()

# Variables for end message
s = 0
l = 0
n = 0

with open(args.input, 'r') as in_file:
    tot = 0
    for line in in_file.readlines():
        if line.startswith('>'):
            tot += 1

# Open the input and output files
with open(args.input, "r") as in_file, open(args.output, "w") as out_file:

    # Parse the fasta sequences
    for record in tqdm(SeqIO.parse(in_file, "fasta"), total=tot):
        # Create list of all levels in the description line
        level_list = record.description.split(';')
        sequence = record.seq
        s += 1

        # Remove empty strings from the list
        level_list = [level for level in level_list if level != '']

        # Get taxonomic id for every level (e.g. kingdom, phylum, etc.) in the description line
        level_taxids = [ncbi.get_name_translator([level]) for level in level_list]

        # Merge the list of dictionaries into one dictionary
        level_taxids = {k: v[0] for d in level_taxids for k, v in d.items()}

        # Determine the lineage of all levels in level_taxids using its taxid and store in a dictionary
        # Note: The lineage is a list of taxids with the most right taxid being the one of the level
        taxid_lineage = {taxid: ncbi.get_lineage(taxid) for taxid in level_taxids.values()}

        # Determine the ranks (i.e. levels) of all taxids in the lineages
        taxid_ranks = {taxid: ncbi.get_rank(lineage) for taxid, lineage in taxid_lineage.items()}

        # Determine the names (i.e. names belonging to the levels) of all taxids in the lineages
        taxid_names = {taxid: ncbi.get_taxid_translator(lineage) for taxid, lineage in taxid_lineage.items()}

        # Check if taxid is taxid of species level
        flag = False

        for taxid, rank in taxid_ranks.items():

            # If taxid is of species level, pick the rank and name list for which the taxid is the key
            if rank[int(taxid)] == 'species':
                flag = True

                species_taxid_rank = rank
                species_taxid_name = taxid_names[taxid]

                # Extract the kingdom, phylum, class, order, family, genus and species
                for level_taxid, level in species_taxid_rank.items():
                    if level == 'kingdom':
                        Kingdom = species_taxid_name[level_taxid]
                    elif level == 'phylum':
                        Phylum = species_taxid_name[level_taxid]
                    elif level == 'class':
                        Class = species_taxid_name[level_taxid]
                    elif level == 'order':
                        Order = species_taxid_name[level_taxid]
                    elif level == 'family':
                        Family = species_taxid_name[level_taxid]
                    elif level == 'genus':
                        Genus = species_taxid_name[level_taxid]
                    elif level == 'species':
                        Species = species_taxid_name[level_taxid]

                break

        if flag == True:
            # Write the new description line and corresponding sequence to output file
            out_file.write(f'>{";".join([Kingdom, Phylum, Class, Order, Family, Genus, Species])}\n{sequence}\n')

        else: 
            # Extract strings with 2 or more words from the level_list
            # Note: species names consist out of 2 or more words
            two_word_string = [s for s in level_list if len(s.split()) >= 2]
            l += 1

            if len(two_word_string) == 0:
                two_word_string = ['Failed auto-species detection'] 
                n+= 1

            # Write a backup description line and corresponding sequence to output file
            out_file.write(f'>{";".join(["NCBI", "failed", "to", "find", "lineage", two_word_string[0]])};\n{sequence}\n')

print(f'\nSequences: {s}\nLineages not found: {l} ({round((l/s)*100,2)}%)\nSpecies not found: {n}  ({round((n/s)*100,2)}%)')