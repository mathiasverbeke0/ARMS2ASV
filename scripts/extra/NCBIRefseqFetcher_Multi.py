####################
## Import modules ##
####################
from Bio import Entrez, SeqIO
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse, sys


##################################
## Parse command line arguments ##
##################################

# Create parser
parser = argparse.ArgumentParser(description='Fetch NCBI reference sequences.')

# Add arguments
parser.add_argument ('-g', '--GOI', required = True, help = 'Specify gene of interest.')
parser.add_argument('-o', '--OUT', required = True, help = 'Specify multifasta output file')
parser.add_argument('-t', '--TAX', nargs='+', help='Specify one or more organisms using their genus and species, or taxonomic level (e.g. Homo). Example: Escherichia_coli Saccharomyces_cerevisiae Homo Mollusca.')
parser.add_argument('-e', '--email', required = True, help = 'Email you use to access NCBI.')

# Parse the arguments
args = parser.parse_args()

def getInfo(id):
    # Search for Refseq gene sequence
    handle = Entrez.efetch(db="nucleotide", id=id, rettype="fasta", retmode="text")
    
    # Extract RefSeq gene sequence
    seq_record = SeqIO.read(handle, "fasta")
    seq = seq_record.seq

    # Search for taxonomy
    handle2 = Entrez.efetch(db="nucleotide", id=id, rettype="gbwithparts", retmode="text")
    
    # Extract taxonomy
    gb_record = SeqIO.read(handle2, "genbank")
    tax_dict = gb_record.annotations["taxonomy"]

    # Prettify taxonomy
    tax_str = ";".join(tax_dict)

    # Write taxonomy and sequence to multifasta file (can be used by DADA2 assignTaxonomy)
    return tax_str, seq

organisms = [' '.join(organism.split('_')).capitalize() for organism in args.TAX]

gene_name = args.GOI.upper()

# Tell NCBI who you are
Entrez.email = args.email

# Loop over the organisms
with open(args.OUT, 'w') as f:
    for organism in organisms:
        print(f'Fetching {gene_name} sequences for {organism}')
        # Search for RefSeq gene ids
        handle = Entrez.esearch(db='nucleotide', term=f'"{organism}"[Organism] AND {gene_name}[gene] AND refseq[filter]', retmax=100000)
        
        # Extract RefSeq gene ids
        record = Entrez.read(handle)
        id_list = record['IdList']
        id_tot = len(id_list)

        if id_tot == 0:
            print('None found')
            continue
        
        with ThreadPoolExecutor(max_workers=3) as executor:
        # Submit tasks to the executor and store the future objects in a list
            futures = [executor.submit(getInfo, id) for id in id_list]

            # Wait for all the tasks to complete and print the results
            for future in tqdm(as_completed(futures), total = id_tot):
                try:
                    tax, seq = future.result()
                    f.write(f'>{tax}\n{seq}\n')
                
                except Exception as e:
                        print('An error occurred: {e}\Halting all subsequent executions. This might take some time.')
                        
                        # Cancel all remaining futures
                        for remaining_future in futures:
                            if not remaining_future.done():
                                remaining_future.cancel()

                        sys.exit()