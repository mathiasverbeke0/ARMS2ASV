from Bio import Entrez, SeqIO
from tqdm import tqdm

# Set search parameters
organisms = ['Annelida', 'Arthropoda', 'Chaetognatha', 'Chordata', 'Cnidaria and Ctenophora', 
                    'Echinodermata', 'Mollusca'] # replace with your organism of interest
gene_name = 'COI' # replace with your gene of interest

# Tell NCBI who you are
Entrez.email = 'mathias.v2000@gmail.com'

with open('testfile.fasta', 'w') as f:
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

        for id in tqdm(id_list, total=id_tot):
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
            f.write(f'>{tax_str}\n{seq}\n')