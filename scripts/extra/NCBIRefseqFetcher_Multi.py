from Bio import Entrez, SeqIO
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed

############################################################################################
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


# Set search parameters
organisms = ['Ascomycota', 'Ochrophyta', 'Rhodophyta', 'Chlorophyta', 'Charophyta', 'Tracheophyta', 'Porifera', 'Cnidaria', 'Ctenophora', 'Platyhelminthes', 'Nemertea', 'Sipuncula', 'Mollusca', 'Annelida', 'Arthropoda', 'Phoronida', 'Bryozoa', 'Brachiopoda', 'Echinodermata', 'Chordata'] # replace with your organism of interest
gene_name = 'COI' # replace with your gene of interest

# Tell NCBI who you are
Entrez.email = input('Provide NCBI email: ')
print('\n')

# Loop over the organisms
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
            print('None found\n')
            continue
        
        with ThreadPoolExecutor(max_workers=3) as executor:
        # Submit tasks to the executor and store the future objects in a list
            futures = [executor.submit(getInfo, id) for id in id_list]

            # Wait for all the tasks to complete and print the results
            for future in tqdm(as_completed(futures), total = id_tot):
                tax, seq = future.result()
                f.write(f'>{tax}\n{seq}\n')
        
        print('\n')