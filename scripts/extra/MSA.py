from Bio.Align.Applications import ClustalOmegaCommandline
import tqdm as tqdm

# Define the input file name
input_file = "/home/guest/Traineeship/ARMSProject/data/temp2/ARMS_SWC_Gbg3_20200206_20200529/06.Seq_Table/COI_ASVS_done.fasta"

# Define the output file name
output_file = "aligned_sequences.fasta"

# Run Clustal Omega multiple sequence alignment command
clustalomega_cline = ClustalOmegaCommandline(infile=input_file, outfile=output_file, verbose=True, auto=True, force = True)
stdout, stderr = clustalomega_cline()