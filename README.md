# Traineeship eDNA Pipeline Repository
This repository holds all the necessary code, scripts, and data that I will use during my traineeship to develop an automated pipeline. The pipeline uses fastq files containing eDNA reads that were collected from several autonomous reef monitoring structures (ARMS). It is designed to generate ASVs, perform taxonomic classification using BOLDigger and locally stored databases, and conduct various ecological analyses.

## Installation and setup

Before you can use this pipeline, you need to install the necessary software and packages. Here are the steps to follow:

1. **Install R and Python**: You will need to have both R (version 4.0 or higher) and Python (version 3.6 or higher) installed on your machine. You can download them from the following links:
   - [R](https://www.r-project.org/)
   - [Python](https://www.python.org/downloads/)

2. **Install BOLDigger**: The pipeline uses the BOLDigger Python package for taxonomic classification. You can install it using pip with the following command:
```bash
pip install --user boldigger_cline cutadapt
```

3. **Clone the repository**: Clone this repository to your local machine using the following command:
```bash
git clone git@github.com:mathiasverbeke0/ARMS2ASV.git
```

4. **Download the reference databases**: Download the desired reference databases and put them in the `./data/databases` directory.

Once you've completed these steps, you're ready to use the pipeline!


## Usage
The DADA2.R script is used to run the pipeline for generating ASVs from eDNA data collected from ARMS using the DADA2 algorithm and performing taxonomic classification using BOLDigger and locally stored databases. The following are the available arguments for the script:

```bash
usage: DADA2.R [-h] -b dir [-d FILENAME] -r {single,multi} [-u USER]
               [-P PASSWORD] [-t] [-p Fwd_Primer Rev_Primer] [-s] [-l Fwd Rev]
               [-m minlen] [-B]
```

For more information on the available arguments and their usage, you can run the command:

```bash
DADA2.R -h
```