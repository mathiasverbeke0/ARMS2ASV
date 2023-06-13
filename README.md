# ReefPipe

ReefPipe is a command line workflow tool developed to address the challenging task of analyzing large amounts of metabarcoding data. It specifically focuses on processing COI and ITS metabarcoding data, automating parallel Amplicon Sequence Variant (ASV) inference and multi-reference taxonomic classification.

## Installation and Setup

To use ReefPipe, follow these setup steps. If you encounter installation issues, consider using the Docker image. It includes all required software components and can be accessed at [ReefPipe-Docker](https://github.com/mathiasverbeke0/ReefPipe-Docker). Using the Docker image eliminates individual installations and potential compatibility problems, as all dependencies are preconfigured.

### R Software

Make sure you have R version 4.1.3 or higher installed on your system. If you don't have it already, you can download and install it from the [official R project website](https://www.r-project.org/).

### Python Software

Make sure you have Python version 3.6 or higher installed on your system. If you don't have it already, you can download and install it from the [official Python website](https://www.python.org/).

### R Packages

R packages are automatically installed by ReefPipe scripts, eliminating the need for manual installation. However, errors might occur during the download of the DADA2<sup>1</sup> package, particularly on certain systems. In such cases, you may need to install additional software specific to your operating system. Please note that we cannot provide troubleshooting support for all operating systems, so it is recommended to address these issues independently.

### Python Modules

After installing Python, open your terminal and run the following command to install the necessary modules:

```bash
pip install boldigger_cline cutadapt argparse biopython tqdm pandas
```

### Clustal Omega Software

For Linux distributions that use the dnf package manager, you can download and install Clustal Omega by executing the following command in your terminal:

```bash
sudo dnf install clustal-omega
```

For Linux distributions that use the apt package manager, you can download and install Clustal Omega by executing the following command in your terminal:

```bash
sudo apt-get install clustalo
```

For Windows and macOS users, a Clustal Omega binary has been incorporated into the ReefPipe script<sup>1</sup>.

### ReefPipe

If you have Git initialized on your system, you can clone or download the ReefPipe code from the GitHub repository by running the following command in your terminal:

```
git clone git@github.com:mathiasverbeke0/ReefPipe.git
```

If you don't have Git initialized on your system, you can download the ReefPipe code as a ZIP file by following these steps:

1. Visit the ReefPipe GitHub repository at https://github.com/mathiasverbeke0/ReefPipe.
2. Click on the green "Code" button.
3. Select "Download ZIP" from the dropdown menu.
4. Once the ZIP file is downloaded, extract its contents to a location of your choice on your system.

## Documentation
For detailed instructions on how to use ReefPipe, refer to the [documentation](./user_guide/).

## References
Please find below the list of references that are relevant to ReefPipe:

* <sup>1</sup> Callahan, B. (s.a.). DADA2 Pipeline Tutorial (1.16). https://benjjneb.github.io/dada2/tutorial.html
* Buchner D, Leese F (2020) BOLDigger – a Python package to identify and organise sequences with the Barcode of Life Data systems. Metabarcoding and Metagenomics 4: e53535. https://doi.org/10.3897/mbmg.4.53535
* Cock PA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B and de Hoon MJL (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25, 1422-1423. https://doi.org/10.1093/bioinformatics/btp163
* Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet. journal, 17(1), 10-12. https://doi.org/10.14806/ej.17.1.200 
* <sup>2</sup> Sievers, F., Wilm, A., Dineen, D., Gibson, T. J., Karplus, K., Li, W., . . . Söding, J. (2011). Fast, scalable generation of high‐quality protein multiple sequence alignments using Clustal Omega. Molecular Systems Biology, 7(1), 539. https://doi.org/10.1038/msb.2011.75 

Please refer to these references for further reading and to understand the underlying concepts and methodologies used in ReefPipe.