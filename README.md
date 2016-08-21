# Calculating-Relative-Frequency-of-Biological-Process-GO
## Installation 
User can either use github interface Download or type the following command in command line:
```bash
git clone https://github.com/nguyenngochuy91/Ancestral-Blocks-Reconstruction
```
The script was written in python3


## Usage
* Step 1: Using info from **"B.Subtilis_Operons_ProOpDB.txt"** , and **"gene_block_names_and_genes.txt**" to write out 2 files.
  1. A text file **“operons_genes.txt”**, this will serve as an input for relative_frequency_bioProcess.py to parse into a dic. An example line is as follow: 'bsub-BSU23440	BSU23420	BSU23430	BSU23440	BSU23402	BSU23401	BSU23410'
  2. A text file **"name_bsucyc_uniprot.txt"**, this will serve as an input for uniprot website to get protein file of all the gene 
  3. Use the command line below:
```bash
./get_genes_name.py -i B.Subtilis_Operons_ProOpDB.txt -g gene_block_names_and_genes.txt -o operons_genes.txt -n name_bsucyc_uniprot.txt
```
* Step 2: 
