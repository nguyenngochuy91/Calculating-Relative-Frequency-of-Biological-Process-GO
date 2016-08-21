# Calculating-Relative-Frequency-of-Biological-Process-GO
## Installation 
User can either use github interface Download or type the following command in command line:
```bash
git clone https://github.com/nguyenngochuy91/Ancestral-Blocks-Reconstruction
```
The script was written in python3


## Usage
* Step 1: Using info from **"B.Subtilis_Operons_ProOpDB.txt"** , and **"gene_block_names_and_genes.txt**" to write out 2 files.
  1. A text file **“operons_genes.txt”**, this will serve as an input for **"relative_frequency_bioProcess.py"** to parse into a dic. A line example: 'bsub-BSU23440	BSU23420	BSU23430	BSU23440	BSU23402	BSU23401	BSU23410'
  2. A text file **"name_bsucyc_uniprot.txt"**, this will serve as an input for uniprot website to get protein file of all the gene. A line example: 'BSUB:BSU23420-MONOMER'
  3. Use the command line below: 
```bash
./get_genes_name.py -i B.Subtilis_Operons_ProOpDB.txt -g gene_block_names_and_genes.txt -o operons_genes.txt -n name_bsucyc_uniprot.txt
```
* Step 2: Getting the protein file from [Uniprot database](http://www.uniprot.org/uploadlists/) website using **"name_bsucyc_uniprot.txt"**
  1. In **"Provide your identifiers"**, upload file **"name_bsucyc_uniprot.txt"** 
  2. In **"Select options"**, choose From **"BioCyc"** to **"UniProtKB"**. Then hit **"Go"**
  3. In UniProtKB results, click on Download tab, choose **"Download all"** with Format **"Text"**
  4. Unzip the downloaded file, save it as **"uniprot.txt"**

* Step 3: Calculating the frequency of biological process of our operons. Using **"relative_frequency_bioProcess.py"** script to write out 4 files.
  1. 2 csv files **"top10_conserved.csv"**,**"bottom10_conserved.csv"** with 1st column is go term, 2nd column is its bioprocess, 3rd column is relative frequency. A line example: 'GO:0006810	P:transport	0.35'
  2. 2 text files **"top10_conserved.txt"**,**"bottom10_conserved.txt"** that has Biological Process Go term on each line. A line example: 'GO:0006810'. This will serve to construct a DART graph for all the Biological Process Go term, with highlight on those that are over represented. (will be provide in Step 4 in the future).
```bash
./relative_frequency_bioProcess.py -i operons_genes.txt -u uniprot.txt -s conservedOperonsSorted.txt 
```
![top10_conserved.csv](https://github.com/nguyenngochuy91/Calculating-Relative-Frequency-of-Biological-Process-GO/blob/master/image/top10_conserved.csv)
