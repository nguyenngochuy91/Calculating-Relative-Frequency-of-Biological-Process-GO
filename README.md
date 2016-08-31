# Calculating-Relative-Frequency-of-Biological-Process-GO
## Installation
1. Download Biopython, recommend using annaconda.
2. Clone or download this  [Ontology Parser](https://github.com/kkoziara/biopython/tree/master/Bio/Ontology) into where Biopython is install (usually in your python/site-packages/Bio)
3. Download go-basic.obo file here [gene ontology file](http://purl.obolibrary.org/obo/go/go-basic.obo)
4. User can either use github interface download or type the following command in command line:
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

* Step 3: Calculating the frequency of biological process of our operons. Using **"relative_frequency_bioProcess.py"** script to write out 5 files.
  1. 2 csv files **"top10_conserved_level2.csv"**,**"bottom10_conserved_level2.csv"** with 1st column is go term, 2nd column is its bioprocess, 3rd column is relative frequency. A line example: 'GO:0006810	P:transport	0.35'
  2. 2 text files **"top10_conserved_level2.txt"**,**"bottom10_conserved_level2.txt"** that has Biological Process Go term on each line. 
  3. 1 csv file **"comparison2.csv"** with 1st is representative_level, 2nd is top10_only, 3rd is bottom10_only, 4th is both.
  4. Following is an example of running the with go term at 2nd level from the root of the Gene Ontology
```bash
./relative_frequency_bioProcess.py -i operons_genes.txt -u uniprot.txt -s conservedOperonsSorted.txt -l 2 -g go-basic.obo 
```

