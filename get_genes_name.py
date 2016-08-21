#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : For each B.Sub operon, get its gene (in bsu-BSU format)
    Start   : 08/18/2016
    End     : 08/20/2016
'''
import argparse

# get the arguments from command line
def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--Operon","-i",help="ProOperon dp file")
    parser.add_argument("--GeneBlock","-g", help="gene blocks that went through filtered")
    parser.add_argument("--Outfile","-o", help="File out for the new result dic")
    parser.add_argument("--GeneName","-n", help="Gene name to parse into uniprot to get protein file")
    args = parser.parse_args()
    return args
###############################################################################
# Quick parsers to get info from files
###############################################################################
# parse the gene_block_names_and_genes.txt file
def parse_gene_block_names_and_genes(myfile):
    infile = open(myfile,'r')
    mydic={}
    for line in infile.readlines():
        line = line.strip()
        items = line.split('\t')
        # print items
        mydic[items[0]] = items[1:]
    return mydic
# parse the B.Subtilis operon txt file
def parse_operon(myfile):
    infile = open(myfile,'r')
    mydic={}
    # ignore the first line
    for line in infile.readlines()[0:]:
        line = line.split('\t')
        mydic[line[1]] = line[2]
    return mydic


# using the 2 dictionary to write out new operon file, where each line starts with 
# starter gene name, follow by what gene in the operon in bsu-BSU format
# in addition, also wrote out gene names in following format BSUB:BSU32510-MONOMER
def return_dic_bsu(dic1,dic2,outfile,GeneName):
    out= open(outfile,'w')
    name = open(GeneName,'w')
    for item in dic1:
        # initiate a list
        out.write(item)
        for gene in dic1[item]: # iterate through gene annotation
            name.write('BSUB:'+dic2[gene]+'-MONOMER'+'\n')
            out.write('\t'+dic2[gene])
        out.write('\n')
    out.close()
    
    
if __name__ == "__main__":
    args       = get_arguments()
    Operon     = args.Operon
    GeneBlock  = args.GeneBlock
    GeneName   = args.GeneName
    Outfile    = args.Outfile
    operon_dic = parse_operon(Operon)
    geneBlock_dic = parse_gene_block_names_and_genes(GeneBlock)
    return_dic_bsu(geneBlock_dic,operon_dic,Outfile,GeneName)
    
