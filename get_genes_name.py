#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : For each E Coli operon, get its gene (in bsu-BSU format)
    Start   : 08/18/2016
    End     : 08/20/2016
'''
import argparse
from Bio import SeqIO
synonym_dic ={'rapZ':'yhbJ','cas1':'ygbT','cas2':'ygbF','yjeE':'tsaE',
              'gutM':'srlM','gutQ':'srlQ','yjeF':'nnr'}
# get the arguments from command line
def get_arguments():
    parser = argparse.ArgumentParser()
    # parser.add_argument("--Operon","-i",help="ProOperon dp file")
    parser.add_argument("--GeneBlock","-g", help="gene blocks that went through filtered 'gene_block_names_and_genes.txt'")
    parser.add_argument("--Outfile","-o", help="File out for the new result dic 'operons_genes.txt'")
    parser.add_argument("--GeneName","-n", help="Gene name to parse into uniprot to get protein file 'name_ecocyc_uniprot.txt'")
    parser.add_argument("--Reference","-r", help="Reference genome 'NC_000913.gbk'")
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

# parse the NC_000913.gbk file
def parse_reference(myfile):
    mydic={}
    record = SeqIO.read(myfile,"genbank")
    for feature in record.features:
        if feature.type == 'CDS':
            for item in feature.qualifiers['db_xref']:
                if item[:2]=='GI':
                    # print feature
                    mydic[feature.qualifiers['gene'][0]] = item[3:]

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
            if gene in synonym_dic:
                gene = synonym_dic[gene]
            try:
                name.write(dic2[gene]+'\n')
                out.write('\t'+dic2[gene])
            except:      
                print gene
                out.write('\t'+dic2[gene])
        out.write('\n')
    out.close()
    
    
if __name__ == "__main__":
    args       = get_arguments()
    # Operon     = args.Operon
    GeneBlock  = args.GeneBlock
    GeneName   = args.GeneName
    Outfile    = args.Outfile
    Reference  = args.Reference
    uniprot_dic = parse_reference(Reference)
    # print 'uniprot_dic',uniprot_dic
    # operon_dic = parse_operon(Operon)
    geneBlock_dic = parse_gene_block_names_and_genes(GeneBlock)
    # print 'geneBlock_dic',geneBlock_dic
    return_dic_bsu(geneBlock_dic,uniprot_dic,Outfile,GeneName)
    
