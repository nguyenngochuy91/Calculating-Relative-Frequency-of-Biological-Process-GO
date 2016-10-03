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
    parser.add_argument("--GeneName","-n", help="Gene name to parse into uniprot to get protein file 'name_bsucyc_uniprot.txt'")
    parser.add_argument("--Reference","-r", help="Reference genome 'NC_000964.gbk'")
    parser.add_argument("--GI_locus","-l", help="Mapping between GI number and the locus in ase can't find by GI (GI_locus.txt)")
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
    GI_dic={}
    locus_dic ={}
    record = SeqIO.read(myfile,"genbank")
    for feature in record.features:
        if feature.type == 'CDS':
            # need the locus tag to query from the downloading file (name is different)
            locus_dic[feature.qualifiers['gene'][0]] = feature.qualifiers['locus_tag'][0]
            for item in feature.qualifiers['db_xref']:
                if item[:2]=='GI':
                    # get the GI number to find from uniprot
                    GI_dic[feature.qualifiers['gene'][0]] = item[3:]
                

    return GI_dic,locus_dic
# using the 2 dictionary to write out new operon file, where each line starts with 
# starter gene name, follow by what gene in the operon in bsu-BSU format
# in addition, also wrote out gene names in following format BSUB:BSU32510-MONOMER
def return_dic_bsu(geneBlock_dic,GI_dic,locus_dic,outfile,GeneName,GI_locus):
    outfile_name= open(outfile,'w')
    GeneName_file = open(GeneName,'w')
    GI_locus_file = open(GI_locus,'w')
    for item in geneBlock_dic:
        # initiate a list
        outfile_name.write(item)
        for gene in geneBlock_dic[item]: # iterate through gene annotation
            if gene in synonym_dic:
                gene = synonym_dic[gene]
            try:
                GeneName_file.write(GI_dic[gene]+'\n')
                outfile_name.write('\t'+locus_dic[gene])
                GI_locus_file.write(GI_dic[gene]+'\t'+locus_dic[gene]+'\n')
            except:      
                print gene
        outfile_name.write('\n')
    outfile_name.close()
    GeneName_file.close()
    GI_locus_file.close()
    
if __name__ == "__main__":
    args       = get_arguments()
    # Operon     = args.Operon
    GeneBlock  = args.GeneBlock
    GeneName   = args.GeneName
    Outfile    = args.Outfile
    Reference  = args.Reference
    GI_locus   = args.GI_locus
    GI_dic,locus_dic = parse_reference(Reference)
    # print 'uniprot_dic',uniprot_dic
    # operon_dic = parse_operon(Operon)
    geneBlock_dic = parse_gene_block_names_and_genes(GeneBlock)
    # print 'geneBlock_dic',geneBlock_dic
    return_dic_bsu(geneBlock_dic,GI_dic,locus_dic,Outfile,GeneName,GI_locus)
    
