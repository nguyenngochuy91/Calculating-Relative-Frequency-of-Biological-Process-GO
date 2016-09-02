#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : For each B.Sub operon, get its gene (in bsu-BSU format)
    Start   : 08/18/2016
    End     : 08/20/2016
'''

from __future__ import division
import csv
import argparse
from Bio.Ontology.IO import OboIO
from Bio.Ontology.Data import OntologyGraph
import time
#gloval variable
is_a ='is_a'
part_of = 'part_of'
MF_root = 'GO:0003674'
BP_root = 'GO:0008150'
CC_root = 'GO:0005575'
genes_with_no_GO_Term = ['BSU23401','BSU22440']
# get the arguments from command line
def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--Operon","-i",help="File to parse into newresult dic(operons_genes.txt)")
    parser.add_argument("--Uniprot","-u", help="Uniprot genes file(uniprot.txt)")
    parser.add_argument("--Score","-s", help="Conservation score of the operon(conservedOperonsSorted.txt)")
    parser.add_argument("--Level","-l", help="At which level (depth) of the BioProcess Ontology that we want to compare")
    parser.add_argument("--GO","-g", help="go-basic.obo file")
    parser.add_argument("--Method","-m", help="Choose either level method or leaf method ('leaf' vs 'level')")
    args = parser.parse_args()
    return args
###############################################################################
# Quick parsers to get info from files
###############################################################################
def return_dic_bsu(operons_genes):
    newdic={}
    infile = open(operons_genes,'r')
    for line in infile.readlines():
        line = line.strip()
        line = line.split('\t')
        newdic[line[0]] = line[1:]
    return newdic

'''@function: reading the uniprot file for each gene and get the go term into
              a dictionary
   @input   : text file
   @output  : dictionary, key is gene name, value is go term
'''   
def getting_go(myfile):
    open_file= open(myfile,'r')
    big_line = open_file.read()
    separated = big_line.split('//')
    new_list=[]
    for item in separated:
        new_list.append(item.split('\n'))
    final_dic={}
    for gene in new_list:
        for item in gene:
            if item[:17] == 'DR   BioCyc; BSUB':
                modified = item.split(' ')
                key = modified[4].replace(';','')
                key= key.split(':')[1]
                key = key.split('-')[0]       
                final_dic[key]=[]
            if item[:7]== 'DR   GO':
                modified = item.split(';')
                go_term = modified[1].strip(' ')+','+modified[2].strip(' ')
                final_dic[key].append(go_term)
    return final_dic
    
'''@function: from the conserved operon sorted text file, create a dictionary              
   @input   : text file
   @output  : dictionary, key is gene name, value is conservation score
''' 
def getting_conservation_score(myfile):
    open_file= open(myfile,'r')
    score = {}
    for line in open_file.readlines()[1:]:
        line = line.split('\t')
        score[line[0]]= line[1]
    return score
    
'''@function: from the conserved operon sorted text file, get 2 list, 1 is 
               top 10 conserved, other 1 is bottom 10 conserved
   @input   : dic
   @output  : 2 list
''' 
def getting_top_bottom(score):
    sorted_dic = sorted(score,key=score.get)
    top10 = sorted_dic[:10]
    bottom10 = sorted_dic[-10:]
    return top10,bottom10
    
    
###############################################################################
# main functions to filter out the biological process and calculating
# relative frequency for go term from either top10 or bottom10 operon,
# and writing function into a csv
###############################################################################

###############################################################################
# ** Method: Leaf
###############################################################################
'''@function: From the go-basic.obo file, parse it into a graph dic, then filter
              into a graph that only contains 'is_a' relationship
   @input   : go_basic.obo
   @output  : OntologyGraph 
''' 





###############################################################################
# ** Method: Level
###############################################################################
## filter out biological process, getting those at level 2 or user input level
###############################################################################    
'''@function: From the go-basic.obo file, parse it into a graph dic, then filter
              into a graph that only contains 'is_a' relationship
   @input   : go_basic.obo
   @output  : OntologyGraph 
''' 
def filter_graph(GO):
    # parsing the go-basic.obo file into a big dic
    go_graph = OboIO.OboReader(open(GO)).read()
    # filtering that go into an only is_a relationship go
    fgraph = OntologyGraph()
    for label,node in go_graph.nodes.items():
        fgraph.update_node(label,node.data) # get the info from the initial graph
        for edge in node.succ:
            # check if edge.data ='is_a'
            if edge.data == is_a:
                # add it into the new filter graph
                fgraph.add_edge(label,edge.to_node.label,edge.data)
    fgraph.synonyms = dict(go_graph.synonyms)
    for rel,typedef in go_graph.typedefs.items():
        if rel == is_a:
            fgraph.typedefs[rel] = typedef
    return fgraph
    
'''@function: From the filter graph with only is_a relationship, get the 
              bioprocesses at level given from the user
   @input   : OntologyGraph (fgraph)
   @output  : list of bioProcess go terms 
''' 
def search_level(fgraph,Level):
    tree = {}
    tree[1]=[BP_root]
    # iterate through the depth of the graph from the root
    for index in range(2,Level+1):
        tree[index]=[]
        for parent in tree[index-1]: # get to the upper level to find the predessor
            node = fgraph.nodes[parent]
            # check the children of this node
            for edge in node.pred:
                # add the children label into tree[index] list
                tree[index].append(edge.to_node.label)
    # print tree
    return tree[Level]
     
'''@function: From the level_bioProcess list, get the biological process of each term
   @input   : OntologyGraph (fgraph),level_bioProcess
   @output  : GO_BioProcess_level_dic
''' 
def retrieve_BP_at_level(fgraph,level_bioProcess):
    GO_BioProcess_level_dic = {}
    for GO in level_bioProcess:
        node =  fgraph.nodes[GO]
        GO_BioProcess_level_dic[GO] = node.data.name
    return GO_BioProcess_level_dic
    
'''@function: using go term dic to pull out the biological process go term for each gene (assumption is that operons have different genes) 
              from assumption, we can say that number of biological processes for all operon is equal
              to total number of biological processes for all genes in final dic
   @input   : dic, key is they gene, value is all its go term
   @output  : 3 parameters: 
              1. dic(key is a gene, value is its biological process) Ex: {'BSU40410':['GO:0000160','GO:0006355']}
              2. dic (key is a specific biological process, value is its count in all operon) Ex: {'GO:0006351':3, 'GO:0006355':4}
              3. dic (key is GO term, value is its biological process) Ex: {'GO:0006351':'P:transcription, DNA-templated'}
''' 
def get_biological_process_and_count(final_dic):
    gene_BioProcess_dic ={} #Ex {'BSU40410':['GO:0000160','GO:0006355']}
    GO_all_count ={} # Ex: {'GO:0006351':3, 'GO:0006355':4}
    GO_BioProcess_dic = {} # Ex: {'GO:0006351':'P:transcription, DNA-templated'}
    for gene in final_dic:
        # initiate the list for biological process
        gene_BioProcess_dic[gene] =[]
        # go through the go term
        for GO in final_dic[gene]:
            GO = GO.split(',')
            GO_term = GO[0]
            purpose = GO[1]
            # check if the go is a biological process:
            if purpose[0] == 'P':
                # increment the total count
                # append this go term into the value of key gene in gene_BioProcess_dic
                gene_BioProcess_dic[gene].append(GO_term)
                # increment, or initiate the count of this GO_term as 1 in the GO_all_count
                if GO_term in GO_all_count:
                    GO_all_count[GO_term] += 1
                else:
                    GO_all_count[GO_term] = 1
                
                # add the go term and its purpose into the GO_BioProcess_dic if the go term is not in the dic already:
                if GO_term not in GO_BioProcess_dic:
                    GO_BioProcess_dic[GO_term] = purpose
                
    return gene_BioProcess_dic,GO_all_count,GO_BioProcess_dic
    
'''@function: For each gene, find their biological process at level 2
   @input   : newdic, gene_BioProcess_dic,level_bioProcess,fgraph
   @output  : gene_BioProcess_at_level_dic
'''
def get_gene_BP_at_level(gene_BioProcess_dic,level_bioProcess,fgraph):
    gene_BioProcess_at_level_dic ={}
    for gene in gene_BioProcess_dic:
        level2_BP=set() # have this as a set since we don't want duplication 
        for BP in gene_BioProcess_dic[gene]:
            if BP in level_bioProcess:
                level2_BP.add(BP)
            else:
                ancestors = fgraph.get_ancestors(BP)
                for ancestor in ancestors:
                    if ancestor in level_bioProcess:
                        level2_BP.add(ancestor) 
        gene_BioProcess_at_level_dic[gene] = level2_BP
    return gene_BioProcess_at_level_dic
    
'''@function: For each operons, find the count of each Biological Process at level 2 in the operon
   @input   : newdic, gene_BioProcess_at_level_dic
   @output  : operon_BP_dic(key: operon name, value: [BP go terms dic count (key: go term, value: count),number of gene in the operon])
''' 
def get_BioProcess_from_operon(newdic,gene_BioProcess_at_level_dic):
    operon_BP_dic = {} # Ex : {'bsub-BSU40410': [{'GO:0000160':4,'GO:0006355':3},5]}
    for operon in newdic: # iterate through the operon names
        length = len(newdic[operon])
        BP_dic = {} # a local dictionary to keep count of each biological process in the operon
        for gene in newdic[operon]: # iterate through the gene in the operon
            if gene in genes_with_no_GO_Term: # those genes were not annotated with any go term
                length -=1
            else:
                for process in gene_BioProcess_at_level_dic[gene]: # iterate through each process in the gene
                    if process in BP_dic:
                        BP_dic[process] +=1
                    else:
                        BP_dic[process] = 1
        operon_BP_dic[operon] = [BP_dic,length]
    return operon_BP_dic
    
'''@function: testing whether there is a go count that is more than number of genes
              in the operon
   @input   : operon_BP_dic
   @output  : True (if pass the test), False (otherwise, give which operon and term that fails)
''' 
def test_operon_BP_dic(operon_BP_dic):
    check = True
    failed = {}
    for operon in operon_BP_dic:
        failed[operon] =[]
        info = operon_BP_dic[operon]
        length = info[1]
        GO = info[0]
        for term in GO:
            if GO[term]>length:
                check = False
                failed[operon].append(term)
    if check: 
        return check
    else:
        return check,failed
                
'''@function: Using the conservation score operon file to filter out irrelevant
              operon
   @input   : operon_BP_dic,score
   @output  : filter_operon_BP_dic
'''             
def filter_by_score(operon_BP_dic,score):
    filter_operon_BP_dic ={}
    for operon in score:
        filter_operon_BP_dic[operon]= operon_BP_dic[operon]
    return filter_operon_BP_dic
    
###############################################################################
## calculating relative frequency for go term from either top10 or bottom10 operon,
## and writing function into a csv
############################################################################### 
'''@function: from filter_operon_BP_dic, get the total number of BP count, and a dictionary
              of each BP count.
   @input   : filter_operon_BP_dic
   @output  : GO_glocal_count (key: biological go term, value: count),global_total_count
'''    
def get_global_count(filter_operon_BP_dic):
    GO_glocal_count ={}
    global_total_count = 0
    for operon in filter_operon_BP_dic:
        BP_dic = filter_operon_BP_dic[operon][0]
        for GO in BP_dic:
            # increment global_total_count
            global_total_count += BP_dic[GO]
            # update the GO_global_count dic
            if GO in GO_glocal_count:
                GO_glocal_count[GO] += BP_dic[GO]
            else:
                GO_glocal_count[GO] = BP_dic[GO]
    return GO_glocal_count,global_total_count
    
'''@function: from a list of operon, get the count of all biological process in each operon,
              as well as the count of each biological process in those operons.
   @input   : list of operon, filter_operon_BP_dic
   @output  : GO_local_count (key: biological go term, value: count),local_total_count
'''
def list_of_10_biological_process_and_count(operon_list,filter_operon_BP_dic):
    GO_local_count = {} # Ex: {'GO:0006351':3, 'GO:0006355':4}
    local_total_count = 0
    for operon in operon_list:
        # get the GO TERM in the operon from filter_operon_BP_dic
        info = filter_operon_BP_dic[operon]
        BP_dic = info[0]
        for GO in BP_dic:
            # increment local_total_count
            local_total_count += BP_dic[GO]
            if GO in GO_local_count:
                GO_local_count[GO] += BP_dic[GO]
            else:
                GO_local_count[GO] = BP_dic[GO]                       
    return GO_local_count,local_total_count
    
'''@function: for each go term that appear in either top 10 or bottom 10, get it relative frequency. The formula is as follow:
              1. Get the count of the go term in top 10, divide it by the total of number go term in all top10. Set this as local frequency
              2. Get the count of the go term in all operon, divide it by the total of number go term in all operon. Set this as global frequency
              3. Divide local by global. get this as the relative frequency
   @input   : GO_local_count,local_total_count, global_total_count,GO_glocal_count
   @output  : relative_frequency_dic (key : bioProcess, value: relative frequency)
'''     
def get_relative_freq(GO_local_count,local_total_count, global_total_count,GO_glocal_count):
    relative_frequency_dic = {} #Ex: {'GO:1900192': .9031}
    for GO_term in GO_local_count:  
        local_freq = GO_local_count[GO_term]/local_total_count
        global_freq = GO_glocal_count[GO_term]/global_total_count
        relative_freq = local_freq/global_freq
        relative_frequency_dic[GO_term] = relative_freq
    return relative_frequency_dic
    
'''@function: given comparison dic and the relative frequency dic, adding item into the 
              comparison dic
   @input   : relative_frequency_dic,comparison_dic
   @output  : comparison_dic
'''        
def helper_add(comparison_dic,relative_frequency_dic,choice):
    for GO in relative_frequency_dic:
        if relative_frequency_dic[GO] <1:
            comparison_dic['under'][choice].add(GO)
        elif relative_frequency_dic[GO] == 1:
            comparison_dic['normal'][choice].add(GO)
        else:
            comparison_dic['over'][choice].add(GO)
    return comparison_dic

    
'''@function: given top10_relative_frequency_dic,bottom10_relative_frequency_dic,
              categorize based on frequency (over/under/normal representative)
   @input   : top10_relative_frequency_dic,bottom10_relative_frequency_dic
   @output  : comparison_dic
'''        
def compare(top10_relative_frequency_dic,bottom10_relative_frequency_dic):
    comparison_dic ={'under' :{'top':set(),'bot':set()},
                     'over'  :{'top':set(),'bot':set()},
                     'normal':{'top':set(),'bot':set()}}
    # update using top10
    comparison_dic = helper_add(comparison_dic,top10_relative_frequency_dic,'top')
    # update using bot10
    comparison_dic = helper_add(comparison_dic,bottom10_relative_frequency_dic,'bot')
    
    # check the one that are the same 
    for level_representative in comparison_dic:
        comparison_dic[level_representative]['same']=  comparison_dic[level_representative]['top'].intersection(comparison_dic[level_representative]['bot'])

    # from the same, remove it from top and bottom
    for level_representative in comparison_dic:
        for GO in comparison_dic[level_representative]['same']:
            comparison_dic[level_representative]['top'].remove(GO)
            comparison_dic[level_representative]['bot'].remove(GO)
    return comparison_dic
    
'''@function: from a list of operon, get the count of all biological process in each operon,
              as well as the count of each biological process in those operons.
   @input   : list of operon,GO_BioProcess_level_dic, gene_BioProcess_dic(key:gene,value: biological process go term)
   @output  : csv file
''' 
def writting_csv(relative_frequency_dic,GO_BioProcess_level_dic,outfile):
    with open(outfile, 'w') as csvfile:
        fieldnames = ['GO_term', 'Biological_process','Relatively_frequency']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for key in sorted(relative_frequency_dic,key=relative_frequency_dic.get):
            
            writer.writerow({'GO_term': key,
                             'Biological_process': GO_BioProcess_level_dic[key],
                             'Relatively_frequency': round(relative_frequency_dic[key],2)})
def to_string(myset,GO_BioProcess_level_dic):
    mystring = ''
    for item in myset:
        mystring += item+':'+GO_BioProcess_level_dic[item]+ '\t'
    return mystring
                    
'''@function: from comparison_dic write out into csv file
   @input   : comparison_dic,GO_BioProcess_level_dic,outfile
   @output  : csv file
''' 
def writting_comparison(comparison_dic,GO_BioProcess_level_dic,outfile):
    with open(outfile, 'w') as csvfile:
        fieldnames = ['representative_level','top10_only','bottom10_only','both']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for key in comparison_dic:
            
            writer.writerow({'representative_level': key,
                             'top10_only': to_string(comparison_dic[key]['top'],GO_BioProcess_level_dic),
                             'bottom10_only': to_string(comparison_dic[key]['top'],GO_BioProcess_level_dic),
                             'both':to_string(comparison_dic[key]['same'],GO_BioProcess_level_dic)})

'''@function: writting out top10 or bottom10 biological process
   @input   : list of go term
   @output  : txt file
''' 
def writting_bioProcess(relative_frequency_dic,outfile):
    out = open(outfile,'w')
    for GO_term in relative_frequency_dic:
        out.write(GO_term+'\n')
    out.close()
###############################################################################
# execute the pipeline
###############################################################################
if __name__ == "__main__":
    # test again
    fgraph = filter_graph('../go-basic.obo')  
    newdic = return_dic_bsu('operons_genes.txt')
    final_dic = getting_go('uniprot.txt')
    score = getting_conservation_score('conservedOperonsSorted.txt')
    gene_BioProcess_dic,GO_all_count,GO_BioProcess_dic = get_biological_process_and_count(final_dic)
    # get the operon_dic that has key as operon name, value is dictionary (key is BP at level, value is count)
    operon_BP_dic =get_BioProcess_from_operon(newdic,gene_BioProcess_dic) 
    # filter our operon_BP_dic using score
    filter_operon_BP_dic = filter_by_score(operon_BP_dic,score)
    # real code
    args    = get_arguments()
    Operon  = args.Operon
    Uniprot = args.Uniprot
    Score   = args.Score
    GO      = args.GO
    Level   = int(args.Level)
    Method  = args.Method
    start   = time.time()
    
    # filter the graph to only get the is_a relationship edge:    
    fgraph = filter_graph(GO)
    
    # getting conservation score for each operon
    score = getting_conservation_score(Score) # ex: 'bsub-BSU40410': '1.6647'
    
    # from the score dictionary, get the top 10 conserved and bottom 10 conserved operon
    top10,bottom10 = getting_top_bottom(score)    
    
    # important dic to know which genes are in an operon
    newdic = return_dic_bsu(Operon) # ex: 'bsub-BSU40410': ['BSU40370',  'BSU40380',  'BSU40360',  'BSU40410',  'BSU40390',  'BSU40400']
    
    # important dic that stores go term for each gene.
    final_dic = getting_go(Uniprot) # ex: 'BSU40390': ['GO:0016021,C:integral component of membrane','GO:0005886,C:plasma membrane']
    
    # getting the biological process Go term, the count for each bioP term as a dic,
    # and the mapping from a go term and its biological process
    
    gene_BioProcess_dic,GO_all_count,GO_BioProcess_dic = get_biological_process_and_count(final_dic)
    if Method.lower() == 'level':
        # given the filter graph, find all the biological process that is at level Level
        level_bioProcess = search_level(fgraph,Level)
        # retrieve the biological process of the term at level
        GO_BioProcess_level_dic = retrieve_BP_at_level(fgraph,level_bioProcess)
        
        # from GO_BioProcess_dic, for each gene, find the BP that is at level specified by user
        gene_BioProcess_dic = get_gene_BP_at_level(gene_BioProcess_dic,level_bioProcess,fgraph)
        # get the operon_dic that has key as operon name, value is dictionary (key is BP at level, value is count)
        operon_BP_dic =get_BioProcess_from_operon(newdic,gene_BioProcess_dic) 
        #Ex:{'bsub-BSU40410': [{'GO:0008152': 1,'GO:0009987': 1,'GO:0050789': 1,'GO:0065007': 1},6]}
        # print test_operon_BP_dic(operon_BP_dic)
    

        # filter our operon_BP_dic using score
        filter_operon_BP_dic = filter_by_score(operon_BP_dic,score)
        # get the GO_global_coutn dic and total number of BP at level
        GO_glocal_count,global_total_count = get_global_count(filter_operon_BP_dic)
        
    elif Method.lower() =='leaf':
        # get the operon_dic that has key as operon name, value is dictionary (key is BP at level, value is count)
        operon_BP_dic =get_BioProcess_from_operon(newdic,gene_BioProcess_dic) 
        # filter our operon_BP_dic using score
        filter_operon_BP_dic = filter_by_score(operon_BP_dic,score)
    # using the top10, bottom10 operon to get the genes in each operon.
    # from these genes, get the go term for each, then calculate the frequency
    
    # get the top10 info:
    top10_GO_local_count,top10_local_total_count = list_of_10_biological_process_and_count(top10,filter_operon_BP_dic)
    # get relative frequency dic for each go term of the top10 operon
    top10_relative_frequency_dic = get_relative_freq(top10_GO_local_count,top10_local_total_count, global_total_count,GO_glocal_count)
    # get the bottom10 info
    bottom10_GO_local_count,bottom10_local_total_count = list_of_10_biological_process_and_count(bottom10,filter_operon_BP_dic)
    # get relative frequency dic for each go term of the bottom10 operon
    bottom10_relative_frequency_dic = get_relative_freq(bottom10_GO_local_count,bottom10_local_total_count, global_total_count,GO_glocal_count)
    
    # comparison between the top10 and bottom10    
    comparison_dic = compare(top10_relative_frequency_dic,bottom10_relative_frequency_dic)
    # print comparison_dic
    # writting the relative into csv file, 1st column is go term, 2nd column is its bioprocess, 3rd column is relative frequency
    writting_csv(top10_relative_frequency_dic,GO_BioProcess_level_dic,'./Result/top10_conserved_level'+args.Level+'.csv')
    writting_csv(bottom10_relative_frequency_dic,GO_BioProcess_level_dic,'./Result/bottom10_conserved_level'+args.Level+'.csv')
    # writting bioProcess term from relative frequency dic into txt file
    writting_bioProcess(top10_relative_frequency_dic,'./Result/top10_conserved_level'+args.Level+'.txt')
    writting_bioProcess(top10_relative_frequency_dic,'./Result/bottom10_conserved_level'+args.Level+'.txt')
    # writting comparison file
    writting_comparison(comparison_dic,GO_BioProcess_level_dic,'./Result/comparison'+args.Level+'.csv')
    stop = time.time()
    print stop-start
