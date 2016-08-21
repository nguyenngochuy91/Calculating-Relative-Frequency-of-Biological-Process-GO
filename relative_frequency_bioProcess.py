#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : For each B.Sub operon, get its gene (in bsu-BSU format)
    Start   : 08/18/2016
    End     : 08/20/2016
'''
import csv
###############################################################################
# Quick parsers to get info from files
###############################################################################
# parse the gene_block_names_and_genes.txt file
def parse_gene_block_names_and_genes(myfile):
    infile = open(myfile)
    mydic={}
    for line in infile.readlines():
        line = line.strip()
        items = line.split('\t')
        # print items
        mydic[items[0]] = items[1:]
    return mydic
# parse the B.Subtilis operon txt file
def parse_operon(myfile):
    infile = open(myfile)
    mydic={}
    # ignore the first line
    for line in infile.readlines()[0:]:
        line = line.split('\t')
        mydic[line[1]] = line[2]
    return mydic


# using the 2 dictionary to write out new operon file, where each line starts with 
# starter gene name, follow by what gene in the operon in bsu-BSU format
def return_dic_bsu(dic1,dic2):
    newdic={}
    for item in dic1:
        # initiate a list
        newdic[item]=[]
        for gene in dic1[item]: # iterate through gene annotation
            newdic[item].append(dic2[gene]) #add the formal format into the newdic
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
    
'''@function: using go term dic to pull out the biological process go term for each gene (assumption is that operons have different genes) 
              from assumption, we can say that number of biological processes for all operon is equal
              to total number of biological processes for all genes in final dic
   @input   : dic, key is they gene, value is all its go term
   @output  : 4 parameters: 
              1. dic(key is a gene, value is its biological process) Ex: {'BSU40410':['GO:0000160','GO:0006355']}
              2. total_count of all biological process go term,
              3. dic (key is a specific biological process, value is its count in all operon) Ex: {'GO:0006351':3, 'GO:0006355':4}
              4. dic (key is GO term, value is its biological process) Ex: {'GO:0006351':'P:transcription, DNA-templated'}
''' 
def get_biological_process_and_count(final_dic):
    count = 0
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
                count +=1 
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
                
    return gene_BioProcess_dic,count,GO_all_count,GO_BioProcess_dic
    
    
'''@function: from a list of operon, get the count of all biological process in each operon,
              as well as the count of each biological process in those operons.
   @input   : list of operon, newdic (key: operon, value: list of genes), gene_BioProcess_dic(key:gene,value: biological process go term)
   @output  : GO_local_count (key: biological go term, value: count),local_total_count
''' 
def list_of_10_biological_process_and_count(operon_list,newdic,gene_BioProcess_dic):
    GO_local_count = {} # Ex: {'GO:0006351':3, 'GO:0006355':4}
    local_total_count = 0
    for operon in operon_list:
        # get the gene in the operon from newdic[operon]
        for gene in newdic[operon]:
            # get the biological process in each gene from gene_BioProcess_dic[gene]
                for bioProcess in gene_BioProcess_dic[gene]:
                    # increment the top10 or bottom10 biological process count
                    local_total_count += 1
                    # initiate the count for that bioProcess in GO_local_count as 1, or increment it 
                    if bioProcess in GO_local_count:
                        GO_local_count[bioProcess] += 1
                    else:
                        GO_local_count[bioProcess] = 1
                        
    return GO_local_count,local_total_count
'''@function: for each go term that appear in either top 10 or bottom 10, get it relative frequency. The formula is as follow:
              1. Get the count of the go term in top 10, divide it by the total of number go term in all top10. Set this as local frequency
              2. Get the count of the go term in all operon, divide it by the total of number go term in all operon. Set this as global frequency
              3. Divide local by global. get this as the relative frequency
   @input   : GO_local_count,local_total_count, total_BioProcess_count,GO_all_count
   @output  : relative_frequency_dic (key : bioProcess, value: relative frequency)
'''     
def get_relative_freq(GO_local_count,local_total_count, total_BioProcess_count,GO_all_count):
    relative_frequency_dic = {} #Ex: {'GO:1900192': .9031}
    for GO_term in GO_local_count:  
        local_freq = GO_local_count[GO_term]/local_total_count
        global_freq = GO_all_count[GO_term]/total_BioProcess_count
        relative_freq = local_freq/global_freq
        relative_frequency_dic[GO_term] = relative_freq
    return relative_frequency_dic
    
'''@function: from a list of operon, get the count of all biological process in each operon,
              as well as the count of each biological process in those operons.
   @input   : list of operon, newdic (key: operon, value: list of genes), gene_BioProcess_dic(key:gene,value: biological process go term)
   @output  : GO_local_count (key: biological go term, value: count),local_total_count
''' 
def writting_csv(relative_frequency_dic,GO_BioProcess_dic,outfile):
    with open(outfile, 'w') as csvfile:
        fieldnames = ['GO_term', 'Biological_process','Relatively_frequency']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for key in sorted(relative_frequency_dic,key=relative_frequency_dic.get):
            
            writer.writerow({'GO_term': key,
                             'Biological_process': GO_BioProcess_dic[key],
                             'Relatively_frequency': round(relative_frequency_dic[key],2)})
        
###############################################################################
# execute the pipeline
###############################################################################
if __name__ == "__main__":
    dic1 = parse_gene_block_names_and_genes('gene_block_names_and_genes.txt')
    dic2 = parse_operon('B.Subtilis_Operons_ProOpDB.txt')
    
    # important dic to know which gene in an operon
    newdic = return_dic_bsu(dic1,dic2) # ex: 'bsub-BSU40410': ['BSU40370',  'BSU40380',  'BSU40360',  'BSU40410',  'BSU40390',  'BSU40400']

    # write into a file in order to query uniprot file for each gene
    outfile = open("name_bsucyc_uniprot","w")
    for operon in newdic:
        for gene in newdic[operon]:
            outfile.write('BSUB:'+gene+'-MONOMER'+'\n')
    outfile.close()
    # the uniprot file after getting from the internet is uniprot.txt (manually)
    # getting the go term for each gene
    
    # important dic that stores go term for each gene.
    final_dic = getting_go('uniprot.txt') # ex: 'BSU40390': ['GO:0016021,C:integral component of membrane','GO:0005886,C:plasma membrane']
    
    
    # getting the biological process Go term, the total bioP count, the count for each bioP term as a dic,
    # and the mapping from a go term and its biological process
    gene_BioProcess_dic,total_BioProcess_count,GO_all_count,GO_BioProcess_dic = get_biological_process_and_count(final_dic)
    # print ("total_BioProcess_count",total_BioProcess_count)
    # print ("GO_BioProcess_dic",GO_BioProcess_dic)
    
    
    # getting conservation score for each operon
    score = getting_conservation_score('conservedOperonsSorted.txt') # ex: 'bsub-BSU40410': '1.6647'
    # from the score dictionary, get the top 10 conserved and bottom 10 conserved operon
    top10,bottom10 = getting_top_bottom(score)
    # print (len(bottom10))
    # using the top10, bottom10 operon to get the genes in each operon.
    # from these genes, get the go term for each, then calculate the frequency
    
    # get the top10 info:
    top10_GO_local_count,top10_local_total_count = list_of_10_biological_process_and_count(top10,newdic,gene_BioProcess_dic)
    # get relative frequency dic for each go term of the top10 operon
    top10_relative_frequency_dic = get_relative_freq(top10_GO_local_count,top10_local_total_count, total_BioProcess_count,GO_all_count)
    # get the bottom10 info
    bottom10_GO_local_count,bottom10_local_total_count = list_of_10_biological_process_and_count(bottom10,newdic,gene_BioProcess_dic)
    # get relative frequency dic for each go term of the bottom10 operon
    bottom10_relative_frequency_dic = get_relative_freq(bottom10_GO_local_count,bottom10_local_total_count, total_BioProcess_count,GO_all_count)
    
    # writting the relative into csv file, 1st column is go term, 2nd column is its bioprocess, 3rd column is relative frequency
    writting_csv(top10_relative_frequency_dic,GO_BioProcess_dic,'top10_conserved.csv')
    writting_csv(top10_relative_frequency_dic,GO_BioProcess_dic,'bottom10_conserved.csv')
