#### Setting up
##Import modules
#import os module for file manipulations
import os
#import regular expression
import re
#import random module
import random

##create a function to get the reverse compliemnt from a sequence
def reverse_compliment(sequence):
    temp = sequence[::-1]# revert the sequence. from last to first
    temp = temp.replace("A", "t").replace("T", "a").replace("G", "c").replace("C", "g")#replace capital DNA with small letter complement of DNA to avoid duouble replacements with captial.
    temp = temp.upper()# capitalize string
    return temp


##create a function to analyse a sat of randomly selected genes from the gff3
def get_random_genes_analysis():
    random.seed()       
    random_genes = random.choices(gff3_genes, k = (len(supplied_genes_seqanalysis.dict_gene_promoter_region)))# select k (length of supplied genes that were in gff3) times genes from gff3_genes
    analysis = seqanalysis(random_genes)#seqanalysis class introduced on line 
    return analysis.dict_promoter_count #return a dictionary of promoters and count from analysis












#### Get inputs
##Get genes file
while True: #function should run until broken
    genes_file = input("What is the name of the genes file. Include file type, e.g .txt\n")#get input
    if genes_file.find(".txt") > -1:#if input has .txt
        try:#try openeing the file and saving in object
            opened_genefile = open(genes_file)
            break
        except:#if fail to open input
            print("Sorry, text file no valid!")#print error message and continue
            continue
    else:#if input does not contain .txt
        print("Sorry file needs to end with .txt!")
        continue 
##get promoters file
while True:
    promoters_file = input("What is the name of the promoters file. Include file type, e.g .txt\n")
    if promoters_file.find(".txt") > -1:
        try:
            opened_promoterfile = open(promoters_file)
            break
        except:
            print("Sorry, text file no valid!")
            continue
    else:
        print("Sorry file needs to end with .txt!")
        continue 
#get name of the folder with gff3 and fasta files.
folder_wgff3_and_fasta = input("What is the name of the folder with gff3 and fasta files.\n")
 











#### Get information from files
##Get information from gene file
genefile_lines = opened_genefile.readlines()
opened_genefile.close()
#Get names of genes we are looking for
wanted_genes = []
for i in genefile_lines:
    wanted_genes.append(i.rstrip())


##get information from promoters file
promoterfile_lines = opened_promoterfile.readlines()
promoters = []
for promoter in promoterfile_lines:
    promoters.append(promoter.rstrip().upper())
opened_promoterfile.close()


##get information from fasta file
for file in os.listdir(folder_wgff3_and_fasta):#for file in directory
    if file.find(".fa") > -1:#if file name has .fa
        opened_fasta = open(folder_wgff3_and_fasta + "/" + file)#open file. need to provide full directory in argument
        fasta_lines = opened_fasta.readlines()#read lines of file to list
        opened_fasta.close() #close opened_fasta
#get fasta sequence
fasta_lines_newline_stripped = []
for i in range(1, len(fasta_lines)):
    fasta_lines_newline_stripped.append(fasta_lines[i].rstrip())
fasta_sequence = "".join(fasta_lines_newline_stripped).upper() 





##Get information from gff3 file
for file in os.listdir(folder_wgff3_and_fasta): #for file in directory
    if file.find(".gff3") > -1:#if file name has .gff3
        opened_gff = open(folder_wgff3_and_fasta + "/" + file) #open file. need to provide full directory in argument
        gff_lines = opened_gff.readlines()#read lines of file to list
        opened_gff.close() #close opened_gff
# Make a dictionary with name of gene as key and a nested list with [start and end positions] and (positive or negative strand) as key
dict_gene_positions_strand = {} 
#get names of genes in gff3 file
gff3_genes = [] #get list of gff3 genes
for line in gff_lines: 
    line_as_list = line.split("\t")
    if not (re.search("^#", line_as_list[0])): #if line does not start with #
        if line_as_list[2] == "gene": 
            temp = line_as_list[8].split(";")
            temp = temp[0].split(":")
            gff3_gene = temp[1]
            gff3_genes.append(gff3_gene)
            position_list = []
            position_list.append(int(line_as_list[3]))#add start position to position list
            position_list.append(int(line_as_list[4]))#add end position
            sublist = []
            sublist.append(position_list)
            sublist.append(line_as_list[6])#add strand sense/antisense   
            dict_gene_positions_strand[gff3_gene] = sublist




#### create a class to analyse information
class seqanalysis(object):

    def __init__(self, genes_to_search): # genes_to_search is a list of the genes being searched against
        self.my_genes = genes_to_search #make my_genes an attribute
        self.dict_gene_promoter_region = {}#make dictionary of gene and the promoter region an attribute. the gene will be the key and promoter region as value
        for gene in self.my_genes:#for gene in interested genes
            if gene in gff3_genes:#if gene in gff3 file
                if dict_gene_positions_strand[gene][1] == "+": # if positive strand
                    self.dict_gene_promoter_region[gene] = fasta_sequence[dict_gene_positions_strand[gene][0][0]- 401 : dict_gene_positions_strand[gene][0][0] - 1]#region is 400nt before start point
                else:#if strand is negative
                    self.dict_gene_promoter_region[gene] = reverse_compliment(fasta_sequence[dict_gene_positions_strand[gene][0][1] : dict_gene_positions_strand[gene][0][1] + 400])#region is the reverse compliment oif 400nt after end point
                    if re.finditer("N{100,}", self.dict_gene_promoter_region[gene]):#if there are a string of 100Ns in promoter region
                        matches = re.finditer("N{100,}", self.dict_gene_promoter_region[gene])
                        new_promoterstart = 0
                        for match in matches: #for each match
                            new_promoterstart = match.end() #change the new promoter start sit to last match.end
                        self.dict_gene_promoter_region[gene] = self.dict_gene_promoter_region[gene][new_promoterstart : ]
        #make a dictionary of promoter counts for each gene
        self.dict_gene_and_promoter_counts = {}
        for gene in self.dict_gene_promoter_region:
            self.dict_gene_and_promoter_counts[gene] = []
            for promoter in promoters:
                cite_to_search = re.compile(promoter)#make site a searchable regex object
                temp = cite_to_search.findall(self.dict_gene_promoter_region[gene])
                list_of_promoter_and_count = []
                list_of_promoter_and_count.append(promoter)
                list_of_promoter_and_count.append((len(temp)))
                self.dict_gene_and_promoter_counts[gene].append(list_of_promoter_and_count)
        
        #make a dict of the counts for all the promoters in selected genes
        self.dict_promoter_count = {}
        for gene in self.dict_gene_and_promoter_counts:
            for i in range(len(self.dict_gene_and_promoter_counts[gene])):
                promoter = self.dict_gene_and_promoter_counts[gene][i][0]
                if promoter in self.dict_promoter_count:
                    self.dict_promoter_count[promoter] = self.dict_promoter_count[promoter] + self.dict_gene_and_promoter_counts[gene][i][1]
                else: 
                    self.dict_promoter_count[promoter] = self.dict_gene_and_promoter_counts[gene][i][1]













#### Analysis to address questions
##analysis supplied genes
supplied_genes_seqanalysis = seqanalysis(wanted_genes)






#list with dictionary of random genes analysis
list_of_dicts = [supplied_genes_seqanalysis.dict_promoter_count]
number_of_samples = 5
for i in range(number_of_samples):
    list_of_dicts.append(get_random_genes_analysis())







#### the print out
for promoter, count in supplied_genes_seqanalysis.dict_promoter_count.items():
    line_to_print = []
    line_to_print.append(promoter)
    for dict in list_of_dicts:
        for key, value in dict.items():
            if promoter == key:
                line_to_print.append(str(value))
    line_to_print = "\t".join(line_to_print)
    print(line_to_print)



