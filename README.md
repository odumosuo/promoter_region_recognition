# Promoter region recognition

## Description
This project searches for known motifs in the upstream promoter regions of a set of coexpressed genes. The file zea_mays_genes.txt
contains a list of maize (Zea mays) gene identifiers for a group of co-regulated maize genes. The file promoters contains a list of experimentally identified transcription factor binding sites (TFBSs) that are active in maize. The script extracts the upstream non-coding sequence for each gene in zea_mays_genes.txt.

## How to run
GFF3 and fasta file should be in one folder. 
Run python script and answer prompts

## Features
Transcription starty site (TSS)  
Transcription factorbinding sites (TFBSs)    
Fasta file  
GTF/GFF file  
Extended regular expression  

## Dataset 
Fasta file containing the maize genome, version AGPv4 was obtained here: ftp.ensemblgenomes.org/pub/release-48/plants/fasta/zea_mays/dna  

GTF/GFF file for the same version to find the TSS for each gene was obtained here:
ftp.ensemblgenomes.org/pub/release-48/plants/gff3/zea_mays
