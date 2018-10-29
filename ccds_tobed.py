
#Reyka Jayasinghe (reyka@wustl.edu)
#Last edited: October 15th, 2018
#This script takes in a CCDS file like below and converts it to a bed like format
#Input CCDS (from UCSC table browser)
#1      NC_000001.10    SAMD11  148398  CCDS2.2 Public  +       861321  879532  [861321-861392, 865534-865715, 866418-866468, 871151-871275, 874419-874508, 874654-874839, 876523-876685, 877515-877630, 877789-877867, 877938-878437, 878632-878756, 879077-879187, 879287-879532]     Identical

#output
#1	861321	861392	SAMD11
#1	865534	865715	SAMD11
#1	866418	866468	SAMD11
#1	871151	871275	SAMD11
#1	874419	874508	SAMD11
#1	874654	874839	SAMD11
#1	876523	876685	SAMD11

#to download ccds 
#clade: Mammal
#genome: Human
#Assembly: Feb. 2009 (GRCh37/hg19)
#group: Genes and Gene Predictions
#track: CCDS
#save as .tsv

#Modules
import sys
import os

ccdsfile=open(sys.argv[1],'r')

for line in ccdsfile:
	(chromosome,nc_accession,gene,gene_id,ccds_id,ccds_status,cds_strand,cds_from,cds_to,cds_locations,match_type) = line.strip().split("\t")
	exons=cds_locations.replace("[","").replace("]","")
	exons_e=exons.strip().split(",")
	for exon in exons_e:
		(start,stop)=exon.strip().split("-")
		fin=chromosome+"\t"+start+"\t"+stop+"\t"+gene
		print(fin)


