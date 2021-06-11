#Reyka Jayasinghe
#Edited 10292018
#Takes in a bed file like below and splits each gene into an individual bed file

#CCDS input bed file below
#1	861321	861392	SAMD11
#1	865534	865715	SAMD11
#1	866418	866468	SAMD11
#1	871151	871275	SAMD11
#1	874419	874508	SAMD11
#1	874654	874839	SAMD11

#Import necessary packages
from collections import defaultdict
import sys
import gzip
import time

start = time.time()

file=sys.argv[1]
BED=open(file,"r")

for line in BED:
	data=line.strip().split("\t") 
	currentgene=data.pop()
	try:
		lastgene
	except NameError:
		###NEW OUT
		outputnew=currentgene+".bed"
		OUT=open(outputnew,'w')
		OUT.write(line)
	else:
		if currentgene == lastgene:
		#add line to current output file
			OUT.write(line)
		else:
		#Open new output file, close last file
			OUT.close()
			outputnew=currentgene+".bed"
			OUT=open(outputnew,'w')
			OUT.write(line)
	lastgene = currentgene

BED.close()
		
