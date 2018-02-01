#Script written to filter data from MC3 MAF
#MC3 MAF Info: https://www.synapse.org/#!Synapse:syn7214402
#You may need to edit column info if modified(Gene Name, Filter Flags, MutationType, Callers)
#by: reykajayasinghe@gmail.com
#Edited on: 01/29/18

import re

og_maf = open("mc3.v0.2.8.CONTROLLED.maf")
#og_maf=open("test")

#Header
#header = og_maf.readline()

for line in og_maf:
	
	data = line.split('\t')
	skip_line = False
#Skip header
	if data[0] == 'Hugo_Symbol':
		skip_line = True

	if skip_line:
		print(line.strip())
		continue

#FILTER 1: Filter flags
#StrandBias, pcadontuse, common_in_exac, oxog, contest, nonpreferredpair, ndp, badseq, broad_PoN_v2
#filters not currently used:PASS, NonExonic, bitgt, gapfiller, native_wga_mix, wga
	filter_list = ['StrandBias','pcadontuse','common_in_exac','oxog','contest','nonpreferredpair','ndp','badseq','broad_PoN_v2']
	
	filters=data[108].split(',')	
	for f in filters:
		if f in filter_list:
			skip_line = True
			break
			
	if skip_line:
		continue

#FILTER 2 (CALLER FILTER)
#For SNPS only keep calls that are derived from two mutation callers
#For InDels only keep calls that are derived from two mutation callers or PINDEL only
	mutationtype = data[9]
	callers = data[110].split('|')
	total_callers = len(callers)
	mut_types = ['ONP','TNP','SNP']
	indel_types = ['INS','DEL']
#SNP Check
	if mutationtype in mut_types:
		if total_callers < 2:
			skip_line = True
	if skip_line:
		continue
#InDel Check
	if mutationtype in indel_types:
		if total_callers < 2:
#			for c in callers:
#				if c != 'PINDEL':
			skip_line = True

	if skip_line:
		continue

#FILTER 3 (VAF) > 0.05 (5%)
	
	t_depth = data[39]
	t_ref_count = data[40]
	t_alt_count = data[41]
	n_depth = data[42]
	n_ref_count = data[43]
	n_alt_count = data[44]
	
 
	VAF = float(t_alt_count)/float(t_depth)

	if VAF < 0.05:
		skip_line = True

	if skip_line:
		continue

	print(line.strip())	

og_maf.close()
