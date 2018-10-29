
#Reyka Jayasinghe (reyka@wustl.edu)
#Edited: October 29th 2018
#Takes in a VCF file with mutliple samples and determines total number of samples with each mutation

#Output
#Needs to be modified. Currently wants input on what ancestry is needed to filter files that are counted.
#This is not integrated into the .gz version currently - only works for .vcf files

#VCF="PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389.vcf.gz"

#Import necessary packages
from collections import defaultdict
import sys
import gzip
import time

VCF=sys.argv[1]
ancestry_group=sys.argv[2]

#CLINICAL data for TCGA PANCANCER DATA
cancertype=defaultdict(dict)
ancestry=defaultdict(dict)
sampleindex=defaultdict(dict)
clin="PanCan_ClinicalData_V4_wAIM_filtered10389.txt"
if not ancestry_group.startswith('NONE'):
	for cline in open(clin,'r'):
	#read in clinical data and parse out group of interest
#bcr_patient_barcode	type	age_at_initial_pathologic_diagnosis	gender	race	ajcc_pathologic_tumor_stage	clinical_stage	histological_type	histological_grade	initial_pathologic_dx_year	menopause_status	birth_days_to	vital_status	tumor_status	last_contact_days_to	death_days_to	cause_of_death	new_tumor_event_type	new_tumor_event_site	new_tumor_event_site_other	new_tumor_event_dx_days_to	treatment_outcome_first_course	margin_status	residual_tumor	OS	OS.time	DSS	DSS.time	DFI	DFI.time	PFI	PFI.time	PFI.1	PFI.time.1	PFI.2	PFI.time.2	Redaction	self_identified_race	self_identified_ethnicity	consensus_call
#TCGA-02-0003	GBM	50	MALE	WHITE	[Not Available]	[Not Applicable]	Untreated primary (de novo) GBM	[Not Available]	2003	[Not Available]	-18341	Dead	WITH TUMOR	144	144	[Not Available]	Recurrence	[Not Applicable][Not Applicable]	40	[Not Available]	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	WHITE	NOT HISPANIC OR LATINO	eur	
	#Grab first entry = sample name and last entry = ancestry
		clindata=cline.strip().split('\t')
		sample=clindata.pop(0)
		c=clindata.pop(0)
		cancertype[sample]=c
		ancestry_call=clindata.pop()
		ancestry[sample]=ancestry_call

###IF GZ file use below
if VCF.endswith(".gz"):
	with gzip.open(VCF,'r') as fin:
		for vline in fin:
			if not vline.startswith(b'#'):
				data=vline.strip().split(b'\t')	
				#first few columns are sample info and the rest are variants per sample
				#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT	SAMPLE1	SAMPLE2	SAMPLE3...
				variant_count=0
				chrom=str(data.pop(0),'utf-8')
				pos=str(data.pop(0),'utf-8')
				ID=str(data.pop(0),'utf-8')
				REF=str(data.pop(0),'utf-8')
				ALT=str(data.pop(0),'utf-8')
				QUAL=str(data.pop(0),'utf-8')
				FILTER=str(data.pop(0),'utf-8')
				INFO=str(data.pop(0),'utf-8')
				FORMAT=str(data.pop(0),'utf-8')
				for sample in data:
					if not sample.startswith((b'./.:.:.:.:.:.:.:.:.:.:.:.:.:.',b'./.:.:.:.:.')):
						variant_count+=1
				print(chrom,pos,REF,ALT,variant_count)
	DATA.close()
else:
	FILE=open(VCF,'r')
	for line in FILE:
		if line.startswith('#CHROM'):
			#This is the header line store sample names and position
			samplelist=line.strip().split('\t')
			del samplelist[:9]
			indexval=0
			for sn in samplelist:
				snshort=sn[0:12]
				sampleindex[indexval]=snshort
				indexval+=1
		if not line.startswith('#'):
			data=line.strip().split('\t')
			variant_count=0
			chrom=data.pop(0)
			pos=data.pop(0)
			ID=data.pop(0)
			REF=data.pop(0)
			ALT=data.pop(0)
			QUAL=data.pop(0)
			FILTER=data.pop(0)
			INFO=data.pop(0)
			FORMAT=data.pop(0)
			vcfindexval=0
			for sample in data:
				vcfindexval+=1
				if not sample.startswith(('./.:.:.:.:.:.:.:.:.:.:.:.:.:.','./.:.:.:.:.')):
					#determine ancestry of individual
					actualval=vcfindexval-1
					samplecurrent=sampleindex[actualval]
					consensus_ancestry_call=ancestry[samplecurrent]
					if consensus_ancestry_call == "eur":
						variant_count+=1
			print(chrom,pos,REF,ALT,variant_count)
	VCF.close()

