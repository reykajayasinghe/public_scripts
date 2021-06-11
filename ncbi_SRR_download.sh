#Reyka Jayasinghe
#06-11-2021

#References
#https://ncbi.github.io/sra-tools/fastq-dump.html
#https://www.biostars.org/p/156909/
#https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/#download-sequence-data-files-usi
#Example SRR: https://www.ncbi.nlm.nih.gov/sra/SRX5914447[accn]

#usage - bash ncbi_SRR_download.sh SRRID

SRR=$1

###########################
##########STEP1############
###########################
#Uncomment and run this command first - this will download the SRA file - takes an hour or two based on internet connection speed
#out=${SRR}.nohupout
#nohup prefetch -t ascp -X 200g ${SRR} > ${out} 2>&1&

###########################
##########STEP2############
###########################
#Uncomment and run this command second - this will unpack fastq files into different files based onr read pair

out=${SRR}_2.nohupout
nohup fastq-dump --split-files --gzip ${SRR}/${SRR}.sra --outdir ${SRR}/.> ${out} 2>&1&

###########################
##########STEP3############
###########################
#Samples will be saved into respective folders by the SRR and fastq will be separated by read file by _1, _2, _3...
#For running downstream scripts like cellranger on scRNA-seq data you will need to rename the read pair files in a specific format - see below
#https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/fastq-input
