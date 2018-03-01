use strict;

my $usage =<<USAGE;
 Usage: $0 <MAF FILE> <OUTPUT SITES DIR> <OUTPUT SAMPLES DIR>
  Example: perl /gscmnt/gc2706/dinglab/medseq/LabCode/Reyka/SpliceInator/SpliceInator/RNA_readcount.pl MAF sitesdat samplesdat

	Creates 4 output directories: sitesdat samplesdat readcountdat and logs.
		sitesdat:5 column format of all variants split by sample and chromosome (sample.chromosome.sites)
		samplesdat:MAF specific to each sample and chromosome (sample.chromsome.txt)
		readcountdat: Final readcount results (sample.chromosome.sites.rc[dna])
		logs: logs of bsub commands

All you need to run is your MAF INPUT. Input bams are already in script along with reference genomes to use for each chromosome.

USAGE
    die $usage unless @ARGV==3;
my (@array1,$sample,$file,$linest);
open(my $MAF,'<',$ARGV[0]) or die "INPUT MAF not found!";
print "Sorting MAF by sample name\n";
my %sample_hash;
my $counter=0;
my $sampledir=$ARGV[2];
my $sitesdir=$ARGV[1];
`mkdir -p $sampledir`;
`mkdir -p $sitesdir`;
`mkdir -p logs`;
#Some cancer types use chr prefix in the bam file and this changes whether or not "chr" is included in the sites file and it changes the rereference fasta used for determining readcount calculations
#This initial section will load in information for each sample derived from samtools view -H to determine whether or not "chr" needs to be used as a prefix.
open(my $RNAbamlist,'<',"/gscmnt/gc2600/dinglab/reyka/Unarchive_Samples/RNA_bampaths_052016"); 
open(my $NOBAM,'>',"NObamfile.log"); 
#TCGA-BT-A0S7-01A-11R-A10U-07 blca /gscmnt/gc8000/info/alignment_data/imported/2892223462/all_sequences.bam @SQ SN:chr1 LN:249250621
#TCGA-FD-A3B3-01A-12R-A206-07 blca /gscmnt/gc13000/info/alignment_data/imported/2892230272/all_sequences.bam @SQ SN:chr1 LN:249250621
#TCGA-GD-A3OP-11A-11R-A220-07 blca /gscmnt/ams1175/info/alignment_data/imported/2892231586/all_sequences.bam @SQ SN:chr1 LN:249250621
my %chrprefixhash;
my %bamhash;
while(my $bamline=<$RNAbamlist>){
	chomp $bamline;
	my @baminfo=split(/ /,$bamline);
	my $bsample=$baminfo[0];
#Skip all normal tissue samples
	if ($bsample=~/TCGA-\w[2]-\w[4]-1.*/){
		next;
	}
	my $sn=substr($bsample,0,16);
	my $snshort=substr($bsample,0,12);
	my $c=$baminfo[1];
	my $bam=$baminfo[2];
	my $chrprefix=$baminfo[4];
	$chrprefixhash{$snshort}=$chrprefix;
	$bamhash{$sn}=$bam;
	$bamhash{$snshort}=$bam;
}
close $RNAbamlist;
##########################
#DNA BAM PATH INFORMATION#
##########################
#     1	Cancer_Type	GBM
#     2	Sample_Name	TCGA-02-0003
#     3	Tumor_Barcode	TCGA-02-0003-01A-01D-1490-08
#     4	Tumor_Bampath	/gscmnt/gc8002/info/build_merged_alignments/merged-alignment-blade8-3-3.gsc.wustl.edu-jeldred-17911-a3bb8b59e2e7452ca2887a7afdb43f7d/a3bb8b59e2e7452ca2887a7afdb43f7d.bam
#     5	Normal_Barcode	TCGA-02-0003-10A-01D-1490-08
#     6	Normal_Bampath	/gscmnt/gc8002/info/build_merged_alignments/merged-alignment-blade11-2-14.gsc.wustl.edu-jwang-16703-34303d0838204d9ab205dcdc441175d4/34303d0838204d9ab205dcdc441175d4.bam

####
my %dnatumor;
my %dnanormal;
open(my $DNAbamlist,'<',"/gscmnt/gc2509/dinglab/reyka/TCGA_Status2/WXS_tumor_normal_bams_052016");
while(my $bamline=<$DNAbamlist>){
    chomp $bamline;
    my @baminfo=split(/\t/,$bamline);
	my $sample=$baminfo[1];
    my $tumor=$baminfo[3];
    my $normal=$baminfo[5];
	$dnatumor{$sample}=$tumor;
	$dnanormal{$sample}=$normal;
}
close $DNAbamlist;
#Going through each line of the MAF and printing all lines pertaining to the same sample name into a different file.
while(<$MAF>){
		$counter++;
#		print " $counter";
		chomp ($linest=$_);
		@array1=split(/\t/,$linest);
		my $sample=$array1[15];
		my $longsn=substr($sample,0,16);
		my $shortsn=substr($sample,0,12);
		my $chr=$array1[4];
		my $cancer=$array1[90];	
####IF MAF has long sample name info TCGA-XX-XXXX-XXX-XXXX.... then check for TCGA-XX-XXXX-XXX in bam hash
		if(length($sample)>12){
			if (exists $bamhash{$longsn}){
#Unless the sample and chromosome is already a key in the hash
				unless ($sample_hash{$sample}{$chr}){
#create a file in sampledat/with sample name
					my $file="$sampledir/$sample.$cancer.$chr.txt";
					$sample_hash{$sample}{$chr}=$file;
					open OUT1,'>', $file or die "Unable to open $file for output: $!\n";	
					close OUT1;
				}	
#if the sample already exists as a key in the hash, then print the line out to the file with the sample name.
				my $filenow=$sample_hash{$sample}{$chr};
				open(OUT,'>>',"$filenow") or die "Unable to open $filenow for output: $!\n";
				print OUT "$_";
				close OUT;
			}
		}
####IF MAF has short sample name info TCGA-XX-XXXX then check for TCGA-XX-XXXX in bam hash
		elsif(length($sample) == 12){
		        if (exists $bamhash{$shortsn}){
#Unless the sample and chromosome is already a key in the hash
      	            unless ($sample_hash{$sample}{$chr}){
#create a file in sampledat/with sample name
                    my $file="$sampledir/$sample.$cancer.$chr.txt";
                    $sample_hash{$sample}{$chr}=$file;
                    open OUT1,'>', $file or die "Unable to open $file for output: $!\n";
                    close OUT1;
                }
#if the sample already exists as a key in the hash, then print the line out to the file with the sample name.
                my $filenow=$sample_hash{$sample}{$chr};
                open(OUT,'>>',"$filenow") or die "Unable to open $filenow for output: $!\n";
                print OUT "$_";
                close OUT;
            }
		}
	else{
		print "$sample $cancer has no bam file\n";
		print $NOBAM "$sample $cancer has no bam file\n";
	}
}
close $MAF;
close $NOBAM;
print "Created $counter files in $sampledir/SAMPLE.cancertype.chr.txt\n";
#to run gmt readcounts, we need only the 5 column format for each sample, this next part goes through each sample and makes a new file with the sample name.sites.
#what makes this differ from WXS readcounts is each chromosome needs to be broken up into a separate sites directory
#listing all files in sitesdirectory
print "Making sites files\n";
my $directory="$sampledir";
opendir(DIR,$directory) or die $!;
while (my $file = readdir(DIR)){
	if ($file=~/(.*)\.([a-z]*)\.([0-9XY]*)\.txt/){
#grabbing 5 column format info and saving to new file with same name and ending with .sites.
		my $sitesout="$sitesdir\/$1.$2.$3.sites";
		my $sitesoutdna="$sitesdir\/$1.$2.$3.sitesdna";
#Some cancer samples don't need to include "chr" prefix in the sites file because the bam and reference fasta dont have the "chr" prefix. IF included, the readcount data will always skip the sites because they won't be found in the bam or reference
		my $sampleshort=substr($1,0,12);	
		my $chrprefixvalue=$chrprefixhash{$sampleshort};
		if ($chrprefixvalue=~/SN:chr/){
			system("cut -f 5,6,7,11,13 $sampledir\/$file |awk '{print \"chr\"\$0}' > $sitesout");
			system("cut -f 5,6,7,11,13 $sampledir\/$file > $sitesoutdna");
		}
		else {
			system("cut -f 5,6,7,11,13 $sampledir\/$file > $sitesout");
			system("cut -f 5,6,7,11,13 $sampledir\/$file > $sitesoutdna");
		}
	}
}
closedir(DIR);

#############COMMAND FOR RUNNING READCOUNT ON RNA###########

my $RNAbamlist="/gscmnt/gc2600/dinglab/reyka/Unarchive_Samples/RNA_bampaths_082015";

`mkdir -p readcountdat`;
print "Submitting jobs for readcount analysis, output directory: readcountdat/\n";
my $directory="$sitesdir";

opendir(DIR,$directory) or die $!;
while (my $sitesfile = readdir(DIR)){
    if ($sitesfile=~/(.*)\.([a-z]*)\.([0-9XY]*)\.sites$/){
		my $sample=$1;
		my $samplefullname;
		my $subsample=substr($sample,0,12);
		my $fullsamplename=substr($sample,0,16);
		my $cancer=$2;
		my $chr=$3;
		my $genome;
		my $chrprefixvalue;
		if (exists $chrprefixhash{$subsample}){
			$chrprefixvalue=$chrprefixhash{$subsample};
		}
		if ($chrprefixvalue=~/SN:chr/){
			$genome="/gscmnt/ams1161/info/model_data/kmeltzst/reference_files/b37_fastawhole/chr$chr";
		}
		else{
			$genome="/gscmnt/ams1102/info/model_data/2869585698/build106942997/all_sequences.fa";
		}
		my $bam=$bamhash{$fullsamplename};
		my $log="logs/$sitesfile.log";
		my $logdna="logs/$sitesfile.dna.log";
		my $outfilerc="readcountdat\/$sitesfile.rc";
		my $outfilercdna="readcountdat\/$sitesfile.dnarc";
		my $command=`bsub -g /rjayasinreadcount -oo $log 'gmt analysis coverage add-readcounts --variant-file=$sitesdir\/$sitesfile --bam-files=$bam --genome-build=$genome --indel-size=1 --output-file=$outfilerc'`;
		chomp $command;
		print "RNA\t$command\n";
		my $normal=$dnanormal{$subsample};
		my $tumor=$dnatumor{$subsample};
		my $sitesfiledna=$sitesfile."dna";
		my $commandwxs=`bsub -g /rjayasinreadcount -oo $logdna 'gmt analysis coverage add-readcounts --variant-file=$sitesdir\/$sitesfiledna --bam-files=$normal,$tumor --header-prefixes=\"nrm,trm\" --genome-build=37 --indel-size-limit=1 --output-file=$outfilercdna'`;
		chomp $commandwxs;
		print "DNA\t$commandwxs\n";

	}
}
closedir(DIR);	
#gmt analysis coverage add-readcounts --variant-file=./sitesdat/TCGA-66-2791.prefix_chr.sites --bam-files=/gscmnt/sata864/info/alignment_data/imported/2887771807/all_sequences.bam --genome-build=/gscmnt/ams1161/info/model_data/kmeltzst/reference_files/b37_fastawhole/chr1 --indel-size-limit=1 --output-file=readcount_mmclella
