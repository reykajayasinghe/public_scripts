use strict;

my $usage =<<USAGE;
 Usage: $0 <MAF FILE> <Copy Number file>
  Example: perl /gscmnt/gc2509/dinglab/reyka/splice_project/scripts/CN_add_to_MAF.pl MAF CN_FILE
	CN_FILE: Generated using parseCNV.pl

USAGE
    die $usage unless @ARGV==2;
########################
#PARSE COPY NUMBER DATA#
########################
###COPY NUMBER DATA INPUT###
#header test_run/10_89692907_89692907_TCGA-02-0055-01A-01D-1490-08_GBM_GBM.CN.data.final
#     1	10
#     2	89692907
#     3	89692907
#     4	TCGA-02-0055-01A-01D-1490-08
#     5	GBM
#     6	10
#     7	6233955
#     8	135225087
#     9	77016
#    10	-0.2949
#    11	TCGA-02-0055-01A-01D-0182-01
#    12	1.63025766173921
my %copynumberdata;
open(my $CN,'<',$ARGV[1]) or die "INPUT COPY NUMBER DATA not found!";
while(<$CN>){
	chomp(my $line=$_);
	my @cnline=split(/\t/,$line);
	my $chr=shift @cnline;
	if ($chr=='23'){
		$chr='X';
	}
	if ($chr=='24'){
		$chr='Y';
	}
	my $start=shift @cnline;
	my $stop=shift @cnline;
	my $sample=shift @cnline;
	my $cancer=shift @cnline;
	my $key=$chr.$start.$stop.$sample.$cancer;
	$copynumberdata{$key}=[@cnline];
#	print "@cnline\n";
}
close $CN;
################
#PARSE MAF FILE#
################

my (@array1,$sample,$file,$linest);
open(my $MAF,'<',$ARGV[0]) or die "INPUT MAF not found!";
my $fileout="$ARGV[0]\.cn";
open(my $OUT,'>',$fileout) or die "INPUT MAF not found!";
open(my $ERROUT,'>',"NoCNdata.maf") or die "INPUT MAF not found!";
my %sample_hash;
my $counter=0;
#Going through each line of MAF and searching for Copy Number Match
while(<$MAF>){
	$counter++;
	chomp ($linest=$_);
	@array1=split(/\t/,$linest);
	$sample=$array1[15];
	my $chr=$array1[4];
	my $cancer=uc($array1[90]);
	my $mafstart=$array1[5];
	my $mafstop=$array1[6];
	my $mafkey=$chr.$mafstart.$mafstop.$sample.$cancer;
	if (exists $copynumberdata{$mafkey}){
		#my @valarray=@{$copynumberdata{$mafkey}};
		my @valarray=join("\t",@{$copynumberdata{$mafkey}});
		print $OUT "$linest\t@valarray\n";
	}
	else{
		print $ERROUT "$linest\n";
		print STDERR "No copy number data:$chr $mafstart $sample $cancer\n";
	}
}
close $MAF;
close $OUT;
close $ERROUT;
