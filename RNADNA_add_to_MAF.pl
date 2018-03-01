use strict;

my $usage =<<USAGE;
 Usage: $0 <MAF FILE> <OUTPUT SITES DIR> <OUTPUT SAMPLES DIR>
  Example: perl RNA_add_to_MAF.pl MAF sitesdat samplesdat

USAGE
    die $usage unless @ARGV==3;
my (@array1,$sample,$file,$linest);
open(my $MAF,'<',$ARGV[0]) or die "INPUT MAF not found!";
my $fileout="$ARGV[0]\.rc";
open(my $OUT,'>',$fileout) or die "INPUT MAF not found!";
open(my $ERROUT,'>',"NoRNAdata.maf") or die "INPUT MAF not found!";
open(my $ERROUTDNA,'>',"NoDNAdata.maf") or die "INPUT MAF not found!";
my %sample_hash;
my $counter=0;
my $sampledir=$ARGV[2];
my $sitesdir=$ARGV[1];
#Going through each line of MAF and searching for RNA Readcount results
while(<$MAF>){
	$counter++;
	chomp ($linest=$_);
	@array1=split(/\t/,$linest);
	$sample=$array1[15];
	my $chr=$array1[4];
	my $cancer=$array1[90];
	my $mafstart=$array1[5];
	my $sitesfile="readcountdat/$sample\.$cancer\.$chr\.sites.rc";
	my $sitesfiledna="readcountdat/$sample\.$cancer\.$chr\.sites.dnarc";
	if (-e $sitesfile){
	if (-e $sitesfiledna){
		open (my $RNA,'<',$sitesfile) or die "Can't open sites file!";
		while(my $line=<$RNA>){
			chomp $line;
			my @rnareadcount=split(/\t/,$line);
			my $start=$rnareadcount[1];
			if ($mafstart==$start){
				shift @rnareadcount;
				shift @rnareadcount;
				shift @rnareadcount;
				shift @rnareadcount;
				shift @rnareadcount;	
				print $OUT "$linest\t$rnareadcount[0]\t$rnareadcount[1]\t$rnareadcount[2]\t";
			}
		}
###DNA RC DATA
#     1	7
#     2	151874731
#     3	151874731
#     4	T
#     5	A
#     6	222
#     7	0
#     8	0.00
#     9	127
#    10	85
#    11	40.09
		open (my $DNA,'<',$sitesfiledna) or die "Can't open sites file!";
		while(my $linedna=<$DNA>){
			chomp $linedna;
			my @dnareadcount=split(/\t/,$linedna);
			my $dnastart=$dnareadcount[1];
			if ($mafstart==$dnastart){
				shift @dnareadcount;
				shift @dnareadcount;
				shift @dnareadcount;
				shift @dnareadcount;
				shift @dnareadcount;
				print $OUT "$dnareadcount[0]\t$dnareadcount[1]\t$dnareadcount[2]\t$dnareadcount[3]\t$dnareadcount[4]\t$dnareadcount[5]\n";
			}
			#close $RNA;
			#close $DNA;
		}
	}
	else{
		print $ERROUTDNA "$linest\n";
	}
	}
	else{
		print $ERROUT "$linest\n";
	}
}
close $MAF;
close $OUT;
close $ERROUT;
