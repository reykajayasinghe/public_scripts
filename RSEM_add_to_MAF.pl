use strict;

my $usage =<<USAGE;
 Usage: $0 <MAF FILE> <RSEM FILE>
  Example: perl RSEM_add_to_MAF.pl MAF RSEMFILE

RSEMFILE:
brca	TCGA-A1-A0SG-01A-11D-A142-09	MRPL46	15_89008087_G_A	228.6023	
brca	TCGA-A1-A0SI-01A-11D-A142-09	BMPR2	2_203329576_C_T	3019.4232	
brca	TCGA-A1-A0SI-01A-11D-A142-09	CHURC1	14_65392796_C_T	1408.9935	
MAF:
Need to have the following in the right columns:
	cancer type in column 91
	sample in column 16
	gene in column 1
	chr start in column 5
	start position in column 6
	ref
	var

USAGE
    die $usage unless @ARGV==2;
my (@array1,$sample,$file,$linest);
open(my $MAF,'<',$ARGV[0]) or die "INPUT MAF not found!";
my $fileout="$ARGV[0]\.rsem";
open(my $OUT,'>',$fileout) or die "INPUT MAF not found!";
open(my $ERROUT,'>',"NoRSEM.maf") or die "INPUT MAF not found!";

my %rsemdata;
open(my $RSEM,'<',$ARGV[1]) or die "INPUT RSEM DATA not found!";
while(<$RSEM>){
    chomp(my $line=$_);
    my @line=split(/\t/,$line);
    my $cancer=shift @line;
    my $sample=shift @line;
    my $gene=shift @line;
	my $key=shift @line;
	my $rsem=shift @line;
	my $subsample=substr($sample,0,12);
    $rsemdata{$key}{$subsample}=$rsem;
}
close $RSEM;

my $counter=0;
my $rsemdata=0;
my $norsemdata=0;
#Going through each line of MAF and searching for RNA Readcount results
while(<$MAF>){
	$counter++;
	chomp ($linest=$_);
	@array1=split(/\t/,$linest);
	my $s=$array1[15];
	my $sn=substr($s,0,12);
	my $chr=$array1[4];
	my $mafstart=$array1[5];
	my $ref=$array1[10];
	my $var=$array1[12];
	my $mafkey=$array1[4]."_".$mafstart."_".$ref."_".$var;
	if (exists $rsemdata{$mafkey}{$sn}){
		$rsemdata++;
		my $value=$rsemdata{$mafkey}{$sn};
		print $OUT "$linest\t$value\n";
	}
	else{
		$norsemdata++;
		print $ERROUT "$linest\n";
	}
}
close $MAF;
close $OUT;
close $ERROUT;

print "Created $fileout\n";
print "SUMMARY:\n";
print "Total mutations: $counter\n";
print "Total no RSEM  : $norsemdata\n";
print "Total with RSEM: $rsemdata\n";

