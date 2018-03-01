#!/gsc/bin/perl
#Mike Wendl's Code, Edited by Reyka Jayasinghe
#Last Edit: May 26th, 2016
#  CLASSIFIES SAMPLES BASED ON EXPRESSION LEVELS
#
#	INPUT:
#	coadread	TCGA-F5-6814-01A-31D-1924-10	ABCB1	7_87179792_C_A	205.6108	404.6906,216.3174,561.3318,40.8068,2363.6194,335.6037,312.1721,68.2505,874.0691,129.5823,19.2308,939.3673,1240.0636,2007.5099,996.8909,1141.9395,948.5078,2262.4693,1746.6226,20.5479,189.1544,5119.7425,234.6227,740.5235,3561.7767,109.0544,94.7752,3780.3030,282.3680,55.7851,2355.6531,1545.8216,838.8498,632.2155,174.0561,111.5257,73.3399,1545.7237,650.0903,1768.6901,286.1089,956.4396,2625.7615,778.3489,914.7856,2745.1172,1885.7843,1116.9291,2460.8968,1830.0626,1497.0060,312.0777,768.4664,987.8424,44.2319,1347.4692,979.9480,316.2502,653.7173,1187.1369,89.5126,158.7302,1900.2532,89.9941,2153.7611,501.5496,356.8399,652.6718,232.8691,325.5379,4024.7844,971.9840,408.8165,195.9799,1093.3382,154.2305,1113.1426,870.7079,18.7154,1262.7045,3745.9195,64.8869,39.0379,34.8177,1544.0759,251.6634,177.4744,1487.8496,2156.7837,939.9915,300.5299,2878.5607,398.5956,2116.9458,716.5775,1252.8965,16.7826,339.1207,426.9360,2222.5000,339.6226,5090.1792,922.5777,3334.7763,1210.8621,964.3520,475.2483,183.2341,2021.7856,456.9872,431.1024,2989.5833,628.7466,16.3351,142.4443,1067.6888,1632.6096,1555.9591
#
#	THIS SCRIPT CAN BE USED TO DETERMINE WHETHER OR NOT A SAMPLE IS EXPRESSED
#	COMPARES THE CASE GENE EXPRESSION VALUE TO CONTROL DISTRIBUTION
#	IF CASE IS GREATER THAN LOWER QUARTILE OF CONTROL DISTRIBUTION:
# 	OUTPUT: 
#	cancer	gene	sample	pvalue	case	lowerquartile	totalcontrols
#	kirp	PPM1D	TCGA-A4-7585-01A-11D-2136-08	0.0104166666666667	961.0129	279.0197	286
#	cesc	PPM1D	TCGA-EK-A3GJ-01A-21D-A20U-09	0.314741035856574	400.9756	251.760175	249
#	hnsc	PPM1D	TCGA-D6-6516-01A-11D-1870-08	0.16632860040568	338.8463	187.882375	49
#  NOTES
#
#  * FDR
#
#    We've used FDR many times and code is widely available, see for example
#
#    research/completed/2011_pathscan/tsp_analysis/process_pathscan_file.pl
#    research/in_progress/blood_clonal_expansion/statistics_blood_normal_test/complete_statistical_pipeline_matched.pl
#
############
#  SET UP  #
############

#__STANDARD PACKAGES
   use strict;
   use warnings;

#__CPAN PACKAGES
   use Statistics::Descriptive;
   use Statistics::RankCorrelation;

#__SPECIAL (NON-STANDARD) TGI PACKAGES
#  use lib "../integrated_analysis/scripts";
 use lib "./";
 use SpliceInator_Quartile ':all';
 use lib "/gscuser/mwendl/submissions/in_progress/splice_site_reyka";  
# use lib "/gscmnt/gc2706/dinglab/medseq/LabCode/Reyka/SpliceInator/SpliceInator"
 use Statistics::CombinePvals;

###################
#  CONFIGURATION  #
###################

#__THRESHOLD P-VALUE CUTOFF FOR ALL GROUPS OF CLASSIFICATION #EXCEPT GENE EXPRESSION - HARD CODED INTO GENE EXPRESSION TEST
   my $threshold_classification = 0.05;

#__INPUT: FULL EXON RPKM EXPRESSION OF CASES AND CONTROLS
   my $path = ".";
##################
##################
##################
my $whole_gene_rsem_file=$ARGV[0];
########## WHOLE GENE RSEM FILE #############
# CANCERTYPE SAMPLE GENE MUTATION CASERSEM CONTROLRSEM(COMMA SEPARATED)
#coadread	TCGA-F5-6814-01A-31D-1924-10	ABCB1	7_87179792_C_A	205.6108	404.6906,216.3174,561.3318,40.8068,2363.6194,335.6037,312.1721,68.2505,874.0691,129.5823,19.2308,939.3673,1240.0636,2007.5099,996.8909,1141.9395,948.5078,2262.4693,1746.6226,20.5479,189.1544,5119.7425,234.6227,740.5235,3561.7767,109.0544,94.7752,3780.3030,282.3680,55.7851,2355.6531,1545.8216,838.8498,632.2155,174.0561,111.5257,73.3399,1545.7237,650.0903,1768.6901,286.1089,956.4396,2625.7615,778.3489,914.7856,2745.1172,1885.7843,1116.9291,2460.8968,1830.0626,1497.0060,312.0777,768.4664,987.8424,44.2319,1347.4692,979.9480,316.2502,653.7173,1187.1369,89.5126,158.7302,1900.2532,89.9941,2153.7611,501.5496,356.8399,652.6718,232.8691,325.5379,4024.7844,971.9840,408.8165,195.9799,1093.3382,154.2305,1113.1426,870.7079,18.7154,1262.7045,3745.9195,64.8869,39.0379,34.8177,1544.0759,251.6634,177.4744,1487.8496,2156.7837,939.9915,300.5299,2878.5607,398.5956,2116.9458,716.5775,1252.8965,16.7826,339.1207,426.9360,2222.5000,339.6226,5090.1792,922.5777,3334.7763,1210.8621,964.3520,475.2483,183.2341,2021.7856,456.9872,431.1024,2989.5833,628.7466,16.3351,142.4443,1067.6888,1632.6096,1555.9591

my $group_i_hits = {};
#__SET MINIMUM NUMBER OF ELEMENTS TO ACTUALLY DO PERMUTATION TEST
   my $min_elements_controls = 50;
   my $min_elements_samples = 10;

####################
#  PRE PROCESSING  #
####################

# hnsc      LRP1B    TCGA-CN-5366-01A-01D-1434-08       NO

   $min_elements_samples = 3 unless $min_elements_samples > 3;

#__READ FULL EXPRESSION FILES
   my $cases = read_full_expression (
      $whole_gene_rsem_file, $group_i_hits
   );


my $controls;

#####################
#  MAIN PROCESSING  #
#####################

#__DATA STRUCTURE FOR THE CATALOG OF P-VALUE RESULTS
   my $catalog_pvals = {};

#######################################################
#  GROUP I: NO OR LOW-LEVEL EXPRESSION OF WHOLE GENE  #
#######################################################

#__DIAGNOSTIC VALUES FOR REYKA --- EXPRESSION LEVELS OF CASE FOR EACH CATEGORY
   my $diagnostics = {}; 

#print "# <GROUP I> NO OR LOW-LEVEL EXPRESSION OF WHOLE GENE\n";
my ($rsem_data);
($rsem_data) = parse_gene_rsem($whole_gene_rsem_file);
#my ($cancertype,$genes,$samples,$caseval,$pval);
#($cancertype,$genes,$samples,$caseval,$pval) = gene_expression($rsem_data);
my $group_i_count;
   ($group_i_hits, $group_i_count,$diagnostics) = gene_expression ($cases,$controls,$catalog_pvals,$threshold_classification,$diagnostics,$rsem_data);
###Bottom part initially in there?  
 $group_i_count = 0;
   foreach my $cancer (keys %{$group_i_hits}) {
      foreach my $gene (keys %{$group_i_hits->{$cancer}}) {
         foreach my $sample (keys %{$group_i_hits->{$cancer}->{$gene}}) {
	    #__RETRIEVE AND OUTPUT
            my ($pval, $test) = @{$group_i_hits->{$cancer}->{$gene}->{$sample}};
            print "$cancer   $gene   $sample   $pval ($test)\n";
           # my ($sample, $gene) = @{$pair};
            $group_i_count++;
           #print "$cancer   $gene   $sample\n";
		 }
      }
   }


