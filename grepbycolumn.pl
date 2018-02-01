#!/bin/perl
#30 July 2014 - Adam D Scott - updated 13 October 2014

use strict;
use warnings;

my $usage = 'perl grepbycolumn.pl file_to_read word column match( 2=exact, 1=contains, 0=mismatch)
';
die $usage , unless @ARGV == 4;

my ( $file_to_read , $word , $col , $match ) = @ARGV;

open ( file_to_read , "<".$file_to_read ) or die "Could not open $file_to_read $!";
while ( <file_to_read> )
{
	chomp;
	my @line = split( "\t" , $_ );
	if ( $match > 0 )
	{#match exact
		if ( $match == 1 )
		{
			if ( $line[$col-1] =~ m/$word/ )
			{
				print $_."\n";
			}
		}
		else
		{
			if ( $line[$col-1] eq $word )
			{
				print $_."\n";
			}
		}
	}
	else
	{#not match
		if ( $line[$col-1] !~ m/$word/ )
		{
			print $_."\n";
		}
	}
}
close file_to_read;
