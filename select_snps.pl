#!/usr/bin/perl
# A program to a) parse dbsnp vcf file b) and select variants whose REF allele matches with the new reference genome assembly 

use strict;
use Data::Dumper;


# make sure the user passed the required arguments
if (scalar @ARGV != 2 ) {
   print STDERR "Usage: select_snps.pl <input vcf file> <output file>\n";
 	  exit(1);
}
chomp(@ARGV);
my ($vcffile, $outfile) = @ARGV;

open(my $ifh, "gunzip -c $vcffile |") || die "ERROR: failed to read vcf file: $!";
open (my $OUTFILE, "> $outfile");

while (my $line = <$ifh>) {
	chomp $line;
	if ($line =~ '#') {
	next;
	} else {
	my $ref =~ /\;Reference\_seq=([ACGT])/;
	my @cols = split (/\t/, $line);
	if ($cols[3] eq $ref) {
		print $OUTFILE "$line\n";
				last;
		      	}
			print $OUTFILE "$line\n";
}

close $ifh;
close $OUTFILE;

exit;	
	