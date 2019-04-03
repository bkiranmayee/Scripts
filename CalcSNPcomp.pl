#!/usr/bin/perl
# CalcSNPcomp.pl 
# Author: Kiranmayee Bakshy
# Date: Jul 24 2018 
#A program to calculate the marker composition data like % GC content, no. of IUPAC code Nucleotides and their number...

# the input file is a tab delimited file with the marker ID, a sequence of 150 bps before  the SNP,
# SNP in the form of [A/B], a sequence of 150 bps after the SNP
 
 # the output file consists of tab delimited file with marker ID, % GC before SNP, %GC after SNP, no. of IUPAC codes before and after the SNP,
 # no. of nucleotides coded as 'N' before and after the SNP
 
 use strict; 
 use warnings; 
 die "usage: CalcSNPcomp.pl <tab delimited ID, seq file>\n" unless @ARGV == 1;
 my $input = "snp_seq";
open (my $OUTFILE, "> snp.seq.comp");

open(my $INFILE, "<$input") || die "Could not open the input file:$input!\n";


while(my $line = <$INFILE>){
	chomp $line;
	my @segs = split(/\t/, $line);
	my @newsegs;
	foreach my $seq(@segs[1,3]) { 
			my $len=length($seq);
			my $GC_content;
			my @bases = split(//,$seq);
			my $numG=0;
			my $numC=0;
			my $numT=0;
			my $numA=0;
			my $numN=0;
			my $numS=0;
			my $numIUPAC=0;
			
				foreach my $bp(@bases) {
					if($bp =~ m/G/i){$numG++};
				    if($bp =~ m/C/i){$numC++};
				   if($bp =~ m/N/i){$numN++};
				   if($bp =~ m/S/i){$numS++};
				   if($bp =~ m/[RYWKMBDHV]/i){$numIUPAC++};
				}
				
				$GC_content =   (($numG+$numC+$numS)/$len)*100;
				push @newsegs, $GC_content,$numN,$numIUPAC;
				
	}
print {$OUTFILE} join("\t", $segs[0], $segs[2],@newsegs) . "\n";
		
}
			

close $INFILE;
close $OUTFILE;
exit;