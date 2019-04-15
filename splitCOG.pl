#!/usr/bin/perl
# This is a script to split grouped COG category counts to individual category counts 

use strict;
use warnings;
use Data::Dumper qw(Dumper);

if (scalar @ARGV != 2 ) {
    print STDERR "Usage: splitCOG.pl <tab-delimited COG counts file> <output filename>\n";
    exit(1);
}

my ($input, $output) = @ARGV;

#my @main_cat = ("D","M","N","O","T","U","V","W","Y","Z","A","B","J","K","L","C","E","F","G","H","I","P","Q","R","S");
open(my $IN, "< $input") || die "Could not open input file: $input!\n";
open(my $OUTFILE, "> $output" );	
print {$OUTFILE} "COG_category\tIA\tPA\tIB\tPB\tIE\tPE\n";

my @cogIdx;
my %cogCls;

my $head = <$IN>;
chomp $head;
my @hsegs = split(/\t/, $head);
shift @hsegs;
@cogIdx = @hsegs;
#print Dumper \@cogIdx;

while(my $line = <$IN>){
	chomp $line;
	my @cols = split(/\t/, $line);
	my $cog = shift @cols;
	for(my $x = 0; $x < scalar(@cols); $x++){
			my $cls = $cogIdx[$x];
			$cogCls{$cog}->{$cls} = $cols[$x];
	}			
}
#print Dumper \%cogCls;

my %newcogCls = map { $_ => $cogCls{$_} } grep { /^[A-Z]$/ } keys %cogCls;
#print Dumper \%newcogCls;

foreach my $cat (sort keys %cogCls)  {
    if ($cat =~ /.*\,.*/) {
		my @subcat = split(/, /, $cat);
		foreach my $i(@subcat){
			foreach my $cls(@cogIdx){
			if ($cogCls{$cat}->{$cls} ne "NA"){
				$newcogCls{$i}->{$cls} += sprintf('%.0f', ($cogCls{$cat}->{$cls})/scalar(@subcat));
			} else {next;
		 }
		}
	}
	} else {next;} 
}
#print Dumper \%cogCls;
#print Dumper \%newcogCls;



foreach my $cat (sort keys %newcogCls) {
	my @val;
	foreach my $cls (@cogIdx) {
	push (@val, $newcogCls{$cat}{$cls});
	}
	printf {$OUTFILE} join("\t", $cat, @val) . "\n";	
}
print "done\n";
close $IN;
close $OUTFILE;
exit;