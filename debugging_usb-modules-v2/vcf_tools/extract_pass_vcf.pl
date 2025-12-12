#!/usr/bin/perl

use strict;
use Getopt::Std;

my %opt;
getopts('n:', \%opt);

my $usage = <<ENDL;
Usage: extract_pass_vcf.pl [-n <max_allowed>] <input.vcf> <filtered.vcf> <filtered_out.vcf>

-n : the maximum number of filters (excluding targetInterval) allowed to hotspots to be whitelisted
ENDL

my $max_filters;
if ($opt{n} eq "") { $max_filters = 2;
} else { $max_filters = $opt{n}; }

if (@ARGV!=3) {
	print STDERR $usage; 
}

open IN, "$ARGV[0]";
open OUT, ">$ARGV[1]";
open EXCLUDE, ">$ARGV[2]";

while (my $line = <IN>) {
	chomp $line;
	#print header lines, and add the new hotspotPASS FILTER key.
	if ($line =~ /^\#\#fileformat\=VCF/) { 
		print OUT $line."\n##FILTER=<ID=hotspotPASS,Description=\"rescued hotspot\">\n"; print EXCLUDE $line."\n##FILTER=<ID=hotspotPASS,Description=\"rescued hotspot\">\n";
	} elsif ($line =~ /^\#/) { 
		print OUT $line."\n"; print EXCLUDE $line."\n";
	}
	else {
		my @line = split /\t/, $line;

		# print anything that has not been filtered out
		# MuSE can have PASS and PASS;Tier[1-5] format
		if ($line[6] eq "." || $line[6] eq "PASS" || $line[6] =~ /^PASS;Tier\d$/ ) {
			print OUT $line."\n";

		# print nothing if the mutation is outside of the target intervals
		} elsif ($line[6] =~ /targetInterval/) { 

		# this part is where the whitelisting of the hotspots happens
		# hotspots are annotated in the INFO field
		} else {
			my $info = $line[7];
			my @filters = split /;/, $line[6];
			my %filters = ();
			foreach my $filter (@filters) {
				$filters{$filter} = 1;
			}

			if ($info =~ /HOTSPOTaa/ || $info =~ /HOTSPOT3Daa/ || $info =~ /HOTSPOTindel/ || $info =~ /HOTSPOTsp/ || $info =~ /HOTSPOTNC/) {
				if (scalar keys %filters <= $max_filters ) { 
					print OUT join("\t", @line[0..6]).";hotspotPASS\t".join("\t", @line[7..$#line])."\n"; # Add hotspotPASS when rescued by hotspot.
				} else { print EXCLUDE $line."\n"; }
			} else { print EXCLUDE $line."\n"; }
		}
	}
}
close IN;
close OUT;
close EXCLUDE;
