use strict;

my $normal_first = 0;

# add FA from AD

while (my $line = <>) {
	chomp $line;
	if ($line =~ /^\#/) { 
		if ($line =~ /^(\#CHROM.+)NORMAL\tTUMOR/) { $line = "$1"."TUMOR\tNORMAL"; $normal_first = 1; }
		print $line; print "\n"; 
		if ($line =~ /^\#\#FORMAT\=\<ID\=SS/) { 
			print "\#\#FORMAT\=\<ID\=FA,Number=2,Type=Float,Description=\"Allele frequency\"\>\n";
		}
	} else {
		my @arr = split /\t/, $line;
		my @format = split /:/, $arr[8];
		my $n_ad = 0;
		for (my $i=0; $i < scalar @format; $i++) {
			if ($format[$i] eq "AD") {
				$n_ad = $i; 
			}
		}
		my @normal = my @tumour = "";
		if ($normal_first == 1) {
			@normal = split /:/, $arr[9];
			@tumour = split /:/, $arr[10];
		} else {
			@tumour = split /:/, $arr[9];
			@normal = split /:/, $arr[10];
		}
		$normal[$n_ad] =~ /^(\d+),(\d+)/;
		my $normal_af = "";
		if ($1 + $2 !=0) { $normal_af = $2/($1 + $2);}

		$tumour[$n_ad] =~ /^(\d+),(\d+)/;
		my $tumour_af = "";
		if ($1 + $2 !=0) { $tumour_af = $2/($1 + $2);}

		push @format, "FA";
		push @normal, $normal_af;
		push @tumour, $tumour_af;

		$arr[8] = join ":", @format;
		$arr[9] = join ":", @tumour;
		$arr[10] = join ":", @normal;

		print join "\t", @arr; print "\n";
	}
}
