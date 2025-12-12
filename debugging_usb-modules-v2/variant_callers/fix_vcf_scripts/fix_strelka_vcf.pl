use strict;

my $normal_first = 0;

while (my $line = <>) {
	chomp $line;
	if ($line =~ /^\#/) { 
		if ($line =~ /^(\#CHROM.+)NORMAL\tTUMOR/) { $line = "$1"."TUMOR\tNORMAL"; $normal_first = 1; }
		print $line; print "\n"; 
		if ($line =~ /^\#\#FORMAT\=\<ID\=TOR/) { 
			print "\#\#FORMAT\=\<ID\=AD,Number=2,Type=Float,Description=\"AD computed from tier 1\"\>\n";
			print "\#\#FORMAT\=\<ID\=FA,Number=2,Type=Float,Description=\"AF computed from tier 1\"\>\n";
		}
	} else {
		my @arr = split /\t/, $line;
		my @format = split /:/, $arr[8];
		my $n_dp = my $n_tir = my $n_alt = 0;
		for (my $i=0; $i < scalar @format; $i++) {
			if ($format[$i] eq "DP") {
				$n_dp = $i; 
			} elsif ($format[$i] eq "TIR") {
				$n_tir = $i;
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
		$normal[$n_tir] =~ /^(\d+),/;
		my $normal_ad = ($normal[$n_dp]-$1).",".$1;
		my $normal_af = "";
		if ($normal[$n_dp]!=0) { $normal_af = $1/$normal[$n_dp];}

		$tumour[$n_tir] =~ /^(\d+),/;
		my $tumour_ad = ($tumour[$n_dp]-$1).",".$1;
		my $tumour_af = "";
		if ($tumour[$n_dp]!=0) { $tumour_af = $1/$tumour[$n_dp];}

		push @format, "AD", "FA";
		push @normal, $normal_ad, $normal_af;
		push @tumour, $tumour_ad, $tumour_af;

		$arr[8] = join ":", @format;
		$arr[9] = join ":", @tumour;
		$arr[10] = join ":", @normal;

		print join "\t", @arr; print "\n";
	}
}
		
	
