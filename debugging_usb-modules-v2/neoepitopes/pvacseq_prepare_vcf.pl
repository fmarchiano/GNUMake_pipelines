while (my $line = <>) {
	chomp $line;
	if ($line =~ /^\#\#FORMAT=<ID=DP/) { 
		print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"."\n";
	} 
	if ($line =~ /^\#\#/) { print $line."\n"; next;}

	my @line = split /\t/, $line;
	$line[8] = "GT:".$line[8];
	$line[9] = "0/1:".$line[9];
	print join "\t", @line[0..9]; print "\n";
}




