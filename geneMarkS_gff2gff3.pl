#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
	perl $0 genome.fasta.gff Locus_tag > genome.fasta.gff3

USAGE
if (@ARGV==0){die $usage}

my $locus_tag = pop @ARGV;

open IN, $ARGV[0] or die $!;
while (<IN>) {
	if (/^#/ or /^\s*$/) { next; }
	elsif (/\tCDS\t/) {
		my @str = split /\t/, $_;
		my $code = $1 if $str[8] =~ m/gene_id=(\d+)/;
		if (length $code == 1) { $code = $locus_tag . "_00000" . $code; }
		if (length $code == 2) { $code = $locus_tag . "_0000" . $code; }
		if (length $code == 3) { $code = $locus_tag . "_000" . $code; }
		if (length $code == 4) { $code = $locus_tag . "_00" . $code; }
		if (length $code == 5) { $code = $locus_tag . "_0" . $code; }

		print "$str[0]\t$str[1]\tgene\t$str[3]\t$str[4]\t\.\t$str[6]\t\.\tID=$code;Name=$code;\n";
		print "$str[0]\t$str[1]\tmRNA\t$str[3]\t$str[4]\t\.\t$str[6]\t\.\tID=$code.mRNA;Name=$code;Parent=$code;\n";
		print "$str[0]\t$str[1]\tCDS\t$str[3]\t$str[4]\t$str[5]\t$str[6]\t$str[7]\tID=$code.cds;Parent=$code.mRNA;\n";
	}
}
close IN;


