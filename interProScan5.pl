#!/usr/bin/perl
use strict;
my $usage = <<USAGE;
Usage:
    perl $0 seq.fasta email_adress outDir threads

USAGE
unless (@ARGV==4){die $usage}

my ($seqFasta, $email, $outDir, $threads) = @ARGV;
mkdir $outDir unless -e $outDir;

open IN, $seqFasta or die $!;
my ($seqID,@seqID,%seq);
while (<IN>) {
	chomp;
	if (/^>(\S+)/) { push @seqID, $1; $seqID = $1 }
	else          { $seq{$seqID} .= $_ }
}

open COM, '>', "$outDir.command" or die $!;

foreach (@seqID) {
	open FASTA, '>', "$outDir/$_.fa" or die $!;
	print FASTA ">$_\n$seq{$_}";
	close FASTA;

	print COM "~/bin/iprscan5_lwp.pl --goterms --email $email --title $_ --jobid $_ --outfile $outDir/$_ $outDir/$_.fa\n";
}

system `ParaFly -c $outDir.command -CPU $threads`;
