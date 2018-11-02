#!/usr/bin/perl
use strict;
my $usage = <<USAGE;
Usage:
	perl $0 SingleCopyOrthologGroups 24[this is CPU threads]

USAGE
die $usage if @ARGV != 2;

my @files = `ls $ARGV[0]`;

open CMDS, '>', "mafft_cmds" or die $!;

foreach (@files) {
	chomp;
	next unless m/\.fasta$/;
	print CMDS "linsi $ARGV[0]/$_ > $ARGV[0]/$_.align\n";
}

close CMDS;

my $parafly_result =`ParaFly -c mafft_cmds -CPU $ARGV[1] -shuffle -v`;

#unlink "clustalo_cmds", "clustalo_cmds.completed" unless -e "FailedCommands";

my @files = `ls $ARGV[0]`;

my %sequences;
my $sequenceId;
my $sequenceLength;

foreach (@files) {
	chomp;
	next unless m/\.align$/;

	open INPUT, "$ARGV[0]/$_" or die $!;

	while (<INPUT>) {
		chomp;
		if (/>(.*?)\|/) { $sequenceId = $1 }
		else         { $sequences{$sequenceId} .= $_; $sequenceLength += length $_ }
	}
	close INPUT;
}

my @sequences = keys %sequences;
@sequences = sort { $a cmp $b } @sequences;

open RESULT, '>', 'allSingleCopyOrthologsAlign.fasta' or die $!;

foreach (@sequences) {
	print RESULT ">$_\n$sequences{$_}\n";
}

print "\nAligned length is $sequenceLength\n";
