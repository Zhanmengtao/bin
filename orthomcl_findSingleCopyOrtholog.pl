#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 groups.txt compliantFasta

USAGE
die $usage if @ARGV != 2;

open GROUPFILE, $ARGV[0] or die $!;
open OUTPUT, '>', "SingleCopyOrthologGroups.txt" or die $!;
mkdir "SingleCopyOrthologGroups" unless -e "SingleCopyOrthologGroups";

my (%fasta, @species);
foreach (<$ARGV[1]/*.fasta>) {
	chomp;
	my $species;
	if (/(\w+).fasta/) { $species = $1; push @species, $species;}

	print STDERR "reading $_\n";
	open IN, $_ or die $!;
	my $protein_id;
	while (<IN>) {
		chomp;
		if (/>(\S+)/) { $protein_id = $1; }
		else { $fasta{$species}{$protein_id} .= $_; }
	}
	close IN;
}

my $species_num = @species;
print STDERR "total $species_num species\n@species\n\n";

my $groupNum;

while (<GROUPFILE>) {
	my $group = $1 if s/^([^:]+?):\s+?//;
	s/\s*$//;
	my @sequences = split /\s+/,$_;
	next if @sequences != $species_num;

	my %species;
	foreach (@sequences) { $species{$1} = $_ if /([^\s\|]+)\|/; }
	next unless keys %species == $species_num;

	$groupNum ++;
	print OUTPUT  "$group: @sequences\n";

	open OUT, '>', "SingleCopyOrthologGroups/$group.fasta" or die $!;

	foreach (@species) {
		my $sequenceId = $species{$_};
		print OUT ">$sequenceId\n$fasta{$_}{$sequenceId}\n";
	}
}

print STDERR "$groupNum groups of Single Copy Orthologs were found!\n";
print STDERR "ALL the files were wrriten to SingleCopyOrthologGroups.txt and directory SingleCopyOrthologGroups\n";
