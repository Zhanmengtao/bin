#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
	perl $0 groups.txt compliantFasta > genes_number_stats.txt

USAGE
if (@ARGV==0) {die $usage}

open IN, $ARGV[0] or die $!;
my (%group_species_gene, @groupID, %speciesID, %all_mcl_genes);
while (<IN>) {
	chomp;
	my $groupID = $1 if s/^([^:]+):\s+//;
	push @groupID, $groupID;
	my @genes = split /\s+/, $_;

	foreach (@genes) {
		if (/(\S+)\|(\S+)/) {
			$speciesID{$1} = 1;
			$group_species_gene{$groupID}{$1}{$2} = 1;
			$all_mcl_genes{$2} = 1;
		}
	}
}
close IN;

unless (-e $ARGV[1]) { die "file $ARGV[1] not exists! $!" }

my @files = `ls $ARGV[1]`;
my %uniqGenes;
foreach my $file (@files) {
	chomp($file);
	if ($file =~ /^(\w+)\.fasta$/) {
		my $species = $1;
		open IN, "$ARGV[1]/$file" or die $!;
		while (<IN>) {
			if (/^>\S+\|(\S+)/) {
				$uniqGenes{$species}{$1} = 1 unless exists $all_mcl_genes{$1};
			}
		}
	}
}

my @speciesID = keys %speciesID;
my $species_number = @speciesID;

foreach my $species (@speciesID) {
	print "$species\n";

	my @genes = keys %{$uniqGenes{$species}};
	foreach (@groupID) {
		if (exists $group_species_gene{$_}{$species}) {
			my @species = keys %{$group_species_gene{$_}};
			next if @species != 1;
			push @genes, keys %{$group_species_gene{$_}{$species}};
		}
	}

	my $genes_num = @genes;

	print "1\t$genes_num\t";
	my $genes_out;
	foreach (@genes) { $genes_out .= "$_ "; }
	$genes_out =~ s/ $//;
	print "$genes_out\n";

	foreach my $genomes_number (2..$species_number) {
		print "$genomes_number\t";
		my @genes = ();
		my $genes_num = 0;
		my $genes_out = "";
		foreach (@groupID) {
			if (exists $group_species_gene{$_}{$species}) {
				my @species = keys %{$group_species_gene{$_}};
				next if @species != $genomes_number;
				push @genes, keys %{$group_species_gene{$_}{$species}};
			}
		}
		$genes_num = @genes;
		print "$genes_num\t";
		foreach (@genes) { $genes_out .= "$_ "; }
		$genes_out =~ s/ $//;
		print "$genes_out\n";
	}
	print "\n";
}
