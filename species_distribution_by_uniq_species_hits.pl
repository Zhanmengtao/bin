#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 blast.tab > species_distribution_by_uniq_species_hits

USAGE
if (@ARGV==0){die $usage}

my %species_to_num;
my ( $geneid, $species );

foreach ( <> ) {
	chomp;
	if  ( /([^\t]+)\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t([^\t]+)/ ) {
		my $queryid = $1;
		my $subject_annotation = $2;

		if ( $geneid ne $queryid ) {
			$species = "";
			$geneid = $queryid;
		}

		my $annoted_species;
		if ( $subject_annotation =~ m/\[(\w+\s+[\w\.]+).*?\]/ ) {
			$annoted_species = $1;

			unless ( $species =~ /$annoted_species/ ) {
				$species_to_num{$annoted_species} ++;
				$species .= "$annoted_species\t";
			}
		}
	}
}

foreach ( sort by_num keys %species_to_num ) {
	print "$_\t$species_to_num{$_}\n";
}

sub by_num {
	$species_to_num{$b} <=> $species_to_num{$a};
}
