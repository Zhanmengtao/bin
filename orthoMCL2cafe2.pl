#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
	perl $0 groups.txt > orthomcl2cafe.tab

USAGE
if (@ARGV==0) { die $usage }

open IN, $ARGV[0] or die $!;
my (%gene_family_size, %taxon, @group_ID);
while (<IN>) {
	if (/^([^:]+?):/) {
		my $group_ID = $1; push @group_ID, $group_ID;
		my @genes = m/(\w\w\w\w)\|/g;
		foreach (@genes) {
			$taxon{$_} = 1;
			$gene_family_size{$group_ID}{$_} ++;
		}
	}
}

my @taxon = sort { $a cmp $b } keys %taxon;
print "Description\tID";
foreach (@taxon) {
	print "\t$_";
}
print "\n";

foreach my $group_ID (@group_ID) {
	my ($out, $num, $number);
	$out .= "\t$group_ID";
	foreach (@taxon) {
                if (exists $gene_family_size{$group_ID}{$_}) {
					$num += $gene_family_size{$group_ID}{$_};
					$out .= "\t$gene_family_size{$group_ID}{$_}";
					$number ++;
                }
                else {
					$out .= "\t0";
                }
	}
	$out .= "\n";
	print $out if $num >= (@taxon / 2) && $number >= 2;
}
