#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 out.cafe species_name > cafe_out.txt 2> cafe_species.txt

USAGE
if (@ARGV==0){die $usage}

open IN, $ARGV[0] or die $!;
my $species = $ARGV[1] or die "No species was given !\n";

# Tree
$_ = <IN>;
print;
# Lambda
$_ = <IN>;
print;

# IDs of nodes
$_ = <IN>;
print "$_\n";;
my %id_node = /(\w+?)<(\d+?)>/g;
my @id_node = sort keys %id_node;

# node order
$_ = <IN>;
my @nodes = /(\d+)/g;

# header
my $out = "\t";
my (%id_node_values, %no_id_node_values);
foreach (@id_node) { $id_node_values{$id_node{$_}} = 1; $out .= "$id_node{$_}\t"; }
foreach (@nodes) { $no_id_node_values{$_} = 1 if ! exists $id_node_values{$_}; }
foreach (sort {$a <=> $b} keys %no_id_node_values) { $out .= "$_\t"; }
$out =~ s/\t$/\n/;
print $out;
my $out = "\t";
foreach (@id_node) { $out .= "$_\t"; }
$out =~ s/\t$/\n/;
print $out;

# Average expansion rate
$_ = <IN>;
$_ = <IN>;
my @aver_expa = m/([\d\.-]+)/g;
my %aver_expa;
foreach (@nodes) { $aver_expa{$_} = shift @aver_expa; }
$out = "Average Expansion\t";
foreach (@id_node) { $out .= "$aver_expa{$id_node{$_}}\t"; }
foreach (sort {$a <=> $b} keys %no_id_node_values) { $out .= "$aver_expa{$_}\t"; }
$out =~ s/\t$/\n/;
print $out;

# Expansion gene family number
$_ = <IN>;
my @expansion = m/(\d+)/g;
my %expansion;
foreach (@nodes) { $expansion{$_} = shift @expansion; }
$out = "Expansion\t";
foreach (@id_node) { $out .= "$expansion{$id_node{$_}}\t"; }
foreach (sort {$a <=> $b} keys %no_id_node_values) { $out .= "$expansion{$_}\t"; }
$out =~ s/\t$/\n/;
print $out;

# No Change
$_ = <IN>;
my @remain = m/(\d+)/g;
my %remain;
foreach (@nodes) { $remain{$_} = shift @remain; }
$out = "Remain\t";
foreach (@id_node) { $out .= "$remain{$id_node{$_}}\t"; }
foreach (sort {$a <=> $b} keys %no_id_node_values) { $out .= "$remain{$_}\t"; }
$out =~ s/\t$/\n/;
print $out;

# Decrease
$_ = <IN>;
my @decrease = m/(\d+)/g;
my %decrease;
foreach (@nodes) { $decrease{$_} = shift @decrease; }
$out = "Decrease\t";
foreach (@id_node) { $out .= "$decrease{$id_node{$_}}\t"; }
foreach (sort {$a <=> $b} keys %no_id_node_values) { $out .= "$decrease{$_}\t"; }
$out =~ s/\t$/\n/;
print $out;

# header
$_ = <>;
my $out = "\nGene family ID\tGene family size\tP-value\t";
foreach (@id_node) { $out .= "$id_node{$_}\t"; }
foreach (sort {$a <=> $b} keys %no_id_node_values) { $out .= "$_\t"; }
$out =~ s/\t$/\n/;
print $out;
my $out = "Gene family ID\tGene family size\tP-value\t";
foreach (@id_node) { $out .= "$_\t"; }
$out =~ s/\t$/\n/;
print $out;

print STDERR "family_id\tdivergence_size\tspecies_size\tP-value\n";
my ($total_expansion_family_number, $expansion_family_number);
# cafe main out
while (<IN>) {
	@_ = split /\t/;
	print "$_[0]\t$_[1]\t$_[2]\t";

	my ($family_id, $family_tree, $p_value, $p_values) = @_;
	my @p_values = $p_values =~ m/([\d\.-]+)/g;
	my %p_values;
	foreach (@nodes) { $p_values{$_} = shift @p_values; }
	my $out = "";
	foreach (@id_node) { $out .= "$p_values{$id_node{$_}}\t"; }
	foreach (sort {$a <=> $b} keys %no_id_node_values) { $out .= "$p_values{$_}\t"; }
	$out =~ s/\t$/\n/;
	print $out;

	if ($family_tree =~ m/$species\_(\d+).*?\)_(\d+)/) {
		if ($p_value <= 0.05) {
			$total_expansion_family_number ++;
			if ($p_values{$id_node{$species}} <= 0.05) {
				$expansion_family_number ++;
				print STDERR "$family_id\t$2\t$1\t$p_value\n";
			}
		}
	}
}
print STDERR "\nTotal Expansion/Contractions family number: $total_expansion_family_number\n";
print STDERR "${species}'s Expansion/Contractions family number: $expansion_family_number\n";
