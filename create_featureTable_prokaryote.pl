#!/usr/bin/perl
use strict;
my $usage = <<USAGE;
Usage:
	perl $0 Locus_Tag dbname Nr.xls *.gff3 > prefix.tbl

USAGE

#2014-06-25	chenllianfu@foxmail.com

if (@ARGV==0) {die $usage}

my $locus_tag = shift @ARGV;
my $dbname = shift @ARGV;
# parsing nr annotation.
my $nrXLS = shift @ARGV;
open IN, $nrXLS or die $!;
my %nr;
while (<IN>) {
	if (/([^\t]+)\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t([^\t]+)/) {
		my ($queryid, $subject_annotation) = ($1, $2);
		if ($subject_annotation=~ s/(.*?) \[.*/$1/) {
			push @{$nr{$queryid}}, $subject_annotation;
		}
	}
}
close IN;

my %nr_best;
foreach my $seqID (keys %nr) {
	my @annot = @{$nr{$seqID}};
	my $annot = $annot[0];
	$annot =~ s/^\s*//;
	$annot =~ s/\s*$//;
	foreach (@annot) {
		next if /hypothetical protein/i;
		next if /predicted protein/i;
		next if /unknown/i;
		next if /homolog/i;
		next if /paralog/i;
		next if /analog/i;
		next if /Conserved domain protein/i;
		next if /uncharacterized protein/i;
		next if /truncat/;
		s/^\s*//;s/\s*$//;
		$annot = $_;
		last;
	}
	if ($annot =~ /hypothetical protein/i) { $annot = "hypothetical protein"; }
	if ($annot =~ /homolog|paralog|analog/i) { $annot = "hypothetical protein"; }
	$annot =~ s/, partial//;
	$nr_best{$seqID} = $annot;
}

# parsing GFF3 files
my %gff3;
foreach (@ARGV) {
	open IN, $_ or die $!;
	while (<IN>) {
		next if /^#/;
		my @fileds = split /\t/, $_;
		my $id;
		$id = $1 if $fileds[8] =~ /ID=([^;]+)/;
		$id =~ s/\.cds//;
		my ($seqID, $feature, $start, $end, $strand) = ($fileds[0], $fileds[2], $fileds[3], $fileds[4], $fileds[6]);
		$gff3{$seqID}{$id}{"feature"} = $feature;
		$gff3{$seqID}{$id}{"start"} = $start;
		$gff3{$seqID}{$id}{"end"} = $end;
		$gff3{$seqID}{$id}{"strand"} = $strand;
	}
}

# out feature table
my @genome_sequeces = keys %gff3;
my %sort_genome_sequeces;
foreach (@genome_sequeces) { if (/(\d+)/) { $sort_genome_sequeces{$_} = $1 } }
@genome_sequeces = sort { $sort_genome_sequeces{$a} <=> $sort_genome_sequeces{$b} } @genome_sequeces;

my ($gene_number, $rRNA16_number, $rRNA23_number, $rRNA5_number);
foreach my $seqID (@genome_sequeces) {
	print ">Feature $seqID\n";
	my @ids = keys %{$gff3{$seqID}};
	@ids = sort { $gff3{$seqID}{$a}{"start"} <=> $gff3{$seqID}{$b}{"start"} or $a cmp $b } @ids;

	foreach my $id (@ids) {
		my $feature = $gff3{$seqID}{$id}{"feature"};
		if ($feature eq "CDS" or $feature eq "tRNA" or $feature eq "rRNA") {
			$gene_number ++;
			my $gene_id;
			if (length $gene_number == 1) { $gene_id = $locus_tag . "_00000" . $gene_number; }
			if (length $gene_number == 2) { $gene_id = $locus_tag . "_0000" . $gene_number; }
			if (length $gene_number == 3) { $gene_id = $locus_tag . "_000" . $gene_number; }
			if (length $gene_number == 4) { $gene_id = $locus_tag . "_00" . $gene_number; }
			if (length $gene_number == 5) { $gene_id = $locus_tag . "_0" . $gene_number; }

			my $strand = $gff3{$seqID}{$id}{"strand"};
			my $start = $gff3{$seqID}{$id}{"start"};
			my $end = $gff3{$seqID}{$id}{"end"};
			if ($strand eq "-") { my $start_aa = $start; $start = $end; $end = $start_aa; }
			print "$start\t$end\tgene\n\t\t\tlocus_tag\t$gene_id\n";

			if ($feature eq "CDS") {
				my $product = "hypothetical protein";
				if ($nr_best{$id}) { $product = $nr_best{$id} }
				my $protein_id = "gnl\|$dbname\|$gene_id";
				print "$start\t$end\tCDS\n\t\t\tproduct\t$product\n\t\t\tprotein_id\t$protein_id\n";
			}
			elsif ($feature eq "tRNA") {
				my $product = $id;
				$product = $1 if $id =~ /(tRNA-\w\w\w)/;
				print "$start\t$end\ttRNA\n\t\t\tproduct\t$product\n";
			}
			elsif ($feature eq "rRNA") {
				my $product = $id;
				$product = "5S ribosomal RNA" if $product =~ m/5s_rRNA/;
				$product = "16S ribosomal RNA" if $product =~ m/16s_rRNA/;
				$product = "23S ribosomal RNA" if $product =~ m/23s_rRNA/;
				print "$start\t$end\trRNA\n\t\t\tproduct\t$product\n";
			}
		}
	}
}
