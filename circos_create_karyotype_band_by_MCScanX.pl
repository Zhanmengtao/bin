#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 data/nc_cp.gff data/nc_cp.collinearity_interspecific 15 1 > karyotype.band_for_circos.txt
	对MCScan的结果进行分析，取第一个物种中共线性基因数目>=15的block作为circos的band信息。
USAGE
if (@ARGV==0){die $usage}

open GFF, $ARGV[0] or die $!;
my %gene_stats;
while (<GFF>) {
	chomp;
	@_ = split /\t/;
	$gene_stats{$_[1]}{"scaffold"} = $_[0];
	$gene_stats{$_[1]}{"start"} = $_[2];
	$gene_stats{$_[1]}{"end"} = $_[3];
}
close GFF;

open MCSCAN, $ARGV[1] or die $!;
my (%mcscan, @mcscan);
while (<MCSCAN>) {
	next if m/^#/;
	if (m/^\s*(\d+)-\s*\d+:\s*(\S+)\s*(\S+)/) {
		push @mcscan, $1 if ! exists $mcscan{$1};
		$mcscan{$1}{"size"} ++;
		$mcscan{$1}{1}{$2} = 1;
		$mcscan{$1}{2}{$3} = 1;
	}
}
close MCSCAN;

foreach my $aligment_code (@mcscan) {
	if ($mcscan{$aligment_code}{"size"} >= $ARGV[2]) {
		my ($scaffold_name, @sites);
		foreach (keys %{$mcscan{$aligment_code}{$ARGV[3]}}) {
			push @sites, $gene_stats{$_}{"start"};
			push @sites, $gene_stats{$_}{"end"};
			$scaffold_name = $gene_stats{$_}{"scaffold"};
		}
		@sites = sort {$a <=> $b} @sites;
		my $start = $sites[0] - 1;
		my $end = $sites[-1] - 1;
		print "band\t$scaffold_name\tmcscan_$aligment_code\tmcscan_$aligment_code\t$start\t$end\tblack\n";
	}
}
