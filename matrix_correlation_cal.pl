#!/usr/bin/perl
use strict;
use Statistics::Basic  qw / correlation  vector /;

my $usage = <<USAGE;
Usage:
    perl $0 matrix.tab

USAGE
if (@ARGV==0){die $usage}

open IN, $ARGV[0] or die "Can not open file $ARGV[0] $!";
$_ = <IN>;
s/^.*?\t//;
s/\s*$//;
my @samples = split /\s+/;

my %data;
while (<IN>) {
	s/^\S+\s+//;
	chomp;
	@_ = split /\t/;
	my @samples1 = @samples;
	foreach (@_) {
		my $sample = shift @samples1;
		push @{$data{$sample}}, $_;
	}
}

my @samples2 = @samples;
my $out = join "\t", @samples;
print "\t$out\n";
foreach (@samples) {
	print "$_\t";
	my @data1 = @{$data{$_}};
	my $v1 = vector(@data1);
	my $out = "";
	foreach (@samples2) {
		my @data2 = @{$data{$_}};
		my $v2 = vector(@data2);
		my $cor = correlation($v1, $v2);
		$out .= "$cor\t";
	}
	$out =~ s/\t$/\n/;
	print $out;
}
