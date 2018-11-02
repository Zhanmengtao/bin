#!/usr/bin/perl
use strict;
my $usage = <<USAGE;
Usage:
    perl $0 blastRemoteResult.tab > blastRemoteResult.sort.tab

USAGE
if (@ARGV == 0){die $usage}

my $head = <>;
print $head;

my %a;
while (<>) {
	if (/([^\t]+)/) { $a{$1} .= $_ }
}

foreach (sort keys %a) {
	my @a = split /\n/, $a{$_};
	my %b;
	foreach (@a) { if (/(\S+)\t[^\t]+$/) { $b{$_} = $1 } }
	@a = sort { $b{$a} <=> $b{$b} } @a;
	my $num = 0;
	foreach (@a) { $num ++; last if $num > 20; print "$_\n" }
}

