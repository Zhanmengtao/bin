#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 key_words

	This program is designed to kill procedures by key words.

USAGE
if (@ARGV==0){die $usage}

my @temp = `ps -aux | grep $ARGV[0]`;
#print @temp;

my $num = 0;
foreach ( @temp ) {
	$_ =~ m/^\S+\s+(\d+)/;
	`kill $1` unless $_ =~ m/kill_procedures_by_name_clf/;
	$num++;
	$_ =~ m/\d+:\d+\s+\d+:\d+\s+(.*?)$/;
	print "$num	$1\n";
}

print "total	$num procedures have been killed.\n";
