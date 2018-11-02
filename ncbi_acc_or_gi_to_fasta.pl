#!/usr/bin/perl
use strict;

my $usage=<<USAGE;
Usage:
	perl $0 accession_or_gi_list parallel_num > target.fasta

USAGE
if (@ARGV ==0) {die $usage}

open OUT, '>', "/tmp/get_seq_by_gcc.pl" or die $!;
print OUT 'use Bio::DB::GenBank;
$db_obj = Bio::DB::GenBank->new();
$seq_obj = $db_obj->get_Seq_by_version($ARGV[0]);
$id = $seq_obj->display_id();
$seq = $seq_obj->seq();
$desc = $seq_obj->desc();
print ">$ARGV[0] \| $id \| $desc\n$seq\n"
';
close OUT;

open IN, $ARGV[0] or die $!;
my %list;
while (<IN>) {
	chomp;
	$list{$_} = 1;
}
close IN;

open COM, '>', "$ARGV[0].commands" or die $!;
mkdir "$ARGV[0].get_seq_by_acc.tmp" unless -e "$ARGV[0].get_seq_by_acc.tmp";
foreach (keys %list) {
	print COM "perl /tmp/get_seq_by_gcc.pl $_ > $ARGV[0].get_seq_by_acc.tmp/$_.fa\n"
}
close COM;

system "ParaFly -c $ARGV[0].commands -CPU $ARGV[1] &> /dev/null";

foreach (keys %list) {
	open IN, "$ARGV[0].get_seq_by_acc.tmp/$_.fa";
	print <IN>;
}

unlink "/tmp/get_seq_by_gcc.pl";
