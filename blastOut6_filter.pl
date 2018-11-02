#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 blast.out6 query.fasta min_identity_percentage min_E-value min_coverage > blast.filter.out6

For example:
    perl $0 blast.out6 query.fasta 20.0 1e-9 0.2 > blast.filter.out6

USAGE
if (@ARGV==0){die $usage}

open IN, $ARGV[0] or die $!;
my %info;
while (<IN>) {
    next if m/^#/;
    my @field = split /\t/, $_;
    next if $field[2] < $ARGV[2];
    next if $field[10] > $ARGV[3];
    $info{$field[0]}{$field[1]} .= $_;
}
close IN;

open IN, $ARGV[1] or die $!;
my (%len, $seq_id);
foreach (<IN>) {
    chomp;
    if (m/^>(\S+)/) { $seq_id = $1; }
    else { $len{$seq_id} += length($_); }
}
close IN;

foreach my $query (sort keys %info) {
    foreach my $subject (sort keys %{$info{$query}}) {
        my $out = $info{$query}{$subject};
        chomp($out);
        my @line = split /\n/, $out;
        my @region;
        foreach (@line) {
            @_ = split /\t/, $_;
            push @region, "$_[6]\t$_[7]";
        }
        my $match_length_query = &match_length(@region);
        my $match_percentage = $match_length_query / $len{$query};
        #print "$out\t$match_length_query\t$len{$query}\t$match_percentage\n";
        print "$out\t$match_length_query\t$len{$query}\t$match_percentage\n" if $match_percentage >= $ARGV[4];
    }
}


sub match_length {
    my @inter_sorted_site;
    foreach (@_) {
        my @aaa = $_ =~ m/(\d+)/g;
        @aaa = sort { $a <=> $b } @aaa;
        push @inter_sorted_site, "$aaa[0]\t$aaa[1]";
    }
    @inter_sorted_site = sort { $a <=> $b } @inter_sorted_site;

    my $out_site_number;
    my $former_region = shift @inter_sorted_site;
    my @aaa = $former_region =~ m/(\d+)/g;
    $out_site_number += ($aaa[1] - $aaa[0] + 1);
    foreach (@inter_sorted_site) {
        my @former_region = $former_region =~ m/(\d+)/g;
        my @present_region = $_ =~ m/(\d+)/g;
        if ($present_region[0] > $former_region[1]) {
            $out_site_number += ($present_region[1] - $present_region[0] + 1);
            $former_region = $_;
        }
        elsif ($present_region[1] > $former_region[1]) {
            $out_site_number += ($present_region[1] - $former_region[1]);
            $former_region = $_;
        }
        else {
            next
        }
    }
    return $out_site_number;
}
