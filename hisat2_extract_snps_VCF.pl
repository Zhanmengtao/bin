#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 snp_indel.vcf > variants.snp_file_for_hisat2

USAGE
if (@ARGV==0){die $usage}

open IN, $ARGV[0] or die $!;
while (<IN>) {
    next if m/^#/;
    my ($chromo_id, $position, $snp_id, $ref, $snp) = split /\t/;
    $position --;
    my $ref_len = length $ref;
    my @alt = split /,/, $snp;
    # 去除为 '*' 的 alt
    my %alt;
    foreach (@alt) {
        $alt{$_} = 1 if $_ ne '*';
    }
    @alt = sort keys %alt;

    if (@alt == 1) {
        my $alt_len = length $alt[0];
        if ($ref_len ==  1 && $alt_len ==  1) {
            print "$snp_id\tsingle\t$chromo_id\t$position\t$alt[0]\n";
        }
        elsif ($alt_len > $ref_len) {
            if ($alt[0] =~ s/^$ref//) {
                $position = $position + $ref_len - 1;
                print "$snp_id\tinsertion\t$chromo_id\t$position\t$alt[0]\n";
            }
        }
        elsif ($alt_len < $ref_len) {
            if ($ref =~ m/^$alt[0]/) {
                $position = $position + $alt_len - 1;
                my $deletion_length = $ref_len - $alt_len;
                print "$snp_id\tdeletion\t$chromo_id\t$position\t$deletion_length\n";
            }
        }
    }
    else {
        my $alt_number = 0;
        foreach my $alt (@alt) {
            my $alt_len = length($alt);
            if ($ref_len ==  1 && $alt_len ==  1) {
                print "$snp_id.$alt_number\tsingle\t$chromo_id\t$position\t$alt\n";
            }
            elsif ($alt_len > $ref_len) {
                if ($alt =~ s/^$ref//) {
                    $position = $position + $ref_len - 1;
                    print "$snp_id.$alt_number\tinsertion\t$chromo_id\t$position\t$alt\n";
                }
            }
            elsif ($alt_len < $ref_len) {
                if ($ref =~ m/^$alt/) {
                    $position = $position + $alt_len - 1;
                    my $deletion_length = $ref_len - $alt_len;
                    print "$snp_id.$alt_number\tdeletion\t$chromo_id\t$position\t$deletion_length\n";
                }
            }
            $alt_number ++;
        }
    }
}
