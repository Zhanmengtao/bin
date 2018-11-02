#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 RNA-Seq.sam > statistics.txt

iterms of statistics：
    total reads
    total Basepairs
    total mapped reads
    fully match
    <=2bp mismatch
    unique match
    multi-positon match
    total unmatched reads

USAGE
if (@ARGV==0) {die $usage}

my (%read, $total_reads, $total_basepairs, $total_mapped_reads, $fully_match, $lt2bp_mismatch, $unique_match, $multi_position_match, $total_unmatched_reads, $best_match_one_pos) = (0,0,0,0,0,0,0,0,0,0);
while (<>) {
    if (/^\@/) { next }

    split;
    my $flag = sprintf("%b",$_[1])+0;
    $flag = "000000" . $flag;
    my $read_number = $1 if $flag =~ m/(\d)\d\d\d\d\d\d$/;
    if (!exists $read{"$_[0]\t$read_number"}) {
        $read{"$_[0]\t$read_number"} = 1;

        $total_reads ++;
        $total_basepairs += length $_[9];

        if ($flag =~ /(\d)\d\d$/ && $1 == 0) {
            $total_mapped_reads ++;
            if (/NM:i:(\d+)/) {
                if ($1 == 0) {
                    $fully_match ++;
                }
                elsif ($1 <= 2) {
                    $lt2bp_mismatch ++;
                }
            }
            if (/ZS:i/ or (/NH:i:(\d+)/ && $1 >= 2)) {
                $multi_position_match ++;
            }
            else {
                $unique_match ++;
            }
            if (/NH:i:(\d+)/ && $1 == 1) {
                $best_match_one_pos ++;
            }
        }
        else {
            $total_unmatched_reads ++;
        }
    }
}
close IN;

# ratio calculation
my ($total_reads_ratio, $total_basepairs_ratio, $total_mapped_reads_ratio, $fully_match_ratio, $lt2bp_mismatch_ratio, $unique_match_ratio, $multi_position_match_ratio, $best_match_one_pos_ratio, $total_unmatched_reads_ratio);
$total_reads_ratio = '100.00%';
$total_basepairs_ratio = '100.00%';
$total_mapped_reads_ratio = ($total_mapped_reads / $total_reads) * 100;
$fully_match_ratio = ($fully_match / $total_reads) * 100;
$lt2bp_mismatch_ratio = ($lt2bp_mismatch / $total_reads) * 100;
$unique_match_ratio = ($unique_match / $total_reads) * 100;
$multi_position_match_ratio = ($multi_position_match / $total_reads) * 100;
$best_match_one_pos_ratio = ($best_match_one_pos / $total_reads) * 100;
$total_unmatched_reads_ratio = ($total_unmatched_reads / $total_reads) * 100;

# output
print "Total reads          \t$total_reads\t$total_reads_ratio\n";
print "Total Basepairs      \t$total_basepairs\t$total_basepairs_ratio\n";
printf "Total Mapped reads   \t$total_mapped_reads\t%2.2f\%\n", $total_mapped_reads_ratio;
printf "Fully match          \t$fully_match\t%2.2f\%\n", $fully_match_ratio;
printf "<=2bp mismatch       \t$lt2bp_mismatch\t%2.2f\%\n", $lt2bp_mismatch_ratio;
printf "Unique match         \t$unique_match\t%2.2f\%\n", $unique_match_ratio;
printf "Multi-position match \t$multi_position_match\t%2.2f\%\n", $multi_position_match_ratio;
printf "One Best Score match \t$best_match_one_pos\t%2.2f\%\n", $best_match_one_pos_ratio;
printf "Total unmatched reads\t$total_unmatched_reads\t%2.2f\%\n", $total_unmatched_reads_ratio;
