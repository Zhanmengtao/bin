#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
	perl $0 blast.xml max_hit_number e_value match_percentage
For example:
	perl $0 nr.xml 20 1e-5 0.2

	if max_hit_number == 0, all hits would be considered.

USAGE
if (@ARGV!=4) {die $usage}

my $max_hit_number = $ARGV[1];
my $E_value = $ARGV[2];
my $Match_percentage = $ARGV[3];

print "Query_id\tQuery_length\tQuery_start\tQuery_end\tSubject_id\tSubject_length\tSubject_start\tSubject_end\tIdentity\tPositive\tGap\tAlign_length\tScore\tE_value\tMatch_percentage\tSubject_annotation\n";

open IN, $ARGV[0] or die $!;
my ($one, $record);
while (<IN>) {
	if (/<Iteration>/) { $record = 1 }
	next unless $record == 1;
	$one .= $_;
	if (/<\/Iteration>/) {
		my $query_id = $1 if $one =~ m#<Iteration_query-def>(.+?)</Iteration_query-def>#;
		my $query_length = $1 if $one =~ m#<Iteration_query-len>(\d+?)</Iteration_query-len>#;

		my @hit = $one =~ m#<Hit>(.*?)</Hit>#gs;
		my $hitNum = 0;

		foreach (@hit) {
			last if  ($max_hit_number && $hitNum > $max_hit_number);

			my $query_start = $1 if $_ =~ m#<Hsp_query-from>(\d+?)</Hsp_query-from>#;
			my $query_end = $1 if $_ =~ m#<Hsp_query-to>(\d+?)</Hsp_query-to>#;
			my $subject_id = $1 if $_ =~ m#<Hit_id>(.+?)</Hit_id>#;
			my $subject_length = $1 if $_ =~ m#<Hit_len>(\d+?)</Hit_len>#;
			my $subject_start = $1 if $_ =~ m#<Hsp_hit-from>(\d+?)</Hsp_hit-from>#;
			my $subject_end = $1 if $_ =~ m#<Hsp_hit-to>(\d+?)</Hsp_hit-to>#;
			my $hsp_identity = $1 if $_ =~ m#<Hsp_identity>(\d+?)</Hsp_identity>#;
			my $hsp_positive = $1 if $_ =~ m#<Hsp_positive>(\d+?)</Hsp_positive>#;
			my $hsp_gaps = $1 if $_ =~ m#<Hsp_gaps>(\d+?)</Hsp_gaps>#;
			my $hsp_align_len = $1 if $_ =~ m#<Hsp_align-len>(\d+?)</Hsp_align-len>#;
			next if $hsp_align_len ==0;
			my $identity = $hsp_identity / $hsp_align_len;
			my $positive = $hsp_positive / $hsp_align_len;
			my $gap = $hsp_gaps / $hsp_align_len;
			my $align_length =  $hsp_align_len;
			my $score = $1 if $_ =~ m#<Hsp_bit-score>(\S+?)</Hsp_bit-score>#;
			my $e_value = $1 if $_ =~ m#<Hsp_evalue>(\S+?)</Hsp_evalue>#;
			my $subject_annotation = $1 if $_ =~ m#<Hit_def>(.*?)</Hit_def>#;

			my (@hsp_query, @hsp_subject);
			my @hsp = $_ =~ m#<Hsp>(.*?)</Hsp>#gs;
			foreach (@hsp) {
				my $query_from = $1 if m#<Hsp_query-from>(\d+?)</Hsp_query-from>#;
				my $query_to = $1 if m#<Hsp_query-to>(\d+?)</Hsp_query-to>#;
				my $hit_from = $1 if m#<Hsp_hit-from>(\d+?)</Hsp_hit-from>#;
				my $hit_to = $1 if m#<Hsp_hit-to>(\d+?)</Hsp_hit-to>#;

				push @hsp_query, "$query_from\t$query_to";
				push @hsp_subject, "$hit_from\t$hit_to";
			}

			my $match_length_query = &match_length(@hsp_query);
			my $match_length_subject = &match_length(@hsp_subject);
			my @lengths = sort {$a <=> $b} ($query_length, $subject_length);
			my @match_lengths = sort {$a <=> $b} ($match_length_query, $match_length_subject);
			my $match_percentage = $match_lengths[0] / $lengths[0];

			if ($e_value <= $E_value && $Match_percentage <= $match_percentage) {
				$hitNum ++;
				printf "$query_id\t$query_length\t$query_start\t$query_end\t$subject_id\t$subject_length\t$subject_start\t$subject_end\t%0.4f\t%0.4f\t%0.4f\t$align_length\t$score\t$e_value\t%0.4f\t$subject_annotation\n",$identity,$positive,$gap,$match_percentage;
			}
		}
		$record = 0;
		$one = "";
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

#sub match_length {
#	my @all_site;
#	foreach (@_) {
#		my @aaa = $_ =~ m/(\d+)/g;
#		@aaa = sort { $a <=> $b } @aaa;
#		foreach ($aaa[0]..$aaa[1]) {
#			$all_site[$_] = 1
#		}
#	}
#	my @length = grep { /1/ } @all_site;
#	my $length = @length;
#	return $length;
#}
