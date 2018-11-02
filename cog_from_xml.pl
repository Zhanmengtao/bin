#!/usr/bin/perl
use strict;
use File::Basename;
use Getopt::Long;

my $usage = <<USAGE;
Usage:
    perl $0 [options] cog.xml

    --E_value <float>    default: 1e-5
    设置Blast结果的E_value阈值。仅选取blast结果中E值>=此设定值的比对结果。

    --coverage <float>    default: 0.2
    设置Blast结果的覆盖率阈值。仅选取blast结果中query序列和hit序列匹配区域的覆盖率>=此设定值的比对结果。

USAGE
if (@ARGV==0){die $usage}

my ($E_value, $coverage);
GetOptions(
        "E_value:s" => \$E_value,
        "coverage:s" => \$coverage,
);
$E_value ||= 1e-5;
$coverage ||= 0.2;

my $cog_script_dir = dirname($0);
my $whogFile = "$cog_script_dir/cog/whog";
my $funFile = "$cog_script_dir/cog/fun.txt";
my (%proteinid_to_cogid, %cog_to_fun, %cog_to_desc, %fun_code_to_category);

open WHOG, '<', $whogFile or die $!;
$/ = "_______";
while (<WHOG>) {
    my $cogID;
    if (/\[(\w+)\]\s+(COG\d+)\s+(.+)\n/) {
        $cog_to_fun{$2} = $1;
        $cog_to_desc{$2} = $3;
        $cogID = $2;
    }
    my @lines = split /\n/;
    foreach (@lines) {
        chomp;
        if (s/^\s+\w+:\s+//) {
            my @proteins = split /\s+/;
            foreach (@proteins) { $proteinid_to_cogid{$_} = $cogID }
        }
    }
}
$/ = "\n";
close WHOG;
print STDERR "Whog File: $whogFile\n";

open FUN, '<', $funFile or die $!;
while (<FUN>) {
    if (/\[(\w)\]\s(.+)/) {
        $fun_code_to_category{$1} = $2
    }
}
close KOG;
print STDERR "Fun File: $funFile\n";

my $name = basename $ARGV[0];
$name =~ s/.xml//;
my $xls = "$name".".xls";
my $geneannot = $name . "_cog_gene_annot.xls";
my $cogclass = $name . "_cog_class_annot.xls";

open XLS, '>', $xls;
open GENEANNOT, '>', $geneannot;

print XLS "Query_id\tQuery_length\tQuery_start\tQuery_end\tSubject_id\tSubject_length\tSubject_start\tSubject_end\tQuery_frame\tIdentity\tPositive\tGap\tAlign_length\tScore\tE_value\tMatch_percentage\tKOG_ID\tSubject_annotation\tFunction_Code\tFunctional-Categories\n";
print GENEANNOT "Query_id\tSubject_id\tScore\tE_value\tMatch_percentage\tKOG_ID\tSubject_annotation\tFunction_Code\tFunctional-Categories\n";

my $one;

open IN, $ARGV[0] or die $!;
print STDERR "Input File: $ARGV[0]\n";
while ( <IN> ) {
    $one .= $_;
    if ( /<\/BlastOutput>/ ) {
        my $Query_id = $1 if $one =~ m#<Iteration_query-def>(.+?)</Iteration_query-def>#;
        my $Query_length = $1 if $one =~ m#<Iteration_query-len>(\S+?)</Iteration_query-len>#;

        my @hit = $one =~ m#<Hit>(.*?)</Hit>#gs;

        foreach ( @hit ) {
            my $Query_start = $1 if $_ =~ m#<Hsp_query-from>(\S+?)</Hsp_query-from>#;
            my $Query_end = $1 if $_ =~ m#<Hsp_query-to>(\S+?)</Hsp_query-to>#;
            my $Subject_id = $1 if $_ =~ m#<Hit_id>(\S+?)</Hit_id>#;
            my $Subject_length = $1 if $_ =~ m#<Hit_len>(\S+?)</Hit_len>#;
            my $Subject_start = $1 if $_ =~ m#<Hsp_hit-from>(\S+?)</Hsp_hit-from>#;
            my $Subject_end = $1 if $_ =~ m#<Hsp_hit-to>(\S+?)</Hsp_hit-to>#;
            my $Query_frame = $1 if $_ =~ m#<Hsp_query-frame>(\S+?)</Hsp_query-frame>#;
            my $Hsp_identity = $1 if $_ =~ m#<Hsp_identity>(\S+?)</Hsp_identity>#;
            my $Hsp_positive = $1 if $_ =~ m#<Hsp_positive>(\S+?)</Hsp_positive>#;
            my $Hsp_gaps = $1 if $_ =~ m#<Hsp_gaps>(\S+?)</Hsp_gaps>#;
            my $Hsp_align_len = $1 if $_ =~ m#<Hsp_align-len>(\S+?)</Hsp_align-len>#;
            my $Identity = $Hsp_identity / $Hsp_align_len;
            my $Positive = $Hsp_positive / $Hsp_align_len;
            my $Gap = $Hsp_gaps / $Hsp_align_len;
            my $Align_length = $Hsp_align_len;
            my $Score = $1 if $_ =~ m#<Hsp_bit-score>(\S+?)</Hsp_bit-score>#;
            my $Evalue = $1 if $_ =~ m#<Hsp_evalue>(\S+?)</Hsp_evalue>#;
            my $Subject_annotation = $1 if $_ =~ m#<Hit_def>(.*?)</Hit_def>#;

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
            my @lengths = sort {$a <=> $b} ($Query_length, $Subject_length);
            my @match_lengths = sort {$a <=> $b} ($match_length_query, $match_length_subject);
            my $match_percentage = $match_lengths[0] / $lengths[0];

            my $cog_id = $proteinid_to_cogid{$Subject_id};
            my $cog_fun = $cog_to_fun{$cog_id};
            my $Subject_annotation = $cog_to_desc{$cog_id};

            my $functional_categories;

            if ( $cog_fun ) {
                my @cog_fun = split //, $cog_fun;
                foreach ( @cog_fun ) {
                    $functional_categories .= $fun_code_to_category{$_};
                    $functional_categories .= " ; ";
                }
            }

            if ( $Evalue <= $E_value && $match_percentage >= $coverage ) {
                printf GENEANNOT "$Query_id\t$Subject_id\t$Score\t$Evalue\t%0.4f\t$cog_id\t$Subject_annotation\t$cog_fun\t$functional_categories\n",$match_percentage;
                printf XLS "$Query_id\t$Query_length\t$Query_start\t$Query_end\t$Subject_id\t$Subject_length\t$Subject_start\t$Subject_end\t$Query_frame\t%0.4f\t%0.4f\t%0.4f\t$Align_length\t$Score\t$Evalue\t%0.4f\t$cog_id\t$Subject_annotation\t$cog_fun\t$functional_categories\n",$Identity,$Positive,$Gap,$match_percentage;
            }
        }
        $one = "";
    }
}
close IN;
close XLS;
close GENEANNOT;
print "Output File 1 : $xls\n";
print "Output File 2 : $geneannot\n";

my %fun_code_to_genes;
open GENEANNOT, '<', $geneannot || die "Can not open the file $geneannot! ($!)";
<GENEANNOT>;
while ( <GENEANNOT> ) {
    if ( /^(.*?)\t.*?\t.*?\t.*?\t.*?\t(.*?)\t.*?\t(.*?)\t/ ) {
        foreach (split //, $3) {
            $fun_code_to_genes{$_}{$1}{$2} = 1;
        }
    }
}
close GENEANNOT;

my @fun_codes = qw( A B C D E F G H I J K L M N O P Q R S T U V W Y Z );
open COGCLASS, '>', $cogclass;
print COGCLASS "Code\tFunctional-Categories\tGene-Number\tGenes\n";
foreach my $function_code ( @fun_codes ) {
    if ( exists $fun_code_to_genes{$function_code} ) {
        my @genes = keys %{$fun_code_to_genes{$function_code}};
        my $gene_number = @genes;
        my $genes_out;
        foreach my $gene_name (@genes) {
            $genes_out .= "$gene_name,";
            foreach (keys %{$fun_code_to_genes{$function_code}{$gene_name}}) {
                $genes_out .= "$_,";
            }
            $genes_out =~ s/,$/ /;
        }
        print COGCLASS "$function_code\t$fun_code_to_category{$function_code}\t$gene_number\t$genes_out\n";
    }
    else {
        print COGCLASS "$function_code\t$fun_code_to_category{$function_code}\t0\n";
    }
}
close COGCLASS;
print "Output File 3 : $cogclass\n";

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
