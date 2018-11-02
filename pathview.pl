#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 query.ko up_gene.list down_gene.list CPUs

    quer.ko format
#Gene ID       KO number
Gene01        K10798
Gene02        K13110

    up_gene.list format
#DEG ID        P_value
Gene01        0.03
Gene02        0.0004

ParaFly and kaas2pathwayEnrichmentAnalysis.pl are required.
R package pathview is required.

USAGE
if (@ARGV==0){die $usage}

my $pwd = `pwd`;
chomp($pwd);
my $kaas_annot = "$pwd/$ARGV[0]";
my $up_genes = "$pwd/$ARGV[1]";
my $down_genes = "$pwd/$ARGV[2]";

my (%up_genes, %down_genes, %ko_num);
open IN, $up_genes or die $!;
while (<IN>) {
    next if m/^\s*$/;
    next if m/^#/;
    chomp;
    @_ = split /\s+/;
    $up_genes{$_[0]} = $_[1];
}
close IN;

open IN, $down_genes or die $!;
while (<IN>) {
    next if m/^\s*$/;
    next if m/^#/;
    chomp;
    @_ = split /\s+/;
    $down_genes{$_[0]} = $_[1];
}
close IN;

open IN, $kaas_annot or die $!;
while (<IN>) {
    if (m/(\S+)\s+(K\d+)/) {
        $ko_num{$2}{"ok"} = 1;
        if (exists $up_genes{$1}) {
            if (exists $ko_num{$2}{"up"}) {
                $ko_num{$2}{"up"} = $up_genes{$1} if $ko_num{$2}{"up"} > $up_genes{$1};
            }
            else {
                $ko_num{$2}{"up"} = $up_genes{$1};
            }
        }
        if (exists $down_genes{$1}) {
            if (exists $ko_num{$2}{"down"}) {
                $ko_num{$2}{"down"} = $ko_num{$2}{"down"} if $ko_num{$2}{"down"} > $ko_num{$2}{"down"};
            }
            else {
                $ko_num{$2}{"down"} = $ko_num{$2}{"down"};
            }
        }
    }
}
close IN;

mkdir "$pwd/pathview" unless -e "$pwd/pathview";
open OUT, '>', "$pwd/pathview/pathview_in.txt" or die $!;
print OUT "\"DCIS_1\"\t\"DCIS_2\"\n";
foreach my $ko_num (sort keys %ko_num) {
    my $out = "$ko_num\t";
    if (exists $ko_num{$ko_num}{"up"} ) {
        if ($ko_num{$ko_num}{"up"} < 0.001) { $out .= "1\t"; }
        elsif ($ko_num{$ko_num}{"up"} < 0.05) { $out .= "0.5\t"; }
    }
    else {
        $out .= "0\t";
    }
    if (exists $ko_num{$ko_num}{"down"}) {
        if ($ko_num{$ko_num}{"down"} < 0.001) { $out .= "-1\n"; }
        elsif ($ko_num{$ko_num}{"down"} < 0.05) { $out .= "-0.5\n"; }
    }
    else {
        $out .= "0\n";
    }
    print OUT $out;
}
close OUT;

# enrichment analysis of up-regulated genes
mkdir "$pwd/up_regulated_genes_enrichment" unless -e "$pwd/up_regulated_genes_enrichment";
chdir "$pwd/up_regulated_genes_enrichment";
symlink "$pwd/ko00001.keg", "ko00001.keg";
my $cmdString = "kaas2pathwayEnrichmentAnalysis.pl $kaas_annot $up_genes";
print STDERR "CMD: $cmdString\n";
(system $cmdString) == 0 or die "CMD Failed: $cmdString";

# enrichment analysis of down_regulated genes
mkdir "$pwd/down_regulated_genes_enrichment" unless -e "$pwd/down_regulated_genes_enrichment";
chdir "$pwd/down_regulated_genes_enrichment";
symlink "$pwd/ko00001.keg", "ko00001.keg";
$cmdString = "kaas2pathwayEnrichmentAnalysis.pl $kaas_annot $down_genes";
print STDERR "CMD: $cmdString\n";
(system $cmdString) == 0 or die "CMD Failed: $cmdString";

my (@up_enrichment_pathway, @down_enrichment_pathway);
open IN, "$pwd/up_regulated_genes_enrichment/enrichment.keggPathway" or die $!;
while (<IN>) {
    if (m/^(ko\d+)/) { push @up_enrichment_pathway, $1; }
}
close IN;
open IN, "$pwd/down_regulated_genes_enrichment/enrichment.keggPathway" or die $!;
while (<IN>) {
    if (m/^(ko\d+)/) { push @down_enrichment_pathway, $1; }
}
close IN;

mkdir "$pwd/pathview/up_enrichment_pathway_maps/" unless -e "$pwd/pathview/up_enrichment_pathway_maps/";
mkdir "$pwd/pathview/up_enrichment_pathway_maps/R_scripts/" unless -e "$pwd/pathview/up_enrichment_pathway_maps/R_scripts/";
chdir "$pwd/pathview/up_enrichment_pathway_maps/";
open CMD, '>', "$pwd/pathview/up_enrichment_pathway_maps/R_scripts.commands"  or die $!;
foreach my $pathway_id (@up_enrichment_pathway) {
    open OUT, '>', "$pwd/pathview/up_enrichment_pathway_maps/R_scripts/$pathway_id.R" or die $!;
    print OUT "library('pathview')\n";
    print OUT "pathview_in <- read.table(\"$pwd/pathview/pathview_in.txt\", header=TRUE)\n";
    print OUT "pv.out <- pathview(gene.data = pathview_in, pathway.id=\"$pathway_id\", species = \"ko\")\n";
    print CMD "cat $pwd/pathview/up_enrichment_pathway_maps/R_scripts/$pathway_id.R | R --vanilla --slave\n";
    close OUT;
}
close CMD;
$cmdString = "ParaFly -c R_scripts.commands -CPU $ARGV[3]";
system $cmdString;

mkdir "$pwd/pathview/down_enrichment_pathway_maps/" unless -e "$pwd/pathview/down_enrichment_pathway_maps/";
mkdir "$pwd/pathview/down_enrichment_pathway_maps/R_scripts/" unless -e "$pwd/pathview/down_enrichment_pathway_maps/R_scripts/";
chdir "$pwd/pathview/down_enrichment_pathway_maps/";
open CMD, '>', "$pwd/pathview/down_enrichment_pathway_maps/R_scripts.commands"  or die $!;
foreach my $pathway_id (@down_enrichment_pathway) {
    open OUT, '>', "$pwd/pathview/down_enrichment_pathway_maps/R_scripts/$pathway_id.R" or die $!;
    print OUT "library('pathview')\n";
    print OUT "pathview_in <- read.table(\"$pwd/pathview/pathview_in.txt\", header=TRUE)\n";
    print OUT "pv.out <- pathview(gene.data = pathview_in, pathway.id=\"$pathway_id\", species = \"ko\")\n";
    print CMD "cat $pwd/pathview/down_enrichment_pathway_maps/R_scripts/$pathway_id.R | R --vanilla --slave\n";
    close OUT;
}
close CMD;
$cmdString = "ParaFly -c R_scripts.commands -CPU $ARGV[3]";
system $cmdString;
