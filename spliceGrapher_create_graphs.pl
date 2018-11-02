#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 GeneModels.gff3 SpiceGrapher_out.gtf SpiceGrapher_filter.Sortedbam Out_dir CPUs

USAGE
if (@ARGV==0){die $usage}

my $output_dir = `pwd`;
chomp($output_dir);
my $bam_file;
my $gtf_file;
if ($ARGV[1] =~ m/^\//) { $gtf_file = $ARGV[1]; }
else { $gtf_file = "$output_dir/$ARGV[1]"; }
if ($ARGV[2] =~ m/^\//) { $bam_file = $ARGV[2]; }
else { $bam_file = "$output_dir/$ARGV[2]"; }
if ($ARGV[3] =~ m/^\//) { $output_dir = $ARGV[3]; }
else { $output_dir = "$output_dir/$ARGV[3]"; }

mkdir $output_dir unless -e $output_dir;
# 对每个基因分别进行画图，需要的文件：
# 1. 画基因的可变剪接图，需要原始的gff3文件
# 2. 画RNA-Seq中有表达的可变剪接，需要SpliceGrapher的结果gtf文件
# 3. 画剪接位点信息和reads深度图，需要sam文件
# 4. 进行画图，需要将上述信息和画图参数整理到一个参素文件中
# 5. 最后，将结果输出到指定文件夹中。
mkdir "$output_dir/gff3_per_gene" unless -e "$output_dir/gff3_per_gene";
mkdir "$output_dir/sam_per_gene" unless -e "$output_dir/sam_per_gene";
mkdir "$output_dir/spliceGrapher_per_gene" unless -e "$output_dir/spliceGrapher_per_gene";
mkdir "$output_dir/plotter_cfg_per_gene" unless -e "$output_dir/plotter_cfg_per_gene";
mkdir "$output_dir/output_pictures" unless -e "$output_dir/output_pictures";

# 1. 得到每个基因的gff3文件
my (%gff3, $gene_id, %gene_region, %gene_scaffold);
open IN, $ARGV[0] or die $!;
while (<IN>) {
    next if m/^#/;
    next if m/^\s*$/;
    if (m/\tgene\t.*ID=(\w+)/) {
        $gene_id = $1;
        $gff3{$gene_id} .= $_;

        @_ = split /\t/;
        $gene_region{$gene_id}{$_[3]} = 1;
        $gene_region{$gene_id}{$_[4]} = 1;
        $gene_scaffold{$gene_id} = $_[0];
    }
    else {
        $gff3{$gene_id} .= $_;
    }
}
close IN;

foreach (keys %gff3) {
    open OUT, ">", "$output_dir/gff3_per_gene/$_.gff3" or die $!;
    print OUT $gff3{$_};
    close OUT;
}

# 2. 通过SpliceGrapher的结果gtf文件，得到有RNA-Seq支持的可变剪接文件
my $cmdString = "gene_model_to_splicegraph.py -m $gtf_file -d $output_dir/spliceGrapher_per_gene";
(system $cmdString) == 0 or die "Failed to execute: $cmdString. $!";
my %gene_id;
open IN, $ARGV[1] or die $!;
while (<IN>) {
    @_ = split /\t/;
    if ($_[8] =~ m/gene_name\s+\"(\w+?)\"/) {
        $gene_id{$1} = 1;
        $gene_region{$1}{$_[3]} = 1;
        $gene_region{$1}{$_[4]} = 1;
    }
}
close IN;

# 3. 得到画图的配置文件
foreach (keys %gff3) {
    my $out;
    if (exists $gene_id{$_}) {
        $out = &get_cfg($_, "yes");
    }
    else {
        $out = &get_cfg($_, "no");
    }
        open OUT, ">", "$output_dir/plotter_cfg_per_gene/$_.cfg" or die $!;
        print OUT $out;
}

# 4. 得到每个基因的画图命令
open OUT, ">", "$output_dir/command.plotter.list";
foreach (keys %gff3) {
    print OUT "plotter.py $output_dir/plotter_cfg_per_gene/$_.cfg\n";
}
close OUT;

# 5. 得到每个基因的sam文件
open OUT, ">", "$output_dir/command.samtools_view.list";
foreach (keys %gene_region) {
    my @region = sort {$a <=> $b} keys %{$gene_region{$_}};
    print OUT "samtools view $bam_file $gene_scaffold{$_}:$region[0]-$region[-1] > $output_dir/sam_per_gene/$_.sam\n";
}
close OUT;
$cmdString = "ParaFly -c $output_dir/command.samtools_view.list -CPU 8";
(system $cmdString) == 0 or die "Failed to execute: $cmdString. $!";

# 6. 运行画图命令
$cmdString = "ParaFly -c $output_dir/command.plotter.list -CPU $ARGV[4]";
(system $cmdString) == 0 or die "Failed to execute: $cmdString. $!";

sub get_cfg {
    my $return;
    my ($gene_id, $yes_or_no) = @_;
    my $GENE_ID = uc($gene_id);
    if ($yes_or_no eq "yes") {
        $return = "[SinglePlotConfig]
legend         = True
output_file    = $output_dir/output_pictures/$gene_id.pdf
shrink_introns = True
width          = 12.0
height         = 8.0

[GeneModelGraph]
plot_type      = gene
gene_name      = $gene_id
relative_size  = 8.0
source_file    = $output_dir/gff3_per_gene/$gene_id.gff3
file_format    = gene_model
title_string   = Gene Model for $gene_id

[SpliceGrapher]
plot_type      = splice_graph
gene_name      = $GENE_ID
relative_size  = 8.0
source_file    = $output_dir/spliceGrapher_per_gene/$GENE_ID.gff
title_string   = SpliceGrapher Prediction for $gene_id

[Junctions]
plot_type      = junctions
labels         = True
relative_size  = 5.0
source_file    = $output_dir/sam_per_gene/$gene_id.sam
title_string   = Splice Junctions of $gene_id

[Reads]
plot_type      = read_depth
relative_size  = 8.0
source_file    = $output_dir/sam_per_gene/$gene_id.sam
title_string   = Read Coverage of $gene_id\n";
    }
    else {
        $return = "[SinglePlotConfig]
legend         = True
output_file    = $output_dir/output_pictures/$gene_id.pdf
shrink_introns = True
width          = 12.0
height         = 8.0

[GeneModelGraph]
plot_type      = gene
gene_name      = $gene_id
relative_size  = 8.0
source_file    = $output_dir/gff3_per_gene/$gene_id.gff3
file_format    = gene_model
title_string   = Gene Model for $gene_id

[Junctions]
plot_type      = junctions
labels         = True
relative_size  = 5.0
source_file    = $output_dir/sam_per_gene/$gene_id.sam
title_string   = Splice Junctions of $gene_id

[Reads]
plot_type      = read_depth
relative_size  = 8.0
source_file    = $output_dir/sam_per_gene/$gene_id.sam
title_string   = Read Coverage of $gene_id\n";
    }
    return $return;
}
