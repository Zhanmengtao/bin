#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 evm_out.gff3 evidence.gff3 Min_coverage_ratio genome.fasta hmmscan_Evalue hmmscan_coverage_ratio CPUs > evm.filter.gff3

For example:
    perl $0 evm_out.gff3 evidence.gff3 0.3 genome.fasta 1e-6 0.3 160 > evm.filter.gff3

    程序用于对EVM得到的gene models进行过滤。示例中的过滤方法为：
    根据转录子或蛋白质序列对基因组的比对信息，得到基因组上可能有基因存在的regions坐标，做成一个数据库。
    然后对EVM得到的gene models和数据库进行比较，若cds的覆盖度 < 0.3，则将该基因翻译成蛋白质序列。
    对蛋白质序列进行Pfam search，设置的hmmscan阈值微 e_value 优于 1e-6，覆盖率 >= 0.3。若Pfam没有比对结果，则过滤该gene model。
    使用hmmscan进行Pfam search的时候，进行了并行化运算，使用了160个并行化运行。默认下使用了/opt/biosoft/hmmer-3.1b2/Pfam-AB.hmm数据库。若需要修改数据库，则需要修改程序内容。
    最后，输出过滤后的gene models。

USAGE
if (@ARGV==0){die $usage}

my $Pfam_db = "/opt/biosoft/hmmer-3.1b2/Pfam-AB.hmm";

open IN, $ARGV[0] or die $!;
$/ = "\n\n";
my (%gene_model, %gene_cds, %cds_length);
while (<IN>) {
    my ($seq_id, $gene_id) = ($1, $2) if m/(\S+?)\tEVM\tgene\t.*?\tID=([^;]+)/;
    $gene_model{$gene_id} = $_;

    my @cds = grep /\tCDS\t/, (split /\n/, $_);
    foreach (@cds) {
        @_ = split /\t/;
        $gene_cds{$gene_id}{"$seq_id\t$_[6]\t$_[3]\t$_[4]\t$_[7]"} = 1;
        $cds_length{$gene_id} += ($_[4] - $_[3] + 1);
    }
}
close IN;

$/ = "\n";
open IN, $ARGV[1] or die $!;
my (%database_collect, %database);
while (<IN>) {
    next if m/^#/;
    next if m/^\s*$/;

    @_ = split /\t/;
    $database_collect{$_[0]}{$_[6]}{"$_[3]\t$_[4]"} = 1;
}
close IN;

foreach my $seq_id (keys %database_collect) {
    foreach my $strand (keys %{$database_collect{$seq_id}}) {
        my @region = sort {$a <=> $b} keys %{$database_collect{$seq_id}{$strand}};
        my ($start, $end) = split /\t/, (shift @region);
        my %region;
        foreach (@region) {
            @_ = split /\t/;
            if ($_[0] > $end) {
                $region{"$start\t$end"} = 1; 
                ($start, $end) = @_;
            }
            elsif ($_[1] > $end) {
                $end = $_[1];
            }
        }
        $region{"$start\t$end"} = 1;

        foreach (keys %region) {
            @_ = split /\t/;
            my $index = int($_[0] / 1000 + 1) * 1000;
            $database{$seq_id}{$strand}{$index}{$_} = 1;
        }
    }
}

my (%for_pfam_search, %validate, $total_gene_number);
foreach my $gene_id (sort keys %gene_model) {
    $total_gene_number ++;
    my $cds_length = $cds_length{$gene_id};
    my @gene_cds = keys %{$gene_cds{$gene_id}};

    my $coverage_ratio = &gene_coverage_ratio($cds_length, @gene_cds);
    #print "$gene_id\t$coverage_ratio\n";
    if ($coverage_ratio >= $ARGV[2]) {
        $validate{$gene_id} = 1;
    }
    else {
        $for_pfam_search{$gene_id} = 1;
    }
}
print STDERR "$total_gene_number genes in total\n";

open IN, $ARGV[3] or die $!;
my (%seq, $seq_id);
while (<IN>) {
    chomp;
    if (/^>(\S+)/) { $seq_id = $1; }
    else { $seq{$seq_id} .= $_; }
}
close IN;

open FASTA, ">", "for_pfam_search.fasta" or die $!;
open COMMAND, ">", "command.hmmscan.list" or die $!;
mkdir "hmmscan.tmp" unless -e "hmmscan.tmp";
my $pfam_search_gene_number;
foreach my $gene_id (sort keys %for_pfam_search) {
    $pfam_search_gene_number ++;
    my @gene_cds = keys %{$gene_cds{$gene_id}};
    my ($seqID, $strand) = split /\t/, $gene_cds[0];
    my %sort;
    foreach (@gene_cds) { @_ = split /\t/; $sort{$_} = $_[2]; }
    my ($cds_seq, $frame);
    
    if ($strand eq "+") {
        @gene_cds = sort {$sort{$a} <=> $sort{$b}} @gene_cds;
        @_ = split /\t/, $gene_cds[0];
        $frame = $_[4];
        foreach (@gene_cds) {
            @_ = split /\t/;
            my $start = $_[2] - 1;
            my $len = $_[3] - $start;
            $cds_seq .= substr($seq{$seqID}, $start, $len);
        }
    }
    elsif ($strand eq "-") {
        @gene_cds = sort {$sort{$b} <=> $sort{$a}} @gene_cds;
        @_ = split /\t/, $gene_cds[0];
        @gene_cds = sort {$sort{$a} <=> $sort{$b}} @gene_cds;
        $frame = $_[4];
        foreach (@gene_cds) {
            @_ = split /\t/;
            my $start = $_[2] - 1;
            my $len = $_[3] - $start;
            $cds_seq .= substr($seq{$seqID}, $start, $len);
        }
        $cds_seq = reverse $cds_seq;
        $cds_seq =~ tr/ATCGatcg/TAGCtagc/;
    }

    $cds_seq =~ s/^\w{$frame}//;
    my $pep_seq = &cds2pep($cds_seq, $gene_id);
    $pep_seq =~ s/\*$//;
    print FASTA ">$gene_id\n$pep_seq\n";

    open OUT, ">", "hmmscan.tmp/$gene_id.fasta" or die $!;
    print OUT ">$gene_id\n$pep_seq\n";
    close OUT;
    print COMMAND "hmmscan -o hmmscan.tmp/$gene_id.txt --cpu 1 -E $ARGV[4] --domE $ARGV[4] --tblout hmmscan.tmp/$gene_id.tbl --domtblout hmmscan.tmp/$gene_id.domtbl $Pfam_db hmmscan.tmp/$gene_id.fasta\n";
}
close FASTA;
close COMMAND;

print STDERR "$pfam_search_gene_number genes were prepared for Pfam search\n";

my $cmdString = "ParaFly -c command.hmmscan.list -CPU $ARGV[6] &> command.hmmscan.log";
system ($cmdString) == 0 or die "Failed to execute: $cmdString\n$!";

my (%num_pfam_ok, $num_pfam_ok);
foreach my $gene_id (keys %for_pfam_search) {
    open IN, "hmmscan.tmp/$gene_id.domtbl" or die $!;
    while (<IN>) {
        next if m/^#/;
        @_ = split /\s+/;
        my $ratio1 = abs($_[15] - $_[16]) / $_[2];
        my $ratio2 = abs($_[17] - $_[18]) / $_[5];
        if ($ratio1 >= $ARGV[5] or $ratio2 >= $ARGV[5]) {
            $validate{$gene_id} = 1;
            $num_pfam_ok{$gene_id} = 1;
        }
    }
    close IN;
}
$num_pfam_ok = keys %num_pfam_ok;
print STDERR "$num_pfam_ok genes were passed through Pfam Search\n";

my $validate_gene_number;
foreach (sort keys %validate) {
    $validate_gene_number ++;
    print $gene_model{$_};
}

my $filterd_gene_number = $pfam_search_gene_number - $num_pfam_ok;
print STDERR "Summary: $filterd_gene_number genes were filtered, $validate_gene_number in the final result\n\n";

sub gene_coverage_ratio {
    my $cds_length = shift;
    my $overlap_length = 0;
    my @cds = @_;
    foreach (@cds) {
        my ($seq_id, $strand, $start, $end) = split /\t/;
        #print "$seq_id, $strand, $start, $end\n";
        my $index1 = int($start / 1000 + 1);
        my $index2 = int($end / 1000 + 1);
        my @region;
        foreach ($index1..$index2) {
            my $index = $_ * 1000;
            #print "$seq_id\t$index\n";
            push @region, keys %{$database{$seq_id}{$strand}{$index}};
        }

        foreach (@region) {
            @_ = split /\t/;
            if ($_[0] <= $end && $_[1] >= $start) {
                if ($_[0] <= $start) {
                    if ($_[1] <= $end) {
                        $overlap_length += ($_[1] - $start + 1);
                    }
                    else {
                        $overlap_length += ($end - $start + 1);
                    }
                }
                else {
                    if ($_[1] <= $end) {
                        $overlap_length += ($_[1] - $_[0] + 1);
                    }
                    else {
                        $overlap_length += ($end - $_[0] + 1);
                    }
                }
            }
        }
    }
    my $ratio = $overlap_length / $cds_length;
    return $ratio;
}

sub cds2pep {
    my %cds2pep = (
        "TTT" => "F",
        "TTC" => "F",
        "TTA" => "L",
        "TTG" => "L",
        "TCT" => "S",
        "TCC" => "S",
        "TCA" => "S",
        "TCG" => "S",
        "TAT" => "Y",
        "TAC" => "Y",
        "TAA" => "*",
        "TAG" => "*",
        "TGT" => "C",
        "TGC" => "C",
        "TGA" => "*",
        "TGG" => "W",
        "CTT" => "L",
        "CTC" => "L",
        "CTA" => "L",
        "CTG" => "L",
        "CCT" => "P",
        "CCC" => "P",
        "CCA" => "P",
        "CCG" => "P",
        "CAT" => "H",
        "CAC" => "H",
        "CAA" => "Q",
        "CAG" => "Q",
        "CGT" => "R",
        "CGC" => "R",
        "CGA" => "R",
        "CGG" => "R",
        "ATT" => "I",
        "ATC" => "I",
        "ATA" => "I",
        "ATG" => "M",
        "ACT" => "T",
        "ACC" => "T",
        "ACA" => "T",
        "ACG" => "T",
        "AAT" => "N",
        "AAC" => "N",
        "AAA" => "K",
        "AAG" => "K",
        "AGT" => "S",
        "AGC" => "S",
        "AGA" => "R",
        "AGG" => "R",
        "GTT" => "V",
        "GTC" => "V",
        "GTA" => "V",
        "GTG" => "V",
        "GCT" => "A",
        "GCC" => "A",
        "GCA" => "A",
        "GCG" => "A",
        "GAT" => "D",
        "GAC" => "D",
        "GAA" => "E",
        "GAG" => "E",
        "GGT" => "G",
        "GGC" => "G",
        "GGA" => "G",
        "GGG" => "G",
    );
    my $seq = shift @_;
    my $gene = shift @_;
    my $pep;
    while ((length $seq) >= 3) {
        $seq =~ s/(\w{3})//;
        $pep .= $cds2pep{$1};
    }
    #print STDERR "Warning: CDS length of $gene is not multiple of 3\n" if (length $seq) > 0;
    #print STDERR "Warning: Stop Codon appear in the middle of $gene\n" if $pep =~ m/\*\w/;
    return $pep;
}
