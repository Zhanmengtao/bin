#!/usr/bin/perl
use strict;
use File::Basename;

my $usage = <<USAGE;
Usage:
    perl $0 geneModels.gff3 proteins.fasta CPUs Max_identity > geneModels.lowIdentity.gff3

For Example:
    perl $0 best_candidates.gff3 best_candidates.fasta 24 0.70 > best_candidates.lowIdentity.gff3

    1. 程序调用了blast构建蛋白质数据库，然后进行了all vs all blast；需要ncbi blast + 和程序同目录下的 blast.pl 程序支持。
    2. 根据blast结果，若identity高于0.70，则去掉序列长度较短的蛋白质。
    3. 根据需要去除的蛋白信息，保留 best_candidates.gff3 中其余的基因信息。


USAGE
if (@ARGV==0){die $usage}

my $path = dirname($0);
unless ($path =~ m/^[\/~]/) { $path = "../$path"; }

open IN, $ARGV[1] or die $!;
my (%length, $id);
while (<IN>) {
    chomp;
    if (m/>(\S+)/) { $id = $1; }
    else { $length{$id} += length; }
}
close IN;

mkdir "remove_redundant_high_identity_genes.tmp" unless -e "remove_redundant_high_identity_genes.tmp";
my $pwd = `pwd`;
chomp($pwd);
chdir "remove_redundant_high_identity_genes.tmp";

my $cmdString = "makeblastdb -in $pwd/$ARGV[1] -dbtype prot -title protein -parse_seqids -out protein -logfile protein.log";
print STDERR "CMD 1/2: $cmdString\n";
unless (-e "makeblastdb.ok") {
    system ($cmdString) == 0 or die "Failed to execute: $!\n";
    print STDERR "Successfully Make Blast Database of target protein sequences\n\n";
    open OK, ">", "makeblastdb.ok";
    close OK;
}
else {
    print STDERR "Skip CMD 1/2, for file makeblastdb.ok exists\n\n";
}

my $cmdString = "$path/blast.pl blastp protein $pwd/$ARGV[1] 1e-10 $ARGV[2] out 6";
print STDERR "CMD 2/2: $cmdString\n";
unless (-e "blast.ok") {
    system ($cmdString) == 0 or die "Failed to execute: $!\n";
    print STDERR "Successfully blast all to all\n\n";
    open OK, ">", "blast.ok";
    close OK;
}
else {
    print STDERR "Skip CMD 2/2, for file blast.ok exists\n\n";
}

open IN, "out.tab" or die $!;
my %delete;
while (<IN>) {
    @_ = split /\t/;
    next if $_[0] eq $_[1];
    next if $_[2] < ($ARGV[3] * 100);
    if ( !exists $delete{$_[0]} && !exists $delete{$_[1]} ) {
       if ($length{$_[0]} >= $length{$_[1]}) { $delete{$_[1]} = 1; }
       elsif ($length{$_[0]} < $length{$_[1]}) { $delete{$_[0]} = 1; }
       else {
           print STDERR "Warning: No such proteins in $pwd/$ARGV[1]\n$_";
       }
   }
}
close IN;

my @redundant = keys %delete;
my $redundant = @redundant;
print STDERR "$redundant removed redundant high identity proteins are:\n";
foreach (@redundant) {
    print STDERR "$_\n";
}

open IN, "$pwd/$ARGV[0]" or die $!;
$/ = "\n\n";
while (<IN>) {
    m/mRNA\t.*?ID=(\S+?);/;
    print unless exists $delete{$1};
}
close IN;
