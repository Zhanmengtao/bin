#!/usr/bin/perl
use strict;
my $usage = <<USAGE;
Usage:
    perl $0

This command will modify 4 files:
    /opt/biosoft/wwwblast/blast.rc
	/opt/biosoft/wwwblast/psiblast.rc
	/opt/biosoft/wwwblast/blast.html
	/opt/biosoft/wwwblast/psiblast.html

USAGE

my (%nucl_db, %prot_db);
foreach (</opt/biosoft/wwwblast/db/*.nsq>) {
	s#/opt/biosoft/wwwblast/db/##;
	s/\.nsq//;
	s/\.\d+$//;
	$nucl_db{$_} = 1;
}
foreach (</opt/biosoft/wwwblast/db/*.psq>) {
	s#/opt/biosoft/wwwblast/db/##;
	s/\.psq//;
	s/\.\d+$//;
	$prot_db{$_} = 1;
}

my $nucl_db = join " ", keys %nucl_db;
my $prot_db = join " ", keys %prot_db;

# /opt/biosoft/wwwblast/blast.rc
open OUT, '>', "/opt/biosoft/wwwblast/blast.rc" or die $!;
print OUT "NumCpuToUse\t4\n\nblastn $nucl_db\nblastp $prot_db\nblastx $prot_db\ntblastn $nucl_db\ntblastx $nucl_db\n";
close OUT;

# /opt/biosoft/wwwblast/psiblast.rc
open OUT, '>', "/opt/biosoft/wwwblast/psiblast.rc" or die $!;
print OUT "NumCpuToUse\t4\n\nblastn $nucl_db\nblastp $prot_db\nblastx $prot_db\ntblastn $nucl_db\ntblastx $nucl_db\n";
close OUT;

# /opt/biosoft/wwwblast/blast.html
my $html_out;
open IN, "/opt/biosoft/wwwblast/blast.html" or die $!;
$html_out = join "", <IN>;
close IN;

my $data_lib;
foreach (sort keys %nucl_db) {
	$data_lib .= "    <option VALUE = \"$_\"> $_\n";
}
foreach (sort keys %prot_db) {
	$data_lib .= "    <option VALUE = \"$_\"> $_\n";
}
$data_lib =~ s/^\n*//;
$data_lib =~ s/\n*$//;
$html_out =~ s#(<select name = \"DATALIB\">).*?(</select>)#$1\n$data_lib\n$2#s;
open OUT, '>', "/opt/biosoft/wwwblast/blast.html" or die $!;
print OUT $html_out;

# /opt/biosoft/wwwblast/psiblast.html
open IN, "/opt/biosoft/wwwblast/psiblast.html" or die $!;
$html_out = join "", <IN>;
close IN;
$html_out =~ s#(<select name = \"DATALIB\">).*?(</select>)#$1\n$data_lib\n$2#s;
open OUT, '>', "/opt/biosoft/wwwblast/psiblast.html" or die $!;
print OUT $html_out;
