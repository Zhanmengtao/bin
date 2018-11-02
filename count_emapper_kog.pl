#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
     perl $0 emapper.code > emapper_KOG_count.txt

USAGE
     if (@ARGV==0) { die $usage }

     open IN, '<', $ARGV[0] or die $!;
     $_= <IN>;
     $_ =~ s/\s//g;
     my $a=($_=~tr/A//);
     my $b=($_=~tr/B//);
     my $c=($_=~tr/C//);
     my $d=($_=~tr/D//);
     my $e=($_=~tr/E//);
     my $f=($_=~tr/F//);
     my $g=($_=~tr/G//);
     my $h=($_=~tr/H//);
     my $i=($_=~tr/I//);
     my $j=($_=~tr/J//);
     my $k=($_=~tr/K//);
     my $l=($_=~tr/L//);
     my $m=($_=~tr/M//);
     my $n=($_=~tr/N//);
     my $o=($_=~tr/O//);
     my $p=($_=~tr/P//);
     my $q=($_=~tr/Q//);
     my $r=($_=~tr/R//);
     my $s=($_=~tr/S//);
     my $t=($_=~tr/T//);
     my $u=($_=~tr/U//);
     my $v=($_=~tr/V//);
     my $w=($_=~tr/W//);
     my $y=($_=~tr/Y//);
     my $z=($_=~tr/Z//);
     print "Code\tFunctiomynal-Categories\tGene-Number\nA\tRNA processing and modification\t$a\nB\tChromatin structure and dynamics\t$b\nC\tEnergy production and conversion\t$c\nD\tCell cycle control, cell division, chromosome partitioning\t$d\nE\tAmino acid transport and metabolism\t$e\nF\tNucleotide transport and metabolism\t$f\nG\tCarbohydrate transport and metabolism\t$g\nH\tCoenzyme transport and metabolism\t$h\nI\tLipid transport and metabolism\t$i\nJ\tTranslation, ribosomal structure and biogenesis\t$j\nK\tTranscription\t$k\nL\tReplication, recombination and repair\t$l\nM\tCell wall/membrane/envelope biogenesis\t$m\nN\tCell motility\t$n\nO\tPosttranslational modification, protein turnover, chaperones\t$o\nP\tInorganic ion transport and metabolism\t$p\nQ\tSecondary metabolites biosynthesis, transport and catabolism\t$q\nR\tGeneral function prediction only\t$r\nS\tFunction unknown\t$s\nT\tSignal transduction mechanisms\t$t\nU\tIntracellular trafficking, secretion, and vesicular transport\t$u\nV\tDefense mechanisms\t$v\nW\tExtracellular structures\t$w\nY\tNuclear structure\t$y\nZ\tCytoskeleton\t$z";
close IN;
exit;
