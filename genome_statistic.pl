#!/usr/bin/perl
use strict;
use File::Basename;
my $usage = <<USAGE;
Usage:
    perl $0 genomeFasta

USAGE
if (@ARGV==0) {die $usage}

foreach ( @ARGV ) {
	my $name = basename( $_ );
	print "\n==>>$name<<==\n\n";
	my ($genome_size, $contig_size, $identifier, %scaffold_length, %scaffold_seq, $N_num, $G_num, $C_num, $T_num, $A_num);
	open GENOME, '<',"$_" or die "Can't open file $_!($!)";
	while ( <GENOME> ) {
		chomp;
		if ( m/^>(\S+)/ ) { $identifier = $1 }
		else {
			my $length = length $_;
			$genome_size += $length;
			$scaffold_length{$identifier} += $length;
			$scaffold_seq{$identifier} .= $_;
			$N_num += ($_ =~ tr/Nn/Nn/);
			$G_num += ($_ =~ tr/Gg/Gg/);
			$C_num += ($_ =~ tr/Cc/Cc/);
			$T_num += ($_ =~ tr/Tt/Tt/);
			$A_num += ($_ =~ tr/Aa/Aa/);
		}
	}
	close GENOME;
	
	my @contig_length;
	foreach (keys %scaffold_seq) {
		my $scaffold_seq = $scaffold_seq{$_};
		#my @contigs = split /N{10,}/, $scaffold_seq;
		my @contigs = split /[Nn]+/, $scaffold_seq;
		foreach (@contigs) {
			my $length = length $_;
			$contig_size += $length;
			push @contig_length, $length
		}
	}
	@contig_length = sort { $b <=> $a } @contig_length;

	my @scaffold_length = values %scaffold_length;
	@scaffold_length = sort { $b <=> $a } @scaffold_length;
	
	my $longest_fragment = $scaffold_length[0];
	my $shortest_fragment = $scaffold_length[$#scaffold_length];
	my $sequence_number = @scaffold_length;
	my $contig_number = @contig_length;
	print "the genome scaffolds number is $sequence_number\nthe genome contigs number is $contig_number\nthe longest length is $longest_fragment\nthe shortest length is $shortest_fragment\n";

	my $rate_of_N = $N_num / $genome_size;
	my $rate_of_GC =  ( $G_num + $C_num ) / ($G_num + $C_num + $T_num + $A_num);
	print "the genome scaffolds size is $genome_size\nthe genome contig size is $contig_size\n";
	print "the rate of N is $rate_of_N\n";
	print "the rate of GC is $rate_of_GC\n";

	my $genome_size50 = $genome_size * 0.5; my $contig_size50 = $contig_size * 0.5;
	my $genome_size90 = $genome_size * 0.9; my $contig_size90 = $contig_size * 0.9;
	my $aclength = 0;
	foreach (@scaffold_length) {
		$aclength += $_;
		if ($aclength >= $genome_size50) {
			print "the scaffold N50 is $_\n";
			last
		}
	}
	$aclength = 0;
	foreach (@contig_length) {
		$aclength += $_;
		if ($aclength >= $contig_size50) {
			print "the contig N50 is $_\n";
			last
		}
	}
	$aclength = 0;
	foreach (@scaffold_length) {
		$aclength += $_;
		if ($aclength >= $genome_size90) {
			print "the scaffold N90 is $_\n";
			last
		}
	}
	$aclength = 0;
	foreach (@contig_length) {
		$aclength += $_;
		if ($aclength >= $contig_size90) {
			print "the contig N90 is $_\n";
			last
		}
	}
	my ($length1000Num,$length2000Num,$length3000Num,$total_length3000,$total_length2000_3000,$total_length1000_2000,$total_length2000,$total_length1000,$length2000_3000Num,$length1000_2000Num);
	foreach (@scaffold_length) {
		if ($_ >= 3000) {
			$length3000Num ++;
			$total_length3000 += $_
		}elsif ($_ >= 2000) {
			$length2000_3000Num ++;
			$total_length2000_3000 += $_
		}elsif ($_ >= 1000) {
			$length1000_2000Num ++;
			$total_length1000_2000 += $_
		}
		$length2000Num = $length3000Num + $length2000_3000Num;
		$length1000Num = $length2000Num + $length1000_2000Num;
		$total_length2000 = $total_length3000 + $total_length2000_3000;
		$total_length1000 = $total_length2000 + $total_length1000_2000;
	}
	print "the number of sequences >= 1kb is $length1000Num\ttotal length is $total_length1000\n";
	print "the number of sequences >= 2kb is $length2000Num\ttotal length is $total_length2000\n";
	print "the number of sequences >= 3kb is $length3000Num\ttotal length is $total_length3000\n";
	print "\n";
}
