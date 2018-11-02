#!/usr/bin/perl
use strict;
use warnings;
use Bio::Tools::Run::RemoteBlast;
my $usage = <<USAGE;
Usage:
    perl $0 blastProgramType blastDB fastaFile

USAGE
if (@ARGV != 3) { die $usage }

my $prog=$ARGV[0];
my $db=$ARGV[1];
my $e_value='1e-3';
my $readmethod='xml';

my @paras = ( '-prog' => $prog,
        '-data' => $db,
        '-expect' => $e_value,
        '-readmethod' => $readmethod );

my $remote_blastxml = Bio::Tools::Run::RemoteBlast->new(@paras);
$remote_blastxml->retrieve_parameter('FORMAT_TYPE', 'XML');

my $v = 1;
my $seqIn = Bio::SeqIO->new(-file=>$ARGV[2], -format => 'fasta' );

my $num = 1;
while (my $query = $seqIn->next_seq()){
	#Blast a sequence against a database:
	
	#Alternatively, you could  pass in a file with many
	#sequences rather than loop through sequence one at a time
	#Remove the loop starting 'while (my $input = $str->next_seq())'
	#and swap the two lines below for an example of that.
	my $report = $remote_blastxml->submit_blast($query);

	print STDERR "\n$num: waiting..." if( $v > 0 );
	$num ++;
	while ( my @rids = $remote_blastxml->each_rid ) {
		foreach my $rid ( @rids ) {
			my $rc = $remote_blastxml->retrieve_blast($rid);
			if( !ref($rc) ) {
				if( $rc < 0 ) {
					$remote_blastxml->remove_rid($rid);
				}
				print STDERR "." if ( $v > 0 );
				sleep 5;
			} else {
				my $result = $rc->next_result();
				#save the output
				my $filename = $result->query_name()."\.out";
				$remote_blastxml->save_output($filename);
				`cat $filename >> blastRemoteResult.xml; rm $filename`;
				$remote_blastxml->remove_rid($rid);
			}
		}
	}
}
