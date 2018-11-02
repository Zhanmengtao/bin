#!/usr/bin/perl -w

use strict;
use LWP::UserAgent;
use Getopt::Std;

###############################################################################
# URL to the Batch CD-Search server
###############################################################################

my $bwrpsb = "http://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi";

###############################################################################
# read list of queries and parameters supplied; queries specified in list piped
# from stdin
###############################################################################

my @queries = <STDIN>;
my $havequery = 0;

###############################################################################
# do some sort of validation and exit if only invalid lines found
###############################################################################

foreach my $line (@queries) {
  if ($line =~ /[a-zA-Z0-9_]+/) {
    $havequery = 1;
  }
}
if ($havequery == 0) {
  die "No valid queries!\n";
}

###############################################################################
# set default values
###############################################################################
my $cdsid = "";
my $cddefl = "false";
my $qdefl = "false";
my $smode = "auto";
my $useid1 = "true";

# Integer, specifies the maximum number of hits to return for each protein in the query list. This parameter only applies if the search mode (smode) is live. If the search mode is set to precalculated or automatic, then the default maxhit number of 500 is applied.
my $maxhit = 250;

# this parameter specifies whether the program will filter out compositionally biased regions from the query sequences.
my $filter = "true";

# Specify the name of the database. Allowable values include: "cdd," "pfam," "smart," "tigrfam," "cog," "kog". This parameter only applies if the search mode (smode) is live. If the search mode is set to precalculated or automatic, then the default CDD database is searched.
my $db = "cdd";

my $evalue = 0.01;

# Specify the data mode desired in the output. Allowable values are:
# "rep" (highest scoring hit, for each region of the query sequence, as shown in the concise results), or
# "std" (best-scoring hit from each source database, for each region of the query sequence, as shown in the standard results), or
# "full" (the complete set of hits in the full results).
# (The value "all" is still allowed and previously was used to display the full result, but now displays the standard result, since Feb. 12, 2014, when the standard mode became available as a new display option.)
my $dmode = "rep";

my $clonly = "false";

# Specify the data type (target data) desired in the output. Allowable values are: "hits" (domain hits), "aligns" (alignment details), or "feats" (features). 
my $tdata = "hits";

###############################################################################
# deal with command line parameters, change default settings if necessary
###############################################################################

our($opt_d, $opt_e, $opt_F, $opt_b, $opt_t, $opt_s, $opt_a, $opt_q);

getopts('d:e:F:b:t:s:a:q');

if ($opt_d) {
  $db = $opt_d;
  print "Databast option set to: $db\n";
}
if ($opt_e) {
  $evalue = $opt_e;
  print "Evalue option set to: $evalue\n";
}
if ($opt_F) {
  if ($opt_F eq "F") {
    $filter = "false"
  } else {
    $filter = "true";
  }
  print "Filter option set to: $filter\n";
}
if ($opt_b) {
  $maxhit = $opt_b;
  print "Maxhit option set to: $maxhit\n";
}
if ($opt_t) {
  $tdata = $opt_t;
  print "Target data option set to: $tdata\n";
}
if ($opt_s) {
  $clonly = "true";
  print "Superfamilies only will be reported\n";
}
if ($opt_a) {
  $dmode = "all";
  print "All hits will be reported\n";
}
if ($opt_q) {
  $qdefl = "true";
  print "Query deflines will be reported\n";
}

###############################################################################
# submitting the search
###############################################################################
my $rid;
{
  my $browser = LWP::UserAgent->new;
  my $response = $browser->post(
    $bwrpsb,
    [
      'useid1' => $useid1,
      'maxhit' => $maxhit,
      'filter' => $filter,
      'db'     => $db,
      'evalue' => $evalue,
      'cddefl' => $cddefl,
      'qdefl'  => $qdefl,
      'dmode'  => $dmode,
      'clonly' => $clonly,
      'tdata'  => "hits",
      ( map {; queries => $_ } @queries )
    ],
  );
  die "Error: ", $response->status_line
    unless $response->is_success;

  if($response->content =~ /^#cdsid\s+([a-zA-Z0-9-]+)/m) {
    $rid =$1;
    print "Search with Request-ID $rid started.\n";
  } else {
    die "Submitting the search failed,\n can't make sense of response: $response->content\n";
  }
}
###############################################################################
# checking for completion, wait 5 seconds between checks
###############################################################################

$|++;
my $done = 0;
my $status = -1;
while ($done == 0) {
  sleep(5);
  my $browser = LWP::UserAgent->new;
  my $response = $browser->post(
    $bwrpsb,
    [
      'tdata' => "hits",
      'cdsid' => $rid
    ],
  );
  die "Error: ", $response->status_line
    unless $response->is_success;

  if ($response->content =~ /^#status\s+([\d])/m) {
    $status = $1;
    if ($status == 0) {
      $done = 1;
      print "Search has been completed, retrieving results ..\n";
    } elsif ($status == 3) {
      print ".";
    } elsif ($status == 1) {
      die "Invalid request ID\n";
    } elsif ($status == 2) {
      die "Invalid input - missing query information or search ID\n";
    } elsif ($status == 4) {
      die "Queue Manager Service error\n";
    } elsif ($status == 5) {
      die "Data corrupted or no longer available\n";
    }
  } else {
    die "Checking search status failed,\ncan't make sense of response: $response->content\n";
  }

}
print "===============================================================================\n\n";

###############################################################################
# retrieve and display results
###############################################################################
{
  my $browser = LWP::UserAgent->new;
  my $response = $browser->post(
    $bwrpsb,
    [
        'tdata'  => $tdata,
        'cddefl' => $cddefl,
        'qdefl'  => $qdefl,
        'dmode'  => $dmode,
        'clonly' => $clonly,
        'cdsid'  => $rid
    ],
  );
  die "Error: ", $response->status_line
    unless $response->is_success;

  print $response->content,"\n";
}
