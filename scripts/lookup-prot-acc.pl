#!/usr/bin/env perl
use LWP::Simple;                
use JSON qw( decode_json );     
use XML::LibXML;

use strict;

my $usage = "perl lookup-prot-acc.pl <list of protein accession.versions>";

my $name_str = "";
if(scalar(@ARGV) != 1) { die $usage; }

my $infile = ($ARGV[0]);

my $chunksize = 50;
my $nattempts = 3;

open(IN, $infile) || die "ERROR unable to open $infile";

my $ct = 0;
my $name_str = "";

while(my $line = <IN>) { 
  chomp $line;
  if(($line !~ m/^\#/) && ($line =~ m/\w/)) { 
    if($name_str ne "") { $name_str .= ","; }
    $name_str .= $line;
    $ct++;
    if($ct == $chunksize) { 
      fetch_chunk($name_str);
      $ct = 0;
      $name_str = "";
    }
  }
}
if($ct > 0) { 
  fetch_chunk($name_str);
}
exit 0;

#######################################
sub fetch_chunk { 
  my $name_str = $_[0];

  my $nattempts = 3;
  my $nseconds = 5;

  my $genbank_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&retmode=xml&id=" . $name_str;
  my $xml_string = get($genbank_url);

  if(! defined $xml_string) { 
    # if NCBI is being hit by a bunch of requests, the get() command
    # may fail in that $got_url may be undefined. If that happens we
    # wait a few seconds ($nseconds) and try again (up to
    # $nattempts) times.
    my $attempt_ctr = 1;
    while((! defined $xml_string) && ($attempt_ctr < $nattempts)) { 
      sleep($nseconds);
      $xml_string = get($genbank_url);
      $attempt_ctr++;
    }
    if(($attempt_ctr >= $nattempts) && (! defined $xml_string)) { 
      die "ERROR trying to fetch sequence data from genbank, reached maximum allowed number of attempts ($attempt_ctr)"; 
    }
  }
  
  #printf("$xml_string\n");
  
  my $xml = XML::LibXML->load_xml(string => $xml_string);

  foreach my $gbseq ($xml->findnodes('//GBSeq')) { 
    my $accver = $gbseq->findvalue('./GBSeq_accession-version');
    
    my $coded_by = $gbseq->findvalue('./GBSeq_feature-table/GBFeature/GBFeature_quals/GBQualifier/GBQualifier_name[text()="coded_by"]/following-sibling::GBQualifier_value');
    if(! defined $coded_by) { 
      die "ERROR did not read coded_by for $accver";
    }
    
    my $transl_table = $gbseq->findvalue('./GBSeq_feature-table/GBFeature/GBFeature_quals/GBQualifier/GBQualifier_name[text()="transl_table"]/following-sibling::GBQualifier_value');
    
    my $taxonomy = $gbseq->findvalue('./GBSeq_taxonomy');
    $taxonomy =~ s/\s+//g;
    printf("accver:$accver coded_by:$coded_by transl_table:%s taxonomy:%s\n", (defined $transl_table) ? $transl_table : "undefined", (defined $taxonomy) ? $taxonomy : "undefined");
  }

  return;
}
