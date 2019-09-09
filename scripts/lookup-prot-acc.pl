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
  my $xml = undef;
  my $xml_string = get($genbank_url);
  my $xml_valid = 0;
  if(defined $xml_string) { 
      # to save memory, remove sequence info from the xml_string since we don't need it
      # remove <GBSeq_sequence> lines
      $xml_string =~ s/[^\n]+\<GBSeq\_sequence\>\w+\<\/GBSeq\_sequence\>\n//g;
      # remove <GBQualifier>\n<GBQualifer_name>translation\nGBQualifier_value\n<\GBQualifier> sets of 4 lines
      $xml_string =~ s/[^\n]+\<GBQualifier\>\n[^\n]+\<GBQualifier\_name\>translation\<\/GBQualifier\_name\>\n[^\n]+\<GBQualifier\_value\>\w+\<\/GBQualifier\_value\>\n[^\n]+\<\/GBQualifier\>\n//g;
      $xml = eval { XML::LibXML->load_xml(string => $xml_string); };
      if($@) { $xml_valid = 0; }
      else   { $xml_valid = 1; }
  }
  if(! $xml_valid) { 
   # the get() command either failed (returned undef) or
    # returned an invalid xml string, either way we
    # wait a few seconds ($nseconds) and try again (up to
    # $nattempts) times BUT we only do this if the ID doesn't look 
    # like a RNAcentral ids. If it does, we do not do more attempts.
    my $attempt_ctr = 1;
    while((! $xml_valid) && ($attempt_ctr < $nattempts)) { 
      sleep($nseconds);
      # printf("Retrying to fetch for $name\n");
      $xml_string = get($genbank_url);
      if(defined $xml_string) { 
        # to save memory, remove sequence info from the xml_string since we don't need it
        # remove <GBSeq_sequence> lines
        $xml_string =~ s/[^\n]+\<GBSeq\_sequence\>\w+\<\/GBSeq\_sequence\>\n//g;
        # remove <GBQualifier>\n<GBQualifer_name>translation\nGBQualifier_value\n<\GBQualifier> sets of 4 lines
        $xml_string =~ s/[^\n]+\<GBQualifier\>\n[^\n]+\<GBQualifier\_name\>translation\<\/GBQualifier\_name\>\n[^\n]+\<GBQualifier\_value\>\w+\<\/GBQualifier\_value\>\n[^\n]+\<\/GBQualifier\>\n//g;
        $xml = eval { XML::LibXML->load_xml(string => $xml_string); };
        if($@) { $xml_valid = 0; }
        else   { $xml_valid = 1; }
      }
      $attempt_ctr++;
    }
    if(($attempt_ctr >= $nattempts) && (! $xml_valid)) { 
      die "ERROR trying to fetch sequence data from genbank, reached maximum allowed number of attempts ($attempt_ctr)"; 
    }
  }
  else { 
    # if we get here: we know that $xml_string is defined and valid
    # and $xml is ready for parsing  }
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
  } # end of 'else' entered if $xml_string is defined

  return;
}
