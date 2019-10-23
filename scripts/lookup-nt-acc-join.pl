#!/usr/bin/env perl
use LWP::Simple;                
use JSON qw( decode_json );     
use XML::LibXML;

require "vadr.pm";

use strict;

my $usage = "perl lookup-nt-acc.pl <info file with coded_by value for all protein accessions>";

my $name_str = "";
if(scalar(@ARGV) != 1) { die $usage; }

my $info_file = ($ARGV[0]);

my $chunksize = 10;
#my $chunksize = 10000;
my $nattempts = 3;
my @source_A = ();

open(IN, $info_file) || die "ERROR unable to open $info_file";

my $ct = 0;
my $name_str = "";
my %codon_start_H = ();
my %coded_by_HA = ();
my ($source, $coords);
my @codon_start_key_A = ();
my %codon_start_key_exists_H = ();

my %line_H = ();
while(my $line = <IN>) { 
  chomp $line;
  #accver:YP_007626734.1 coded_by:NC_020752.1:5328..6872 transl_table:2 taxonomy:Eukaryota;Metazoa;Chordata;Craniata;Vertebrata;Euteleostomi;Mammalia;Eutheria;Laurasiatheria;Cetartiodactyla;Ruminantia;Pecora;Bovidae;Bovinae;Tragelaphus  
  my @el_A = split(/\s+/, $line);
  my $coded_by_line = $el_A[1];
  my $strand = "+";
  my $coords = undef;

  # first extract the source accession
  if($coded_by_line =~ /coded_by:join\(([^\:]+)\:/) { 
    #coded_by:join(NC_044902.1:16688..17577,NC_044902.1:18432..19113)
    $source = $1;
  }
  elsif($coded_by_line =~ /coded_by:complement\(join\(([^\:]+)\:/) { 
    $source = $1;
  }
  elsif($coded_by_line =~ /coded_by:complement\(([^\:]+)\:/) { 
    $source = $1;
  }
  elsif($coded_by_line =~ /coded_by:([^\:]+)\:/) { 
    $source = $1;
  }
  else { 
    die "ERROR couldn't figure out accession in coded_by value on line $line"; 
  }

  my $no_source_coded_by_line = $coded_by_line;
  $no_source_coded_by_line =~ s/\Q$source\E\://g; 

  if($no_source_coded_by_line =~ /coded_by:complement\((\S+)\)/) { 
    ($coords) = ($1);
    $strand = "-";
  }
  elsif($no_source_coded_by_line =~ /coded_by:(\S+)/) { 
    ($coords) = ($1);
    $strand = "+";
  }
  else { 
    die "ERROR unable to parse no_source_coded_by_line: $no_source_coded_by_line in line\n$line\n"; 
  }
  if($coords =~ /\<?(\d+)\.\.\>?(\d+)/) { 
    my($start, $stop) = ($1, $2);
    if($strand eq "-") { 
      my $tmp = $start;
      $start = $stop;
      $stop = $tmp;
    }
  }
  else { 
    die "ERROR unable to parse (2) line $line\n";
  }
  if($name_str ne "") { $name_str .= ","; }
  $name_str .= $source;
  $ct++;
  if(! defined $coords) { die "ERROR coords undef for source $source\n"; }
  my $coded_by_value = "";
  if($strand eq "-") { 
    $coded_by_value = "complement($coords)"; 
  }
  else { 
    $coded_by_value = $coords; 
  }
  if(! defined $coded_by_HA{$source}) { 
    @{$coded_by_HA{$source}} = ();
  }
  push(@{$coded_by_HA{$source}}, $coded_by_value);
  my $codon_start_key = $source . ":::" . $coded_by_value;
  $line_H{$codon_start_key} = $line;
  if(defined $codon_start_key_exists_H{$codon_start_key}) { 
    die "ERROR two instances for codon start key: $codon_start_key"; 
  }
  push(@codon_start_key_A, $codon_start_key);
  if($ct == $chunksize) { 
    fetch_chunk($name_str, "nuccore", \%codon_start_H, \%coded_by_HA);
    $ct = 0;
    $name_str = "";
  }
}
if($ct > 0) { 
  fetch_chunk($name_str, "nuccore", \%codon_start_H, \%coded_by_HA);
}

foreach my $codon_start_key (@codon_start_key_A) { 
  if(! defined $codon_start_H{$codon_start_key}) { die "ERROR no codon start for $source"; }
  print $line_H{$codon_start_key} . " codon_start=" . $codon_start_H{$codon_start_key} . "\n";
}
exit 0;

#######################################
sub fetch_chunk { 
  my ($name_str, $db, $codon_start_HR, $coded_by_HAR) = (@_);

  my $nattempts = 3;
  my $nseconds = 5;

  my $genbank_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=$db&retmode=xml&id=" . $name_str;

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
      #printf("parsing $accver\n");
      if(! defined $coded_by_HAR->{$accver}) { die "ERROR no coded_by for $accver"; }
      my $found_location = 0;
      
      foreach my $coded_by (@{$coded_by_HAR->{$accver}}) { 
        my $codon_start_key = $accver . ":::" . $coded_by; 
        foreach my $feature ($gbseq->findnodes('//GBFeature')) { 
          my $feature_key = $feature->findvalue('./GBFeature_key');
          if($feature_key eq "CDS") { 
            my $location    = $feature->findvalue('./GBFeature_location');
            if($location eq $coded_by) { 
              $found_location = 1;
              my $codon_start = $feature->findvalue('./GBFeature_quals/GBQualifier/GBQualifier_name[text()="codon_start"]/following-sibling::GBQualifier_value');
              #if(defined $codon_start) { printf("DEFINED $accver $location $codon_start\n"); }
              if(! defined $codon_start) { 
                #printf("UNDEFINED $accver $location\n");
                $codon_start = 1; 
              }
              $codon_start_H{$codon_start_key} = $codon_start;
            }
          }
        }
        if(! $found_location) { 
          printf STDERR "ERROR did not find CDS feature in $accver with location $coded_by\n"; 
          $codon_start_H{$codon_start_key} = "?";
        }
      }
    }
  } # end of 'else' entered if $xml_string is defined

  return;
}
