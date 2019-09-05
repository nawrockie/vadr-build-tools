#!/usr/bin/env perl
use LWP::Simple;                
use JSON qw( decode_json );     
use XML::LibXML;

use strict;

my $usage = "perl lookup-nt-acc.pl <info file with coded_by value for all protein accessions>";

my $name_str = "";
if(scalar(@ARGV) != 1) { die $usage; }

my $info_file = ($ARGV[0]);

my $chunksize = 50;
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
  if($coded_by_line =~ /coded_by:complement\(([^\:]+)\:(\S+)\)/) { 
    ($source, $coords) = ($1, $2);
    $strand = "-";
  }
  elsif($coded_by_line =~ /coded_by:([^\:]+)\:(\S+)/) { 
    ($source, $coords) = ($1, $2);
    $strand = "+";
  }
  else { 
    die "ERROR unable to parse line $line\n"; 
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
  #exit 0;

  my $xml = XML::LibXML->load_xml(string => $xml_string);

  foreach my $gbseq ($xml->findnodes('//GBSeq')) { 
    my $accver = $gbseq->findvalue('./GBSeq_accession-version');
    if(! defined $coded_by_HAR->{$accver}) { die "ERROR no coded_by for $accver"; }
    # printf("$accver $coded_by\n");
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

  return;
}
