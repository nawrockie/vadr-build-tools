#!/usr/bin/env perl
# EPN, Tue Jan  1 06:39:31 2019
use warnings;
use strict;
use Getopt::Long;

my $usage;
$usage  = "perl fasta-ncbi-idfetch-name-long-to-short.pl <fasta file>\n";

while(my $line = <>) { 
  if($line =~ /^>(\S+)(\s*.*)/) { 
    chomp $line;
    my ($name, $desc) = ($1, $2);
    if($name =~ /^gi\|\d+\|\S+\|(\S+\.\d+)\|.*/) { 
      $name = $1;
    }
    printf ">%s%s\n", $name, $desc;
  }
  else { 
    print $line;
  }
}
exit 0;
