use strict;

my $usage = "perl split-by-tax.pl <info file> <KEEP string (or NONE)> <comma separated SKIP string(s) or NONE> <tax level> <extra string to add to output file name (NONE for none)> <output root>";
if(scalar(@ARGV) != 6) { die $usage; }
my ($info_file, $keep_str, $skip_str, $tax_level, $extra_string, $file_root) = (@ARGV);

my @skip_A = ();
if($skip_str ne "NONE") { 
  @skip_A = split(",", $skip_str);
}
if(($keep_str ne "NONE") && ($skip_str ne "NONE")) { 
  die "ERROR neither <KEEP string> nor <SKIP string> are NONE";
}
# they can both equal NONE

open(IN, $info_file) || die "ERROR unable to open $info_file";
my %bytax_HA = ();
while(my $line = <IN>) { 
  chomp $line;
  my @el_A = split(/\s+/, $line);
  my $taxline = $el_A[3];
  if($taxline =~ /taxonomy:(\S+)/) { 
    my $tax = $1;
    my $include_this_line = 0;
    if(($keep_str eq "NONE") && ($skip_str eq "NONE")) { 
      $include_this_line = 1;
    }
    if($keep_str ne "NONE") { 
      if($tax =~ /$keep_str/) { $include_this_line = 1; }
    }
    if(scalar(@skip_A) > 0) { 
      $include_this_line = 1;
      foreach my $skip_tok (@skip_A) { 
        if($tax =~ /$skip_tok/) { 
          $include_this_line = 0; 
        }
      }
    }
    # printf("include_this_line: $include_this_line\n");
    if($include_this_line) { 
      my @tax_A = split(";", $tax);
      my $taxkey = (scalar(@tax_A) < $tax_level) ? $tax_A[(scalar(@tax_A)-1)] : $tax_A[($tax_level-1)];
      if(! defined $bytax_HA{$taxkey}) { 
        @{$bytax_HA{$taxkey}} = (); 
      }
      push(@{$bytax_HA{$taxkey}}, $line); 
    }
  }
}
foreach my $taxkey (sort keys %bytax_HA) { 
  my $taxkey_lc = $taxkey;
  $taxkey_lc =~ tr/A-Z/a-z/;
  my $tmp_root = sprintf("%s.%s%s", $file_root, (($extra_string eq "NONE") ? "" : $extra_string), $taxkey_lc);
  my $file_name = $tmp_root . ".info";
  printf("$tmp_root\n");
  open(OUT, ">", $file_name) || die "ERROR unable to open $file_name for writing";
  foreach my $line (@{$bytax_HA{$taxkey}}) { 
    print OUT $line . "\n";
  }
  close OUT;
}

