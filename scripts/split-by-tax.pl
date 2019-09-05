my $usage = "perl sort-by-tax.pl <info file> <tax level to split on> <extra string to add to output file name (NONE for none)> <output root>";
if(scalar(@ARGV) != 4) { die $usage; }
my ($info_file, $tax_level, $extra_string, $file_root) = (@ARGV);

open(IN, $info_file) || die "ERROR unable to open $info_file";
while($line = <IN>) { 
  chomp $line;
  my @el_A = split(/\s+/, $line);
  my $taxline = $el_A[3];
  if($taxline =~ /taxonomy:(\S+)/) { 
    $tax = $1;
    @tax_A = split(";", $tax);
    if(scalar(@tax_A) < $tax_level) { die "ERROR user specified tax level $tax_level but only " . scalar(@tax_A) . " tax tokens exist on line:\n$line\n"; }
    my $taxkey = $tax_A[($tax_level-1)];
    if(! defined $bytax_HA{$taxkey}) { 
      @{$bytax_HA{$taxkey}} = (); 
    }
    push(@{$bytax_HA{$taxkey}}, $line); 
  }
}
foreach $taxkey (sort keys %bytax_HA) { 
  $taxkey_lc = $taxkey;
  $taxkey_lc =~ tr/A-Z/a-z/;
  my $tmp_root = sprintf("%s.%s%s", $file_root, (($extra_string eq "NONE") ? "" : $extra_string), $taxkey_lc);
  $file_name = $tmp_root . ".info";
  printf("$tmp_root\n");
  open(OUT, ">", $file_name) || die "ERROR unable to open $file_name for writing";
  foreach $line (@{$bytax_HA{$taxkey}}) { 
    print OUT $line . "\n";
  }
  close OUT;
}

