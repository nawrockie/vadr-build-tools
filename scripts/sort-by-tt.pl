my $usage = "perl sort-by-tt.pl <info file> <output root>\n";
if(scalar(@ARGV) != 2) { die $usage; }
my ($info_file, $root) = (@ARGV);

printf("#Translation table values:\n");
open(IN, $info_file) || die "ERROR unable to open $info_file";
while($line = <IN>) { 
  chomp $line;
  my @el_A = split(/\s+/, $line);
  my $ttline = $el_A[2];
  if($ttline =~ /transl_table:(\d+)/) { 
    $tt = $1;
    $out_file = $root . ".tt" . $tt . ".info";
    if(! exists $open_H{$tt}) { 
      open($FH_H{$tt}, ">", $out_file) || die "ERROR unable to open $out_file for writing";
      $open_H{$tt} = 1;
      printf("$tt\n");
    }
    $fh = $FH_H{$tt};
    print $fh ($line . "\n");
  }
}
foreach $fh (keys %FH_H) { 
  close $fh;
}
