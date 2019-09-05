my $usage = "perl parse-by-source.pl <info file>";
if(scalar(@ARGV) != 1) { die $usage; }
my ($info_file) = (@ARGV);

my $nsimple = 0;
my $nhard = 0;
open(IN, $info_file) || die "ERROR unable to open $info_file";
while($line = <IN>) { 
  chomp $line;
  my @el_A = split(/\s+/, $line);
  my $coded_by_line = $el_A[1];
  my $strand = "+";
  if($coded_by_line =~ /coded_by:complement\(([^\:]+)\:(\S+)\)/) { 
    ($source, $coords) = ($1, $2);
    $strand = "+";
  }
  elsif($coded_by_line =~ /coded_by:([^\:]+)\:(\S+)/) { 
    ($source, $coords) = ($1, $2);
    $strand = "-";
  }
  else { 
    die "ERROR unable to parse line $line\n"; 
  }
  print $source . "\n";
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
}
#print("nsimple: $nsimple\n");
#print("nhard:   $nhard\n");

