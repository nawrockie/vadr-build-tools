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
  print $source . "\n";
}
#print("nsimple: $nsimple\n");
#print("nhard:   $nhard\n");

