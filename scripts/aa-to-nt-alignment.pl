use strict;

my $usage = "perl aa-to-nt-alignment.pl <aa afa alignment file> <aa fa unaligned> <nt fa unaligned file> <map file>\n";
if(scalar(@ARGV) != 4) { die $usage; }

my ($aa_aln_file, $aa_unaln_file, $nt_unaln_file, $map_file) = (@ARGV);

open(MAP, $map_file) || die "ERROR unable to open map file";
my $line;
my %nt_seq_H = ();
my %nt_len_H = ();
my %aa2nt_map_H = ();
my %nt2aa_map_H = ();

while($line = <MAP>) { 
  chomp $line;
  my @el_A = split(/\s+/, $line);
  if(scalar(@el_A) != 2) { die "ERROR unable to parse map line $line"; }
  my ($aa_name, $nt_name) = ($el_A[0], $el_A[1]);
  if((defined $aa2nt_map_H{$aa_name}) && ($aa2nt_map_H{$aa_name} ne $nt_name)) { die "ERROR read two different map lines for protein $aa_name"; }
  if((defined $nt2aa_map_H{$aa_name}) && ($nt2aa_map_H{$nt_name} ne $aa_name)) { die "ERROR read two different map lines for nt $nt_name"; }
  $aa2nt_map_H{$aa_name} = $nt_name;
  $nt2aa_map_H{$nt_name} = $aa_name;
}
close(MAP);

# parse unaligned nt file
open(NT, $nt_unaln_file) || die "ERROR unable to open nt file";
my $nt_name = undef;
while($line = <NT>) { 
  chomp $line;
  if($line =~ /^\>(\S+)/) { 
    $nt_name = $1;
  }
  elsif($line =~ m/\w/) { 
    if(! defined $nt_seq_H{$nt_name}) { 
      $nt_seq_H{$nt_name} = ""; 
      $nt_len_H{$nt_name} = 0; 
    }
    $nt_seq_H{$nt_name} .= $line;
    $nt_len_H{$nt_name} += length($line);
  }
}
close(NT);

# parse unaligned aa file
open(AA, $aa_unaln_file) || die "ERROR unable to open aa file";
my $aa_name = undef;
my %aa_seq_H = ();
my %aa_len_H = ();
while($line = <AA>) { 
  chomp $line;
  if($line =~ /^\>(\S+)/) { 
    $aa_name = $1;
  }
  elsif($line =~ m/\w/) { 
    if(! defined $aa_seq_H{$aa_name}) { 
      $aa_seq_H{$aa_name} = ""; 
      $aa_len_H{$aa_name} = 0; 
    }
    $aa_seq_H{$aa_name} .= $line;
    $aa_len_H{$aa_name} += length($line);
  }
}
close(AA);

# make sure that length of nt and aa seq make sense
# if aa_len is M, nt_len should be 3*M, (3*M)+1, (3*M)+2 or (3*M)+3
# if it does not we won't output that nucleotide seq
# because we'll be reverse translating it incorrectly
my %skip_H = (); # names of nt seqs to skip because lengths don't make sense
foreach my $nt_name (sort keys %nt_seq_H) { 
  if(! defined $nt2aa_map_H{$nt_name}) { 
    die "ERROR no mapping aa seq for nt seq $nt_name"; 
  }
  my $aa_name = $nt2aa_map_H{$nt_name};
  my $aa_len = $aa_len_H{$aa_name};
  my $nt_len = $nt_len_H{$nt_name};
  if(($nt_len != (3 * $aa_len)) && ($nt_len != ((3 * $aa_len) + 3))) { 
    if   ($nt_len == ((3 * $aa_len)+1))  { ; } # okay, we'll just not output the final nt
    elsif($nt_len == ((3 * $aa_len)+2))  { ; } # okay, we'll just not output the final 2 nt
    elsif($nt_len < (3 * $aa_len))       { $skip_H{$nt_name} = 1; print STDERR "SKIPPING nt_len < 3*M     (nt $nt_name: $nt_len aa $aa_name: $aa_len)\n"; }
    elsif($nt_len > ((3 * $aa_len) + 3)) { $skip_H{$nt_name} = 1; print STDERR "SKIPPING nt_len > 3*(M+1) (nt $nt_name: $nt_len aa $aa_name: $aa_len)\n"; }
    else { die "ERROR uncovered case (nt: $nt_len aa: $aa_len)\n"; }   
  }
}

open(AA, $aa_aln_file) || die "ERROR unable to open aa file";
my $aa_name = undef;
my $nt_seq = undef;
my @nt_seq_A = ();
my $nt_idx; 
my $nt_len;
my $aa_name;
my $nseq = 0;
my $add_3_gaps = 0;
my $skip_flag = 0;
while($line = <AA>) { 
  chomp $line;
  if($line =~ /^\>(\S+)/) { 
    # first finish off prev seq if there is one:
    if(! $skip_flag) { 
      $add_3_gaps = 1;
      if(($nseq > 0) && ($nt_idx < $nt_len)) { # check and add the stop codon if 
        if(($nt_len - $nt_idx) > 3) { 
          die "ERROR more than 3 nt left when adding stop for $nt_name nt_idx: $nt_idx nt_len: $nt_len\n"; 
        } 
        elsif(($nt_idx + 3) == ($nt_len)) { # only add full stops (we could change this to add partial ones too) 
          my $stop_codon = $nt_seq_A[$nt_idx] . $nt_seq_A[($nt_idx+1)] . $nt_seq_A[($nt_idx+2)];
          print $stop_codon . "\n";
          $add_3_gaps = 0;
        }
      }
      if(($nseq > 0) && ($add_3_gaps)) { print "---\n"; }
    }

    $aa_name = $1;
    if(! defined $aa2nt_map_H{$aa_name}) { die "ERROR read $aa_name but no nt seq exists in map"; }
    $nt_name = $aa2nt_map_H{$aa_name};
    if(! defined $nt_seq_H{$nt_name}) { die "ERROR read $aa_name with mapping nt $nt_name, but no seq exists from nt file"; }
    $nt_seq = $nt_seq_H{$nt_name};
    @nt_seq_A = split("", $nt_seq);
    $nt_idx = 0;
    $nt_len = scalar(@nt_seq_A);
    $skip_flag = (defined $skip_H{$nt_name}) ? 1 : 0;
    if(! $skip_flag) { 
      print(">$nt_name $aa_name\n");
      $nseq++;
    }
  }
  elsif($line ne "") { 
    if(! $skip_flag) { 
      my @aa_seq_A = split("", $line); 
      my $alen = scalar(@aa_seq_A);
      for(my $i = 0; $i < $alen; $i++) { 
        if(($aa_seq_A[$i] eq "-") || ($aa_seq_A[$i] eq ".")) { 
          print "---"; 
        }
        else { 
          if(($nt_idx+2) >= $nt_len) { die "ERROR went off the end of the sequence for $nt_seq (idx:" . ($nt_idx+2) . " >= len:$nt_len)\n"; }
          print $nt_seq_A[$nt_idx] . $nt_seq_A[($nt_idx+1)] . $nt_seq_A[($nt_idx+2)];
          $nt_idx += 3;
        }
      }
      print "\n";
    }
  }
}
if(! $skip_flag) { 
  $add_3_gaps = 1;
  if(($nseq > 0) && ($nt_idx < $nt_len)) { # check and add the stop codon
    if(($nt_len - $nt_idx) > 3) { 
      die "ERROR more than 3 nt left when adding stop for $nt_name nt_idx: $nt_idx nt_len: $nt_len\n"; 
    } 
    elsif(($nt_idx + 3) == ($nt_len)) { # only add full stops (we could change this to add partial ones too) 
      my $stop_codon = $nt_seq_A[$nt_idx] . $nt_seq_A[($nt_idx+1)] . $nt_seq_A[($nt_idx+2)];
      print $stop_codon . "\n";
      $add_3_gaps = 0;
    }
  }
  if(($nseq > 0) && ($add_3_gaps)) { print "---\n"; }
}
close(AA);

