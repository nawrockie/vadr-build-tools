#!/usr/bin/env perl
use strict;

my $usage = "perl fetch-given-info.pl <info file> <fasta file to fetch from> <fasta file to create>";
if(scalar(@ARGV) != 3) { die $usage; }
my ($info_file, $source_fa_file, $output_fa_file) = (@ARGV);

if(! exists($ENV{"VADRBUILDTOOLSDIR"})) { 
  die "ERROR, the environment variable VADRBUILDTOOLSDIR is not set";
}
if(! (-d $ENV{"VADRBUILDTOOLSDIR"})) { 
  die "ERROR, the directory specified by your environment variable VADRBUILDTOOLSDIR does not exist.\n"; 
}    
if(! exists($ENV{"VADREASELDIR"})) { 
  die "ERROR, the environment variable VADREASELDIR is not set";
}
if(! (-d $ENV{"VADREASELDIR"})) { 
  die "ERROR, the directory specified by your environment variable VADREASELDIR does not exist.\n"; 
}    

my $scripts_dir = $ENV{"VADRBUILDTOOLSDIR"} . "/scripts";
my $easel_dir = $ENV{"VADREASELDIR"};

# remove output file if it exists
if(-e $output_fa_file) { unlink $output_fa_file; }

my $cmd = undef;

# for each file that has more than one segment, we have to fetch it independently
# for each file that only has one segment we can create a single esl-sfetch input file
# and fetch all of those seqs all at once at the end

my $one_seg_sfetch_file = "tmp.sfetch";
open(SFETCH, ">", $one_seg_sfetch_file) || die "ERROR unable to open $one_seg_sfetch_file for writing";

my $n2fetch = 0;
open(IN, $info_file) || die "ERROR unable to open $info_file";
while(my $line = <IN>) { 
  chomp $line;
  my @el_A = split(/\s+/, $line);
  my $coded_by_line = $el_A[1];
  my $strand = "";

  # first extract the source accession
  my $source;
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

  my $ncbi_location = $coded_by_line;
  $ncbi_location =~ s/\Q$source\E\://g; 
  $ncbi_location =~ s/^coded_by\://;
  
  my $vadr_coords = vdr_CoordsFromLocation($ncbi_location, undef);
  
  my @start_A  = ();
  my @stop_A   = ();
  my @strand_A = ();
  vdr_FeatureStartStopStrandArrays($vadr_coords, \@start_A, \@stop_A, \@strand_A, undef);
  my $nsgm = scalar(@start_A);
  # create name from $vadr_coords string
  my $vadr_coords_no_strand = $vadr_coords;
  $vadr_coords_no_strand =~ s/\:\+//g;
  $vadr_coords_no_strand =~ s/\:\-//g;
  my $newname = $source . "/" . $vadr_coords_no_strand;

  my $sfetch_str = sprintf("%s %d %d %s\n", $newname, $start_A[0], $stop_A[0], $source);
  if($nsgm == 1) { 
    printf SFETCH $sfetch_str;
    $n2fetch++;
  }
  else { 
    # multi-segment, fetch first segment with seq name, then fetch each other one stripping name away
    # to concatenate each segment into one sequence
    my $tmp_sfetch = $source . ".$vadr_coords_no_strand.tmp.sfetch";
    open(TMP, ">", $tmp_sfetch) || die "ERROR unable to open $tmp_sfetch for writing";
    printf TMP $sfetch_str; 
    close(TMP);
    $cmd = "$easel_dir/esl-sfetch -Cf $source_fa_file $tmp_sfetch >> $output_fa_file\n";
    RunCommand($cmd, 1);
    
    $cmd = "rm $tmp_sfetch";
    RunCommand($cmd, 1);

    for(my $i = 1; $i < $nsgm; $i++) { 
      $cmd = "$easel_dir/esl-sfetch -c $start_A[$i]..$stop_A[$i] $source_fa_file $source | grep -v ^\\>  >> $output_fa_file";
      RunCommand($cmd, 1);
    }
  }
}
close(SFETCH);

if($n2fetch > 0) { 
  # fetch all single segment seqs
  $cmd = "$easel_dir/esl-sfetch -Cf $source_fa_file $one_seg_sfetch_file >> $output_fa_file";
  RunCommand($cmd, 1);
}
#unlink $one_seg_sfetch_file;

#################################################################
# SUBROUTINES from VADR 0.991
#################################################################
# Subroutine: vdr_CoordsFromLocation
# Incept:     EPN, Wed Mar 13 14:17:08 2019
# 
# Purpose:    Convert a GenBank file 'location' value to 
#             a coords string in the format:
#             <start1>-<stop2>:<strand1>,<start2>-<stop2>:<strand2>,...,<startN>-<stopN>:<strandN>
# 
#             Any carrots before start/stop positions in the 
#             location string are removed.
#             See vdr_CoordsFromLocationWithCarrots() for 
#             a similar subroutine that keeps carrots.
#
#             This function has to call itself recursively in some
#             cases.
# 
# Arguments:
#   $location: GenBank file location string
#   $FH_HR:    REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
# 
# Dies:       if unable to parse $location
#
# Ref: GenBank release notes (release 230.0) as of this writing
#      and
#      https://www.ncbi.nlm.nih.gov/genbank/samplerecord/
#################################################################
sub vdr_CoordsFromLocation { 
  my $sub_name = "vdr_CoordsFromLocation";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($location, $FH_HR) = @_;

  # Examples we can parse: 
  # $location                          return value
  # ---------------------------------  -----------------
  # 1..200                             1..200:+
  # <1..200                            1..200:+
  # 100..>200                          100..200:+
  # <1..>200                           1..200:+
  # complement(1..200)                 200..1:-
  # join(1..200,300..400)              1..200:+,300..400:+
  # complement(join(1..200,300..400))  400..300:-,200..1:-
  # join(1..200,complement(300..400))  1..200:+,400..300:- ! NOT SURE IF THIS IS CORRECT !
  # join(complement(300..400),1..200)  400..300:-,1..200:+ ! NOT SURE IF THIS IS CORRECT !

  my $ret_val = "";
  if($location =~ /^join\((.+)\)$/) { 
    my $location_to_join = $1;
    $ret_val = vdr_CoordsFromLocation($location_to_join, $FH_HR);
  }
  elsif($location =~ /^complement\((.+)\)$/) { 
    my $location_to_complement = $1;
    my $coords_to_complement = vdr_CoordsFromLocation($location_to_complement, $FH_HR);
    $ret_val = vdr_CoordsComplement($coords_to_complement, $FH_HR);
  }
  elsif($location =~ /\,/) { 
    # not wrapped in join() or complement(), but multiple segments
    foreach my $location_el (split(",", $location)) { 
      if($ret_val ne "") { $ret_val .= ","; }
      $ret_val .= vdr_CoordsFromLocation($location_el, $FH_HR);
    }
  }
  elsif($location =~ /^\<?(\d+)\.\.\>?(\d+)$/) { 
    $ret_val = $1 . ".." . $2 . ":+"; # a recursive call due to the complement() may complement this
  }
  elsif($location =~ /^\<?(\d+)$/) { # single nucleotide
    $ret_val = $1 . ".." . $1 . ":+"; # a recursive call due to the complement() may complement this
  }
  elsif($location =~ /^\>?(\d+)$/) { # single nucleotide
    $ret_val = $1 . ".." . $1 . ":+"; # a recursive call due to the complement() may complement this
  }
  else { 
    die "ERROR in $sub_name, unable to parse location token $location";
  }

  return $ret_val;
}

#################################################################
# Subroutine: vdr_CoordsFromLocationWithCarrots
# Incept:     EPN, Fri Apr 12 11:51:16 2019
# 
# Purpose:    Convert a GenBank file 'location' value to 
#             a coords string in the format:
#             <start1>-<stop2>:<strand1>,<start2>-<stop2>:<strand2>,...,<startN>-<stopN>:<strandN>
#             
#             <startN>: may begin with "<" carrot.
#             <stopN>: may begin with ">" carrot.
#
#             vdr_CoordsFromLocation() does the same thing but 
#             removes carrots.
#
#             This function has to call itself recursively in some
#             cases.
# 
# Arguments:
#   $location: GenBank file location string
#   $FH_HR:    REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
# 
# Dies:       if unable to parse $location
#
# Ref: GenBank release notes (release 230.0) as of this writing
#      and
#      https://www.ncbi.nlm.nih.gov/genbank/samplerecord/
#################################################################
sub vdr_CoordsFromLocationWithCarrots { 
  my $sub_name = "vdr_CoordsFromLocationWithCarrots";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($location, $FH_HR) = @_;

  # Examples we can parse: 
  # $location                          return value
  # ---------------------------------  -----------------
  # 1..200                             1..200:+
  # <1..200                            <1..200:+
  # 100..>200                          100..>200:+
  # <1..>200                           <1..>200:+
  # complement(1..200)                 200..1:-
  # join(1..200,300..400)              1..200:+,300..400:+
  # complement(join(1..200,300..400))  400..300:-,200..1:-
  # join(1..200,complement(300..400))  1..200:+,400..300:- ! NOT SURE IF THIS IS CORRECT !
  # join(complement(300..400),1..200)  400..300:-,1..200:+ ! NOT SURE IF THIS IS CORRECT !

  my $ret_val = "";
  if($location =~ /^join\((.+)\)$/) { 
    my $location_to_join = $1;
    $ret_val = vdr_CoordsFromLocationWithCarrots($location_to_join, $FH_HR);
  }
  elsif($location =~ /^complement\((.+)\)$/) { 
    my $location_to_complement = $1;
    my $coords_to_complement = vdr_CoordsFromLocationWithCarrots($location_to_complement, $FH_HR);
    $ret_val = vdr_CoordsComplementWithCarrots($coords_to_complement, $FH_HR);
  }
  elsif($location =~ /\,/) { 
    # not wrapped in join() or complement(), but multiple segments
    foreach my $location_el (split(",", $location)) { 
      if($ret_val ne "") { $ret_val .= ","; }
      $ret_val .= vdr_CoordsFromLocationWithCarrots($location_el, $FH_HR);
    }
  }
  elsif($location =~ /^(\<?\d+\.\.\>?\d+)$/) { 
    $ret_val = $1 . ":+"; # a recursive call due to the complement() may complement this
  }
  elsif($location =~ /^(\<?\d+)$/) { # single nucleotide
    $ret_val = $1 . ".." . $1 . ":+"; # a recursive call due to the complement() may complement this
  }
  elsif($location =~ /^(\>?\d+)$/) { # single nucleotide
    $ret_val = $1 . ".." . $1 . ":+"; # a recursive call due to the complement() may complement this
  }
  else { 
    die "ERROR in $sub_name, unable to parse location token $location";
  }

  return $ret_val;
}

#################################################################
# Subroutine: vdr_CoordsComplement
# Incept:     EPN, Wed Mar 13 15:00:24 2019
# 
# Purpose:    Complement a coords string by complementing all
#             elements within it. Removes carrots "<" and ">"
#             before start and stop positions, if they exist.
#             See vdr_CoordsComplementWithCarrots() to keep
#             carrots.
# 
# Arguments:
#   $coords:   coords string to complement
#   $FH_HR:    REF to hash of file handles, including "log" and "cmd"
#
# Returns:    complemented $coords
# 
# Dies:       if unable to parse $coords, or any segment in $coords
#             is already on the negative strand.
#
#################################################################
sub vdr_CoordsComplement { 
  my $sub_name = "vdr_CoordsComplement";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($coords, $FH_HR) = @_;

  # Examples we can parse: 
  # $coords                  return value
  # -----------------------  -----------------
  # 1-200:+                  200-1:-
  # 1-200:+,300-400:+        400-300:-,200-1:-

  my $ret_val = "";
  my @el_A = split(",", $coords);
  for(my $i = scalar(@el_A)-1; $i >= 0; $i--) { 
    if($el_A[$i] =~ /^\<?(\d+)\.\.\>?(\d+)\:\+/) { 
      my ($start, $stop) = ($1, $2, $3, $4);
      if($ret_val ne "") { $ret_val .= ","; }
      $ret_val .= $stop . ".." . $start . ":-";
    }
    else { 
      die "ERROR in $sub_name, unable to parse coords token $coords";
    }
  }

  # printf("\tin $sub_name, coords: $coords ret_val: $ret_val\n");

  return $ret_val;
}

#################################################################
# Subroutine: vdr_CoordsComplementWithCarrots
# Incept:     EPN, Fri Apr 12 11:49:06 2019
# 
# Purpose:    Complement a coords string by complementing all
#             elements within it, and reverse any carrots. 
#             Just like vdr_CoordsComplement() but keeps
#             and complements carrots.
# 
# Arguments:
#   $coords:   coords string to complement
#   $FH_HR:    REF to hash of file handles, including "log" and "cmd"
#
# Returns:    complemented $coords
# 
# Dies:       if unable to parse $coords, or any segment in $coords
#             is already on the negative strand.
#
#################################################################
sub vdr_CoordsComplementWithCarrots { 
  my $sub_name = "vdr_CoordsComplementWithCarrots";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($coords, $FH_HR) = @_;

  # Examples we can parse: 
  # $coords                  return value
  # -----------------------  -----------------
  # 1-200:+                  200-1:-
  # 1-200:+,300-400:+        400-300:-,200-1:-

  my $ret_val = "";
  my @el_A = split(",", $coords);
  for(my $i = scalar(@el_A)-1; $i >= 0; $i--) { 
    if($el_A[$i] =~ /^(\<?)(\d+)\.\.(\>?)(\d+)\:\+/) { 
      my ($start_carrot, $start, $stop_carrot, $stop) = ($1, $2, $3, $4);
      if($start_carrot eq "<") { $start_carrot = ">"; }
      if($stop_carrot  eq ">") { $stop_carrot  = "<"; }
      if($ret_val ne "") { $ret_val .= ","; }
      $ret_val .= $stop_carrot . $stop . ".." . $start_carrot . $start . ":-";
    }
    else { 
      die "ERROR in $sub_name, unable to parse coords token $coords";
    }
  }

  # printf("\tin $sub_name, coords: $coords ret_val: $ret_val\n");

  return $ret_val;
}

#################################################################
# Subroutine: vdr_FeatureStartStopStrandArrays()
# Incept:     EPN, Sat Mar  9 05:50:10 2019
#
# Synopsis: Given a comma separated coords string, parse it, 
#           validate it, and fill @{$start_AR}, @{$stop_AR} and
#           @{$strand_AR} based on it.
# 
# Arguments:
#  $coords:       coordinate string
#  $start_AR:     REF to start position array to fill here, FILLED here, can be undef
#  $stop_AR:      REF to stop position array to fill here, FILLED here, can be undef
#  $strand_AR:    REF to strand array to fill here with "+" or "-", FILLED here, can be undef
#  $FH_HR:        REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies: if unable to parse $coords
#
#################################################################
sub vdr_FeatureStartStopStrandArrays {
  my $sub_name = "vdr_FeatureStartStopStrandArrays";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($coords, $start_AR, $stop_AR, $strand_AR, $FH_HR) = @_;
  if(! defined $coords) { 
    die "ERROR in $sub_name, coords is undefined";
  }

  my @start_A  = ();
  my @stop_A   = ();
  my @strand_A = ();
  my ($start, $stop, $strand, $sgm_idx);
  my @coords_A  = split(",", $coords);
  my $nsgm = scalar(@coords_A);
  for($sgm_idx = 0; $sgm_idx < $nsgm; $sgm_idx++) { 
    ($start, $stop, $strand) = vdr_CoordsTokenParse($coords_A[$sgm_idx], $FH_HR);
    # vdr_CoordsTokenParse() will fail if unable to parse $coords_A[$sgm_idx]
    push(@start_A,  $start);
    push(@stop_A,   $stop);
    push(@strand_A, $strand); 
  }

  if(defined $start_AR)  { @{$start_AR}   = @start_A;  }
  if(defined $stop_AR)   { @{$stop_AR}    = @stop_A;   }
  if(defined $strand_AR) { @{$strand_AR}  = @strand_A;  }

  return;
}
#################################################################
# Subroutine: vdr_CoordsTokenParse()
# Incept:     EPN, Tue Mar 26 06:15:09 2019
#
# Synopsis: Given a single coords token, validate it, 
#           and return its start, stop, strand values.
# 
# Arguments:
#  $coords_tok:   coordinate token
#  $FH_HR:        REF to hash of file handles, including "log" and "cmd"
#
# Returns:    3 values:
#             $start:  start position
#             $stop:   stop position
#             $strand: strand
#
# Dies: if unable to parse $coords
#
#################################################################
sub vdr_CoordsTokenParse {
  my $sub_name = "vdr_CoordsTokenParse";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($coords_tok, $FH_HR) = @_;
  if(! defined $coords_tok) { 
    die "ERROR in $sub_name, coords is undefined";
  }
  if($coords_tok =~ /^\<?(\d+)\.\.\>?(\d+)\:([\+\-])$/) { 
    return ($1, $2, $3);
  }
  die "ERROR in $sub_name, unable to parse coords token $coords_tok";

  return; # NEVER REACHED
}
#################################################################
# Subroutine: RunCommand()
# Incept:     EPN, Thu Feb 11 13:32:34 2016 [dnaorg.pm]
#
# Purpose:     Runs a command using system() and exits in error 
#              if the command fails. If $be_verbose, outputs
#              the command to stdout.
#
# Arguments:
#   $cmd:         command to run, with a "system" command;
#   $be_verbose:  '1' to output command to stdout before we run it, '0' not to
#
# Returns:    void
#
# Dies:       if $cmd fails
#################################################################
sub RunCommand {
  my $sub_name = "RunCommand()";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($cmd, $be_verbose) = (@_);
  
  if($be_verbose) { 
    print ("Running cmd: $cmd\n"); 
  }

  system($cmd);

  if($? != 0) { 
    die "ERROR in $sub_name, the following command failed:\n$cmd\n";
  }

  return;
}

