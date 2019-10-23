#!/usr/bin/env perl
use strict;

my $usage = "perl vb-step2-taxinfo2muscle-qsub.pl <protein fasta file> <tax-split file> <output root>\n";
if(scalar(@ARGV) != 3) { die $usage; }

my ($fa_file, $taxsplit_file, $root) = (@ARGV);
my $cmd;
my $line;

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

# make sure you have the required files created as output by earlier vb-step*pl scripts
my $short_fa_file = $root . ".aa.short.fa";
my $info2_file    = $root . ".info2.txt";
my $tt_file       = $root . ".tt.list";
my @reqd_files_A  = ($short_fa_file, $short_fa_file.".ssi", $info2_file, $tt_file);
foreach my $reqd_file (@reqd_files_A) { 
  if(! -e $reqd_file) { die "ERROR required file $reqd_file does not exist. Did you (succesfully) run vb-step1-fasta2taxinfo.pl?"; }
  if(! -s $reqd_file) { die "ERROR required file $reqd_file exists but is empty. Did you (succesfully) run vb-step1-fasta2taxinfo.pl?"; }
}

#########################################################
# parse the taxsplit file
# top of that file
## translation_table taxonomy_level prefix
## tt2: level 7 (mammalia, etc.) 
#2 7 vertebrata-
my %taxsplit_level_H = ();
my %taxsplit_prefix_H = ();
open(IN, $taxsplit_file) || die "ERROR unable to open $taxsplit_file for reading";
while($line = <IN>) { 
  chomp $line;
  if($line !~ m/^#/) { 
    my @el_A = split(/\s+/, $line);
    if(scalar(@el_A) != 3) { die "ERROR unable to parse line in $taxsplit_file:\n$line\n"; }
    my ($tt, $level, $prefix) = (@el_A);
    $taxsplit_level_H{$tt} = $level;
    $taxsplit_prefix_H{$tt} = $prefix;
  }
}
close(IN);
#####################################################

#########################################################
# determine count of seqs
my $nseq3 = `wc -l $info2_file`;
chomp $nseq3;
#########################################################

#########################################################
# split files up by translation_table value, all seqs for each model in VADR
# must currently use the same genetic code (translation_table).
# get translation table values:
my @tt_A;
open(IN, $tt_file) || die "ERROR unable to open $tt_file for reading";
while($line = <IN>) { 
  chomp $line;
  push(@tt_A, $line);
}

my @tt_info_file_A = ();
my @tt_aa_fa_file_A = ();
my $ntt = scalar(@tt_A);
my $i;
my $model_root_file = $root . ".model.list";
# create new file $model_root_file;
open(OUT, ">", $model_root_file) || die "ERROR unable to open $model_root_file for writing"; 
close(OUT);
for($i = 0; $i < $ntt; $i++) { 
  my $tt = $tt_A[$i];
  my $tt_info_file  = $root . ".tt" . $tt . ".info";
  my $tt_aa_fa_file  = $root . ".tt" . $tt . ".aa.fa";
  push(@tt_info_file_A, $tt_info_file);
  push(@tt_aa_fa_file_A, $tt_aa_fa_file);

  #########################################################
  # fetch sequences in for each tt into a new file
  $cmd = "cat $tt_info_file | awk '{ print \$1 }' | sed 's/accver://' | sort | uniq | $easel_dir/esl-sfetch -f $short_fa_file - > $tt_aa_fa_file";
  RunCommand($cmd, 1);
  #########################################################

  #########################################################
  # split each tt group by taxonomy
  if(! exists $taxsplit_level_H{$tt})  { die "ERROR did not read taxsplit info from $taxsplit_file for translation table $tt"; }
  if(! exists $taxsplit_prefix_H{$tt}) { die "ERROR did not read taxsplit info from $taxsplit_file for translation table $tt"; }
  my $level  = $taxsplit_level_H{$tt};
  my $prefix = $taxsplit_prefix_H{$tt};
  $cmd = "perl $scripts_dir/split-by-tax.pl $tt_info_file $level $prefix $root.tt$tt >> $model_root_file ";
  RunCommand($cmd, 1);
}

# parse the model_root file
my @mdl_A = ();
my @mdl_info_file_A = ();
open(IN, $model_root_file) || die "ERROR unable to open $model_root_file for reading";
while($line = <IN>) { 
  chomp $line;
  $line =~ s/^$root\.//;
  push(@mdl_A, $line);
  push(@mdl_info_file_A, $root . "." . $line . ".info");
}

# for each model, create a new fasta file for protein seqs for that model
my $nmdl = scalar(@mdl_A);
my $m;
my @mdl_aa_fa_file_A = ();
my @mdl_aa_aln_file_A = ();
my $muscle_qsub_file = $root . ".muscle.qsub";
my $nseq4 = 0;

open(MUSCLE, ">", $muscle_qsub_file) || die "ERROR unable to open $muscle_qsub_file for writing";
for($m = 0; $m < $nmdl; $m++) { 
  my $mdl = $mdl_A[$m];
  my $mdl_info_file = $mdl_info_file_A[$m];
  # determine how many sequences are in the mdl_info_file
  my $tmp_nseq = `wc -l $mdl_info_file`;
  chomp $tmp_nseq;
  $nseq4 += $tmp_nseq;

  my $mdl_aa_fa_file = $root . "." . $mdl . ".aa.fa"; 
  my $mdl_aa_aln_file = $root . "." . $mdl . ".aa.afa"; 
  my $mdl_muscle_file = $root . "." . $mdl . ".muscle";

  push(@mdl_aa_fa_file_A, $mdl_aa_fa_file);
  push(@mdl_aa_aln_file_A, $mdl_aa_aln_file);

  $cmd = "cat $mdl_info_file | awk '{ print \$1 }' | sed 's/accver://' | sort | uniq | $easel_dir/esl-sfetch -f $short_fa_file - > $mdl_aa_fa_file";
  RunCommand($cmd, 1);

  print MUSCLE ("qsub -N $mdl -b y -v SGE_FACILITIES -P unified -S /bin/bash -cwd -V -j n -o /dev/null -e $mdl.muscle.err -l m_mem_free=8G,h_rt=2880000,mem_free=8G,h_vmem=8G -m n \"/usr/local/muscle/3.7/bin/muscle -maxmb 4096 -in $mdl_aa_fa_file -out $mdl_aa_aln_file > $mdl_muscle_file\"\n");
}
close(MUSCLE);

if($nseq4 != $nseq3) { 
  die "ERROR, the $nmdl translation_table/taxonomic_group models only comprise $nseq4 of the $nseq3 files in $info2_file\nYou need to specify additional taxonomic groups in $taxsplit_file";
}
printf("\nScript to submit $nmdl muscle jobs to the farm is in:\n$muscle_qsub_file\n");
printf("\nRun that script, wait for all jobs to finish, then run:\n\$VADRBUILDTOOLSDIR/vb-step3-muscle-alns2cmbuild-qsub.pl $model_root_file\n");
  
#########################################################

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

