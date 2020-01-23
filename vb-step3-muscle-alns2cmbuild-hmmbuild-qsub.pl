#!/usr/bin/env perl
use strict;

my $usage = "perl vb-step3-muscle-alns2cmbuild-qsub.pl <model list>\n";
if(scalar(@ARGV) != 1) { die $usage; }

my $version = "0.02";

my ($model_root_file) = (@ARGV);
my $root = $model_root_file;
if($root !~ m/\.model\.list$/) { 
  die "ERROR unable to parse model list file name $model_root_file"; 
}
$root =~ s/\.model\.list$//;

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
if(! exists($ENV{"VADRINFERNALDEVDIR"})) { 
  die "ERROR, the environment variable VADRINFERNALDEVDIR is not set";
}
if(! (-d $ENV{"VADRINFERNALDEVDIR"})) { 
  die "ERROR, the directory specified by your environment variable VADRINFERNALDEVDIR does not exist.\n"; 
}    
if(! exists($ENV{"VADRHMMERDIR"})) { 
  die "ERROR, the environment variable VADRHMMERDIR is not set";
}
if(! (-d $ENV{"VADRHMMERDIR"})) { 
  die "ERROR, the directory specified by your environment variable VADRHMMERDIR does not exist.\n"; 
}    

my $scripts_dir  = $ENV{"VADRBUILDTOOLSDIR"} . "/scripts";
my $easel_dir    = $ENV{"VADREASELDIR"};
my $infernal_dir = $ENV{"VADRINFERNALDEVDIR"};
my $hmmer_dir    = $ENV{"VADRHMMERDIR"};
my $cmd;

# parse the model_root file
my @mdl_A = ();
open(IN, $model_root_file) || die "ERROR unable to open $model_root_file for reading";
while(my $line = <IN>) { 
  chomp $line;
  my $mdl = $line;
  $mdl =~ s/^$root\.//;
  push(@mdl_A, $mdl);

  my @reqd_files_A = ();
  my $info_file         = $root . "." . $mdl . ".info";
  my $aa_fa_file        = $root . "." . $mdl . ".aa.fa"; 
  my $aa_aln_file       = $root . "." . $mdl . ".aa.afa"; 
  my @reqd_files_A = ($info_file, $aa_fa_file, $aa_aln_file);
  foreach my $reqd_file (@reqd_files_A) { 
    if(! -e $reqd_file) { die "ERROR required file $reqd_file does not exist. Did you (succesfully) run vb-step2-taxinfo2muscle-qsub.pl?"; }
    if(! -s $reqd_file) { die "ERROR required file $reqd_file exists but is empty. Did you (succesfully) run vb-step2-taxinfo2muscle-qsub.pl?"; }
  }
}

my $build_qsub_file  = $root . ".build.qsub";
open(BUILD,  ">", $build_qsub_file)  || die "ERROR unable to open $build_qsub_file for writing";

# for each model: 
my $nmdl = scalar(@mdl_A);
for(my $m = 0; $m < $nmdl; $m++) { 
  my $mdl = $mdl_A[$m];

  # required input files (we checked that these existed above when we created @mdl_A)
  my $info_file         = $root . "." . $mdl . ".info";
  my $aa_fa_file        = $root . "." . $mdl . ".aa.fa"; 
  my $aa_aln_file       = $root . "." . $mdl . ".aa.afa"; 

  # output files that we'll create
  my $source_list       = $root . "." . $mdl . ".source.list";
  my $source_nt_fa_file = $root . "." . $mdl . ".source.nt.fa";
  my $nt_fa_file        = $root . "." . $mdl . ".nt.fa";
  my $map_file          = $root . "." . $mdl . ".map";
  my $nt_aln_file       = $root . "." . $mdl . ".nt.afa";
  my $nt_stk_file       = $root . "." . $mdl . ".nt.stk";
  my $alistat_file      = $root . "." . $mdl . ".alistat";

  # create list of 'source' nt accessions
  $cmd = "perl $scripts_dir/parse-by-source.pl $info_file | sort | uniq > $source_list";
  RunCommand($cmd, 1);
  
  # fetch the source sequences
  $cmd = "/netopt/ncbi_tools64/bin/idfetch -t 5 -c 1 -G $source_list | perl $scripts_dir/fasta-ncbi-idfetch-name-long-to-short.pl > $source_nt_fa_file";
  RunCommand($cmd, 1);

  # index the source sequence files
  $cmd = "$easel_dir/esl-sfetch --index $source_nt_fa_file";
  RunCommand($cmd, 1);
  
  # fetch the subseqs that encode the proteins from their sources
  $cmd = "perl $scripts_dir/fetch-given-info.pl $info_file $source_nt_fa_file $nt_fa_file $map_file";
  RunCommand($cmd, 1);

  # create a nucleotide alignment by translating the protein alignment in place
  $cmd = "perl $scripts_dir/aa-to-nt-alignment.pl $aa_aln_file $aa_fa_file $nt_fa_file $map_file > $nt_aln_file";
  RunCommand($cmd, 1);
  
  # reformat the nucleotide alignment to stockholm format
  $cmd = "$easel_dir/esl-reformat stockholm $nt_aln_file > $nt_stk_file";
  RunCommand($cmd, 1);

  # get stats on the alignment
  $cmd = "$easel_dir/esl-alistat $nt_stk_file > $alistat_file";
  RunCommand($cmd, 1);

  # make qsub commands for building the models
  my $cm_name = $root . "." . $mdl;
  my @ere_opt_A  = ("1.0");
  my @ere_name_A = ("1p0");
  my $nere = scalar(@ere_opt_A);

  for(my $i = 0; $i < $nere; $i++) { 
    my $cm_root = $cm_name . "." . $ere_name_A[$i];
    my $cm_file_name = $cm_root . ".vadr.cm";
    my $cmbuild_file_name = $cm_root . ".vadr.cmbuild";
    my $cmbuild_cmd = "$infernal_dir/cmbuild -F -n $cm_name --emaxseq 100000 --noss --ere $ere_opt_A[$i] $cm_file_name $nt_stk_file > $cmbuild_file_name";
    printf BUILD ("qsub -N $cm_root -b y -v SGE_FACILITIES -P unified -S /bin/bash -cwd -V -j n -o /dev/null -e $cm_root.err -l m_mem_free=8G,h_rt=2880000,mem_free=8G,h_vmem=8G -m n \"$cmbuild_cmd\"\n");
  }

  # build HMM
  my $hmm_name = $cm_name;
  my $hmm_root = $hmm_name;
  my $hmm_file_name = $hmm_root . ".vadr.hmm";
  my $hmmbuild_file_name = $hmm_root . ".vadr.hmmbuild";
  my $hmmbuild_cmd = "$hmmer_dir/hmmbuild -n $hmm_name $hmm_file_name $aa_aln_file > $hmmbuild_file_name";
  printf BUILD ("qsub -N $hmm_root.hmm -b y -v SGE_FACILITIES -P unified -S /bin/bash -cwd -V -j n -o /dev/null -e $hmm_root.hmm.err -l m_mem_free=8G,h_rt=2880000,mem_free=8G,h_vmem=8G -m n \"$hmmbuild_cmd\"\n");
}
close(BUILD);
printf("\nScript to submit $nmdl cmbuild and hmmbuild jobs to the farm is in:\n$build_qsub_file\n");
printf("\nRun that script, wait for all jobs to finish, then run:\n");
printf("perl \$VADRBUILDTOOLSDIR/vb-step4-create-vadr-files.pl $model_root_file <name of vadr model dir to create> <gene value (use _ for space) <product value (use _ for space)>\n");


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

