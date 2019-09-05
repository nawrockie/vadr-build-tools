#!/usr/bin/env perl
use strict;

my $usage = "perl vb-step2.pl <model list>\n";
if(scalar(@ARGV) != 1) { die $usage; }

my ($model_root_file) = (@ARGV);
my $root = $model_root_file;
if($root !~ m/\.model\.list$/) { 
  die "ERROR unable to parse model list file name $model_root_file"; 
}
$root =~ s/\.model\.list$//;

my $scripts_dir = "./scripts";
my $easel_dir = "/usr/local/infernal/1.1.2/bin";

my $cmd;
my $line;

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

# for each model: 
my $nmdl = scalar(@mdl_A);
my $m;
my @mdl_aa_fa_file_A = ();
my @mdl_aa_aln_file_A = ();

my $cmbuild_qsub_file = $root . ".cmbuild.qsub";
open(CMBUILD, ">", $cmbuild_qsub_file) || die "ERROR unable to open $cmbuild_qsub_file for writing";

for($m = 0; $m < $nmdl; $m++) { 
  my $mdl = $mdl_A[$m];
  my $mdl_info_file = $mdl_info_file_A[$m];
  my $mdl_aa_fa_file = $root . "." . $mdl . ".aa.fa"; 
  my $mdl_aa_aln_file = $root . "." . $mdl . ".aa.afa"; 
  push(@mdl_aa_fa_file_A, $mdl_aa_fa_file);
  push(@mdl_aa_aln_file_A, $mdl_aa_aln_file);

  my $source_list = $root . "." . $mdl . ".source.list";
  my $info_file   = $root . "." . $mdl . ".info";
  my $source_nt_fa_file = $root . "." . $mdl . ".source.nt.fa";
  my $nt_fa_file = $root . "." . $mdl . ".nt.fa";
  my $nt_aln_file = $root . "." . $mdl . ".nt.afa";
  my $nt_stk_file = $root . "." . $mdl . ".nt.stk";
  my $aa_fa_file = $root . "." . $mdl . ".aa.fa"; 
  my $aa_aln_file = $root . "." . $mdl . ".aa.afa"; 
  my $map_file = $root . "." . $mdl . ".map";
  my $alistat_file = $root . "." . $mdl . ".alistat";

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
  $cmd = "perl $scripts_dir/parse-to-sfetch.pl $info_file | $easel_dir/esl-sfetch -Cf $source_nt_fa_file - > $nt_fa_file";
  RunCommand($cmd, 1);

  # make a map file that maps protein accessions to nucleotide accessions
  $cmd = "perl $scripts_dir/make-map.pl $info_file > $map_file";
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
  my @ere_opt_A  = ("0.6", "0.7", "0.8", "0.9", "1.0", "1.1", "1.2", "1.3", "1.4");
  my @ere_name_A = ("0p6", "0p7", "0p8", "0p9", "1p0", "1p1", "1p2", "1p3", "1p4");
  my $nere = scalar(@ere_opt_A);

  for(my $i = 0; $i < $nere; $i++) { 
    my $cm_root = $cm_name . "." . $ere_name_A[$i];
    my $cm_file_name = $cm_root . ".vadr.cm";
    my $cmbuild_file_name = $cm_root . ".vadr.cmbuild";
    my $cmbuild_cmd = "cmbuild -F -n $cm_name --noss --ere $ere_opt_A[$i] $cm_file_name $nt_stk_file > $cmbuild_file_name";
    printf CMBUILD ("qsub -N $cm_root -b y -v SGE_FACILITIES -P unified -S /bin/bash -cwd -V -j n -o /dev/null -e $cm_root.err -l m_mem_free=8G,h_rt=2880000,mem_free=8G,h_vmem=8G -m n \"$cmbuild_cmd\"\n");
  }
}
close(CMBUILD);
printf("\nScript to submit $nmdl cmbuild jobs to the farm is in:\n$cmbuild_qsub_file\n");
printf("\nRun that script, wait for all jobs to finish, then run vb-step3.pl\n");



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

