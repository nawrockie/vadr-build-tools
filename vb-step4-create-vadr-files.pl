#!/usr/bin/env perl
use strict;

my $usage = "perl vb-step4-create-vadr-files.pl <'auto' or 'all'> <model list> <name of vadr model dir to create> <gene value (use _ for space)> <product value (use _ for space)>\n";
if(scalar(@ARGV) != 5) { die $usage; }

my $version = "0.03";

my $auto_or_all_opt = undef;
my ($auto_or_all_opt, $model_root_file) = (@ARGV);
my $do_auto = 0;
my $do_all  = 0;

my ($auto_or_all_opt, $model_root_file, $vadr_model_dir, $gene, $product) = (@ARGV);
my $root = $model_root_file;
if($root !~ m/\.model\.list$/) { 
  die "ERROR unable to parse model list file name $model_root_file"; 
}
$root =~ s/\.model\.list$//;

if   ($auto_or_all_opt eq "auto") { $do_auto = 1; }
elsif($auto_or_all_opt eq "all")  { $do_all  = 1; }
else { die "ERROR can't parse first commandline arg, should be 'auto' or 'all'"; }

my $out_root = $root; 
if($do_all) { $out_root .= ".all"; }
else        { $out_root .= ".auto"; }

$gene =~ s/\_/ /g;
$product =~ s/\_/ /g;

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
if(! exists($ENV{"VADRBLASTDIR"})) { 
  die "ERROR, the environment variable VADRBLASTDIR is not set";
}
if(! (-d $ENV{"VADRBLASTDIR"})) { 
  die "ERROR, the directory specified by your environment variable VADRBLASTDIR does not exist.\n"; 
}    
if(! exists($ENV{"VADRINFERNALDIR"})) { 
  die "ERROR, the environment variable VADRINFERNALDIR is not set";
}
if(! (-d $ENV{"VADRINFERNALDIR"})) { 
  die "ERROR, the directory specifiedy by your environment variable VADRINFERNALDIR is not set";
}
if(! exists($ENV{"VADRHMMERDIR"})) { 
  die "ERROR, the environment variable VADRHMMERDIR is not set";
}
if(! (-d $ENV{"VADRHMMERDIR"})) { 
  die "ERROR, the directory specified by your environment variable VADRHMMERDIR does not exist.\n"; 
}    

my $scripts_dir = $ENV{"VADRBUILDTOOLSDIR"} . "/scripts";
my $easel_dir = $ENV{"VADREASELDIR"};
my $blast_dir = $ENV{"VADRBLASTDIR"};
my $hmmer_dir = $ENV{"VADRHMMERDIR"};
my $infernal_dir = $ENV{"VADRINFERNALDIR"};

if(! exists($ENV{"VADRINSTALLDIR"})) { die "ERROR the environment variable VADRINSTALLDIR is not set"; }

my $cmpress_path = $infernal_dir . "/cmpress";
my $cmemit_path = $infernal_dir . "/cmemit";
my $hmmpress_path = $hmmer_dir . "/hmmpress";
my $reformat_path = $easel_dir . "/esl-reformat";
my $makeblastdb_path = $blast_dir . "/makeblastdb";

if(! -s $hmmpress_path)  { die "ERROR hmmpress does not exist at $hmmpress_path"; }
if(! -s $cmpress_path)   { die "ERROR cmpress does not exist at $cmpress_path"; }
if(! -s $cmemit_path)    { die "ERROR cmemit does not exist at $cmemit_path"; }
if(! -s $reformat_path)  { die "ERROR reformat does not exist at $reformat_path"; }

my $cmd;
my $line;

# create the output directory
if(-e $vadr_model_dir) { die "ERROR file or dir named $vadr_model_dir already exists, remove it"; }
my $cmd = "mkdir $vadr_model_dir"; 
system("$cmd");

# parse the model_root file
my @mdl_A = ();
my @mdl_info_file_A = ();
open(IN, $model_root_file) || die "ERROR unable to open $model_root_file for reading";
while($line = <IN>) { 
  chomp $line;
  my $mdl = $line;
  $mdl =~ s/^$root\.//;
  push(@mdl_A, $mdl);
  my $mdl_info_file = $root . "." . $mdl . ".info";
  push(@mdl_info_file_A, $mdl_info_file);

  # make sure all the files we require exist and are non-empty
  my $cm_name = $out_root . "." . $mdl;
  my $cm_file_name = $cm_name . ".1p0" . ".vadr.cm";
  my $mdl_aa_fa_file = $root . "." . $mdl . ".aa.fa"; 

  my @hmm_ere_opt_A  = ("0.1", "0.2", "0.3", "0.4", "0.5", "0.59", "0.7", "0.8", "0.9", "1.0", "1.1", "1.2");
  my @hmm_ere_name_A = ("0p1", "0p2", "0p3", "0p4", "0p5", "0p59", "0p7", "0p8", "0p9", "1p0", "1p1", "1p2");
  my $nhmm_ere = scalar(@hmm_ere_opt_A);

  my @hmm_file_name_A = ();
  my $i;
  for($i = 0; $i < $nhmm_ere; $i++) {
    my $hmm_root = $out_root . "." . $mdl . "." . $hmm_ere_name_A[$i];
    my $hmm_file_name = $hmm_root . ".vadr.hmm";
    push(@hmm_file_name_A, $hmm_file_name);
  }

  my @reqd_files_A = ($mdl_info_file, $cm_file_name, $mdl_aa_fa_file, @hmm_file_name_A);
  foreach my $reqd_file (@reqd_files_A) { 
    if(! -e $reqd_file) { die "ERROR required file $reqd_file does not exist. Did you (succesfully) run vb-step2-taxinfo2muscle-qsub.pl?"; }
    if(! -s $reqd_file) { die "ERROR required file $reqd_file exists but is empty. Did you (succesfully) run vb-step2-taxinfo2muscle-qsub.pl?"; }
  }
}

# for each model: 
my $nmdl = scalar(@mdl_A);
my $m;

my $vadr_minfo_file = $root . ".minfo";
open(MINFO, ">", $vadr_minfo_file) || die "ERROR unable to open $vadr_minfo_file for writing";

for($m = 0; $m < $nmdl; $m++) { 
  my $mdl = $mdl_A[$m];
  my $tt;
  my $subgroup;
  my $group = $root;
  if($mdl =~ m/^tt(\d+)\.(\S+)$/) { 
    ($tt, $subgroup) = ($1, $2);
  }
  else { 
    die "ERROR unable to parse mdl $mdl\n";
  }
  my $mdl_info_file = $mdl_info_file_A[$m];
  my $mdl_aa_fa_file = $root . "." . $mdl . ".aa.fa"; 

  my $cm_and_hmm_name  = $out_root . "." . $mdl;
  my $cm_file_name  = $cm_and_hmm_name . ".1p0" . ".vadr.cm";
  my $hmm_file_name = $cm_and_hmm_name . "vadr.hmm";

  my $blast_db_dst_file = $root . "." . $mdl . ".vadr.protein.fa";

  # make the vadr .minfo file
  my $clen = `grep CLEN $cm_file_name`;
  if($clen =~ s/^CLEN\s+//) {
    chomp $clen;
  }
  print MINFO ("MODEL $root.$mdl group:\"$group\" subgroup:\"$subgroup\" transl_table:\"$tt\" length:\"$clen\" blastdb:\"$blast_db_dst_file\" hmmfile:\"$hmm_file_name\"\n");
  print MINFO ("FEATURE $root.$mdl type:\"CDS\" coords:\"1..$clen\:+\" gene:\"$gene\" product:\"$product\"\n");

  my $cmd = "cp $mdl_aa_fa_file $blast_db_dst_file";
  system("$cmd");

  # make the blast db
  my $cmd = "$blast_dir/makeblastdb -in $blast_db_dst_file -dbtype prot > /dev/null";
  system("$cmd");

  # move to the new dir
  my $cmd = "mv $blast_db_dst_file* $vadr_model_dir";
  system("$cmd");

  # copy alignments to the new dir
  my $nt_stk_src_file = $out_root . "." . $mdl . ".nt.stk";
  my $aa_stk_src_file = $out_root . "." . $mdl . ".hmmbuild.aa.stk";
  my $nt_stk_dst_file = $root . "." . $mdl . ".nt.stk";
  my $aa_stk_dst_file = $root . "." . $mdl . ".hmmbuild.aa.stk";
  my $cmd = "cp $nt_stk_src_file $vadr_model_dir/$nt_stk_dst_file";
  system("$cmd");
  my $cmd = "cp $aa_stk_src_file $vadr_model_dir/$aa_stk_dst_file";
  system("$cmd");
}
close(MINFO);
my $cmd = "mv $vadr_minfo_file $vadr_model_dir";
system("$cmd");

# concatenate all the CM files together and press them
# cmemit the consensus sequence, convert to dna and make a blastn db 
my @cm_ere_opt_A  = ("1.0");
my @cm_ere_name_A = ("1p0");
my $ncm_ere = scalar(@cm_ere_opt_A);
for(my $i = 0; $i < $ncm_ere; $i++) { 
  my $cat_cmd = "cat ";
  my $cm_ere = $cm_ere_name_A[$i];
  #my $big_cm_file_name = $root . "." . $cm_ere . ".cm";
  #my $blastn_db_file_name = $root . "." . $cm_ere . ".fa";
  my $big_cm_file_name = $root . ".cm";
  my $blastn_db_file_name = $root . ".fa";
  for($m = 0; $m < $nmdl; $m++) { 
    my $cm_file_name = $out_root . "." . $mdl_A[$m] . "." . $cm_ere . ".vadr.cm";
    $cat_cmd .= " " . $cm_file_name;
  }
  $cat_cmd .= " > $big_cm_file_name";
  RunCommand($cat_cmd, 1);

  $cmd = "$cmpress_path $big_cm_file_name > /dev/null";
  RunCommand($cmd, 1);

  $cmd = "$cmemit_path -c $big_cm_file_name | sed 's/-hmmconsensus//' | sed 's/-cmconsensus//' | $reformat_path -d fasta - > $blastn_db_file_name";
  RunCommand($cmd, 1);

  # run makeblastdb
  $cmd = "$makeblastdb_path -in $blastn_db_file_name -dbtype nucl";
  RunCommand($cmd, 1);
  
  # move to the new dir
  $cmd = "mv $big_cm_file_name* $vadr_model_dir";
  RunCommand($cmd, 1);

  $cmd = "mv $blastn_db_file_name* $vadr_model_dir";
  RunCommand($cmd, 1);

}

# concatenate all the HMM files together and press them
my @hmm_ere_opt_A  = ("0.1", "0.2", "0.3", "0.4", "0.5", "0.59", "0.7", "0.8", "0.9", "1.0", "1.1", "1.2");
my @hmm_ere_name_A = ("0p1", "0p2", "0p3", "0p4", "0p5", "0p59", "0p7", "0p8", "0p9", "1p0", "1p1", "1p2");
my $nhmm_ere = scalar(@hmm_ere_opt_A);
for(my $i = 0; $i < $nhmm_ere; $i++) { 
  my $cat_cmd = "cat ";
  my $hmm_ere = $hmm_ere_name_A[$i];
  my $big_hmm_file_name = $root . "." . $hmm_ere . ".hmm";
  for($m = 0; $m < $nmdl; $m++) { 
    my $hmm_file_name = $out_root . "." . $mdl_A[$m] . "." . $hmm_ere . ".vadr.hmm";
    $cat_cmd .= " " . $hmm_file_name;
  }
  $cat_cmd .= " > $big_hmm_file_name";
  RunCommand($cat_cmd, 1);

  $cmd = "$hmmpress_path $big_hmm_file_name > /dev/null";
  RunCommand($cmd, 1);

  # move to the new dir
  $cmd = "mv $big_hmm_file_name* $vadr_model_dir";
  RunCommand($cmd, 1);
}

# tar and gzip
my $tarball = $vadr_model_dir . ".tar";
$cmd = "tar -cf $tarball $vadr_model_dir";
RunCommand($cmd, 1);

$cmd = "gzip $tarball";
RunCommand($cmd, 1);

printf("\nTarball with models is here:\n$tarball.gz\n\n");



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

