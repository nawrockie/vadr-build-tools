#!/usr/bin/env perl
use strict;

my $usage = "perl vb-step3.pl <model list> <name of vadr model dir to create> <gene value (use _ for space)> <product value (use _ for space)>\n";
if(scalar(@ARGV) != 4) { die $usage; }

my ($model_root_file, $vadr_model_dir, $gene, $product) = (@ARGV);
my $root = $model_root_file;
if($root !~ m/\.model\.list$/) { 
  die "ERROR unable to parse model list file name $model_root_file"; 
}
$root =~ s/\.model\.list$//;

$gene =~ s/\_/ /g;
$product =~ s/\_/ /g;

my $scripts_dir = "./scripts";
my $easel_dir = "/usr/local/infernal/1.1.2/bin";
if(! exists($ENV{"VADRINSTALLDIR"})) { die "ERROR the environment variable VADRINSTALLDIR is not set"; }
my $cmpress_path = $ENV{"VADRINSTALLDIR"} . "/infernal-dev/src/cmpress";
if(! -s $cmpress_path) { die "ERROR cmpress does not exist at $cmpress_path"; }

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
  $line =~ s/^$root\.//;
  push(@mdl_A, $line);
  push(@mdl_info_file_A, $root . "." . $line . ".info");
}

# for each model: 
my $nmdl = scalar(@mdl_A);
my $m;
my @mdl_aa_fa_file_A = ();
my @mdl_aa_aln_file_A = ();

my $vadr_minfo_file = "vadr." . $root . ".minfo";
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
  my $mdl_aa_aln_file = $root . "." . $mdl . ".aa.afa"; 
  push(@mdl_aa_fa_file_A, $mdl_aa_fa_file);
  push(@mdl_aa_aln_file_A, $mdl_aa_aln_file);

  my $cm_name = $root . "." . $mdl;
  my $cm_file_name = $cm_name . ".1p0" . ".vadr.cm";

  my $blast_db_dst_file = $root . "." . $mdl . ".vadr.protein.fa";

  # make the vadr .minfo file
  my $clen = `grep CLEN $cm_file_name`;
  if($clen =~ s/^CLEN\s+//) {
    chomp $clen;
  }
  print MINFO ("MODEL $root.$mdl group:\"$group\" subgroup:\"$subgroup\" transl_table:\"$tt\" length:\"$clen\" blastdb:\"$blast_db_dst_file\"\n");
  print MINFO ("FEATURE $root.$mdl type:\"CDS\" coords:\"1..$clen\:+\" gene:\"$gene\" product:\"$product\"\n");

  my $cmd = "cp $mdl_aa_fa_file $blast_db_dst_file";
  system("$cmd");

  # make the blast db
  my $cmd = "/usr/bin/makeblastdb -in $blast_db_dst_file -dbtype prot > /dev/null";
  system("$cmd");

  # move to the new dir
  my $cmd = "mv $blast_db_dst_file* $vadr_model_dir";
  system("$cmd");
}
close(MINFO);
my $cmd = "mv $vadr_minfo_file $vadr_model_dir";
system("$cmd");

# concatenate all the CM files together and press them
# make qsub commands for building the models
my @ere_opt_A  = ("0.6", "0.7", "0.8", "0.9", "1.0", "1.1", "1.2", "1.3", "1.4");
my @ere_name_A = ("0p6", "0p7", "0p8", "0p9", "1p0", "1p1", "1p2", "1p3", "1p4");
my $nere = scalar(@ere_opt_A);
for(my $i = 0; $i < $nere; $i++) { 
  my $cat_cmd = "cat ";
  my $ere = $ere_name_A[$i];
  my $big_cm_file_name = "vadr." . $root . "." . $ere . ".cm";
  for($m = 0; $m < $nmdl; $m++) { 
    my $cm_name = 
    my $cm_file_name = $root . "." . $mdl_A[$m] . "." . $ere . ".vadr.cm";
    $cat_cmd .= " " . $cm_file_name;
  }
  $cat_cmd .= " > $big_cm_file_name";
  RunCommand($cat_cmd, 1);

  $cmd = "$cmpress_path $big_cm_file_name > /dev/null";
  RunCommand($cmd, 1);

  # move to the new dir
  $cmd = "mv $big_cm_file_name* $vadr_model_dir";
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

