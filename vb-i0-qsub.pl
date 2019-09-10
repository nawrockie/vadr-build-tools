#!/usr/bin/env perl
use strict;

my $usage = "perl vb-i0-qsub.pl <model list>  <name of vadr model dir to use>\n";
if(scalar(@ARGV) != 2) { die $usage; }

my ($model_root_file, $vadr_model_dir) = (@ARGV);

my $root = $model_root_file;
if($root !~ m/\.model\.list$/) { 
  die "ERROR unable to parse model list file name $model_root_file"; 
}
$root =~ s/\.model\.list$//;

my $blast_db_dir    = $vadr_model_dir;
my $vadr_minfo_file = $vadr_model_dir . "/vadr." . $root . ".minfo";
my $vadr_cm_file    = $vadr_model_dir . "/vadr." . $root . ".1p0.cm";

# parse the model_root file
my @mdl_A = ();
my @mdl_info_file_A = ();
open(IN, $model_root_file) || die "ERROR unable to open $model_root_file for reading";
while(my $line = <IN>) { 
  chomp $line;
  my $fafile = $line . ".nt.fa";
  my $outdir = "va-i0-$line";
  my $outfile = "va-i0-$line.out";
  my $jobname= "i0." . $line;
  my $errfile= "va-i0-$line.err";
  my $cmd = "v-annotate.pl -f -p --keep -i $vadr_minfo_file -b $blast_db_dir -m $vadr_cm_file $fafile $outdir > $outfile";
  printf("qsub -N $jobname -b y -v SGE_FACILITIES -P unified -S /bin/bash -cwd -V -j n -o /dev/null -e $errfile -l m_mem_free=8G,h_rt=2880000,mem_free=8G,h_vmem=16G -m n \"$cmd\"\n");
}
close(IN);
