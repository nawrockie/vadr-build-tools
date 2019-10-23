#!/usr/bin/env perl
use strict;

my $usage = "perl vb-protein-to-nt.pl <protein fasta file> <output root <root> (output nucleotide fasta file will be <root>.nt.fa)>\n";
if(scalar(@ARGV) != 2) { die $usage; }

my ($fa_file, $root) = (@ARGV);
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

#########################################################
# convert the seqs in the fasta file to short names (acc.version)
my $short_fa_file = $root . ".aa.short.fa";
$cmd = "perl $scripts_dir/fasta-ncbi-idfetch-name-long-to-short.pl $fa_file > $short_fa_file";
RunCommand($cmd, 1);
#########################################################

#########################################################
# index the short fa file
$cmd = "$easel_dir/esl-sfetch --index $short_fa_file > /dev/null";
RunCommand($cmd, 1);
#########################################################

#########################################################
# make list of seqs (short names)
my $short_name_file = $root . ".aa.short.list";
$cmd = "$easel_dir/esl-seqstat -a $short_fa_file | grep ^\= | awk '{ print \$2 }' > $short_name_file";
RunCommand($cmd, 1);
# determine count of seqs
my $nseq = `wc -l $short_name_file`;
chomp $nseq;
#########################################################

#########################################################
# get info on the protein seqs by looking them up with eutils
my $info_file = $root . ".info";
$cmd = "perl $scripts_dir/lookup-prot-acc.pl $short_name_file > $info_file";
RunCommand($cmd, 1);
#########################################################

#########################################################
# HOPEFULLY TEMPORARY
# remove any sequence with a join 
#my $info2_file = $root . ".info2.txt";
#$cmd = "grep -v join $info_file > $info2_file";
#RunCommand($cmd, 1);
#my $nseq2 = `wc -l $info2_file`;
#chomp $nseq2;
#########################################################

#########################################################
# fetch info on codon_start using edirect for nucleotide accessions
my $info2_file = $root . ".info2.txt";
#$cmd = "perl $scripts_dir/lookup-nt-acc.pl $info_file > $info2_file";
$cmd = "perl $scripts_dir/lookup-nt-acc-join.pl $info_file > $info2_file";
RunCommand($cmd, 1);
#########################################################

my $source_list       = $root . ".source.list";
my $source_nt_fa_file = $root . ".source.nt.fa";
my $nt_fa_file        = $root . ".nt.fa";
my $map_file          = $root . ".map";

# create list of 'source' nt accessions
#$cmd = "perl $scripts_dir/parse-by-source.pl $info2_file | sort | uniq > $source_list";
$cmd = "perl $scripts_dir/parse-by-source-join.pl $info2_file | sort | uniq > $source_list";
RunCommand($cmd, 1);

# fetch the source sequences
$cmd = "/netopt/ncbi_tools64/bin/idfetch -t 5 -c 1 -G $source_list | perl $scripts_dir/fasta-ncbi-idfetch-name-long-to-short.pl > $source_nt_fa_file";
RunCommand($cmd, 1);

# index the source sequence files
$cmd = "$easel_dir/esl-sfetch --index $source_nt_fa_file";
RunCommand($cmd, 1);
  
# fetch the subseqs that encode the proteins from their sources
#$cmd = "perl $scripts_dir/parse-to-sfetch.pl $info_file | $easel_dir/esl-sfetch -Cf $source_nt_fa_file - > $nt_fa_file";
$cmd = "perl $scripts_dir/fetch-given-info.pl $info_file $source_nt_fa_file $nt_fa_file";
RunCommand($cmd, 1);

exit 0;

# make a map file that maps protein accessions to nucleotide accessions
$cmd = "perl $scripts_dir/make-map.pl $info_file > $map_file";
RunCommand($cmd, 1);

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

