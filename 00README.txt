EPN, Mon Sep  9 14:23:33 2019
vadr-build-tools v0.02
https://github.com/nawrockie/vadr-build-tools

Outline of vadr-model-build procedure for marker genes.

============================
DOWNLOADING VADR-BUILD-TOOLS
============================
The first step is to clone the github repo that includes the scripts
for building a VADR database with this command:

git clone https://github.com/nawrockie/vadr-build-tools.git

Then move into the directory that gets created called 'vadr-build-tools'.

=============================
SETTING ENVIRONMENT VARIABLES
=============================
Before using these scripts you need to be set up to run VADR. Email
eric.nawrocki@nih.gov for instructions. 

AND

you need to set two additional environment variables in your .bashrc or
.chsrc file:

If you are using the bash shell, add the following
lines to the '.bashrc' file in your home directory:

export VADRBUILDTOOLSDIR=<path to current directory (created by git clone command above)>
export VADRINFERNALDEVDIR=/panfs/pan1/dnaorg/virseqannot/code/vadr-install/infernal/binaries
export VADRHMMERDIR=~nawrocke/src/hmmer-3.3/src

After adding the export lines to your .bashrc file, source that file
to update your current environment with the command:

source ~/.bashrc

--
If you are using the C shell, add the following
lines to the '.cshrc' file in your home directory:

setenv VADRBUILDTOOLSDIR <path to current directory (created by git clone command above)>
setenv VADRINFERNALDEVDIR /panfs/pan1/dnaorg/virseqannot/code/vadr-install/infernal/binaries
setenv VADRHMMERDIR ~nawrocke/src/hmmer-3.3/src

After adding the setenv lines to your .cshrc file, source that file
to update your current environment with the command:

source ~/.cshrc

(To determine which shell you use, type: 'echo $SHELL')


------------------------------------------------------------

1. Collect set of possibly trustable sequences.
   cox1: Susan's 7K manually curated dataset
   cytb: RefSeq query (9175 seqs)
         (cytb [ti] OR cytochrome b [ti]) NOT WGS [filter] NOT
         chromosome [ti] NOT supercontig [ti] NOT contig [ti] NOT
         scaffold [ti] NOT linkage [ti] NOT mRNA [filter] AND
         "mitochondrion"[Filter] NOT plastid [filter] NOT chloroplast
         [filter] NOT plastid [ti] NOT chloroplast [ti] AND RefSeq
         [filter] NOT srcdb pdb[prop] NOT uncultured NOT unverified

------------------------------------------------------------

2. Run vb-step1-fasta2taxinfo.pl to get taxonomic information and
   translation table information for all sequences.

   Example command: 
   $ perl $VADRBUILDTOOLSDIR/vb-step1-fasta2taxinfo.pl cytb.9175.fa cytb > cytb.step1.out

------------------------------------------------------------

3. Create a .tax file which specifies the resolution of taxonomic
   level you want per translation table when creating your VADR
   models.

   First, look at how many translation tables there were in your
   input:
   $ cat cytb.tt.list
   2
   5
   4
   9
   1
   33
   13
   14
   22
   3
   11
   24

   Then for each translation table, look at how many taxonomic groups
   there are at different taxonomic levels (an example of a taxonomic
   level is 'phyla', which is level 3).

   Example, list groups and size of groups for tax level 3 (phyla) for translation table 5

   $ cat cytb.tt5.info | awk '{ print $4 }' | sort | awk -F\; '{ printf("%s;%s;%s;\n", $1, $2, $3); }' | sort | uniq -c
      5 taxonomy:Eukaryota;Metazoa;Chaetognatha;
      9 taxonomy:Eukaryota;Metazoa;Chordata;
   2608 taxonomy:Eukaryota;Metazoa;Ecdysozoa;
      2 taxonomy:Eukaryota;Metazoa;Gnathostomulida;
      2 taxonomy:Eukaryota;Metazoa;Hemichordata;
    572 taxonomy:Eukaryota;Metazoa;Lophotrochozoa;
      1 taxonomy:Eukaryota;Metazoa;Platyhelminthes;
      1 taxonomy:Eukaryota;Metazoa;Porifera;
     10 taxonomy:Eukaryota;Metazoa;Xenacoelomorpha;

   To look at level 2, use this command
   $ cat cytb.tt5.info | awk '{ print $4 }' | sort | awk -F\; '{ printf("%s;%s;\n", $1, $2); }' | sort | uniq -c
   3207 taxonomy:Eukaryota;Metazoa;

   Pick the level you want, and then update a file called 'cytb.tax'
   with this information.

   The 'cytb.tax' needs 1 line per translation table, with exactly 3
   tokens, separated by whitespace. Lines that start with "#" are
   comment lines and are ignored.

   The first token is the translation table, the second is the
   taxonomic level you want for that translation table, and the third
   token is the prefix you want to put before the taxonomic group
   name, usually this will be NONE for no prefix, but for vertebrate
   classes, you may want this to be 'vertebrata-' for example.

   So for tt5 we would add this line:
   5 3 NONE

   A file with the levels I chose for all translation tables for cytb
   is:
   $ cat tax-files/cytb.tax
   # translation_table taxonomy_level prefix
   # tt2: level 7 (mammalia, etc.) 
   2 7 vertebrata-
   # tt5: level 3 (porifera, etc.)
   5 3 NONE
   # tt4: level 3 (cnidaria, ctenophora, placozoa, porifera...
   4 3 NONE
   # tt9: level 3 (platy, hemichordate, echinodermata)
   9 3 NONE
   # tt1: level 2 (Amoebozoa, Viridplantae)
   1 2 NONE
   # tt13: level 4 (tunicata ONLY)
   13 4 NONE
   # tt3: level 3 (dikarya ONLY)
   3 3 NONE
   # tt22: level 2 (viridiplantae)
   22 2 NONE
   # tt11: level 2 (jakobida, nucleariideaandFonticulagroup)
   11 2 NONE
   # tt14: level 3 (platyhelminthes ONLY)
   14 3 NONE
   # tt33: level 3 (hemichordate ONLY)
   33 3 NONE
   # tt24: level 3 (hemichordate ONLY)
   24 3 NONE

------------------------------------------------------------

4. When you are finished creating the .tax file, run
   vb-step2-taxinfo2muscle-qsub.pl to create a script that will submit
   jobs to the compute farm to generate protein sequence alignments
   using muscle or each taxonomic group/translation table set.

   Example command:
   $ perl $VADRBUILDTOOLSDIR/vb-step2-taxinfo2muscle-qsub.pl cytb.9175.fa tax-files/cytb.tax cytb > step2.out

   When that finishes, look at the end of step2.out:
   
   $ tail -n 5 step2.out
   ~~~~~~~~~~~~~~~~~~
   Script to submit 49 muscle jobs to the farm is in:
   cytb.muscle.qsub

   Run that script, wait for all jobs to finish, then run:
   $VADRBUILDTOOLSDIR/vb-step3-muscle-alns2cmbuild-qsub.pl cytb.model.list
   ~~~~~~~~~~~~~~~~~~

   Run the script that submits the muscle jobs:
   $ sh cytb.muscle.qsub 
   
   Then wait for all of those jobs to finish. They are finished when
   the command 'qstat' returns no text.

   Here is an example of a 'qstat' command output when there is still
   a single muscle alignment job running:
   $ qstat
   job-ID     prior   name       user         state submit/start at     queue                          jclass                         slots ja-task-ID 
   ------------------------------------------------------------------------------------------------------------------------------------------------
   226980328 0.54677 tt5.ecdyso nawrocke     r     09/10/2019 11:37:28 unified@sge531.be-md.ncbi.nlm.                                    1        

------------------------------------------------------------

5. When the muscle jobs are all finished. Run
   vb-step3-muscle-alns2cmbuild-qsub.pl to create a script that
   will submit the 'cmbuild' jobs to create the CM files for VADR.
   
   Example command:
   $ perl $VADRBUILDTOOLSDIR/vb-step3-muscle-alns2cmbuild-hmmbuild-qsub.pl auto cytb.model.list > step3.out 

   You may see some output like this:
   SKIPPING nt_len < 3*M     (nt NC_012447.1/3518-4638: 1121 aa YP_002790679.1: 374)
   SKIPPING nt_len < 3*M     (nt NC_035825.1/1286-1: 1286 aa YP_009424481.1: 429)

   Indicating that those sequences were not included in the nucleotide
   alignments. 

   When that finishes, look at the end of step3.out:
   
   $ tail -n 5 step3.out
   ~~~~~~~~~~~~~~~~~~
   Script to submit 49 cmbuild and hmmbuild jobs to the farm is in:
   cytb.build.qsub

   Run that script, wait for all jobs to finish, then run:
   perl $VADRBUILDTOOLSDIR/vb-step4-create-vadr-files.pl auto cytb.model.list <name of vadr model dir to create> <gene value (use _ for space) <product value (use _ for space)>
   ~~~~~~~~~~~~~~~~~~

   Run the script that submits the cmbuild jobs:
   $ sh cytb.cmbuild.qsub 
   
   Then wait for all of those jobs to finish. They are finished when
   the command 'qstat' returns no text.

------------------------------------------------------------

6. When the cmbuild jobs are all finished. Run
   vb-step4-create-vadr-files.pl
   to make the final directory with all the files required to run
   VADR: 

   Example command:
   $ perl $VADRBUILDTOOLSDIR/vb-step4-create-vadr-files.pl auto cytb.model.list vadr-cytb-i0-models-0.991.1-dir CYTB cytochrome_b > step4.out 
   
   This will create a directory called vadr-cytb-i0-models-0.991.1-dir
   with all the files you need to run VADR on cytb sequences. 

------------------------------------------------------------

7. Test run of v-annotate.pl:

   Pick one of the sequence files that ends with .nt.fa to use for a
   test run of v-annotate.pl to make sure the model files were created
   correctly. You probably want to pick one of the .nt.fa files with
   not too many sequences (less than 50) so the command below doesn't
   take too long. (To find out how many sequences are in a file:
   $VADREASELDIR/esl-seqstat <fasta-file>.)

   For the cytb example, we'll use cytb.tt5.chordata.nt.fa which has
   about 10 sequences:

   Example command:
   $ v-annotate.pl -f --keep -i vadr-cytb-i0-models-0.991.1-dir/vadr.cytb.minfo -b vadr-cytb-i0-models-0.991.1-dir -m vadr-cytb-i0-models-0.991.1-dir/vadr.cytb.1p0.cm --xlongest cytb.tt5.chordata.nt.fa va-test
   
------------------------------------------------------------
   
8. Run all i0 sequences used to build the i0 models against the i0 models:

   Sequences that FAIL should probably be removed for the iteration 1
   (i1) models.

   You can use the vb-i0-qsub.pl script to create a script that will
   submit these jobs to the farm: 

   Example command:
   $ perl $VADRBUILDTOOLSDIR/vb-i0-qsub.pl cytb.model.list vadr-cytb-i0-models-0.991.1-dir > i0.qsub.sh

   Then 
   $ sh i0.qsub.sh

------------------------------------------------------------

