EPN, Tue Jun 11 2019 

Creating VADR models from an unaligned fasta file of homologous
protein sequences:

------------------------
Step 1. Create intermediate files and a qsub script to submit MUSCLE
alignment jobs to the compute farm:

Usage:
> perl vb-step1.pl <fastafile> <tax-split-file> <output-root e.g. cox1>

Example usage:
> perl vb-step1.pl cox1seqv3 cox1.tax cox1

For an example of the tax-split-file see cox1.tax
That file may need to be manually updated to include taxonomic groups
that cover all the taxonomic groups in <fastafile>.

------------------------
Step 2. Create nucleotide alignments and intermediate files and a qsub
script to submit cmbuild jobs to the compute farm:

Usage:
> perl vb-step2.pl <output-root>.model.list

Example usage:
> perl vb-step2.pl cox1.model.list 

------------------------
Step 3. Create the VADR CM, model info and blast db files: 

Usage:
> perl vb-step3.pl <output-root>.model.list <name of output directory to create> <gene name> <product name>

Example usage:
> perl vb-step3.pl cox1.model.list vadr-cox1-models-0.971.1-dir COX1 cytochrome_c_oxidase_subunit_I

--------------------------

Step 3 will create a file: vadr-cox1-models-0.98.1-dir.tar.gz 
which when unpacked will create a directory called
vadr-cox1-models-0.98.1-dir that includes the models for use with
v-annotate.pl

To run v-annotate.pl:
   > v-annotate.pl \
     -f \ 
     -i <PATH-TO-VADR-COX1-MODELS>/vadr.cox1.minfo \
     -b <PATH-TO-VADR-COX1-MODELS> \ 
     -m <PATH-TO-VADR-COX1-MODELS>/vadr.cox1.1p0.cm \
     <PATH-TO-INPUT-FASTA-FILE> \
     <PATH-TO-DESIRED-OUTPUT-DIR>

If you add the '-p' option, v-annotate.pl will run in 'parallel' mode
and split up the input fasta files into chunks and run across multiple
CPUs on the compute farm.

----------------------------
Relevant output files:
   suffix           description
   ------           -----------
   .vadr.pass.list  list of passing sequences (those with zero errors)
   .vadr.fail.list  list of failing sequences (those with at least one error)
   .vadr.alt.list   list of all errors (alerts) that caused a seq to fail
   .vadr.pass.ft    feature table for passing seqs
   .vadr.fail.ft    feature table for failing seqs
   .vadr.ftr.tbl    one-line-per-seq file with feature coords,
                    pass/fail, and list of alerts
   .vadr.*.tbl      other easy to parse one-line-per files with
                    stats on alerts (errors), etc.


