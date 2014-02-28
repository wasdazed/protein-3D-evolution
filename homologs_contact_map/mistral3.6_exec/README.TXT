AUGUST 2009

              MISTRAL: multiple structural alignment
              ======================================


                          NOTICE

The provided executable (available for Linux and Mac-OSX architectures)
is the computational engine implemented in the MISTRAL web-server
discussed in:

"MISTRAL: a tool for energy-based multiple structural alignment of proteins
Cristian Micheletti; Henri Orland
Bioinformatics 2009", vol. 25, p. 2663-2669


The above article should be cited in publications where the provided code is used.

This software is copyrighted and provided free of charge for academic
use and with no liability on the part of the author.


===============================================


Instructions:
-------------

Make sure that the MISTRAL executable, mistral.x, is in your PATH, or
in the local directory where you have a number of PDB structures to
align.

The Mistral alignment (with default parameters) of the structures is
simply carried out by invoking:

mistral.x name_pdb_file1 name_pdb_file2 name_pdb_file3 etc...

For example, assuming that in your current directory you have the
mistral executable and two PDB structures such as 1nh0.pdb and
1aid.pdb, you can do their pairwise alignment by invoking:

./mistral.x 1nh0.pdb 1aid.pdb

The main quantitative results of the alignment are provided on the
screen as follows:

Pair:  1nh0.pdb  1aid.pdb  lengths: 204  198  energy: -586.709  n_matches: 192  RMSD:    1.0  sid: 182  z_score:    8.8  p-value: 7.30e-06

Further details are provided as output files. In particular:

1) alignment_list.dat: contains the list of corresponding amino acids
   returned by the alignment. At the bottom of the file there is a
   SUMMARY providing the RMSD, percentage of sequence identity
   etc. for each pairwise alignment involving the pivot protein (the
   first one in the table) and every other protein.

2) aligned_structures.pdb: is a PDB file containing the aligned
   proteins in their optimal orientation. Proteins are treated as
   "models", see REMARK field to associate the model number to the
   original PDB filename.  In each model, aligned residues are
   distinguished by a non-zero entry in the "B-factor" (or
   temperature) column.

3) The individual aligned proteins, in their optimal orientation, are
   also provided as separate PDB files named as "aligned_XXX", where
   XXX is the original PDB filename. For the specific example that we
   are considering the two files are, therefore, aligned_1nh0.pdb and
   aligned_1aid.pdb .  In each model, aligned residues are
   distinguished by a non-zero entry in the "B-factor" (or
   temperature) column.



Advanced usage
--------------

You can override the default parameters (described in the MISTRAL
article) and default name of output files by invoking mistral with
suitable "switches". A complete list of the options can be obtained by
invoking mistral with no argument, which produces the following
output:

Usage: ./mistral.x pdb_0 pbd_1 ... pdb_n [-en_int_range energy_interaction_range] [-seed_tol seed_tolerance] [-al_tol alignment_tolerance] [-min_seg_len minimum_aligned_segment_length] [-no_seq_inversion] [-oPDB output_alignment_pdb_file] [-olist output_alignment_list_aa_correspondences]

where: 

-en_int_range       is the interaction range (in Angstroms) for the energy
                    function to be minimised. Default value: 8.0 .

-seed_tol           is the seed tolerance (in Angstroms) used to search for
                    matching (alignable) segments in a protein pair. Default value: 3.5 .

-al_tol             is the maximum tolerance (in Angstroms) for the distances of
                    matching (alignable) amino acid in a protein pair. Default value: 4.0 .

-min_seg_len        minimum length (number of amino acids) of the segments of aligned
                    residues in a protein pair. Default value: 4.

-no_seq_inversion   if this option is present, then no matching are
                    considered for segments that have inverse sequence
                    directionality in two proteins.

-oPDB               name of the "multi-model" file containing the aligned
                    proteins. Default: aligned_structures.pdb .

-olist              name of the file containing the table of corresponding amino
                    acids taking part to the protein core. Default: alignment_list.dat .


For example, to align the above mentioned proteins, using an energy
interaction of 10 A, an alignment tolerance of 5 A, with matching
segments involving at lest 10 amino acids, and disallowing sequence
inversions in matching segments, you would invoke the command:


./mistral.x 1nh0.pdb 1aid.pdb -en_int_range 10  -al_tol 5 -min_seg_len 10 -no_seq_inversion 


For further information, or to report bugs, please contact C. Micheletti, michelet@sissa.it .
