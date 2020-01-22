# GFF-to-GenBank
Convert GFF file to GenBank file format while extracting the sequences between the annotated regions

Original script written by Brad Chapman and Hsiao Yi - https://github.com/chapmanb/bcbb/blob/master/gff/Scripts/gff/gff_to_genbank.py

Edited by Kartik Chundru to crop the FASTA sequence between the first and last annotated regions.

Usage:
    gff_to_genbank_edit.py \<GFF annotation file\>  \<FASTA sequence file\>
  
  Requires:
  * BioPython - https://biopython.org/
  * BCBio - https://github.com/chapmanb/bcbb
