#!/usr/bin/env python
"""Convert a GFF and associated FASTA file into GenBank format.

Original script written by Brad Chapman and Hsiao Yi - https://github.com/chapmanb/bcbb/blob/master/gff/Scripts/gff/gff_to_genbank.py
Edited by Kartik Chundru to crop the FASTA sequence between the first and last annotated regions.

Usage:
    gff_to_genbank.py <GFF annotation file> <FASTA sequence file>
"""
from __future__ import print_function

import sys
import os

from Bio import SeqIO
from Bio import Seq
from Bio import SeqFeature

from BCBio import GFF

def main(gff_file, fasta_file):
    out_file = "%s.gb" % os.path.splitext(gff_file)[0]
    fasta_input = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    gff_iter = GFF.parse(gff_file, fasta_input)
    SeqIO.write(_check_gff(_fix_ncbi_id(_extract_regions(gff_iter))), out_file, "genbank")

def _fix_ncbi_id(fasta_iter):
    """GenBank identifiers can only be 16 characters; try to shorten NCBI.
    """
    for rec in fasta_iter:
        if len(rec.name) > 16 and rec.name.find("|") > 0:
            new_id = [x for x in rec.name.split("|") if x][-1]
            print("Warning: shortening NCBI name %s to %s" % (rec.id, new_id))
            rec.id = new_id
            rec.name = new_id
        """ Edited by KC Jan 2020. The following loop shortens the feature name length. as opposed to the fragment name length above.
        """
        for i in range(len(rec.features)):
            if len(rec.features[i].type)>15:
                rec.features[i].type=rec.features[i].type[0:15]
        yield rec

def _check_gff(gff_iterator):
    """Check GFF files before feeding to SeqIO to be sure they have sequences.
    """
    for rec in gff_iterator:
        if isinstance(rec.seq, Seq.UnknownSeq):
            print("Warning: FASTA sequence not found for '%s' in GFF file" % (
                    rec.id))
        yield _flatten_features(rec)

def _extract_regions(gff_iterator):
    """Function added by KC Jan 2020. This Extracts regions from the first annotated position to the last annotated position, and updates the locations to correspond to the location in the sequence.
    """
    for rec in gff_iterator:
        pos=[]
        loc=min([i.location.start for i in rec.features])
        endloc=max([i.location.end for i in rec.features])
        for i in range(len(rec.features)):
            pos+=range(int(rec.features[i].location.start),int(rec.features[i].location.end))
            rec.features[i].location=SeqFeature.FeatureLocation(SeqFeature.ExactPosition(rec.features[i].location.start-loc),SeqFeature.ExactPosition(rec.features[i].location.end-loc), strand=rec.features[i].strand)
            for j in range(len(rec.features[i].sub_features)):
                rec.features[i].sub_features[j].location=SeqFeature.FeatureLocation(SeqFeature.ExactPosition(rec.features[i].sub_features[j].location.start-loc),SeqFeature.ExactPosition(rec.features[i].sub_features[j].location.end-loc), strand=rec.features[i].sub_features[j].strand)
        rec.seq=rec.seq[loc:endloc]
        rec.annotations["molecule_type"] = "DNA"
        yield rec

def _flatten_features(rec):
    """Make sub_features in an input rec flat for output.

    GenBank does not handle nested features, so we want to make
    everything top level.
    """
    out = []
    for f in rec.features:
        cur = [f]
        while len(cur) > 0:
            nextf = []
            for curf in cur:
                out.append(curf)
                if len(curf.sub_features) > 0:
                    nextf.extend(curf.sub_features)
            cur = nextf
    rec.features = out
    return rec

if __name__ == "__main__":
    main(*sys.argv[1:])
