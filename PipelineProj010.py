################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id:  $
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#################################################################################
"""
==============================================
PipelineProj010.py - custom tasks for proj010
==============================================

:Author: Jethro Johnson
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This module file contains generic functions for use in the script
pipeline_proj010_lncrna.py

"""

import gzip, os, re
import pysam
import random
import pickle
import shutil
import tempfile
import copy
import numpy as np
import collections
import itertools
import pandas as pd
import pandas.rpy.common as com

import rpy2.robjects as robjects
import rpy2.robjects.pandas2ri as pandas2ri
from rpy2.robjects import r as R
from rpy2.robjects import numpy2ri as rpyn
from rpy2.robjects.packages import importr

from scipy import stats
from bx.bbi.bigwig_file import BigWigFile
from pybedtools.contrib.bigwig import bam_to_bigwig

import CGAT.IOTools as IOTools
import CGAT.GTF as GTF
import CGAT.Bed as Bed
import CGATPipelines.Pipeline as P
import CGAT.Experiment as E
import CGAT.Expression as Expression
import CGAT.IndexedGenome as IndexedGenome
import CGATPipelines.PipelineUtilities as PU
from CGATPipelines.Pipeline import cluster_runnable

import PipelineProj010QC as QC

#################################################################################
#################################################################################
#################################################################################
## section: miscellaneous
#################################################################################
# Things that are too clunky to sit in the main pipeline
#################################################################################


def postProcessClosestBed( in_bed, outfile, cage = True ):
    """
    Resolve intervals that are equidistant
    """
    outf = IOTools.openFile( outfile, "w" )
    if cage:
        outf.write( "lnc_contig\t"
                    "lnc_start\t"
                    "lnc_end\t"
                    "lnc_id\t"
                    "lnc_score\t"
                    "lnc_strand\t"
                    "cage_contig\t"
                    "cage_start\t"
                    "cage_end\t"
                    "cage_id\t"
                    "cage_score\t"
                    "cage_strand\t"
                    "cage_thickStart\t"
                    "cage_thickEnd\t"
                    "cage_rgb\t"
                    "distance\n" )
    else:
        outf.write( "lnc_contig\t"
                    "lnc_start\t"
                    "lnc_end\t"
                    "lnc_id\t"
                    "lnc_score\t"
                    "lnc_strand\t"
                    "refcoding_contig\t"
                    "refcoding_start\t"
                    "refcoding_end\t"
                    "refcoding_id\t"
                    "refcoding_score\t"
                    "refcoding_strand\t"
                    "distance\n" )

    def _get_distance( tss, start, end ):
        tss, start, end = map( int, ( tss, start, end ) ) 
        assert end > start, "Error in bedfile"
        if tss > end: 
            distance = tss - end
        elif start > tss:
            distance = start - (tss + 1)
        else: 
            distance = 0
        return distance

    # iterate through bedfile, add each entry to dictionary with gene_id as key
    # for cage intervals resolve closest by selecting interval with closest thickSart
    # for refcoding tss resolve closest by taking first interval
    near_peaks = collections.defaultdict( list ) 
    for line in IOTools.openFile( in_bed ):
        near_peaks[ line.split()[3] ].append( line )

    for key, values in near_peaks.iteritems():
        if len( values ) == 1:
            outf.write( near_peaks[ key ][0] )
            pass
        elif not cage:
            outf.write( near_peaks[ key ][0] )
        else:
            dst = [ x.split()[-1] for x in values ]
            assert len( set( dst ) ) == 1, "Duplicate entries with different distances"
            tss = values[0].split()[1]
            
            closest = None
            closest_dist = None
            dists = []
            for value in enumerate( values ):
                value_start = value[1].split()[12]
                value_end = value[1].split()[13]
                distance = _get_distance( tss, value_start, value_end )
                dists.append( distance )
                if closest_dist:
                    if distance < closest_dist:
                        closest = value[0]
                        closest_dist = distance
                else:
                    closest = value[0]
                    closest_dist = distance
            E.info( "Resolving equidistant peaks for %s: %s"
                    " Keeping %s, which corresponds to %s" %( key, 
                                                              str( dists ), 
                                                              str( closest_dist ), 
                                                              near_peaks[key][closest].split()[9] ) )
            outf.write( near_peaks[ key ][ closest ] )
    outf.close()         

def summarize_eRNAs( lncRNAs, infiles, outfile ):
    """
    Takes sorted set of lncRNAs and makes empty dictionary.
    Takes bedfiles containing gene/transcript, eRNAs/pRNAs for robust/permissive
    /original datasets and outputs a single file that summarizes lncRNA status
    for all datasets.    
    """
    # set up empty dictionary
    summary_dict = {}
    for gtf in GTF.flat_gene_iterator( GTF.iterator( IOTools.openFile( lncRNAs ) ) ):
        gene_id = gtf[0].gene_id
        assert gene_id not in summary_dict.keys(), "Duplicate Gene IDs"
        summary_dict[ gene_id ] = { "gene": "ambiguous", 
                                    "gene_permissive": "ambiguous", 
                                    "gene_robust": "ambiguous",
                                    "transcript": "ambiguous", 
                                    "transcript_permissive": "ambiguous", 
                                    "transcript_robust": "ambiguous" }
    # set up counting dictionary
    count_dict = { "gene": { "eRNA": 0, "pRNA": 0, "conflicted": 0 },
                   "gene_permissive": { "eRNA": 0, "pRNA": 0, "conflicted": 0 },
                   "gene_robust": { "eRNA": 0, "pRNA": 0, "conflicted": 0 },
                   "transcript": { "eRNA": 0, "pRNA": 0, "conflicted": 0 },
                   "transcript_permissive": { "eRNA": 0, "pRNA": 0, "conflicted": 0 },
                   "transcript_robust": { "eRNA": 0, "pRNA": 0, "conflicted": 0 } }

    # iterate through infiles
    for infile in infiles:
        E.info( "Processing file: %s" % os.path.basename( infile ) )
        # specify category based on filename
        col_name = None
        category = None
        inf_name = os.path.basename( infile ).split("_")
        if inf_name[1] in [ "eRNA", "pRNA" ]:
            col_name = inf_name[0]
            category = inf_name[1]
        elif inf_name[1] in [ "permissive", "robust" ]:
            assert inf_name[2] in [ "eRNA", "pRNA" ], "Unrecognized file name"
            col_name = "_".join( inf_name[0:2] )
            category = inf_name[2]
        else: 
            raise Exception( "Unrecognised file name %s" % os.path.basename(infile) )

        # HACK... transcripts have different TSS in bedfile, hence are not merged.
        gene_id_list = [ bed.fields[0] for bed in Bed.iterator( IOTools.openFile( infile ) ) ]
        if col_name.startswith( "gene" ):
            assert len( gene_id_list ) == len( set( gene_id_list ) ), "Ambiguous TSS for genes" 
        gene_id_list = list( set( gene_id_list ) )

        # iterate through gene_ids and add information
        # for each interval to summary_dict
        for gene_id in gene_id_list:
            # E.info( "Processing gene_id %s" % gene_id )
            # it is not possible for field to be conflicted until it has been
            # parsed twice, no entry should be parsed more than twice.
            assert summary_dict[ gene_id ][ col_name ] != "conflict", "Duplicate Gene IDs in eRNA files"
            # reset ambiguous fields, based on current infile
            if summary_dict[ gene_id ][ col_name ] == "ambiguous":
                summary_dict[ gene_id ][ col_name ] = category
                count_dict[ col_name ][ category ] += 1
            # check no field is reassigned the same lncRNA status twice
            elif summary_dict[ gene_id ][ col_name ] in [ "eRNA", "pRNA" ]:
                if summary_dict[ gene_id ][ col_name ] == category:
                    raise ValueError( "lncRNA is categorised twice as %s. "
                                      "Failed on %s, %s for file %s" 
                                      % ( category, 
                                          gene_id, 
                                          col_name, 
                                          os.path.basename( infile ) ) )
                else:
                    # remove count from original classification
                    to_remove = summary_dict[ gene_id ][ col_name ] 
                    count_dict[ col_name ][ to_remove ] -= 1
                    summary_dict[ gene_id ][ col_name ] = "conflicted"
                    count_dict[ col_name ][ "conflicted" ] += 1
            else:
                raise Exception( "Unrecognized field entry %s for gene_id"
                                 " %s, file %s" % ( summary_dict[gene_id][col_name],
                                                    gene_id, 
                                                    col_name ))

    # write outfile
    headers = [ "gene", "gene_permissive", "gene_robust", 
                "transcript", "transcript_permissive", "transcript_robust" ]
    outf = IOTools.openFile( outfile, "w" )
    outf.write( "gene_id\t" + "\t".join( headers ) + "\n" )
    for gene_id, fields in summary_dict.iteritems():
        out_fields = []
        for field in headers:
            out_fields.append( fields[ field ] )
        outf.write( gene_id + "\t" + "\t".join( out_fields ) + "\n" )
    outf.close()

    # write summary outfile
    headers = [ "eRNA", "pRNA", "conflicted" ]
    outf_summary = P.snip( outfile, ".tsv.gz" ) + "_summary.tsv.gz" 
    outf_summary = IOTools.openFile( outf_summary, "w" )
    outf_summary.write("data_set\t" + "\t".join( headers ) + "\n" )
    for data_set, fields in count_dict.iteritems():
        out_fields = []
        for field in headers: 
            out_fields.append( str( fields[ field ] ) )
        outf_summary.write( data_set + "\t" + "\t".join( out_fields ) + "\n" )
    outf_summary.close()


def compareOverlaps(genes_gtf, 
                    tf_bed, 
                    genome_file, 
                    genome_file_ungapped,
                    isochore_file,
                    outf_stub, 
                    annotation="lncRNA" ):
    """
    A wrapper for various commandline statements associated with comparing
    the overlap between a geneset (genes_gtf) and a bed file of transcription
    factor binding sites (tf_bed). The steps are as follows:
    i) slop geneset (1kb upstream) and convert to bedfile
    ii) run diff_bed.py to get overlap statistics
    iii) run GAT to see if tf binding across gene intervals is enriched vs 
        background
    iv) create boolean table detailing which gene_ids have tfbs in gene body
    """

    # create outfiles
    genes_bed_slop = outf_stub + "_" + annotation + "_1kb.bed.gz"
    outf_diff = outf_stub + "_overlap.tsv"
    outf_gat = outf_stub + "_bg_gat.tsv"
    outf_tab = outf_stub + "_coverage.tsv.gz"
    outf_bed = outf_stub + "_coverage.bed.tsv.gz"

    # Create intervals from gtf gene models
    statement = ( "zcat %(genes_gtf)s |"
                  " python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=merge-transcripts"
                  "  --log=%(outf_stub)s.log |"
                  " bedtools slop "
                  "  -i stdin"
                  "  -l 1000"
                  "  -r 0" # necessary to state -r if -l given
                  "  -g %(genome_file)s |"
                  " python %(scriptsdir)s/gff2bed.py"
                  "  --is-gtf"
                  "  --log=%(outf_stub)s.log |"
                  " gzip > %(genes_bed_slop)s" )
    P.run()
    E.info( "Completed slop for sample %s" % genes_gtf )

    # run diff_bed.py to compare overlap between tfbs and gene models
    statement = ( "python %(scriptsdir)s/diff_bed.py"
                  " --pattern-identifier='.*/(.+).bed.gz'"
                  " --log=%(outf_stub)s.log"
                  " %(genes_bed_slop)s %(tf_bed)s"
                  " > %(outf_diff)s" )
    P.run()
    E.info( "Completed diff_bed comparison for sample %s" % genes_gtf )
    
    # run gat to find significance of overlap between tfbs and 
    # geneset vs background
    # gat annotations file requires 4th column to specify annotation
    tmpf = P.getTempFilename("/ifs/scratch")
    statement = ( "zcat %(genes_bed_slop)s |"
                  " awk 'BEGIN {OFS=\"\\t\"} {print $1,$2,$3,\"%(annotation)s\"}'"
                  " > %(tmpf)s" )
    P.run()

    job_options = "-l mem_free=5G"
    statement = ( "gat-run.py"
                  "  --segments=%(tf_bed)s"
                  "  --annotations=%(tmpf)s" # replacement for gene_bed_slop
                  "  --workspace-bed-file=%(genome_file_ungapped)s"
                  "  --isochore-file=%(isochore_file)s"
                  "  --ignore-segment-tracks"
                  "  --truncate-segments-to-workspace"
                  "  --num-samples=10000"
                  "  -v5"
                  "  --log=%(outf_stub)s.log"
                  " > %(outf_gat)s" )
    P.run()
    E.info( "Completed GAT run for sample %s" % genes_gtf )
    os.unlink( tmpf )

    # create table of booleans specifying whether a gene intersects a tfbs
    # run bedtools intersect to retrieve gene ids intersecting tfbs
    statement = ( "bedtools intersect"
                  " -a %(genes_bed_slop)s"
                  " -b %(tf_bed)s"
                  " -wa"
                  " -u |"
                  " gzip > %(outf_bed)s" )
    P.run()

    # create list of gene_ids with tf binding
    bound_genes = []
    for line in IOTools.openFile( outf_bed ):
        gene_id = line.split()[3]
        bound_genes.append( gene_id )

    # write table
    outf = IOTools.openFile( outf_tab, "w" )
    outf.write( "gene_id\t%s\n" % os.path.basename(outf_stub) )
    for gtfs in GTF.flat_gene_iterator( 
            GTF.iterator( IOTools.openFile( genes_gtf ) ) ):
        gene_id = gtfs[0].gene_id
        if gene_id in bound_genes:
            outf.write( gene_id + "\tTRUE\n" )
        else: 
            outf.write( gene_id + "\tFALSE\n" )
    outf.close()       
    E.info( "Created overlap table for sample %s" % genes_gtf )


# extract chromatin de regions from DESeq output in csvdb
def getDEIntervals(database, table_name, direction, fold_change):
    """
    Use PU fetch to return interavals at required thresholds from DESeq results
    table.
    """
    statement = ( "SELECT test_id, control_name, treatment_name, l2fold, pvalue"
                  " FROM %(table_name)s"
                  " WHERE l2fold %(direction)s %(fold_change)s" % locals() )
    df = PU.fetch_DataFrame( statement, database = database )
    
    return  df
        



#################################################################################
#################################################################################
#################################################################################
## section: manipulate bam files
#################################################################################

@cluster_runnable
def normalizeBams(bamfile, normalization_file, outfile):
    # create dictionary of normalization factors. 
    norm_factors = {}
    for line in IOTools.openFile(normalization_file).readlines()[1:]:
        line = line.split()
        norm_factors[line[1]] = float(line[0])

    # get the smallest normalization factor (this is the smallest bamfile). 
    norm_min = min(norm_factors.values())
    # calculate this as a proportion of all normalization factors
    for key, value in norm_factors.iteritems():
        norm_factors[key] = norm_min/value

    # pull out the relevant normalization factor
    bam_id = P.snip(os.path.basename(bamfile), ".bam")
    norm_factor = norm_factors[bam_id]

    print bam_id, norm_factor
#    # iterate through the bamfile
#    pysam_in = pysam.Samfile(bamfile, "rb")
#    pysam_out = pysam.Samfile(outfile, "wb", template=pysam_in)

    # write out the specified proportion of reads
#    for read in pysam_in.fetch():
#        if random.random() <= norm_factor:
#            pysam_out.write(read)

#    pysam_in.close()
#    pysam_out.close()
#    pysam.index(outfile)


def normalize( infile, larger_nreads, outfile, smaller_nreads ):
    threshold = float(smaller_nreads) / float(larger_nreads)

    pysam_in = pysam.Samfile( infile, "rb" )
    pysam_out = pysam.Samfile( outfile, "wb", template = pysam_in )

    for read in pysam_in.fetch():
         if random.random() <= threshold:
             pysam_out.write( read )

    pysam_in.close()
    pysam_out.close()
    pysam.index( outfile )

    return outfile

@cluster_runnable
def normalizeBamfiles( sample_file, input_file, sample_outfile, input_outfile, submit=True ):
    sample_sam = pysam.Samfile( sample_file, "rb" )
    sample_nreads = 0
    for read in sample_sam:
        sample_nreads += 1

    input_sam = pysam.Samfile( input_file, "rb" )
    input_nreads = 0
    for read in input_sam:
        input_nreads += 1

    if input_nreads > sample_nreads: 
        P.info( "%s bam has %s reads, %s bam has %s reads"
                % ( input_file, input_nreads, sample_file, sample_nreads ) )
        P.info( "%s being downsampled to match %s" 
                % ( input_file, sample_file ) )
        input_outfile = normalize( input_file,
                                   input_nreads, 
                                   input_outfile, 
                                   sample_nreads )
        shutil.copyfile( sample_file, sample_outfile )
        pysam.index( sample_outfile )
        
#        return sample_outfile, input_outfile
        
    elif sample_nreads > input_nreads:
        P.info( "%s bam has %s reads, %s bam has %s reads"
                % ( sample_file, sample_nreads, input_file, input_nreads ) )
        P.info( "%s being downsampled to match %s" 
                % ( sample_file, input_file ) )
        sample_outfile = normalize( sample_file, 
                                    sample_nreads, 
                                    sample_outfile, 
                                    input_nreads )
        shutil.copyfile( input_file, input_outfile )
        pysam.index( input_outfile )
        
#        return sample_outfile, input_outfile

    else: 
        E.info ( "WARNING: Both bamfiles are the same size!!" )
        shutil.copyfile( sample_file, sample_outfile )
        pysam.index( sample_outfile )
        shutil.copyfile( input_file, input_outfile  )
        pysam.index( input_outfile )

#        return sample_outfile, input_outfile


# merge bams using samtools merge
def mergeBam( infile_list, outfile ):
    out_stub = P.snip( outfile, ".bam" )
    to_cluster = True
    job_options = "-l mem_free=5G"
    statement = ( "samtools merge - %(infile_list)s"
                  " | samtools sort - %(out_stub)s"
                  " 2>%(outfile)s.log;"
                  " checkpoint;"
                  " samtools index %(outfile)s"
                  " 2>%(outfile)s.bai.log" )
    P.run()


# randomly split reads from bamfile into two bamfiles of approx. equal size
def splitBam( infile, outfiles ):
    # rewrite as a dictionary to take multiple splits
    outfile_00, outfile_01 = outfiles
    pysam_in = pysam.Samfile( infile, "rb" )
    pysam_out_00 = pysam.Samfile( outfile_00, "wb", template = pysam_in )
    pysam_out_01 = pysam.Samfile( outfile_01, "wb", template = pysam_in )

    for read in pysam_in.fetch():
        if random.random() <= 0.5:
            pysam_out_00.write( read )
        else:
            pysam_out_01.write( read )

    pysam_out_00.close()
    pysam_out_01.close()
    pysam.index( outfile_00 )
    pysam.index( outfile_01 )


# randomly split reads from bamfile into a specified number of bamfiles
def multiSplitBam( infiles, outfiles, params ):
    infile = infiles
    outfile_stub = outfiles
    n_outfiles = int(params[0])

    pysam_in = pysam.Samfile( infile, "rb" )

    outfile_handles = []
    outfile_names = []

    # create list of upper bounds for intervals
    intervals = []
    lower = 0
    for i in range( n_outfiles ):
        upper = lower + 1.0/n_outfiles
        intervals.append( upper )
        # add an outfile handle to list of outfile handles
        outf = outfile_stub + "_" + str(i).zfill(2) + ".bam" 
        outfile_names.append( outf )
        outfile_handles.append( pysam.Samfile( outf, "wb", template = pysam_in ) )
        lower = upper 

    # iterate through reads in samfile and write them to an outfile at random
    for read in pysam_in.fetch():
        r_num = random.random()
        for i in range( len( intervals ) ):
            if r_num < intervals[i]:
                outfile_handles[i].write( read ) 
                break
            else: continue

    # close outfiles
    for i in range( n_outfiles ):
        outfile_handles[i].close()

    # index outfiles
    for split_sam in outfile_names:
        pysam.index( split_sam )


# split bam file into separate files for + and - strands
def splitBamByStrand( infile, outfiles ):
    out_pos, out_neg = outfiles
    out_pos = P.snip( out_pos, ".bam" )
    out_neg = P.snip( out_neg, ".bam" )
    tmpf_pos = P.getTempFilename( "/scratch" )
    tmpf_neg = P.getTempFilename( "/scratch" )

    pysam_in = pysam.Samfile( infile, "rb" )
    pysam_out_pos = pysam.Samfile( tmpf_pos, "wb", template = pysam_in )
    pysam_out_neg = pysam.Samfile( tmpf_neg, "wb", template = pysam_in )

    for read in pysam_in.fetch():
        if read.is_read1 and not read.is_reverse:
            pysam_out_pos.write( read )
        elif read.is_read2 and read.is_reverse:
            pysam_out_pos.write( read )
        elif read.is_read1 and read.is_reverse:
            pysam_out_neg.write( read )
        elif read.is_read2 and not read.is_reverse:
            pysam_out_neg.write( read )

    pysam_out_pos.close()
    pysam_out_neg.close()

    pysam.sort( tmpf_pos, out_pos )
    pysam.sort( tmpf_neg, out_neg )

    # check this 
    pysam.index( out_pos + ".bam" )
    pysam.index( out_neg + ".bam" )

    os.unlink( tmpf_pos )
    os.unlink( tmpf_neg )


# converts bamfile to bigwig, via bedgraph
# contigs - list of contig sizes (PARAMS_ANNOTATIONS["interface_contigs"])
# scale provides option for normalizing the bigwigs
def bam2bigwig( bamfile, bigwigfile, contigs, scale=1.0 ):
    tmpf_1 = P.getTempFilename( "./bigwigs" )
    tmpf_2 = P.getTempFilename( "./bigwigs" )

    to_cluster = True
    job_options = "-l mem_free=10G"
    statement = ( "bedtools genomecov"
                  "  -split"
                  "  -scale %(scale)s"
                  "  -bga"
                  "  -ibam %(bamfile)s"
                  "  -g %(contigs)s"
                  "  > %(tmpf_2)s"
                  " bedGraphToBigWig"
                  "  %(tmpf_2)s"
                  "  %(contigs)s"
                  "  %(bigwigfile)s"
                  "  2> %(bigwigfile)s.log" )
    P.run()

    os.unlink( tmpf_1 )
    os.unlink( tmpf_2 )

#                  "  2> %(bigwigfile)s.log;"
#                  " sort  -k1g,1g -k2n,2n"
#                  "  %(tmpf_1)s"
#                  "  > %(tmpf_2)s"
#                  "  2> %(bigwigfile)s.log;"

"""
    # write sam header to a tempfile
    tmpf = P.getTempFilename(".")
    statement = ( "samtools view -H %(infile)s > %(tmpf)s.header" )
    P.run()

    # split bam into two files with original sam header
    bamfile = pysam.Samfile( infile, "rb" )
    nreads = 0
    for entry in bamfile: nreads += 1
    nreads_split = (nreads + 1)/num_split
    to_cluster = True
    statement = ( "samtools view %(infile)s"
                  " | shuf"
                  " split -d -l %(nreads_split)s - %(out_stub)s;"
                  " checkpoint;"

                  " cat %(tmpf)s.header > %(tmpf)s00;"
                  " cat %(out_stub)s00 >> %(tmpf)s00;"
                  " samtools view -bT %(genome_dir)s/%(genome)s.fa %(tmpf)s00 "
                    " > %(out_stub)s00.bam;"
                  " samtools index %(out_stub)s00.bam;"

                  " cat %(tmpf)s.header > %(tmpf)s01;"
                  " cat %(out_stub)s01 >> %(tmpf)s01;"
                  " samtools view -bT %(genome_dir)s/%(genome)s.fa %(tmpf)s01 "
                    " > %(out_stub)s01.bam;"
                  " samtools index %(out_stub)s01.bam;"

                  " rm %(tmpf)s.header %(tmpf)s00 %(tmpf)s01 %(out_stub)s00 %(out_stub)s01" )

    P.run()
"""

#################################################################################
#################################################################################
#################################################################################
## section: manipulate gtf files
#################################################################################
# NB. GTF.iterator object attributes are zero based half open, which is pythonic
#     as opposed to the original gtf which is 1-based inclusive (i.e. closed).
# Bed format is zero-based half open.
# Therefore accessing gtf attributes and writing them to bed format does not 
# require any numeric conversion. 

def defineTranscripts( exon_list ):
    """
    Receives a list of gtf objects pertaining to multiple transcripts as 
    supplied by GTF.flat_gene_iterator.
    Returns an ordered dictionary with transcript_id as key and (start, end, 
    number_of_intervals) as value. Dictionary is ordered by transcript start. 
    """
    transcript_dict = {}

    for exon in exon_list:
        t_id = exon.transcript_id
        if t_id not in transcript_dict:
            transcript_dict[ t_id ] = [ exon.start, exon.end, 1 ] 
        else:
            transcript_dict[ t_id ][2] += 1
            if transcript_dict[ t_id ][0] > exon.start:
                transcript_dict[ t_id ][0] = exon.start
            if transcript_dict[ t_id ][1] < exon.end:
                transcript_dict[ t_id ][1] = exon.end
    
    return collections.OrderedDict( sorted( transcript_dict.items(), 
                                            key = lambda t: t[1][0] ) )


def findTSSFromGTF( infile, outfile ):
    inf = IOTools.openFile( infile, "r" )
    outf = IOTools.openFile( outfile, "w" )

    for exons in GTF.transcript_iterator( GTF.iterator( inf ) ):
        chrom = exons[0].contig
        start = str(exons[0].start)
        end = str(exons[0].start + 1)
        gene_id = exons[0].gene_id
        strand = exons[0].strand
        outf.write( chrom + "\t" 
                    + start + "\t" 
                    + end + "\t" 
                    + gene_id + "\t"
                    + ".\t"
                    + strand + "\n" )
    outf.close()


# splits transcripts in gtf based on a chosen field, using dictionary
#  output suffix as key and matching regex as value
# if retain == true, then all transcripts not matching regex will be output
# if dict contains value "other", failed matches are output using given key as
#  suffix, otherwise failed matches are output using suffix 'unassigned'
def splitGTFByFieldValues( infile, out_stub, field, dictionary, retain = False ):
    inf = IOTools.openFile( infile, "r" )
    suffix = ".gtf" + infile.split( ".gtf" )[1]

    other = "unassigned"
    for key, value in dictionary.iteritems():
        if key.lower() == "other":
            other = dictionary.pop( key ) 
            break                       # can't change dict size during iteration
    other = IOTools.openFile( out_stub + "_" + other + suffix, "w" )

    file_dict = {}
    for key in dictionary.iterkeys():
        file_dict[ key ] = IOTools.openFile( out_stub + "_" + key + suffix, "w" )

    for transcript in GTF.transcript_iterator( GTF.iterator( inf ) ):
        match = False
        for key, value in dictionary.iteritems():
            if re.search( value, transcript[0].__getattribute__( field ) ):
                flat_transcript = IOTools.flatten( transcript )
                for interval in flat_transcript:
                    file_dict[ key ].write( str( interval ) + "\n" )
                match = True
            else: continue
        if retain and not match:
            flat_transcript = IOTools.flatten( transcript )
            for interval in flat_transcript:
                other.write( str( interval ) + "\n" )

    for value in file_dict.itervalues():
        value.close()
    other.close()

    other =  out_stub + "_unassigned" + suffix
    if os.path.exists( other ):
        if not retain:
            os.unlink( other )




# intersect two gtf files
# remove any gene_id in gtf_a that intersects with gtf_b
def intersectGTFs( lncRNA_gtf, 
                   reference_gtf, 
                   filtered_gtf, 
                   rejected_gtf, 
                   control_gtf = None,
                   control_biotypes = None ):
    """
    Intersects two gtf files using bedtools intersect. Removes any transcript for
    which an interval in gtf_a intersects an interval in gtf_b. 
    Has the option of providing a third - control - gtf (plus list of biotypes). 
    If provided, any transcript to be removed is checked to see if it intersects
    with the control. Those that intersect transcripts of specified biotype in 
    control file are rescued.
    This is a hack designed to rescue transcripts that are labelled as protein 
    coding in refseq annotations but as lincRNAs in ensembl annotations.
    """

    # intersect the lncRNA_gtf with the reference_gtf - any locus that intersects 
    # is output to tmpfile
    tmpf = P.getTempFilename( "." )
    statement = ( "bedtools intersect"
                  " -a %(lncRNA_gtf)s"
                  " -b %(reference_gtf)s"
                  " -u"
                  " -s"
                  " > %(tmpf)s" )
    P.run()

    # iterate through tmpfile and add intersecting gene_ids to gene_ids_to_reject
    gene_ids_to_reject = []
    for gtf in GTF.iterator( IOTools.openFile( tmpf ) ):
        gene_ids_to_reject.append( gtf.gene_id )
    gene_ids_to_reject = set( gene_ids_to_reject )
    os.unlink( tmpf )

    # iterate through the lncRNA_gtf, moving gene_ids_to_reject to temp_rejected
    inf = IOTools.openFile( lncRNA_gtf )
    temp_filtered = P.getTempFilename( "." )
    temp_rejected = P.getTempFilename( "." )

    temp_filtered_out = IOTools.openFile( temp_filtered, "w" )
    temp_rejected_out = IOTools.openFile( temp_rejected, "w" )
    for gtf in GTF.flat_gene_iterator( GTF.iterator( inf ) ):
        if gtf[0].gene_id in gene_ids_to_reject:
            for exon in IOTools.flatten( gtf ):
                temp_rejected_out.write( str( exon ) + "\n" )
        else:
            for exon in IOTools.flatten( gtf ):
                temp_filtered_out.write( str( exon ) + "\n" )
    temp_rejected_out.close()

    # if a control gtf is provided, then provide a chance to save temp_rejected
    # i) intersect intervals in temp_rejected with intervals in control_gtf
    # ii) if intersecting control intervals are of a specified biotype, then 
    # intervals in temp_rejected are saved (i.e. written to rescued and filtered )
    if control_gtf:
        assert isinstance( control_biotypes, (list,) ), "need to supply biotypes as list"
        to_rej_gtf = IOTools.openFile( temp_rejected )
        control_gtf = IOTools.openFile( control_gtf )
        def_rejected = P.getTempFilename( "." )
        rescued_gtf = re.sub( "filtered", "rescued", filtered_gtf ) 
        to_rescue = []

        # create index of intervals in control file
        control_index = GTF.readAndIndex( GTF.iterator( control_gtf ) )

        # iterate through the currently rejected transcripts
        for gtf in GTF.transcript_iterator( GTF.iterator( to_rej_gtf ) ):
            for exon in gtf:
                if control_index.contains( exon.contig, exon.start, exon.end ):
                    # check if any intersecting control intervals are on the same strand
                    for interval in control_index.get( exon.contig, exon.start, exon.end ):
                        if exon.strand == interval[2].strand:
                            # for those that are, check biotype
                            if interval[2].source in control_biotypes:
                                to_rescue.append( exon.gene_id )
                            else: 
                                continue
                        else: 
                            continue
                else:
                    continue
        to_rej_gtf.close()
        control_gtf.close()

        # create a set gene_ids to save
        to_rescue = set( to_rescue )
        # iterate through the temp_rejected_gtfs,
        #  if they are in to_save write them to rescued_gtf and to filtered_gtf
        #  if not then write them to def_rejected
        def_rejected_out = IOTools.openFile( def_rejected, "w" )
        rescued_gtf_out = IOTools.openFile( rescued_gtf, "w" )
        for gtf in GTF.flat_gene_iterator( GTF.iterator( IOTools.openFile( temp_rejected ) ) ):
            if gtf[0].gene_id in to_rescue:
                for exon in IOTools.flatten( gtf ):
                    temp_filtered_out.write( str( exon ) + "\n" )
                    rescued_gtf_out.write( str( exon ) + "\n" )
            else:
                for exon in IOTools.flatten( gtf ):
                    def_rejected_out.write( str( exon ) + "\n" )


        os.unlink( temp_rejected )
        def_rejected_out.close()
        rescued_gtf_out.close()

    else: 
        def_rejected = temp_rejected
    
    temp_filtered_out.close()

    # sort the outfiles by gene
    statement = ( "cat %(temp_filtered)s |"
                  " python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=sort --sort-order=gene+transcript"
                  "  --log %(filtered_gtf)s.log |"
                  " gzip"
                  " > %(filtered_gtf)s;"
                  " cat %(def_rejected)s |"
                  " python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=sort --sort-order=gene+transcript"
                  "  --log %(filtered_gtf)s.log |"
                  " gzip"
                  " > %(rejected_gtf)s;" )
    P.run()

    os.unlink( temp_filtered )
    os.unlink( def_rejected )

def collateChunks( gtf_chunks, distance ):
    """
    Receives an iterator_sorted_chunks object and a specified distance. If the 
    adjacent chunks are less than specified distance apart, they are returned
    as a list in the resulting generator. 
    """
    # take first interval in generator 
    last = gtf_chunks.next()
    to_join = []
    to_join.extend( last )
    
    for gtf in gtf_chunks:
        # check distance between adjacent intervals in generator
        # NB. last will always be the last element in to_join, regardless of genomic position
        dist = gtf[0].start - max([ x.end for x in to_join ])
        # if adjacent intervals are on same contig & strand check order
        if gtf[0].contig == last[0].contig and gtf[0].strand == last[0].strand:
            assert gtf[0].start >= last[0].start, "gene features are not sorted!"
        # if adjacent intervals are not suitable for merging, yield list
        if gtf[0].contig != last[0].contig or gtf[0].strand != last[0].strand or dist > distance:
            yield to_join
            to_join = []
        # set current interval to last and add to list
        last = gtf
        to_join.extend( last )

    yield to_join
    raise StopIteration

def calcIntergenicDist( gtf_chunks, ignore_overlap = False ):
    """
    Receives an iterator_sorted_chunks object and calulates the distance between
    adjacent chunks... returns as list. When adjacent gene models overlap the 
    one with the last end co-ordinate is retained if ignore_overlap, 
    then overlapping intervals are ignored.
    """
    # take first interval in generator
    last = gtf_chunks.next()     
    distances = {}
    contains = collections.defaultdict( list )
    
    for gtf in gtf_chunks:
        # if adjacent intervals are not on the same contig then continue
        if gtf[0].contig != last[0].contig or gtf[0].strand != last[0].strand:
            last = gtf
            continue
        # sanity check
        assert gtf[0].start >= last[0].start, "gene features are not sorted!"
        # check distance between adjacent intervals in generator
        dist = gtf[0].start - last[-1].end
        
        # if distance is negative, then check which has largest end co-ordinate
        # if gtf is contained within last, write to a dictionary
        if dist < 0:
            if gtf[-1].end < last[-1].end:
                contains[ last[0].gene_id ].append( gtf[0].gene_id )
                # print( last[0].gene_id, last[0].start, last[-1].end, 
                #        gtf[0].gene_id, gtf[0].start, gtf[-1].end )
                if ignore_overlap:
                    continue
            else: 
                last = gtf
                if ignore_overlap:
                    continue
        
        distances[ gtf[0].gene_id ] = dist
                                                                     
        last = gtf

    return  distances, contains

def splitGTFOnStrand( gtf_file ):
    """
    Splits transcripts on strand, returning two tempfile names
    """
    gtf_plus = P.getTempFile( "." )
    gtf_minus = P.getTempFile( "." )
    for gtf in GTF.transcript_iterator( GTF.iterator( IOTools.openFile( gtf_file ) ) ):
        if gtf[0].strand == "+":
            for exon in gtf: 
                gtf_plus.write( str( exon ) + "\n" )
        elif gtf[0].strand == "-":
            for exon in gtf:
                gtf_minus.write( str( exon ) + "\n" )
        else:
            raise ValueError( "Unrecognised strand"
                              " orientation in gtffile: %s" % gtf[0].strand )
    gtf_plus.close()
    gtf_minus.close()
    
    return gtf_plus.name, gtf_minus.name

def returnSourceIntervals( in_gtf, source, interval = "gene" ):
    """
    Iterate through a gtf file and return a list of gene_ids/transcript_ids for
    which source field ($2) matches the given string. 
    """
    assert interval in [ "gene", "transcript" ], "interval must either be 'gene' or 'transcript'"

    # iterate through gtf and return list of ids with specified source. 
    id_set = set()
    for gtf in GTF.iterator( IOTools.openFile( in_gtf ) ):
        if gtf.source == source:
            if interval == "gene":
                id_set.add( gtf.gene_id )
            else: 
                id_set.add( gtf.transcript_id )
        else:
            continue

    return list( id_set )
        

def getCommonNames( infile, outfile ):
    """
    Iterate over all flat gene models in a gtf file, retreive the gene name, 
    the gene id and the old id... write them to an outfile to be loaded into 
    csvdb
    If no common name exists, then write gene_id to gene_name field.
    """
    inf = IOTools.openFile( infile )
    outf = IOTools.openFile( outfile, "w" )
    outf.write( "gene_id\told_gene_id\tgene_name\n" )

    # there are gene_ids with duplicate gene names... fix this. 
    gene_names = {}
    n_duplicates = 0

    for gtf in GTF.flat_gene_iterator( GTF.iterator( inf ) ):
        gene_id = gtf[0].gene_id
        old_id = gtf[0].oId
        # for lncRNAs, keep the gene_id as gene_name
        if re.match( "LNC", gene_id ):
            gene_name = gene_id
        else:
            try: 
                gene_name = gtf[0].gene_name
                if gene_name in gene_names.keys():
                    gene_names[gene_name] += 1
                    gene_name = gene_name + "_" + str( gene_names[gene_name] )
                    n_duplicates += 1
                else:
                    gene_names[gene_name] = 0

            except AttributeError:
                gene_name = gene_id
        outf.write( "\t".join( [ gene_id, old_id, gene_name ] ) + "\n" )

    outf.close()
    E.info( "There are %i duplicate gene names in the"
            " protein coding gene set" % n_duplicates )



def gtfToBed12(infile, outfile):
    """
    Convert a gtf file of transcripts to bed12 format 
    where name additional name field is transcript id
    """
    outfile = IOTools.openFile(outfile, "w")

    for all_exons in GTF.transcript_iterator(
            GTF.iterator(IOTools.openFile(infile, "r"))):
        chrom = all_exons[0].contig
        # GTF.iterator returns start co-ordinates as zero-based
        start = str(all_exons[0].start)
        end = str(all_exons[len(all_exons) - 1].end)
        name_1 = all_exons[0].gene_id 
        name_2 = all_exons[0].transcript_id
        score = "0"
        strand = all_exons[0].strand
        thickStart = start
        thickEnd = end
        colourRGB = "0"
        blockCount = str(len(all_exons))

        sizes = []
        starts = []
        for exon in all_exons:
            blockSize = str(exon.end - (exon.start))
            sizes.append(blockSize)
            # start + blockStart should return a zero-based co-ordinate for the
            # exon start
            blockStart = str((exon.start) - int(start))
            starts.append(str(blockStart))
        blockSizes = ','.join(sizes)
        blockStarts = ','.join(starts)

        outfile.write(chrom + "\t"
                      + start + "\t"
                      + end + "\t"
                      + name_1 + "\t"
                      + name_2 + "\t"
                      + score + "\t"
                      + strand + "\t"
                      + thickStart + "\t"
                      + thickEnd + "\t"
                      + colourRGB + "\t"
                      + blockCount + "\t"
                      + blockSizes + "\t"
                      + blockStarts + "\n")

    outfile.close()


#################################################################################
#################################################################################
#################################################################################
## section: manipulate bed files
#################################################################################
def splitBedOnStrand( bed_file ):
    """
    Splits bed file on strand, returning two tempfilenames
    """
    bed_plus = P.getTempFile( "." )
    bed_minus = P.getTempFile( "." )
    for bed in Bed.iterator( IOTools.openFile( bed_file ) ):
        if bed.strand == "+":
            bed_plus.write( str( bed ) + "\n" )
        elif bed.strand == "-":
            bed_minus.write( str( bed ) + "\n" )
        else:
            raise ValueError( "Unrecognised strand"
                              " orientation in bedfile: %s" % bed.strand )
    bed_plus.close()
    bed_minus.close()

    return bed_plus.name, bed_minus.name


def resetTSS( in_df, 
              fantom_bed, 
              tss_bed, 
              outfile, 
              distance = 500,
              reporter= "gene", 
              stringency = "permissive" ): 
    """
    Receives i) dataframe of distance between lncRNA and nearest CAGE peak, 
    ii) fantom bedfile containing details of cage peaks, iii) bedfile of lncRNA
    TSS. For any pairs less than 'distance' apart, the cage peak summit replaces
    the lnc TSS in outfile.
    """
    ## Select only peaks within specified distance
    # for clarity, specify distance and cage id column names
    dist_col = "_".join( ( reporter, stringency, "distance" ) )
    cage_id = "_".join( ( reporter, stringency, "cage_id" ) )
    # create df of features <500bp from CAGE peak
    # Pandas prints a random error when faced with single element boolean queries!
    df = in_df.query( "%s < %s" % ( dist_col, str( distance ) )  )
    # # retain only columns of interest
    df = df[ [ cage_id, "lnc_id" ] ]

    cage_ids = df[cage_id].tolist()
    ## create dictionary containing lnc_id and new TSS coordinates
    lnc_dict = {}
    x = 0
    for line in IOTools.openFile( fantom_bed ):
        line = line.split()
        if line[3] in cage_ids:
            x += 1
            # get all entries in df for cage_id (there
            match = df.ix[ df[ cage_id ] == line[3] ]
            # reset index... necessary for iterating through rows
            match = match.reset_index()
            if reporter == "gene":
                assert len( match.index.values ) == 1, \
                    "CAGE peak matches multiple lncRNAs!!"
            # in the event that one cage peak matches more than one tss (likely 
            # for transcripts), add each lnc_id to dictionary separately.
            for index, row in match.iterrows():
                lnc_dict[ match.iloc[index]["lnc_id"] ] = ( line[6], line[7] )

    assert len( df.index.values ) == len( [ x for x in lnc_dict.values() if x ] ), \
        "Mismatch between # of matches in dataframe and number of entries in dictionary"

    ## iterate through TSS bedfile and replace co-ordinates for entries in dictionary
    # write new co-ordinates to outfile
    outf = IOTools.openFile( outfile, "w" )
    x = 0
    for bed in Bed.iterator( IOTools.openFile( tss_bed ) ):
        assert re.match( "LNC", bed.fields[0] ), "The entries in TSS bed aren't lncRNAs"
        if bed.fields[0] in lnc_dict.keys():
            x += 1
            bed.start = lnc_dict[ bed.fields[0] ][0]
            bed.end = lnc_dict[ bed.fields[0] ][1]
            bed.fields[1] = 1
            outf.write( str( bed ) + "\n" )
        else:
            outf.write( str( bed ) + "\n" )
    outf.close()
    assert len( df.index.values ) == x, "Entries in table that are not in lnc_tss file"


def collateK4CoverageBed( infiles, outfile, params ):
    """
    Pull the coverage across intervals from me1 and me3 coverageBed files,
    output the result as a table
    """
    # retrieve params
    pseudocount, min_cov = [ int(x) for x in params ]

    # assign infiles
    for file in infiles:
        if re.search( "_K4me3_", file ):
            me3_bed = file
        elif re.search( "_K4me1_", file ):
            me1_bed = file
        else:
            raise IOError( "Unrecognised condition in infile: %s" % file )

    assert os.path.exists( me3_bed ), "Missing H3K4me3 file"
    assert os.path.exists( me1_bed ), "Missing H3K4me1 file"
    
    # collate coverage information
    coverage_dict = collections.defaultdict( list )
    for bed in Bed.iterator( IOTools.openFile( me1_bed ) ):
        coverage_dict[ bed.fields[0] ].append( bed.fields[3] )
    for bed in Bed.iterator( IOTools.openFile( me3_bed ) ):
        coverage_dict[ bed.fields[0] ].append( bed.fields[3] )

    # write outfile
    outf = IOTools.openFile( outfile, "w" )
    outf.write( "lnc_id\tme1_cov\tme3_cov\tratio\tpass\n" )
    for key, values in coverage_dict.iteritems():
        # calculate me1/me3 ratio
        me1 = float( values[0] )
        me3 = float( values[1] )
        ratio = ( me1 + pseudocount )/( me3 + pseudocount )
        suff_cov = None
        if max( [ me1, me3 ] ) >= min_cov:
            suff_cov = 1
        else:
            suff_cov = 0
        outf.write( "\t".join( map( str, [ key, 
                                           me1, 
                                           me3, 
                                           ratio, 
                                           suff_cov ] ) ) + "\n" )
    outf.close()


def findOverlap( in_query, in_target, out_overlap ):
    """
    Run bedtools interesect and report the number of intervals in query that do
    and do not intersect target. Assumes strandedness.
    """
    # specify outfiles
    out_no_overlap = P.snip( out_overlap, ".overlap.bed" ) + ".no_overlap.bed" 

    # get the number of intervals overlapping
    statement = ( "bedtools intersect"
                  " -a %(in_query)s"
                  " -b %(in_target)s"
                  " -wa"
                  " -u"
#                  " -s"
                  " > %(out_overlap)s" )
    P.run()

    statement = ( "bedtools intersect"
                  " -a %(in_query)s"
                  " -b %(in_target)s"
                  " -wa"
                  " -v"
#                  " -s"
                  " > %(out_no_overlap)s" )
    P.run()


def findIntersect( infile_1, infile_2, outfile ):
    """
    Output a bedfile of the intersecting gene_ids in two bedfiles
    """
    first = [ line.split()[3] for line in IOTools.openFile( infile_1 ).readlines() ]
    second =  [ line.split()[3] for line in IOTools.openFile( infile_2 ).readlines() ]
    
#    assert len(first) == len(set(first)), "Duplicates in %s" % infile_1
#    assert len(second) == len(set(second)), "Duplicates in %s" % infile_2
    
    intersection = list(set(first) and set(second))

    outf = IOTools.openFile( outfile, "w" )
    for line in IOTools.openFile( infile_1 ).readlines():
        if line.split()[3] in intersection:
            outf.write( line )
    outf.close()




#################################################################################
#################################################################################
#################################################################################
## section: Coverage Counter
#################################################################################

def assessReadOrientation( read, strand, library_type, include_mate = True ):
    """
    This function takes an instance of pysam read object (returned by 
    samfile.fetch()) together with strand of a gene feature. According to library
    type specified (fr-firststrand, fr-secondstrand, unstranded), it will return 
    true or false depending on whether the read maps to same strand as feature.
    Additionally if include_mate == False, will return False if mate is second in
    pair. 
    """

    if library_type not in [ "fr-firststrand", "fr-secondstrand", "unstranded" ]:
        raise ValueError( "Unrecognised library type: %s" % library_type )
    if strand not in [ "+", "-", "." ]:
        raise ValueError( "Feature strand must be '+', '-', of '.' in gtf " )

    if read.is_read2 and not include_mate:
        condition == False
    elif library_type == "fr-firststrand": 
        if strand == "+":
            condition = ( read.is_read1 and read.is_reverse ) or ( read.is_read2 and not read.is_reverse )
        elif strand == "-":
            condition = ( read.is_read2 and read.is_reverse ) or ( read.is_read1 and not read.is_reverse )
        else:
            raise ValueError( "Feature strand must be either '+' or '-' for standed libraries" )
    elif library_type == "fr-secondstrand":
        if strand == "+":
            condition = ( read.is_read1 and not read.is_reverse ) or ( read.is_read2 and read.is_reverse )
        elif strand == "-":
            condition = ( read.is_read2 and not read.is_reverse ) or ( read.is_read1 and read.is_reverse )
        else:
            raise ValueError( "Feature strand must be either '+' or '-' for standed libraries" )
    else:
        condition = True
            
    return condition
       

def coverageCounter( samfile, 
                     contig, 
                     start, 
                     end, 
                     strand, 
                     library_type = "fr-firststrand", 
                     include_mate = True ):

    """
    Takes bamfile interval co-ordinates (contig, start, end, strand) and returns 
    i) pptn of bases in interval that are covered by one or more reads ,
    ii) number of reads mapping to the interval,
    iii) a list of number of reads covering each base in interval. 
    Accounts for stranded data (library_type fr-firststrand, fr-secondstrand,
    unstranded). Allows coverage to be counted using both reads in pair
    (include_mate == True), or only first read in pair (include_mate == False).
    N.B. for spliced reads, coverage is only counted for the portion read mapping
    to interval. However, spliced read are still counted as mapping to interval
    (nreads += 1), meaning they may inflate read counts if the returned nreads
    values are summed for all exons within a gene model.
    """

    samfile = pysam.Samfile( samfile, "rb" )
    interval_width = end - start
    nreads = 0
    
    # create a counter for the interval
    counts = []
    for n in range( 0, interval_width ):
        counts.append( 0 )
        
    for read in samfile.fetch( contig, start, end ):
        # only consider reads on the correct strand
        condition = assessReadOrientation( read, strand, library_type, include_mate )
#        if stranded:
#            if strand ==  "+":
#                condition = ( read.is_read1 and read.is_reverse ) or ( read.is_read2 and not read.is_reverse )
#            else: 
#                condition = ( read.is_read2 and read.is_reverse ) or ( read.is_read1 and not read.is_reverse )
#        else: 
#            condition = True
        
        if condition:
            # check to see if reads are spliced
            if "N" in read.cigarstring:
                # iterate through the cigar string
                # create a list of tuples with start end of read fragments
                fragments = []
                fragment_start = read.pos
                for operation in read.cigar:
                    if operation[0] == 0:
                        fragments.append( (fragment_start, 
                                           fragment_start + operation[1] ) )
                        fragment_start += operation[1]
                    elif operation[0] == 3:
                        fragment_start += operation[1]
                    else: continue
                
                # iterate through list of fragments
                # discard those that fall outside interval
                for fragment in fragments:
                    if fragment[0] >= end: 
                        break
                    else:
                        rstart = max( 0, fragment[0] - start )
                        rend = min( interval_width, fragment[1] - start )
                        for i in range( rstart, rend ):
                            counts[i] += 1
                nreads += 1
            
            # iterate count the entire read
            else:
                nreads += 1
                rstart = max( 0, read.pos - start ) 
                rend = min( interval_width, read.aend - start )
                for i in range( rstart, rend ):
                    counts[i] += 1
        else: continue    
    
    coverage = 0
    for i in counts:
        if i != 0:
            coverage += 1
    coverage = float(coverage)/float( len( counts ) )
    
    return coverage, nreads, counts


def findExpressedGenesFromCoverage( gtf_file, 
                                    bamfiles, 
                                    coverage_threshold, 
                                    replicate_threshold,
                                    sum_stat = "mean" ):

    """
    Receives a gtf file and list of bamfiles, returns a pandas Series with mean 
    coverage for each feature in gtf file. 
    Also returns a pandas Series specifying whether each feature is 'expressed', 
    based on a specified coverage threshold and minumum number of replicates in 
    which this threshold must be reached. 
    N.B. coverage just refers to pptn of bases in feature that are covered by 
    one or more reads, it doesn't account for depth. 
    """

    # create nested dict for storing per-sample coverage
    coverage_dict = collections.defaultdict( dict )
    # create dict for storing number of samples above threshold
    threshold_dict = {}

    # loop through bamfiles
    for bamfile in bamfiles:
        gene_ids = []
        sample_id = P.snip( os.path.basename( bamfile ), ".bam" )
        print( "\t\tProcessing bamfile for %s" % sample_id )

        # obtain coverage for each feature in gtf
        n = 0
        for gtfs in GTF.flat_gene_iterator( GTF.iterator( IOTools.openFile( gtf_file ) ) ):
            n += 1
            gene_counts = []
            gene_id = gtfs[0].gene_id
            assert gene_id not in gene_ids, "Gene_id appears in gtf multiple times! "
            gene_ids.append( gene_id )
            contig = gtfs[0].contig
            strand = gtfs[0].strand
            # iterate through exons
            for gtf in gtfs:
                start = gtf.start
                end = gtf.end 

                # use coverage counter to get stranded coverage for exon
                coverage, nreads, counts = coverageCounter( bamfile, 
                                                            contig, 
                                                            start, 
                                                            end, 
                                                            strand )
                # extend transcript coverage count by exon
                gene_counts.extend( counts )
            
            # calculate proportion of merged transcript that is covered
            coverage = 0
            for i in gene_counts:
                if i != 0:
                    coverage += 1
            coverage = float(coverage)/float( len( gene_counts ) )
            # add coverage value to 
            coverage_dict[ gene_id ][sample_id] = coverage

            # if merged transcript coverage is above expression threshold
            if coverage >= coverage_threshold:
                if gene_id not in threshold_dict:
                    threshold_dict[ gene_id ] = 1
                else: threshold_dict[ gene_id ] += 1
            else:
                if gene_id not in threshold_dict:
                    threshold_dict[ gene_id ] = 0

    # create a dictionary of mean coverage
    mean_cov_dict = {}
    for gene in gene_ids:
        if sum_stat == "median":
            mean_cov = np.median( coverage_dict[ gene ].values() )
        elif sum_stat == "mean":
            mean_cov = np.mean( coverage_dict[ gene ].values() )
        else:
            raise ValueError( "Unrecognised summary statistic" )
        mean_cov_dict[ gene ] = mean_cov 
        # convert the threshold dict to boolean
        if threshold_dict[ gene ] >= replicate_threshold:
            threshold_dict[ gene ] = 1
        else:
            threshold_dict[ gene ] = 0

    series_coverage = pd.Series( mean_cov_dict, index = gene_ids )
    series_expressed = pd.Series( threshold_dict, index = gene_ids )

    return series_coverage, series_expressed

#################################################################################
#################################################################################
#################################################################################
## section: lncRNA filters
#################################################################################

def filterOnIntergenicCoverage( lncRNA_gtf, 
                                tts_bed,
                                bamfiles,
                                window, 
                                max_coverage, 
                                max_occurrence, 
                                filtered_gtf, 
                                rejected_gtf, 
                                rejected_tsv, 
                                outdir ):

    """
    Searches for high intergenic coverage between protein coding genes and 
    downstream lncRNAs. lncRNAs within a specified distance from protein coding
    tts are checked for high intergenic coverage using each bamfile supplied.
    Any lncRNA with > than max_coverage in > max_occurrence samples will be 
    removed from the final set. 
    Coverage is calculated as the proportion of bases that have one or more 
    reads mapping to them. Regions encompassed by spliced reads are not counted 
    as covered. 
    """
    
    # index regions downstream of protein coding gene tts
    tts_index = IndexedGenome.IndexedGenome()

    for bed in Bed.iterator( IOTools.openFile( tts_bed ) ):
        if bed.fields[2] == "+":
            tts_index.add( bed.contig, 
                           bed.end, 
                           bed.end + window, 
                           str( bed ).split() )
        else:
            tts_index.add( bed.contig, 
                           bed.start - window, 
                           bed.start, 
                           str( bed ).split() )

    # for all lncRNA within <window> dist downstream of a protein coding tts
    # output a bed entry of the shortest interval between lncRNA start and tts
    test_bed = tempfile.NamedTemporaryFile( dir=outdir, 
                                            delete=False, 
                                            suffix="_test_intervals" )

    inf = IOTools.openFile( lncRNA_gtf )
    for gtf in GTF.merged_gene_iterator( GTF.iterator( inf ) ):
        if tts_index.contains( gtf.contig, gtf.start, gtf.end ):
            interval_dist = 0
            out_fields = []
            for interval in tts_index.get( gtf.contig, gtf.start, gtf.end ):
                if interval[2][5] == gtf.strand == "+":
                    if interval_dist == 0:
                        interval_dist = gtf.start - interval[0]
                        out_fields = [ gtf.contig, 
                                       interval[0], 
                                       gtf.start + 1, 
                                       gtf.gene_id + "__" + interval[2][3], 
                                       "0", 
                                       gtf.strand ] 
                    elif interval_dist > gtf.start - interval[0]:
                        interval_dist = gtf.start - interval[0]
                        out_fields = [ gtf.contig, 
                                       interval[0],
                                       gtf.start + 1, 
                                       gtf.gene_id + "__" + interval[2][3], 
                                       "0", 
                                       gtf.strand ]
                    else: continue
                elif interval[2][5] == gtf.strand == "-":
                    if interval_dist == 0:
                        interval_dist = interval[1] - gtf.end
                        out_fields = [ gtf.contig, 
                                       gtf.end - 1, 
                                       interval[1],
                                       gtf.gene_id + "__" + interval[2][3], 
                                       "0",
                                       gtf.strand ]
                    elif interval_dist > interval[1] - gtf.end:
                        interval_dist = interval[1] - gtf.end
                        out_fields = [ gtf.contig, 
                                       gtf.end -1, 
                                       interval[1],
                                       gtf.gene_id + "__" + interval[2][3], 
                                       "0",
                                       gtf.strand ]
                    else: continue
                else: continue

            if out_fields:
                if out_fields[1] > out_fields[2]:
                    lnc_id, pc_id = out_fields[3].split("__")
                    E.warn( "Warning, sense overlap between %s and %s"
                           " - interval skipped" % (lnc_id, pc_id) )
                else:
                    output_string = "\t".join( [ str( x ) for x in out_fields ] )
                    test_bed.write( output_string + "\n" )
                    
        else: continue
                
    test_bed.close()

    # for each bamfile, output a bedfile containing coverage in $5
    # output a histogram of per-base read coverage
    for bamfile in bamfiles:
        out_bed = P.snip( os.path.basename( bamfile ), ".bam" ) + "_coverage.bed.gz" 
        out_bed = IOTools.openFile( os.path.join( outdir, out_bed ), "w" )
        out_hist = P.snip( os.path.basename( bamfile ), ".bam" ) + "_coverage.map.gz"
        out_hist = IOTools.openFile( os.path.join( outdir, out_hist ), "w" )

        # loop through the test_bed file containing intervals
        # is it necessary to close and re-open the test intervals?
        for bed in Bed.iterator( IOTools.openFile( test_bed.name ) ):
            strand = bed.fields[2]
            coverage, nreads, counts = coverageCounter( bamfile, 
                                                        bed.contig, 
                                                        bed.start, 
                                                        bed.end, 
                                                        strand )
            out_list = [bed.contig, 
                        bed.start, 
                        bed.end, 
                        bed.fields[0], 
                        coverage, 
                        bed.fields[2] ]
            string = "\t".join( [ str(x) for x in out_list ] )           
            out_bed.write( string + "\n" )

            lnc_id, pc_id = bed.fields[0].split( "__" )
            out_hist.write( lnc_id + "\t" + ",".join( [ str(x) for x in counts ] ) )

        out_bed.close()
        out_hist.close()

    # for each coverage bedfile, add interval with > max_coverage to dictionary
    high_coverage = {}
    for coverage_file in os.listdir( outdir ):
        if not coverage_file.endswith( "coverage.bed" ): 
            continue
        else:
            coverage_file = os.path.join( outdir, coverage_file )
            sample = P.snip( os.path.basename( coverage_file ), "_coverage.bed" )
            for bed in Bed.iterator( IOTools.openFile( coverage_file ) ):
                lncRNA_id, protein_coding_id = bed.fields[0].split( "__" )
                if float( bed.fields[1] ) >= max_coverage:
                    if lncRNA_id not in high_coverage:
                        high_coverage[ lncRNA_id ] = [ 1, 
                                                       [ protein_coding_id ], 
                                                       [ sample ] ]
                    else:
                        high_coverage[ lncRNA_id ][0] += 1
                        high_coverage[ lncRNA_id ][1].append( protein_coding_id )
                        high_coverage[ lncRNA_id ][2].append( sample )
                else: continue
    E.info( "%i lncRNA loci have high upstream read coverage in one or more"
            " samples" % len( [ x for x in high_coverage.iterkeys() ] ) )

    # create a list of lncRNA gene_ids to be rejected 
    to_reject = []
    rejected_tsv = IOTools.openFile( rejected_tsv, "w" )
    rejected_tsv.write( "lncrna_id\tupstr_gene_id\t#_samples\tsample_ids\n" )

    for key, value in high_coverage.iteritems():
        if value[0] >= int( max_occurrence ):
            to_reject.append( key )
            out_list = [ key, 
                         ",".join( set( value[1] ) ), 
                         value[0], 
                         ",".join( value[2] ) ]  
            string = "\t".join( [ str(x) for x in out_list ] )
            rejected_tsv.write( string + "\n" )
        else: continue

    rejected_tsv.close()

    # remove lncRNAs with > max_coverage in > max_occurrance samples
    rejected_gtf = IOTools.openFile( rejected_gtf, "w" )
    filtered_gtf = IOTools.openFile( filtered_gtf, "w" )

    c_dict = { "rejected": 0, "filtered": 0 }

    inf = IOTools.openFile( lncRNA_gtf )
    for gtf in GTF.flat_gene_iterator( GTF.iterator( inf ) ):
        if gtf[0].gene_id in to_reject:
            outf = rejected_gtf
            c_dict[ "rejected" ] += 1
        else:
            outf = filtered_gtf
            c_dict[ "filtered" ] += 1
        for exon in IOTools.flatten( gtf ):
            outf.write( str( exon ) + "\n" )

    rejected_gtf.close()
    filtered_gtf.close()

    E.info( "%i lncRNA loci were removed due to high read coverage between"
            " lncRNA tss and tts of upstream protein coding gene\n %i lncRNA"
            " in filtered gtf" % ( c_dict["rejected"], c_dict["filtered"] ) )

#    os.unlink( test_bed )
#    return outdir


def filterOnSharedSplicedReads( params ):

    """
    Search for spliced reads in lncRNA exons that are shared with a set of 
    reference exons. For all spliced reads intersecting lncRNA exon, if start/end
    of spliced read falls within exon boundaries, check to see if the opposite
    end falls within an indexed set of protein coding gene exons. 
    Doesn't account for instances where lncRNA and protein coding exons overlap.
    """

    # function is set up to run with P.submit, meaning that all parameters have
    # to be passed as strings. 
    lncRNA_gtf, reference_gtf, bamfiles, max_occurrence, filtered_gtf, rejected_gtf, rejected_tsv = params
    bamfiles = bamfiles.split( "__" )
    max_occurrence = int( max_occurrence )
    
    # create an index of exons in reference set
    ref_index = IndexedGenome.IndexedGenome()
    for gtf in GTF.iterator( IOTools.openFile( reference_gtf ) ):
        ref_index.add( gtf.contig, 
                       gtf.start, 
                       gtf.end, 
                       ( gtf.gene_id, 
                         gtf.strand ) )
    
    # create dictionary of lncRNAs that share spliced reads with protein coding 
    # genes.
    ## Question - what happens to unmapped reads wrt read.is_reverse?
    overlap = {}
    for bamfile in bamfiles:
        sample = P.snip( os.path.basename( bamfile ), ".bam" )
        samfile = pysam.Samfile( bamfile, "rb" )
        for gtf in GTF.iterator( IOTools.openFile( lncRNA_gtf ) ):
            for read in samfile.fetch( gtf.contig, gtf.start, gtf.end ):
                if ( gtf.strand == "+" and ( ( read.is_read1 and read.is_reverse ) or ( read.is_read2 and not read.is_reverse ) ) ) or ( gtf.strand == "-" and ( ( read.is_read2 and read.is_reverse ) or ( read.is_read1 and not read.is_reverse ) ) ):
                    if "N" in read.cigarstring:
                        if gtf.start < read.pos < gtf.end and ref_index.contains( gtf.contig, read.aend -1, read.aend ):
                            for interval in ref_index.get( gtf.contig, 
                                                           read.aend - 1, 
                                                           read.aend + 1 ):
                                if interval[2][1] == gtf.strand:
                                    if gtf.gene_id not in overlap:
                                        # first list is protein_coding gene_id, 
                                        # second list is sample
                                        overlap[ gtf.gene_id ] = [ 1, 
                                                                   [ interval[2][0] ], 
                                                                   [ sample ] ]
                                    else:
                                        overlap[ gtf.gene_id ][0] += 1
                                        overlap[ gtf.gene_id ][1].append( interval[2][0] )
                                        overlap[ gtf.gene_id ][2].append( sample )
                                    break
                        elif gtf.start < read.aend < gtf.end and ref_index.contains( gtf.contig, read.pos, read.pos + 1 ):
                            for interval in ref_index.get( gtf.contig, read.pos, read.pos + 1 ):
                                if interval[2][1] == gtf.strand:
                                    if gtf.gene_id not in overlap:
                                        overlap[ gtf.gene_id ] = [ 1, 
                                                                   [ interval[2][0] ], 
                                                                   [ sample ] ]
                                    else:
                                        overlap[ gtf.gene_id ][0] += 1
                                        overlap[ gtf.gene_id ][1].append( interval[2][0] )
                                        overlap[ gtf.gene_id ][2].append( sample )
                                    break

    # write gene_ids of lncRNA that share spliced reads to a tsv file
    to_reject = []

    lncRNA_rejected_tsv = IOTools.openFile( rejected_tsv, "w" )
    lncRNA_rejected_tsv.write( "lncRNA_id\t#_shared_spliced_reads\t"
                               "#_shared_genes\t#_samples\t"
                               "shared_gene_ids\tsample_ids\n" )    
    for key, value in overlap.iteritems():
        gene_id = key
        sj_num = value[0]
        sj_num_genes = len( set( value[1] ) )
        sj_num_samples = len( set( value[2] ) )
        sj_genes = ",".join( x for x in set( value[1] ) )
        sj_samples = ",".join( x for x in set( value[2] ) )
        if sj_num_samples >= max_occurrence:
            to_reject.append( gene_id )
        out_list = [ gene_id, sj_num, sj_num_genes, sj_num_samples, sj_genes, sj_samples ]
        out_string = "\t".join( str(x) for x in out_list )
        lncRNA_rejected_tsv.write( out_string + "\n" )

    lncRNA_rejected_tsv.close()

    # write filtered and rejected lncRNAs to different gtf files
    lncRNA_rejected = IOTools.openFile( rejected_gtf, "w" )
    lncRNA_filtered = IOTools.openFile( filtered_gtf, "w" )


    for gtf in GTF.flat_gene_iterator( GTF.iterator( IOTools.openFile( lncRNA_gtf ) ) ):
        if gtf[0].gene_id in to_reject:
            outf = lncRNA_rejected
        else:
            outf = lncRNA_filtered
        for exon in gtf:
            outf.write( str( exon ) + "\n" )

    lncRNA_rejected.close()
    lncRNA_filtered.close()



def write_to_temp(tempfile, interval_list, transcript, check_strand=True):
    if check_strand:
        for interval in interval_list:
            if interval[2][6] == transcript[0].strand:
                tempfile.write(transcript[0].gene_id + "\t"
                               + str(interval[0]) + "\t"
                               + str(interval[1]) + "\t"
                               + str(transcript[0]) + "\t"
                               + "\t".join(interval[2]) + "\n")
    else:
        for interval in interval_list:
            tempfile.write(transcript[0].gene_id + "\t"
                           + str(interval[0]) + "\t"
                           + str(interval[1]) + "\t"
                           + str(transcript[0]) + "\t"
                           + "\t".join(interval[2]) + "\n")


def ClassifyLncRNAGenes_detail(lncRNA_gtf,
                               reference_gtf,
                               outfile,
                               upstr_dist=5,
                               dstr_dist=5,
                               wdir=".", 
                               prioritize_flanks=True,
                               distinguish_antisense_introns=False):
    """
    This re-write of reClassifyLncRNAGenes(). However, it prioritizes a lncRNA
    as antisense_upstream/downstream rather than antisense. 
    It also separately classifies antisense, and antisense intronic loci.
    To be in any category requires a bp overlap of  1 or more... 
    """

    # index exons in the reference gene-set
    ref_index = IndexedGenome.IndexedGenome()
    for exon in GTF.iterator(IOTools.openFile(reference_gtf)):
        ref_index.add(exon.contig, exon.start, exon.end, str(exon).split())

    # create index for all other intervals to be classified
    intron = IndexedGenome.IndexedGenome()
    plus_up = IndexedGenome.IndexedGenome()
    plus_down = IndexedGenome.IndexedGenome()
    minus_up = IndexedGenome.IndexedGenome()
    minus_down = IndexedGenome.IndexedGenome()

    # iterate over reference transcripts and create intervals in memory
    ref_file = IOTools.openFile(reference_gtf)
    for transcript in GTF.transcript_iterator(GTF.iterator(ref_file)):
        start = transcript[0].end
        for i in range(1, len(transcript)):
            intron.add(transcript[i].contig,
                       start,
                       transcript[i].start,
                       str(transcript[i]).split())
            start = transcript[i].end

        # create up and downstream intervals on plus strand
        if transcript[0].strand == "+":
            plus_up.add(transcript[0].contig,
                        transcript[0].start - (upstr_dist * 1000),
                        transcript[0].start,
                        str(transcript[0]).split())
            plus_down.add(transcript[0].contig,
                          transcript[len(transcript) - 1].end,
                          transcript[
                              len(transcript) - 1].end + (dstr_dist * 1000),
                          str(transcript[len(transcript) - 1]).split())

        # create up and downstream intervals on minus strand
        elif transcript[0].strand == "-":
            minus_up.add(transcript[0].contig,
                         transcript[len(transcript) - 1].end,
                         transcript[len(transcript) - 1].end +
                         (upstr_dist * 1000),
                         str(transcript[len(transcript) - 1]).split())
            minus_down.add(transcript[0].contig,
                           transcript[0].start - (dstr_dist * 1000),
                           transcript[0].start,
                           str(transcript[0]).split())
        else:
            E.warn("WARNING: no strand specified for %s" %
                   transcript[0].transcript_id)

    # create single representative transcript for each lncRNA gene_id
    merged_lncRNA_gtf = P.getTempFilename(wdir)
    to_cluster = False
    statement = ("zcat %(lncRNA_gtf)s |"
                 " python %(scriptsdir)s/gtf2gtf.py"
                 "  --method=sort --sort-order=gene"
                 "  --log=%(outfile)s.log |"
                 " python %(scriptsdir)s/gtf2gtf.py"
                 "  --method=merge-exons"
                 "  --log=%(outfile)s.log"
                 " > %(merged_lncRNA_gtf)s")
    P.run()

    # create a temp directory containing the indexed intervals used to classify
    # the lncRNA transcripts created (for debugging purposes)
    # create a temporary count of # of gene_models in each category
    tempdir = P.getTempDir(wdir)
    E.info("intersecting intervals are being written to %s"
           % os.path.abspath(tempdir))
    temp_file_names = ["sense",
                       "sense_intronic",
                       "sense_overlap",
                       "antisense",
                       "sense_downstream",
                       "sense_upstream",
                       "antisense_downstream",
                       "antisense_upstream",
                       "intergenic"]
    temp_files = {}
    temp_count = {}
    for handle in temp_file_names:
        temp_count[handle] = 0
        temp_files[handle] = IOTools.openFile(os.path.join(tempdir,
                                                           handle), "w")

    # iterate through the representative (i.e. merged) lncRNA transcripts
    # each lncRNA transcript is classified only once.
    # In situations where a lncRNA fits > 1 classification, priority is:
    # (sense > antisense)
    # & (overlap_exons > overlap_introns > downstream > upstream > intergenic)
    lnc_file = IOTools.openFile(merged_lncRNA_gtf)
    gene_class = {}  # dictionary of gene_id : classification
    input_transcripts = 0  # keep track of # transcripts in lncRNA_gtf
    for transcript in GTF.transcript_iterator(GTF.iterator(lnc_file)):
        input_transcripts += 1
        gene_id = transcript[0].gene_id
        strand = transcript[0].strand

        # create lists of indexed intervals that intersect transcript exons
        overlap_list = []
        intron_list = []
        plus_down_list = []
        minus_down_list = []
        plus_up_list = []
        minus_up_list = []
        for exon in transcript:
            if exon.contig in ref_index.mIndex.keys():
                overlap_list.extend([x for x in list(ref_index.get(exon.contig,
                                                                   exon.start,
                                                                   exon.end))])
            else:
                E.warn("Contig %s not in reference exon index "
                       "failed to retrieve intervals for %s" % (exon.contig,
                                                                exon.gene_id))
            if exon.contig in intron.mIndex.keys():
                intron_list.extend([x for x in list(intron.get(exon.contig,
                                                               exon.start,
                                                               exon.end))])
            else:
                E.warn("Contig %s not in reference intron index, "
                       "failed to retrieve intervals for %s" % (exon.contig,
                                                                exon.gene_id))
            if exon.contig in plus_down.mIndex.keys():
                plus_down_list.extend([x for x in list(plus_down.get(
                    exon.contig,
                    exon.start,
                    exon.end))])
            else:
                E.warn("Contig %s not in plus downstream index, "
                       "failed to retrieve intervals for %s" % (exon.contig,
                                                                exon.gene_id))
            if exon.contig in minus_down.mIndex.keys():
                minus_down_list.extend([x for x in list(minus_down.get(
                    exon.contig,
                    exon.start,
                    exon.end))])
            else:
                E.warn("Contig %s not in minus downstream index, "
                       "failed to retrieve intervals for %s" % (exon.contig,
                                                                exon.gene_id))
            if exon.contig in plus_up.mIndex.keys():
                plus_up_list.extend([x for x in list(plus_up.get(exon.contig,
                                                                 exon.start,
                                                                 exon.end))])
            else:
                E.warn("Contig %s not in plus upstream index, "
                       "failed to retrieve intervals for %s" % (exon.contig,
                                                                exon.gene_id))
            if exon.contig in minus_up.mIndex.keys():
                minus_up_list.extend([x for x in list(minus_up.get(exon.contig,
                                                                   exon.start,
                                                                   exon.end))])
            else:
                E.warn("Contig %s not in minus upstream index, "
                       "failed to retrieve intervals for %s" % (exon.contig,
                                                                exon.gene_id))

        # check if any exon in lncRNA intersects an reference exon
        if overlap_list:
            # if the intersecting exons are on the same strand,
            # classify lncRNA as sense.
            if strand in [x[2][6] for x in overlap_list]:
                gene_class[gene_id] = "sense"
                write_to_temp(temp_files[gene_class[gene_id]],
                              overlap_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # otherwise check if lncRNA has sense overlap with a reference
            # intron
            elif intron_list and strand in [x[2][6] for x in intron_list]:
                last = len(transcript) - 1
                start_list = [(x[0], x[1]) for x in list(intron.get(
                    transcript[0].contig,
                    transcript[0].start,
                    transcript[0].start + 1)) if x[2][6] == strand]
                end_list = [(x[0], x[1]) for x in list(intron.get(
                    transcript[last].contig,
                    transcript[last].end,
                    transcript[last].end + 1)) if x[2][6] == strand]
                # if start and end of transcript are within the same sense
                # introns, then lncRNA is classified as 'sense_intronic'
                if set(start_list) == set(end_list):
                    gene_class[gene_id] = "sense_intronic"
                    write_to_temp(temp_files[gene_class[gene_id]],
                                  intron_list,
                                  transcript)
                    temp_count[gene_class[gene_id]] += 1

                # if start/end are within different sense introns,
                # then lncRNA is classified as 'sense overlap'
                else:
                    gene_class[gene_id] = "sense_overlap"
                    write_to_temp(temp_files[gene_class[gene_id]],
                                  intron_list,
                                  transcript)
                    temp_count[gene_class[gene_id]] += 1

            # ...check if lncRNA is sense downstream on the plus strand...
            elif plus_down_list and strand in [x[2][6] for x in
                                               plus_down_list]:
                gene_class[gene_id] = "sense_downstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              plus_down_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # ...check if lncRNA is sense downstream on the minus strand...
            elif minus_down_list and strand in [x[2][6] for x in
                                                minus_down_list]:
                gene_class[gene_id] = "sense_downstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              minus_down_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # ...check if lncRNA is sense upstream on the plus strand...
            elif plus_up_list and strand in [x[2][6] for x in plus_up_list]:
                gene_class[gene_id] = "sense_upstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              plus_up_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # ...check if lncRNA is sense upstream on the minus strand...
            elif minus_up_list and strand in [x[2][6] for x in minus_up_list]:
                gene_class[gene_id] = "sense_upstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              minus_up_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # ...check if lncRNA is antisense downstream on the plus strand...
            elif plus_down_list and prioritize_flanks:
                gene_class[gene_id] = "antisense_downstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              plus_down_list,
                              transcript,
                              check_strand=False)
                temp_count[gene_class[gene_id]] += 1

            # ...check if lncRNA is antisense downstream on the minus strand...
            elif minus_down_list and prioritize_flanks:
                gene_class[gene_id] = "antisense_downstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              minus_down_list,
                              transcript,
                              check_strand=False)
                temp_count[gene_class[gene_id]] += 1

            # ...check if lncRNA is antisense upstream on the plus strand...
            elif plus_up_list and prioritize_flanks:
                gene_class[gene_id] = "antisense_upstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              plus_up_list,
                              transcript,
                              check_strand=False)
                temp_count[gene_class[gene_id]] += 1

            # ...check if lncRNA is antisense upstream on the minus strand...
            elif minus_up_list and prioritize_flanks:
                gene_class[gene_id] = "antisense_upstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              minus_up_list,
                              transcript,
                              check_strand=False)
                temp_count[gene_class[gene_id]] += 1

            # ...if none of the above... classify as antisense
            else:
                gene_class[gene_id] = "antisense"
                write_to_temp(temp_files[gene_class[gene_id]],
                              overlap_list,
                              transcript,
                              check_strand=False)
                temp_count[gene_class[gene_id]] += 1

        # if lncRNA doesn't intersect a reference exon,
        # check if it overlaps a reference intron
        elif intron_list:
            if strand in [x[2][6] for x in intron_list]:
                last = len(transcript) - 1
                start_list = [(x[0], x[1]) for x in list(intron.get(
                    transcript[0].contig,
                    transcript[0].start,
                    transcript[0].start + 1)) if x[2][6] == strand]
                end_list = [(x[0], x[1]) for x in list(intron.get(
                    transcript[last].contig,
                    transcript[last].end,
                    transcript[last].end + 1)) if x[2][6] == strand]
                # if start and end of transcript are within the same sense
                # introns, then lncRNA is classified as 'sense_intronic'
                if set(start_list) == set(end_list):
                    gene_class[gene_id] = "sense_intronic"
                    write_to_temp(temp_files[gene_class[gene_id]],
                                  intron_list,
                                  transcript)
                    temp_count[gene_class[gene_id]] += 1

                # if start/end are within different sense introns,
                # then lncRNA is classified as 'sense overlap'
                else:
                    gene_class[gene_id] = "sense_overlap"
                    write_to_temp(temp_files[gene_class[gene_id]],
                                  intron_list,
                                  transcript)
                    temp_count[gene_class[gene_id]] += 1

            # ...check if lncRNA is sense downstream on the plus strand...
            elif plus_down_list and strand in [x[2][6] for x in
                                               plus_down_list]:
                gene_class[gene_id] = "sense_downstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              plus_down_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # ...check if lncRNA is sense downstream on the minus strand...
            elif minus_down_list and strand in [x[2][6] for x in
                                                minus_down_list]:
                gene_class[gene_id] = "sense_downstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              minus_down_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # ...check if lncRNA is sense upstream on the plus strand...
            elif plus_up_list and strand in [x[2][6] for x in plus_up_list]:
                gene_class[gene_id] = "sense_upstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              plus_up_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # ...check if lncRNA is sense upstream in the minus strand...
            elif minus_up_list and strand in [x[2][6] for x in minus_up_list]:
                gene_class[gene_id] = "sense_upstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              minus_up_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # ...if none of the above, lncRNAs intersecting introns on
            # the opposite strand are classified as antisense
            else:
                if distinguish_antisense_introns:
                    gene_class[gene_id] = "antisense_intronic"
                else:
                    gene_class[gene_id] = "antisense"
                write_to_temp(temp_files[gene_class[gene_id]],
                              intron_list,
                              transcript,
                              check_strand=False)
                temp_count[gene_class[gene_id]] += 1

        # if lncRNA doesn't intersect reference introns or exons...
        # check if it's downstream on the plus strand...
        elif plus_down_list:
            # ... check if lncRNA is sense downstream...
            if strand in [x[2][6] for x in plus_down_list]:
                gene_class[gene_id] = "sense_downstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              plus_down_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # ...check if lncRNA is sense downstream on the minus strand...
            elif minus_down_list and strand in [x[2][6] for x in
                                                minus_down_list]:
                gene_class[gene_id] = "sense_downstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              minus_down_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # ...check if lncRNA is sense usptream on the plus strand...
            elif plus_up_list and strand in [x[2][6] for x in plus_up_list]:
                gene_class[gene_id] = "sense_upstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              plus_up_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # ...check if lncRNA is sense upstream on the pluse strand...
            elif minus_up_list and strand in [x[2][6] for x in minus_up_list]:
                gene_class[gene_id] = "sense_upstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              minus_up_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # if none of the above, lncRNA is classified as
            # antisense_downstream
            else:
                gene_class[gene_id] = "antisense_downstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              plus_down_list,
                              transcript,
                              check_strand=False)
                temp_count[gene_class[gene_id]] += 1

        # check if lncRNA is downstream on the minus strand...
        elif minus_down_list:
            # check if lncRNA is sense downstream
            if strand in [x[2][6] for x in minus_down_list]:
                gene_class[gene_id] = "sense_downstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              minus_down_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # ...check if lncRNA is sense upstream on the plus strand...
            elif plus_up_list and strand in [x[2][6] for x in plus_up_list]:
                gene_class[gene_id] = "sense_upstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              plus_up_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # ...check if lncRNA is sense upstream on the minus strand...
            elif minus_up_list and strand in [x[2][6] for x in minus_up_list]:
                gene_class[gene_id] = "sense_upstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              minus_up_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # if none of the above, lncRNA is classified as
            # antisense_downstream
            else:
                gene_class[gene_id] = "antisense_downstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              minus_down_list,
                              transcript,
                              check_strand=False)
                temp_count[gene_class[gene_id]] += 1

        # check if lncRNA is upstream on the plus strand...
        elif plus_up_list:
            # check if lncRNA is sense upstream...
            if strand in [x[2][6] for x in plus_up_list]:
                gene_class[gene_id] = "sense_upstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              plus_up_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # ...check if lncRNA is sense upstream on the plus strand...
            elif minus_up_list and strand in [x[2][6] for x in minus_up_list]:
                gene_class[gene_id] = "sense_upstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              minus_up_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # if none of the above, lncRNA is classified as
            # antisense upstream
            else:
                gene_class[gene_id] = "antisense_upstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              plus_up_list,
                              transcript,
                              check_strand=False)
                temp_count[gene_class[gene_id]] += 1

        # check if lncRNA is upstream on the minus strand...
        elif minus_up_list:
            # check if lncRNA is sense upstream...
            if strand in [x[2][6] for x in minus_up_list]:
                gene_class[gene_id] = "sense_upstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              minus_up_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # otherwise classify as antisense upstream
            else:
                gene_class[gene_id] = "antisense_upstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              minus_up_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

        # lncRNA that do not fall into any of the above categories
        # are classified as intergenic
        else:
            gene_class[gene_id] = "intergenic"
            temp_files[gene_class[gene_id]].write(str(transcript[0]) + "\n")
            temp_count[gene_class[gene_id]] += 1

    # check that all the numbers add up
    E.info("Number of lncRNA loci falling into each category are as follows:")
    for key, value in temp_count.iteritems():
        print(key + "\t" + str(value))
    total_classified = sum(temp_count.values())
    E.info("Total number of lncRNA loci classified: %i" % total_classified)
    E.info("Total number of lncRNA loci in input gtf: %i" % input_transcripts)

    # sanity check:
    assert total_classified == input_transcripts, (
        "Not all lncRNAs in input gtf were successfully classified")

    # close the tempfiles
    for handle in temp_file_names:
        temp_files[handle].close()

    # write the genes plus their classification to the outfile
    outf = IOTools.openFile(outfile, "w")
    for gtf in GTF.iterator(IOTools.openFile(lncRNA_gtf)):
        if gtf.gene_id in gene_class:
            gtf.source = gene_class[gtf.gene_id]
            outf.write(str(gtf) + "\n")
        else:
            E.info("Warning the gene_id %s is not classified" % gtf.gene_id)
    outf.close()

    os.unlink(merged_lncRNA_gtf)

    return tempdir


#################################################################################
#################################################################################
#################################################################################
## section: Cuffdiff expression calculator
#################################################################################

def calculateSummaryCuffdiffFPKM( cuffdiff_fpkms, 
                                  outfile, 
                                  id_dict,
                                  stat="mean", 
                                  num_samples = False,
                                  to_exclude = False ):
    """
    Designed to be passed a genes.read_group_tracking file from cuffdiff. 
    Calculates the median or mean for each gene in each cell type
    To exclude includes a list of samples to be ignored when calculating summary - 
    uses WGCNA syntax for naming of samples (tissue.condition.replicate)
    """
    # open a logfile
    logfile = outfile + ".log" 
    outf_log = IOTools.openFile( logfile, "w" )
    n_fail = 0

    # generate a list of samples to be excluded
    ids_to_remove = []
    if to_exclude:
        for sample in to_exclude:
            # samples in list of those to exclude have WGCNA nomenclature
            sample = re.sub( "\.", "-", sample )
            # id_dict has cuffdiff_id as key, sample id as value
            for cuffdiff_id, sample_id in id_dict.iteritems():
                if sample == sample_id:
                    ids_to_remove.append(cuffdiff_id)
                else:
                    continue
    
    # add each fpkm value to a dictionary with the gene/cell_type as key
    fpkm_dict = collections.defaultdict( list )
    for line in IOTools.openFile( cuffdiff_fpkms ).readlines():
        if line.startswith( "tracking_id" ):
            outf_log.write( line )
            continue
        line = line.split()
        if to_exclude:
            cuffdiff_id = line[1] + "_" + line[2]
            if cuffdiff_id in ids_to_remove:
                outf_log.write( "\t".join( line ) + "\n" )
                continue
        if line[8] != "OK":
            # E.warn( "Skipping %s in sample %s %s because it has"
            #         " cuffdiff status %s" % (line[0], line[1], line[2], line[8] ) )
            outf_log.write( "\t".join( line ) + "\n" )
            n_fail += 1
            gene_sample_id = line[0] + "__" + line[1]
            fpkm = ( np.nan )
            fpkm_dict[ gene_sample_id ].append( fpkm )
        else:
            # add remaining sample fpkm values to fpkm_dict
            gene_sample_id = line[0] + "__" + line[1]
            fpkm = float( line[6] )
            fpkm_dict[ gene_sample_id ].append( fpkm )

    outf_log.write( "\nDeleted gene_ids:" )

    # now create a nested dictionary with gene_ids the key and
    # cell_type: summary_fpkm as value
    summary_dict = collections.defaultdict( dict )
    for key, value in fpkm_dict.iteritems():
        gene_id, cell_type = key.split( "__" )
        if stat == "mean":
            summary_fpkm = np.nanmean( value )
        elif stat == "median":
            summary_fpkm = np.nanmedian( value )
        else:
            raise ValueError( "Unrecognised summary statistic" )
        summary_dict[ gene_id ][ cell_type ] = summary_fpkm


    # Note: this is obsolete because the number of samples will always be
    # max now that missing samples are assigned np.nan (above)
    if num_samples:
        to_remove = []
        for key, value in summary_dict.iteritems():
            if len( value.keys() ) < num_samples:
                E.warn( "Removing gene %s from fpkms because cuffdiff failed"
                        " to calculate fpkms for any samples in one or more"
                        " cell type -- see logfile %s" % ( key, logfile ) )
                outf_log.write( key + str( value ) )
                to_remove.append(  key )
        for x in to_remove:
            del summary_dict[ x ]

    # write the summary_fpkms to a flat file using OrderedDict
    outf = IOTools.openFile( outfile, "w" )
    start = True
    for gene, nested_dict  in summary_dict.iteritems():
        # from the collections docs...
        value = collections.OrderedDict( sorted( nested_dict.items(), 
                                                 key=lambda t: t[0] ) )
        if start:
            outf.write( "gene_id" + "\t" + "\t".join( value.keys() ) + "\n" )
            outf.write( gene + "\t" 
                        + "\t".join( [ str(x) for x in value.values() ] ) + "\n" )
            start = False
        else:
            outf.write( gene + "\t" 
                        + "\t".join( [ str(x) for x in value.values() ] ) + "\n" )
    outf.close()
    outf_log.close()
    E.info( "There were %s sample values which were either FAIL or HIDATA" % str( n_fail ) )


def mapCuffdiffIDs( infile, suffix=".bam" ):
    """ 
    Takes a cuffdiff read_groups.info.gz file and returns a dictionary mapping the 
    sample_id to the condition_replicate values
    """
    id_dict = {}
    header = True
    for line in IOTools.openFile( infile ).readlines():
        if header:
            header = False
            continue
        line = line.split()
        sample_id = P.snip( os.path.basename(line[0]), ".bam" )
        cuffdiff_id = "_".join( line[1:3] )
        id_dict[ cuffdiff_id ] = sample_id
    
    return id_dict


def extractPerSampleCuffdiffFPKM( cuffdiff_fpkms, outfile, id_dict = None ):
    """
    Designed to be passed a genes.read_group_tracking file from cuffdiff
    Outputs the fpkm value for each sample in the file 
    Cuffdiff IDs are mapped to original sample ids using dictionary, if provided
    """

    # extract fpkm values from cuffdiff output, checking that cuffdiff completed
    fpkm_dict = collections.defaultdict( dict )
    for line in IOTools.openFile( cuffdiff_fpkms ).readlines():
        if line.startswith( "tracking_id" ): continue
        line = line.split()
        gene_id = line[0]
        sample_id = line[1] + "_" + line[2]
        if id_dict:
            sample_id = id_dict[ sample_id ]
        if line[8] == "OK":
            fpkm = line[6]
        else:
            fpkm = "NULL"
        fpkm_dict[ gene_id ][ sample_id ] = fpkm
    
    # write the fpkm values to a flatfile
    outf = IOTools.openFile( outfile, "w" )
    start = True
    for gene, nested_dict in fpkm_dict.iteritems():
        value = collections.OrderedDict( sorted( nested_dict.items(), 
                                                 key=lambda t: t[0] ) )
        if start:
            outf.write( "gene_id" + "\t" + "\t".join( value.keys() ) + "\n" )
            outf.write( gene + "\t" + "\t".join( [ str(x) for x in value.values() ] ) + "\n" )
            start = False
        else:
            outf.write( gene + "\t" + "\t".join( [ str(x) for x in value.values() ] ) + "\n" )

    outf.close()


def extractPerSampleCuffdiffFPKM_stacked( cuffdiff_fpkms, outfile, id_dict ):
    """
    Designed to be passed a genes.read_group_tracking file from cuffdiff
    Modifies the replicate column to match hard-coded replicates.
    """

    # obsolete... as this assumption was invalid
    # # define mapping on the assumption that order of cuffdiff replicates equals
    # # the order of the replicate numbers in the original samples.
    # mapping = { "pro": [ "R1", "R2", "R3", "R4", "R5" ], 
    #             "pre": [ "R1", "R2", "R3", "R4", "R5" ],
    #             "immature": [ "R1", "R2", "R3", "R4", "R5" ],
    #             "mature": [ "R1", "R2", "R3", "R4", "R5" ],
    #             "follicular": [ "R11", "R12", "R13", "R14", "R15" ],
    #             "marginal": [ "R11", "R12", "R13", "R14", "R15" ],
    #             "b1a": [ "R1", "R2", "R3", "R4", "R5" ],
    #             "germinal": [ "R11", "R12", "R13", "R14", "R15" ] }


    # parse the input file, output tracking_id, condition, replicate, raw_frags, FPKM, status
    outf = IOTools.openFile( outfile, "w" )
    outf.write( "gene_id\t" # naming must match gtf2table output
                "cell_type\t"
                "replicate\t"
                "cuffdiff_raw_frags\t"
                "cuffdiff_FPKM\t"
                "cuffdiff_status\n" )
    header = True
    for line in IOTools.openFile( cuffdiff_fpkms ).readlines():
        if header:
            header = False
            continue
        line = line.split()
        if id_dict:
            tracking_id = id_dict[ "_".join(line[1:3]) ]
            condition = tracking_id.split("-")[1]
            replicate = tracking_id.split("-")[2]
        else:
            tracking_id = line[0]
            condition = line[1]
            replicate = line[2]
        raw_frags = line[3]
        fpkm = line[6]
        status = line[8]
        line_out = [ tracking_id, condition, replicate, raw_frags, fpkm, status ]
        outf.write( "\t".join( line_out ) + "\n" )
    outf.close()


def stackCoverage_readpairs( ref_gtf, infiles, tables, outfile ):
    """
    Take gene_id, antisense & sense mapping read counts and concatenate for 
    each sample.
    """
    # generate dictionary of strand for all gene_ids
    strand_dict = {}
    for gtfs in GTF.flat_gene_iterator( GTF.iterator( IOTools.openFile( ref_gtf ) ) ):
        gene_id = gtfs[0].gene_id
        strand = gtfs[0].strand
        strand_dict[gene_id] = strand

    # generate tempfile containing stacked coverage data for either strand
    to_cluster = False
    tmpf = P.getTempFilename(".")
    statement = ( "python %(scriptsdir)s/combine_tables.py"
                  "  --cat=CAT"
                  "  --log=%(outfile)s.log"
                  " %(tables)s |"
                  " cut -f1,2,3,9"
                  " > %(tmpf)s" )
    P.run()

    # open outfile
    outf = IOTools.openFile( outfile, "w" )
    outf.write( "gene_id\t"
                "strand\t"
                "cell_type\t"
                "replicate\t"
                "sense\t"
                "antisense\t"
                "pptn_sense\t"
                "pptn_plus\n" )

    # post process file
    header = True
    for line in IOTools.openFile( tmpf ).readlines():
        if header:
            header = False
            continue
        filename, gene_id, antisense, sense = line.split()
        sample_id = os.path.basename( filename ).split("_")[0]
        tissue, cell_type, replicate = sample_id.split("-")

        assert gene_id in strand_dict.keys(), "Missing gene_id: %s" % gene_id
        strand = strand_dict[ gene_id ]

        try: 
            pptn_sense = float(sense) / ( float(sense) + float(antisense) )
            if strand == "+":
                pptn_plus = pptn_sense
            elif strand == "-":
                pptn_plus = 1 - pptn_sense
            else:
                raise ValueError( "Unrecognised strand: %s" % strand )
        except ZeroDivisionError:
            # zero coverage on either strand
            pptn_sense = "NULL"
            pptn_plus = "NULL"

        line_out = [ gene_id, 
                     strand, 
                     cell_type, 
                     replicate, 
                     sense, 
                     antisense, 
                     pptn_sense, 
                     pptn_plus ]
        line_out = [ str(x) for x in line_out ]
        outf.write( "\t".join( line_out ) + "\n" )

    outf.close()
    os.unlink( tmpf )


#################################################################################
#################################################################################
#################################################################################
## section: Filter pandas dataframes
#################################################################################
# Mostly functions that receive an array and return a boolean, dependant on some 
# condition. They are designed to be used when subsetting pandas dataframes with
# pd.apply().
#################################################################################
# Utility functions...
# calculate robust parameters from an array (i.e. ignore NaN entries)
# N.B. these functions will be obsolete from numpy 1.8 onwards (see np.nanmean)
def robust_mean( array, keep_dims = False ):
    return np.mean( [ x for x in array if not np.isnan(x) ], keepdims = keep_dims )
 

def robust_std( array, keep_dims = False ):
    return np.std( [ x for x in array if not np.isnan(x) ], keepdims = keep_dims )

def robust_median( array, keep_dims = False ):
    return np.median( [ x for x in array if not np.isnan(x) ] )

def maskNA( series, 
              missing = ("na", "Nan", None, "", "NA", "NULL", np.nan), 
              dtype = np.float ):
    mask = [ i in missing for i in series ]
    return np.array( [ series[i] for i in range(len(series)) if not mask[i]], dtype = dtype )

# convert proportion thresholds to integer values
def set_threshold( threshold, length ):
    """
    Recieves a threshold value and an array length. 
    If threshold value is less than 1.0, then it is treated as a proportion; the 
    minimum number of elements in array required to meet this proportion is 
    calculated and returned as an integer value. 
    If the threshold value is greater than 1, it is returned as an integer value.
    """
    threshold = float( threshold )
    if threshold < 1.0:
        threshold = np.ceil( threshold * length )
    return int( threshold )

#################################################################################
def record_outliers( df, outfile, std_threshold = 5 ):
    """
    Takes a dataframe and calculates summary statistics across rows.
    Subsets those rows in which one or more sample values lie greater than 
    'std_threshold' standard deviations from the mean. Writes the subset 
    dataframe to an outfile. 
    """

    # determine whether an array has outliers
    def contains_outliers( array, std_threshold ):
        mean = robust_mean( array )
        std = robust_std( array )
        contains_outlier = False
        for i in array:
            if np.absolute( i - mean )/ std >= std_threshold:
                contains_outlier = True
                break
            else: continue
        return contains_outlier

    # select those rows containing arrays with one or more outlier
    series = df.apply( contains_outliers, axis = 1, args=( std_threshold, ) ) 
    df_out = df[ series ]
    df_out.to_csv( outfile, sep ="\t", na_rep = "NaN" )

    ## bonus output
    ## write out series used to subset the datafame
    #series_out = outfile + "_series"
    #series_out = IOTools.openFile( series_out, "w" )
    #series_out.write( str( series ) )
    #series_out.close()

    ## write zscores
    #z_out = outfile + "_zscores"
    #df_zscore = df.apply( stats.mstats.zscore, axis = 1 )
    #df_zscore.to_csv( z_out, sep= "\t", na_rep = "NaN" )


def mask_outliers( array, std_thresh = 5 ):
    """
    Receives either numpy 1darray or pandas series. 
    Returns an array in which any value greater than 'std_threshold' standard 
    deviations from the mean is replaced with 'nan'.
    Note: Ignores 'nan' values in the input series...
     ...numpy's nanmean/nanstd are available in version 1.8 onwards
     ...scipy.stats.mstats.zscore doesn't have option to ignore 'nan'
    """
    mean = robust_mean( array )
    std = robust_std( array )

    def filter_outliers( x, mean, std, std_thresh ):
        if std_thresh and np.absolute( x - mean )/ std >= std_thresh:
            return float( "NaN" )
        else:
            return x

    v_filter_outliers = np.vectorize( filter_outliers )
    
    return v_filter_outliers( array, mean, std, std_thresh )


def drop_nan( array, data_thresh = 1, rejected = False):
    """
    Receives an array. Returns True if more than 'data_thresh' values are not
    'NaN'. Where 'data_thresh' is an integer value.
    Note: pandas.DataFrame.dropna does not have the option to retain inverse
    """
    has_data = 0
    for i in array:
        if not np.isnan(i):
            has_data += 1
        else: 
            continue

    if has_data >= data_thresh:
        keep = True
    else:
        keep = False
        
    if rejected:
        keep = not keep

    return keep


# return true if greater than x samples have fpkm > y
def drop_fpkm( array, min_fpkm, fpkm_thresh, rejected = False ):
    """
    Receives an array of FPKM values.
    Returns True if there are greater than 'fpkm_thresh' entries in an array
    for which FPKM value is above 'min_fpkm'. Where 'fpkm_thresh' is an integer
    value.
    """
    high_fpkm = 0
    for i in array:
        if i >= min_fpkm:
            high_fpkm += 1
        else: 
            continue

    if high_fpkm >= fpkm_thresh:
        keep = True
    else: 
        keep = False

    if rejected: 
        keep = not keep

    return keep
    

def drop_cv( array, cv_thresh, rejected = False ):
    """
    Returns True if the coefficient of variation of an array is greater than the
    specified threshold
    Note: Ignores nan values in the input array.
    """
    mean = robust_mean( array )
    std = robust_std( array )
    cv = std/mean

    if cv >= cv_thresh:
        keep = True
    else:
        keep = False
    
    if rejected:
        keep = not keep

    return keep


def calc_cv( array ):
    """Calculates the coefficient of variation for an array"""
    mean = robust_mean( array )
    std = robust_std( array )
    return std/mean


def calc_zscore( array ):
    """Calculates standard deviation, ignoring NaN values"""
    mean = robust_mean( array )
    centered = array - mean
    stddev = robust_std( centered )
    zscore = centered/stddev

    return zscore


def returnHighestValue( array, min_thresh ):
    """
    Receives an array of values 
    Returns a boolean array, where max value is true if greater than min_thresh,
    else false, and all other values are False
    """

    if max(abs(array)) >= min_thresh:
        bool_array = abs(array).isin([max(abs(array)),])
        assert len( array[bool_array] ) == 1, "Two modules have same correlation score"
        
        return bool_array
    else:
        bool_array = array.copy()
        for i in range(len(bool_array)):
            bool_array.iloc[i] = False
        
        return bool_array    

################################################################################
## section: alter cell value based on row/column statistics
################################################################################
def zeroAdjust( series ): 
    """
    Recieves a series, if the minimum value is negative, zero adjusts data.
    np.nan values are returned as nan.
    """
    min_val = min( series)
    if min_val < 0.0:
        def zadj( x, min_val ):
            return x + abs( min_val )
        v_zadj = np.vectorize( zadj )

        return v_zadj( series, min_val )
    else:
        return( array )


#################################################################################
## section: Summarize Axis functions
#################################################################################
def summarizeAxis( df, function, subset=None, axis=1 ):
    """
    A wrapper for pandas apply, which allows selection of a subset of indices
    from the  axis to be summarized.
    Applies a pre-defined function to a pandas dataframe, returning a series.
    subset -> a subset of index labels for the axis to be condensed
    ??How to pass on kwargs??
    """
    assert axis in [0,1], "axis must be either 0 (columns) or 1 (rows)"
    # axis to be summarized needs to be rows for subsetting to work
    if axis == 0: df = df.T
    
    # subsetting selects columns
    if subset:   
        try:
            df = df[subset.split(",")] 
        except IOError:
            print( "Unable to subset dataframe index (%s) using values given"
                   " %s" % ( ",".join( [x for x in df.columns] ), subset ) )

    return df.apply( function, axis=1 )


def calcSpecificity( series, mask_na = True ):
    """
    Calculates the specificity of gene expression, following Yani et al. 2005
    Bioinformatics 21, 650-659.   
    Note this gives ZeroDivisionError if all values in series are 0, but that 
    return of ZeroDivisionError doesn't cause pd.apply to fail... it just reports
    nan.
    """
    if mask_na:
        series = maskNA( series )
    max_val = max( series )

    return sum( [ 1 - x/max_val for x in series ] ) / ( len( series ) -1 )

#################################################################################
#################################################################################
#################################################################################
## section: Return venn coverage stats from pandas datafames
#################################################################################

def getNDVenn( df, stages, threshold, greater_than=True ):
    """
    Generates output compatible with R package VennDiagram
    Takes i) a pandas data frame with variables as rows and samples as columns, 
    ii) a list of colnames to compare, iii) a threshold specifying the value at 
    which variables are considered true or false. Returns a pandas.series 
    specifying the number of values that are true for every combination.
    """
    coverage_values = pd.Series()
    
    # check names
    for name in stages:
        assert name in df.axes[1], "Name '%s' not in dataframe colnames" % name
    
    # generate combinations
    combinations = []
    for i in range( 1, len( stages ) + 1 ):
        for j in itertools.combinations( sorted( stages ), i ):
            combinations.append( j )
            
    # generate combination headers
    combination_headers = [ "_".join( x ) for x in combinations ]
    
    # create function to generate combination specific statements
    def get_statement( combination, threshold, greater=greater_than ):
        statement = []
        for i in combination:
            if greater:
                 statement.append( "(" + str(i) + ">" + str(threshold) + ")" )
            else:
                 statement.append( "(" + str(i) + "<" + str(threshold) + ")" )
        return "&".join( statement )

    for combination in enumerate( combinations ):
        df_t = df.query( get_statement( combination[1], threshold ) )
        coverage = len( df_t.index )
        header = combination_headers[ combination[0] ]
        coverage_values.set_value( header, coverage )
        
    return coverage_values


def plotNDVenn( df, combinations, threshold, outf_stub, biotype ):
    """
    For each of the lists provided in combinations, 
    at the given threshold,retrieve the overlap between list elements from df,
    then plot the appropriate venn diagram
    """
    # iterate through each B cell combination, and run get NDVenn to retrieve 
    # series to be plotted.
    R('''suppressMessages(library('VennDiagram'))''')
    for combination in combinations:
        # generate outfile name
        outf_suffix = "_".join( combination ) + ".pdf"
        outf = outf_stub + "_" + biotype + "_" + outf_suffix
        outf_2 = IOTools.openFile( P.snip( outf, ".pdf" ) + ".tsv", "w" )

        # collect overlap using getNDVenn()
        overlap = getNDVenn( df, combination, threshold )
        outf_2.write("For biotype & combination:\n%s %s\n%s" % (biotype, str(combination), str(overlap)) )
        outf_2.close()

        # I can't work out how to save robjects from python... pushing to R env
        R('''rm(list=ls())''')
        R.assign( "ov", robjects.IntVector( overlap.tolist() ) )
        R.assign("cbn", robjects.Vector( overlap.index.tolist() ) )

        if len(combination) == 2:
            # plot pairwise venn
            R('''pdf("%(outf)s")''' % locals() )
            R('''draw.pairwise.venn(area1=ov[1], area2=ov[2], cross.area=ov[3], category=cbn[1:2])''')
            R('''dev.off()''')
        elif len( combination ) == 3:
            # plot three-way venn
            R('''pdf("%(outf)s")''' % locals() )
            R('''draw.triple.venn(area1=ov[1], area2=ov[2], area3=ov[3], n12=ov[4],
                                  n13=ov[5], n23=ov[6], n123=ov[7], category=cbn[1:3])''')
            R('''dev.off()''')
        elif len( combination ) == 4:
            # plot four-way venn
            R('''pdf("%(outf)s")''' % locals() )
            R('''draw.quad.venn(area1=ov[1], area2=ov[2], area3=ov[3], area4=ov[4], 
                                n12=ov[5], n13=ov[6], n14=ov[7], 
                                n23=ov[8], n24=ov[9], n34=ov[10],
                                n123=ov[11], n124=ov[12], n134=ov[13], n234=ov[14], n1234=ov[15],
                                category=cbn[1:4])''')
            R('''dev.off()''')
        else:
            raise Exception( "Not capable of plotting a %i-way venn" % len(combination) )

    R('''rm(list=ls())''')


def plotVennDiagrams( count_table, 
                      biotype_table,
                      outf_stub, 
                      biotype, 
                      combinations, 
                      threshold ):
    """
    Retrieve count data matching the specified biotype from csvdb.
    Pass resulting dataframe to plotNDVenn
    """
    # Retrieve relevant data...
    if biotype == "lncRNA_all":
        statement = ( "SELECT a.*"
                      " FROM %(count_table)s AS a"
                      " INNER JOIN %(biotype_table)s AS b"
                      " ON a.gene_id = b.gene_id"
                      " WHERE b.biotype != 'protein_coding'" % locals() )
    else:
        statement = ( "SELECT a.*"
                      " FROM %(count_table)s AS a"
                      " INNER JOIN %(biotype_table)s AS b"
                      " ON a.gene_id = b.gene_id"
                      " WHERE b.biotype = '%(biotype)s'" % locals() )

    df = PU.fetch_DataFrame( statement )
    df = df.set_index( "gene_id" )

    # Hack: need to make the order of cell types in dataframe match ontogeny
    new_names = []
    for i in df.columns.tolist():
        if i == "pro":
            new_names.append("apro")
        elif i == "pre":
            new_names.append("cpre")
        elif i == "immature":
            new_names.append("dimmature")
        elif i == "mature":
            new_names.append("bmature")
        elif i == "b1a":
            new_names.append("eb1a")
        elif i == "follicular":
            new_names.append("gfollicular")
        elif i == "marginal":
            new_names.append("hmarginal")
        elif i == "germinal":
            new_names.append("fgerminal")
        else:
            raise Exception( "Unrecognized header %s" " ".join( df.columns.tolist() ) )

    df.columns = new_names
    E.info("New header names in df: %s" % " ".join(df.columns.tolist()))

    plotNDVenn( df, combinations, threshold, outf_stub, biotype )


#################################################################################
#################################################################################
#################################################################################
## section: Miscellaneous statistical functions
#################################################################################

def filterMasked( xvals, 
                  yvals, 
                  missing = ("na", "Nan", None, "", "NA", "NULL"), 
                  dtype = np.float, 
                  mask_nan = False ):
    """
    Converts two lists of paired values into numpy arrays,
    removing any element that has values missing in one or both lists.
    """
    # create boolean list identifying values missing in first and second array
    if mask_nan:
        mask = [ i in missing or np.isnan(i) for i in xvals ] 
        ymask = [ i in missing or np.isnan(i) for i in yvals ]
    else: 
        mask = [ i in missing for i in xvals ]
        ymask = [ i in missing for i in yvals ]
    # create a consensus boolean list identifying samples missing in one or both arrays
    for i in range( len( xvals ) ):
        if ymask[i]: mask[i] = True
    return ( np.array( [xvals[i] for i in range(len(xvals)) if not mask[i]], dtype = dtype  ),
             np.array( [yvals[i] for i in range(len(yvals)) if not mask[i]], dtype = dtype) )

"""
def calculateCorrelations(infile, 
                          outf_1, 
                          outf_2, 
                          method="pearson", 
                          min_obs=5,
                          regex1="ENSMUSG",
                          regex2="LNC", 
                          header=True ):
"""

def calculateCorrelations( params ):
    """
    This function takes an infile and subdivides rows based on first column 
    matching one of two supplied regular expressions (rows that contain neither
    regex will cause an error). Correlation stats are then calculated for each
    pairwise comparison between the two subdivided groups.
    Two outfiles are generated, the first containing summary statistics, the
    second contianing a matrix of correlation coefficients.
    """

    # P.submit takes only a string... 
    infile, outf_1, outf_2, method, min_obs, regex1, regex2, header = params
    if header == "True":
        header = True
    else:
        header = False
    min_obs = int( min_obs )
    
    # assign genes to separate dictionaries based on supplied regex
    coding_dict = {}
    lncRNA_dict = {}
    for line in IOTools.openFile( infile ).readlines():
        if header:
            header = False
            continue
        line= line.split()
        key = line[0]
        values = line[1:]
        if re.search( regex1, key ):
            coding_dict[ key ] = values
        elif re.search( regex2, key ):
            lncRNA_dict[ key ] = values
        else:
            raise ValueError( "Infile contains row that doesn't"
                              " match either regex" )

    # create nested dict to store results summary
    nested_dict = collections.defaultdict( dict )

    # open file to contain correlation statistics
    outf1 = IOTools.openFile( outf_1, "w" )
    outf1.write( "Gene_id_1\tGene_id_2\tpvalue\tmethod\t"
                 "no_observations\tcoefficient\n" )

    # iterate through all pairwise combinations of entries in files,
    # calculate correlation for each combination
    for combination in itertools.product( lncRNA_dict.keys(), 
                                          coding_dict.keys() ):
        lnc_key, coding_key = combination
        lnc_values = lncRNA_dict[ lnc_key ]
        coding_values = coding_dict[ coding_key ]

        lnc_values, coding_values = filterMasked( lnc_values, coding_values )
        
        assert len( lnc_values ) == len( coding_values )

        if method == "pearson":
            result = stats.pearsonr( lnc_values, coding_values )
        elif method == "spearman":
            result = stats.spearmanr( lnc_values, coding_values )
        else:
            raise ValueError( "Unknown method %s" % method ) 

        coefficient, pvalue = result
        nobservations = len( lnc_values )

        # write stats to outfile of correlation statistics
        output = [ lnc_key, 
                   coding_key, 
                   pvalue, 
                   method, 
                   nobservations, 
                   coefficient ] 
        output_string = "\t".join( [ str(x) for x in output ] )
        outf1.write( output_string + "\n" )

        # add the correlation coefficient to summary dictionary, 
        # unless there were < min_observations
        if nobservations < min_obs:
            nested_dict[ lnc_key ][ coding_key ] = "NA"
            E.warn( "Samples %s vs %s only had %s"
                    " observations:" % (lnc_key, coding_key, nobservations) )
        else:
            nested_dict[ lnc_key ][ coding_key ] = coefficient

    # create a sorted list of nested keys (coding_ids) to use as headers
    header_list = sorted( nested_dict[ nested_dict.keys()[0] ].keys() )
    headers = "\t".join( [ x for x in header_list ] )
    # open file to contain matrix of correlation coefficients
    outf2 = IOTools.openFile( outf_2, "w" )
    outf2.write( "Gene_id\t" + headers + "\n" )

    # iterate through nested dict and output the coefficient values
    # in the correct order
    for lnc_id in nested_dict.keys():
        coeff_list = []
        for coding_id in header_list:
            coeff_list.append( nested_dict[ lnc_id ][ coding_id ] )
        coeff_string = "\t".join( [ str(x) for x in coeff_list ] )
        outf2.write( lnc_id + "\t" + coeff_string + "\n" )

    outf1.close()
    outf2.close()


def calculateCorrelations_df( in_df1, in_df2, min_n, outfile, outfile_rej ):
    """
    This function receives two dataframes, in which genes are rows and samples
    are columns. It asserts that all the columns in df1 are also present in df2
    (but not the reciprocal). 
    Writes an outfile in which df1_genes are columns, df2_genes are rows and
    the values are the pearson correlation.
    Samples with values missing for one or other gene are masked. Pairwise 
    comparisons with too many missing values (min_n) are rejected and written
    to outfile as NULL.
    """

    corr_dict = collections.defaultdict( dict )
    outf = IOTools.openFile( outfile, "w" )
    outf_rej = IOTools.openFile( outfile_rej, "w" ) 
    outf_rej.write( "df1_id\tdf2_id\tn_samples\n" )
    
    # first dataframe has samples as columns and modules as rows
    # assert that all the sample_ids in in_df1 are also in in_df2
    for sample_id in in_df1.columns:
        assert sample_id in in_df2.columns, "Sample %s missing from second dataframe" % sample_id

    # iterate through modules
    for module in in_df1.index:
        eigen_ser = in_df1.loc[ module ]
        # create a list of eigen_vals
        eigen_vals = []
        for sample_id in eigen_ser.index:
            eigen_vals.append( eigen_ser.loc[ sample_id ] )


        # second dataframe has samples as columns and gene_ids as rows
        # iterate through gene_ids
        for gene in in_df2.index:
            gene_ser = in_df2.loc[ gene ]
            gene_vals = []
            # create list of gene_vals
            for sample_id in eigen_ser.index:
                gene_vals.append( gene_ser.loc[ sample_id ] )

            # mask NULL values, rejecting those where there are too few samples
            eigen_v, gene_v = filterMasked( eigen_vals, gene_vals, mask_nan = True )
            if len( eigen_v ) < min_n:
                outf_rej.write( gene + "\t" + str( len( eigen_v ) ) + "\n" )
                corr_dict[ gene ][ module ] = "NULL" 
                continue
            else:
                result = stats.pearsonr( eigen_v, gene_v )
                corr_dict[ gene ][ module ] = result[0]

    # write the dictionary of correlations to outfile
    module_ids = [ x for x in in_df1.index ]
    outf.write( "gene_id\t" + "\t".join( module_ids ) + "\n" )
    for gene_id in corr_dict.iterkeys():
        outf.write( gene_id + "\t" )
        correlations = []
        for module_id in module_ids:
            correlations.append( corr_dict[ gene_id ][ module_id ] )
        outf.write( "\t".join( [ str(x) for x in correlations ] ) + "\n" )

    outf.close()
    outf_rej.close()

#################################################################################
#################################################################################
#################################################################################
## section: R functions
#################################################################################
R('''
suppressMessages( library( DESeq ) )
''')

def calcNormalisedCountData( infile, outfile ):
    """
    This function takes a file of read counts, with sample_ids as columns and
    gene_ids as rows. 
    Uses DESeq's estimateSizeFactorsForMatrix to normalise the count data on the
    basis of then
    outputs the normalsed data in the same tabular form. 
    """

    R('''
normalise.counts <- function( infile, outfile ){
  # read table as a dataframe
  table <- read.delim( infile, header=TRUE, sep="\t", row.names=1, stringsAsFactors=FALSE )
  # get the relative library sizes using DESeq
  sizes <- estimateSizeFactorsForMatrix(table)
  # create normalised table
  table.n <- table / do.call(rbind, rep(list(sizes), length(rownames(table))))
  
  # add rownames and colnames to the table
  rownames(table.n) <- rownames(table)
  colnames(table.n) <- colnames(table)

  # open outfile as gzip file
  outfile <- gzfile(outfile, "w")
  # write normalised table to outfile
  write.table(table.n, file=outfile, append=FALSE, quot e=FALSE, col.names=NA, row.names=TRUE, sep="\t")
  close(outfile)
}
''')
    norm_count_data = R['normalise.counts']
    norm_count_data( infile, outfile )


def calcRLENormalisedCountData( pddf, meth = "RLE"):
    """
    This function takes a pandas dataframe where genes are rows and samples 
    are columns and returns a pandas dataframe with normalised count values.
    Normalization is done using EdgeR RLE.
    """
    edgeR = importr("edgeR")

    # convert dataframe to R object
    df = com.convert_to_r_dataframe( pddf )
    # normalize data 
    n_factors = edgeR.calcNormFactors( df, method= meth )
    # it seems that edgeR is returning an ndarray anyway
    assert isinstance( n_factors, np.ndarray ), "edgeR is not returning an np.ndarray"
    # convert resulting vector back to a numpy array
    # n_factors = rpyn.ri2numpy( n_factors )
    # scale the dataframe
    df = pddf.mul( n_factors, axis=1 )
    
    return n_factors, df
    

def calcPoissonResiduals( infile, outfiles ):
    """
    This function takes a file of read counts, with sample_ids as columns and
    gene_ids as rows. 
    It fits a poisson glm to data for each gene_id and returns the standardized
    residuals in the same tabular form (see Crawley p.520).
    Data are normalised to account for different library sizes using DESeq's
    estimateSizeFactorsForMatrix() function.
    """

    outfile_1, outfile_2 = outfiles

    R('''
poisson.residuals <- function( infile, outfile.1, outfile.2 ){
  # read the infile as a dataframe
  table <- read.delim( infile, header=TRUE, sep="\t", row.names=1, stringsAsFactors=FALSE )
  # extract the cell types from the column names of this table
  cell.type <- as.factor(do.call(rbind, strsplit(colnames(table), ".", fixed=TRUE))[,2])
  
  # create a data frame of count data normalised using DESeq
  sizes <- estimateSizeFactorsForMatrix(table)
  # do.call performs rbind on all the lists produced by rep
  # rep(X, Y) returns X, Y times
  # divide the original table by the table of scaling factors
  table.n <- table / do.call(rbind, rep(list(sizes), length(rownames(table))))
  
  # create an empty data frame to store residuals
  residuals <- data.frame()
  # create an empty data frame to store the fitted values
  fitted.val <- data.frame()
  # create an empty character vector to contain gene_ids
  gene.ids <- character()
  
  # iterate through rows in table
  for (i in 1:length(rownames(table.n))) {
    # extract the gene_id from table
    gene.id <- rownames(table.n)[i]
    # append gene_id to list of gene_ids
    gene.ids <- append(gene.ids, gene.id)
    # extract the count data from table
    counts <- as.numeric(as.vector(table.n[i,]))
    # fit glm using sample condition as categorical predictor
    model <- glm(counts~cell.type)
    # calculate the standardised residuals for poisson errors
    std.res <- (model$y - model$fitted.values)/sqrt(model$fitted.values)
    # append standardised residuals to the data frame
    residuals <- rbind(residuals, std.res)

    # append the fitted values to the data frame
    fitted.val <- rbind(fitted.val, model$fitted.values)
  }
  
  # apply the column names to data frame of residuals and fitted values
  colnames(fitted.val) <- colnames(table)
  colnames(residuals) <- colnames(table)
  # apply the row gene_ids as rownames
  rownames(fitted.val) <- gene.ids
  rownames(residuals) <- gene.ids

  # open outfiles as gzfiles
  outfile.1 <- gzfile(outfile.1)
  outfile.2 <- gzfile(outfile.2)
  
  # write the resulting dataframe to an outfile
  # NB. col.names=NA creates a blank column name for columns when row.names=TRUE
  write.table(fitted.val, file=outfile.1, append=FALSE, quote=FALSE, col.names=NA, row.names=TRUE, sep="\t")
  write.table(residuals, file=outfile.2, append=FALSE, quote=FALSE, col.name=NA, row.names=TRUE, sep="\t")

  # close files
  close(outfile.1)
  close(outfile.2)
}
''' )
    poisson_residuals = R['poisson.residuals']
    poisson_residuals( infile, outfile_1, outfile_2 )


def calcExpressionThresholds( table_name, regex, threshold, outfile ):
    """
    This function was replaced with the python function threshold expression
    This function receives a table of expression values and, given a particular
    expression threshold, returns a boolean. Table is passed as tablename for
    sqlite database.
    Can't work out how to return rpy2 dataframe to python 
    (pandas.rpy.common.load_data)
    Writing flatfile from within r function instead. 
    """

    # define R function to threshold expression
    R('''
      threshold.expression <- function( df, threshold, outfile ){
        for( i in 1:ncol(df) ){
          df[,i] = df[,i] >= as.numeric( threshold ) }
        df = df*1
        outfile = gzfile( outfile )
        write.table( df, file=outfile, append=F, quote=F, col.names=T, row.names=T, sep="\t" )
        close( outfile )
        }
      ''')

    # use pipeline utilities to retrieve csvdb data as data frame
    statement = ( "SELECT * FROM %(table_name)s"
                  " WHERE gene_id LIKE '%%(regex)s%" % locals() )
    df = PU.fetch_DataFrame( statement )

    # convert dataframe to an R object
    df = pd.rpy.common.convert_to_r_dataframe(df)

    # run R function from python
    threshold_expression = R['threshold.expression']  
    threshold_expression( df, threshold, outfile )

#################################################################################
#################################################################################
#################################################################################
## section: access BioMart
#################################################################################
# Functions for accessing and retrieving data from biomart using python package
# biomart
#################################################################################
# 
# database = None
# dataset = "mm9_phase1_samples" or "mm9_phase1_peaks2"
# attributes = [] list of fields in order to be returned
# filters = {} nested dictionary of filter: values
# header = boolean
def getBiomartDataset( host = "http://fantom.gsc.riken.jp/5/biomart/",
                       dataset = "mm9_phase1_samples", 
                       attributes = [ "library_id", "ff_id", "sample_name" ], 
                       filters = {}, 
                       output_headers = False ):
    if output_headers:
        output_headers = 1
    else: 
        output_headers = 0

    # connect to host
    server = BiomartServer( host )

    # establish dataset
    dataset = server.datasets[ dataset ]

    # return data as tsv
    return dataset.search( params = { 'attributes': attributes,
                                      'filters': filters },
                           header = output_headers  )

    
#################################################################################
#################################################################################
#################################################################################
## section: WCGNA analysis (R)
#################################################################################
# These functions are wrappers for the fundamental steps of WGCNA analysis.
# However, instead of just calling the underlying R function, they also plot the
# related diagnostic plots and save pertinent data for subsequent pipeline steps.
#################################################################################
def createRkwargs( dictionary ):
    """
    An attempt to handle passing kwargs from python to R. String arguments still
    need to be double enclosed in quotations. 
    """
    for key, val in dictionary.iteritems():
        if isinstance( val, ( bool, ) ):
            if val: 
                dictionary[ key ] = "TRUE"
            else:
                dictionary[ key ] = "FALSE"
        else:
            continue
    kw_list = ", " + ", ".join( [ x + "=" + str(y) for x, y in dictionary.iteritems() ] )

    return kw_list


def wgcnaRunGoodSampleGenes( in_tsv, out_filtered, out_rejected ):
    """
    Reads a tsv file where genes are rows and samples are columns,
    runs goodSampleGenes, removes genes with expression that is too low, writes 
    two outfiles, one containing filtered genes, the other containing removed 
    genes. 
    """

    # read the infile as a dataframe
    in_tsv = os.path.abspath( in_tsv )
    R( '''datExpr <- read.delim( "%(in_tsv)s",
                                   header=TRUE,
                                   sep="\t", 
                                   row.names=1, 
                                   stringsAsFactors=FALSE ) ''' % locals() )
    E.info( "Refcoding FPKMs successfully loaded into R" )
    R( '''head( datExpr )''' )

    # transpose, so samples are rows and genes are columns
    R( '''datExpr <- t(datExpr)''' )

    # run goodSampleGenes
    E.info( "Running WGCNA goodSampleGenes:" )
    R( '''gsg <- goodSamplesGenes( datExpr, verbose = 5 )''' )

    # check to see if all data passed WGCNA quality thresholds
    all_ok = R('''toString(gsg$allOK)''')[0]
    if all_ok == "TRUE":
        E.info( "All genes/samples are suitable for WCGNA" )
        # touch rejected outfile
        P.touch( out_rejected )
    else:
        E.info( "Some genes/samples are unsuitable for WGCNA\n"
                "Filtering rejected data to: %s" % ( out_rejected ) )
        # removing rejected samples/genes
        R( '''datExpr.rejected <- datExpr[!gsg$goodSamples, !gsg$goodGenes]''' )
        R( '''datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]''' )
        # writing rejected samples to out_rejected
        R( '''out.rejected <- gzfile( "%(out_rejected)s", "w" ) ''' % locals() )
        R( '''write.table( datExpr.rejected, file=out.rejected, sep=",", quote=FALSE )''' % locals() )
        R( '''close( out.rejected )''' )

    # write the filtered and rejected data to flatfiles... misses title field for rownames
    R( '''out.filtered <- gzfile( "%(out_filtered)s", "w" ) ''' % locals() )
    R( '''write.table( datExpr, file=out.filtered, sep=",", quote=FALSE )''' % locals() )
    R( '''close( out.filtered )''' )

    # Debug
    #out_rdata = P.snip( out_rejected, ".gz" ) + ".rds"
    #R( '''save( datExpr, gsg, file="%(out_rdata)s" )''' % locals() )

def wgcnaClusterSamples( in_tsv, out_png, out_tree, out_dat ):
    """
    Reads a csv file where samples are rows and genes are columns, 
    Following WGCNA analysis, uses dist and flash clust to create a sampleTree.
    Plots the sample tree and saves it as an RDS file
    """    
    # read in dataframe from tsv file
    in_tsv = os.path.abspath( in_tsv )
    R( '''datExpr <- read.delim( "%(in_tsv)s",
                                  header=TRUE,
                                  sep=",", 
                                  row.names=1, 
                                  stringsAsFactors=FALSE ) ''' % locals() )
    # create sample dendrogram
    R( '''sampleTree <- flashClust( dist( datExpr ), method = "average" )''' )

    # plot sample dendrogram
    R( '''png( "%(out_png)s" )''' % locals() )
    R( '''plot( sampleTree, 
                main = "Sample clustering to detect outliers",
                cex.lab = 1.5,
                cex.axis = 1.5,
                cex.main = 2 )''' )

    R( '''dev.off()''' )
    
    # save the sample dendrogram as an .rds file
    R( '''saveRDS( sampleTree, file="%(out_tree)s" )''' % locals() )
    # save the expression dataframe as an .rds file
    R( '''saveRDS( datExpr, file="%(out_dat)s" )''' % locals() )

def wgcnaRemoveOutlierSamples( in_tree, 
                               in_dat,
                               out_dat,
                               cut_height = False, 
                               min_size = 10 ):
    """
    Receives two rds files (saved R objects),  one containing a flashClust
    dendrogram, the other containing the dataframe from which it was derived. 
    If a cut height is passed, then uses WGCNA cutreeStatic to cluster data and
    saves data for samples belonging to first cluster, before re-clustering them
    and plotting the resulting dendrogram.
    If no cut height is passed then the input data are saved and no plot is
    created.
    """

    in_tree = os.path.abspath( in_tree )
    in_dat = os.path.abspath( in_dat )

    # load the R data
    R( '''datExpr <- readRDS( "%(in_dat)s" )''' % locals() )
    R( '''sampleTree <- readRDS( "%(in_tree)s" )''' % locals() )

    if cut_height:
        # create filenames for plots of the cut and re-clustered dendrograms
        out_cut_height = P.snip( in_dat, "_dat.rds" ) + "_cutHeight.png"
        out_tree = P.snip( in_dat, "_dat.rds" ) + "_outliersRemoved.png"

        # create plot showing where tree is to be cut
        R( '''png( "%(out_cut_height)s" )''' % locals() )
        R( '''plot( sampleTree, 
                    main = "Sample clustering prior to outlier removal",
                    cex.lab = 1.5,
                    cex.axis = 1.5,
                    cex.main = 2 )''' )
        R( '''abline( h = %(cut_height)s, col = "red" )''' % locals() )
        R( '''dev.off()''' )
        
        # cluster samples based on supplied cut height
        R( '''clust <- cutreeStatic( sampleTree, 
                                     cutHeight = %(cut_height)s,
                                     minSize = %(min_size)s )''' % locals() )
        # identify samples forming part of the principal cluster
        R( '''keepSamples <- (clust == 1)''' )
        # retrieve samples belonging to the principal cluster
        R( '''datExpr <- datExpr[keepSamples,]''' )
        # re-cluster the samples
        R( '''sampleTree1 <- flashClust( dist( datExpr ), method = "average" )''' )
        # replot the resulting dendrogram
        R( '''png( "%(out_tree)s" )''' % locals() )
        R( '''plot( sampleTree1, 
                    main = "Sample clustering following outlier removal",
                    cex.lab = 1.5,
                    cex.axis = 1.5,
                    cex.main = 2 )''' )
        R( '''dev.off()''' )

    # write the expression dataframe as an rds file
    R( '''saveRDS( datExpr, file="%(out_dat)s" )''' % locals() )


def wgcnaPlotSoftThreshold( in_dat, 
                            powers, 
                            threshold, 
                            out_scaleInd, 
                            out_meanCon ):
    """
    Recieves a dataframe as an .rds file, plus a list of potential thresholds.
    Outputs the two suggested plots for picking soft-threshold used to 
    approximate a scale free network. 
    """
    # load dataframe
    R( '''datExpr <- readRDS( "%(in_dat)s" )''' % locals() )

    # create R vector from list of powers
    powers = robjects.IntVector( powers )
    R( '''powers <- %s''' % powers.r_repr() )

    # run WGCNA pickSoftThreshold
    R( '''sft <- pickSoftThreshold( datExpr,
                                    powerVector = powers, 
                                    verbose = 5 )''' % locals() )

    # define plot functions
    R( ''' plot.top.fit <- function( sft, cut.off, powers, outfile="./top_fit.png" ){
             sizeGrWindow( 9, 5 )
             cex1 = 0.9
             png( outfile )
             plot( sft$fitIndices[,1], 
                   -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
                   xlab="Soft threshold (power)",
                   ylab="Scale free topology model fit, signed R^2", 
                   type="n", 
                   main= paste("Scale independence") )
             text( sft$fitIndices[,1], 
                   -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
                   labels=powers,
                   cex=0.9,
                   col="red" )
             abline( h=cut.off, col="red")
             dev.off()
            }''' )
    R( '''plot.connectivity <- function( sft, powers, outfile="./mean_con.png" ){
            sizeGrWindow( 9, 5 )
            cex1 = 0.9
            png( outfile )
            plot(sft$fitIndices[,1], 
                 sft$fitIndices[,5],
                 xlab="Soft threshold (power)", 
                 ylab="Mean Connectivity", 
                 type="n", 
                 main=paste("Mean connectivity"))
            text( sft$fitIndices[,1], 
                  sft$fitIndices[,5], 
                  labels=powers, 
                  cex=cex1, 
                  col="red")
            dev.off()
            }''' )

    # plot figures
    R( '''plot.top.fit( sft, %(threshold)s, powers, outfile="%(out_scaleInd)s" )''' % locals() )
    R( '''plot.connectivity( sft, powers, outfile="%(out_meanCon)s" )''' % locals() )


def wgcnaCalcCoexpressionNetwork( in_dat, out_adj, out_tom, out_scaledCon, power ):
    """
    Receives a dataframe as an .rds file, using power provided, calculates an 
    adjacency matrix, a topological overlap matrix, and network concepts.
    """
    # load dataframe from infile
    R( '''datExpr <- readRDS( "%(in_dat)s" )''' % locals() )
    # calculate adjacency matrix
    R( '''adj <- adjacency( datExpr, power = %(power)s )''' % locals() )
    R( '''saveRDS( adj, file="%(out_adj)s" )''' % locals() )
    # calculate TOM
    R( '''TOM <- TOMsimilarity( adj )''' )
    R( '''saveRDS( TOM, file="%(out_tom)s" )''' % locals() )

    # calculate Fundamental Network Concepts... is taking days to complete
    #R( '''NC <- networkConcepts( datExpr, power = %(power)s )''' % locals() )
    #R( '''saveRDS( NC, file="%(out_nc)s" )''' % locals() )

    # in the interim... calcuated and plot the scaled connectivity distribution. 
    # see connectivity calculation in fundamentalNetworkConcepts()
    R( '''conn <- apply( adj, 2, sum )''' )
    R( '''scaledConn <- conn/max(conn, na.rm=T)''' )
    R( '''png( "%(out_scaledCon)s" )''' % locals() )
    R( '''hist( scaledConn, breaks=1000 )''' )
    R( '''dev.off()''' )

#def wgcnaClusterTOM( inf, outf, params ):
def wgcnaClusterTOM( in_dat, out_geneTree, out_dendro ):
    """
    Receives a TOM as an .rds file, converts it to dissTOM, uses flashClust to 
    cluster and saves cluster. Plots resulting dendrogram and visualized TOM.
    """
    # setting up to work with P.submit doesn't work due to plotting. 
    # R( '''suppressMessages(library(WGCNA))''' )
    # R( '''allowWGCNAThreads(nThreads=8)''' )
    # in_dat, out_geneTree, out_dendro = params

    # load TOM from infile, convert to dissTOM
    R( ''' TOM <- readRDS( "%(in_dat)s" )''' % locals() )
    R( '''dissTOM <- 1-TOM''' )

    # cluster TOM and save
    R( '''geneTree <- flashClust( as.dist( dissTOM ), method="average" )''' )
    R( ''' saveRDS( geneTree, file="%(out_geneTree)s" )''' % locals() )
    E.info( "Gene clustering successfully completed, saved to: %s" % out_geneTree  )

    # plot the resulting geneTree as dendrogram
    R( '''sizeGrWindow(12, 9)''' )
    R( '''png( "%(out_dendro)s" )''' % locals() )
    R( '''plot( geneTree, 
                xlab = "", 
                sub = "", 
                main = "Gene clustering on TOM-based dissimilarity", 
                labels = FALSE, 
                hang = 0.04 )''' )
    R('''dev.off()''' )
    E.info( "Plotted the resulting cluster as dendrogram: %s" % out_dendro )

    # visualise the TOM using TOMplot... is causing a segfault
    # E.info( "Beginning to plot TOM visualization" )
    # R( '''png( "%(out_TOMplot)s" )''' % locals() )
    # R( '''TOMplot( dissTOM, 
    #                geneTree )''' ) 
    # R( '''dev.off()''' )
    # E.info( "TOM visualization successfully completed, saved to: %s" % out_TOMplot )


def wgcnaCutreeDynamic( in_tree, 
                        in_TOM, 
                        method,
                        out_modules, 
                        out_tree_png,
                        deep_split, 
                        min_cluster_size, 
                        **kwargs ):
    """
    Receives a geneTree (hclust object), TOM, and original dataframe.
    Runs WGCNA cutreeDynamic. Outputs a plot of modules under dendrogram 
    and a character vector in which colours represent module assignment for 
    each gene (grey = unassigned ).
    Requires method (tree/hybrid), deep_split ([0-3]/Bool), and min_module_size
    to be passed explicitly, other options may be passed to cutreeDynamic 
    via kwargs.
    """
    # sort kwargs
    if kwargs:
        kw_list = createRkwargs( kwargs )
    else:
        kw_list = ""

    # handle boolean... there seems to be no easy way to assign boolean robjects
    if isinstance( deep_split, (bool,) ):
        if deep_split:
            deep_split = "TRUE"
        else:
            deep_split = "FALSE"

    # load the geneTree()
    R( '''geneTree <- readRDS( "%(in_tree)s" )''' % locals() )
    R( '''TOM <- readRDS( "%(in_TOM)s" )''' % locals() )
    R( '''dissTOM <- 1 - TOM''' )

    # run cutreeDynamic
    R( ''' dynamicMods <- cutreeDynamic( dendro = geneTree, 
                                         distM = dissTOM, 
                                         method = "%(method)s",
                                         deepSplit = %(deep_split)s, 
                                         minClusterSize = %(min_cluster_size)s 
                                         %(kw_list)s )''' % locals() )
    
    # convert numeric module labels to colours. 
    R( '''dynamicColours <- labels2colors( dynamicMods )''' )
    R( '''saveRDS( dynamicColours, file = "%(out_modules)s" )''' % locals() )

    # plot modules below gene dendrogram
    R( '''png( "%(out_tree_png)s" )''' % locals() )
    R( '''plotDendroAndColors( geneTree, 
                               dynamicColours, 
                               groupLabels = "DynamicTreeCut",
                               dendroLabels = FALSE, 
                               hang = 0.03,
                               addGuide = TRUE, 
                               guideHang = 0.05, 
                               main = "Gene dendrogram and module colours" )''' )                                    
    R( '''dev.off()''' )


def wgcnaModuleEigengenes( in_dat, 
                           in_modules, 
                           out_eigen, 
                           out_eigenplot, 
                           merge_dist = False,
                           **kwargs ):
    """
    Receives a data frame of expression values and a character vector assigning
    genes to colour modules. Runs WGCNA moduleEigenGenes to return the first PC
    for each module, saves the resulting list and plots eigengene dendrogram.
    Additional parameters can be passed to moduleEigenGenes via **kwargs.
    """
    # assert MEList$allOK == True
    
    if kwargs:
        kw_list = createRkwargs( kwargs )
    else:
        kw_list = ""

    # load expression dataframe and module assignment vector
    R( '''datExpr <- readRDS( "%(in_dat)s" )''' % locals() )
    R( '''dynamicColours <- readRDS( "%(in_modules)s" )''' % locals() )
    
    # determine module eigengenes
    R( '''MEList <- moduleEigengenes( datExpr, 
                    colors = dynamicColours
                    %(kw_list)s )''' % locals() )
    R( '''saveRDS( MEList, file = "%(out_eigen)s" )''' % locals() )

    # extract and cluster the eigengenes
    R( '''MEs <- MEList$eigengenes''' )
    R( '''MEDiss <- 1-cor(MEs)''' )
    R( '''METree <- flashClust( as.dist( MEDiss ), method = "average" )''' )

    # plot
    R( '''sizeGrWindow(7, 6)''' )
    R( '''png( "%(out_eigenplot)s" )''' % locals() )
    R( ''' plot( METree, 
                 main = "Clustering module eigengenes", 
                 xlab = "", 
                 sub = "" )''' )
    if merge_dist:
        R( '''abline( h=%(merge_dist)s, col="red" )''' % locals() )
    R( '''dev.off()''' )


def wgcnaMergeCloseModules( in_dat, 
                            in_modules, 
                            out_modules, 
                            out_eigen, 
                            merge_dist, 
                            **kwargs ):
    """
    Receives a dataframe of expression values, a character vector assigning genes
    to colour modules, and a merge_distance specifying the value 
    (1 - pearson correlation) at which to merge modules, based on similarity of 
    their eigengenes.
    Outputs a new character vector assigning genes to merged modules and a list 
    containing information about the merged module eigengenes.
    Additional parameters can be passed to mergeCloseModules via **kwargs.
    NB. There is no way to pass additional parameters for the function to 
    calculate merged module eigenvectors. 
    """
    if kwargs:
        kw_list = createRkwargs( kwargs )
    else:
        kw_list = ""

    # load expression dataframe and module assignment vector
    R( '''datExpr <- readRDS( "%(in_dat)s" )''' % locals() )
    R( '''dynamicColours <- readRDS( "%(in_modules)s" )''' % locals() )

    # merge close modules and save new module assignment vector
    R( '''merge <- mergeCloseModules( datExpr, 
                                      dynamicColours, 
                                      cutHeight = %(merge_dist)s, 
                                      verbose = 5
                                      %(kw_list)s )''' % locals() )
    assert R( '''toString(merge$allOK)''' )[0] == "TRUE"
    R( '''mergedColours <- merge$colors''' )
    R( '''saveRDS( mergedColours, file="%(out_modules)s" )''' % locals() )

    # calculate module eigengenes for the new vectors and save
    R( '''MEList <- moduleEigengenes( datExpr, colors = mergedColours )''' )
    assert R( '''toString(MEList$allOK)''' )[0] == "TRUE"
    R( '''saveRDS( MEList, file = "%(out_eigen)s" )''' % locals() )


def wgcnaPlotDendroAndColors( in_geneTree, in_modules, out_png, **kwargs ):
    """
    Receives a geneTree (hclust object) and character vector assigning genes to 
    colour modules. Runs plotDendroAndColors. Additional options can be passsed 
    via **kwargs.
    """
    if kwargs:
        kw_list = createRkwargs( kwargs )
    else:
        kw_list = ""

    # load geneTree and module assignment vector
    R( '''geneTree <- readRDS( "%(in_geneTree)s" )''' % locals() )
    if isinstance( in_modules, basestring ):
        R( '''mergedColours <- readRDS( "%(in_modules)s" )''' % locals() )
    elif isinstance( in_modules, (list, tuple) ):
        merged_modules = []
        for i in range( len( in_modules ) ):
            current_module = in_modules[i]
            merged_module = "mergedColours"  + str(i)
            merged_modules.append( merged_module )
            R( ''' %(merged_module)s <- readRDS( "%(current_module)s" )  ''' % locals() )
        merged_modules =  "cbind(" + ",".join( merged_modules ) + ")"
        R( '''mergedColours <- %(merged_modules)s ''' % locals() )
    else:
        raise ValueError( "Unknown type for in_modules: %s" % in_modules )         

    # plot dendrogram
    R( '''sizeGrWindow(12, 6)''' )
    R( '''png( filename = "%(out_png)s", width = 960, height = 480 )''' % locals() )
    R( '''plotDendroAndColors( geneTree, 
                               mergedColours, 
                               dendroLabels = FALSE, 
                               hang= 0.03,
                               addGuide = TRUE,
                               guideHang = 0.05
                               %(kw_list)s )''' % locals() )
    R( '''dev.off()''' )


def wgcnaExtractModuleAssignment( in_modules, in_dat, outfile ):
    # load the datafiles
    R( '''geneModules <- readRDS( "%(in_modules)s" )''' % locals() )
    R( '''datExpr <- read.delim( "%(in_dat)s",
                                  header = TRUE, 
                                  sep = ",", 
                                  row.names=1, 
                                  stringsAsFactors=FALSE ) ''' % locals() )

    # set geneIDs and tables
    R( '''gene.ids <- names( datExpr )''' )
    R( '''module.table <- table( geneModules )''' )
    R( '''gene.table <- gene.ids''' )
    
    # iteratively add module assignment
    R( '''for( mod in names( module.table ) ){
            gene.table <- cbind( gene.table, geneModules == mod )
            }''' )
    
    # assign headers
    R( '''colnames( gene.table ) <- append( "gene_id", names( module.table ) ) ''' )

    # write outfile
    R( '''gene.df <- as.data.frame( gene.table )''' )
    R( '''write.table( gene.df, 
                       file = "%(outfile)s", 
                       quote=FALSE,
                       sep="\t",
                       row.names=FALSE )''' % locals() )


def wgcnaPlotModuleSpecificHeatmaps( in_dat, in_modules, in_tree, out_dir ):
    # load required library and set plot options
    R( '''library( gplots )''' )
    R( '''cols <- colorRampPalette(c("white", "darkred" ), space="rgb")(100)''' )

    # load sampleTree, module assignment vector, expression dataframe
    R( '''sampleTree <- readRDS( "%(in_tree)s" ) ''' % locals() )
    R( '''sampleDendro <- as.dendrogram( sampleTree )''' )
    R( '''geneModules <- readRDS( "%(in_modules)s" )''' % locals() )
    R( '''datExpr <- read.delim( "%(in_dat)s",
                                  header = TRUE, 
                                  sep = ",", 
                                  row.names=1, 
                                  stringsAsFactors=FALSE ) ''' % locals() )

    # create a table of module assignment
    R( '''module.table <- table( geneModules )''' )

    # loop through modules, creating a heatmap for each one
    R( '''
    for(i in 1:(length(module.table))){
        # create outfile name
        x <- names( module.table )[i]
        outfile <- paste("%(out_dir)s/clusteredHeatmap2Module_", x, ".png", sep="" )
  
        # create boolean that identifies genes beloning to a module. 
        module.genes <- ( geneModules == x )
        # extract genes belonging to module from the  original expression dataset
        module.geneset <- datExpr[, module.genes]

        ## transform geneset
        ## samples are rows and genes are columns
        # iterate over columns
        # zero adjust the gene expression expression data, based on the minimum value
        #module.geneset <- sweep( module.geneset, 2, apply( module.geneset, 2, function(n) min(n)*-1), "+")
        module.geneset <- sweep( module.geneset, 2, apply( module.geneset, 2, function(n) sqrt(min(n)^2)), "+")
        # standardize the lncRNA expression by dividing by the maximum value
        module.geneset <- sweep( module.geneset, 2, apply( module.geneset, 2, max), "/")
        # make sure the rows for the module geneset are in the same order as the rows for the expression data
        module.geneset <- module.geneset[rownames(datExpr),]
        # create a matrix from the module geneset... 
        module.geneset.matrix <- as.matrix( module.geneset )
  
        # plot
        outfile.title <- paste( "Module", x, "heatmap (", length(colnames(module.geneset)), "Genes)" )
        png( outfile, width=1440, height=960 )
        
        heatmap.2( module.geneset.matrix, 
                   Rowv = sampleDendro,
                   Colv = FALSE,
                   col=cols,
                   dendrogram="row",
                   trace = "none", 
                   #na.rm = TRUE,
                   margins=c(8,8),
                   cexRow = 1.5,
                   breaks = 101,
                   labCol = FALSE, 
                   key = FALSE, 
                   main = outfile.title )
  
 
        dev.off()
     }''' % locals() )


def wgcnaPlotModuleEigengenes( in_eig, in_dat, out_dir, to_remove ):
    """
    """
    # Retrieve module eigengenes
    R( '''module.eigen <- readRDS( "%(in_eig)s" )''' % locals() )
    R( '''MEs <- module.eigen$eigengenes''' )

    # Generate list of sample names for eigengenes
    R( '''datExpr <- read.delim( "%(in_dat)s",
                                  header = TRUE, 
                                  sep = ",", 
                                  row.names=1, 
                                  stringsAsFactors=FALSE ) ''' % locals() )
    R( ''' eigen.names <- rownames(datExpr)''' )
    for sample in to_remove:
        R( '''eigen.names <- eigen.names[eigen.names != "%s"]''' % sample )

    # plot eigengenes using barplot
    R( '''
    for( i in 1:length(names(MEs)) ){
        mod.name <- names(MEs)[i]
        outfile <- paste("%(out_dir)s/moduleEigenGene_", mod.name, ".png", sep="" )
        main <- paste( "Module", mod.name, "eigengene expression\n(", round(module.eigen$varExplained[[i]], 2), ")" )
        png( outfile )
        barplot(MEs[[i]], names.arg=eigen.names, las=2, main=main)
        dev.off()
        } ''' % locals() )


def extractEigengenes( in_csv, in_eigen, to_remove, outfile ):
    """
    Receives i) an expression dataframe in which columns are genes and rows are
    samples ii) an RDS file containing output from extractmoduleEigengenes()
    iii) a list of samples removed before calculating eigengenes
    Creates a tsv file with eigengenes as rows and samples as columns
    com.load doesn't work...
    """

    R('''rm(list=ls())''')
    # read in dataframe from csv
    in_csv = os.path.abspath( in_csv )
    R( '''datExpr <- read.delim( "%(in_csv)s",
                                 header=TRUE, 
                                 sep=",", 
                                 row.names=1,
                                 stringsAsFactors=FALSE )''' % locals() )
    # fetch the sample_ids from the expression dataframe
    R( '''sample.ids <- row.names( datExpr )''' )

    # com.load_data doesn't work because of issues with numpy2ri.explicit
    sample_ids = com.load_data( "sample.ids" )    

    # exclude ids for removed samples
    eigen_ids = []
    for i in sample_ids:
        if i in to_remove:
            continue
        else:
            eigen_ids.append( i )

    # read in output from moduleEigengenes()
    in_eigen = os.path.abspath( in_eigen )
    R( '''MEList <- readRDS( "%(in_eigen)s" )''' % locals() )
    R( '''eigen.genes <- MEList$eigengenes''' )

    # Fails to import...
    eigen_genes = com.load_data( "eigen.genes" )
    eigen_genes = eigen_genes.transpose()

    # check number of elements in eigenvector equals number of sample_ids
    assert len( eigen_ids ) == len( eigen_genes.columns ), "# sample_ids doesn't match # columns"

    # write outfile
    eigen_genes.to_csv( outfile, 
                        sep = "\t", 
                        na_rep = "NA", 
                        header = eigen_ids, 
                        index_label="module_id" )


def extractEigengenes_R( in_csv, in_eigen, to_remove, outfile ):
    """
    Receives i) an expression dataframe in which columns are genes and rows are
    samples ii) an RDS file containing output from extractmoduleEigengenes()
    iii) a list of samples removed before calculating eigengenes
    Creates a tsv file with eigengenes as rows and samples as columns
    """

    R('''rm(list=ls())''')
    # read in dataframe from csv
    in_csv = os.path.abspath( in_csv )
    R( '''datExpr <- read.delim( "%(in_csv)s",
                                 header=TRUE, 
                                 sep=",", 
                                 row.names=1,
                                 stringsAsFactors=FALSE )''' % locals() )

    # fetch the sample_ids from the expression dataframe
    R( '''sample.ids <- row.names( datExpr )''' )
    # assign the to_remove list to r environment
    R.assign("to_remove", robjects.StrVector(to_remove))
    # remove sample ids for absent samples
    R("""sample.ids <- sample.ids[!(sample.ids %in% to_remove)]""")

    # read in output from moduleEigengenes()
    in_eigen = os.path.abspath( in_eigen )
    R( '''MEList <- readRDS( "%(in_eigen)s" )''' % locals() )
    R( '''eigen.genes <- t(MEList$eigengenes)''' )
    R( '''colnames(eigen.genes) <- sample.ids''' )

    # write to tempfile
    tmpf = P.getTempFilename("/ifs/scratch")
    R('''write.table(eigen.genes, file="%s", append=FALSE, quote=FALSE, sep="\t")''' % tmpf)

    # fix table column names in pandas... 
    df = pd.read_table(tmpf, index_col=0)
    df.to_csv(outfile, sep="\t", index_label="module_id", na_rep="NA")

    os.unlink(tmpf)

def calcEigenRINCorrelation( eigen_df, rin_vals, out_dir, outfile ):
    
    eigen_rin_cor_dict = collections.defaultdict( list )
    p_thresh = 0.05/len( eigen_df.index )

    for module in eigen_df.index:
        rin_list = []
        eig_list = []
        module_eigen_vals = eigen_df.loc[ module ] 

        # create lists of eign values and rin values
        for sample in module_eigen_vals.index:
            if sample in rin_vals.index:
                rin_list.append( rin_vals.loc[ sample ] )
                eig_list.append( module_eigen_vals.loc[ sample ] )
            else:
                raise Exception( "The sample %s does not have an associated RIN value" % sample )

        eig, rin = filterMasked( eig_list, rin_list, mask_nan = True )
        result = stats.pearsonr( eig, rin )
        eigen_rin_cor_dict[ module ].append( result[0] ) # rval
        eigen_rin_cor_dict[ module ].append( result[1] ) # pval
        eigen_rin_cor_dict[ module ].append( p_thresh ) # p_adj
        if result[1] > p_thresh:
            eigen_rin_cor_dict[ module ].append( "PASS" )
        else:
            eigen_rin_cor_dict[ module ].append( "FAIL" )

            # plot failed correlations
            failed_plot = out_dir + "/rinCor_" + module + ".png" 
            R.assign('eig', robjects.vectors.FloatVector( eig ) )
            R.assign('rin', robjects.vectors.FloatVector( rin ) )

            R( '''png( "%(failed_plot)s" )''' % locals() )
            R( '''plot( eig, 
                        rin, 
                        main="%(module)s", 
                        pch=16, 
                        cex.lab=1.5, 
                        col="darkred", 
                        xlab="expression", 
                        ylab="RIN #" )''' % locals() )
            R( '''dev.off()''' )

    # write outfile 
    outf = IOTools.openFile( outfile, "w" )
    outf.write( "module_id\tpearson_r\tp_val\tp_threshold\tresult\n" )
    for key, value in eigen_rin_cor_dict.iteritems():
        outf.write( key + "\t" + "\t".join( [ str(x) for x in value ] ) + "\n" )
    outf.close()



def plotZscoreHeatmaps( infiles, outfiles ):
    """
    Plot heatmaps of zscores.
    """
    in_modules = P.snip(os.path.basename( infiles[0] ), ".tsv.gz" )
    in_expression = P.snip(os.path.basename( infiles[1] ), ".load" )
    out_stub = os.path.dirname( infiles[0] )
    out_stub = re.sub( "_assignment", "_heatmaps_zscore", out_stub ) 

    # get a list of the modules in no particular order... 
    statement = ("SELECT * FROM %(in_modules)s LIMIT 1" % locals())
    df = PU.fetch_DataFrame( statement )
    df = df.set_index( "gene_id" )
    modules = list( df.columns.values )
    E.info( "Retieved module IDs for %s:\n%s" % ( in_modules, " ".join( modules ) ) )

    # load R library
    R( '''suppressMessages(library(gplots))''' )

    # iterate through modules and plot a heatmap of zscores for genes in each module
    for module in modules:
        E.info( "Plotting module %s" % module )
        # specify outfile
        outfile1 = "clusteredHeatmapModule_" + module + ".png"
        outfile2 = "clusteredHeatmapModule_" + module + ".eps"
        outfile1 = os.path.join( out_stub, outfile1 )
        outfile2 = os.path.join( out_stub, outfile2 )

        # fetch expression data for module
        statement = ( "SELECT"
                      " a.gene_id,"
                      " a.pro,"
                      " a.pre,"
                      " a.immature," 
                      " a.mature,"
                      " a.follicular," 
                      " a.marginal,"
                      " a.b1a,"
                      " a.germinal"
                      " FROM %(in_expression)s AS a"
                      " INNER JOIN %(in_modules)s AS b"
                      " ON a.gene_id = b.gene_id"
                      " WHERE b.%(module)s = 1" % locals() )
        df = PU.fetch_DataFrame( statement )
        df = df.set_index( "gene_id" )
        E.info( "Retrieved expression data for module %s" % module )
        tmpf = P.getTempFilename(".")
        df.to_csv( tmpf, sep="\t", na_rep="NaN" )
        E.info( "Written module expression data to %s" % tmpf )

        # set header
        length = len( df.index )
        header = module + "\n(" + str(length) + " genes)" 

        # pass expression data to R... isn't working.. writing to flatfile instead
        # df = com.convert_to_r_dataframe( df )
        # R.assign( "df", df )

        # transpose matrix, calculate zscores (rows/genes), and plot heatmap
        E.info( "Reading table into R" )
        R( '''df <- read.table( "%s", header=T, row.names=1, stringsAsFactors=F)''' % tmpf )
        R( '''df <- t(df)''' )
        E.info( "Calculating zscores" )
        R( '''df <- scale(df, center=TRUE, scale=TRUE)''' )
        # plot
        E.info( "Writing to outfile %s" %  outfile1 )
        R( '''png( "%s" )''' % outfile1 )
        R( '''heatmap.2(df, Rowv=F, Colv=T, dendrogram='column', scale='none',
                        col=bluered(256), trace="none", labCol="none", density.info="none",
                        main='%s', col.main='%s' )''' % ( header, module ) )
        R( '''dev.off()''' )
        E.info( "Successfully plotted module %s" % module )
        R( '''setEPS()''' )
        R( '''postscript("%s")''' % outfile2 )
        R( '''heatmap.2(df, Rowv=F, Colv=T, dendrogram='column', scale='none',
                         col=bluered(256), trace="none", labCol="none", density.info="none" )''' )
        R( '''dev.off()''' )
        E.info( "clearing global namespace" )
        R( '''rm(list=ls())''' )
        os.unlink(tmpf)

def wgcnaPlotModuleDendrograms( df_expr, pos_cor_lnc, neg_cor_lnc, outfile, module ):
    """
    Takes a dataframe of lncRNA and gene expression values for a particular wgcna
    module, plus lists of which lncRNAs are positively/negatively correlated with
    the module eigengene.
    Plots a dendrogram genes in module
    """
    outfile_orig = P.snip( outfile, ".tsv" ) + "_original.png"
    outfile_orig_eps = P.snip( outfile_orig, ".png" ) + ".eps"
    outfile_abs = P.snip( outfile, ".tsv" ) + "_absolute.png"
    outfile_abs_eps = P.snip( outfile_abs, ".png" ) + ".eps"
    # Because I'm paranoid
    R( '''rm(list=ls())''' )

    # hack... it's necessary to alter the width of the outfile based on the
    # number of genes in the module.
    n_genes = float( len( IOTools.openFile( outfile ).readlines() ) - 1 )
    # assuming that 31 genes fits into standard width of 480
    out_width = str(15.5 * n_genes)

    # convert to r objects
    #    R.assign( "df.exp", com.convert_to_r_dataframe( df_expr ) )
    #    R('''df.exp <- as.matrix(df.exp)''')
    #    R('''df.exp <- as.numeric(df.exp)''')
    #    R('''print(class(df.exp))''')
    #    R('''print(typeof(df.exp))''')
    #    R('''print(names(df.exp))''')
    #    R('''print(row.names(df.exp))''')
    #    R('''print(df.exp[,1:5])''')
    #    exit
    # convert to r objects
    R.assign( "pos.cor", robjects.vectors.StrVector( pos_cor_lnc ) )
    R.assign( "neg.cor", robjects.vectors.StrVector( neg_cor_lnc ) )
    
    # having to read in data frame because of issuse with assigning... 
    E.info( "READING IN TABLE %s" % outfile )
    R('''df.exp <- read.table('%s', header=TRUE, row.names=1, sep="\t", stringsAsFactors=F, na.strings="NULL")''' % outfile )

    # remove genes with NA values
    R('''df.exp <- df.exp[complete.cases(df.exp),]''')

    # write function for labelling colours in R
    # see http://rstudio-pubs-static.s3.amazonaws.com/1876_df0bf890dd54461f98719b461d987c3d.html
    R("""
      colLab <- function(n){
        if(is.leaf(n)){
          a <- attributes(n)
        if( a$label %in% pos.cor ){labCol = "red"}
        else if( a$label %in% neg.cor ){labCol = "blue"}
        else {labCol = "black"}
        attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
        }
        n
      }
    """)

    # calculate pearson correlations using the expression dataframe
    R('''m.cor <- cor(t(df.exp), method="pearson")''')
    # convert to absolute pearson
    R('''m.cor.abs <- abs(m.cor)''')
    # convert to distance matrix
    R('''m.dist <- as.dist(1 - m.cor)''')
    R('''m.dist.abs <- as.dist(1 - m.cor.abs)''')
    # cluster
    R('''clust <- hclust(m.dist)''')
    R('''clust.abs <- hclust(m.dist.abs)''')
    # convert to dendrogram
    R('''dend <- as.dendrogram(clust)''')
    R('''dend.abs <- as.dendrogram(clust.abs)''')
    # apply function
    R('''dend.1 <- dendrapply(dend, colLab)''')
    R('''dend.1.abs <- dendrapply(dend.abs, colLab)''')
 
    # plot
    R('''png('%s', width=%s)''' % ( outfile_orig, out_width ) )
    R('''plot(dend.1, cex=0.6, main='%s')''' % module)
    R('''dev.off()''')
    R('''setEPS()''')
    R('''postscript('%s', width=%s)''' % ( outfile_orig_eps, out_width ) )
    R('''plot(dend.1, cex=0.6)''')
    R('''dev.off()''')

    R('''png('%s', width=%s)''' % ( outfile_abs, out_width ) )
    R('''plot(dend.1.abs, cex=0.6, main='%s')''' % module)
    R('''dev.off()''')
    R('''setEPS()''')
    R('''postscript('%s', width=%s)''' % ( outfile_abs_eps, out_width ) )
    R('''plot(dend.1.abs, cex=0.6)''')
    R('''dev.off()''')


def test_nn_module_membership( df_lncRNA, df_refcoding ):
    """
    Receives two pandas dataframes: 
    df_lncRNA contains lnc_id, nn_gene_id, chromatin_status, genomic_location,
    lnc module affiliation (As True/False)
    df_refcoding contains gene_id and gene module affiliation (As 1/0)

    

    NB. latter df has prefix ME on all module names
    """
    # get a list of module ids from the refcoding dataframe... 
    # (colnames are unicode)
    module_ids = [ str(x) for x in df_refcoding.columns ]
    
    # generate a dictionary to store results
    module_dict = collections.defaultdict( dict )
    
    for module in module_ids:
        # modules have ME prefix in df_lncRNA
        module_p = "ME" + module
    
        # set up dictionary to store list of nn_genes in the same module vs
        # different module compared to lncRNA
        module_dict[module]["Same"] = []
        module_dict[module]["Different"] = []
        module_dict[module]["Excluded"] = []
    
        # get list of nn_genes for lncRNAs affiliated with module
        nn_list = [ str(x) for x in list( 
            df_lncRNA["nearest_protein_coding"][df_lncRNA[module_p]=="True"] ) ]
        
        for nn_gene in nn_list:
            try:
                same_module = df_refcoding.loc[nn_gene][module] 
                if same_module:
                    module_dict[module]["Same"].append( nn_gene )
                else:
                    module_dict[module]["Different"].append( nn_gene )
            except KeyError:
                module_dict[module]["Excluded"].append( nn_gene )
                
    # return module_dict
    
    # calculate the probability of a gene being in a particular module
    geneset_size = sum( [sum(df_refcoding[x]) for x in module_ids] )
    module_probs = {}
    for module in module_ids:
        module_size = sum(df_refcoding[module])
        module_probs[module] = float(module_size)/float(geneset_size)

    # create a dictionary containing:
    # i) number of times a nearest neighbour gene appears in module (successes)
    # ii) number of lncRNAs correlated with a module (trials)
    # iii) probability of a gene being in a module (probs)
    prop_dict = {}
    # create lists for R stats
    s = []
    t = []
    p = []
    for module, values in module_dict.iteritems():
        successes = float( len(values["Same"]) )
        s.append( successes )
        trails = float( len(values["Same"]) + len(values["Different"]) )
        t.append( trials )
        probs = module_probs[module]
        p.append( probs )
        prop_dict[module] = ( successes, trials, probabilities )

    # Run test... doing this in R because of subsetting
    # convert data to R objects
    R.assign("suc", robjects.FloatVector( s ) )
    R.assign("tr", robjects.FloatVector( t ) )
    R.assign("prob", robjects.FloatVector( p ) )
    R("""suc.adj = suc[tr>0]""")
    R("""tr.adj = tr[tr>0]""")
    R("""prob.adj = prob[tr>0]""")
    R("""result = prop.test(x=suc.adj, n=tr.adj, p=prob.adj)""")
    R("""cat(result, file=%(outfile)s)""")

    return module_dict, prop_dict

#################################################################################
#################################################################################
#################################################################################
## section: eRNA analysis 
#################################################################################
#
#
#
#################################################################################

def plotMatrixAndControlPeakShape( infiles, outfile, params=None ):
    """
    Receive matrix and control file output from bam2peakshape, plot heatmaps
    using heatmap2, merge heatmaps using ImageMagick
    """
    in_matrix, in_control = infiles
    # write heatmap2 function for plotting matrix and control heatmaps
    R( '''library( gplots )''' )
    R( '''library( RColorBrewer )''' )

    R('''
      plot.heatMaps <- function( df.matrix, df.control, outfile1, outfile2 ){
        cols= brewer.pal(9, "Blues")
        png( outfile1, width=480, height=960 )
        heatmap.2( as.matrix(df.matrix), 
                   trace="none", 
                   dendrogram="none", 
                   Rowv=FALSE, 
                   Colv=FALSE, 
                   labRow="", 
                   labCol="",
                   col=cols, 
                   breaks=seq(0, 100, 11) )
        dev.off()
        cols= brewer.pal(9, "Greys")
        png( outfile2, width=480, height=960 )
        heatmap.2( as.matrix(df.control), 
                   trace="none", 
                   dendrogram="none", 
                   Rowv=FALSE, 
                   Colv=FALSE, 
                   labRow="", 
                   labCol="",
                   col=cols, 
                   breaks=seq(0, 100, 11) )
        dev.off()
        }
      ''')

    # generate R dataframes for matrix and control 
    df_matrix = pd.read_table( in_matrix, 
                               compression = "gzip", 
                               index_col = None,
                               header = 0 )
    if params:
        lnc_ids = pickle.load( open( params[0] ) )
        df_matrix = df_matrix[df_matrix['name'].isin( lnc_ids )]
    df_matrix = df_matrix.set_index('name')
    df_matrix = pd.rpy.common.convert_to_r_dataframe( df_matrix )
    df_control = pd.read_table( in_control, 
                                compression = "gzip", 
                                index_col = None,
                                header = 0 )
    if params:
        lnc_ids = pickle.load( open( params[0] ) )
        df_control = df_control[df_control['name'].isin( lnc_ids)]
    df_control = df_control.set_index('name')
    df_control = pd.rpy.common.convert_to_r_dataframe( df_control )

    # run R plot function from python
    outf1 = P.getTempFilename(".")
    outf2 = P.getTempFilename(".")
    plot_heatMaps = R['plot.heatMaps']  
    plot_heatMaps( df_matrix, df_control, outf1, outf2 )

    # merge the two heatmaps on the command line
    statement = ( "montage"
                  " -geometry '1x1+0+0<'"
                  " %(outf1)s"
                  " %(outf2)s"
                  " %(outfile)s" )
    P.run()


def plotMatrixPeakShape( infiles, outfile, params ):
    """
    A refinement on previous function, receives two peakshape files, 
    outputs a merged plot of coverage ranked by me1/me3 ratio
    """
    # specify infiles
    in_me3, in_me1 = infiles
    lnc_ids, ratios = [ pickle.load( open( x ) ) for x in params ]
    
    # subset infiles and convert to R objects
    ratios = robjects.vectors.StrVector( ratios )

    df_me3 = pd.read_table( in_me3, 
                            compression = "gzip", 
                            index_col = 0, 
                            header = 0 )
    df_me3 = df_me3.loc[ lnc_ids, : ]
    df_me3 = pd.rpy.common.convert_to_r_dataframe( df_me3 )

    df_me1 = pd.read_table( in_me1, 
                            compression = "gzip", 
                            index_col = 0, 
                            header = 0 )
    df_me1 = df_me1.loc[ lnc_ids, : ]
    df_me1 = pd.rpy.common.convert_to_r_dataframe( df_me1 )

    # write heatmap2 function for plotting matrix and control heatmaps
    R( '''library( gplots )''' )
    R( '''library( RColorBrewer )''' )

    R('''
      plot.heatMaps <- function( df.me3, df.me1, ratios, outfile1, outfile2 ){
        cols= brewer.pal(9, "Blues")
        png( outfile1, width=480, height=960 )
        heatmap.2( as.matrix(df.me3), 
                   trace="none", 
                   dendrogram="none", 
                   Rowv=FALSE, 
                   Colv=FALSE, 
                   labRow=ratios,
                   cexRow=1.25,
                   labCol="",
                   col=cols, 
                   breaks=seq(0, 100, 11),
                   main="H3K4me3" )
        dev.off()
        cols= brewer.pal(5, "Greys")
        png( outfile2, width=480, height=960 )
        heatmap.2( as.matrix(df.me1), 
                   trace="none", 
                   dendrogram="none", 
                   Rowv=FALSE, 
                   Colv=FALSE, 
                   labRow=ratios,
                   cexRow=1.25,
                   labCol="",
                   col=cols, 
                   breaks=seq(0, 50, 10), 
                   main="H3K4me1" )
        dev.off()
        }
      ''')

    # run R plot function from python
    outf1 = P.getTempFilename(".")
    outf2 = P.getTempFilename(".")
    plot_heatMaps = R['plot.heatMaps']  
    plot_heatMaps( df_me3, df_me1, ratios, outf1, outf2 )

    # Save Image of workspace so that heatmap plots can be tweaked if need be
    out_workspace = P.snip(outfile, ".png") + ".RData"
    R('''save.image("%s")''' % out_workspace)

    # merge the two heatmaps on the command line
    statement = ( "montage"
                  " -geometry '1x1+0+0<'"
                  " %(outf1)s"
                  " %(outf2)s"
                  " %(outfile)s" )
    P.run()


def findHighestTSS( in_bed, in_table, out_dir, out_suffix ):
    # generate list of cell types. 
    statement = ( "SELECT DISTINCT cell_type FROM %(in_table)s" % locals() )
    cell_types = PU.fetch( statement )
    cell_types = [ str( x[0] ) for x in cell_types ]  

    # for each cell type, select generate list of lnc_ids
    for cell_type in cell_types:
        statement = ( "SELECT lnc_id"
                      " FROM %(in_table)s" 
                      " WHERE cell_type == '%(cell_type)s'" % locals() )
        lnc_ids = PU.fetch( statement )
        lnc_ids = [ str( x[0] ) for x in lnc_ids ]  

        # iterate through bedfile and write out any tss that appear in lnc_id list
        outf = "".join( [ out_dir, cell_type, out_suffix ] )
        outf = IOTools.openFile( outf, "w" )
        for bed in Bed.iterator( IOTools.openFile( in_bed ) ):
            if bed.fields[0] in lnc_ids:
                outf.write( str( bed ) + "\n" )
            else: 
                continue
        outf.close()


def calculateGeneSetEnrichment( coverage_table, 
                                gene_sets, 
                                outdir,
                                out_suffix,
                                transcript=False ):
    """   
    """
    # select cell types
    statement = ( "SELECT DISTINCT cell_type FROM %(coverage_table)s" % locals() )
    # print( "\n" + statement )
    cell_types = PU.fetch( statement )
    cell_types = [ str( x[0] ) for x in cell_types ]  

    # for genes gene_id is in column lnc_id
    if transcript:
        column = "gene_id"
    else:
        column = "lnc_id"
    # select gene sets
    df_gene_sets = pd.read_table( gene_sets,
                                  compression = "gzip", 
                                  index_col = None, 
                                  header = 0 )
    tmpf1_name = P.getTempFilename(".")
    tmpf1 = open( tmpf1_name, "wb" )
    pickle.dump( df_gene_sets, tmpf1 )
    tmpf1.close()

    # select gene level statistics
    for cell_type in cell_types: 
        # specify outfile
        outfile = cell_type + "_" + out_suffix
        outfile = os.path.join( outdir, outfile )
        # create a dataframe of gene level statistics where rownames are gene_ids
        statement = ( "SELECT %(column)s,ratio"
                      " FROM %(coverage_table)s"
                      " WHERE cell_type='%(cell_type)s' AND pass=1 " % locals() )
        # print( statement + "\n" )
        df_gene_stats = PU.fetch_DataFrame( statement )
        df_gene_stats = df_gene_stats.set_index( column )
        
        tmpfn_name = P.getTempFilename(".")
        tmpfn = open( tmpfn_name, "wb" )
        pickle.dump( df_gene_stats, tmpfn )
        tmpfn.close()

        # Submit wrapper for R GSEA function
        P.submit( "/ifs/devel/projects/proj010/PipelineProj010",
                  "runGSEA",
                  infiles = [ tmpf1_name, tmpfn_name ],
                  outfiles = outfile )
        os.unlink( tmpfn_name )
    os.unlink( tmpf1_name )

def runGSEA( infiles, outfile):
    """
    Run GSEA from r library piano. Pvalue adjustment and nPerm are hardcoded.
    Output is written to table.
    """
    # create outfile for RDS of gsea results object
    outfile1 = P.snip( outfile, "summary_table.tsv" ) + "result.rds"

    # open the two dataframes
    df_gene_sets, df_gene_stats = infiles
    df_gene_sets = pickle.load( open( df_gene_sets ) ) 
    df_gene_stats = pickle.load( open( df_gene_stats ) )
    # turn them into R matrix/dataframe
    df_gene_sets = pd.rpy.common.convert_to_r_dataframe( df_gene_sets ) 
    df_gene_stats = pd.rpy.common.convert_to_r_dataframe( df_gene_stats )

    # Not sure how rpy2 handles module specific classes... do everything else in R
    R( '''library(piano)''' )
    R.assign( "df.gene.sets", df_gene_sets )
    R.assign( "df.gene.stats", df_gene_stats )
    
    # log2 transform dataframe... so foldchange is centred on zero
    R( '''df.gene.stats <- log2(df.gene.stats)''' )
    # establish gene sets... loadGSC requires columns in specific order
    R( '''df.gene.sets <- cbind( df.gene.sets$gene_id, df.gene.sets$gene_set )''' )
    R( '''gene.sets <- loadGSC(df.gene.sets)''' )
    # run GSEA, using a reduced number of permutations. 
    R( '''gsaRes <- runGSA( df.gene.stats, geneSetStat="gsea", adjMethod="none", gsc=gene.sets, nPerm=101 )''' )
    R( '''GSAsummaryTable( gsaRes, save=TRUE, file="%(outfile)s" )''' % locals() )
    R( '''saveRDS( gsaRes, file="%(outfile1)s" )''' % locals() )
    


#################################################################################
#################################################################################
#################################################################################
## section: GAT TFBS overlap testing
#################################################################################

def assignIntervalsToModules( module_table, interval_bed, outfile, bed6=True ):
    """
    Take a table of gene wgcna module assignment, a bedfile of intervals (bed3)
    Return a file that is either bed4 or bed6 (ie stranded) where the 4th column
    is the module assignment for the interval.
    """

    # create dictionary of module assignment
    table = P.snip( os.path.basename( module_table ), ".load" )
    module_dict = {}

    statement = ( "SELECT * FROM %s" % table )
    df = PU.fetch_DataFrame( statement )
    df.set_index("gene_id", inplace=True )

    gene_ids = [str(x) for x in list(df.index)]
    module_ids = df.columns
    for gene_id in gene_ids:
        modules = module_ids[df.loc[gene_id]==1]
        assert len(modules) == 1, "Gene assigned to multiple modules %s" % gene_id
        module_dict[gene_id] = str(modules[0])

    # write bed intervals to outfile, including module assignment
    outf = IOTools.openFile( outfile, "w" )
    for interval in IOTools.openFile( interval_bed ).readlines():
        contig, start, end, gene_id, null, strand = interval.split()
        if gene_id in module_dict.keys():
            if bed6:
                outf.write( "\t".join( [ contig, 
                                         start, 
                                         end, 
                                         module_dict[gene_id], 
                                         ".", 
                                         strand ] )
                            + "\n" )
            else:
                outf.write( "\t".join( [contig, 
                                        start, 
                                        end, 
                                        module_dict[gene_id]] )
                            + "\n" )
        else:
            continue
    outf.close()


################################################################################
# Section: DESeq2 related tasks
################################################################################

def filterCountTable( count_table, 
                      design_table, 
                      biotype_table, 
                      outfile, 
                      biotypes_to_remove = ["antisense", "antisense_downstream"] ):
    """
    Recapitulates filtering of count table in differential expression
    pipeline.
    Removes anything that is antisense
    """
    # clear global namespace... DE pipeline keeps things in namespace
    R( '''rm(list=ls())''' )
    # filter the data according to DE pipeline functions
    Expression.loadTagData( count_table, design_table )
    Expression.filterTagData() # Does everything in global namespace!
    # fetch dataframe into pandas and discard antisense lncRNAs
    # Gives:
    # ValueError: All parameters must be of type Sexp_Type,or Python int/long, float, bool, or None
    # df = com.load_data("countsTable")
    # df.to_csv( outfile, sep="\t", index=False )

    # Write to tempfile. 
    tmpf = P.getTempFilename("/ifs/scratch")
    R('''write.table(countsTable, "%s", quote=FALSE, sep="\t")''' % tmpf )
    R( '''rm(list=ls())''' )

    # get a list of antisense loci
    antisense_list = []
    for locus in IOTools.openFile( biotype_table ).readlines():
        locus = locus.split()
        if locus[2] in biotypes_to_remove:
            antisense_list.append( locus[0] )
        else:
            continue

    # read tempfile and filter out antisense
    outf = IOTools.openFile( outfile, "w" )
    outf_rejected = IOTools.openFile( P.snip(outfile, ".tsv.gz") + ".rejected.tsv.gz", "w" )
    header = True
    for line in IOTools.openFile( tmpf ).readlines():
        if header:
            outf.write( "gene_id\t" + line )
            outf_rejected.write( "gene_id\t" + line )
            header = False
        elif line.split()[0] in antisense_list:
            outf_rejected.write( line )
        else: 
            outf.write( line )
    outf.close()
    outf_rejected.close()

    os.unlink( tmpf )


        
def padjDESeq2Results( infiles ):
    """
    Take the different results, perform multiple testing correction (BH)
    and return a single flatfile containing all results
    WARNING rbind with duplicate rownames causes problems...
    """
    for inf in infiles:
        if re.search("Mat_", os.path.basename(inf)):
            mature_results = os.path.abspath(inf)
        elif re.search("Pro_", os.path.basename(inf)):
            pro_results = os.path.abspath(inf)
        else:
            raise Exception("Don't recognise infile: %s" % inf)

    R('''rm(list=ls())''')
    # add Mat column
    R('''df.mat <- read.table("%(mature_results)s", 
                              header=TRUE,
                              sep=",",
                              row.names=1, 
                              stringsAsFactors=FALSE)''' % locals())
    # R('''df.mat <- df.mat[,-which(names(df.mat) %in% c("padj"))]''')
    R('''df.mat <- cbind(df.mat, Bcell=rep("Mature", length(rownames(df.mat))))''')
    R('''df.mat$gene_id <- rownames(df.mat)''')
    # add Pro column
    R('''df.pro <- read.table("%(pro_results)s", 
                              header=TRUE,
                              sep=",",
                              row.names=1, 
                              stringsAsFactors=FALSE)''' % locals())
    # R('''df.pro <- df.pro[,-which(names(df.pro) %in% c("padj"))]''')
    R('''df.pro <- cbind(df.pro, Bcell=rep("Pro", length(rownames(df.pro))))''')
    R('''df.pro$gene_id <- rownames(df.pro)''')
    # recalculate padj
    R('''df <- rbind(df.mat, df.pro)''')
    R('''df <- cbind(df, padj_global=p.adjust(df$pvalue, method="BH"))''')
    # write to flatfile
    tmpf = P.getTempFilename(".")
    print( tmpf )
    R('''write.table(df, file="%(tmpf)s", quote=FALSE, sep="\t", row.names=FALSE)''' % locals())
    R('''rm(list=ls())''')

    return tmpf

            
def transformCountData( count_file, design_file, outfile, transform="rlog" ):
    """
    Read in count table & design file. 
    Check the transformation to be applied.
    Create a DESeqDataSetObject
    Apply DESeq transformation function
    Plot diagnostic plots suitable for each. 
    """
    if not transform in [ "rlog", "vst", "fpm", "logfpm" ]:
        raise ValueError ( "Unrecognised transformation: %s" % transform )

    # clear R global namespace
    R('''rm(list=ls())''')
    R('''suppressMessages(library("DESeq2"))''')
    R('''suppressMessages(library("vsn"))''')
    
    R('''countData <- read.table( "%(count_file)s",
                                  header=TRUE,
                                  row.names=1 )''' % locals() )
    E.info( "Loaded count data" )
    R('''colData <- read.table( "%(design_file)s",
                                header=TRUE,
                                row.names=1 )''' % locals() )
    E.info( "Loaded design file" )
    R('''dds <- DESeqDataSetFromMatrix(countData = countData,
                                       colData = colData, 
                                       design = ~ condition)''')
    if transform == "rlog":
        R('''sumExp <- rlog(dds, blind=TRUE)''')
        R('''countM <- assay(sumExp)''')
    elif transform == "vst":
        R('''sumExp <- varianceStabilizingTransformation(dds, blind=TRUE)''')
        R('''countM <- assay(sumExp)''')
    elif transform == "fpm":
        R('''countM <- fpm(dds, robust=TRUE)''')
    elif transform == "logfpm":
        R('''countM <- log(fpm(dds, robust=TRUE) + 1)''')
    else:
        raise ValueError( "No transform function for %s" % transform )
    E.info( "Completed %s transformation of count data" % transform )

    # Save transformed count matrix (adding header for row names)
    tmpf = P.getTempFilename("/ifs/scratch")
    R('''write.table(countM, file="%(tmpf)s", quote=FALSE, sep="\t")''' % locals())
    header = True
    outf = IOTools.openFile( outfile, "w" )
    for line in IOTools.openFile( tmpf ).readlines():
        if header:
            line = ["gene_id",] + line.split()
            outf.write( "\t".join(line) + "\n" )
            header = False
        else:
            outf.write( line )
    outf.close()
    os.unlink( tmpf )
    E.info( "Written transformed dataframe to flat file" )

    ## plot diagnostic plots
    # plotting meanSDPlot
    outf_1 = P.snip( outfile, ".tsv.gz" ) + "_meanSDPlot.png"
    R('''notAllZero <- (rowSums(counts(dds))>0)''')
    R('''png("%(outf_1)s")''' % locals() )
    R('''meanSdPlot(countM[notAllZero,], ylim=c(0,2.5))''')
    R('''dev.off()''')
    E.info( "Plotted meanSDPlot for %s transformed data" % transform )

    # plotting hierachical clustering of samples
    outf_2 = P.snip( outfile, ".tsv.gz" ) + "_hclust.png"
    R('''png("%(outf_2)s")''' % locals() )    
    R('''plot(hclust(dist(t(countM))))''')
    R('''dev.off()''')
    E.info( "Plotted hclust dendrogram for %s transformed data" % transform )

    # plotting PCA of vst and rlog transformed data
    if transform in [ "rlog", "vst" ]:
        outf_3 = P.snip( outfile, ".tsv.gz" ) + "_pca500.png"
        R('''png("%(outf_3)s")''' % locals() )    
        R('''plot(plotPCA(sumExp, intgroup="condition"))''')
        R('''dev.off()''')
        E.info( "Plotted pca of top 500 genes for %s transformed data" % transform )
    else:
        pass
    R('''rm(list=ls())''')


def summarizeCountTables( infile, outfile, summary="mean", index_id = "gene_id", rnd = False ):
    """
    Take a table with gene_id as first column and subsequent tissue-condition-replicate
    columns. Take mean/median summary for each tissue-condition and write to outfile. 
    Assumes infile is gzipped.
    If round, then will round values to nearest integer (necessary for count data).
    """

    # read in count table
    df = pd.read_table( infile, 
                        compression = "gzip", 
                        na_values = ["NA", "NaN", "NULL", "nan"], 
                        index_col = 0, 
                        header = 0 )
    
    # extract column headers and remove index id
    sample_ids = list( df.columns )
    assert df.index.name == index_id, "First column is not %s" % index_id

    # generate set of tissue-condition combinations
    sample_ids = set([ "-".join(x.split("-")[:2]) for x in sample_ids ])
    sample_ids = iter( sorted( sample_ids ) )

    # create summary dataframe
    first_id = sample_ids.next()
    df_2 = df.filter( regex = first_id )
    if summary == "mean":
        df_2  = df_2.mean( axis = 1 )
    elif summary == "median":
        df_2 = df_2.median( axis = 1 )
    else:
        raise Exception( "Unrecognised summary value %s" % summary )
    df_2 = pd.DataFrame( df_2, columns = [first_id] )

    # add summaries for remaining tissue-condition combinations
    for tissue_condition in sample_ids:
        E.info( "Processing sample %s" % tissue_condition )
        df_3 = df.filter( regex = tissue_condition )
        if summary == "mean":
            df_3  = df_3.mean( axis = 1 )
        elif summary == "median":
            df_3 = df_3.median( axis = 1 )
        # df.join requires series to have name
        df_3.name = tissue_condition
        df_2 = df_2.join(df_3)

    # round count tables to integer values
    if rnd:
        df_2 = df_2.apply( pd.Series.round )
        df_2 = df_2.apply( lambda x: x.astype(int) )

    # write table
    outf = IOTools.openFile( outfile, "w" )
    df_2.to_csv( outf, sep = "\t", na_rep = "NA", index_label = index_id )


def transformCuffdiffFPKMs( fpkm_file, outfile, constant = 1 ):
    """
    Read table of fpkms into R and output transformed data.
    """
    R('''rm(list=ls())''')
    R('''fpkmData <- read.table( "%(fpkm_file)s",
                                 header=TRUE,
                                 na.strings = c("NA", "NULL"),
                                 row.names=1 ) ''' % locals() )
    R('''fpkmData <- log( fpkmData + %s )''' % str(constant) )
    
    # Save transformed fpkm matrix (adding header for row names)
    tmpf = P.getTempFilename("/ifs/scratch")
    R('''write.table(fpkmData, file="%(tmpf)s", quote=FALSE, sep="\t")''' % locals())
    header = True
    outf = IOTools.openFile( outfile, "w" )
    for line in IOTools.openFile( tmpf ).readlines():
        if header:
            line = ["gene_id",] + [ re.sub("\.", "-", x) for x in line.split() ]
            outf.write( "\t".join(line) + "\n" )
            header = False
        else:
            outf.write( line )
    outf.close()
    os.unlink( tmpf )
    E.info( "Written transformed dataframe to flat file" )


def plotCuffdiffFPKMDispersions( count_file, outfile ):
    R('''rm(list=ls())''')
    R('''suppressMessages(library("vsn"))''')

    R('''countData <- read.table( "%(count_file)s",
                                  header=TRUE,
                                  row.names=1, 
                                  na.string=c("NULL", "NA"))''' % locals() )
    R('''countM <- as.matrix(countData)''')
    R('''png("%(outfile)s")''' % locals() )
    R('''notAllZero <- (rowSums(countM, na.rm = TRUE)>0)''')
    R('''meanSdPlot(log2(countM[notAllZero,] + 1), ylim=c(0,2.5))''')
    R('''dev.off()''')
    E.info( "Plotted meanSDPlot for cuffdiff data" )

    outfile2 = P.snip( outfile, "_meanSDPlot.png" ) + "_hclust.png"
    R('''png("%(outfile2)s")''' % locals() )
    R('''plot(hclust(dist(t(log2(countM + 1)))))''')
    R('''dev.off()''')
    E.info( "Plotted hclust dendrogram for cuffdiff data" )
    R('''rm(list=ls())''')



def fetchPairwiseResults( dds, numerator, denominator, outf_stub, variable="condition" ):
    """
    Receives an RDS containing DESeqDataSet object, a numerator and a denominator, 
    and an outfile_stub. Writes a flat file of the results of Wald test, plus an
    MAPlot.
    """
    R('''rm(list=ls())''')
    R('''suppressMessages(library("DESeq2"))''')

    R('''dds <- readRDS("%(dds)s")''' % locals())
    R('''res <- results( dds, contrast = c( "%(variable)s", "%(numerator)s", "%(denominator)s" ) )''' % locals() )
    R('''png( paste("%(outf_stub)s", "_MAplot.png", sep="") )''' % locals() )
    R('''plotMA(res)''')
    R('''dev.off()''')

    # write data to outfile
    tmpf = P.getTempFilename("/ifs/scratch")
    R('''res_df <- as.data.frame( res )''')
    R('''res_df["numerator"] = rep("%(numerator)s", length(rownames(res_df)))''' % locals())
    R('''res_df["denominator"] = rep("%(denominator)s", length(rownames(res_df)))''' % locals())
    R('''write.table(res_df, file="%(tmpf)s", quote=FALSE, sep="\t")''' % locals())
    header = True
    outf = IOTools.openFile( outf_stub + ".tsv.gz", "w" )
    for line in IOTools.openFile( tmpf ).readlines():
        if header:
            line = ["gene_id",] + line.split()
            outf.write( "\t".join(line) + "\n" )
            header = False
        else:
            outf.write( line )
    outf.close()

    R('''rm(list=ls())''')

    # # pull results into pandas for editing. 
    # df = com.load_data("res_df")
    # df["numerator"] = [ numerator, ]*len( df.index )
    # df["denominator"] = [ denominator, ]*len( df.index )

    # # write results to outfile 
    # outf = outf_stub + ".tsv" 
    # df.to_csv( outfile, sep="\t", index_label="gene_id" )


def correctPairwiseComparisons( infile, outfile ):
    """
    Using BH correction.
    """
    outfile = P.snip( outfile, ".gz" )
    R('''rm(list=ls())''')
    R('''df <- read.table( "%(infile)s", 
                           header=TRUE, 
                           na.string=c("NA", "NaN", "NULL"), 
                           sep="\t")''' % locals() )
    R('''pvals = df[["pvalue"]]''')
    R('''pval_adj = p.adjust(pvals, method="BH")''')
    R('''df["padj_global"] <- pval_adj''')
    R('''write.table(df, file="%(outfile)s", quote=FALSE, sep="\t", row.names=FALSE)''' % locals())

    to_cluster = False
    statement = ("gzip %(outfile)s")
    P.run()

    # Causes error 'rpy2.rinterface.RRuntimeError(Error in (function (data = NA, dim = length(data), dimnames = NULL)  : 
    #                                                   'dims' cannot be of length 0
    # R('''rm(list=ls())''')
    # # fetch pvalues
    # df = pd.read_table( infile, 
    #                     compression = "gzip", 
    #                     na_values = ["NA", "NaN", "NULL"] )
    # E.info( "Read in dataframe" )

    # pvals = df['pvalue'].tolist()
    # E.info( "Extracted %i pvalues" % len(pvals) )

    # # perform BH correction
    # rstats = importr("stats")
    # E.info( "Imported stats to python env" )
    # padj = rstats.p_adjust( pvals, method="BH" )
    # E.info( "Applied BH correction to pvalues" )
    # padj = rpyn.ri2numpy(padj)
    # E.info( "Converted pvals to numpy array" )

    # # add qvalues and write out dataframe
    # df["padj_global"] = padj
    # df.to_csv( outfile, sep="\t", na_rep="NA" )


@cluster_runnable
def plotHeatmapsggplot( infile, outfile, limit ):
    """
    Takes a matrix as tsv, plots as a heat map
    with upper limit specified as limit
    """

    R('''rm(list=ls())''')
    R('''suppressMessages(library("ggplot2"))''')
    R('''suppressMessages(library("reshape"))''')
    R('''df <- read.table("%(infile)s", header=TRUE, stringsAsFactors=FALSE, sep="\t", row.names=1)''' % locals())
    R('''df <- as.matrix(df)''')
    R('''df <- melt(df)''')
    R('''p <- ggplot(df, aes(x=X2, y=X1, fill=value))''')
    R('''p <- p + geom_tile()''')
    R('''p <- p + geom_text(aes(fill=df$value, label=round(df$value, digits=2), alpha=0.3))''')
    R('''p <- p + ylim(rev(c("pro","pre","immature","mature","follicular","marginal","b1a","germinal")))''')
    R('''p <- p + xlim(c("pro","pre","immature","mature","follicular","marginal","b1a","germinal"))''')
    R('''p <- p + scale_fill_gradientn(limits=c(0, %s), colours=c("blue", "white", "red"))''' % str(limit) )
    R('''p <- p + ylab( "Gene Set" )''')
    R('''p <- p + xlab( "H3K4me1:H3K4me3 ratio" )''')
    R('''p <- p + theme_bw() + theme(line=element_blank(), panel.border=element_blank())''')
    R('''p <- p + theme(axis.text=element_text(size=15, angle=45, hjust=1))''')
    R('''png("%(outfile)s")''' % locals())
    R('''print(p)''')
    R('''dev.off()''')


#################################################################################
#################################################################################
#################################################################################
## section: Conservation
#################################################################################
# Things relating to estmating conservation scores.
#################################################################################

def getPhastConsScores( gtf_file, bigwig_file, outfile ):
    """
    The resulting arrays assume exons are sorted by order. 
    It is not necessary for this assumption to hold for the summary.
    """
    # bx.bbi.bigwig_file.BigWigFile
    bw = BigWigFile( open(bigwig_file) )

    # populate dictionary with array of phastcons scores
    score_dict = collections.defaultdict( dict )
    for gtf in GTF.iterator( IOTools.openFile( gtf_file ) ):
        # can't .extend() np.ndarray
        pc_score = bw.get_as_array(gtf.contig, gtf.start, gtf.end).tolist()
        if gtf.gene_id in score_dict:
            if gtf.transcript_id in score_dict[gtf.gene_id]:
                score_dict[gtf.gene_id][gtf.transcript_id].extend( pc_score )
            else:
                score_dict[gtf.gene_id][gtf.transcript_id] = pc_score
        else:
            score_dict[gtf.gene_id][gtf.transcript_id] = pc_score

    # pickle the array in case it is useful downstream. 
    outfile_2 = P.snip( outfile, ".tsv.gz" ) + ".p"
    outf2 = open( outfile_2, "wb" )
    pickle.dump( score_dict, outf2 )
    outf2.close()

    # collapse dict into summary stats and write to outfile
    # NB. to copy collections.dfdct(dict) requires copy.deepcopy()
    outf = IOTools.openFile( outfile, "w" )
    outf.write( "gene_id\t"
                "transcript_id\t"
                "mean_phastcons\t"
                "max_phastcons\t"
                "min_phastcons\t"
                "stddev_phastcons\t"
                "transcript_length\n" )

    for k1, v1 in score_dict.iteritems():
        for k2, v2 in v1.iteritems():
            out_list = [ k1, 
                         k2, 
                         np.mean(v2), 
                         np.max(v2), 
                         np.min(v2), 
                         np.std(v2), 
                         len(v2) ]
            outf.write( "\t".join([ str(x) for x in out_list ] ) + "\n" )
    outf.close()

#################################################################################
#################################################################################
#################################################################################
## section: Plot GViz Tracks
#################################################################################
# This rips Ian's code from CirclularCandidates.py and turns it into a callable
# module (making it worse in the process).
#################################################################################


import rpy2.robjects as ro

Gviz = importr("Gviz")
GenomicFeatures = importr("GenomicFeatures")
AnnotationDbi = importr("AnnotationDbi")

class GenomePlot(object):
    DataTracks = {}
    DataTracks_type = "h"
    # GTF file containing annotations is hardcoded!
    annotations = "/ifs/projects/proj010/analysis_lncrna051/lncRNA_refcoding.gtf.gz"

    gene_track_options = {}
    plot_options = {}
    extend = 100
    height = 7
    width = 14

    def populateDataTracks(self, tracks):
        """
        Assumes track dictionary is passed to call method.
        Simply returns the dictionary.
        """
        return tracks

    def getDataTracks(self, track):
        """
        Return a list of tracks to be passed to GViz plotTracks
        """
        return [data_track["obj"] for data_track in self.DataTracks]

    def __init__(self, *args, **kwargs):
        # Create a TranscriptDb object from lncRNA_refcoding.gtf.gz
        GenomicFeatures = importr("GenomicFeatures")
        txdb = GenomicFeatures.makeTranscriptDbFromGFF(self.annotations,
                                                       format = "gtf",
                                                       exonRankAttribute="exon_number")

    def __call__(self, intervals, tracks, outfile_stub):
        """
        i) Takes nested list of intervals to be plotted [[gene_id, contig, start, end],]
        ii) Takes a nested dict of the files to create tracks for plus type:
        DataTracks = {"First Bam": {"file": "a.bam", "type": "h"},}
        Outputs a single plot for each interval.
        """

        # populate tracks for GViz plot using child class method... 
        # ...or just return input dictionary
        self.DataTracks = self.populateDataTracks(tracks)

        # for each track, create a GViz DataTrack object.
        for track in self.DataTracks:
            # add Gviz DataTrack object to nested dictionary
            self.DataTracks[track]["obj"] = Gviz.DataTrack(
                self.DataTracks[track]["file"],
                # pass type specified in tracks dict, failing that add self.DataTracks_type
                type=self.DataTracks[track].get("type", self.DataTracks_type),
                name=track,
                # passing **dict, which unpacks dict and passes elements as keyword arguments
                **self.DataTracks[track].get("import options", {}))
        E.info( "Created annotation DB object" )

        for i in intervals:
            E.info( "Processing gene_id %s" % i[0] )
            chrom, start, end = i
            # Generate GeneRegionTrack from GTF. 
            gene_track = Gviz.GeneRegionTrack(self.txdb,
                                              chromosome=chrom,
                                              start=start,
                                              end=end,
                                              **self.gene_track_options)

            # Return a list of DataTracks
            data_tracks = self.getDataTracks(track)

            # Create a generic axis track
            axisTrack = Gviz.GenomeAxisTrack()

            # Combine all tracks
            all_tracks = [axisTrack,gene_track] + data_tracks

            # Generate outfile name
            filename = outfile_stub + "_" + gene_id + ".png"

            # Open outfile
            R.png(filename,
                  units="in",
                  res=200,
                  height=self.height,
                  width=self.width)
        
            # plot
            Gviz.plotTracks(all_tracks, main=track, **self.plot_options)
            R["dev.off"]()

#@cluster_runnable
def runGOSeq(infile, gene_sizes, outfile):
    """
    Run GOSeq from python
    Infile contains two columns, first containing gene_ids, second containing
    module assignment. 
    Outputs a table containing GO results, with an additional column specifying
    module id
    """
    R("""rm(list=ls())""")

    grdevices = importr("grDevices")
    goseq = importr("goseq")

    # load file containing module assignment
    df = pd.read_table(infile)
    df.set_index("gene_id", inplace=True)
    # there's only one column per infile
    module = df.columns[0]
    background = map(str, df.index.tolist())

    # calculated gene lengths based on interval data in csvdb
    gene_lengths = []
    for gene in background:
        locus = gene_sizes.loc[gene]
        size = locus[0]
        size = size.split(":")[1]
        size = int(size.split("-")[1]) - int(size.split("-")[0])
        gene_lengths.append(size)
    gene_lengths = robjects.IntVector(gene_lengths)

    # capture plot output by GOSeq
    out_plot = P.snip(infile, ".tsv") + ".png"
    grdevices.png(file=out_plot)

    # run goseq
    E.info("Running GOseq for module: %s" % module)
    genes = df[module]
    genes = robjects.IntVector(genes)
    genes.names = background

    pwf = goseq.nullp(DEgenes=genes, genome="mm10", id="ensGene", bias_data=gene_lengths)#, plot_fit=False)
    GO_wall = goseq.goseq(pwf, "mm10", "ensGene")

    df_out = pandas2ri.ri2pandas(GO_wall)
    df_out["module"] = len(df_out.index)*[module,]
    df_out.to_csv(IOTools.openFile(outfile, "w"), sep="\t", index=False)
    grdevices.dev_off()

    E.info("Completed GOSeq for module: %s" % module)


def getLncRNA_refcoding_miRNA_overlap(lnc_module,
                                      lnc_class,
                                      overlap_tab,
                                      lnc_type="pRNA",
                                      permutations=1000):
    """
    For a particular class of lncRNAs, get the pptn of miRNA sites that overlap
    between lncRNA and all refcoding genes in same wgcna module. 
    (To overcome module size differences, for each lncRNA, get the mean overlap
    with genes in module)
    Permute lnc class assignment to get a null distribution for overlaps. 
    Return the true mean overlap values, the mean overlap values for all 
    permutations, the mean of the mean overlap values for each permutation. 
    """

    def _get_overlap(lnc, module, df):
        df_tmp = df[df["module"] == module]
        overlap = map(float, df_tmp[lnc].tolist())

        return np.mean(overlap)

    # Get the table summarizing lncRNA refcoding miRNA overlap
    df = pd.read_table(overlap_tab, index_col = 0, compression="gzip")

    # Get table of lncRNAs and their associated module
    statement = ("SELECT a.*, b.class_empirical AS classification"
                 " FROM %(lnc_module)s AS a"
                 " INNER JOIN %(lnc_class)s AS b"
                 " ON a.lnc_id = b.gene_id" % locals())
    df_lnc = PU.fetch_DataFrame(statement)

    # There are lncRNAs that are missing from targetscan output, they are all
    # small repeat regions, which presumably don't have any overlaping miRNA
    # targets. They need to be removed from the df_lnc table... 
    from_miRNA_table = set(df.columns.tolist()[1:])
    from_classification_table = set(map(str, df_lnc["lnc_id"].tolist()))
    assert not from_miRNA_table - from_classification_table, \
        "There are lncs in the miRNA output that are missing from chrom class table"
    lost = from_classification_table - from_miRNA_table
    for i in lost:
        E.warn("WARNING: lncRNA %s is missing from targetScan output" % i)
        df_lnc = df_lnc[df_lnc.lnc_id != i]

    ## Calculating TRUE overlap distribution
    true_overlap = []
    # Subset the lncRNAs of interest
    df_pRNA = df_lnc[df_lnc["classification"] == lnc_type]
    # Get lncRNA list
    lncs = map(str, df_pRNA["lnc_id"].tolist())
    # Get module list, removing the ME prefix
    modules = map(str, [x[2:] for x in df_pRNA["module"].tolist()])
    # get lnc module assignment
    lnc_modules = zip(lncs, modules)

    for lnc, module in lnc_modules:
        true_overlap.append(_get_overlap(lnc, module, df))

    ## Calculating the null distribution
    # Get null dataframes
    null_lnc_class = map(str, df_lnc["classification"].tolist())
    df_lnc_null = df_lnc.copy()
    null_dist = []
    null_mean = []

    # Get null distributions
    i = 0
    while i < 100:
        i += 1
        # Create null lncRNA dataframe
        random.shuffle(null_lnc_class)
        df_lnc_null["classification"] = null_lnc_class
    
        # subset null lncRNAs of interest
        df_pRNA_null  = df_lnc_null[df_lnc_null["classification"] == lnc_type]
        # Get null lncRNA list
        lncs_null = map(str, df_pRNA_null["lnc_id"].tolist())
        # Get null module list, removing the ME prefix
        modules_null = map(str, [x[2:] for x in df_pRNA_null["module"].tolist()])
        # get null lnc module assignment
        lnc_modules_null = zip(lncs_null, modules_null)

    
        # calculate null overlap
        null_overlap = []
        for lnc_null, module_null in lnc_modules_null:
            null_overlap.append(_get_overlap(lnc_null, module_null, df))
    
        null_dist.extend(null_overlap)
        null_mean.append(np.mean(null_overlap))

    # create dataframe of true and null 
    level = ["True",]*len(true_overlap) + ["Null",]*len(null_dist)
    overlap=true_overlap + null_dist
    df_bp = pd.DataFrame([level, overlap])
    df_bp = df_bp.transpose()
    df_bp.columns = ["level", "overlap"]

    return true_overlap, null_dist, null_mean, df_bp


@cluster_runnable
def calcLncRNARefcodingOverlap(infiles, outfile):
    """
    For each lncRNA, write out overlap with each refcoding gene. 
    """
    lncRNA, refcoding_dict = infiles
    lnc_id = P.snip(os.path.basename(lncRNA), ".tsv")

    lnc_miRNAs = IOTools.openFile(lncRNA).readline().split()
    pc_dict = pickle.load(open(refcoding_dict, "rb"))


    outf = IOTools.openFile(outfile, "w")
    outf.write("gene_id" + "\t" + lnc_id + "\n")

    n = 0
    for gene in pc_dict.iterkeys():
        n += 1
        if n % 100 == 0 and n > 100:
            E.info("%i genes processed" % n)
        shared = 0
        pc_miRNAs = pc_dict[gene]
        total = len(lnc_miRNAs) + len(pc_miRNAs)
        for miRNA in lnc_miRNAs:
            if miRNA in pc_miRNAs:
                shared += 1
            else:
                continue
        # number of pc miRNA sites also found in lncRNA
        for miRNA in pc_miRNAs:
            if miRNA in lnc_miRNAs:
                shared += 1
            else:
                continue
        # the pptn of total miRNA sites shared between lncRNA and pc gene
        pptn = float(shared)/float(total)
        outf.write(gene + "\t" + str(pptn) + "\n")

    outf.close()

###############################################################################
## Section: Creating UCSC Track Hub
###############################################################################

@cluster_runnable
def bamToBigWig(infile, outfile):
    """
    Run pybedtools bam_to_bigwig on the cluster. 
    """
    bam_to_bigwig(bam=infile, genome="mm10", output=outfile, scale=True)


import CGATPipelines.PipelineUCSC as PipelineUCSC


PROJECT_ROOT = '/ifs/projects'

# gets set by importing script
#PARAMS = {}


def publish_tracks(export_files,
                   params,
                   prefix="",
                   project_id=None,
                   project_name=None,
                   hub_id="ucsc"):
    '''publish a UCSC Track Hub.


    '''
    PARAMS=params

    if not prefix:
        prefix = PARAMS.get("report_prefix", "")

    web_dir = PARAMS["web_dir"]
    if project_id is None:
        project_id = P.getProjectId()
    if project_name is None:
        project_name = P.getProjectName()

    src_export = os.path.abspath("export")
    dest_report = prefix + "report"
    dest_export = prefix + "export"

    hubdir = os.path.join(PARAMS["web_dir"], hub_id)

    if not os.path.exists(hubdir):
        E.info("creating %s" % hubdir)
        os.mkdir(hubdir)

    # write the UCSC hub file
    hubfile = os.path.join(hubdir, "hub.txt")
    genomesfile = os.path.join(hubdir, "genomes.txt")
    trackdir = os.path.join(hubdir, PARAMS["genome"])
    trackfile = os.path.join(hubdir, PARAMS["genome"], "trackDb.txt")
    trackrelpath = os.path.join(PARAMS["genome"], "trackDb.txt")

    if os.path.exists(hubfile):
        with IOTools.openFile(hubfile) as infile:
            hubdata = PipelineUCSC.readUCSCFile(infile)
    else:
        hubdata = [('hub', project_name),
                   ('shortLabel', project_name),
                   ('longLabel', "Data for CGAT collaboration %s" % project_name),
                   ('genomesFile', "genomes.txt"),
                   ('email', 'andreas.heger@gmail.com')]

    E.info("writing to %s" % hubfile)
    with IOTools.openFile(hubfile, "w") as outfile:
        PipelineUCSC.writeUCSCFile(outfile, hubdata)

    # create the genomes.txt file - append to it if necessary.
    if os.path.exists(genomesfile):
        with IOTools.openFile(genomesfile) as infile:
            genomes = PipelineUCSC.readUCSCFile(infile)
    else:
        genomes = []

    if ("genome", PARAMS["genome"]) not in genomes:
        genomes.append(("genome", PARAMS["genome"]))
        genomes.append(("trackDb", trackrelpath))

    E.info("writing to %s" % genomesfile)
    with IOTools.openFile(genomesfile, "w") as outfile:
        PipelineUCSC.writeUCSCFile(outfile, genomes)

    # create the track data
    if not os.path.exists(trackdir):
        os.mkdir(trackdir)

    if os.path.exists(trackfile):
        E.debug('reading existing tracks from %s' % trackfile)
        with IOTools.openFile(trackfile) as infile:
            tracks = PipelineUCSC.readTrackFile(infile)
    else:
        tracks = []

    tracks = collections.OrderedDict(tracks)

    # a dictionary of the track, shortLabel, longlabel for each track
    track_map = {"lncrna_loci": ("LncRNA_Loci", "LncRNA Loci", "LncRNA Loci"),
                 "rnaseq_pro": ("RNASeq_PRO", "RNASeq_PRO", "RNASeq read alignments for pro B cells"),
                 "rnaseq_pre": ("RNASeq_PRE", "RNASeq_PRE", "RNASeq read alignments for pre B cells"),
                 "rnaseq_immature": ("RNASeq_IMM", "RNASeq_IMM", "RNASeq read alignments for immature B cells"),
                 "rnaseq_mature": ("RNASeq_MAT", "RNASeq_MAT", "RNASeq read alignments for mature B cells"),
                 "rnaseq_follicular": ("RNASeq_FO", "RNASeq_FO", "RNASeq read alignments for follicular B cells"),
                 "rnaseq_marginal": ("RNASeq_MZ", "RNASeq_MZ", "RNASeq read alignments for marginal zone B cells"),
                 "rnaseq_b1a": ("RNASeq_B1A", "RNASeq_B1A", "RNASeq read alignments for B1A B cells"),
                 "rnaseq_germinal": ("RNASeq_GC", "RNASeq_GC", "RNASeq read alignments for germinal center B cells"),
                 "chip_pro": ("ChIPSeq_PRO", "ChIPSeq_PRO", "ChIPSeq read alignments for chromatin marks H3K4me1 and H3K4me3 for pro B cells"),
                 "chip_pre": ("ChIPSeq_PRE", "ChIPSeq_PRE", "ChIPSeq read alignments for chromatin marks H3K4me1 and H3K4me3 for pre B cells"),
                 "chip_immature": ("ChIPSeq_IMM", "ChIPSeq_IMM", "ChIPSeq read alignments for chromatin marks H3K4me1 and H3K4me3 for immature B cells"),
                 "chip_mature": ("ChIPSeq_MAT", "ChIPSeq_MAT", "ChIPSeq read alignments for chromatin marks H3K4me1 and H3K4me3 for mature B cells"),
                 "chip_follicular": ("ChIPSeq_FO", "ChIPSeq_FO", "ChIPSeq read alignments for chromatin marks H3K4me1 and H3K4me3 for follicular B cells"),
                 "chip_marginal": ("ChIPSeq_MZ", "ChIPSeq_MZ", "ChIPSeq read alignments for chromatin marks H3K4me1 and H3K4me3 for marginal zone B cells"),
                 "chip_b1a": ("ChIPSeq_B1A", "ChIPSeq_B1A", "ChIPSeq read alignments for chromatin marks H3K4me1 and H3K4me3 for B1A B cells"),
                 "chip_germinal": ("ChIPSeq_GC", "ChIPSeq_GC", "ChIPSeq read alignments for chromatin marks H3K4me1 and H3K4me3 for germinal center B cells"),
                 }

    cell_map = {"pro": "PRO", "pre": "PRE", "immature": "IMM", "mature": "MAT", "follicular": "FO",
                "marginal": "MZ", "b1a": "B1A", "germinal": "GC"}

    # Iterate over export_files. 
    for track_group in export_files.keys():
        # LncRNA loci
        if track_group == "lncrna_loci":
            tracks[track_group] = (("track", track_map[track_group][0]),
                                   ("compositeTrack", "on"),
                                   ("shortLabel", track_map[track_group][1]),
                                   ("longLabel", track_map[track_group][2]),
                                   ("type", "bigBed 12"),
                                   ("visibility", "full"),
                                   ("allButtonPair", "on"))

            for track in export_files[track_group]:
                trackname = os.path.basename(track).split(".")[0]
                dest = os.path.join(trackdir, prefix + os.path.basename(track))
                dest = os.path.abspath(dest)
                # create a symlink
                if not os.path.exists(dest):
                    try:
                        os.symlink(os.path.abspath(track), dest)
                    except OSError, msg:
                        E.warn("could not create symlink from %s to %s: %s" %
                               (os.path.abspath(track), dest, msg))

                if trackname == "lncRNA_gene_models":
                    tracks[trackname] = (("\ttrack", trackname),
                                         ("\tparent", "LncRNA_Loci on"),
                                         ("\tbigDataUrl", os.path.basename(dest)),
                                         ("\tshortLabel", "Gene models"),
                                         ("\tlongLabel", "Collapsed gene models"),
                                         ("\ttype", "bigBed 12"))

                else:
                    assert trackname == "lncRNA_transcript_models"
                    tracks[trackname] = (("\ttrack", trackname),
                                         ("\tparent", "LncRNA_Loci off"),
                                         ("\tbigDataUrl", os.path.basename(dest)),
                                         ("\tshortLabel", "Transcript models"),
                                         ("\tlongLabel", "Transcript models"),
                                         ("\ttype", "bigBed 12"))

        # ChIPSeq peaks
        elif track_group.startswith("H3K4me"):
            chromatin_mark = track_group.split("_")[0]
            group_name = chromatin_mark + "_Peaks"

            if chromatin_mark.endswith("me1"):
                colour = "169,168,168"
            else:
                colour =  "20,78,138"

            # create composite stanza
            tracks[track_group] = (("track", group_name),
                                   ("compositeTrack", "on"),
                                   ("longLabel", "Called peaks for %s" % chromatin_mark),
                                   ("shortLabel", group_name),
                                   ("type", "bigBed 3"),
                                   ("visibility", "dense"))

            # generate stanza for each chipseq alignment
            for track in export_files[track_group]:
                trackname = os.path.basename(track).split(".")[0]
                dest = os.path.join(trackdir, prefix + os.path.basename(track))
                dest = os.path.abspath(dest)
                # create a symlink
                if not os.path.exists(dest):
                    try:
                        os.symlink(os.path.abspath(track), dest)
                    except OSError, msg:
                        E.warn("could not create symlink from %s to %s: %s" %
                               (os.path.abspath(track), dest, msg))

                # only show the merged intervals by default.
                if trackname.startswith("H3K4me"):
                    parent = group_name + " on"
                else:
                    parent = group_name + " off"
                    # tweak trackname
                    cell, mark = trackname.split("-")
                    trackname = mark + "_" + cell_map[cell]

                tracks[trackname] = (("\ttrack", trackname + "_peaks"),
                                     ("\tparent", parent),
                                     ("\tbigDataUrl", os.path.basename(dest)),
                                     ("\tshortLabel", trackname),
                                     ("\tlongLabel", trackname),
                                     ("\ttype", "bigBed 3"),
                                     ("\tcolor", colour))

        # RNASeq alignments
        elif track_group.startswith("rnaseq"):
            group_name, short_label, long_label = track_map[track_group]
            
            if re.search("B1A", group_name):
                vis = "full"
            else:
                vis = "hide"

            # create composite stanza
            tracks[track_group] = (("track", group_name),
                                   ("compositeTrack", "on"),
                                   ("longLabel", long_label),
                                   ("shortLabel", short_label),
                                   ("smoothingWindow", "4"),
                                   ("autoScale", "on"),
                                   ("type", "bigWig"),
                                   ("visibility", vis),
                                   ("allButtonPair", "on"),
                                   ("color", "209,0,28"))

            # generate stanza for each rnaseq alignment
            for track in export_files[track_group]:
                trackname = os.path.basename(track).split(".")[0]
                dest = os.path.join(trackdir, prefix + os.path.basename(track))
                dest = os.path.abspath(dest)
                # create a symlink
                if not os.path.exists(dest):
                    try:
                        os.symlink(os.path.abspath(track), dest)
                    except OSError, msg:
                        E.warn("could not create symlink from %s to %s: %s" %
                               (os.path.abspath(track), dest, msg))

                
                # only show the merged files for a single cell type (B1A) by defauld
                if trackname.endswith("R0"):
                    parent = group_name + " on"
                else:
                    parent = group_name + " off"
                tracks[trackname] = (("\ttrack", trackname),
                                     ("\tparent", parent),
                                     ("\tbigDataUrl", os.path.basename(dest)),
                                     ("\tshortLabel", trackname),
                                     ("\tlongLabel", trackname),
                                     ("\ttype", "bigWig"))
                                  
        elif track_group.startswith("chip"):
            group_name, short_label, long_label = track_map[track_group]

            if re.search("B1A", group_name):
                vis = "full"
            else:
                vis = "hide"
                
            # create multiWig stanza
            tracks[track_group] = (("track", group_name),
                                   ("longLabel", long_label),
                                   ("shortLabel", short_label),
                                   ("container", "multiWig"),
                                   ("aggregate", "transparentOverlay"),
                                   ("showSubtrackColorOnUi", "on"),
                                   ("smoothingWindow", "4"),
                                   ("autoScale", "on"),
                                   ("type", "bigWig"),
                                   ("visibility", vis))

            for track in export_files[track_group]:
                trackname = os.path.basename(track).split(".")[0]
                dest = os.path.join(trackdir, prefix + os.path.basename(track))
                dest = os.path.abspath(dest)
                # create a symlink
                if not os.path.exists(dest):
                    try:
                        os.symlink(os.path.abspath(track), dest)
                    except OSError, msg:
                        E.warn("could not create symlink from %s to %s: %s" %
                               (os.path.abspath(track), dest, msg))

            
                # get sample specific details
                cell = trackname.split("-")[0]
                chromatin_mark = trackname.split("-")[1]
                new_name = chromatin_mark + "_" + cell_map[cell]
                long_label = "ChIPSeq read alignment for %s" % chromatin_mark
                if chromatin_mark == "H3K4me1":
                    colour = "169,168,168"
                else:
                    colour =  "20,78,138"

                tracks[trackname] = (("\ttrack", new_name),
                                     ("\tparent", group_name),
                                     ("\tbigDataUrl", os.path.basename(dest)),
                                     ("\tshortLabel", new_name),
                                     ("\tlongLabel", long_label),
                                     ("\tcolor", colour),
                                     ("\ttype", "bigWig"))
                            
    with IOTools.openFile(trackfile, "w") as outfile:
        for track, trackdata in tracks.iteritems():
            for line in trackdata:
                outfile.write(" ".join((line[0], line[1])) + "\n")
            outfile.write("\n")
