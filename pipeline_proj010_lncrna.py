################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_snps.py 2870 2010-03-03 10:20:29Z andreas $
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
===========================
pipeline_proj010_lncrna.py
===========================

:Author: Jethro Johnson
:Release: $Id$
:Date: |today|
:Tags: Python

A pipeline for the analysis of lncRNAs in B cells

Overview
========

This pipeline contains code used to produce figures and results
presented in the study by Brazao et al. 'LncRNAs in early B cell
development and activation'.

Usage
=====

Running this pipeline requires installation of the CGAT code collection
(github.com/CGATOxford/cgat) and the CGAT pipeline collection 
(github.com/CGATOxford/CGATPipelines). 

Parameters for specific analyses are specified in the accompanying
configuration file (./pipeline_proj010_lncrna/pipeline.ini)

Code
====

"""

from ruffus import *
from ruffus.combinatorics import *
from scipy import stats

import sys, glob, gzip, os, itertools, re, math, types, collections, time
import optparse, shutil
import sqlite3
import pysam
import random
import tempfile
import numpy as np
import scipy
import pickle
import collections
import pandas as pd
#import pandas.rpy.common as com
import rpy2.robjects as robjects
import rpy2.robjects.pandas2ri as pandas2ri
from rpy2.robjects import r as R
from rpy2.robjects import numpy2ri as rpyn
from rpy2.robjects.packages import importr

import CGAT.IOTools as IOTools
import CGAT.Experiment as E
import CGAT.Database as Database
import CGAT.Expression as Expression
import CGAT.GTF as GTF
import CGAT.Bed as Bed
import CGAT.FastaIterator as FastaIterator
import CGATPipelines.PipelineLncRNA as PipelineLncRNA
import CGATPipelines.PipelineUtilities as PU
import CGATPipelines.PipelineRnaseq as PipelineRnaseq
import CGATPipelines.PipelinePublishing as PipelinePublishing

import PipelineProj010QC as QC
import PipelineProj010 as P10

#################################################################################
#################################################################################
#################################################################################
# TO DO:

#################################################################################
#################################################################################
#################################################################################
## Pipeline configuration
#################################################################################

# load options from the config file
import CGATPipelines.Pipeline as P
P.getParameters( 
    ["%s/pipeline.ini" % __file__[:-len(".py")],
     "../pipeline.ini",
     "pipeline.ini" ] )

PARAMS = P.PARAMS
PARAMS_ANNOTATIONS = P.peekParameters( PARAMS["annotations_dir"],
                                       "pipeline_annotations.py" )

#################################################################################
#################################################################################
## Helper functions mapping tracks to conditions, etc
#################################################################################

import CGATPipelines.PipelineTracks as PipelineTracks

Sample = PipelineTracks.Sample3

# define tracks based on all sample bamfiles that are not input or index
TRACKS = PipelineTracks.Tracks( Sample ).loadFromDirectory( 
    os.listdir( PARAMS[ "location_bamfiles"] ),
    "(\S+).star.bam", 
    exclude = ["(\S+).bai"] )

# list tissue, condition, experiment
TISSUES = PipelineTracks.Aggregate( TRACKS, labels = ("tissue", ) )
CONDITIONS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", ) )
EXPERIMENTS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", "tissue" ) )


#################################################################################
#################################################################################
#################################################################################
def connect():
    '''connect to database.

    Use this method to connect to additional databases.

    Returns a database connection.
    '''

    dbh = sqlite3.connect( PARAMS["database"] )
    statement = '''ATTACH DATABASE '%s' as annotations''' % (PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute( statement )
    cc.close()

    return dbh

#################################################################################
#### METASECTION #### Process & Filter LncRNAs from discovery Pipeline #### 
# * For historical reasons single-exon and multi-exon lncRNA loci are processed 
# separately. This has been left in the pipeline in case any of the filtering 
# thresholds need to be set differently for the two categories. 
# * There is an issue with the way that cufflinks originally assigns gene_ids to 
# lncRNAs. If a transcript that joins two loci is filtered out following naming, 
# then these loci will continue to have the same gene_id, even if they no longer 
# overlap. For this reason: i) non-overlapping loci are renamed at the start of
# the processing pipeline; ii) are re-classified as single exon or multi-exon;
# are subsequently divided so that me and se lncRNA can be filtered separately. 
#################################################################################
#################################################################################
#################################################################################
# Section: resolve lncRNA naming issues
#################################################################################
#################################################################################
## subsection: split lncRNA gtf on the basis of exon_status_locus
#################################################################################
# create separate gtfs for single and multi-exon loci
@follows( mkdir( "./filter_se_lncrna" ), mkdir( "./filter_me_lncrna" ) )
@split( os.path.join( PARAMS[ "location_lncrnafiles" ], "lncrna_final.gtf.gz" ), 
            regex( "(.+)/(.+)_final.gtf.gz" ),
            [r"./filter_se_lncrna/\2_se.gtf.gz", 
             r"./filter_me_lncrna/\2_me.gtf.gz"] )
def splitMEAndSELncRNA( infile, outfiles ):
    """
    Separates lncRNA loci on the basis of gtf attribute exon_status_locus
    """
    outf_se, outf_me = outfiles
    tmpf_se = P.getTempFile( os.path.dirname( outf_se ) )
    tmpf_me = P.getTempFile( os.path.dirname( outf_me ) )

    # separate intervals in gtf file on the basis of exon_status_locus
    for gtf in GTF.iterator( IOTools.openFile( infile ) ):
        if gtf.exon_status_locus == "s":
            tmpf_se.write( "%s\n" % str( gtf ) )
        elif gtf.exon_status_locus == "m":
            tmpf_me.write( "%s\n" % str( gtf ) )
        else: raise ValueError( "Unrecognized gtf attribute %s" 
                                % gtf.exon_status_locus )
    tmpf_se.close()
    tmpf_me.close()
            
    # sort resulting gtfs by gene
    for gtf in [ ( tmpf_se.name, outf_se ), ( tmpf_me.name, outf_me ) ]:
        inf, outf = gtf
        statement = ( "cat %(inf)s |"
                      " python %(scriptsdir)s/gtf2gtf.py"
                      "  --method=sort --sort-order=gene+transcript"
                      "  --log=%(outf)s.log |"
                      " gzip >  %(outf)s" )
        P.run()

    os.unlink( tmpf_se.name )
    os.unlink( tmpf_me.name )

#################################################################################
## subsection: rename multi-exon lncRNA loci
#################################################################################
@transform( splitMEAndSELncRNA, 
            regex( "(.+)/(.+)_me.gtf.gz" ), 
            add_inputs( r"filter_se_lncrna/\2_se.gtf.gz" ), 
            r"\1/\2_me_renamed01.gtf.gz" )
def renameMELncRNA( infiles, outfile ):
    """
    Assigns multi-exon lncRNAs new gene_ids, old gene_ids are passed to gene_oId.
    Gene models containing space not covered by either an exon or intron are split
    Requires input gtf to be sorted by gene+transcript for transcript_iterator
    to work
    """
    # specify infiles
    me_lnc_gtf, se_lnc_gtf = infiles
    tmp_gtf = P.getTempFilename( "." )
    tmp_gtf2 = P.getTempFilename( "." )
    
    # open outfiles
    outf = IOTools.openFile( tmp_gtf2, "w" )
    outf_split = P.snip( outfile, ".gtf.gz" ) + "_split.tsv.gz"
    outf_split = IOTools.openFile( outf_split, "w" )

    # generate template for new naming conventions
    MEG = 0
    MET = 0
    # number of se transcripts
    n_se_t = len( [ x for x in IOTools.openFile( se_lnc_gtf ).readlines() ] )
    # number of me transcripts
    n_me_t = len( [ x for x in GTF.transcript_iterator( GTF.iterator(IOTools.openFile( me_lnc_gtf ) ) ) ] )
    len_z = len( str( n_se_t + n_me_t  ) )
    def _getGeneID( meg ):
        return "LNCGme" + str( meg ).zfill( len_z )

    def _getTranscriptID( met ):
        return "LNCTme" + str( met ).zfill( len_z )

    # sort me lncRNA by gene_id 
    statement = ( "zcat %(me_lnc_gtf)s |"
                  " python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=sort --sort-order=gene+transcript"
                  "  --log=%(outfile)s.log"
                  "   > %(tmp_gtf)s" )
    P.run()

    # create dictionary to record new gene_id for each transcript
    transcript_gene_ids = {}
    split_gene_ids = []

    # iterate through me lnc gtf, retrieve transcript co-ordinates for each gene
    for flat_gene in GTF.flat_gene_iterator( GTF.iterator( IOTools.openFile ( tmp_gtf ) ) ):
        # return ordered dictionary with k = transcript_id, v = [start, stop, len]
        transcript_dict = P10.defineTranscripts( flat_gene )
        trans_list = transcript_dict.items()

        # set the current gene_id
        MEG += 1
        curr_gene_id = _getGeneID( MEG )
        # set current transcript_id and coordinates as first in flat gene list
        curr_tran_id = trans_list[0][0]
        curr_tran_start = trans_list[0][1][0]
        curr_tran_end = trans_list[0][1][1]

        # assign current gene_id to current transcript
        transcript_gene_ids[ curr_tran_id ] = curr_gene_id

        # iterate through subsequent transcripts in list
        for next_transcript in trans_list[1:]:
            next_tran_id = next_transcript[0]
            next_tran_start = next_transcript[1][0]
            next_tran_end = next_transcript[1][1]

            assert curr_tran_start <= next_tran_start, "transcripts are not sorted"
            assert curr_tran_start < curr_tran_end
    
            # if next transcript overlaps current transcript, 
            # assign them the same gene_id
            if curr_tran_end >= next_tran_start:
                transcript_gene_ids[ next_tran_id ] = curr_gene_id
                # check if current transcript encompasses next transcript
                if curr_tran_end >= next_tran_end:
                    continue
                else:
                    curr_tran_id = next_tran_id
                    curr_tran_start = next_tran_start
                    curr_tran_end = next_tran_end
            # if the next transcript doesn't overlap current transcript, 
            # assign subsequent transcripts a new gene_id
            else:
                split_ids = []
                split_ids.append( curr_gene_id )
                MEG += 1
                curr_gene_id = _getGeneID( MEG )
                split_ids.append( curr_gene_id )
                split_gene_ids.append( split_ids )                
                transcript_gene_ids[ next_tran_id ] = curr_gene_id
                curr_tran_id = next_tran_id
                curr_tran_start = next_tran_start
                curr_tran_end = next_tran_end

    # iterate through transcripts
    for gtf in GTF.transcript_iterator( GTF.iterator( IOTools.openFile( tmp_gtf ) ) ):
        # set new transcript id and retrieve gene_id
        MET += 1
        transcript_id = _getTranscriptID( MET )
        gene_id = transcript_gene_ids[ gtf[0].transcript_id ]
    
        # iterate through exons
        for exon in gtf:
            # set new gene_id
            exon.setAttribute( "gene_oId", exon.gene_id )
            exon.setAttribute( "gene_id", gene_id )
            # set new transcript_id
            exon.setAttribute( "transcript_oId", exon.transcript_id )
            exon.setAttribute( "transcript_id", transcript_id )
            outf.write( str( exon ) + "\n" )
        
    # write file of newly split gene_ids
    for gene_ids in split_gene_ids:
        outf_split.write( "\t".join( gene_ids ) + "\n" )

    outf.close()
    outf_split.close()

    # sort the renamed gtf by gene_id
    statement = ( "cat %(tmp_gtf2)s |" 
                  " python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=sort --sort-order=gene+transcript"
                  "  --log=%(outfile)s.log |"
                  " gzip  > %(outfile)s" )
    P.run()

    os.unlink( tmp_gtf )
    os.unlink( tmp_gtf2 )

    E.info( "%s multi exon lncRNA renamed\n"
            "%s multi exon transcript renamed\n"
            "%s multi exon genes arising from split gene models"
            % ( str( MEG ), str( MET ), str( len( split_gene_ids ) ) ) ) 


@transform( renameMELncRNA, 
            regex( "(.+)/(.+)_me_renamed01.gtf.gz" ), 
            r"\1/\2_me_renamed02.gtf.gz" )
def correctRenamedMELncRNA( infile, outfile ):
    """
    Some of the multi-exon gene models split by renameMELncRNA result in loci
    that contain only single exon transcripts. These single-exon loci need to 
    be renamed. If single exon loci contain more than one overlapping 
    transcript, they're renamed to sm. Otherwise they're renamed to se.
    Note: gene_iterator calls transcript_iterator, which requires intervals
    to be sorted by gene+transcript
    """
    tmpf1 = P.getTempFilename(".")
    tmpf2 = P.getTempFilename(".")
    outf = IOTools.openFile(tmpf2, "w")
    outf2 = IOTools.openFile( re.sub( "02", "03", outfile ), "w" )

    # sort infile
    statement = ( "zcat %(infile)s |"
                  " python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=sort --sort-order=gene+transcript"
                  "  --log=%(outfile)s.log"
                  " > %(tmpf1)s" )
    P.run()

    for gene in GTF.gene_iterator( GTF.iterator( IOTools.openFile( tmpf1 ) ) ):
        # summarise exon_status of all transcripts in a gene
        exon_status_transcripts = []
        for transcript in gene:
            assert transcript[0].exon_status in [ "m", "s" ], \
            "Unknown exon_status %s" % transcript[0].exon_status
            exon_status_transcripts.append( transcript[0].exon_status ) 
        # ignore genes with multi-exon transcripts
        if "m" in exon_status_transcripts:
            for transcript in gene:
                for exon in transcript:
                    assert exon.exon_status_locus == "m"
                    outf.write( str( exon ) + "\n" )
            continue
        else:
            # reset gene_id and exon_status_locus for single-exon genes
            if len( exon_status_transcripts ) > 1:
                gene_rename = "sm"
            else:
                gene_rename = "se"
            for transcript in gene: 
                for exon in transcript:
                    new_gene_id = re.sub( "me", gene_rename, exon.gene_id )
                    new_tran_id = re.sub( "me", "se", exon.transcript_id )
                    exon.setAttribute( "gene_id", new_gene_id )
                    exon.setAttribute( "transcript_id", new_tran_id )
                    exon.setAttribute( "exon_status_locus", "s" )
                    outf.write( str( exon ) + "\n" )
                    outf2.write( str( exon ) + "\n" )
    outf.close()
    outf2.close()

    statement = ( "cat %(tmpf2)s |"
                  " python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=sort --sort-order=gene+transcript"
                  "  --log=%(outfile)s.log |"
                  " gzip > %(outfile)s" )
    P.run()

    os.unlink(tmpf1)
    os.unlink(tmpf2)

#################################################################################
## subsection: rename single-exon lncRNA loci
#################################################################################
@transform( splitMEAndSELncRNA, 
            suffix( "_se.gtf.gz" ),
            add_inputs( renameMELncRNA ),
            "_se_renamed01.gtf.gz" )
def renameSELncRNA( infiles, outfile ):
    # specify infiles
    se_lnc_gtf, me_lnc_gtf = infiles
    tmp_gtf = P.getTempFilename( "." )

    # generate template for new naming conventions
    # establish how many me_lncRNA genes and transcripts there are
    SEG = len( [ x for x in GTF.flat_gene_iterator( GTF.iterator( IOTools.openFile( me_lnc_gtf ) ) ) ] )
    SET = len( [ x for x in GTF.transcript_iterator( GTF.iterator( IOTools.openFile( me_lnc_gtf ) ) ) ] )
    # number of se transcripts
    n_se_t = len( [ x for x in IOTools.openFile( se_lnc_gtf ).readlines() ] )
    # number of me transcripts
    n_me_t = len( [ x for x in GTF.transcript_iterator( GTF.iterator(IOTools.openFile( me_lnc_gtf ) ) ) ] )
    len_z = len( str( n_se_t + n_me_t  ) )
    def _getGeneID( SEG ):
        return "LNCGse" + str( SEG ).zfill( len_z )

    def _getTranscriptID( SET ):
        return "LNCTse" + str( SET ).zfill( len_z )

    # sort se lncRNA by gene id
    statement = ( "zcat %(se_lnc_gtf)s |"
                  " python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=sort --sort-order=gene+transcript"
                  "  --log=%(outfile)s.log"
                  "   > %(tmp_gtf)s" )
    P.run()

    # open outfiles
    outf = IOTools.openFile( outfile, "w" )
    outf_split = P.snip( outfile, ".gtf.gz" ) + "_split.tsv.gz"
    outf_split = IOTools.openFile( outf_split, "w" )


    # create dictionary to record new gene_id for each transcript
    transcript_gene_ids = {}
    split_gene_ids = []

    # iterate through me lnc gtf, retrieve transcript co-ordinates for each gene
    for flat_gene in GTF.flat_gene_iterator( GTF.iterator( IOTools.openFile ( tmp_gtf ) ) ):
        # return ordered dictionary with k = transcript_id, v = [start, stop, len]
        transcript_dict = P10.defineTranscripts( flat_gene )
        trans_list = transcript_dict.items()

        # set the current gene_id
        SEG += 1
        curr_gene_id = _getGeneID( SEG )
        # set current transcript_id and coordinates as first in flat gene list
        curr_tran_id = trans_list[0][0]
        curr_tran_start = trans_list[0][1][0]
        curr_tran_end = trans_list[0][1][1]
    
        # assign current gene_id to current transcript
        transcript_gene_ids[ curr_tran_id ] = curr_gene_id

        # iterate through subsequent transcripts in list
        for next_transcript in trans_list[1:]:
            next_tran_id = next_transcript[0]
            next_tran_start = next_transcript[1][0]
            next_tran_end = next_transcript[1][1]

            assert curr_tran_start <= next_tran_start, "transcripts are not sorted"
            assert curr_tran_start < curr_tran_end
    
            # if next transcript overlaps current transcript, 
            # assign them the same gene_id
            if curr_tran_end >= next_tran_start:
                transcript_gene_ids[ next_tran_id ] = curr_gene_id
                # check if current transcript encompasses next transcript
                if curr_tran_end >= next_tran_end:
                    continue
                else:
                    curr_tran_id = next_tran_id
                    curr_tran_start = next_tran_start
                    curr_tran_end = next_tran_end
            # if the next transcript doesn't overlap current transcript, 
            # assign subsequent transcripts a new gene_id
            else:
                split_ids = []
                split_ids.append( curr_gene_id )
                SEG += 1
                curr_gene_id = _getGeneID( SEG )
                split_ids.append( curr_gene_id )
                split_gene_ids.append( split_ids )                
                transcript_gene_ids[ next_tran_id ] = curr_gene_id
                curr_tran_id = next_tran_id
                curr_tran_start = next_tran_start
                curr_tran_end = next_tran_end

    # iterate through transcripts
    for gtf in GTF.transcript_iterator( GTF.iterator( IOTools.openFile( tmp_gtf ) ) ):
        # set new transcript id and retrieve gene_id
        SET += 1
        transcript_id = _getTranscriptID( SET )
        gene_id = transcript_gene_ids[ gtf[0].transcript_id ]
    
        # iterate through exons
        for exon in gtf:
            # set new gene_id
            exon.setAttribute( "gene_oId", exon.gene_id )
            exon.setAttribute( "gene_id", gene_id )
            # set new transcript_id
            exon.setAttribute( "transcript_oId", exon.transcript_id )
            exon.setAttribute( "transcript_id", transcript_id )
            outf.write( str( exon ) + "\n" )

        
    # write file of newly split gene_ids
    for gene_ids in split_gene_ids:
        outf_split.write( "\t".join( gene_ids ) + "\n" )

    outf.close()
    outf_split.close()
    os.unlink( tmp_gtf )

    ng = SEG - len( [ x for x in GTF.flat_gene_iterator( GTF.iterator( IOTools.openFile( me_lnc_gtf ) ) ) ] )
    nt = SET - len( [ x for x in GTF.transcript_iterator( GTF.iterator( IOTools.openFile( me_lnc_gtf ) ) ) ] )
    E.info( "%s multi exon lncRNA renamed\n"
            "%s multi exon transcript renamed\n"
            "%s multi exon genes arising from split gene models"
            % ( str( ng ), str( nt ), str( len( split_gene_ids ) ) ) ) 


@transform( renameSELncRNA, 
            regex( "(.+)/(.+)_se_renamed01.gtf.gz" ), 
            r"\1/\2_se_renamed02.gtf.gz" )
def correctRenamedSELncRNA( infile, outfile ):
    """
    Some of the single-exon gene models contain more than one overlapping 
    transcript. The gene_ids for these loci are renamed from se to sm. 
    """
    tmpf1 = P.getTempFilename( "." )
    outf1 = IOTools.openFile( outfile, "w" )
    outf2 = IOTools.openFile( re.sub( "02", "03", outfile ), "w" )

    # sort infile
    statement = ( "zcat %(infile)s |"
                  " python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=sort --sort-order=gene+transcript"
                  "  --log=%(outfile)s.log"
                  " > %(tmpf1)s" )
    P.run()

    for gene in GTF.gene_iterator( GTF.iterator( IOTools.openFile( tmpf1 ) ) ):
        n_transcripts = 0
        for transcript in gene:
            n_transcripts += 1
            assert len( transcript ) == 1, "File contains multi-exon transcripts"
            assert transcript[0].exon_status_locus == "s"
        if n_transcripts > 1:
            for transcript in gene:
                for exon in transcript:
                    new_gene_id = re.sub( "se", "sm", exon.gene_id )
                    exon.setAttribute( "gene_id", new_gene_id )
                    outf1.write( str( exon ) + "\n" )
                    outf2.write( str( exon ) + "\n" ) 
        else:
            for transcript in gene:
                for exon in transcript:
                    outf1.write( str( exon ) + "\n" )

    outf1.close()
    outf2.close()

#################################################################################
## subsection: re-split multi-exon and single-exon lncRNA 
#################################################################################
@follows( mkdir( "./filter_lncrna" ) )
@merge( [ correctRenamedMELncRNA, correctRenamedSELncRNA ],
        "./filter_lncrna/lncrna_combined_pre_filtering.gtf.gz")
def mergeRenamedLncRNA( infiles, outfile ):
    """
    Following renaming, some loci in the multi-exon gtf contain only
    single-exon transcripts. 
    """
    inf = " ".join( infiles )
    statement = ( "zcat %(inf)s |"
                  " python %(scriptsdir)s/gtf2gtf.py "
                  "  --method=sort --sort-order=gene+transcript"
                  "  --log=%(outfile)s.log |"
                  " gzip > %(outfile)s" )
    P.run()


@split( mergeRenamedLncRNA, 
        regex( "(.+)/(.+)_combined_pre_filtering.gtf.gz" ),
        [ r"./filter_me_lncrna/\2_me_filtered_0.gtf.gz", 
          r"./filter_se_lncrna/\2_se_filtered_0.gtf.gz" ] )
def splitRenamedLncRNA( infile, outfiles ):
    """
    Split multi-exon and single-exon lncRNA on the basis of gene_id
    """
    me_out, se_out = outfiles
    me_out = IOTools.openFile( me_out, "w" )
    se_out = IOTools.openFile( se_out, "w" )

    for exon in GTF.iterator( IOTools.openFile( infile ) ):
        if re.match( "LNCGm", exon.gene_id ):
            me_out.write( str( exon ) + "\n" )
        elif re.match( "LNCGs", exon.gene_id ):
            se_out.write( str( exon ) + "\n" )
        else:
            raise ValueError( "Unrecognised gene_id: %s" % exon.gene_id )

    me_out.close()
    se_out.close()


@follows( splitRenamedLncRNA )
def renameLncRNA(): pass

#################################################################################
#################################################################################
#################################################################################
# Section: filter multi-exon lncRNA loci
#################################################################################
# Filter anything that is within 5kb sense upstream/downstream of pc gene model
# liftover human 5' 3'UTR and remove anything within these regions
# remove anything within 20(?)kb downstream of pc gene, with >60% intergenic cov
#################################################################################
#################################################################################
## subsection: classify lncRNA loci in relation to refcoding/reference genesets
#################################################################################
@transform( splitRenamedLncRNA,
            suffix( "_me_filtered_0.gtf.gz" ),
            add_inputs( os.path.join( PARAMS["location_transcriptfiles"],
                                      "reference.gtf.gz" ) ), 
            r"_me_classified.gtf.gz" )
def classifyVSReference( infiles, outfile ):
    lncRNA_gtf, reference_gtf = infiles

    dist_upstr = PARAMS[ "filter_distance_upstream" ]
    dist_dstr = PARAMS[ "filter_distance_downstream" ]
    reference_biotypes = PARAMS[ "filter_classify_vs" ].split( "," )

    # preprocess reference file
    tmpf_ref = P.getTempFile( "./filter_me_lncrna" )
    tmpf_ref_name = tmpf_ref.name
    for gtf in GTF.iterator( IOTools.openFile( reference_gtf ) ):
        if gtf.source in reference_biotypes:
            if gtf.feature == "exon": 
                tmpf_ref.write( str( gtf ) + "\n" )
            else: continue
        else: continue
    tmpf_ref.close()

    tempdir = PipelineLncRNA.reClassifyLncRNAGenes( lncRNA_gtf,
                                                    tmpf_ref_name, 
                                                    outfile, 
                                                    dist_upstr,
                                                    dist_dstr, 
                                                    wdir = "./filter_me_lncrna")
    to_cluster=False
    tmpf_ref_name = os.path.basename( tmpf_ref_name )
    statement = ( "mv ./filter_me_lncrna/%(tmpf_ref_name)s "
                  "%(tempdir)s/%(tmpf_ref_name)s;"
                  " printf 'tempfiles in this dir have been kept in case it is"
                  " necessary to check that the classification is working"
                  " correctly. \\n%(tempdir)s - contains a record the intervals"
                  " used to classify each lncRNA locus together with the gtf"
                  " entries for the relevant lncRNA and reference exons"
                  " \\n%(tempdir)s/%(tmpf_ref_name)s - contains all the gtf"
                  " entries used to classify lncRNAs\\n'"
                  " >> ./filter_me_lncrna/README" )
    P.run()
    
#    os.unlink( tmpf_ref_name )
#    shutil.rmtree( tempdir )


@transform( classifyVSReference, 
            suffix( "_classified.gtf.gz" ), 
            "_filtered_1.gtf.gz" )
def filterBasedOnClassification( infile, outfile ):
    rejected = P.snip( outfile, "_filtered_1.gtf.gz" ) + "_rejected_1.gtf.gz"
    inf = IOTools.openFile( infile )
    outf = IOTools.openFile( outfile, "w" )
    outf_reject = IOTools.openFile( rejected, "w" )
    to_reject = PARAMS[ "filter_reject" ].split( "," )
    
    rejected = 0
    
    for gtf in GTF.transcript_iterator( GTF.iterator( inf ) ):
        if gtf[0].source in to_reject:
            rejected += 1
            for exon in IOTools.flatten( gtf ):
                outf_reject.write( str( exon ) + "\n" )
        else:
            for exon in IOTools.flatten( gtf ):
                outf.write( str( exon ) + "\n" )

    E.info( "%i lncRNA transcripts rejected based on their"
            " classification" % rejected  )

    outf.close()
    outf_reject.close()

#################################################################################
## subsection: filter lncRNA that intersect undesired ensembl biotypes
#################################################################################
# This is as sanity check as intervals overlapping these biotypes should have 
# been removed as part of the lncRNA pipeline (catches some pseudogenes).

@transform( filterBasedOnClassification, 
            suffix( "_filtered_1.gtf.gz" ),
            add_inputs( os.path.join( PARAMS["location_transcriptfiles"],
                                      "reference.gtf.gz" ), 
                        os.path.join( PARAMS[ "annotations_dir" ], 
                                      PARAMS_ANNOTATIONS[ "interface_pseudogenes_gtf" ] ),
                        os.path.join( PARAMS[ "annotations_dir" ],
                                      PARAMS_ANNOTATIONS[ "interface_numts_gtf" ] ) ),
            "_filtered_2.gtf.gz" )
def filterAgainstBiotypes( infiles, outfile ):
    """
    Is a sanity check step... uses bedtools intersect (-s) to see if there are
    transcripts in lncRNA gtf that overlap gene features specified
    It also currently removes things that intersect pseudogenes, unless they 
    also intersect lincRNA or antisense transcripts.
    """

    lncRNA_gtf, reference_gtf, pseudogene_gtf, numt_gtf  = infiles
    rejected = P.snip( outfile, "_filtered_2.gtf.gz" ) + "_rejected_2.gtf.gz"
    tmpf_1_name  = ( "./filter_me_lncrna/biotypes.gtf.gz" )
    tmpf_1 = IOTools.openFile( tmpf_1_name, "w" )

    E.info( "Biotypes to be filtered against have been"
            " written to %s" % tmpf_1_name )

    # pull out gtf entries matching the specified biotypes
    biotypes = PARAMS[ "filter_intersect_vs" ].split( "," )
    for gtf in GTF.iterator( IOTools.openFile( reference_gtf ) ):
        if gtf.source in biotypes:
            tmpf_1.write( str( gtf ) + "\n" )
        else: continue
    tmpf_1.close()

    # concatenate biotypes with numts and pseudogenes
    statement = ( "zcat %(pseudogene_gtf)s %(numt_gtf)s"
                  " >> %(tmpf_1_name)s" )
    to_cluster=False
    P.run()

    # pass the biotypes to rescue - any filtered transcript that intersect an 
    # ensembl interval with these annotations is rescued.
    control_biotypes = PARAMS[ "filter_rescue_biotypes" ].split(",") 

    # run intersect gtfs, but rescuing anything that intersects with antisense 
    # or lincRNA annotated transcripts in the ENSEMBL reference geneset.
    P10.intersectGTFs( lncRNA_gtf, 
                       tmpf_1_name, 
                       outfile, 
                       rejected, 
                       control_gtf = reference_gtf, 
                       control_biotypes = control_biotypes )

    to_cluster=False
    statement = ( "printf '\\n %(tmpf_1_name)s contains the biotypes for which"
                  " overlapping lncRNAs were removed following"
                  " filterAgainstBiotypes' >> ./filter_me_lncrna/README" )
    P.run()

#    os.unlink( tmpf_1_name )

#################################################################################
## subsection: filter lncRNA that intersect specified RefSeq reference biotypes
#################################################################################
# This removes any lncRNAs that intersect refseq mRNA transcripts. 
# There is no such thing as a gene_id in refseq, only transcript ids. 
# Transcripts with accession prefix NM_ are mRNA, while those with prefix NR_ are
# ncRNA (see http://www.ncbi.nlm.nih.gov/books/NBK21091/)
# Currently filtered against anything with NM prefix in refseq refGene table 
# downloaded from UCSC tableBrowser. 
# All processing of the file is done in the pipeline

@transform( filterAgainstBiotypes, 
            suffix( "_filtered_2.gtf.gz" ), 
            add_inputs( PARAMS[ "filter_refseq_gtf" ], 
                        os.path.join( PARAMS["location_transcriptfiles"],
                                      "reference.gtf.gz" ) ),
            "_filtered_3.gtf.gz" )
def filterAgainstRefseqTranscripts( infiles, outfile ):
    lncRNA_gtf, refseq_gtf, control_gtf = infiles
    rejected = P.snip( outfile, "_filtered_3.gtf.gz" ) + "_rejected_3.gtf.gz"
    refseq_filtered_gtf = P.getTempFile( "." )
    refseq_sorted_gtf = P.getTempFilename( "." )

    # filter refseq gtf file so that it contains only mRNA exons
    for gtf in GTF.iterator( IOTools.openFile( refseq_gtf ) ):
        if gtf.feature == "exon":
            if re.match( "NM_", gtf.gene_id ):
                refseq_filtered_gtf.write( str(gtf) + "\n" )
            else: continue
        else: continue
    refseq_filtered_gtf.close()
    refseq_filtered_gtf = refseq_filtered_gtf.name

    statement = ( "cat %(refseq_filtered_gtf)s |"
                  " python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=sort --sort-order=transcript"
                  "  --log=%(outfile)s.log"
                  " > %(refseq_sorted_gtf)s" )
    P.run()

    # pass the 
    control_biotypes = PARAMS[ "filter_rescue_biotypes" ].split(",") 

    P10.intersectGTFs( lncRNA_gtf, 
                       refseq_filtered_gtf, 
                       outfile, 
                       rejected, 
                       control_gtf, 
                       control_biotypes )

    os.unlink( refseq_filtered_gtf )
    os.unlink( refseq_sorted_gtf )

#################################################################################
## subsection: filter lncRNA in hg19 UTRs following liftover #NOT IMPLEMENTED#
#################################################################################
# This filter step is currently missed out of the pipeline because the liftOver 
# of the UTRs between hg19 and mm10 is really poor... very few UTRs lifted over. 
@transform( filterAgainstBiotypes,  
            suffix( "_filtered_2_gtf.gz" ), 
            "_filtered_3_gtf.gz" )
def filterAgainstHg19UTR( infile, outfile ):
    tempdir = P.getTempDir( "./filter_me_lncrna" )
    cds_gtf= PARAMS[ "hg19_cds" ]
    exon_gtf= PARAMS[ "hg19_exon" ]
    utr_gtf= P.getTempFilename( tempdir )
    cds_bed = P.getTempFilename( tempdir )
    exon_bed = P.getTempFilename( tempdir )
    utr_bed = P.getTempFilename( tempdir )
    chain_file = PARAMS[ "hg19_chain" ]
    utr_mm10_bed = P.getTempFilename( tempdir )

    liftOver_failed = P.snip( outfile, "gtf.gz" ) + "liftover_failed.bed"
    rejected = P.snip( outfile, "_filtered_3.gtf.gz" ) + "_rejected_3.gtf.gz"
    utr_mm10_gtf = P.snip( outfile, "_filtered_3.gtf.gz" ) + "_hg19utr_3.gtf.gz"

    statement =  ( "zcat %(cds_gtf)s |" 
                   "  python %(scriptsdir)s/gtf2gtf.py"
                   "   --method=sort --sort-order=gene"
                   "   --log=%(outfile)s.log |"
                   "  python %(scriptsdir)s/gtf2gtf.py"
                   "   --method=merge-transcripts"
                   "   --log=%(outfile)s.log |"
                   "  python %(scriptsdir)s/gff2bed.py"
                   "   --is-gtf"
                   "   --log=%(outfile)s.log"
                   "  > %(cds_bed)s;"                   
                   " zcat %(exon_gtf)s" 
                   "  python %(scriptsdir)s/gtf2gtf.py"
                   "   --method=sort --sort-order=gene"
                   "   --log=%(outfile)s.log |"
                   "  python %(scriptsdir)s/gtf2gtf.py"
                   "   --method=merge-transcripts"
                   "   --log=%(outfile)s.log |"
                   "  python %(scriptsdir)s/gff2bed.py"
                   "   --is-gtf"
                   "   --log=%(outfile)s.log"
                   "  > %(exon_bed)s; "
                   " bedtools subtract"
                   "   -a %(exon_bed)s"
                   "   -b %(cds_bed)s"
                   "   -s |"
                   "  bedtools sort"
                   "  > %(utr_bed)s;"
                   " liftOver"
                   "   %(utr_bed)s"
                   "   %(chain_file)s"
                   "   %(utr_mm10_bed)s"
                   "   %(liftover_failed)s;"
                   " zcat %(utr_mm10_bed)s"
                   "  python %(scriptsdir)s/bed2gff.py"
                   "   --as-gtf"
                   "   --log%(outfile)s.log |"
                   "  python %(scriptsdir/gtf2gtf.py)s"
                   "   --method=sort --sort-order=gene"
                   "   --log=%(outfile)s.log"
                   "  > utr_mm10_gtf" )
    P.run()

    P10.intersectGTFs( infile, utr_mm10_gtf, outfile, rejected )

    shutil.rmtree( os.path.abspath( tempdir ) )

#################################################################################
## subsection: filter on the basis of upstream intergenic coverage
#################################################################################
# need to check that co-ordinates are correct when switching bed/gtf co-ordinates
@follows( mkdir( "./filter_me_lncrna/me_intergenic_coverage" ) )
@split( filterAgainstRefseqTranscripts,  
        regex( "(.+)/(.+)_filtered_3.gtf.gz" ), 
        add_inputs( os.path.join( PARAMS[ "annotations_dir" ], 
                                  PARAMS_ANNOTATIONS[ "interface_tts_gene_bed" ] ) ),
        [ r"\1/\2_filtered_4.gtf.gz", 
          r"\1/me_intergenic_coverage/*coverage.bed.gz" ] )

def filterAgainstIntergenicCoverage( infiles, outfiles ):
    """
    """

    # infiles...
    lncRNA_gtf, tts_bed = infiles

    # outfile
    outfile = outfiles[0]
    # bamfiles
    bamfiles = []
    for entry in os.listdir( PARAMS[ "location_bamfiles" ] ):
        if entry.endswith( ".star.bam" ):
            bamfiles.append( os.path.join( PARAMS[ "location_bamfiles" ], entry ) )

    for bamfile in bamfiles:
        assert os.path.exists( bamfile + ".bai" ), "bamfile missing index %s" % bamfile

    # parameters
    window = int( PARAMS[ "filter_downstream" ] * 1000 )
    max_coverage = float( PARAMS[ "filter_coverage" ] )
    max_occurrence = int( PARAMS[ "filter_threshold" ] )

    # outfiles
    filtered_gtf = outfile
    rejected_gtf = P.snip( outfile, "_filtered_4.gtf.gz" ) + "_rejected_4.gtf.gz"
    rejected_tsv = P.snip( outfile, "_filtered_4.gtf.gz" ) + "_rejected_4.tsv.gz"

    P10.filterOnIntergenicCoverage( lncRNA_gtf,
                                    tts_bed, 
                                    bamfiles,
                                    window,
                                    max_coverage,
                                    max_occurrence, 
                                    filtered_gtf, 
                                    rejected_gtf, 
                                    rejected_tsv, 
                                    outdir = "./filter_me_lncrna/me_intergenic_coverage" )

    
@transform( filterAgainstIntergenicCoverage, 
            suffix( "_coverage.bed.gz" ), 
            "_coverage_vs_dist.tsv.gz" )
def calcCoverageVSDist( infile, outfile ):
    outf = IOTools.openFile( outfile, "w" )
    for line in IOTools.openFile( infile ).readlines():
        line = line.split()
        lnc_id = line[3].split("__")[0]
        dist = str( int(line[2]) - int(line[1]) )
        outf.write( "\t".join( [ lnc_id, dist, line[4] ] ) + "\n" )
    outf.close()


@merge( calcCoverageVSDist, 
        "./filter_me_lncrna/me_lncrna_intergenic_coverage.tsv.gz" )
def combineCoverageVSDist( infiles, outfile ):
    outf = P.snip( outfile, ".gz" )
    file_names = " ".join( [ x for x in infiles ] ) 

    # combine tables will not write headers when --no-titles is specified
    # writing headers to file first
    headers = [ P.snip( os.path.basename(x), 
                        ".star_coverage_vs_dist.tsv.gz" ) for x in infiles ]
    headers =  "lnc_id\tdistance\t" + "\t".join( headers )
    out = IOTools.openFile( outf, "w" )
    out.write( headers + "\n" )
    out.close()

    # using sed to split the join columns
    to_cluster = False
    statement = ( "python %(scriptsdir)s/combine_tables.py"
                  "  --columns=1,2"
                  "  --no-titles"
#                  " --header-names=%(headers)s"
                  "  --log=%(outfile)s.log"
                  "  %(file_names)s |"
                  " sed 's/-/\\t/'" 
                  " >> %(outf)s;"
                  " gzip %(outf)s" )
    P.run()


@transform( combineCoverageVSDist, regex( "(.+)/(.+).tsv.gz" ), r"\2.load" )
def loadCoverageVSDist( infile, outfile ):
    P.load( infile, outfile )

            
@follows( calcCoverageVSDist )
@transform( filterAgainstIntergenicCoverage, 
            suffix( "_coverage.bed.gz" ), 
            "_coverage.hist.gz" ) 
def calcCoverageHist( infile, outfile ):
    window = str( int( PARAMS[ "filter_downstream" ] ) * 1000 )
    header = P.snip( os.path.basename( infile ), "_coverage.bed.gz" )

    statement = ( "zcat %(infile)s |"
                  " awk '$3-$2 <= %(window)s {print $5}' |"
                  " python %(scriptsdir)s/data2histogram.py"
                  "  --header-names=%(header)s"
                  "  --range=0,1"
                  "  --bin-size=0.01"
                  "  --log=%(outfile)s.log |"
                  " gzip > %(outfile)s" )
    P.run()


@merge( calcCoverageHist, 
        "./filter_me_lncrna/me_lncrna_intergenic_coverage.hist.gz" )
def combineCoverageHist( infiles, outfile ):
    file_names = " ".join( [ x for x in infiles ] ) 
    to_cluster = False
    statement = ( "python %(scriptsdir)s/combine_tables.py" 
                  "  --columns=1"
                  "  --log=%(outfile)s.log"
                  "  %(file_names)s |"
                  " gzip"
                  " > %(outfile)s" )
    P.run()


@transform( combineCoverageHist, regex( "(.+)/(.+).gz" ), r"\2.load" )
def loadCoverageHist( infile, outfile ):
    P.load( infile, outfile )

#################################################################################
## subsection: filter on the basis of shared splice junctions
#################################################################################
@transform( filterAgainstIntergenicCoverage,
            suffix( "_filtered_4.gtf.gz" ), 
            add_inputs( os.path.join( PARAMS[ "location_transcriptfiles" ], 
                                      "reference.gtf.gz" ) ),
            "_filtered_5.gtf.gz" )
def filterOnSharedSpliceJunctions( infiles, outfile ):

    # infiles
    lncRNA_gtf, reference_gtf = infiles
    
    # preprocess reference file
    reference_biotypes = PARAMS[ "filter_classify_vs" ].split( "," )
    tmpf_ref = P.getTempFile( "./filter_me_lncrna" )
    tmpf_ref_name = tmpf_ref.name
    for gtf in GTF.iterator( IOTools.openFile( reference_gtf ) ):
        if gtf.source in reference_biotypes:
            if gtf.feature == "exon":
                tmpf_ref.write( str( gtf ) + "\n" )
            else: continue
        else: continue
    tmpf_ref.close()

    # creating a string of bamfiles 
    bamfiles = []
    for entry in os.listdir( PARAMS[ "location_bamfiles" ] ):
        if entry.endswith( ".star.bam" ):
            bamfiles.append( os.path.join( PARAMS[ "location_bamfiles" ], entry ) )
    bamfiles = "__".join( bamfiles )

    # parameters
    max_occurrence = str( PARAMS[ "filter_sj_threshold" ] )

    # outfiles
    filtered_gtf = outfile
    rejected_gtf = P.snip( outfile, "_filtered_5.gtf.gz" ) + "_rejected_5.gtf.gz"
    rejected_tsv = P.snip( outfile, "_filtered_5.gtf.gz" ) + "_rejected_5.tsv.gz"

    # submitting job to cluster
    params = [ lncRNA_gtf, 
               tmpf_ref_name, 
               bamfiles, 
               max_occurrence, 
               filtered_gtf, 
               rejected_gtf, 
               rejected_tsv ]
    P.submit( "/ifs/devel/projects/proj010/PipelineProj010", 
              "filterOnSharedSplicedReads", 
              params, 
              jobOptions = " -l mem_free=25G" )

#################################################################################
## subsection: filter unstranded gene models
#################################################################################
@transform( filterOnSharedSpliceJunctions, 
            suffix( "_filtered_5.gtf.gz" ),
            "_filtered_6.gtf.gz" )
def filterMEUnstranded( infile, outfile ):
    lncRNA_filtered = IOTools.openFile( outfile, "w" )
    lncRNA_rejected = P.snip( outfile,  "filtered_6.gtf.gz" ) + "rejected_6.gtf.gz" 
    lncRNA_rejected = IOTools.openFile( lncRNA_rejected, "w" )

    tmp_gtf = P.getTempFilename( "." )
    statement = ( "zcat %(infile)s |" 
                  "python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=sort --sort-order=gene+transcript"
                  "  --log=%(outfile)s.log"
                  " > %(tmp_gtf)s" )
    P.run()

    for gtfs in GTF.flat_gene_iterator( GTF.iterator( IOTools.openFile( tmp_gtf ) ) ):
        exon_strands = [ x.strand for x in gtfs ] 
        if "." in exon_strands:
            outf = lncRNA_rejected
        else: 
            outf = lncRNA_filtered

        for exon in gtfs:
            outf.write( str( exon ) + "\n" )

    lncRNA_filtered.close()
    lncRNA_rejected.close()
    os.unlink( tmp_gtf )


#################################################################################
## subsection: filter anything classified as downstream of a coding gene
#################################################################################
@transform( filterMEUnstranded, 
            suffix( "_filtered_6.gtf.gz" ), 
            "_filtered_7.gtf.gz" )
def filterMESenseDownstream( infile, outfile ):
    """
    Remove anything classified as sense_downstream
    """
    inf = IOTools.openFile( infile )
    outf = IOTools.openFile( outfile, "w" )
    outf_rej = IOTools.openFile( re.sub( "filtered", "rejected", outfile ), "w" )

    for gtf in GTF.iterator( inf ):
        if gtf.source == "sense_downstream":
            outf_rej.write( str( gtf ) + "\n" )
        else:
            outf.write( str( gtf ) + "\n" )
    
    inf.close()
    outf.close()
    outf_rej.close()


@follows( filterMESenseDownstream, 
          loadCoverageHist )
def filterMELncRNA(): pass

#################################################################################
#################################################################################
#################################################################################
# Section: filter single-exon lncRNA loci
#################################################################################
# Filter anything that is within 5kb sense upstream/downstream of pc gene model
# liftover human 5' 3'UTR and remove anything within these regions
# remove anything within 20(?)kb downstream of pc gene, with >60% intergenic cov
#################################################################################
## subsection: classify se lncRNA relative to the reference geneset
#################################################################################
@transform( splitRenamedLncRNA,
            suffix( "_se_filtered_0.gtf.gz" ), 
            add_inputs( os.path.join( PARAMS["location_transcriptfiles"],
                                      "reference.gtf.gz" ) ), 
            "_se_classified.gtf.gz" )
def classifySELncRNAVSReference( infiles, outfile ):
    lncRNA_gtf, reference_gtf = infiles
    tmp_gtf = P.getTempFilename( "./filter_se_lncrna" )
    tmp_gtf = tmp_gtf + ".gz"

    to_cluster = False
    statement = ( " zcat %(lncRNA_gtf)s |"
                  " python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=sort --sort-order=gene+transcript"
                  "  --log=%(outfile)s.log |"
                  " gzip"
                  " > %(tmp_gtf)s" )
    P.run()

    dist_upstr = PARAMS[ "filter_distance_upstream" ]
    dist_dstr = PARAMS[ "filter_distance_downstream" ]
    reference_biotypes = PARAMS[ "filter_classify_vs" ].split( "," )

    # preprocess reference file
    tmpf_ref = P.getTempFile( "./filter_se_lncrna" )
    tmpf_ref_name = tmpf_ref.name
    for gtf in GTF.iterator( IOTools.openFile( reference_gtf ) ):
        if gtf.source in reference_biotypes:
            if gtf.feature == "exon": 
                tmpf_ref.write( str( gtf ) + "\n" )
            else: continue
        else: continue
    tmpf_ref.close()

    tempdir = PipelineLncRNA.reClassifyLncRNAGenes( tmp_gtf,
                                                    tmpf_ref_name, 
                                                    outfile, 
                                                    dist_upstr,
                                                    dist_dstr, 
                                                    wdir = "./filter_se_lncrna")
    shutil.rmtree( tempdir )
    os.unlink( tmpf_ref_name )
    os.unlink( tmp_gtf )


@transform( classifySELncRNAVSReference, 
            suffix( "_classified.gtf.gz" ), 
            "_filtered_1.gtf.gz" )
def filterSELncRNABasedOnClassification( infile, outfile ):
    rejected = P.snip( outfile, "_filtered_1.gtf.gz" ) + "_rejected_1.gtf.gz"
    inf = IOTools.openFile( infile )
    outf = IOTools.openFile( outfile, "w" )
    outf_reject = IOTools.openFile( rejected, "w" )
    to_reject = PARAMS[ "filter_reject_se" ].split( "," )
    
    rejected = 0
    
    for gtf in GTF.transcript_iterator( GTF.iterator( inf ) ):
        if gtf[0].source in to_reject:
            rejected += 1
            for exon in IOTools.flatten( gtf ):
                outf_reject.write( str( exon ) + "\n" )
        else:
            for exon in IOTools.flatten( gtf ):
                outf.write( str( exon ) + "\n" )

    E.info( "%i lncRNA transcripts rejected based on their"
            " classification" % rejected  )

    outf.close()
    outf_reject.close()

#################################################################################
## subsection: filter lncRNA that intersect undesired ensembl biotypes
#################################################################################
# A sanity check to make sure that all transcripts intersecting known coding loci 
# have been removed in the lncRNA pipeline.
@transform( filterSELncRNABasedOnClassification, 
            suffix( "_filtered_1.gtf.gz" ),
            add_inputs( os.path.join( PARAMS["location_transcriptfiles"],
                                      "reference.gtf.gz" ), 
                        os.path.join( PARAMS[ "annotations_dir" ], 
                                      PARAMS_ANNOTATIONS[ "interface_pseudogenes_gtf" ] ),
                        os.path.join( PARAMS[ "annotations_dir" ],
                                      PARAMS_ANNOTATIONS[ "interface_numts_gtf" ] ) ),
            "_filtered_2.gtf.gz" )
def filterSELncRNAAgainstBiotypes( infiles, outfile ):
    """
    Again, a sanity check... using bedtools intersect (-s) to remove transcripts 
    in lncRNA_gtf that overlap specified features. Unless they also intersect 
    intervals in control file annotated as lncRNA or antisense.
    """

    lncRNA_gtf, reference_gtf, pseudogene_gtf, numt_gtf  = infiles
    rejected = P.snip( outfile, "_filtered_2.gtf.gz" ) + "_rejected_2.gtf.gz"
    tmpf_1_name  = ( "./filter_se_lncrna/biotypes.gtf.gz" )
    tmpf_1 = IOTools.openFile( tmpf_1_name, "w" )

    E.info( "Biotypes to filter SE lncRNA against have been"
            " written to %s" % tmpf_1_name )

    # pull out gtf entries matching the specified biotypes
    biotypes = PARAMS[ "filter_intersect_vs" ].split( "," )
    for gtf in GTF.iterator( IOTools.openFile( reference_gtf ) ):
        if gtf.source in biotypes:
            tmpf_1.write( str( gtf ) + "\n" )
        else: continue
    tmpf_1.close()

    # concatenate biotypes with numts and pseudogenes
    statement = ( "zcat %(pseudogene_gtf)s %(numt_gtf)s"
                  " >> %(tmpf_1_name)s" )
    to_cluster=False
    P.run()

    # pass the 
    control_biotypes = PARAMS[ "filter_rescue_biotypes" ].split(",") 

    # run intersect GTFs
    P10.intersectGTFs( lncRNA_gtf, 
                       tmpf_1_name, 
                       outfile, 
                       rejected, 
                       control_gtf = reference_gtf, 
                       control_biotypes = control_biotypes )

    os.unlink( tmpf_1_name )

#################################################################################
## subsection: filter lncRNA that intersect specified RefSeq reference biotypes
#################################################################################
# This removes any lncRNAs that intersect refseq mRNA transcripts. 
# There is no such thing as a gene_id in refseq, only transcript ids. 
# Transcripts with accession prefix NM_ are mRNA, while those with prefix NR_ are
# ncRNA (see http://www.ncbi.nlm.nih.gov/books/NBK21091/)
# Currently filtered against anything with NM prefix in refseq refGene table 
# downloaded from UCSC tableBrowser. 
# All processing of the file is done in the pipeline

@transform( filterSELncRNAAgainstBiotypes, 
            suffix( "_filtered_2.gtf.gz" ), 
            add_inputs( PARAMS[ "filter_refseq_gtf" ], 
                        os.path.join( PARAMS["location_transcriptfiles"], 
                                      "reference.gtf.gz" ) ),
            "_filtered_3.gtf.gz" )
def filterSELncRNAAgainstRefseqTranscripts( infiles, outfile ):
    lncRNA_gtf, refseq_gtf, control_gtf  = infiles
    rejected = P.snip( outfile, "_filtered_3.gtf.gz" ) + "_rejected_3.gtf.gz"
    refseq_filtered_gtf = P.getTempFile( "." )
    refseq_sorted_gtf = P.getTempFilename( "." )

    # filter refseq gtf file so that it contains only mRNA exons
    for gtf in GTF.iterator( IOTools.openFile( refseq_gtf ) ):
        if gtf.feature == "exon":
            if re.match( "NM_", gtf.gene_id ):
                refseq_filtered_gtf.write( str(gtf) + "\n" )
            else: continue
        else: continue
    refseq_filtered_gtf.close()
    refseq_filtered_gtf = refseq_filtered_gtf.name

    statement = ( "cat %(refseq_filtered_gtf)s |"
                  " python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=sort --sort-order=transcript"
                  "  --log=%(outfile)s.log"
                  " > %(refseq_sorted_gtf)s" )
    P.run()

    # pass the control biotypes
    control_biotypes = PARAMS[ "filter_rescue_biotypes" ].split(",") 

    P10.intersectGTFs( lncRNA_gtf, 
                       refseq_filtered_gtf, 
                       outfile, 
                       rejected, 
                       control_gtf, 
                       control_biotypes )

    os.unlink( refseq_filtered_gtf )
    os.unlink( refseq_sorted_gtf )

#################################################################################
## subsection: filter lncRNA in hg19 UTRs following liftover # NOT IMPLEMENTED#
#################################################################################
# This filter step is currently missed out of the pipeline because the liftOver 
# of the UTRs between hg19 and mm10 is really poor... very few UTRs lifted over. 
@transform( filterSELncRNAAgainstBiotypes,  
            suffix( "_filtered_2_gtf.gz" ), 
            "_filtered_3_gtf.gz" )
def filterSELncRNAAgainstHg19UTR( infile, outfile ):
    tempdir = P.getTempDir( "./filter_me_lncrna" )
    cds_gtf= PARAMS[ "hg19_cds" ]
    exon_gtf= PARAMS[ "hg19_exon" ]
    utr_gtf= P.getTempFilename( tempdir )
    cds_bed = P.getTempFilename( tempdir )
    exon_bed = P.getTempFilename( tempdir )
    utr_bed = P.getTempFilename( tempdir )
    chain_file = PARAMS[ "hg19_chain" ]
    utr_mm10_bed = P.getTempFilename( tempdir )

    liftOver_failed = P.snip( outfile, "gtf.gz" ) + "liftover_failed.bed"
    rejected = P.snip( outfile, "_filtered_3.gtf.gz" ) + "_rejected_3.gtf.gz"
    utr_mm10_gtf = P.snip( outfile, "_filtered_3.gtf.gz" ) + "_hg19utr_3.gtf.gz"

    statement =  ( "zcat %(cds_gtf)s |" 
                   "  python %(scriptsdir)s/gtf2gtf.py"
                   "   --method=sort --sort-order=gene"
                   "   --log=%(outfile)s.log |"
                   "  python %(scriptsdir)s/gtf2gtf.py"
                   "   --method=merge-transcripts"
                   "   --log=%(outfile)s.log |"
                   "  python %(scriptsdir)s/gff2bed.py"
                   "   --is-gtf"
                   "   --log=%(outfile)s.log"
                   "  > %(cds_bed)s;"                   
                   " zcat %(exon_gtf)s" 
                   "  python %(scriptsdir)s/gtf2gtf.py"
                   "   --method=sort --sort-order=gene"
                   "   --log=%(outfile)s.log |"
                   "  python %(scriptsdir)s/gtf2gtf.py"
                   "   --method=merge-transcripts"
                   "   --log=%(outfile)s.log |"
                   "  python %(scriptsdir)s/gff2bed.py"
                   "   --is-gtf"
                   "   --log=%(outfile)s.log"
                   "  > %(exon_bed)s; "
                   " bedtools subtract"
                   "   -a %(exon_bed)s"
                   "   -b %(cds_bed)s"
                   "   -s |"
                   "  bedtools sort"
                   "  > %(utr_bed)s;"
                   " liftOver"
                   "   %(utr_bed)s"
                   "   %(chain_file)s"
                   "   %(utr_mm10_bed)s"
                   "   %(liftover_failed)s;"
                   " zcat %(utr_mm10_bed)s"
                   "  python %(scriptsdir)s/bed2gff.py"
                   "   --as-gtf"
                   "   --log%(outfile)s.log |"
                   "  python %(scriptsdir/gtf2gtf.py)s"
                   "   --method=sort --sort-order=gene"
                   "   --log=%(outfile)s.log"
                   "  > utr_mm10_gtf" )
    P.run()

    P10.intersectGTFs( infile, utr_mm10_gtf, outfile, rejected )

    shutil.rmtree( os.path.abspath( tempdir ) )

#################################################################################
## subsection: filter on the basis of upstream intergenic coverage
#################################################################################
# need to check that co-ordinates are correct when switching bed/gtf co-ordinates
@follows( mkdir( "./filter_se_lncrna/se_intergenic_coverage" ) )
@split( filterSELncRNAAgainstRefseqTranscripts,  
        regex( "(.+)/(.+)_filtered_3.gtf.gz" ), 
            add_inputs( os.path.join( PARAMS[ "annotations_dir" ], 
                                      PARAMS_ANNOTATIONS[ "interface_tts_gene_bed" ] ) ),
        [ r"\1/\2_filtered_4.gtf.gz",  
          r"\1/se_intergenic_coverage/*coverage.bed.gz" ] )
def filterSELncRNAAgainstIntergenicCoverage( infiles, outfiles ):

    # infiles...
    lncRNA_gtf, tts_bed = infiles
    #outfile
    outfile = outfiles[0]
    # bamfiles
    bamfiles = []
    for entry in os.listdir( PARAMS[ "location_bamfiles" ] ):
        if entry.endswith( ".star.bam" ):
            bamfiles.append( os.path.join( PARAMS[ "location_bamfiles" ], entry ) )

    # parameters
    window = int( PARAMS[ "filter_downstream" ] * 1000 )
    max_coverage = float( PARAMS[ "filter_coverage_se" ] )
    max_occurrence = int( PARAMS[ "filter_threshold_se" ] )

    # outfiles
    filtered_gtf = outfile
    rejected_gtf = P.snip( outfile, "_filtered_4.gtf.gz" ) + "_rejected_4.gtf.gz"
    rejected_tsv = P.snip( outfile, "_filtered_4.gtf.gz" ) + "_rejected_4.tsv.gz"

    P10.filterOnIntergenicCoverage( lncRNA_gtf,
                                    tts_bed, 
                                    bamfiles,
                                    window,
                                    max_coverage,
                                    max_occurrence, 
                                    filtered_gtf, 
                                    rejected_gtf, 
                                    rejected_tsv, 
                                    outdir = "./filter_se_lncrna/se_intergenic_coverage" )


@transform( filterSELncRNAAgainstIntergenicCoverage, 
            suffix( "_coverage.bed.gz" ), 
            "_coverage_vs_dist.tsv.gz" )
def calcSELncRNACoverageVSDist( infile, outfile ):
    outf = IOTools.openFile( outfile, "w" )
    for line in IOTools.openFile( infile ).readlines():
        line = line.split()
        lnc_id = line[3].split("__")[0]
        dist = str( int(line[2]) - int(line[1]) )
        outf.write( "\t".join( [ lnc_id, dist, line[4] ] ) + "\n" )
    outf.close()


@merge( calcCoverageVSDist, 
        "./filter_se_lncrna/se_lncrna_intergenic_coverage.tsv.gz" )
def combineSELncRNACoverageVSDist( infiles, outfile ):
    outf = P.snip( outfile, ".gz" )
    file_names = " ".join( [ x for x in infiles ] ) 

    # combine tables will not write headers when --no-titles is specified
    # writing headers to file first
    headers = [ P.snip( os.path.basename(x), 
                        ".star_coverage_vs_dist.tsv.gz" ) for x in infiles ]
    headers =  "lnc_id\tdistance\t" + "\t".join( headers )
    out = IOTools.openFile( outf, "w" )
    out.write( headers + "\n" )
    out.close()

    # using sed to split the join columns
    to_cluster = False
    statement = ( "python %(scriptsdir)s/combine_tables.py"
                  "  --columns=1,2"
                  "  --no-titles"
#                  " --header-names=%(headers)s"
                  "  --log=%(outfile)s.log"
                  "  %(file_names)s |"
                  " sed 's/-/\\t/'" 
                  " >> %(outf)s;"
                  " gzip %(outf)s" )
    P.run()


@transform( combineSELncRNACoverageVSDist, regex( "(.+)/(.+).tsv.gz" ), r"\2.load" )
def loadSELncRNACoverageVSDist( infile, outfile ):
    P.load( infile, outfile )


@follows( calcSELncRNACoverageVSDist )
@transform( filterSELncRNAAgainstIntergenicCoverage, 
            suffix( "_coverage.bed.gz" ), 
            "_coverage.hist.gz" ) 
def calcSELncRNACoverageHist( infile, outfile ):
    window = str( int( PARAMS[ "filter_downstream" ] ) * 1000 )
    header = P.snip( os.path.basename( infile ), "_coverage.bed.gz" )

    statement = ( "zcat %(infile)s |"
                  " awk '$3-$2 <= %(window)s {print $5}' |"
                  " python %(scriptsdir)s/data2histogram.py"
                  "  --header-names=%(header)s"
                  "  --range=0,1"
                  "  --bin-size=0.01"
                  "  --log=%(outfile)s.log |"
                  " gzip > %(outfile)s" )
    P.run()


@merge( calcSELncRNACoverageHist, 
        "./filter_se_lncrna/se_lncrna_intergenic_coverage.hist.gz" )
def combineSELncRNACoverageHist( infiles, outfile ):
    file_names = " ".join( [ x for x in infiles ] ) 
    to_cluster = False
    statement = ( "python %(scriptsdir)s/combine_tables.py" 
                  "  --columns=1"
                  "  --log=%(outfile)s.log"
                  "  %(file_names)s |"
                  " gzip"
                  " > %(outfile)s" )
    P.run()


@transform( combineSELncRNACoverageHist, regex( "(.+)/(.+).gz" ), r"\2.load" )
def loadSELncRNACoverageHist( infile, outfile ):
    P.load( infile, outfile )

#################################################################################
## subsection: filter on the basis of shared splice junctions
#################################################################################
@transform( filterSELncRNAAgainstIntergenicCoverage,
            suffix( "_filtered_4.gtf.gz" ), 
            add_inputs( os.path.join( PARAMS[ "location_transcriptfiles" ], 
                                      "reference.gtf.gz" ) ),
            "_filtered_5.gtf.gz" )
def filterSELncRNAOnSharedSpliceJunctions( infiles, outfile ):

    # infiles
    lncRNA_gtf, reference_gtf = infiles
    
    # preprocess reference file
    reference_biotypes = PARAMS[ "filter_classify_vs" ].split( "," )
    tmpf_ref = P.getTempFile( "./filter_me_lncrna" )
    tmpf_ref_name = tmpf_ref.name
    for gtf in GTF.iterator( IOTools.openFile( reference_gtf ) ):
        if gtf.source in reference_biotypes:
            if gtf.feature == "exon":
                tmpf_ref.write( str( gtf ) + "\n" )
            else: continue
        else: continue
    tmpf_ref.close()

    # creating a string of bamfiles that can be passed via P.submit
    bamfiles = []
    for entry in os.listdir( PARAMS[ "location_bamfiles" ] ):
        if entry.endswith( ".star.bam" ):
            bamfiles.append( os.path.join( PARAMS[ "location_bamfiles" ], entry ) )
    bamfiles = "__".join( bamfiles )

    # parameters
    max_occurrence = str( PARAMS[ "filter_sj_threshold_se" ] )

    # outfiles
    filtered_gtf = outfile
    rejected_gtf = P.snip( outfile, "_filtered_5.gtf.gz" ) + "_rejected_5.gtf.gz"
    rejected_tsv = P.snip( outfile, "_filtered_5.gtf.gz" ) + "_rejected_5.tsv.gz"

    # submitting job to cluster
    params = [ lncRNA_gtf, 
               tmpf_ref_name, 
               bamfiles, 
               max_occurrence, 
               filtered_gtf, 
               rejected_gtf, 
               rejected_tsv ]
    P.submit( "/ifs/devel/projects/proj010/PipelineProj010", 
              "filterOnSharedSplicedReads", 
              params, 
              jobOptions = " -l mem_free=25G" )

#################################################################################
## subsection: filter unstranded gene models
#################################################################################
@transform( filterSELncRNAOnSharedSpliceJunctions, 
            suffix( "_filtered_5.gtf.gz" ),
            "_filtered_6.gtf.gz" )
def filterSEUnstranded( infile, outfile ):
    """
    A small proportion c.19 selncRNA have strand recorded as "."
    These are filtered out at this point because they screw up subsequent
    pipeline steps. An alternative would be to split duplicate these loci,
    recording a separate transcript on the forward and reverse strand. 
    """
    lncRNA_filtered = IOTools.openFile( outfile, "w" )
    lncRNA_rejected = P.snip( outfile,  "filtered_6.gtf.gz" ) + "rejected_6.gtf.gz" 
    lncRNA_rejected = IOTools.openFile( lncRNA_rejected, "w" )

    tmp_gtf = P.getTempFilename( "." )
    statement = ( "zcat %(infile)s |" 
                  "python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=sort --sort-order=gene+transcript"
                  "  --log=%(outfile)s.log"
                  " > %(tmp_gtf)s" )
    P.run()

    for gtfs in GTF.flat_gene_iterator( GTF.iterator( IOTools.openFile( tmp_gtf ) ) ):
        exon_strands = [ x.strand for x in gtfs ] 
        if "." in exon_strands:
            outf = lncRNA_rejected
        else: 
            outf = lncRNA_filtered

        for exon in gtfs:
            outf.write( str( exon ) + "\n" )

    lncRNA_filtered.close()
    lncRNA_rejected.close()
    os.unlink( tmp_gtf )

#################################################################################
## subsection: filter anything classified as downstream of a coding gene
#################################################################################
@transform( filterSEUnstranded, 
            suffix( "_filtered_6.gtf.gz" ), 
            "_filtered_7.gtf.gz" )
def filterSESenseDownstream( infile, outfile ):
    """
    Remove anything classified as sense_downstream
    """
    inf = IOTools.openFile( infile )
    outf = IOTools.openFile( outfile, "w" )
    outf_rej = IOTools.openFile( re.sub( "filtered", "rejected", outfile ), "w" )

    for gtf in GTF.iterator( inf ):
        if gtf.source == "sense_downstream":
            outf_rej.write( str( gtf ) + "\n" )
        else:
            outf.write( str( gtf ) + "\n" )
    
    inf.close()
    outf.close()
    outf_rej.close()

#################################################################################
## subsection: classify se lncRNA relative to the me lncRNA geneset
#################################################################################
@transform( filterSESenseDownstream, 
            suffix( "_filtered_7.gtf.gz" ),
            add_inputs( filterMESenseDownstream ), 
            "_classified_vs_melncRNA.gtf.gz" )
def classifySElncRNAVSMElncRNA( infiles, outfile ):
    se_lnc_gtf, me_lnc_gtf = infiles
    tmpf = P.getTempFilename( "./filter_se_lncrna" )
    tmpf_log = tmpf + ".log"

    # specify distance at which selncRNA are classified as upstr/dstr of melncNRA
    dist_upstr = PARAMS[ "filter_distance_upstream" ]
    dist_dstr = PARAMS[ "filter_distance_downstream" ]

    tempdir = PipelineLncRNA.reClassifyLncRNAGenes( se_lnc_gtf, 
                                                    me_lnc_gtf, 
                                                    tmpf, 
                                                    dist_upstr, 
                                                    dist_dstr,
                                                    wdir = "./filter_se_lncrna" )
    statement = ( "cat %(tmpf)s |" 
                  " python %(scriptsdir)s/gtf2gtf.py"
                  "   --method=sort --sort-order=gene+transcript"
                  "   --log=%(outfile)s.log |"
                  " gzip > %(outfile)s" )
    P.run()

    shutil.rmtree( tempdir )
    os.unlink( tmpf )
    os.unlink( tmpf_log )


@transform( filterSESenseDownstream,
            suffix( "_filtered_7.gtf.gz" ), 
            add_inputs( classifySElncRNAVSMElncRNA ),
            "_filtered_8.gtf.gz" )
def filterSElncRNABasedOnClassificationVSMElncRNA( infiles, outfile ):
    se_lnc_gtf, se_lnc_class = infiles
    outf = IOTools.openFile( outfile, "w" )
    outf_rej = IOTools.openFile( re.sub( "filtered", "rejected", outfile ), "w" )
    
    # select se lncRNA classifications to be rejected
    rejected_class = PARAMS[ "filter_reject_se_vs_me" ].split( "," )
    to_reject = []
    rejected = 0

    # iterate through gtf containing se lncRNA classifications,
    # pull out those that are to be rejected.
    for gtfs in GTF.flat_gene_iterator( GTF.iterator( IOTools.openFile( se_lnc_class ) ) ):
        if len( gtfs ) > 1:
            assert re.search( "sm", gtfs[0].gene_id ), \
            "There are unmerged se lncRNA with the same gene_id"
        if gtfs[0].source in rejected_class:
            to_reject.append( gtfs[0].gene_id )
        else:
            continue

    # filter se lncRNAs
    for gtf in GTF.iterator( IOTools.openFile( se_lnc_gtf ) ):
        if gtf.gene_id in to_reject:
            rejected += 1
            outf_rej.write( str( gtf ) + "\n" )
        else:
            outf.write( str( gtf ) + "\n" )

    E.info( "%i lncRNA transcripts rejected based on their"
            " classification" % rejected  )

    outf.close()
    outf_rej.close()


#          loadSELncRNACoverageHist )
@follows( filterSElncRNABasedOnClassificationVSMElncRNA ) 
def filterSELncRNA(): pass

#################################################################################
#################################################################################
#################################################################################
# Section: merge multi-exon and single-exon lncRNAs
#################################################################################
#################################################################################
## subsection: merge multi-exon and single-exon lncRNAs
#################################################################################
@follows( mkdir( "./filter_lncrna" ) )
@merge( [ filterMESenseDownstream, 
          filterSElncRNABasedOnClassificationVSMElncRNA ], 
        "./filter_lncrna/lncRNA_combined.gtf.gz" )
def combineMEandSELncRNA( infiles, outfile ):
    me_lncRNA, se_lncRNA = infiles
    tmpf = P.getTempFilename( "." )

    to_cluster=False
    statement = ( "zcat"
                  "  %(me_lncRNA)s"
                  "  %(se_lncRNA)s"
                  " > %(tmpf)s;"
                  " cat %(tmpf)s |"
                  " python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=sort --sort-order=gene+transcript"
                  "  --log=%(outfile)s.log |"
                  " gzip > %(outfile)s" )
    P.run()

    os.unlink( tmpf ) 

@jobs_limit( 1 )
@transform( [ combineMEandSELncRNA, 
              os.path.join( PARAMS[ "location_transcriptfiles" ], 
                            "refcoding.gtf.gz" ) ], 
            regex( "(.+)/(.+).gtf.gz" ), 
            r"./filter_lncrna/\2_intergenic_dist.tsv.gz" )
def calcIntergenicDistances( infile, outfile ):
    """
    Calculate the distribution of intergenic distances between adjacent gene 
    models.
    """
    inf = IOTools.openFile( infile )
    outf = IOTools.openFile( outfile, "w" )

    # create a generator of sorted chunks, and calc distance between chunks
    gtf_chunks = GTF.iterator_sorted_chunks( GTF.flat_gene_iterator( GTF.iterator( inf ) ), 
                                             sort_by = "contig-strand-start-end" )
    distances_dict, contained = P10.calcIntergenicDist( gtf_chunks, 
                                                        ignore_overlap = True )

    # write distances to outfiles
    distances = distances_dict.values()
    for i in distances: 
        outf.write( str(i) + "\n" )
    outf.close()

    # plot log densiy plot of distance
    outf_den = P.snip( outfile, ".tsv.gz" ) + "density.png"
    distances = robjects.IntVector( distances )
    R.assign( "distances", distances )
    R( '''png( "%s" )''' % outf_den )
    R( '''plot( density( log( distances), na.rm=TRUE ) )''' )
    R( '''dev.off()''' )
    
    # plot histogram of 10th percentile
    outf_hist = P.snip( outfile, ".tsv.gz" ) + "p10_hist.png"
    dist10 = [ x for x in distances 
               if x < stats.scoreatpercentile( distances, per = 10 ) ]
    dist10 = robjects.IntVector( dist10 )
    R.assign( "dist10", dist10 )
    R( '''png( "%s" )''' % outf_hist )
    R( '''hist(dist10, breaks=1000)''' )
    R( '''dev.off()''' )


# debug:
# @follows( mkdir( "filter_lncrna_debug" ) )
# @transform( combineMEandSELncRNA, 
#             regex( "(.+)/(.+)_combined.gtf.gz" ),
#             add_inputs(filterMESenseDownstream ),
#             r"filter_lncrna_debug/\2_merged.gtf.gz" )
@transform( combineMEandSELncRNA, 
           suffix( "_combined.gtf.gz" ),
           add_inputs( filterMESenseDownstream ), 
           "_merged.gtf.gz" )
def mergeMEandSELncRNA( infiles, outfile ):
    """
    On the assumption that se lncRNA located close to me lncRNA are unresolved
    exons or extensions, rather than independent loci:
    Merges all gene models that are less than a specified distance from one
    another on the same strand. For single-exon genes, the original transcript is
    removed from the outfile. For multi-exon genes the gene_id is renamed (ms),
    the merged transcript is added to the outfile and all existing transcripts
    are also appended to the outfile.
    """

    all_gtf, me_gtf = infiles

    # create gtf of introns for multi-exon lncRNAs
    me_introns_gtf = P.getTempFilename( "./filter_lncrna" )
    to_cluster=False
    statement = ( "zcat %(me_gtf)s |"
                  " python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=exons2introns"
                  "  --log=%(outfile)s.log"
                  " > %(me_introns_gtf)s" )
    P.run()

    # plot histogram of intron size distribution (legacy option)
    out_hist = P.snip( outfile, "_merged.gtf.gz" ) + "_me_intronSize.hist.png"
    me_introns = []
    for gtf in GTF.iterator( IOTools.openFile( me_introns_gtf ) ):
        me_introns.append( gtf.end - gtf.start )
    intron_dist = stats.scoreatpercentile( me_introns, per = 50 )
    intron_vec = robjects.IntVector( me_introns ) 
    R.assign( "intron.vec", intron_vec  )
    per = int( intron_dist )
    R.assign( "intron.per", per )
    R( '''png( "%(out_hist)s" )''' % locals() )
    R( '''hist( intron.vec, breaks=1000, xlim=c(0,10000) )''' )
    R( '''abline( v=intron.per, col="red" )''' )
    R( '''dev.off()''' )
  
    # ascertain the number of gene models in the infile
    n_genes_in = 0
    for gtf in GTF.flat_gene_iterator( GTF.iterator( IOTools.openFile( all_gtf) ) ):
        n_genes_in += 1

    # create a generator object that yields list of exons with the same gene_id,
    # sorted by contig, strand, start
    gtf_chunks = GTF.iterator_sorted_chunks( GTF.flat_gene_iterator( GTF.iterator( IOTools.openFile( all_gtf ) ) ),
                                              sort_by = "contig-strand-start-end" )

    # create a generator object that yields list of exons which are chunks, 
    # joined if the start of a chunk is less than merge_dist away from end
    # of the previous chunk.
    merge_dist = PARAMS[ "filter_merge_dist" ]
    merged_chunks = P10.collateChunks( gtf_chunks, merge_dist  )

    # iterate through chunks, if they contain more than one gene_id, then rename
    # write to outfile. 
    n_genes_unchanged = 0
    n_genes_merged = 0
    n_merged_genes = 0
    tmpf = P.getTempFilename( "./filter_lncrna" )
    outf = IOTools.openFile( tmpf, "w" )
    outf2 = IOTools.openFile( re.sub( "gtf", "tsv", outfile ), "w" ) 
    for chunk in merged_chunks:
        gene_ids = set( [exon.gene_id for exon in chunk ] )
        n_gene_ids = len( gene_ids )
        if n_gene_ids > 1:
            n_genes_merged += n_gene_ids
            n_merged_genes += 1
            exon_status_locus = None
            new_gene_id = None
            # check if any of the merged gene models are multi-exon
            if True in [ re.search( "me", gene_id ) != None for gene_id in gene_ids ]:
                for gene_id in gene_ids:
                    # set new gene_id based on the first multi-exon id encountered
                    if re.search( "me", gene_id ):
                        new_gene_id = re.sub( "me", "mm", gene_id )
                        exon_status_locus = "m"
                        break
            # check if any of the merged gene models are single-exon
            elif True in [ re.search( "se", gene_id ) != None for gene_id in gene_ids ]:
                for gene_id in gene_ids:
                    # set new gene_id based on the first single-exon id encountered
                    if re.search( "se", gene_id ):
                        new_gene_id = re.sub( "se", "sm", gene_id )
                        exon_status_locus = "s"
                        break
            elif True in [ re.search( "sm", gene_id ) != None for gene_id in gene_ids ]:
                for gene_id in gene_ids:
                    # set new gene_id based on the first merged single-exon id encountered
                    if re.search( "sm", gene_id ):
                        new_gene_id = gene_id
                        exon_status_locus = "s"
                        break
            else:
                raise ValueError( "Unrecognized gene_id list %s. Ids need to be"
                                   "either LNCGme* or LNCGse*" % ",".join( gene_ids ) )

            assert new_gene_id, "Error in code: no new gene_id assigned!"

            # write merged gene ids to tsv
            gene_ids = ",".join( gene_ids )
            outf2.write( gene_ids + "\t" + new_gene_id + "\n" )

            # write all exons to outfile
            for exon in chunk:
                exon.setAttribute( "gene_oId2", exon.gene_id )
                exon.setAttribute( "gene_id", new_gene_id )
                exon.setAttribute( "exon_status_locus", exon_status_locus )
                outf.write( str( exon ) + "\n" )
        else:
            n_genes_unchanged += 1
            for exon in chunk:
                outf.write( str( exon ) + "\n" )
    
    E.info( "%s genes are present prior to merging" % str(n_genes_in) )
    E.info( "%s genes are merged to form %s genes" % ( str(n_genes_merged), str(n_merged_genes) ) )
    E.info( "%s genes are left un-merged" % n_genes_unchanged )
    assert n_genes_in == n_genes_unchanged + n_genes_merged, "Missing genes after merging!"
    outf.close()
    outf2.close()

    to_cluster = False
    statement = ( "cat %(tmpf)s |"
                  " python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=sort --sort-order=gene+transcript"
                  "  --log=%(outfile)s.log"
                  " | gzip > %(outfile)s" )
    P.run()
    os.unlink( me_introns_gtf )
    os.unlink( tmpf )


#################################################################################
#################################################################################
#################################################################################
# Section: Filter lncRNAs on coding potential
#################################################################################
#################################################################################
## subsection: calculate lncRNA CPC score
#################################################################################
@follows( mkdir( "filter_lncrna_coding_potential") )
@transform( mergeMEandSELncRNA, 
            regex( "(?:.+)/(.+)_merged.gtf.gz" ), 
            r"filter_lncrna_coding_potential/\1.fasta" )
def buildLncRNAFasta( infile, outfile ):
    """
    Create fasta file prior to running CPC
    """
    genome = os.path.join( PARAMS["genome_dir"], PARAMS["genome"] + ".fasta" )
    statement = ( "zcat %(infile)s |"
                  " python %(scriptsdir)s/gff2fasta.py"
                  "  --genome-file=%(genome)s"
                  "  --log=%(outfile)s.log"
                  "  --is-gtf"
                  " > %(outfile)s" )
    P.run()


@transform( buildLncRNAFasta, suffix( ".fasta" ), "_cpc.result" )
def runCPC( infile, outfile ):
    # farm.py is called from within cpc.sh
    assert P.which("farm.py"), "farm.py needs to be in $PATH for cpc to run"

    # Default cpc parameters don't work with later versions of blast
    E.info("Running cpc with blast version:%s" % P.which("blastx"))

    result_evidence = P.snip( outfile, ".result" ) + ".evidence"
    working_dir = os.path.join( os.path.dirname( infile ), "lncRNA_cpc" )
    statement = ("%(scriptsdir)s/cpc.sh"
                 " %(infile)s"
                 " %(outfile)s"
                 " %(working_dir)s"
                 " %(result_evidence)s")
    P.run()


@transform( runCPC, regex( "(?:.+)/(.+).result" ), r"\1.load" )
def loadCPC( infile, outfile ):
    options = ( "--header-names=transcript_id,feature,Coding,CP_score"
                " --add-index=transcript_id" )
    P.load( infile, outfile, options = options )

#################################################################################
## subsection: calculate ENSEMBL lincrna CPC score
#################################################################################
@follows( mkdir( "filter_lncrna_coding_potential", loadCPC ) )
@transform( os.path.join( PARAMS["location_transcriptfiles" ], 
                          "reference.gtf.gz" ),
            regex( "(.+)/reference.gtf.gz" ),
            r"./filter_lncrna_coding_potential/ensembl_lincRNA.gtf.gz" )
def extractEnsemblLincRNAs_CPC( infile, outfile ):
    outf = IOTools.openFile( outfile, "w" )

    for gtf in GTF.iterator( IOTools.openFile( infile ) ):
        if gtf.source == "lincRNA":
            outf.write( str( gtf ) + "\n" )
        else:
            continue
    outf.close()


@transform( extractEnsemblLincRNAs_CPC,
            suffix( ".gtf.gz" ), 
            ".fasta" )
def buildEnsemblLincRNAFasta( infile, outfile ):
    """
    Create fasta file prior to running CPC
    """
    genome = os.path.join( PARAMS["genome_dir"], PARAMS["genome"] + ".fasta" )
    statement = ( "zcat %(infile)s |"
                  " python %(scriptsdir)s/gff2fasta.py"
                  "  --genome-file=%(genome)s"
                  "  --log=%(outfile)s.log"
                  "  --is-gtf"
                  " > %(outfile)s" )
    P.run()


@transform( buildEnsemblLincRNAFasta, suffix( ".fasta" ), "_cpc.result" )
def runEnsemblLincRNACPC( infile, outfile ):
    """
    Is run with ncbiblast-2.2.25+/bin/blastx as the blast parameters have changed
    """

    # farm.py is called from within cpc.sh
    assert P.which("farm.py"), "farm.py needs to be in $PATH for cpc to run"

    # Default cpc parameters don't work with later versions of blast
    E.info("Running cpc with blast version:%s" % P.which("blastx"))

    result_evidence = P.snip( outfile, ".result" ) + ".evidence"
    working_dir = os.path.join( os.path.dirname( infile ), "ensembl_cpc" )
    statement = ("%(scriptsdir)s/cpc.sh"
                 " %(infile)s"
                 " %(outfile)s"
                 " %(working_dir)s"
                 " %(result_evidence)s")
    P.run()


@transform( runEnsemblLincRNACPC, regex( "(?:.+)/(.+).result" ), r"\1.load" )
def loadEnsemblLincRNACPC( infile, outfile ):
    options = ( "--header-names=transcript_id,feature,Coding,CP_score"
                " --add-index=transcript_id" )
    P.load( infile, outfile, options = options )

    
#################################################################################
## subsection: calculate lncRNA phyloCSF score
#################################################################################
@follows( mkdir( "filter_lncrna_phyloCSF" ) )
@merge( os.path.join( PARAMS[ "phyloCSF_location_axt" ], "*.axt.gz" ), 
        "./filter_lncrna_phyloCSF/filtered_alignment.maf.gz" )
def createMAFAlignment( infiles, outfile ):
    """
    Takes all .axt files in the input directory, filters them to remove
    files based on supplied regular expressions, converts to a single maf file
    using axtToMaf, filters maf alignments under a specified length.
    """
    outfile = P.snip( outfile, ".gz" )
    axt_dir = PARAMS[ "phyloCSF_location_axt" ] 
    to_ignore = re.compile( PARAMS[ "phyloCSF_ignore" ] )

    axt_files = []
    for axt_file in os.listdir( axt_dir ):
        if axt_file.endswith( "net.axt.gz" ) and not to_ignore.search( axt_file ):
            axt_files.append( os.path.join( axt_dir,  axt_file ) )
    axt_files = (" ").join( sorted( axt_files ) )
    
    E.info( "axt files from which MAF alignment will be created: %s" % axt_files )

    target_genome = PARAMS[ "phyloCSF_target_genome" ]
    target_contigs = os.path.join( PARAMS[ "annotations_dir" ], 
                                   PARAMS_ANNOTATIONS[ "interface_contigs" ] )
    query_genome = PARAMS[ "phyloCSF_query_genome" ] 
    query_contigs = os.path.join( PARAMS[ "phyloCSF_query_assembly" ],
                                  PARAMS_ANNOTATIONS[ "interface_contigs" ] )

    tmpf1 = P.getTempFilename( "./filter_lncrna_phyloCSF" )
    tmpf2 = P.getTempFilename( "./filter_lncrna_phyloCSF" )
    to_cluster = False
    # concatenate axt files, then remove headers
    statement = ( "zcat %(axt_files)s"
                  " > %(tmpf1)s;"
                  " checkpoint;"
                  " axtToMaf "
                  "  -tPrefix=%(target_genome)s."
                  "  -qPrefix=%(query_genome)s."
                  "  %(tmpf1)s"
                  "  %(target_contigs)s"
                  "  %(query_contigs)s"
                  "  %(tmpf2)s" )
    P.run()

    E.info( "Temporary axt file created %s" % os.path.abspath( tmpf1 ) )
    E.info( "Temporary maf file created %s" % os.path.abspath( tmpf2 ) )    

    removed = P.snip( outfile, ".maf" ) + "_removed.maf"
    to_cluster = False
    filtered = PipelineLncRNA.filterMAF( tmpf2, 
                                         outfile, 
                                         removed,
                                         PARAMS[ "phyloCSF_filter_alignments" ] )
    E.info( "%s blocks were ignored in MAF alignment"
            " because length of target alignment was too short" % filtered[0] )
    E.info( "%s blocks were output to filtered MAF alignment" % filtered[1] )

    os.unlink( tmpf1 )
    os.unlink( tmpf2 )
    to_cluster = False
    statement = ( "gzip %(outfile)s;"
                  " gzip %(removed)s" )
    P.run()


@follows( mkdir( "filter_lncrna_phyloCSF" ) )
@transform( mergeMEandSELncRNA, 
            regex( "(.+)/(.+).gtf.gz" ), 
            r"./filter_lncrna_phyloCSF/\2.bed.gz" )
def convertGTFToBed12( infile, outfile ):
    """
    Transform the lncrna_final.gtf.gz into lncrna_final.bed
    """
    PipelineLncRNA.gtfToBed12( infile, outfile, "transcript" )


@merge( [ convertGTFToBed12, createMAFAlignment ], 
        "./filter_lncrna_phyloCSF/lncRNA_aligned.fasta.gz" )
def extractFastaAlignments( infiles, outfile ):
    """
    Recieves a MAF file containing pairwise alignments and a gtf12 file
    containing intervals. Outputs a single fasta file containing aligned
    sequence for each interval.
    """
    bed_file, maf_file = infiles
    maf_tmp = P.getTempFilename( "./filter_lncrna_phyloCSF" )
    to_cluster = False
    statement = ( "gunzip -c %(maf_file)s > %(maf_tmp)s" )
    P.run()

    target_genome = PARAMS[ "genome" ]
    query_genome = PARAMS[ "phyloCSF_query_genome" ] 

    genome_file = os.path.join( PARAMS[ "genome_dir" ], PARAMS[ "genome" ] )

    gene_models = PipelineLncRNA.extractMAFGeneBlocks( bed_file, 
                                                       maf_tmp, 
                                                       genome_file,
                                                       outfile, 
                                                       target_genome, 
                                                       query_genome, 
                                                       keep_gaps = False )
    E.info( "%i gene_models extracted" % gene_models )
    os.unlink( maf_tmp )


@follows( mkdir( "./filter_lncrna_phyloCSF/fasta_files" ) )
@split( extractFastaAlignments, 
        "./filter_lncrna_phyloCSF/fasta_files/*.fasta" )
def splitFasta( infile, outfiles ):
    out_dir = "./filter_lncrna_phyloCSF/fasta_files"

    name_dict = {}
    for mapping in PARAMS[ "phyloCSF_map_species_names" ].split(","):
        pair = mapping.split( ":" )
        key = ">" + pair[0]
        value = ">" + pair[1]
        name_dict[key] = value
    E.info( "Name mapping: %s" % name_dict )

    PipelineLncRNA.splitAlignedFasta( infile, out_dir, name_dict )


@transform( splitFasta, suffix( ".fasta" ), ".phyloCSF" )
def runPhyloCSF( infile, outfile ):
    phylogeny = PARAMS[ "phyloCSF_phylogeny" ]
    n_frames = int( PARAMS[ "phyloCSF_n_frames" ] )
    if PARAMS[ "phyloCSF_options" ]:
        options = PARAMS[ "phyloCSF_options" ]
    else: 
        options = ""

    species = []
    for mapping in PARAMS[ "phyloCSF_map_species_names" ].split(","):
        species.append( mapping.split( ":" )[1] )
    species = ",".join( species )

    to_cluster=True
    statement = ( "PhyloCSF %(phylogeny)s"
                  "  %(infile)s"
                  "  --frames=%(n_frames)s"
                  "  --species=%(species)s"
                  " %(options)s"
                  " > %(outfile)s" )
    P.run()


@merge( runPhyloCSF, "./filter_lncrna_phyloCSF/lncRNA_phyloCSF.tsv.gz" )
def mergePhyloCSF( infiles, outfile ):
    file_names = " ".join( [ x for x in infiles ] ) 
    statement = '''
                python %(scriptsdir)s/combine_tables.py
                 --no-titles
                 --cat=CAT
                 --missing-value=0
                 --log=%(outfile)s.log
                 %(file_names)s
                | gzip > %(outfile)s
                 '''
    P.run()


@transform( mergePhyloCSF, suffix( ".tsv.gz" ), "_scores.tsv.gz" )
def postProcessPhyloCSF( infile, outfile ):
    """
    Process the collated phyloCSF output to return the gene_id, transcript_id
    phyloCSF_score (maxScore(decibans)), and start & stop of the scored 
    interval.
    """
    outf = IOTools.openFile( outfile, "w" )
    outf.write( "\t".join( [ "gene_id", 
                             "transcript_id", 
                             "phyloCSF_score", 
                             "start", 
                             "stop" ] ) + "\n" )

    for line in IOTools.openFile( infile ).readlines():
        line = line.split()
        ids = P.snip(os.path.basename( line[0] ), ".phyloCSF" )
        gene_id, transcript_id = ids.split( "__" )
        outf.write( "\t".join( [ gene_id, 
                                 transcript_id, 
                                 line[3], 
                                 line[4], 
                                 line[5] ] ) + "\n" )
    outf.close()


@transform( postProcessPhyloCSF,
            regex( "(.+)/(.+).tsv.gz" ),
            r"\2.load" )
def loadPhyloCSF( infile, outfile ):
    P.load( infile, outfile, options="--add-index=gene_id --add-index=transcript_id" )


#################################################################################
## subsection: calculate ENSEMBL lincrna phyloCSF score
#################################################################################
@follows( mkdir( "./filter_lncrna_phyloCSF" ), mergePhyloCSF )
@transform( os.path.join( PARAMS["location_transcriptfiles" ], 
                          "reference.gtf.gz" ),
            regex( "(.+)/reference.gtf.gz" ),
            r"./filter_lncrna_phyloCSF/ensembl_lincRNA.gtf.gz" )
def extractEnsemblLincRNAs_phyloCSF( infile, outfile ):
    outf = IOTools.openFile( outfile, "w" )

    for gtf in GTF.iterator( IOTools.openFile( infile ) ):
        if gtf.source == "lincRNA":
            outf.write( str( gtf ) + "\n" )
        else:
            continue
    outf.close()


@transform( extractEnsemblLincRNAs_phyloCSF,
            suffix( ".gtf.gz" ), 
            r".bed.gz" )
def convertEnsemblGTFToBed12( infile, outfile ):
    """
    Transform the lncrna_final.gtf.gz into lncrna_final.bed
    """
    PipelineLncRNA.gtfToBed12( infile, outfile, "transcript" )


@merge( [ convertEnsemblGTFToBed12, createMAFAlignment ], 
        "./filter_lncrna_phyloCSF/ensembl_lincRNA_aligned.fasta.gz" )
def extractEnsemblFastaAlignments( infiles, outfile ):
    """
    Recieves a MAF file containing pairwise alignments and a gtf12 file
    containing intervals. Outputs a single fasta file containing aligned
    sequence for each interval.
    """
    bed_file, maf_file = infiles
    maf_tmp = P.getTempFilename( "./filter_lncrna_phyloCSF" )
    to_cluster = False
    statement = ( "gunzip -c %(maf_file)s > %(maf_tmp)s" )
    P.run()

    target_genome = PARAMS[ "genome" ]
    query_genome = PARAMS[ "phyloCSF_query_genome" ] 

    genome_file = os.path.join( PARAMS[ "genome_dir" ], PARAMS[ "genome" ] )

    gene_models = PipelineLncRNA.extractMAFGeneBlocks( bed_file, 
                                                       maf_tmp, 
                                                       genome_file,
                                                       outfile, 
                                                       target_genome, 
                                                       query_genome, 
                                                       keep_gaps = False )
    E.info( "%i gene_models extracted" % gene_models )
    os.unlink( maf_tmp )


@follows( mkdir( "./filter_lncrna_phyloCSF/ensembl_fasta_files" ) )
@split( extractEnsemblFastaAlignments, 
        "./filter_lncrna_phyloCSF/ensembl_fasta_files/*.fasta" )
def splitEnsemblFasta( infile, outfiles ):
    out_dir = "./filter_lncrna_phyloCSF/ensembl_fasta_files"

    name_dict = {}
    for mapping in PARAMS[ "phyloCSF_map_species_names" ].split(","):
        pair = mapping.split( ":" )
        key = ">" + pair[0]
        value = ">" + pair[1]
        name_dict[key] = value
    E.info( "Name mapping: %s" % name_dict )

    PipelineLncRNA.splitAlignedFasta( infile, out_dir, name_dict )


@transform( splitEnsemblFasta, suffix( ".fasta" ), ".phyloCSF" )
def runEnsemblPhyloCSF( infile, outfile ):
    phylogeny = PARAMS[ "phyloCSF_phylogeny" ]
    n_frames = int( PARAMS[ "phyloCSF_n_frames" ] )
    if PARAMS[ "phyloCSF_options" ]:
        options = PARAMS[ "phyloCSF_options" ]
    else: 
        options = ""

    species = []
    for mapping in PARAMS[ "phyloCSF_map_species_names" ].split(","):
        species.append( mapping.split( ":" )[1] )
    species = ",".join( species )

    to_cluster=True
    statement = ( "PhyloCSF %(phylogeny)s"
                  "  %(infile)s"
                  "  --frames=%(n_frames)s"
                  "  --species=%(species)s"
                  " %(options)s"
                  " > %(outfile)s" )
    P.run()


@merge( runEnsemblPhyloCSF, "./filter_lncrna_phyloCSF/ensembl_lincRNA_phyloCSF_scores.tsv.gz" )
def mergeEnsemblPhyloCSF( infiles, outfile ):
    """
    Iterate over the ensembl linc transcript phyloCSF scores and output 
    gene_id transcript_id score, start, stop
    Doing this with combine_tables.py causes problems:
      OSError: [Errno 7] Argument list too long
    """
    outf = IOTools.openFile(outfile, "w")
    outf.write("gene_id\ttranscript_id\tphyloCSF_score\tstart\tstop\n")

    for inf in infiles:
        line = IOTools.openFile(inf).readline()
        line = line.split()
        ids = P.snip(os.path.basename(line[0]), ".fasta")
        gene_id, transcript_id = ids.split("__")
        score, start, stop = line[-3:]
        outf.write("\t".join([gene_id, transcript_id, score, start, stop]) + "\n")

    outf.close()


@transform( mergeEnsemblPhyloCSF,
            regex( "(.+)/(.+).tsv.gz" ),
            r"\2.load" )
def loadEnsemblPhyloCSF( infile, outfile ):
    P.load( infile, outfile, options="--add-index=gene_id --add-index=transcript_id" )


#################################################################################
## subsection: filter lncRNAs on CPC and phyloCSF score
#################################################################################
@transform( mergeMEandSELncRNA, 
            suffix( "_merged.gtf.gz" ),
            add_inputs( loadCPC, loadPhyloCSF ),
            "_cp_filtered.gtf.gz" )
def filterlncRNAsOnCodingPotential( infiles, outfile ):
    """
    """
    lnc_gtf, CPC_table, pCSF_table = infiles
    CPC_table = P.snip( os.path.basename(CPC_table), ".load" )
    pCSF_table = P.snip( os.path.basename(pCSF_table), ".load" )

    # fetch a dataframe containing the CP and PhyloCSF scores for each transcript
    statement = ( "SELECT"
                  "  a.gene_id,"
                  "  a.transcript_id,"
                  "  a.phyloCSF_score,"
                  "  b.CP_score"
                  " FROM %(pCSF_table)s AS a"
                  " INNER JOIN %(CPC_table)s AS b"
                  " ON a.transcript_id = b.transcript_id" % locals() )
    df = PU.fetch_DataFrame( statement )

    # retrieve the filtering thresholds for each transcript. 
    CP_threshold = int(PARAMS["filter_coding_potential_threshold"])

    pCSF_threshold = int(PARAMS["filter_phylocsf_threshold"])
    
    if PARAMS["filter_cp_level"] == "transcript":
        # if filtering transcripts that fail both CP thresholds
        df_pass = df.query( "CP_score < %i | phyloCSF_score < %i"
                            % ( CP_threshold, pCSF_threshold ) )
        # retrieve a list of gene_ids that pass this filter threshold. 
        genes_pass = [ x for x in df_pass["gene_id"] ]
        genes_pass = set( genes_pass )
    elif PARAMS["filter_cp_level"] == "locus": 
        # if filtering loci that fail both CP thresholds
        df_fail_cp = df.query( "CP_score > %i" % CP_threshold )
        df_fail_pCSF = df.query( "phyloCSF_score > %i" % pCSF_threshold )

        genes_fail_cp = set( [ x for x in df_fail_cp["gene_id"] ] )
        genes_fail_pCSF = set( [ x for x in df_fail_pCSF["gene_id"] ] )
    
        genes_fail = genes_fail_cp & genes_fail_pCSF
        # genes_pass = [ x for x in df["gene_id"] if not in genes_fail ]
        genes_pass = []
        # iterate through series
        for gene in df["gene_id"]:
            if gene in genes_fail:
                continue
            else:
                genes_pass.append( gene )
        genes_pass = set( genes_pass )
    else: 
        raise ValueError( "Unrecognised filter option: %s" % PARAMS["filter_cp_level"] )

    # write tables 
    out_filtered_tsv = P.snip(outfile, ".gtf.gz") + ".tsv"
    out_rejected_tsv = P.snip(outfile, "_filtered.gtf.gz") + "_rejected.tsv"
    df_out_pass = df[df["gene_id"].isin(genes_pass)]
    df_out_fail = df[~df["gene_id"].isin(genes_pass)]
    df_out_pass.to_csv(out_filtered_tsv, sep="\t")
    df_out_fail.to_csv(out_rejected_tsv, sep="\t")
    
    # iterate through lncRNA gtf and discard those that fail this filtering
    outf_filtered = IOTools.openFile( outfile, "w" )
    outf_rejected = P.snip( outfile, "_filtered.gtf.gz" ) + "_rejected.gtf.gz"
    outf_rejected = IOTools.openFile( outf_rejected, "w" )

    for gtf in GTF.iterator( IOTools.openFile( lnc_gtf ) ):
        if gtf.gene_id in genes_pass:
            outf_filtered.write( str(gtf) + "\n" )
        else:
            outf_rejected.write( str(gtf) + "\n" )

    outf_filtered.close()
    outf_rejected.close()

#################################################################################
## subsection: finish merging and filtering lncRNAs
#################################################################################
# This is set up so that 'Filter lncRNAs on coding potential' doesn't
# necessarily have to be run.

if PARAMS["filter_filter_on_coding_potential"]:
    in_func = filterlncRNAsOnCodingPotential
else:
    in_func = mergeMEandSELncRNA

@transform( in_func, 
            regex( "(.+)/(.+).gtf.gz" ), 
            r"\1/lncRNA_filtered.gtf.gz" )
def mergeAndFilterLncRNAs( infile, outfile ):
    shutil.copyfile( infile, outfile )


#################################################################################
#################################################################################
#################################################################################
# Section: Re-classify merged and filtered lncRNAs
#################################################################################
#################################################################################
## subsection: classify merged & filtered lncRNAs based on genomic location
#################################################################################
## prior to this some merged loci were a combination of different classifications
@transform( mergeAndFilterLncRNAs, 
            suffix( "_filtered.gtf.gz" ), 
            add_inputs( os.path.join( PARAMS["location_transcriptfiles"], 
                                      "reference.gtf.gz" ) ),
            "_final.gtf.gz" )
def classifyMergedLncRNAs( infiles, outfile ):
    lncRNA_gtf, reference_gtf = infiles

    dist_upstr = PARAMS[ "filter_distance_upstream" ]
    dist_dstr = PARAMS[ "filter_distance_downstream" ]
    reference_biotypes = PARAMS[ "filter_classify_vs" ].split( "," )

    # preprocess reference file
    tmpf_ref = P.getTempFile( "/ifs/scratch" )
    tmpf_ref_name = tmpf_ref.name
    for gtf in GTF.iterator( IOTools.openFile( reference_gtf ) ):
        if gtf.source in reference_biotypes:
            if gtf.feature == "exon": 
                tmpf_ref.write( str( gtf ) + "\n" )
            else: continue
        else: continue
    tmpf_ref.close()

    tmp_gtf = P.getTempFilename( "/ifs/scratch" )
    tempdir = PipelineLncRNA.reClassifyLncRNAGenes( lncRNA_gtf,
                                                    tmpf_ref_name, 
                                                    tmp_gtf, 
                                                    dist_upstr,
                                                    dist_dstr, 
                                                    wdir = "")

    shutil.rmtree( tempdir )
    os.unlink( tmpf_ref_name )

    statement = ( "cat %(tmp_gtf)s |"
                  " python %(scriptsdir)s/gtf2gtf.py"
                  " --method=sort --sort-order=gene+transcript"
                  " --log=%(outfile)s.log |"
                  " gzip > %(outfile)s" )
    P.run()

#################################################################################
## subsection: merge the all of the lncRNAs with protein coding geneset
#################################################################################
@merge( [ classifyMergedLncRNAs, 
          os.path.join( PARAMS[ "location_transcriptfiles" ], 
                        "refcoding.gtf.gz" ) ],
        "./lncRNA_refcoding.gtf.gz" )
def combineRefcodingAndLncRNA( infiles, outfile ):
    gtf = " ".join( infiles )
    tmpf = P.getTempFilename( "." )

    statement = ( "zcat %(gtf)s"
                  "  > %(tmpf)s;"
                  " cat %(tmpf)s |"
                  " python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=sort --sort-order=gene"
                  "  --log=%(outfile)s.log |"
                  " gzip >"
                  " %(outfile)s" )
    P.run()
    os.unlink( tmpf )

          
@follows( combineRefcodingAndLncRNA )
def filterLncRNA():
    pass


#################################################################################
#################################################################################
#################################################################################
# Section: Load various useful statistics
#################################################################################
#################################################################################
## subsection: load summary of relationship between transcript_id & gene_id
#################################################################################
@transform( classifyMergedLncRNAs, 
            suffix( "_final.gtf.gz" ), 
            "_gene_vs_transcript_id.tsv.gz" )
def mapTranscriptToGeneID( infile, outfile ):
    interval_ids = []
    for gtf in GTF.iterator( IOTools.openFile( infile ) ):
        interval_ids.append( ( gtf.gene_id, gtf.transcript_id, gtf.source ) )
    outf = IOTools.openFile( outfile, "w" )
    outf.write( "gene_id\ttranscript_id\tbiotype\n" )
    for pair in set( interval_ids ):
        outf.write( pair[0] + "\t" + pair[1] + "\t" + pair[2] + "\n" )
    outf.close()    


@transform( mapTranscriptToGeneID, 
            regex( "(.+)/(.+).tsv.gz" ),
            r"./\2.load" )
def loadMappedIDs( infile, outfile ):
    P.load( infile, outfile, options="--add-index=gene_id --add-index=transcript_id" )            


#################################################################################
## subsection: extract lncRNA locations
#################################################################################
@follows(mkdir("characterize_lncrna_location"))
@transform(combineRefcodingAndLncRNA,
           regex(".*/lncRNA_refcoding.gtf.gz"), 
           r"characterize_lncrna_location/lncRNA_refcoding_coordinates.tsv.gz")
def extractLocusCoordinates(infile, outfile):
    """
    For the merged lncRNA file, extract co-ordinates as contig:start-end
    """
    tmpf = P.getTempFilename("/ifs/scratch")
    to_cluster = False
    statement = ("zcat %(infile)s |"
                 " python %(scriptsdir)s/gtf2gtf.py"
                 "  --method=sort"
                 "  --sort-order=gene"
                 "  --log=%(outfile)s.log |"
                 " python %(scriptsdir)s/gtf2gtf.py"
                 "  --method=merge-transcripts"
                 "  --log=%(outfile)s.log"
                 " > %(tmpf)s")
    P.run()

    outf = IOTools.openFile(outfile, "w")
    outf.write("gene_id\tLocation\tStrand\n")
    for gtf in GTF.iterator(IOTools.openFile(tmpf)):
        gene_id = gtf.gene_id
        contig = gtf.contig
        start = gtf.start
        end = gtf.end
        strand = gtf.strand
        interval = contig + ":" + str(start) + "-" + str(end)
        outf.write( "\t".join([gene_id, interval, strand]) + "\n" )

    outf.close()


@transform(extractLocusCoordinates, regex(".+/(.+).tsv.gz"), r"\1.load")
def loadLocusCoordinates(infile, outfile):
    P.load(infile, outfile)


#################################################################################
## subsection: load a table of lncRNA and refcoding biotypes
#################################################################################
@transform( combineRefcodingAndLncRNA, 
            suffix( ".gtf.gz" ),
            "_biotypes.load" )
def summarizelncRNARefcodingBiotypes( infile, outfile ):
    """
    Load a csvdb table that contains the source field for each lncRNA and 
    refcoding gene model.
    """
    tmpfile = P.getTempFilename( "/ifs/scratch" )

    tmpf = IOTools.openFile( tmpfile, "w" )
    tmpf.write( "gene_id\tbiotype\n" )
    for gene in GTF.flat_gene_iterator( GTF.iterator( IOTools.openFile( infile ) ) ):
        assert len(set([gtf.source for gtf in gene])) == 1, \
            "Multiple biotypes for gene %s" % gene[0].gene_id
        line_out = [ gene[0].gene_id, gene[0].source ]
        tmpf.write( "\t".join( line_out ) + "\n" )
    tmpf.close()

    P.load( tmpfile, 
            outfile, 
            options="--add-index=gene_id --table=lncRNA_refcoding_biotypes" )
    
    os.unlink( tmpfile )


#################################################################################
# Section: load a table of refcoding common names
#################################################################################
@transform( combineRefcodingAndLncRNA, 
            suffix( ".gtf.gz" ), 
            "_gene_names.tsv.gz" )
def getRefcodingCommonNames( infile, outfile ):
    """
    Parse gtf file and retrieve common names. 
    Note there are 74 duplicate common names... these have been given a suffix.
    """
    P10.getCommonNames( infile, outfile )


@transform( getRefcodingCommonNames, 
            suffix( ".tsv.gz" ), 
            ".load" )
def loadRefcodingCommonNames( infile, outfile ):
    P.load( infile, outfile, options="--add-index=gene_id" )


#################################################################################
## subsection: calculate distance to nearest lncRNA on the antisense strand
#################################################################################
@follows(mkdir("characterize_antisense_distance"))
@transform(classifyMergedLncRNAs,
           regex(".+/(.+)_final.gtf.gz"),
          r"characterize_antisense_distance/\1_plus.gtf.gz" )
def splitLncRNAFile(infile, outfile):
    """
    Split the lncRNAs into those on the plus strand and those on the minus. 
    """
    lnc_plus = outfile
    lnc_minus = P.snip(lnc_plus, "_plus.gtf.gz") + "_minus.gtf.gz"

    lncRNA_plus = IOTools.openFile(lnc_plus, "w")
    lncRNA_minus = IOTools.openFile(lnc_minus, "w")
    for gtf in GTF.flat_gene_iterator(GTF.iterator(IOTools.openFile(infile))):
        if gtf[0].strand == "+":
            lncRNA_plus.write("\n".join([str(x) for x in gtf]) + "\n")
        elif gtf[0].strand == "-":
            lncRNA_minus.write("\n".join([str(x) for x in gtf]) + "\n")
        else:
            raise ValueError("Unrecognised strand: %s" % gtf[0].strand )
            
    lncRNA_plus.close()
    lncRNA_minus.close()


@transform(splitLncRNAFile, suffix(".gtf.gz"), "_nearest.tsv")
def findNearestPlus(infile, outfile):
    """
    Find the nearest antisense lncRNA to lncRNAs on the plus strand 
    """
    lnc_plus = infile
    lnc_minus = P.snip(lnc_plus, "_plus.gtf.gz") + "_minus.gtf.gz"

    statement = ("zcat %(lnc_plus)s |"
                 " python %(scriptsdir)s/gtf2table.py"
                 "  --counter=distance-genes"
                 "  --log=%(outfile)s.log"
                 "  --gff-file=%(lnc_minus)s |"
                 " cut -f1,2,3 | sort -k3n,3n > %(outfile)s")
    to_cluster = False
    P.run()


@transform(splitLncRNAFile, suffix("_plus.gtf.gz"), "_minus_nearest.tsv")
def findNearestMinus(infile, outfile):
    """
    Find the nearest antisense lncRNA to lncRNAs on the minus strand 
    """
    lnc_plus = infile
    lnc_minus = P.snip(lnc_plus, "_plus.gtf.gz") + "_minus.gtf.gz"

    statement = ("zcat %(lnc_minus)s |"
                 " python %(scriptsdir)s/gtf2table.py"
                 "  --counter=distance-genes"
                 "  --log=%(outfile)s.log"
                 "  --gff-file=%(lnc_plus)s |"
                 " cut -f1,2,3 | sort -k3n,3n > %(outfile)s")
    to_cluster = False
    P.run()


@merge([findNearestPlus, findNearestMinus], 
       "characterize_antisense_distance/lncRNA_nearest.tsv")
def collateNearestAntisenseLncRNA(infiles, outfile):
    """
    Combine tables
    """

    infiles = " ".join(infiles)
    statement = ("python %(scriptsdir)s/combine_tables.py"
                 " --cat=strand"
                 " --log=%(outfile)s.log"
                 " %(infiles)s > %(outfile)s")
    to_cluster=False
    P.run()




#################################################################################
#################################################################################
#### METASECTION #### LncRNA Descriptive Statistics ####
#################################################################################
#################################################################################
# Section:  lncRNA location summary
#################################################################################
#################################################################################
#################################################################################
## subsection: find nearest gene features for lncRNA dataset
#################################################################################
@follows( mkdir( "./location_lncrna_neighbours" ) )
@transform( classifyMergedLncRNAs, 
            regex( "(.+)/(.+)_final.gtf.gz" ),
            r"./location_lncrna_neighbours/\2_nearest_tss.tsv.gz" )
def findNearestTSS( infile, outfile ):

    antn_file = os.path.join( PARAMS["annotations_dir"], 
                              PARAMS_ANNOTATIONS["interface_tss_gene_bed"] )

    statement = ( "zcat %(infile)s |"
                  " python %(scriptsdir)s/gtf2table.py"
                  "  --counter=distance-tss"
                  "  --log=%(outfile)s.log"
                  "  --gff-file=%(antn_file)s"
                  "  --filename-format='bed'"
                  "  --genome-file=%(genome_dir)s/%(genome)s |"
                  " gzip > %(outfile)s" )
    P.run()


#################################################################################
## subsection: find nearest protein coding genes
#################################################################################
@follows( mkdir( "./location_lncrna_neighbours" ) )
@transform( classifyMergedLncRNAs, 
            regex( "(.+)/(.+)_final.gtf.gz" ), 
            add_inputs( os.path.join( PARAMS[ "location_transcriptfiles" ], 
                                      "refcoding.gtf.gz" ) ),
            r"./location_lncrna_neighbours/\2_nearest_refcoding_gene.tsv.gz" )
def findNearestProteinCodingGenes( infiles, outfile ):
    lncRNA_gtf, refcoding_gtf = infiles

    to_cluster = False
    statement = ( "zcat %(lncRNA_gtf)s |"
                  " python %(scriptsdir)s/gtf2table.py"
                  "  --counter=distance-genes"
                  "  --gff-file=%(refcoding_gtf)s"
                  "  --log=%(outfile)s.log |"
                  " gzip > %(outfile)s" )
    P.run()


@transform( findNearestProteinCodingGenes, 
            regex( "(.+)/(.+).tsv.gz" ), 
            r"./\2.load" )
def loadNearestProteinCodingGenes( infile, outfile ):
    P.load( infile, outfile, options="--add-index=gene_id" )


#################################################################################
## subsection: find nearest protein coding genes accounting for strand
#################################################################################
@follows( mkdir( "./location_lncrna_neighbours" ) )
@transform( classifyMergedLncRNAs, 
            regex( "(.+)/(.+)_final.gtf.gz" ), 
            add_inputs( os.path.join( PARAMS[ "location_transcriptfiles" ], 
                                      "refcoding.gtf.gz" ) ),
            r"location_lncrna_neighbours/\2_nearest_refcoding_gene_stranded.tsv.gz" )
def findNearestStrandedProteinCodingGenes( infiles, outfile ):
    lnc_gtf, refcoding_gtf = infiles

    # split refcoding and lncRNA gtfs by strand. 
    # run gtf2table on each of the four combinations. 
    # create nested dict of lnc_id : { id5_sense: dist5_sense, id3_sense: ... } 
    # output nested dict as table that also contains closest_id, closest_dist, closest_strand

    # get tempfiles
    tmpdir = P.getTempDir( "." )
    lnc_gtf_plus = P.getTempFilename( tmpdir )
    lnc_gtf_minus = P.getTempFilename( tmpdir )
    refcoding_gtf_plus = P.getTempFilename( tmpdir ) + ".gtf"
    refcoding_gtf_minus = P.getTempFilename( tmpdir ) + ".gtf"

    # split files
    for infiles in ( [ lnc_gtf, lnc_gtf_plus, lnc_gtf_minus ], 
                 [ refcoding_gtf, refcoding_gtf_plus, refcoding_gtf_minus ] ):
        inf, out_plus, out_minus = infiles
        out_plus = IOTools.openFile( out_plus, "w" )
        out_minus = IOTools.openFile( out_minus , "w" )

        for gtf in GTF.iterator( IOTools.openFile( inf ) ):
            if gtf.strand == "+":
                out_plus.write( str( gtf ) + "\n" )
            elif gtf.strand == "-":
                out_minus.write( str( gtf ) + "\n" )
            else:
                E.warn( "LncRNA %s has no strand assigned to it" % gtf.gene_id )
        out_plus.close()
        out_minus.close()

    # get tempfiles
    lnc_plus_sense = P.getTempFilename( tmpdir )
    lnc_plus_antisense = P.getTempFilename( tmpdir )
    lnc_minus_sense = P.getTempFilename( tmpdir )
    lnc_minus_antisense = P.getTempFilename( tmpdir )

    # calculate nearest sense and antisense genes for lncRNAs using gtf2table
    for infiles in ( [ lnc_gtf_plus, refcoding_gtf_plus, lnc_plus_sense ], 
                     [ lnc_gtf_plus, refcoding_gtf_minus, lnc_plus_antisense ], 
                     [ lnc_gtf_minus, refcoding_gtf_minus, lnc_minus_sense ], 
                     [ lnc_gtf_minus, refcoding_gtf_plus, lnc_minus_antisense] ):
        inf_lnc, inf_refcoding, outf = infiles
        to_cluster = True
        statement = ( "cat %(inf_lnc)s |"
                      " python %(scriptsdir)s/gtf2table.py"
                      "  --counter=distance-genes"
                      "  --gff-file=%(inf_refcoding)s"
                      "  --log=%(outfile)s.log "
                      " > %(outf)s" )
        P.run()

    # create nested dictionary
    closest_dict = collections.defaultdict( dict )

    # collect the nearest 5' and 3' genes from gtf2table output
    for inf in ( [ lnc_plus_sense, "_sense" ], 
                 [ lnc_minus_sense, "_sense" ], 
                 [ lnc_plus_antisense, "_antisense" ], 
                 [ lnc_minus_antisense, "_antisense" ] ):
        closest_id = "closest_id" + inf[1]
        closest_dist = "closest_dist" + inf[1]
        id5 = "id5" + inf[1]
        id3 = "id3" + inf[1]
        dist5 = "dist5" + inf[1]
        dist3 = "dist3" + inf[1] 

        header = True 
        for line in IOTools.openFile( inf[0] ).readlines():
            if header:
                header = False
                continue
            line = line.split()
            gene_id = line[0]
            # id5=line[4], dist5=line[5], id3=line[7], dist3=line[8]
            closest_dict[ gene_id ][ closest_id ] = line[1]
            closest_dict[ gene_id ][ closest_dist ] = line[2]
            closest_dict[ gene_id ][ id5 ] = line[4]
            closest_dict[ gene_id ][ dist5 ] = line[5]
            closest_dict[ gene_id ][ id3 ] = line[7] 
            closest_dict[ gene_id ][ dist3 ] = line[8]

    # write outfile
    outf = IOTools.openFile( outfile, "w" )
    outf.write( "gene_id\t"
                "closest_id\tclosest_dist\t"
                "id5_sense\tdist5_sense\t"
                "id5_antisense\tdist5_antisense\t"
                "id3_sense\tdist3_sense\t"
                "id3_antisense\tdist3_antisense\n" )

    # write 5' and 3' gene distances to outfile
    for key, value in closest_dict.iteritems():
        closest_sense =  int( closest_dict[ key ][ "closest_dist_sense" ] )
        closest_antisense = int( closest_dict[ key ][ "closest_dist_antisense" ] )
        if closest_sense < closest_antisense:
            closest_dist = closest_sense
            closest_id = closest_dict[ key ][ "closest_id_sense" ]
        elif closest_antisense < closest_sense:
            closest_dist = closest_antisense
            closest_id = closest_dict[ key ][ "closest_id_antisense" ]
        else: 
            closest_dist = closest_sense
            closest_id = closest_dict[ key ][ "closest_id_sense" ]           
            P.warn( "LncRNA %s is equidistant between sense and antisense genes."
                    " (Sense gene recorded as closest) "
                    " Sense gene_id: %s Sense dist: %s"
                    " Antisense gene_id: %s Antisense_dist: %s" % ( key,
                                                                    closest_dict[ key ][ "closest_id_sense" ], 
                                                                    closest_sense, 
                                                                    closest_dict[ key ][ "closest_id_antisense" ], 
                                                                    closest_antisense ) )
        out_list = [ key, 
                     closest_id, 
                     closest_dist, 
                     closest_dict[ key ][ "id5_sense" ], 
                     closest_dict[ key ][ "dist5_sense" ], 
                     closest_dict[ key ][ "id5_antisense" ], 
                     closest_dict[ key ][ "dist5_antisense" ], 
                     closest_dict[ key ][ "id3_sense" ], 
                     closest_dict[ key ][ "dist3_sense" ], 
                     closest_dict[ key ][ "id3_antisense" ], 
                     closest_dict[ key ][ "dist3_antisense" ] ]

        out_str = "\t".join( map( str, out_list ) ) + "\n" 
        outf.write( out_str )

    outf.close()
    shutil.rmtree( tmpdir )


@transform( findNearestStrandedProteinCodingGenes, 
            regex( "(.+)/(.+).tsv.gz" ), 
            r"./\2.load" )
def loadNearestStrandedProteinCodingGenes( infile, outfile ):
    P.load( infile, outfile, options="--add-index=gene_id" )


# #################################################################################
# ## subsection: find lncRNA overlap with repeat regions
# #################################################################################
# @follows( findNearestTSS )
# @transform( mergeMEandSELncRNA, 
#             regex( "(.+)/(.+)_merged.gtf.gz" ),
#             r"./characterize_lncrna_location/\2_repeat_overlap.tsv.gz" )
# def findOverlapWithRepeatRegions( infile, outfile ): 

#     antn_file = os.path.join( PARAMS["annotations_dir"], 
#                               PARAMS_ANNOTATIONS["interface_repeats_gff"] )

#     statement = ( "zcat %(infile)s |"
#                   " python %(scriptsdir)s/gtf2table.py"
#                   "  --counter=overlap"
#                   "  --log=%(outfile)s.log"
#                   "  --gff-file=%(antn_file)s"
#                   "  --genome-file=%(genome_dir)s/%(genome)s |"
#                   " gzip > %(outfile)s" )
#     P.run()
# #          findOverlapWithRepeatRegions, 



@follows( loadNearestProteinCodingGenes, 
          loadNearestStrandedProteinCodingGenes )
def runLncRNALocationStats():
    pass


#################################################################################
# Section: lncRNA expression summary
#################################################################################
#################################################################################
## subsection: generate raw and transformed count data using FeatureCounts/DESeq2
#################################################################################

@follows( mkdir( "./expression_featureCounts" ) )
@transform( os.path.join( PARAMS["location_bamfiles_filtered"], "*.bam"), 
            regex( "(.+)/(.+).bam" ), 
            add_inputs( combineRefcodingAndLncRNA ), 
            r"expression_featureCounts/\2_vs_lncRNA_refcoding.tsv.gz" )
def buildFeatureCounts( infiles, outfile ):
    """
    Runs PipelineRnaseq.runFeatureCounts() 
    Assuming library is fr-firststrand (i.e. featureCounts strand=2)
    Discard multimapping reads (as suggested by SA at CSAMA)
    Minimum mapping quality (-Q) set to 10
    """
    bamfile, annotations = infiles
    PipelineRnaseq.runFeatureCounts(
        annotations,
        bamfile,
        outfile,
        nthreads=4,
        strand=2,
        options='-Q10' )


@collate( buildFeatureCounts, 
          regex( "(.+)/(.+)_vs_(.+).tsv.gz" ),
          r"\1/\3_raw_counts.tsv.gz" )
def summarizeFeatureCounts( infiles, outfile ):
    """
    Collate count data into a single file to be read as matrix of counts
    (Lifted out of pipeline_rnaseqdiffexpression.py, with minor changes)
    """
    infiles = " ".join( infiles )
    statement = ( "python %(scriptsdir)s/combine_tables.py"
                  " --columns=1"
                  " --take=7"
                  " --use-file-prefix"
                  " --regex-filename='(.+)_vs.+.tsv.gz'"
                  " --log=%(outfile)s.log"
                  " %(infiles)s"
                  " | sed 's/Geneid/gene_id/'"
                  " | gzip > %(outfile)s" )
    P.run()


@transform( summarizeFeatureCounts, 
            suffix( "_raw_counts.tsv.gz" ),
            "_mean_raw_counts.tsv.gz" )
def extractMeanRawCounts( infile, outfile ):
    P10.summarizeCountTables( infile, outfile, rnd = True )


@transform( summarizeFeatureCounts, 
            suffix( "_raw_counts.tsv.gz" ),
            "_median_raw_counts.tsv.gz" )
def extractMedianRawCounts( infile, outfile ):
    P10.summarizeCountTables( infile, outfile, summary = "median", rnd = True )


@follows( mkdir("expression_DESeq_transformed_data") )
@merge( os.path.join( PARAMS["location_bamfiles_filtered"], "*.bam" ),
        "./expression_DESeq_transformed_data/colData.tsv" )
def generateDESeq2DesignFile( infiles, outfile ):
    """
    Generate a single file that can be read as the colData for a 
    DESeqDataSet object. Has sample_id as $1 and condition as $2
    """
  
    sample_ids = [ P.snip(os.path.basename( infile ), ".bam" ) for infile in infiles ]
    conditions = [ re.split( "-", sample_id )[1] for sample_id in sample_ids ]
    outf = IOTools.openFile( outfile, "w" )
    outf.write( "sample\tcondition\n" )
    for line in zip( sample_ids, conditions ):
        outf.write( "\t".join( line ) + "\n" )

    outf.close()


@jobs_limit( 1, "RGlobalEnv" )
@transform( summarizeFeatureCounts, 
            regex( "(?:.+)/(.+)_raw_counts.tsv.gz" ), 
            add_inputs( generateDESeq2DesignFile ), 
            r"expression_DESeq_transformed_data/\1_rlog_trans_counts.tsv.gz" )
def rlogTransformData( infiles, outfile ):
    """
    Create regularized log transformation of count data, together with 
    sundry diagnostic plots. 
    """
    count_file, design_file = infiles
    P10.transformCountData( count_file, design_file, outfile, transform="rlog" )


@jobs_limit( 1, "RGlobalEnv" )
@transform( rlogTransformData, 
            suffix("_counts.tsv.gz" ), 
            "_mean_counts.tsv.gz" )
def extractMeanRlogTransformedCountData( infile, outfile ):
    """
    Create regularized log transformation of count data, together with 
    sundry diagnostic plots. 
    """
    P10.summarizeCountTables( infile, outfile, summary = "mean" )


@jobs_limit( 1, "RGlobalEnv" )
@transform( rlogTransformData, 
            suffix("_counts.tsv.gz" ), 
            "_median_counts.tsv.gz" )
def extractMedianRlogTransformedCountData( infile, outfile ):
    """
    Create regularized log transformation of count data, together with 
    sundry diagnostic plots. 
    """
    P10.summarizeCountTables( infile, outfile, summary = "median" )


@transform( [ extractMeanRlogTransformedCountData, 
              extractMedianRlogTransformedCountData ],
            regex( "(.+)/(.+).tsv.gz" ),
            r"\2.load" )
def loadSummaryRlogTransformedCountData( infile, outfile ):
    """
    Is necessary for extracting module specific heatmaps for wgcna analysis
    """
    table_name = P.snip(os.path.basename(outfile), ".load")
    statement   = ( "zcat %(infile)s |"
                    " sed 's/Bcell-//g' |"
                    " python %(scriptsdir)s/csv2db.py"
                    "  --add-index=gene_id"
                    "  --table=%(table_name)s"
                    " > %(outfile)s" )
    P.run()


@follows( mkdir("expression_DESeq_transformed_data") )
@merge( os.path.join( PARAMS["location_bamfiles_filtered"], "*.bam" ),
        "./expression_DESeq_transformed_data/colData_summary.tsv" )
def generateDESeq2SummaryDesignFile( infiles, outfile ):
    """
    Generate a single file that can be read as the colData for a 
    DESeqDataSet object. Has sample_id as $1 and condition as $2
    """
    sample_ids = [ P.snip(os.path.basename( infile ), ".bam" ) for infile in infiles ]
    sample_ids  = list(set( [ "-".join( x.split("-")[:2] ) for x in sample_ids ] ) )
    print( sample_ids )
    conditions = [ re.split( "-", sample_id )[1] for sample_id in sample_ids ]
    print( conditions )
    outf = IOTools.openFile( outfile, "w" )
    outf.write( "sample\tcondition\n" )
    for line in zip( sample_ids, conditions ):
        outf.write( "\t".join( line ) + "\n" )

    outf.close()


@jobs_limit( 1, "RGlobalEnv" )
@transform( summarizeFeatureCounts, 
            regex( "(?:.+)/(.+)_raw_counts.tsv.gz" ), 
            add_inputs( generateDESeq2DesignFile ), 
            r"expression_DESeq_transformed_data/\1_vst_trans_counts.tsv.gz" )
def vstTransformData( infiles, outfile ):
    """
    Create variance stabilized transformation of count data, together
    with sundry diagnostic plots. 
    """
    count_file, design_file = infiles
    P10.transformCountData( count_file, design_file, outfile, transform="vst" )


@jobs_limit( 1, "RGlobalEnv" )
@transform( summarizeFeatureCounts, 
            regex( "(?:.+)/(.+)_raw_counts.tsv.gz" ), 
            add_inputs( generateDESeq2DesignFile ), 
            r"expression_DESeq_transformed_data/\1_fpm_trans_counts.tsv.gz" )
def fpmTransformData( infiles, outfile ):
    """
    Create fragments per million transformation of count data, 
    using DESeq's robust library normalization method.
    """
    count_file, design_file = infiles
    P10.transformCountData( count_file, design_file, outfile, transform="fpm" )

@jobs_limit( 1, "RGlobalEnv" )
@transform( fpmTransformData, 
            suffix("_counts.tsv.gz" ), 
            "_mean_counts.tsv.gz" )
def extractMeanfpmTransformedCountData( infile, outfile ):
    """
    Create regularized log transformation of count data, together with 
    sundry diagnostic plots. 
    """
    P10.summarizeCountTables( infile, outfile, summary = "mean" )


@jobs_limit( 1, "RGlobalEnv" )
@transform( fpmTransformData, 
            suffix("_counts.tsv.gz" ), 
            "_median_counts.tsv.gz" )
def extractMedianfpmTransformedCountData( infile, outfile ):
    """
    Create regularized log transformation of count data, together with 
    sundry diagnostic plots. 
    """
    P10.summarizeCountTables( infile, outfile, summary = "median" )



@jobs_limit( 1, "RGlobalEnv" )
@transform( summarizeFeatureCounts, 
            regex( "(?:.+)/(.+)_raw_counts.tsv.gz" ), 
            add_inputs( generateDESeq2DesignFile ), 
            r"expression_DESeq_transformed_data/\1_fpm_l_trans_counts.tsv.gz" )
def logfpmTransformData( infiles, outfile ):
    """
    Create fragments per million transformation of count data, 
    using DESeq's robust library normalization method, then log(x + 1) transform.
    """
    count_file, design_file = infiles
    P10.transformCountData( count_file, design_file, outfile, transform="logfpm" )


#################################################################################
## subsection: calculate upper-quartile normalised FPKM values using cuffdiff
#################################################################################
@follows( mkdir("expression_cuffdiff_fpkms") )
@merge( os.path.join( PARAMS["location_bamfiles_filtered"], "*.bam" ),
        "./expression_cuffdiff_fpkms/design.tsv" )
def generateCuffdiffDesignFile( infiles, outfile ):
    """
    Generate a simple design file for use with Expression.py functions
    """
    sample_ids = [ P.snip(os.path.basename( infile ), ".bam" ) for infile in infiles ]
    include = ["1"]*len( sample_ids )
    pair = ["0"]*len( sample_ids )
    conditions = [ re.split( "-", sample_id )[1] for sample_id in sample_ids ]
    outf = IOTools.openFile( outfile, "w" )

    outf.write( "track\tinclude\tcondition\tpair\n" )
    for line in zip( sample_ids, include, conditions, pair ):
        outf.write( "\t".join( line ) + "\n" )
    outf.close()


@merge( [ os.path.join( PARAMS["location_bamfiles_filtered"], "*.bam" ), 
          generateCuffdiffDesignFile,
          combineRefcodingAndLncRNA,
          os.path.join( PARAMS["annotations_dir"],
                        PARAMS_ANNOTATIONS["interface_geneset_all_gtf"] ) ],
        "expression_cuffdiff_fpkms/cuffdiff_results.tsv.gz" )
def runCuffdiff( infiles, outfile ):
    """
    Run Cuffdiff using Expression.py to obtain upper-quartile normalised fpkm 
    values
    """
    geneset_all = infiles.pop()
    geneset_file = infiles.pop()
    design_file = infiles.pop()
    bamfiles = infiles

    # build mask gtf to discount reads mapping to rRNA and chrM
    tmp_maskfile = P.getTempFilename("/ifs/scratch")
    tmpf = IOTools.openFile( tmp_maskfile, "w" )
    for gtf in GTF.iterator( IOTools.openFile( geneset_all ) ):
        if re.search( "rRNA", gtf.source ):
            tmpf.write( str(gtf) + "\n" )
        elif re.search( "chrM", gtf.source ):
            tmpf.write( str(gtf) + "\n" )
        else: 
            continue
    tmpf.close()

    # set cuffdiff options
    options = ( "--upper-quartile-norm"
                " --max-bundle-frags 2000000"
                " --library-type fr-firststrand" )

    # run cuffdiff from Expression.py
    Expression.runCuffdiff( bamfiles=bamfiles,
                            design_file=design_file, 
                            geneset_file=geneset_file,
                            outfile=outfile,
                            threads=4,
                            cuffdiff_options=options,
                            fdr=0.05,
                            mask_file=tmp_maskfile )
    os.unlink( tmp_maskfile )


@transform( runCuffdiff, 
            suffix( "_results.tsv.gz" ),
            "_id_mapping.p")
def mapCuffdiffIDs( infile, outfile ):
    """
    Create dictionary mapping cuffdiff IDs to sample IDs
    """
    mapping_file = os.path.join( infile + ".dir", 
                                 "read_groups.info.gz" )
    id_dict = P10.mapCuffdiffIDs( mapping_file )
    pickle.dump( id_dict, open( outfile, "wb" ) )


@transform( runCuffdiff,
            regex( "(.+)/cuffdiff_results.tsv.gz" ),
            add_inputs( mapCuffdiffIDs ),
            r"\1/lncRNA_refcoding_cuffdiff_fpkms.tsv.gz"  )
def extractPerSampleFPKMs( infiles, outfile ):
    infile, id_file = infiles
    infile = os.path.join( infile + ".dir", 
                           "genes.read_group_tracking.gz" )
    id_dict = pickle.load( open( id_file , "rb" ) )
    P10.extractPerSampleCuffdiffFPKM( infile, outfile, id_dict )


@transform( runCuffdiff,
            regex( "(.+)/cuffdiff_results.tsv.gz" ),
            add_inputs( mapCuffdiffIDs ),
            r"\1/lncRNA_refcoding_cuffdiff_fpkms.tsv.gz"  )
def extractPerSampleFPKMs_stacked( infiles, outfile ):
    infile, id_file = infiles
    infile = os.path.join( infile + ".dir", 
                           "genes.read_group_tracking.gz" )
    id_dict = pickle.load( open( id_file , "rb" ) )
    P10.extractPerSampleCuffdiffFPKM_stacked( infile, outfile, id_dict )


@jobs_limit( 1, "RGlobalEnv" )
@transform( extractPerSampleFPKMs,
            suffix( "_fpkms.tsv.gz" ), 
            "_l_fpkms.tsv.gz" )
def logTransformCuffdiffFPKMs( infile, outfile ):
    """
    log(x + 1) transformation of fpkm values
    """
    P10.transformCuffdiffFPKMs( infile, outfile )


@jobs_limit( 1, "RGlobalEnv" )
@transform( extractPerSampleFPKMs,
            suffix( "_fpkms.tsv.gz" ), 
            "_ls_fpkms.tsv.gz" )
def logTransformCuffdiffFPKMs_smallconstant( infile, outfile ):
    """
    log(x + 1e-4) transformation of fpkm values
    """
    P10.transformCuffdiffFPKMs( infile, outfile, constant = 0.00001 )


@transform( runCuffdiff,
            regex( "(.+)/cuffdiff_results.tsv.gz" ),
            add_inputs( mapCuffdiffIDs ),
            r"\1/lncRNA_refcoding_cuffdiff_mean_fpkms.tsv.gz"  )
def extractMeanFPKMs( infiles, outfile ):
    """
    There are 54 gene values that are either HIDATA or FAIL
    There is one gene ENSMUSG00000024610 that has no data for fol/mat (HIDATA)
    """
    infile, id_file = infiles
    infile = os.path.join( infile + ".dir", 
                           "genes.read_group_tracking.gz" )
    P10.calculateSummaryCuffdiffFPKM( infile, 
                                      outfile,
                                      id_file,
                                      stat="mean", 
                                      num_samples = 8, 
                                      to_exclude = False )


@transform( runCuffdiff,
            regex( "(.+)/cuffdiff_results.tsv.gz" ),
            add_inputs( mapCuffdiffIDs ),
            r"\1/lncRNA_refcoding_cuffdiff_median_fpkms.tsv.gz"  )
def extractMedianFPKMs( infiles, outfile ):
    """
    There are 54 gene values that are either HIDATA or FAIL
    There is one gene ENSMUSG00000024610 that has no data for fol/mat (HIDATA)
    """
    infile, id_file = infiles
    infile = os.path.join( infile + ".dir", 
                           "genes.read_group_tracking.gz" )
    P10.calculateSummaryCuffdiffFPKM( infile, 
                                      outfile,
                                      id_file,
                                      stat="median", 
                                      num_samples = 8,
                                      to_exclude = False )


@transform( [ logTransformCuffdiffFPKMs,
              logTransformCuffdiffFPKMs_smallconstant ],
            regex( "(.+)/lncRNA_refcoding_cuffdiff_(.+)_fpkms.tsv.gz" ),
            r"\1/lncRNA_refcoding_cuffdiff_mean_\2_fpkms.tsv.gz"  )
def extractMeanLogTransformedFPKMs( infile, outfile ):
    """
    Use summarizeCountTables to get mean of transformed datatable
    """
    P10.summarizeCountTables( infile, outfile, summary = "mean" )


@transform( [ logTransformCuffdiffFPKMs,
              logTransformCuffdiffFPKMs_smallconstant ],
            regex( "(.+)/lncRNA_refcoding_cuffdiff_(.+)_fpkms.tsv.gz" ),
            r"\1/lncRNA_refcoding_cuffdiff_median_\2_fpkms.tsv.gz"  )
def extractMedianLogTransformedFPKMs( infile, outfile ):
    """
    Use summarizeCountTables to get median of transformed datatable
    """
    P10.summarizeCountTables( infile, outfile, summary = "median" )


@jobs_limit( 1, "RGlobalEnv" )
@transform( [ extractPerSampleFPKMs, 
              extractMeanFPKMs,
              extractMedianFPKMs,
              logTransformCuffdiffFPKMs,
              extractMeanLogTransformedFPKMs,
              extractMedianLogTransformedFPKMs,
              logTransformCuffdiffFPKMs_smallconstant ],
            suffix( ".tsv.gz" ),
            "_meanSDPlot.png" )
def plotFPKMDispersions( infile, outfile ):
    """
    Use R library vsn to plot the relationship between standard deviation 
    and expression for transformed and untransformed data.
    """
    P10.plotCuffdiffFPKMDispersions( infile, outfile )


@transform( [ summarizeFeatureCounts, 
              rlogTransformData,
              fpmTransformData,
              extractPerSampleFPKMs,
              extractMedianFPKMs,
              extractMeanFPKMs,
              extractMedianLogTransformedFPKMs ],
            regex( "(?:.+)/(.+).tsv.gz" ), 
            r"\1.load" )
def loadReadCounts( infile, outfile ):
    P.load( infile, outfile, options="--add-index=gene_id" )


#################################################################################
## subsection: plot refcoding and lncRNA ordinations
#################################################################################
@jobs_limit( 1, "RGlobalEnv" )
@transform( [ summarizeFeatureCounts,
              rlogTransformData, 
              vstTransformData,
              fpmTransformData,
              logfpmTransformData,
              extractPerSampleFPKMs,
              logTransformCuffdiffFPKMs ],
            regex( "(.+)_(counts|fpkms).tsv.gz" ), 
            r"\1_pca.rds" )
def plotOrdinationsBasedOnExpression( infile, outfile ):
    """
    Run r script transcriptome_pca_plots.R, which will output PCA of supplied
    dataframe, along with sundry diagnostic plots. 
    Takes infile, dendrogram cut height, outfile_stub as arguments.
    """
    outf_stub = P.snip( outfile, "_pca.rds" )

    E.info( "Plotting PCA plots for %s" % infile )
    to_cluster = False
    statement = ( "Rscript /ifs/devel/projects/proj010/transcriptome_pca_plots.R"
                  " %(infile)s"
                  " 0.6"
                  " %(outf_stub)s" )
    P.run()
    E.info( "Completed PCA plots for %s" % infile )


@jobs_limit( 1, "RGlobalEnv" )
@follows( mkdir( "expression_ordination_plots" ) )
@split( [ summarizeFeatureCounts, 
          rlogTransformData, 
          extractPerSampleFPKMs, 
          logTransformCuffdiffFPKMs],
        regex( "(.+)/(.+)_(counts|fpkms).tsv.gz" ),
        add_inputs( combineRefcodingAndLncRNA ),
        r"expression_ordination_plots/\2_intergenic_lncRNA_pca.rds" )
def lncRNAOrdinationPlots( infiles, outfile ):
    """
    Plot pca plots of lncRNAs separate from protein coding genes.
    Also separates intergenic and non-intergenic lncRNAs into different plots.
    """

    count_file, gtf_file = infiles
    
    E.info( "Plotting ordinations for %s" % count_file )
    # create lists of gene_ids for all_lncRNA, intergenic_lncRNA, and refcoding
    lnc_all, lnc_intergenic, refcoding = set(), set(), set()
    for gtf in GTF.iterator( IOTools.openFile( gtf_file ) ):
        if re.match("ENSMUSG", gtf.gene_id):
            assert gtf.source == "protein_coding", "Non coding ensembl gene_id"
            refcoding.add( gtf.gene_id )
        elif re.match("LNCG", gtf.gene_id):
            assert gtf.source in [ "antisense", 
                                   "antisense_downstream", 
                                   "antisense_upstream", 
                                   "intergenic", 
                                   "sense_upstream" ], "Unrecognised lnc classification"
            lnc_all.add( gtf.gene_id )
            if gtf.source == "intergenic":
                lnc_intergenic.add( gtf.gene_id )
        else:
            raise Exception( "Unrecognised gene_id %s" % gtf.gene_id )

    # create temporary count files and add relevant count data
    tmpf_ref = P.getTempFilename("/ifs/scratch")
    tmpf_lnc_all = P.getTempFilename("/ifs/scratch")
    tmpf_lnc_intergenic = P.getTempFilename("/ifs/scratch")
    tmpf_1 = IOTools.openFile( tmpf_ref, "w" )
    tmpf_2 = IOTools.openFile( tmpf_lnc_all, "w" )
    tmpf_3 = IOTools.openFile( tmpf_lnc_intergenic, "w" )

    header = True
    for line in IOTools.openFile( count_file ).readlines():
        if header:
            tmpf_1.write( line )
            tmpf_2.write( line )
            tmpf_3.write( line )
            header = False
            continue
        gene_id = line.split()[0]
        if gene_id in refcoding:
            tmpf_1.write( line )
        elif gene_id in lnc_all:
            tmpf_2.write( line )
            if gene_id in lnc_intergenic:
                tmpf_3.write( line )
        else:
            raise Exception( "Gene id in count data is absent from the "
                             "lncRNA_refcoding.gtf" % gene_id )
    tmpf_1.close()
    tmpf_2.close()
    tmpf_3.close()


    # plot ordinations for refcoding genes
    # using P.snip on a glob expression for outfile works, I checked.
    outf_stub = P.snip( outfile, "intergenic_lncRNA_pca.rds" ) + "protein_coding"
    statement = ( "Rscript /ifs/devel/projects/proj010/transcriptome_pca_plots.R"
                  " %(tmpf_ref)s"
                  " 0.6"
                  " %(outf_stub)s" )
    to_cluster = False
    P.run()

    # plot ordinations for all lncRNAs
    outf_stub = P.snip( outfile, "intergenic_lncRNA_pca.rds" ) + "all_lncRNA"
    statement = ( "Rscript /ifs/devel/projects/proj010/transcriptome_pca_plots.R"
                  " %(tmpf_lnc_all)s"
                  " 0.6"
                  " %(outf_stub)s" )
    to_cluster = False
    P.run()

    # plot ordinations for all lncRNAs
    outf_stub = P.snip( outfile, "_pca.rds" ) 
    statement = ( "Rscript /ifs/devel/projects/proj010/transcriptome_pca_plots.R"
                  " %(tmpf_lnc_intergenic)s"
                  " 0.6"
                  " %(outf_stub)s" )
    to_cluster = False
    P.run()

    os.unlink( tmpf_ref )
    os.unlink( tmpf_lnc_all )
    os.unlink( tmpf_lnc_intergenic )
    E.info( "Completed plotting ordinations for %s" % count_file )    


#################################################################################
## subsection: plot refcoding and lncRNA venn diagrams of cuffdiff fpkms
#################################################################################
@jobs_limit( 1, "RGlobalEnv" )
@follows( mkdir("expression_venn_plots") )
@subdivide( loadReadCounts, 
            regex( "(.+)_median_fpkms.load" ),
            add_inputs( summarizelncRNARefcodingBiotypes ),
            r"expression_venn_plots/lncRNA_refcoding_venn_*.pdf" )
def plotVennForCuffdiffData( infiles, outfile ):
    """
    Plot venn diagrams of the number of genes that are co-expressed between
    cell types at a given fpkm threshold. 
    """
    count_table, biotype_table = [ P.snip(os.path.basename(x), ".load") 
                                   for x in infiles ]
    outfile_stub = "expression_venn_plots/lncRNA_refcoding_venn"

    # plot venn diagrams for i) protein coding genes ii) all lncRNA 
    # iii) intergenic lncRNA
    biotypes = ["protein_coding", "intergenic", "lncRNA_all"] 

    # plot venn diagrams for overlap between various B cell combinations
    combinations = [ [ "apro", "cpre", "dimmature", "bmature" ],
                     [ "gfollicular", "hmarginal", "eb1a" ], 
                     [ "gfollicular", "hmarginal", "fgerminal" ], 
                     [ "gfollicular", "fgerminal" ],
                     [ "gfollicular", "hmarginal", "eb1a", "fgerminal" ],
                     [ "apro", "cpre", "gfollicular", "fgerminal" ] ]


    for biotype in biotypes:
        P10.plotVennDiagrams( count_table, 
                              biotype_table, 
                              outfile_stub, 
                              biotype, 
                              combinations,
                              1.0 )


#################################################################################
## subsection: summarize expression values
#################################################################################
@follows(mkdir("expression_heatmaps"))
@transform( [ extractPerSampleFPKMs,
              logTransformCuffdiffFPKMs,
              logTransformCuffdiffFPKMs_smallconstant,
              extractMeanFPKMs,
              extractMedianFPKMs,
              extractMeanLogTransformedFPKMs,
              extractMedianLogTransformedFPKMs,
              fpmTransformData,
              extractMeanfpmTransformedCountData, 
              extractMedianfpmTransformedCountData,
              rlogTransformData, 
              extractMeanRlogTransformedCountData, 
              extractMedianRlogTransformedCountData ],
            regex( "(.+)/(.+)s.tsv.gz" ), 
            r"expression_heatmaps/\2_zscore.tsv.gz" )
def calcExpressionZScores( infile, outfile ):
    """
    Use pandas to calculate zscores
    """
    outf = P.snip( outfile, ".gz" )
    df = pd.read_table( infile, 
                        compression = "gzip", 
                        index_col = 0, 
                        header = 0, 
                        na_filter = True, 
                        na_values = ["NULL", "NA", "nan", "NaN"] )

    # calc zscore acoss rows
    df = df.apply( P10.calc_zscore, 
                   axis = 1 )

    df.to_csv( outf, sep = "\t", na_rep = "NaN" )

    # different transformation methods output different titles...
    to_cluster=False
    statement = ( "sed -i 's/Bcell-//g' %(outf)s; gzip %(outf)s" )
    P.run()


@subdivide( calcExpressionZScores, 
        regex( "(.+)/lncRNA_refcoding_(.+)(mean|median)(.+).tsv.gz" ),
        r"\1/*_\2\3\4_heatmap.png" )
def plotExpressionZScoreHeatmaps( infile, outfiles ):
    """
    Plot heatmap.2 heatmaps for per-gene zscores across samples
    """
    out_dir = os.path.dirname( infile )
    outf_stub = os.path.basename( infile )[len("lncRNA_refcoding_"):]
    outf_stub = P.snip( outf_stub, ".tsv.gz" ) + "_heatmap.png"


    job_options = "-l mem_free=10G"
    statement = ( "Rscript /ifs/devel/projects/proj010/transcriptome_heatmaps.R"
                  " %(infile)s"
                  " %(out_dir)s"
                  " %(outf_stub)s" )
    P.run()

#################################################################################
## subsection: calculate tissue specificity for each gene
#################################################################################
# extract pandas dataframe of mean expression values from csvdb
@follows( mkdir( "./expression_tissue_specificity" ) )
@transform( [ extractMeanFPKMs, 
              extractMedianFPKMs, 
              extractMedianLogTransformedFPKMs,
              extractMeanLogTransformedFPKMs,
              extractMeanRlogTransformedCountData,
              extractMedianRlogTransformedCountData ],
            regex( "(.+)/(.+).tsv.gz" ), 
            r"./expression_tissue_specificity/\2_tissueSpecificity.tsv.gz" )
def calcTissueSpecificity( infile, outfile ):
#    out_name = P.snip( os.path.basename( outfile ), ".tsv.gz" )
    out_name = "specificity"
    out_stub = P.snip( outfile, ".gz" )

    # open infile as pandas dataframe
    df = pd.read_table( infile, 
                        compression = "gzip", 
                        index_col = 0, 
                        header = 0, 
                        na_filter = True, 
                        na_values = ["NULL", "NA", "nan", "NaN"] )

    # reject genes that have data for fewer than 5 cell stages
    min_data_limit = 5
    df = df[ df.apply( P10.drop_nan, 
                       axis = 1, 
                       data_thresh = min_data_limit, 
                       rejected = False ) ]

    # calculate specificity
    specificity = df.apply( P10.calcSpecificity, axis = 1 )
    specificity.name = out_name
    
    # write to outfile
    specificity.to_csv( out_stub, 
                        index = True, 
                        sep = "\t", 
                        na_rep = "NA", 
                        header = True )
    to_cluster = False
    statement = ( "gzip %(out_stub)s" )
    P.run()


@collate( calcTissueSpecificity, 
          regex( "(.+)/(.+)_tissueSpecificity.tsv.gz" ),
          r"\1/lncRNA_refcoding_tissue_specificity.tsv.gz" )
def combineTissueSpecificity( infiles, outfile ):
    """
    Merge all tissue specificity estimates
    """
    infiles = " ".join( infiles )
    statement = ( "python %(scriptsdir)s/combine_tables.py"
                  "  --columns=1"
                  "  --log %(outfile)s.log"
                  "  --add-file-prefix"
                  "  --regex-filename '.*lncRNA_refcoding_(.+)_tissueSpecificity.tsv.gz'"
                  " %(infiles)s |"
                  " gzip > %(outfile)s" )
    P.run()


@transform( combineTissueSpecificity,
            regex( "(?:.+)/(.+).tsv.gz" ),
            r"./\1.load" )
def loadTissueSpecificity( infile, outfile ):
    P.load( infile, outfile, options="--add-index=gene_id" )


#################################################################################
## subsection: calculate zero adjusted tissue specificity for each gene
#################################################################################
# extract pandas dataframe of mean expression values from csvdb
# Its worth noting that the rlog transformation of data with all zero values
# returns all zero values, but rlog tranformation of a gene with one read returns
# negative values.
# For this reason I've removed all genes for which there are no read counts. 
# These should not really be in consideration for this study anyway. 

@follows( mkdir( "./expression_tissue_specificity" ) )
# @transform( [ extractMeanFPKMs, 
#               extractMedianFPKMs, 
#               extractMedianLogTransformedFPKMs,
#               extractMeanLogTransformedFPKMs,
#               extractMeanRlogTransformedCountData,
#               extractMedianRlogTransformedCountData ],
@transform( extractMedianRlogTransformedCountData,
            regex( "(.+)/(.+).tsv.gz" ), 
            r"./expression_tissue_specificity/\2_zadj_tissueSpecificity.tsv.gz" )
def calcTissueSpecificity_zadj( infile, outfile ):
    """
    Zero adjust transformed datavalues... see conversation thread
    https://support.bioconductor.org/p/59369/
    WARNING: The discarding genes with all zeros is based on the observation 
    that rlog tranformed data with small # of read counts are negative, while 
    rlog transformed data with all zeros are zero.
    """
#    out_name = P.snip( os.path.basename( outfile ), ".tsv.gz" )
    out_name = "specificity"

    # open infile as pandas dataframe
    df = pd.read_table( infile, 
                        compression = "gzip", 
                        index_col = 0, 
                        header = 0, 
                        na_filter = True, 
                        na_values = ["NULL", "NA", "nan", "NaN"] )

    # reject genes that have data for fewer than 5 cell stages
    min_data_limit = 5
    df = df[ df.apply( P10.drop_nan, 
                       axis = 1, 
                       data_thresh = min_data_limit, 
                       rejected = False ) ]

    # discard genes that have no expression across any cell stage. 
    df = df.loc[(df != 0).any(axis=1)]

    # If the dataframe contains negative values (due to transformation), 
    # then zero-adjust all the data by the minimum value.
    # First find the minimum of the minimum value in each column...
    df_min = min(df.min(0))
    if df_min < 0.0:
        E.warn("%s contains negative values... adjusting by %s" % (infile, str(df_min)))
        df = df.apply(lambda x: x + abs(df_min))

    # calculate specificity
    specificity = df.apply( P10.calcSpecificity, axis = 1 )
    specificity.name = out_name
    
    # write to outfile
    specificity.to_csv( IOTools.openFile(outfile, "w"), 
                        index = True, 
                        sep = "\t", 
                        na_rep = "NA", 
                        header = True )


@collate( calcTissueSpecificity_zadj, 
          regex( "(.+)/(.+)_zadj_tissueSpecificity.tsv.gz" ),
          r"\1/lncRNA_refcoding_zadj_tissue_specificity.tsv.gz" )
def combineTissueSpecificity_zadj( infiles, outfile ):
    """
    Merge all tissue specificity estimates
    """
    infiles = " ".join( infiles )
    statement = ( "python %(scriptsdir)s/combine_tables.py"
                  "  --columns=1"
                  "  --log %(outfile)s.log"
                  "  --add-file-prefix"
                  "  --regex-filename '.*lncRNA_refcoding_(.+)_zadj_tissueSpecificity.tsv.gz'"
                  " %(infiles)s |"
                  " gzip > %(outfile)s" )
    P.run()


@transform( combineTissueSpecificity_zadj,
            regex( "(?:.+)/(.+).tsv.gz" ),
            r"./\1.load" )
def loadTissueSpecificity_zadj( infile, outfile ):
    P.load( infile, outfile, options="--add-index=gene_id" )


#################################################################################
## subsection: run DESeq2 for raw read count data. 
#################################################################################
@follows( mkdir( "expression_DESeq_pairwise_tests" ) )
@transform( summarizeFeatureCounts, 
            regex( "(.+)/(.+)_raw_counts.tsv.gz" ),
            add_inputs( generateDESeq2DesignFile ),
            r"expression_DESeq_pairwise_tests/\2_results.rds" )
def runDESeq2( infiles, outfile ):
    """
    Runs DESeq2 with single factor design (cell type) and outputs an RDS 
    containing the DESeqDataSet object for querying downstream.
    """
    count_table = os.path.abspath( infiles[0] )
    design_file = os.path.abspath( infiles[1] )
    statement = ( "Rscript /ifs/devel/projects/proj010/lncRNA_refcoding_deseq2.R"
                   " %(count_table)s"
                   " %(design_file)s"
                   " %(outfile)s" % locals() )
    P.run()

@jobs_limit( 1 )
@split( runDESeq2, 
        regex( "(.+)/(.+)_results.rds" ), 
        add_inputs( generateDESeq2DesignFile ),
        r"\1/\2_results_*.tsv.gz" )
def extractDESeqPairwiseComparisons( infiles, outfile ):
    """
    Iterate through pairwise combinations of Bcells and extract results
    of wald test between pairs. 
    """
    result_file, design_file = infiles   
    # get groups from infile
    groups = set()
    for line in IOTools.openFile( design_file ).readlines():
        if line.startswith("sample"):
            continue
        groups.add( line.split()[1] )
    groups = sorted ( [ x for x in groups ] )

    for combination in itertools.combinations( groups, 2 ):
        numerator, denominator = combination
        E.info( "Fetching result for %s vs %s" % ( numerator, denominator ) )
        outf_stub = P.snip( result_file, ".rds" ) + "_" + "_vs_".join( combination )
        P10.fetchPairwiseResults( result_file, numerator, denominator, outf_stub )
        

@jobs_limit( 1 )
@collate( extractDESeqPairwiseComparisons, 
          regex( "(.+)/(.+)_results_(.+).tsv.gz" ), 
          r"\1/\2_DESeq_results.tsv.gz" )
def combineDESeqPairwiseComparisons( infiles, outfile ):
    """
    """
    infiles = " ".join( infiles )
    to_cluster=False
    tmpf = P.getTempFilename( "/ifs/scratch" )
    statement = ( "python %(scriptsdir)s/combine_tables.py"
                  " --cat=comparison"
                  " --missing-value=NA"
                  " --log=%(outfile)s.log"
                  " %(infiles)s"
                  " > %(tmpf)s" )
    P.run()

    outf = IOTools.openFile( outfile, "w" )
    for line in IOTools.openFile( tmpf ).readlines():
        if line.startswith("comparison"):
            outf.write( line )
            continue
        line = line.split()
        line[0] = P.snip( os.path.basename(line[0]), ".tsv.gz" )
        outf.write( "\t".join( line ) + "\n" )
    outf.close()

    os.unlink( tmpf )


@jobs_limit( 1 )
@transform( combineDESeqPairwiseComparisons,
            suffix( ".tsv.gz" ),
            "_padj.tsv.gz" )
def correctDESeqPairwiseComparisons( infile, outfile ):
    """
    Apply a global correction for pvalues across all pairwise comparisons.
    """
    P10.correctPairwiseComparisons( infile, outfile )


@transform( correctDESeqPairwiseComparisons, 
            regex( "(.+)/(.+).tsv.gz" ),
            r"\2.load" )
def loadDESeqPairwiseComparisons( infile, outfile ):
    P.load( infile, outfile, options="--add-index=gene_id" )


@follows( loadDESeqPairwiseComparisons )
def runDESeq2ForCountData():
    pass


#################################################################################
## subsection: bin genes based on overall expression value
#################################################################################
# This is done separately for all transformed data because rlog transformation
# can alter the order of relative expression, not sure about FPKM calculations
@follows( mkdir( "expression_bins" ) )
@transform( [ summarizeFeatureCounts, 
              rlogTransformData, 
              extractPerSampleFPKMs,
              logTransformCuffdiffFPKMs ],
            regex( "(.+)/(.+).tsv.gz" ),
            r"expression_bins/\2_median.tsv.gz" )
def calcMedianAbsoluteExpression( infile, outfile ):
    """
    Write a table containing the median expression value across all samples
    """
    df = pd.read_table( infile, 
                        header = 0, 
                        index_col = 0, 
                        sep = "\t", 
                        compression = "gzip" )
    # remove genes with zero expression across all samples
    # sanity check...  no lncRNAs with sum expression zero
    df_empty = df[df.sum(axis=1) == 0]
    non_expressed = [str(x) for x in df_empty.index]
    E.info("There are %i genes with zero expression acros samples" % len(non_expressed))
    non_expressed = [x for x in non_expressed if x.startswith("LNC")]
    assert len(non_expressed) == 0, "LncRNAs with zero expression across all samples"

    df = df[df.sum(axis=1) != 0]
  
    # get median
    df = df.apply( P10.robust_median, axis = 1 )
    df.index.name = "gene_id"
    df.columns = ["median",]
    outf = IOTools.openFile( outfile, "w" )
    df.to_csv( outf, sep = "\t", header = ["median",] )

 
@subdivide( calcMedianAbsoluteExpression,
            regex( "(.+).tsv.gz" ),
            r"\1.bin*" )
def binMedianAbsoluteExpression( infile, outfiles ):
    """
    Bin lncRNA and refcoding genes into 20 percentile bins, based on expression.
    """
    outfile_stub = P.snip( infile, ".tsv.gz" ) + "\.bin"
    outfile_log = P.snip( outfile_stub, ".bin" ) + ".log"
    statement = ( "zcat %(infile)s |"
                  " python %(scriptsdir)s/data2bins.py"
                  "  --log %(outfile_log)s"
                  "  --column 2"
                  "  --num-bins 4"
                  "  --output-filename-pattern %(outfile_stub)s%%i" )
    P.run()


@collate( binMedianAbsoluteExpression,
          regex( "(.+).bin[0-9]$" ),
          r"\1.binned.tsv.gz" )
def combineMedianAbsoluteExpressionBins( infiles, outfile ):
    """
    Merge sample median expression bins into a single file.
    """
    infiles = " ".join( infiles )
    statement = ( "python %(scriptsdir)s/combine_tables.py"
                  "  --cat Bin"
                  "  --log %(outfile)s.log"
                  "  --regex-filename='.+\.(bin[0-9])'"
                  " %(infiles)s |"
                  " gzip > %(outfile)s" )
    P.run()
    

@collate( combineMedianAbsoluteExpressionBins,
         regex( "(.+)/lncRNA_refcoding(.+)_median.binned.tsv.gz" ),
          r"\1/lncRNA_refcoding_median_expression_bins.tsv.gz" )
def combineMedianExpressionBinsAcrossCountMethods( infiles, outfile ):
   """  
   Merge all binned expression values into one table
   """
   infiles = " ".join( infiles )
   statement = ( "python %(scriptsdir)s/combine_tables.py"
                 " --columns=2"
                 " --log %(outfile)s.log"
                 " --add-file-prefix"
                 " --regex-filename '.*lncRNA_refcoding_(.+)_median.binned.tsv.gz'"
                 " %(infiles)s |"
                 " gzip > %(outfile)s" )
   P.run()


@transform( combineMedianExpressionBinsAcrossCountMethods,
            regex( ".+/(.+).tsv.gz" ),
            r"\1.load" )
def loadMedianExpressionBins( infile, outfile ):
    P.load( infile, outfile, options="--add-index=gene_id" )

################################################################################
## subsection: create pairwise DE summaries
################################################################################
@jobs_limit( 1 )
@split( loadDESeqPairwiseComparisons,
        regex("(lncRNA)_(refcoding)_DESeq_results_padj.load"),
        add_inputs( generateDESeq2DesignFile,
                    summarizelncRNARefcodingBiotypes ),
        [ r"expression_DESeq_pairwise_tests/\1_pairwise_test_summary_abs_count",
          r"expression_DESeq_pairwise_tests/\1_pairwise_test_summary_dir_count",
          r"expression_DESeq_pairwise_tests/\1_pairwise_test_summary_abs_percent",
          r"expression_DESeq_pairwise_tests/\1_pairwise_test_summary_dir_percent",
          r"expression_DESeq_pairwise_tests/\2_pairwise_test_summary_abs_count",
          r"expression_DESeq_pairwise_tests/\2_pairwise_test_summary_dir_count",
          r"expression_DESeq_pairwise_tests/\2_pairwise_test_summary_abs_percent",
          r"expression_DESeq_pairwise_tests/\2_pairwise_test_summary_dir_percent" ] )
def summarizePairwiseDEResults(infiles, outfile):
    de_table, design_table, biotype = infiles
    outdir = "expression_DESeq_pairwise_tests"
    outstub = "_pairwise_test_summary"

    # get conditions from design file
    groups = set()
    for line in IOTools.openFile( design_table ).readlines():
        if line.startswith("sample"):
            continue
        groups.add(line.split()[1])
    groups = sorted( [x for x in groups] )


    # extract all DE data from table. 
    de_table = P.snip( os.path.basename( de_table ), ".load" )
    statement = ("SELECT gene_id,log2FoldChange,numerator,denominator,padj_global"
                 " FROM %s" % de_table)
    df_orig = PU.fetch_DataFrame( statement )

    # For the protein coding and lncRNA genes separately...
    biotype_table = P.snip(os.path.basename(biotype), ".load")
    for gene_type in [ "refcoding", "lncRNA" ]:
        if gene_type == "refcoding":
            statement = "SELECT gene_id FROM %s WHERE biotype = '%s'" % (biotype_table, "protein_coding")
        else:
            statement = "SELECT gene_id FROM %s WHERE biotype != '%s'" % (biotype_table, "protein_coding")
        biotype_gene_ids = [str(x[0]) for x in PU.fetch(statement)]

        # subset the required gene_type data from the original dataframe
        df = df_orig[df_orig["gene_id"].isin( biotype_gene_ids )]
        n_genes = len(set(df["gene_id"]))

        # create dataframe to be output...
        columns = ["pro", 
                   "pre", 
                   "immature", 
                   "mature", 
                   "follicular", 
                   "marginal", 
                   "b1a", 
                   "germinal"]
        df_abs = pd.DataFrame(columns=columns, index=columns)

        # fill output dataframe in non-directional manner
        for group in itertools.combinations( groups, 2 ):
            numerator, denominator = group
            de_genes = df[ (df.padj_global < 0.01) 
                           & (df.numerator==numerator) 
                           & (df.denominator==denominator) 
                           & (abs(df.log2FoldChange) > 1) ]
            assert len(de_genes.index) != 0
            assert len( df[ (df.padj_global < 0.01) 
                            & (df.numerator==denominator) 
                            & (df.denominator==numerator) 
                            & (abs(df.log2FoldChange) > 1) ] ) == 0
            de_genes = len(de_genes.index)
            df_abs.loc[numerator, denominator] = de_genes
            df_abs.loc[denominator, numerator] = de_genes
        # fill in centre
        for i in groups:
            df_abs.loc[i, i] = 0

        # write undirectional outfile of counts
        outf = os.path.join( outdir, gene_type + outstub + "_abs_count" )
        df_abs.to_csv(outf, index=True, index_label="cell_type", sep="\t")
        
        # write undirectional outfile of percentages
        df_per = (df_abs/n_genes*100)
        outf = os.path.join( outdir, gene_type + outstub + "_abs_percent" )
        df_per.to_csv(outf, index=True, index_label="cell_type", sep="\t")

        
        # fill output dataframe in a directional manner
        df_abs_dir = pd.DataFrame(columns=columns, index=columns)
        for group in itertools.combinations( groups, 2 ):
            numerator, denominator = group
            de_genes_up = df[ (df.padj_global < 0.01) 
                              & (df.numerator==numerator) 
                              & (df.denominator==denominator) 
                              & (df.log2FoldChange > 1) ]
            de_genes_down = df[ (df.padj_global < 0.01) 
                                & (df.numerator==numerator) 
                                & (df.denominator==denominator) 
                                & (df.log2FoldChange < -1) ]
            # assert len(de_genes_up.index) != 0
            de_genes_up = len(de_genes_up.index)
            de_genes_down = len(de_genes_down.index)
            df_abs_dir.loc[numerator, denominator] = de_genes_up
            df_abs_dir.loc[denominator, numerator] = de_genes_down
        # fill centre
        for i in groups:
            df_abs_dir.loc[i, i] = 0

        # write directional outfile of counts
        outf = os.path.join( outdir, gene_type + outstub + "_dir_count" )
        df_abs_dir.to_csv(outf, index=True, index_label="cell_type", sep="\t")

        # write directional outfile of counts
        df_per_dir = (df_abs_dir/n_genes*100)
        outf = os.path.join( outdir, gene_type + outstub + "_dir_percent" )
        df_per_dir.to_csv(outf, index=True, index_label="cell_type", sep="\t")


#################################################################################
## subsection: create pairwise DE summaries for expression bins
#################################################################################
@jobs_limit( 1 )
@subdivide( binMedianAbsoluteExpression,
            regex( "(.+)/(lncRNA)_(refcoding)_rlog_trans_counts_median.bin([0-9])" ),
            add_inputs( loadDESeqPairwiseComparisons,
                        generateDESeq2DesignFile,
                        summarizelncRNARefcodingBiotypes ),
            [ r"expression_DESeq_pairwise_tests/\2_pairwise_test_summary_bin\4_abs_count",
              r"expression_DESeq_pairwise_tests/\2_pairwise_test_summary_bin\4_dir_count",
              r"expression_DESeq_pairwise_tests/\2_pairwise_test_summary_bin\4_abs_percent",
              r"expression_DESeq_pairwise_tests/\2_pairwise_test_summary_bin\4_dir_percent",
              r"expression_DESeq_pairwise_tests/\3_pairwise_test_summary_bin\4_abs_count",
              r"expression_DESeq_pairwise_tests/\3_pairwise_test_summary_bin\4_dir_count",
              r"expression_DESeq_pairwise_tests/\3_pairwise_test_summary_bin\4_abs_percent",
              r"expression_DESeq_pairwise_tests/\3_pairwise_test_summary_bin\4_dir_percent" ] )
def summarizePairwiseDEResultsAcrossBins( infiles, outfile ):
    """
    Take the 4 bins based on cuffdiff FPKM value. 
    For each, extract the genes that are differentially expressed, 
    in that bin. Done for i) protein coding genes and ii) all lncRNAs. 
    """
    bin_genes, de_table, design_table, biotype = infiles
    bin_number = bin_genes[-1]
    outdir = "expression_DESeq_pairwise_tests"
    outstub = "_pairwise_test_summary_bin" + bin_number

    # get conditions from design file
    groups = set()
    for line in IOTools.openFile( design_table ).readlines():
        if line.startswith("sample"):
            continue
        groups.add(line.split()[1])
    groups = sorted( [x for x in groups] )

    # extract all DE data from table. 
    de_table = P.snip( os.path.basename( de_table ), ".load" )
    statement = ("SELECT gene_id,log2FoldChange,numerator,denominator,padj_global"
                 " FROM %s" % de_table)
    df_orig = PU.fetch_DataFrame( statement )

    # subset the dataframe by the values in the bin file
    binned_gene_ids = []
    for line in IOTools.openFile( bin_genes ):
        if line.startswith( "gene_id" ): continue
        line = line.split()
        # ignore zero values... taken out bc zero vals were removed when binning
        # if int(line[1]) > 0:
        #     binned_gene_ids.append( line[0] )
        binned_gene_ids.append(line[0])
    df_orig = df_orig[df_orig["gene_id"].isin( binned_gene_ids )]

    # For the protein coding and lncRNA genes separately...
    biotype_table = P.snip(os.path.basename(biotype), ".load")
    for gene_type in [ "refcoding", "lncRNA" ]:
        if gene_type == "refcoding":
            statement = "SELECT gene_id FROM %s WHERE biotype = '%s'" % (biotype_table, "protein_coding")
        else:
            statement = "SELECT gene_id FROM %s WHERE biotype != '%s'" % (biotype_table, "protein_coding")
        biotype_gene_ids = [str(x[0]) for x in PU.fetch(statement)]

        # subset the required gene_type data from the original dataframe
        df = df_orig[df_orig["gene_id"].isin( biotype_gene_ids )]
        n_genes = len(set(df["gene_id"]))

        # print "Bin: %s" % bin_genes
        # print gene_type
        # print "ngenes %i" % n_genes
        # print "\n"

        # create dataframe to be output...
        columns = ["pro", 
                   "pre", 
                   "immature", 
                   "mature", 
                   "follicular", 
                   "marginal", 
                   "b1a", 
                   "germinal"]
        df_abs = pd.DataFrame(columns=columns, index=columns)

        # fill output dataframe in non-directional manner
        for group in itertools.combinations( groups, 2 ):
            numerator, denominator = group
            de_genes = df[ (df.padj_global < 0.01) 
                           & (df.numerator==numerator) 
                           & (df.denominator==denominator) 
                           & (abs(df.log2FoldChange) > 1) ]
            assert len(de_genes.index) != 0
            assert len( df[ (df.padj_global < 0.01) 
                            & (df.numerator==denominator) 
                            & (df.denominator==numerator) 
                            & (abs(df.log2FoldChange) > 1) ] ) == 0
            de_genes = len(de_genes.index)
            df_abs.loc[numerator, denominator] = de_genes
            df_abs.loc[denominator, numerator] = de_genes
        # fill in centre
        for i in groups:
            df_abs.loc[i, i] = 0

        # write undirectional outfile of counts
        outf = os.path.join( outdir, gene_type + outstub + "_abs_count" )
        df_abs.to_csv(outf, index=True, index_label="cell_type", sep="\t")
        
        # calculate total number of genes in bin file
        E.info( "Bin_file %s has %s %s genes" % ( str(bin_number), 
                                                  str(n_genes), 
                                                  gene_type) )

        # write undirectional outfile of percentages
        df_per = (df_abs/n_genes*100)
        outf = os.path.join( outdir, gene_type + outstub + "_abs_percent" )
        df_per.to_csv(outf, index=True, index_label="cell_type", sep="\t")

        
        # fill output dataframe in a directional manner
        df_abs_dir = pd.DataFrame(columns=columns, index=columns)
        for group in itertools.combinations( groups, 2 ):
            numerator, denominator = group
            de_genes_up = df[ (df.padj_global < 0.01) 
                              & (df.numerator==numerator) 
                              & (df.denominator==denominator) 
                              & (df.log2FoldChange > 1) ]
            de_genes_down = df[ (df.padj_global < 0.01) 
                                & (df.numerator==numerator) 
                                & (df.denominator==denominator) 
                                & (df.log2FoldChange < -1) ]
            # assert len(de_genes_up.index) != 0
            de_genes_up = len(de_genes_up.index)
            de_genes_down = len(de_genes_down.index)
            df_abs_dir.loc[numerator, denominator] = de_genes_up
            df_abs_dir.loc[denominator, numerator] = de_genes_down
        # fill centre
        for i in groups:
            df_abs_dir.loc[i, i] = 0

        # write directional outfile of counts
        outf = os.path.join( outdir, gene_type + outstub + "_dir_count" )
        df_abs_dir.to_csv(outf, index=True, index_label="cell_type", sep="\t")

        # write directional outfile of counts
        df_per_dir = (df_abs_dir/n_genes*100)
        outf = os.path.join( outdir, gene_type + outstub + "_dir_percent" )
        df_per_dir.to_csv(outf, index=True, index_label="cell_type", sep="\t")


@transform(summarizePairwiseDEResultsAcrossBins,
           regex( "(.+)_(count|percent)" ),
#           regex("(.+)/lncRNA_pairwise_test_summary_bin0_abs_count" ),
           r"\1_\2.png" )
#           r"\1/lncRNA_pairwise_test_summary_bin0_abs_count.png" )
def plotPairwiseDEResultsAcrossBins( infile, outfile ):
    """
    Read in matrix and plot as heatmap using ggplot
    """
    if re.search( os.path.basename("count"), infile ):
        lim = 3000
    else:
        lim = 100
    P10.plotHeatmapsggplot( infile, outfile, lim, submit=True )


#################################################################################
## subsection: calculate corrleation between lncRNAs protein_coding genes
#################################################################################
@follows( mkdir( "./correlation_lncRNA_vs_protein" ) )
@split( [ extractPerSampleFPKMs, logTransformCuffdiffFPKMs, rlogTransformData ],
        regex( "(.+)/(.+)_(fpkms|counts).tsv.gz" ), 
        [ r"./correlation_lncRNA_vs_protein/\2_\3_per_sample_pr_correlation_stats.tsv.gz", 
          r"./correlation_lncRNA_vs_protein/\2_\3_per_sample_pr_correlation_summary.tsv.gz" ] )
def calculateLncRNAProteinCodingPearsonCorrelations( infile, outfiles ):
    outf_1, outf_2 = outfiles
    method = "pearson"
    min_obs = str(10)
    params = [ infile, outf_1, outf_2, method, min_obs, "ENSMUSG", "LNC", "True" ]
    log = P.snip( outf_1, "_stats.tsv.gz" ) + ".log"
    P.submit( "/ifs/devel/projects/proj010/PipelineProj010", 
              "calculateCorrelations",
              params = params,
              logfile = log,
              jobOptions = "-l mem_free=30G" )


@transform( calculateLncRNAProteinCodingPearsonCorrelations, 
            regex( "(.+)/(.+)_stats.tsv.gz" ),
            r"\2_stats.load" )
def loadLncRNAProteinCodingPearsonCorrelations( infile, outfile ):
    to_cluster = True
    job_options = "-l mem_free=30G"
    P.load( infile, outfile, options= "-i Gene_id_1 -i Gene_id_2" )


@transform( calculateLncRNAProteinCodingPearsonCorrelations, 
            regex( "(.+)/(.+)_stats.tsv.gz" ),
            add_inputs( findNearestProteinCodingGenes ),
            r"\1/\2_nn_stats.tsv.gz" )
def findLncRNAProteinCodingNNPearsonCorrelations( infiles, outfile ):
    corr_tsv, nn_tsv = infiles

    neighbours_dict = {}
    for line in IOTools.openFile( nn_tsv ).readlines():
        if line.startswith("gene_id"):
            line = line.split()
            assert line[1] == "closest_id", "Infile not in expected order"
            assert line[4] == "id5", "Infile not in expected order"
            assert line[7] == "id3", "Infile not in expected order"
            continue
        else:
            line = line.split()
            lnc_id = line[0]
            # get nearest, nearest5, nearest3
            neighbours = [ [line[1],], [line[4],], [line[7],] ]
            # append "na" as correlation val for missing neigbours
            for i in neighbours:
                if i[0] == "na":
                    i.append("na")
            neighbours_dict[ lnc_id ] = neighbours

    for line in IOTools.openFile( corr_tsv ).readlines():
        if line.startswith( "Gene" ):
            line = line.split()
            assert line[0] == "Gene_id_1",  "Infile not in expected order" 
            assert line[1] == "Gene_id_2",  "Infile not in expected order"
            assert line[5] == "coefficient", "Infile not in expected order"
            continue
        else:
            line = line.split()
            lnc_id = line[0]
            pc_id = line[1]
            coeff = line[5]
            
            # get nested list of neighbours
            neighbours = neighbours_dict[ lnc_id ]
            if pc_id in [x[0] for x in neighbours]:
                for i in neighbours:
                    if i[0] == pc_id:
                        i.append(coeff)
                neighbours_dict[ lnc_id ] = neighbours
            else:
                continue

    outf = IOTools.openFile( outfile, "w" )
    outf.write( "lnc_id\tidClosest\tcorrClosest\t"
                "id5\tcorr5\tid3\tcorr3\tcorrMax\n" )
    for lncRNA, correlations in neighbours_dict.iteritems():
        idClosest = correlations[0][0]
        try:
            corrClosest = correlations[0][1]
        except IndexError:
            print "\n"
            print lncRNA
            print correlations
            print "\n"
            break
        id5 = correlations[1][0]
        corr5 = correlations[1][1]
        id3 = correlations[2][0]
        corr3 = correlations[2][1]
        corrs = [x[1] for x in correlations if x[1] != "na"]
        # correlation is also NA if refcoding has 0 expression across all samples
        if len(corrs) == 0:
            corrMax = "na"
        else:
            corrMax = str(max(map(float, corrs)))
        
        outf.write( "\t".join([lncRNA,
                               idClosest, 
                               corrClosest, 
                               id5, 
                               corr5, 
                               id3, 
                               corr3, 
                               corrMax]) + "\n" )
    outf.close()

    # outf = IOTools.openFile( outfile, "w" )
    # outf.write( "lnc_id\tidClosest\tcorrClosest\t"
    #             "id5\tcorr5\tid3\tcorr3\tcorrMax\n" )
    # for lncRNA, correlations in neighbours_dict.iteritems():
    #     idClosest = correlations[0][0]
    #     corrClosest = correlations[0][1]
    #     id5 = correlations[1][0]
    #     corr5 = correlations[1][1]
    #     id3 = correlations[2][0]
    #     corr3 = correlations[2][1]
    #     corrs = [x[1] for x in correlations if x[1] != "na"]
    #     corrMax = str(max(map(float, corrs)))

    #     outf.write( "\t".join([lncRNA,
    #                            idClosest, 
    #                            corrClosest, 
    #                            id5, 
    #                            corr5, 
    #                            id3, 
    #                            corr3, 
    #                            corrMax]) + "\n" )
    # outf.close()


@transform( findLncRNAProteinCodingNNPearsonCorrelations, 
            regex( "(.+)/(.+)_stats.tsv.gz" ),
            r"\2_stats.load" )
def loadLncRNAProteinCodingNNPearsonCorrelations( infile, outfile ):
    job_options = "-l mem_free=30G"
    P.load( infile, outfile, options = "-i Gene_id_1" )


@follows( mkdir( "./correlation_lncRNA_vs_protein" ) )
@split( [ extractMedianFPKMs, extractMedianLogTransformedFPKMs, extractMedianRlogTransformedCountData ],
        regex( "(.+)/(.+)_median_(fpkms|counts|ls_fpkms|l_fpkms).tsv.gz" ), 
        [ r"./correlation_lncRNA_vs_protein/\2_\3_median_pr_correlation_stats.tsv.gz", 
          r"./correlation_lncRNA_vs_protein/\2_\3_median_pr_correlation_summary.tsv.gz" ] )
def calculateMedianLncRNAProteinCodingPearsonCorrelations( infile, outfiles ):
    outf_1, outf_2 = outfiles
    method = "pearson"
    min_obs = str(5)
    log = P.snip( outf_1, "_stats.tsv.gz" ) + ".log"
    params = [ infile, outf_1, outf_2, method, min_obs, "ENSMUSG", "LNC", "True" ]
    P.submit( "/ifs/devel/projects/proj010/PipelineProj010", 
              "calculateCorrelations",
              params=params,
              logfile=log,
              jobOptions = "-l mem_free=30G" )


@transform( calculateMedianLncRNAProteinCodingPearsonCorrelations, 
            regex( "(.+)/(.+)_stats.tsv.gz" ), 
            r"\2_stats.load" )
def loadMedianLncRNAProteinCodingPearsonCorrelations( infile, outfile ):
    P.load( infile, outfile, options= "-i Gene_id_1 -i Gene_id_2" )


@transform( calculateMedianLncRNAProteinCodingPearsonCorrelations, 
            regex( "(.+)/(.+)_stats.tsv.gz" ),
            add_inputs( findNearestProteinCodingGenes ),
            r"\1/\2_nn_stats.tsv.gz" )
def findMedianLncRNAProteinCodingNNPearsonCorrelations( infiles, outfile ):
    corr_tsv, nn_tsv = infiles
    outf = IOTools.openFile( outfile, "w" )

    neighbours_dict = {}
    for line in IOTools.openFile( nn_tsv ).readlines():
        if line.startswith("gene_id"):
            line = line.split()
            assert line[1] == "closest_id", "Infile not in expected order"
            assert line[4] == "id5", "Infile not in expected order"
            assert line[7] == "id3", "Infile not in expected order"
            continue
        else:
            line = line.split()
            lnc_id = line[0]
            # get nearest, nearest5, nearest3
            neighbours = [ [line[1],], [line[4],], [line[7],] ]
            # append "na" as correlation val for missing neigbours
            for i in neighbours:
                if i[0] == "na":
                    i.append("na")
            neighbours_dict[ lnc_id ] = neighbours

    for line in IOTools.openFile( corr_tsv ).readlines():
        if line.startswith( "Gene" ):
            line = line.split()
            assert line[0] == "Gene_id_1",  "Infile not in expected order" 
            assert line[1] == "Gene_id_2",  "Infile not in expected order"
            assert line[5] == "coefficient", "Infile not in expected order"
            continue
        else:
            line = line.split()
            lnc_id = line[0]
            pc_id = line[1]
            coeff = line[5]
            
            # get nested list of neighbours
            neighbours = neighbours_dict[ lnc_id ]
            if pc_id in [x[0] for x in neighbours]:
                for i in neighbours:
                    if i[0] == pc_id:
                        i.append(coeff)
                neighbours_dict[ lnc_id ] = neighbours
            else:
                continue

    outf = IOTools.openFile( outfile, "w" )
    outf.write( "lnc_id\tidClosest\tcorrClosest\t"
                "id5\tcorr5\tid3\tcorr3\tcorrMax\n" )
    for lncRNA, correlations in neighbours_dict.iteritems():
        idClosest = correlations[0][0]
        try:
            corrClosest = correlations[0][1]
        except IndexError:
            print "\n"
            print lncRNA
            print correlations
            print "\n"
            break
        id5 = correlations[1][0]
        corr5 = correlations[1][1]
        id3 = correlations[2][0]
        corr3 = correlations[2][1]
        corrs = [x[1] for x in correlations if x[1] != "na"]
        # correlation is also NA if refcoding has 0 expression across all samples
        if len(corrs) == 0:
            corrMax = "na"
        else:
            corrMax = str(max(map(float, corrs)))
        
        outf.write( "\t".join([lncRNA,
                               idClosest, 
                               corrClosest, 
                               id5, 
                               corr5, 
                               id3, 
                               corr3, 
                               corrMax]) + "\n" )
    outf.close()


@transform( findMedianLncRNAProteinCodingNNPearsonCorrelations, 
            regex( "(.+)/(.+)_stats.tsv.gz" ),
            r"\2_stats.load" )
def loadMedianLncRNAProteinCodingNNPearsonCorrelations( infile, outfile ):
    P.load( infile, outfile, options = "-i Gene_id_1 -i Gene_id_2" )
@follows(findLncRNAProteinCodingNNPearsonCorrelations,
         findMedianLncRNAProteinCodingNNPearsonCorrelations)
def calculateNearestNeighbourCorrelations():
    pass


@follows(loadLncRNAProteinCodingNNPearsonCorrelations,
         loadMedianLncRNAProteinCodingNNPearsonCorrelations)
def loadNearestNeighbourCorrelations():
    pass

#################################################################################
# Section: lncRNA conservation summary
#################################################################################
#################################################################################
## subsection: calculate phastcons score summaries for exons.
#################################################################################
@follows( mkdir("conservation_phastCons") )
@transform( combineRefcodingAndLncRNA, 
            regex( "(.+)/(.+).gtf.gz" ),
            add_inputs( PARAMS["location_phastcons"] ),
            r"conservation_phastCons/\2.phastCons_summary.tsv.gz" )
def extractTranscriptPhastConsSummary( infiles, outfile ):
    """
    Retrieve summary of PhastCons scores across transcripts
    """
    gtf_file, bigwig_file = infiles
    P10.getPhastConsScores( gtf_file, bigwig_file, outfile )


@transform( extractTranscriptPhastConsSummary,
            regex( "(.+)/(.+).tsv.gz" ),
            r"\2.load" )
def loadPhastConsSummary( infile, outfile ):
    P.load( infile, outfile, options="--add-index=gene_id" )

#################################################################################
# Section: additional lncRNA statistics
#################################################################################
# determineRefcodingAndLncRNAStatus - summarizes exon and coding status of genes
# calculateLncRNASummaryStats - produces summary stats for the four exon status
# calculateLncRNAGeneLengthStats - produces length stats for merged gene models
# calculateGeneCompositionStats - produces cpg composition for all transcripts
#################################################################################
## subsection: summarize lncRNA genomic positions
#################################################################################

@follows( mkdir( "./characterize_lncrna_refcoding" ) )
@transform( classifyMergedLncRNAs, 
            regex( "(.+)/(.+)_final.gtf.gz" ),
            r"./characterize_lncrna_refcoding/\2_classification.tsv.gz" )
def summarizeLncRNAGenomicPosition( infile, outfile ):
    """
    Iterate through flat gene models, assert all intervals in merged gene model
    have same source field, output tsv containing gene_id, classification.
    """
    outf = IOTools.openFile( outfile, "w" )
    outf.write( "gene_id\tclassification\n" )

    # output from mergeMEandSELncRNA is sorted by gene
    N_gtf = 0
    for gtfs in GTF.flat_gene_iterator( GTF.iterator( IOTools.openFile( infile ) ) ):
        N_gtf += 1
        source_ids = [ x.source for x in gtfs ]
        gene_id = gtfs[0].gene_id
        assert len( set( source_ids ) ) == 1, "Multiple classifications for %s:\t%s" % ( gene_id, " ".join(source_ids) )  
        classification = source_ids[0]
        outf.write( gene_id + "\t" + classification + "\n" )
    
    outf.close()
    E.info( "There are %i entries in gtf file" % N_gtf )


@transform( summarizeLncRNAGenomicPosition, 
            regex( "(.+)/(.+).tsv.gz" ), 
            r"\2.load" )
def loadLncRNAGenomicPositionSummary( infile, outfile ):
    P.load( infile, outfile, options="--add-index=gene_id" )


#################################################################################
## subsection: 
#################################################################################
@follows( mkdir( "./characterize_lncrna_refcoding" ) )
@transform( classifyMergedLncRNAs, 
            regex( "(.+)/(.+)_final.gtf.gz" ), 
            r"./characterize_lncrna_refcoding/\2_interval_stats.tsv.gz" )
def countLncRNAFeatures( infile, outfile ):
    """
    Run gff2stats on merged lncRNA set. 
    """
    statement = ( "zcat %(infile)s |"
                  " python %(scriptsdir)s/gff2stats.py"
                  "  --is-gtf"
                  "  --log=%(outfile)s.log |"
                  " gzip > %(outfile)s" )
    P.run()

@follows( mkdir( "./characterize_lncrna_refcoding" ) )
@transform( combineRefcodingAndLncRNA, 
            regex( "(.+)/(.+).gtf.gz" ),
            r"./characterize_lncrna_refcoding/\2_status.tsv.gz" )
def determineRefcodingAndLncRNAStatus( infile, outfile ):
    """
    A miscellaneous function for collating gene characteristics from the 
    pooled geneset. 
    Amongst other things, checks strand consistency within gene models. 
    """   
    outf = IOTools.openFile( outfile, "w" )
    outf.write( "gene_id\t"
                "gene_name\t"
                "strand\t"
                "coding\t"
                "n_transcripts\t"
                "n_exon_max\t"
                "n_exon_min\t"
                "n_exon_mean\t"
                "exon_status_lnc\n" )

    tmpf = P.getTempFilename('.')
    statement = ("zcat %(infile)s |"
                 " python %(scriptsdir)s/gtf2gtf.py"
                 " --method=sort --sort-order=gene+transcript"
                 " --log=%(outfile)s.log"
                 " > %(tmpf)s")
    to_cluster = False
    P.run()

    for gene in GTF.gene_iterator( GTF.iterator( IOTools.openFile( tmpf ) ) ):
        n_trans = 0
        n_exons = []
        strand = []
        gene_id = gene[0][0].gene_id
        if "gene_name" in gene[0][0].asDict().keys():
            gene_name = gene[0][0].gene_name
        else:
            gene_name = "NA"
        if gene[0][0].gene_id.startswith( "LNCGme" ):
            coding = "LNCme"
        elif gene[0][0].gene_id.startswith( "LNCGse" ):
            coding = "LNCse"
        else:
            coding = "Coding"
        if "exon_status" in gene[0][0].asDict().keys():
            exon_status_lnc = gene[0][0].asDict()[ "exon_status" ]
        else:
            exon_status_lnc = "pc"
        for transcript in gene:
            n_trans += 1
            n_exons.append( len( [ x for x in transcript ] ) )
            strand.extend( [ x.strand for x in transcript ] ) 
        if re.search( "sm", gene[0][0].gene_id ):
            n_exon_max = n_trans
            n_exon_min = 1
            n_exon_mean = "NA"
        elif  re.search( "mm", gene[0][0].gene_id ):
            n_exon_max = "NA"
            n_exon_min = "NA"
            n_exon_mean = "NA"
        else:
            n_exon_max = max( n_exons )
            n_exon_min = min( n_exons )
            n_exon_mean = np.mean( n_exons )
        strand = ",".join( set( strand ) )
        outf.write( "\t".join( map( str, [ gene_id, 
                                           gene_name, 
                                           strand, 
                                           coding, 
                                           n_trans, 
                                           n_exon_max, 
                                           n_exon_min, 
                                           n_exon_mean, 
                                           exon_status_lnc ] ) ) )
        outf.write( "\n" )
    outf.close()


# @transform( determineRefcodingAndLncRNAStatus, 
#             regex( "(.+)/(.+).tsv.gz" ), 
#             r"./\2.load" )
# def loadRefcodingAndLncRNAStatus( infile, outfile ):
#     P.load( infile, outfile, options="--add-index=gene_id" )


# @follows( mkdir( "./characterize_lncrna_refcoding" ) )
# @transform( combineRefcodingAndLncRNA, 
#             regex( "(.+)/(.+).gtf.gz" ),
#             r"./characterize_lncrna_refcoding/\2_summary_stats.tsv" )
# def calculateLncRNASummaryStats( infile, outfile ):
#     """
#     Produces summary stats for the lncRNAs remaining after filtering. Uses all
#     the counters present in PipelineLncRNA and counts for multi-exon and single
#     exon lncRNAs separately. Many of these counters are redundant for one or other
#     classification, but they are used as a sanity check. 
#     """
#     # open outfile and write headers
#     outf = IOTools.openFile( outfile, "w" )
#     outf.write( "\t".join( [ "lncrna_type", 
#                              "no_transcripts" ,
#                              "no_genes",
#                              "no_exons_per_transcript",
#                              "no_exons_per_gene",
#                              "no_single_exon_transcripts",
#                              "no_multi_exon_transcripts", 
#                              "no_single_exon_genes", 
#                              "no_multi_exon_genes\n" ] ) )

#     # split me and se lncrna into two sorted tempfiles
#     me_tmp = P.getTempFilename(".")
#     se_tmp = P.getTempFilename(".")
#     mm_tmp = P.getTempFilename(".")
#     sm_tmp = P.getTempFilename(".")
#     pc_tmp = P.getTempFilename(".")
#     statement = '''
#                 zcat %(infile)s |
#                  awk '$0 ~ /gene_id \"LNCGme/' | 
#                  python %(scriptsdir)s/gtf2gtf.py
#                   --method=sort --sort-order=gene
#                   --log=%(outfile)s.log
#                  > %(me_tmp)s;
#                 checkpoint;
#                 zcat %(infile)s |
#                  awk '$0 ~ /gene_id \"LNCGse/' | 
#                  python %(scriptsdir)s/gtf2gtf.py
#                   --method=sort --sort-order=gene
#                   --log=%(outfile)s.log
#                  > %(se_tmp)s
#                 checkpoint;
#                 zcat %(infile)s |
#                  awk '$0 ~ /gene_id \"LNCGsm/' | 
#                  python %(scriptsdir)s/gtf2gtf.py
#                   --method=sort --sort-order=gene
#                   --log=%(outfile)s.log
#                  > %(sm_tmp)s;
#                 checkpoint;
#                 zcat %(infile)s |
#                  awk '$0 ~ /gene_id \"LNCGmm/' |  
#                  python %(scriptsdir)s/gtf2gtf.py
#                   --method=sort --sort-order=gene
#                   --log=%(outfile)s.log
#                  > %(mm_tmp)s;
#                 checkpoint;
#                 zcat %(infile)s |
#                  awk '$0 ~ /gene_id \"ENSMUS/' |  
#                  python %(scriptsdir)s/gtf2gtf.py
#                   --method=sort --sort-order=gene
#                   --log=%(outfile)s.log
#                  > %(pc_tmp)s;

#                 ''' 
#     to_cluster = False
#     P.run()

#     # print "me_tmp : %s" % me_tmp
#     # print "se_tmp : %s" % se_tmp
#     # print "sm_tmp : %s" % sm_tmp
#     # print "mm_tmp : %s" % mm_tmp
#     # print "pc_tmp : %s" % pc_tmp

#     # calculate summary stats for all counters in PipelineLncRNA.py
#     # could add mean number of transcripts per gene.
#     outf.write("\t".join( map( str, [ "multi_exon",
#                                       PipelineLncRNA.CounterTranscripts( me_tmp ).count(),
#                                       PipelineLncRNA.CounterGenes( me_tmp ).count(),
#                                       PipelineLncRNA.CounterExonsPerTranscript( me_tmp ).count(),
#                                       PipelineLncRNA.CounterExonsPerGene( me_tmp ).count(),
#                                       PipelineLncRNA.CounterSingleExonTranscripts( me_tmp ).count(),
#                                       PipelineLncRNA.CounterMultiExonTranscripts( me_tmp ).count(),
#                                       PipelineLncRNA.CounterSingleExonGenes( me_tmp ).count(),
#                                       PipelineLncRNA.CounterMultiExonGenes( me_tmp ).count() ] ) ) )
#     outf.write( "\n" )
#     outf.write("\t".join( map( str, [ "single_exon",
#                                       PipelineLncRNA.CounterTranscripts( se_tmp ).count(),
#                                       PipelineLncRNA.CounterGenes( se_tmp ).count(),
#                                       PipelineLncRNA.CounterExonsPerTranscript( se_tmp ).count(),
#                                       PipelineLncRNA.CounterExonsPerGene( se_tmp ).count(),
#                                       PipelineLncRNA.CounterSingleExonTranscripts( se_tmp ).count(),
#                                       PipelineLncRNA.CounterMultiExonTranscripts( se_tmp ).count(),
#                                       PipelineLncRNA.CounterSingleExonGenes( se_tmp ).count(),
#                                       PipelineLncRNA.CounterMultiExonGenes( se_tmp ).count() ] ) ) )
#     outf.write( "\n" )
#     outf.write("\t".join( map( str, [ "single_exon_merged",
#                                       PipelineLncRNA.CounterTranscripts( sm_tmp ).count(),
#                                       PipelineLncRNA.CounterGenes( sm_tmp ).count(),
#                                       PipelineLncRNA.CounterExonsPerTranscript( sm_tmp ).count(),
#                                       PipelineLncRNA.CounterExonsPerGene( sm_tmp ).count(),
#                                       PipelineLncRNA.CounterSingleExonTranscripts( sm_tmp ).count(),
#                                       PipelineLncRNA.CounterMultiExonTranscripts( sm_tmp ).count(),
#                                       PipelineLncRNA.CounterSingleExonGenes( sm_tmp ).count(),
#                                       PipelineLncRNA.CounterMultiExonGenes( sm_tmp ).count() ] ) ) )
#     outf.write( "\n" )
#     outf.write("\t".join( map( str, [ "multi_exon_merged",
#                                       PipelineLncRNA.CounterTranscripts( mm_tmp ).count(),
#                                       PipelineLncRNA.CounterGenes( mm_tmp ).count(),
#                                       PipelineLncRNA.CounterExonsPerTranscript( mm_tmp ).count(),
#                                       PipelineLncRNA.CounterExonsPerGene( mm_tmp ).count(),
#                                       PipelineLncRNA.CounterSingleExonTranscripts( mm_tmp ).count(),
#                                       PipelineLncRNA.CounterMultiExonTranscripts( mm_tmp ).count(),
#                                       PipelineLncRNA.CounterSingleExonGenes( mm_tmp ).count(),
#                                       PipelineLncRNA.CounterMultiExonGenes( mm_tmp ).count() ] ) ) )               
#     outf.write( "\n" )
#     outf.write("\t".join( map( str, [ "protein_coding",
#                                       PipelineLncRNA.CounterTranscripts( pc_tmp ).count(),
#                                       PipelineLncRNA.CounterGenes( pc_tmp ).count(),
#                                       PipelineLncRNA.CounterExonsPerTranscript( pc_tmp ).count(),
#                                       PipelineLncRNA.CounterExonsPerGene( pc_tmp ).count(),
#                                       PipelineLncRNA.CounterSingleExonTranscripts( pc_tmp ).count(),
#                                       PipelineLncRNA.CounterMultiExonTranscripts( pc_tmp ).count(),
#                                       PipelineLncRNA.CounterSingleExonGenes( pc_tmp ).count(),
#                                       PipelineLncRNA.CounterMultiExonGenes( pc_tmp ).count() ] ) ) )               
#     os.unlink( se_tmp )
#     os.unlink( me_tmp )
#     os.unlink( sm_tmp )
#     os.unlink( mm_tmp )
#     os.unlink( pc_tmp )
#     outf.close()


# @transform( calculateLncRNASummaryStats, 
#             regex( r"(.+)/(.+).tsv" ),
#             r"./\2.load" )
# def loadLncRNASummaryStats( infile, outfile ):
#     P.load( infile, outfile )

@follows( mkdir( "./characterize_lncrna_refcoding" ) )
@transform( combineRefcodingAndLncRNA, 
            regex( "(.+)/(.+).gtf.gz" ), 
            add_inputs( os.path.join( PARAMS[ "genome_dir" ], 
                                      PARAMS[ "genome" ] + ".fasta" ) ),
            r"./characterize_lncrna_refcoding/\2_geneLength_stats.tsv.gz" )
def calculateLncRNAGeneLengthStats( infiles, outfile ):
    """
    This function runs gtf2table --counter=length on the lncRNAs
    _gtf2table.getSegments() returns all intervals in a gtf via GTF.asRanges()
    (NB. getExons() actually returns exon CDS UTR, not just exons)
    GTF.asRanges() returns list of tuples containing interval ( start, end )
    asRanges() uses flat_gene_iterator which returns all intervals assoc. with
    a particular gene_id
    Stats.Summary() returns nval, min, max, mean, median, stdev, sum, q1, q3
    on returned list. 
    Therefore, to get stats for consensus gene_model, it's necessary to generate
    one using gtf2gtf.py --method=merge-exons
    """
    infile, genome = infiles
    tmpf = P.getTempFilename('.')
    statement = ( "zcat %(infile)s |"
                  # " python %(scriptsdir)s/gtf2gtf.py"
                  # "  --method=merge-exons"
                  # "  --log=%(outfile)s.log |"
                  " python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=sort --sort-order=gene+transcript"
                  "  --log=%(outfile)s.log |"
                  " python %(scriptsdir)s/gtf2table.py"
                  "  --counter=length"
                  "  --reporter=genes"
                  "  --genome-file=%(genome)s"
                  "  --log=%(outfile)s.log |"
                  " gzip "
                  " > %(outfile)s;"
                  " checkpoint;" )
    to_cluster = False
    P.run()


@follows( mkdir( "./characterize_lncrna_refcoding" ) )
@transform( combineRefcodingAndLncRNA, 
            regex( "(.+)/(.+).gtf.gz" ), 
            add_inputs( os.path.join( PARAMS[ "genome_dir" ], 
                                      PARAMS[ "genome" ] + ".fasta" ) ),
            r"./characterize_lncrna_refcoding/\2_transcriptLength_stats.tsv.gz" )
def calculateLncRNATranscriptLengthStats( infiles, outfile ):
    """
    This function runs gtf2table --counter=length on the lncRNAs
    """
    infile, genome = infiles
    tmpf = P.getTempFilename('.')
    statement = ( "zcat %(infile)s |"
                  # " python %(scriptsdir)s/gtf2gtf.py"
                  # "  --method=merge-exons"
                  # "  --log=%(outfile)s.log |"
                  " python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=sort --sort-order=gene+transcript"
                  "  --log=%(outfile)s.log |"
                  " python %(scriptsdir)s/gtf2table.py"
                  "  --counter=length"
                  "  --reporter=transcripts"
                  "  --genome-file=%(genome)s"
                  "  --log=%(outfile)s.log |"
                  " gzip "
                  " > %(outfile)s;"
                  " checkpoint;" )
    to_cluster = False
    P.run()


# @transform( calculateLncRNAGeneLengthStats,
#             regex( "(.+)/(.+).tsv.gz" ),
#             r"./\2.load" )
# def loadLncRNAGeneLengthStats( infile, outfile ):
#     P.load( infile, outfile, options="--add-index="  ) # need to index by gene_id


# # gtf2table --counter composition-na composition-cpg length splice
# @transform( combineRefcodingAndLncRNA, 
#             regex( "(.+)/(.+).gtf.gz" ), 
#             add_inputs( os.path.join( PARAMS[ "genome_dir" ], 
#                                       PARAMS[ "genome" ] + ".fasta" ) ),
#             r"./characterize_lncrna_refcoding/\2_composition_stats.tsv.gz" )
# def calculateLncRNAGeneCompositionStats( infiles, outfile ):
#     infile, genome = infiles
#     statement = ( "zcat %(infile)s |"
#                   " python %(scriptsdir)s/gtf2table.py"
#                   "  --counter=composition-na,composition-cpg"
#                   "  --reporter=genes"
#                   "  --genome-file=%(genome)s"
#                   "  --log=%(outfile)s.log |"
#                   " gzip > %(outfile)s" )
#     P.run()


# @transform( calculateLncRNAGeneCompositionStats, 
#             regex( "(.+)/(.+).tsv.gz" ),
#             r"./\2.load" )
# def loadLncRNACompositionStats( infile, outfile ):
#     P.load( infile, outfile, options="--add-index=" ) # need to index by transcript_id


# @follows( loadRefcodingAndLncRNAStatus,
#           loadLncRNASummaryStats, 
#           loadLncRNAGeneLengthStats, 
#           loadLncRNACompositionStats )
# def runLncRNARefcodingSummaryStats(): pass


#################################################################################
#################################################################################
#### METASECTION #### LncRNA Overlap with Previous Genesets ####
#################################################################################
#################################################################################
#################################################################################
# Section:  lncRNA overlap with Hu et al. T cell data & Ensembl lincRNAs
#################################################################################
#################################################################################
## WARNING: There are three duplicate ids in the Tcell data
# lincR-Azi2-5-1672K
# lincR-Crem-3-260K
# lincR-Il6st-5AS-24K
##
@follows( mkdir( "characterize_geneset_overlap_me" ) )
@transform( os.path.join( PARAMS["location_external_datafiles"], 
                          "Hu_et_al",
                          "GSE48138_LincRNA-cluster-coordinates.gtf.gz" ),
            regex( "(.+)/GSE48138_LincRNA-cluster-coordinates.gtf.gz" ),
            add_inputs( "/ifs/mirror/ucsc/mm9/liftOver/mm9ToMm10.over.chain.gz" ),
            r"characterize_geneset_overlap_me/tcell_cluster_cordinates.bed.gz" )
def liftOverTcellLncRNAs( infiles, outfile ):
    """
    Convert the T cell cluster co-ordinates produced by Hu et al. into bed format
    then liftover to mm10.
    """
    gtf_file, chain_file = infiles
    tmpf_bed = P.getTempFilename( "/ifs/scratch" )
    tmpf_lifted = P.getTempFilename( "/ifs/scratch" )
    unlifted = P.snip( outfile, ".gz" ) + "_unlifted" 

    to_cluster = False
    # the infile has a header line...
    statement = ( "zcat %(gtf_file)s |"
                  " sed 1d |"
                  " python %(scriptsdir)s/gff2bed.py"
                  "  --is-gtf"
                  "  --log=%(outfile)s.log"
                  " > %(tmpf_bed)s;"
                  " liftOver"
                  "  %(tmpf_bed)s"
                  "  %(chain_file)s"
                  "  %(tmpf_lifted)s"
                  "  %(unlifted)s;"
                  " cat %(tmpf_lifted)s |"
                  " sort"
                  "  -V"
                  "  -k1,1"
                  "  -k2,2 |"
                  " gzip > %(outfile)s;"
                  " gzip %(unlifted)s" )
    P.run()
    
    os.unlink( tmpf_bed )
    os.unlink( tmpf_lifted )


@follows( mkdir( "characterize_geneset_overlap_me" ) )
@transform( os.path.join( PARAMS["location_transcriptfiles"], 
                          "reference.gtf.gz" ),
            regex( "(.+)/reference.gtf.gz" ),
            r"characterize_geneset_overlap_me/ensembl_lincRNAs.bed.gz" )
def getEnsemblLincRNAs( infile, outfile ):
    """
    Extract ensembl lincRNAs, generate a single interval spanning the entire
    gene model, convert to bed format.
    """
    
    statement = ( "zcat %(infile)s |"
                  " awk '$2 == \"lincRNA\"' |"
                  " python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=sort --sort-order=gene"
                  "  --log=%(outfile)s.log |"
                  " python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=merge-transcripts"
                  "  --log=%(outfile)s.log |"
                  " python %(scriptsdir)s/gff2bed.py"
                  "  --is-gtf"
                  "  --log=%(outfile)s.log |"
                  " gzip > %(outfile)s" )
    P.run()


@follows( mkdir( "characterize_geneset_overlap_me" ) )
@transform( classifyMergedLncRNAs, 
            regex( "(.+)/(.+)_final.gtf.gz" ),
            r"./characterize_geneset_overlap_me/\2_intergenic.bed.gz")
def getIntergenicLncRNAs( infile, outfile ):
    """
    Extract all lncRNAs classified as intergenic, generate a single interval 
    spanning the entire gene model, convert to bed format.
    """
    statement = ( "zcat %(infile)s |"
                  " awk '$2 == \"intergenic\"' |"
                  " python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=merge-transcripts"
                  "  --log=%(outfile)s.log |"
                  " python %(scriptsdir)s/gff2bed.py"
                  "  --is-gtf"
                  "  --log=%(outfile)s.log |"
                  " grep LNCGm |" # Filter out se lncRNA only
                  " gzip > %(outfile)s" )
    P.run()
 

@split( getEnsemblLincRNAs,
        regex( "(.+)/(.+)_lincRNAs.bed.gz" ), 
        add_inputs( liftOverTcellLncRNAs ), 
        r"\1/\2_vs_Tcell*" )
def intersectEnsemblVsTcell( infiles, outfiles ):
    """
    Find the number of Ensembl intervals that intersect and do not intersect
    T cell lncRNA loci
    """
    ensembl_bed, tcell_bed = infiles

    out_overlap = os.path.join( os.path.dirname( ensembl_bed ), 
                                "ensembl_vs_Tcell.overlap.bed" )
    P10.findOverlap( ensembl_bed, tcell_bed, out_overlap )

    out_overlap_2 = os.path.join( os.path.dirname( tcell_bed ), 
                                  "Tcell_vs_ensembl.overlap.bed" )
    P10.findOverlap( tcell_bed, ensembl_bed, out_overlap_2 )

    
@split( getIntergenicLncRNAs,
        regex( "(.+)/(.+)_intergenic.bed.gz" ), 
        add_inputs( liftOverTcellLncRNAs ), 
        r"\1/\2_vs_Tcell*" )
def intersectLncRNAVsTcell( infiles, outfiles ):
    """
    Find the number of lncRNAs that intersect T cell lncRNA loci
    """
    lncRNA_bed, tcell_bed = infiles

    out_overlap = os.path.join( os.path.dirname( lncRNA_bed ), 
                                "lncRNA_vs_Tcell.overlap.bed" )
    P10.findOverlap( lncRNA_bed, tcell_bed, out_overlap )

    out_overlap_2 = os.path.join( os.path.dirname( tcell_bed ), 
                                  "Tcell_vs_lncRNA.overlap.bed" )
    P10.findOverlap( tcell_bed, lncRNA_bed, out_overlap_2 )


@split( getIntergenicLncRNAs,
        regex( "(.+)/(.+)_intergenic.bed.gz" ), 
        add_inputs( getEnsemblLincRNAs ), 
        r"\1/\2_vs_Ensembl*" )
def intersectLncRNAVsEnsembl( infiles, outfiles ):
    """
    Find the number of lncRNAs that intersect T cell lncRNA loci
    """
    lncRNA_bed, ensembl_bed = infiles

    out_overlap = os.path.join( os.path.dirname( lncRNA_bed ), 
                                "lncRNA_vs_Ensembl.overlap.bed" )
    P10.findOverlap( lncRNA_bed, ensembl_bed, out_overlap )

    out_overlap_2 = os.path.join( os.path.dirname( ensembl_bed ), 
                                  "Ensembl_vs_lncRNA.overlap.bed" )
    P10.findOverlap( ensembl_bed, lncRNA_bed, out_overlap_2 )


@collate( [ intersectLncRNAVsEnsembl, intersectLncRNAVsTcell ],
          regex( "(.+)/lncRNA_vs_(.+)\.overlap.bed" ), # escape necessary!
          r"\1/lncRNA_vs_both.overlap.bed" )
def getLncRNAOverlappingBoth( infiles, outfile ):
    """
    Find lncRNA ids that overlap both Ensembl and Tcell loci
    """
    overlap_ensembl, overlap_tcell = infiles
    ensembl = [ line.split()[3] for line in IOTools.openFile( overlap_ensembl ).readlines() ]
    tcell =  [ line.split()[3] for line in IOTools.openFile( overlap_tcell ).readlines() ]
    
    assert len(ensembl) == len(set(ensembl)), "Duplicates in %s" % overlap_ensembl
    assert len(tcell) == len(set(tcell)), "Duplicates in %s" % overlap_tcell
    
    intersect = list(set(ensembl) and set(tcell))

    outf = IOTools.openFile( outfile, "w" )
    for line in IOTools.openFile( overlap_ensembl ).readlines():
        if line.split()[3] in intersect:
            outf.write( line )
    outf.close()


@follows( intersectEnsemblVsTcell, intersectLncRNAVsEnsembl, intersectLncRNAVsTcell )
@merge( [ "characterize_geneset_overlap_me/Tcell_vs_ensembl.no_overlap.bed",
          "characterize_geneset_overlap_me/Tcell_vs_lncRNA.no_overlap.bed" ], 
        "characterize_geneset_overlap_me/Tcell_vs_both.no_overlap.bed" )
def getTcellUnique( infiles, outfile ):
    P10.findIntersect( infiles[0], infiles[1], outfile )


@follows( intersectEnsemblVsTcell, intersectLncRNAVsEnsembl, intersectLncRNAVsTcell )
@merge( [ "characterize_geneset_overlap_me/Ensembl_vs_lncRNA.no_overlap.bed",
          "characterize_geneset_overlap_me/ensembl_vs_Tcell.no_overlap.bed" ], 
        "characterize_geneset_overlap_me/Ensembl_vs_both.no_overlap.bed" )
def getEnsemblUnique( infiles, outfile ):
    P10.findIntersect( infiles[0], infiles[1], outfile )


@follows( intersectEnsemblVsTcell, intersectLncRNAVsEnsembl, intersectLncRNAVsTcell )
@merge( [ "characterize_geneset_overlap_me/lncRNA_vs_Ensembl.no_overlap.bed",
          "characterize_geneset_overlap_me/lncRNA_vs_Tcell.no_overlap.bed" ], 
        "characterize_geneset_overlap_me/lncRNA_vs_both.no_overlap.bed" )
def getLncRNAUnique( infiles, outfile ):
    P10.findIntersect( infiles[0], infiles[1], outfile )

@follows( getTcellUnique, getEnsemblUnique, getLncRNAUnique, getLncRNAOverlappingBoth )
@merge( [ "characterize_geneset_overlap_me/ensembl_lincRNAs.bed.gz", 
          "characterize_geneset_overlap_me/tcell_cluster_cordinates.bed.gz",
          "characterize_geneset_overlap_me/lncRNA_intergenic.bed.gz",
          "characterize_geneset_overlap_me/lncRNA_vs_Ensembl.overlap.bed",
          "characterize_geneset_overlap_me/lncRNA_vs_Tcell.overlap.bed",
          "characterize_geneset_overlap_me/ensembl_vs_Tcell.overlap.bed",
          "characterize_geneset_overlap_me/lncRNA_vs_both.overlap.bed" ],
        "characterize_geneset_overlap_me/summary.tsv" )
def getOverlapSummary( infiles, outfile ):
    """
    Return the intersect between three genesets as required by the R
    package VennDiagram. Intersects return the number of p10 lncRNAs overlapping in
    each case, i.e. if one p10 lncRNA overlaps two tcell lncRNAs... the overlap
    is recorded as one. Ensembl then has priority over t cell.
    """
    # area1 the size of the first set
    lncRNA =  len(IOTools.openFile(infiles[2]).readlines())
    # area2 the size of the second set
    ensembl = len(IOTools.openFile(infiles[0]).readlines())
    # area3 the size of the third set
    tcell =  len(IOTools.openFile(infiles[1]).readlines())

    # n12 size of the intersect between first and second set
    lncRNA_ensembl =  len(IOTools.openFile(infiles[3]).readlines())
    # n23 size of the intersect between second and third set
    ensembl_tcell =  len(IOTools.openFile(infiles[5]).readlines())
    # n13 size of the intersect between first and third set
    lncRNA_tcell =  len(IOTools.openFile(infiles[4]).readlines())

    # n123 size of intersect between all three sets
    center =  len(IOTools.openFile(infiles[6]).readlines())

    outf = IOTools.openFile( outfile, "w" )

    outf.write( "ensembl_total" + "\t" + str(ensembl) + "\n" +
                "tcell_total" + "\t" + str(tcell) + "\n" +
                "lncRNA_total" + "\t" + str(lncRNA) + "\n" +
                "ensembl_tcell" + "\t" + str(ensembl_tcell) + "\n" +
                "lncRNA_tcell" + "\t" + str(lncRNA_tcell) + "\n" +
                "lcnRNA_ensembl" + "\t" + str(lncRNA_ensembl) + "\n" +
                "All" + "\t" + str(center) + "\n" )
    outf.close()


#################################################################################
# Section:  lncRNA overlap with Hu et al. T cell data & Ensembl 78 lincRNAs
#################################################################################
@follows( mkdir( "characterize_geneset_overlap_latest" ) )
@transform( classifyMergedLncRNAs, 
            regex( "(.+)/(.+)_final.gtf.gz" ),
            r"./characterize_geneset_overlap_latest/\2_intergenic.bed.gz")
def getIntergenicLncRNAs_latest( infile, outfile ):
    """
    Extract all lncRNAs classified as intergenic, generate a single interval 
    spanning the entire gene model, convert to bed format.
    """
    statement = ( "zcat %(infile)s |"
                  " awk '$2 == \"intergenic\"' |"
                  " python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=merge-transcripts"
                  "  --log=%(outfile)s.log |"
                  " python %(scriptsdir)s/gff2bed.py"
                  "  --is-gtf"
                  "  --log=%(outfile)s.log |"
                  " gzip > %(outfile)s" )
    P.run()


@follows(mkdir("characterize_geneset_overlap_latest"))
@split( PARAMS["annotations_lincrna_latest"],
        regex( "(.+)/(.+).bed.gz" ), 
        add_inputs( liftOverTcellLncRNAs ), 
        r"characterize_geneset_overlap_latest/ensembl_vs_Tcell*" )
def intersectEnsemblVsTcell_latest( infiles, outfiles ):
    """
    Find the number of Ensembl intervals that intersect and do not intersect
    T cell lncRNA loci
    """
    ensembl_bed, tcell_bed = infiles

    out_overlap = os.path.join( "characterize_geneset_overlap_latest", 
                                "ensembl_vs_Tcell.overlap.bed" )
    P10.findOverlap( ensembl_bed, tcell_bed, out_overlap )

    out_overlap_2 = os.path.join( "characterize_geneset_overlap_latest", 
                                  "Tcell_vs_ensembl.overlap.bed" )
    P10.findOverlap( tcell_bed, ensembl_bed, out_overlap_2 )


@follows(mkdir("characterize_geneset_overlap_latest"))
@split( getIntergenicLncRNAs_latest,
        regex( "(.+)/(.+)_intergenic.bed.gz" ), 
        add_inputs( liftOverTcellLncRNAs ), 
        r"characterize_geneset_overlap_latest/\2_vs_Tcell*" )
def intersectLncRNAVsTcell_latest( infiles, outfiles ):
    """
    Find the number of lncRNAs that intersect T cell lncRNA loci
    """
    lncRNA_bed, tcell_bed = infiles

    out_overlap = os.path.join( "characterize_geneset_overlap_latest", 
                                "lncRNA_vs_Tcell.overlap.bed" )
    P10.findOverlap( lncRNA_bed, tcell_bed, out_overlap )

    out_overlap_2 = os.path.join( "characterize_geneset_overlap_latest", 
                                  "Tcell_vs_lncRNA.overlap.bed" )
    P10.findOverlap( tcell_bed, lncRNA_bed, out_overlap_2 )


@follows(mkdir("characterize_geneset_overlap_latest"))
@split( getIntergenicLncRNAs_latest,
        regex( "(.+)/(.+)_intergenic.bed.gz" ), 
        add_inputs( PARAMS["annotations_lincrna_latest"] ), 
        r"characterize_geneset_overlap_latest/\2_vs_Ensembl*" )
def intersectLncRNAVsEnsembl_latest( infiles, outfiles ):
    """
    Find the number of lncRNAs that intersect T cell lncRNA loci
    """
    lncRNA_bed, ensembl_bed = infiles

    out_overlap = os.path.join( "characterize_geneset_overlap_latest", 
                                "lncRNA_vs_Ensembl.overlap.bed" )
    P10.findOverlap( lncRNA_bed, ensembl_bed, out_overlap )

    out_overlap_2 = os.path.join( "characterize_geneset_overlap_latest", 
                                  "Ensembl_vs_lncRNA.overlap.bed" )
    P10.findOverlap( ensembl_bed, lncRNA_bed, out_overlap_2 )


@collate( [ intersectLncRNAVsEnsembl_latest, intersectLncRNAVsTcell_latest ],
          regex( "(.+)/lncRNA_vs_(.+)\.overlap.bed" ), # escape necessary!
          r"\1/lncRNA_vs_both.overlap.bed" )
def getLncRNAOverlappingBoth_latest( infiles, outfile ):
    """
    Find lncRNA ids that overlap both Ensembl and Tcell loci
    """
    overlap_ensembl, overlap_tcell = infiles
    ensembl = [ line.split()[3] for line in IOTools.openFile( overlap_ensembl ).readlines() ]
    tcell =  [ line.split()[3] for line in IOTools.openFile( overlap_tcell ).readlines() ]
    
    assert len(ensembl) == len(set(ensembl)), "Duplicates in %s" % overlap_ensembl
    assert len(tcell) == len(set(tcell)), "Duplicates in %s" % overlap_tcell
    
    intersect = list(set(ensembl) and set(tcell))

    outf = IOTools.openFile( outfile, "w" )
    for line in IOTools.openFile( overlap_ensembl ).readlines():
        if line.split()[3] in intersect:
            outf.write( line )
    outf.close()


@follows( intersectEnsemblVsTcell_latest, 
          intersectLncRNAVsEnsembl_latest,
          intersectLncRNAVsTcell_latest )
@merge( [ "characterize_geneset_overlap_latest/Tcell_vs_ensembl.no_overlap.bed",
          "characterize_geneset_overlap_latest/Tcell_vs_lncRNA.no_overlap.bed" ], 
        "characterize_geneset_overlap_latest/Tcell_vs_both.no_overlap.bed" )
def getTcellUnique_latest( infiles, outfile ):
    P10.findIntersect( infiles[0], infiles[1], outfile )


@follows( intersectEnsemblVsTcell_latest, 
          intersectLncRNAVsEnsembl_latest, 
          intersectLncRNAVsTcell_latest )
@merge( [ "characterize_geneset_overlap_latest/Ensembl_vs_lncRNA.no_overlap.bed",
          "characterize_geneset_overlap_latest/ensembl_vs_Tcell.no_overlap.bed" ], 
        "characterize_geneset_overlap_latest/Ensembl_vs_both.no_overlap.bed" )
def getEnsemblUnique_latest( infiles, outfile ):
    P10.findIntersect( infiles[0], infiles[1], outfile )


@follows( intersectEnsemblVsTcell_latest, 
          intersectLncRNAVsEnsembl_latest, 
          intersectLncRNAVsTcell_latest )
@merge( [ "characterize_geneset_overlap_latest/lncRNA_vs_Ensembl.no_overlap.bed",
          "characterize_geneset_overlap_latest/lncRNA_vs_Tcell.no_overlap.bed" ], 
        "characterize_geneset_overlap_latest/lncRNA_vs_both.no_overlap.bed" )
def getLncRNAUnique_latest( infiles, outfile ):
    P10.findIntersect( infiles[0], infiles[1], outfile )


@follows( getTcellUnique_latest, 
          getEnsemblUnique_latest, 
          getLncRNAUnique_latest, 
          getLncRNAOverlappingBoth_latest )
@merge( [ PARAMS["annotations_lincrna_latest"], 
          "characterize_geneset_overlap_me/tcell_cluster_cordinates.bed.gz",
          "characterize_geneset_overlap_latest/lncRNA_intergenic.bed.gz",
          "characterize_geneset_overlap_latest/lncRNA_vs_Ensembl.overlap.bed",
          "characterize_geneset_overlap_latest/lncRNA_vs_Tcell.overlap.bed",
          "characterize_geneset_overlap_latest/ensembl_vs_Tcell.overlap.bed",
          "characterize_geneset_overlap_latest/lncRNA_vs_both.overlap.bed" ],
        "characterize_geneset_overlap_latest/summary.tsv" )
def getOverlapSummary_latest( infiles, outfile ):
    """
    Return the intersect between three genesets as required by the R
    package VennDiagram. Intersects return the number of p10 lncRNAs overlapping in
    each case, i.e. if one p10 lncRNA overlaps two tcell lncRNAs... the overlap
    is recorded as one. Ensembl then has priority over t cell.
    """
    # area1 the size of the first set
    lncRNA =  len(IOTools.openFile(infiles[2]).readlines())
    # area2 the size of the second set
    ensembl = len(IOTools.openFile(infiles[0]).readlines())
    # area3 the size of the third set
    tcell =  len(IOTools.openFile(infiles[1]).readlines())

    # n12 size of the intersect between first and second set
    lncRNA_ensembl =  len(IOTools.openFile(infiles[3]).readlines())
    # n23 size of the intersect between second and third set
    ensembl_tcell =  len(IOTools.openFile(infiles[5]).readlines())
    # n13 size of the intersect between first and third set
    lncRNA_tcell =  len(IOTools.openFile(infiles[4]).readlines())

    # n123 size of intersect between all three sets
    center =  len(IOTools.openFile(infiles[6]).readlines())

    outf = IOTools.openFile( outfile, "w" )

    outf.write( "ensembl_total" + "\t" + str(ensembl) + "\n" +
                "tcell_total" + "\t" + str(tcell) + "\n" +
                "lncRNA_total" + "\t" + str(lncRNA) + "\n" +
                "ensembl_tcell" + "\t" + str(ensembl_tcell) + "\n" +
                "lncRNA_tcell" + "\t" + str(lncRNA_tcell) + "\n" +
                "lcnRNA_ensembl" + "\t" + str(lncRNA_ensembl) + "\n" +
                "All" + "\t" + str(center) + "\n" )
    outf.close()


#################################################################################
# Section:  lncRNA overlap with Hu et al. T cell data & Mouse Encode Contigs
#################################################################################


@follows(mkdir("characterize_geneset_overlap_encode"))
@transform(os.path.join(PARAMS["location_external_datafiles"],
                        "mouse_encode",
                        "Supp_Table_S4_cshl_cell_longPolyA_projected_contigs.bed.gz"),
           regex("(.+)/Supp_Table_S4_cshl_cell_longPolyA_projected_contigs.bed.gz"),
           add_inputs( "/ifs/mirror/ucsc/mm9/liftOver/mm9ToMm10.over.chain.gz" ),
           "characterize_geneset_overlap_encode/mencode_cluster_coordinates.bed.gz")
def liftOverMouseEncodeIntervals( infiles, outfile ):
    """
    Lift over the half million intervals identified by the mouse encode project.
    """
    bed_file, chain_file = infiles
    tmpf_bed = P.getTempFilename( "/ifs/scratch" )
    tmpf_lifted = P.getTempFilename( "/ifs/scratch" )
    unlifted = P.snip( outfile, ".gz" ) + "_unlifted" 

    to_cluster = False
    # the infile has a header line...
    statement = ( "zcat %(bed_file)s"
                  " > %(tmpf_bed)s;"
                  " liftOver"
                  "  %(tmpf_bed)s"
                  "  %(chain_file)s"
                  "  %(tmpf_lifted)s"
                  "  %(unlifted)s;"
                  " cat %(tmpf_lifted)s |"
                  " sort"
                  "  -V"
                  "  -k1,1"
                  "  -k2,2 |"
                  " gzip > %(outfile)s;"
                  " gzip %(unlifted)s" )
    P.run()
    
    os.unlink( tmpf_bed )
    os.unlink( tmpf_lifted )


@transform( liftOverMouseEncodeIntervals, 
            suffix( ".bed.gz" ),
            add_inputs( os.path.join( PARAMS["location_transcriptfiles"],
                                      "reference.gtf.gz" ), 
                        os.path.join( PARAMS[ "annotations_dir" ], 
                                      PARAMS_ANNOTATIONS[ "interface_pseudogenes_gtf" ] ),
                        os.path.join( PARAMS[ "annotations_dir" ],
                                      PARAMS_ANNOTATIONS[ "interface_numts_gtf" ] ) ),
            "_filtered_1.bed.gz" )
def filterMouseEncodeAgainstEnsemblAnnotations( infiles, outfile ):
    """
    Following filterAgainstBiotypes, this task filters mouse encode intervals against
    annotations that intersect ensembl annotations. 
    It also filters things 5kb upstr or dstr of them. 
    """

    encode_bed, reference_gtf, pseudogene_gtf, numt_gtf  = infiles
    rejected = P.snip( outfile, "_filtered_1.bed.gz" ) + "_rejected_1.gtf.gz"
    tmpf_1_name  = ( "./characterize_geneset_overlap_encode/biotypes.gtf" )
    tmpf_2_name  = ( "./characterize_geneset_overlap_encode/biotypes_1.gtf" )
    tmpf_1 = IOTools.openFile( tmpf_1_name, "w" )

    # pull out gtf entries matching the specified biotypes
    biotypes = PARAMS[ "filter_intersect_vs" ].split( "," )
    for gtf in GTF.iterator( IOTools.openFile( reference_gtf ) ):
        if gtf.source in biotypes:
            tmpf_1.write( str( gtf ) + "\n" )
        else: continue
    tmpf_1.close()
    E.info( "Biotypes to be removed have been"
            " written to %s" % tmpf_1_name )
    # add known pseudogenes
    statement = ("zcat %(pseudogene_gtf)s >> %(tmpf_1_name)s")
    to_cluster = False
    P.run()

    # merge gtf intervals into a single interval and slop by 5kb.
    genome_file = os.path.join(PARAMS["annotations_dir"], 
                               PARAMS_ANNOTATIONS["interface_contigs"])
    to_cluster = False
    statement = ("cat %(tmpf_1_name)s |"
                 " python %(scriptsdir)s/gtf2gtf.py"
                 " --method=sort --sort-order=gene"
                 " --log=%(outfile)s.log |"
                 " python %(scriptsdir)s/gtf2gtf.py"
                 " --method=merge-transcripts"
                 " --log=%(outfile)s.log |"
                 " python %(scriptsdir)s/gff2bed.py"
                 " --log=%(outfile)s.log"
                 " --is-gtf |"
                 " bedtools slop"
                 "  -i stdin"
                 "  -b 5000"
                 "  -g %(genome_file)s"
                 " > %(tmpf_2_name)s")
    P.run()

    # remove any encode interval that intersects these intervals. 
    tmpf_3_name = P.getTempFilename("/ifs/scratch")
    to_cluster = False
    statement = ("bedtools intersect"
                 " -a %(encode_bed)s"
                 " -b %(tmpf_2_name)s"
                 " -v"
                 " > %(tmpf_3_name)s")
    P.run()

    # for good measure, remove anything from encode that could be numt
    tmpf_4_name = P.getTempFilename("/ifs/scratch")
    to_cluster = False
    statement = ("zcat %(numt_gtf)s |"
                 " python %(scriptsdir)s/gff2bed.py"
                 " --is-gtf"
                 " --log=/dev/null"
                 " > %(tmpf_4_name)s;"
                 " bedtools intersect"
                 "  -a %(tmpf_3_name)s"
                 "  -b %(tmpf_4_name)s"
                 "  -v |"
                 " gzip > %(outfile)s" )
    P.run()

    os.unlink(tmpf_1_name)
    os.unlink(tmpf_2_name)
    os.unlink(tmpf_3_name)
    os.unlink(tmpf_4_name)


@transform(classifyMergedLncRNAs,
           regex("(.+)/(.+)_final.gtf.gz"),
           add_inputs(liftOverMouseEncodeIntervals),
           r"characterize_geneset_overlap_encode/\2_mouse_encode_overlap.tsv")
def getLncRNAvsEncode_Overlap(infiles, outfile):
    """
    Use gtf2table to count the overlap between lncRNAs and mouse encode
    expressed transcripts. 
    N.B. 
    nover1 is the number of exons in lncRNA that overlap annotation, 
    nover2 is the number of exons in encode that overlap lncRNA
    nover is the number of bases overlapped. 
    pover1 is the percentage of lncRNA that is overlapped by annotation
    pover2 is the percentage of annotations overlapped by lncRNA. 
    """
    lnc_gtf, encode_bed = infiles
    

    # calculate overlap between intergenic lncRNAs and encode intervals. 
    statement = ("zcat %(lnc_gtf)s |"
                 " python %(scriptsdir)s/gtf2table.py"
                 "  --reporter=genes"
                 "  --counter=overlap"
                 "  --gff-file=%(encode_bed)s"
                 "  --filename-format=bed"
                 "  --log=%(outfile)s.log"
                 "  > %(outfile)s")
    P.run()
    

@transform(getLncRNAvsEncode_Overlap,
           regex(".+/(.+).tsv"),
           r"\1.load")
def loadLncRNAvsEncode_Overlap(infile, outfile):
    P.load(infile, outfile, options="--add-index=gene_id")

#################################################################################
#################################################################################
#################################################################################
# Section: Compare lncRNA tss with FANTOM5 Data
#################################################################################
#################################################################################
# connect to FANTOM BioMart and retrieve i) sample details, ii) robust peaks
# summarise distance between robust peaks and ensembl gene_ids
# summarise distance between robust peaks and lncRNAs
#################################################################################
## subsection: liftOver FANTOM data from mm9 to mm10
#################################################################################
@follows( mkdir( "characterize_fantom5_overlap" ) )
@transform( os.path.join( PARAMS["location_external_datafiles"], 
                          "FANTOM5", 
                          "mm9.cage_peak_coord_robust.bed.gz" ),
            regex( "(.+)/mm9.cage_peak_coord_robust.bed.gz" ),
            add_inputs( "/ifs/mirror/ucsc/mm9/liftOver/mm9ToMm10.over.chain.gz" ),
            r"./characterize_fantom5_overlap/robust_cage_peaks_mm10.bed.gz" )
def liftOverRobustPeaks( infiles, outfile ):
    """
    Lift over the robust CAGE peak intervals from mm9 to mm10
    """
    peakfile, chainfile = infiles
    tmpfile = P.getTempFilename(".")
    unlifted = P.snip( re.sub( "mm10", "unlifted", outfile ), ".gz" )

    statement = ( "liftOver"
                  "  %(peakfile)s"
                  "  %(chainfile)s"
                  "  %(tmpfile)s"
                  "  %(unlifted)s;"
                  " cat %(tmpfile)s |"
                  " sort"
                  "  -V"
                  "  -k1,1"
                  "  -k2,2 |"
                  " gzip > %(outfile)s;"
                  " gzip %(unlifted)s" )
    P.run()

    os.unlink( tmpfile )
    

@follows( mkdir( "characterize_fantom5_overlap" ) )
@transform( os.path.join( PARAMS["location_external_datafiles"], 
                          "FANTOM5", 
                          "mm9.cage_peak_coord_permissive.bed.gz" ),
            regex( "(.+)/mm9.cage_peak_coord_permissive.bed.gz" ),
            add_inputs( "/ifs/mirror/ucsc/mm9/liftOver/mm9ToMm10.over.chain.gz" ),
            r"./characterize_fantom5_overlap/permissive_cage_peaks_mm10.bed.gz" )
def liftOverPermissivePeaks( infiles, outfile ):
    """
    Lift over the permissive CAGE peak intervals from mm9 to mm10
    """
    peakfile, chainfile = infiles
    tmpfile = P.getTempFilename(".")
    unlifted = P.snip( re.sub( "mm10", "unlifted", outfile ), ".gz" )

    statement = ( "liftOver"
                  "  %(peakfile)s"
                  "  %(chainfile)s"
                  "  %(tmpfile)s"
                  "  %(unlifted)s;"
                  " cat %(tmpfile)s |"
                  " sort"
                  "  -V"
                  "  -k1,1"
                  "  -k2,2 |"
                  " gzip > %(outfile)s;"
                  " gzip %(unlifted)s" )
    P.run()

    os.unlink( tmpfile )


@follows( mkdir( "characterize_fantom5_overlap" ) )
@transform( os.path.join( PARAMS["location_external_datafiles"], 
                          "FANTOM5", 
                          "mm9.cage_peak_*_permissive.bed.gz" ),
            regex( "(.+)/mm9.cage_peak_((?!coord).+)_permissive.bed.gz" ),
            add_inputs( "/ifs/mirror/ucsc/mm9/liftOver/mm9ToMm10.over.chain.gz" ),
            r"./characterize_fantom5_overlap/\2_cage_peaks_mm10.bed.gz" )
def liftOverEnsemblCagePeaks( infiles, outfile ):
    """
    Lift over the permissive CAGE peak intervals annotated by Wilfried as either:
    i) assocated with an ENSEMBL id
    ii) not associated with an ENSEMBL id
    """
    peakfile, chainfile = infiles
    tmpfile = P.getTempFilename(".")
    unlifted = P.snip( re.sub( "mm10", "unlifted", outfile ), ".gz" )

    statement = ( "liftOver"
                  "  %(peakfile)s"
                  "  %(chainfile)s"
                  "  %(tmpfile)s"
                  "  %(unlifted)s;"
                  " cat %(tmpfile)s |"
                  " sort"
                  "  -V"
                  "  -k1,1"
                  "  -k2,2 |"
                  " gzip > %(outfile)s;"
                  " gzip %(unlifted)s" )
    P.run()

    os.unlink( tmpfile )


@follows( liftOverRobustPeaks, 
          liftOverPermissivePeaks, 
          liftOverEnsemblCagePeaks )
def liftOverCAGEPeaks():
    pass


@collate( liftOverEnsemblCagePeaks, 
          regex( "(.+)/(.+)_cage_peaks_mm10.bed.gz" ), 
          add_inputs( os.path.join( PARAMS["location_transcriptfiles"], 
                                    "reference.gtf.gz" ) ),
          r"\1/ensembl_nc_cage_peaks_mm10.bed.gz" )
def extractNonCodingCagePeaks( infiles, outfile ):
    """
    Select known non-coding and unannotated CAGE intervals based on the 
    CAGE annotations supplied by Wilfried
    """
    cage_files = [ x[0] for x in infiles ]
    ref_gtf = infiles[0][1]
    control_biotypes = PARAMS["filter_rescue_biotypes"].split(",")

    non_coding_ensembl = set()
    for gtf in GTF.iterator( IOTools.openFile( ref_gtf ) ):
        if gtf.source in control_biotypes:
            non_coding_ensembl.add( gtf.gene_id )

    outf = IOTools.openFile( outfile, "w" )
    for infile in cage_files:
        if re.search( "unannotated", infile ):
            for line in IOTools.openFile( infile ):
                outf.write( line )
        elif re.search( "ensembl", infile ):
            n_nc = 0
            for line in IOTools.openFile( infile ):
                if line.split()[3] in non_coding_ensembl:
                    n_nc += 1
                    outf.write( line )
                else:
                    continue
        else:
            raise Exception( "Unrecognised infile %s" % infile )

    E.info( "There are %i cage intervals associated with"
            " known non-coding intervals" % n_nc )
    outf.close()


@merge( [ extractNonCodingCagePeaks, liftOverPermissivePeaks ],
        "./characterize_fantom5_overlap/permissivenc_cage_peaks_mm10.bed.gz" )
def collateNonCodingCagePeaks_permissive( infiles, outfile ):
    """
    Intersect the original permissive cage peak intervals, with those
    identified by Wilfried as being either i) unannotated, or ii) non-coding. 
    Write out the interval that appears in the original permissive cage
    peak file.
    """
    
    for infile in infiles:
        if re.match( "permissive", os.path.basename(infile) ):
            cage_peaks = infile
        else:
            nc_peaks = infile

    statement = ( "bedtools intersect"
                  " -a %(cage_peaks)s"
                  " -b %(nc_peaks)s"
                  " -wa" # write interval in a for each overlap
                  " -u"  # write interval in a only once
                  " -s |"  # force strandedness
                  " gzip > %(outfile)s" )
    P.run()

    n_nc = len( IOTools.openFile( nc_peaks ).readlines() )
    n_out = len( IOTools.openFile( outfile ).readlines() )

    assert n_nc == n_out, "The number of intervals in the noncoding cage \
    file doesn't match the number of intervals in the outfile" 


@merge( [ extractNonCodingCagePeaks, liftOverRobustPeaks ],
        "./characterize_fantom5_overlap/robustnc_cage_peaks_mm10.bed.gz" )
def collateNonCodingCagePeaks_robust( infiles, outfile ):
    """
    Intersect the original robust cage peak intervals, with those identified
    by Wilfried as being either i) unannotated, or ii) non-coding. Write out
    the interval that appears in the original permissive cage peak file.

    """
    for infile in infiles:
        if re.match( "robust", os.path.basename(infile) ):
            cage_peaks = infile
        else:
            nc_peaks = infile

    statement = ( "bedtools intersect"
                  " -a %(cage_peaks)s"
                  " -b %(nc_peaks)s"
                  " -wa" # write interval in a for each overlap
                  " -u"  # write interval in a only once
                  " -s |"  # force strandedness
                  " gzip > %(outfile)s" )
    P.run()


@follows( collateNonCodingCagePeaks_permissive,
          collateNonCodingCagePeaks_robust )
def collateNonCodingCagePeaks():
    pass


#################################################################################
## subsection: generate bed/gtf files of lncRNA transcription start sites
#################################################################################                 
@follows( mkdir( "characterize_fantom5_overlap" ) )
@split( classifyMergedLncRNAs, 
        regex( "(.+)/(.+)_final.gtf.gz" ), 
        [ r"characterize_fantom5_overlap/\2_gene_tss.gtf.gz", 
          r"characterize_fantom5_overlap/\2_gene_tss.bed.gz" ] )
def findLncRNAGeneTSS( infile, outfiles ): 
    """
    Returns gff containing the transcription start site for each gene model
    """
    out_gtf, out_bed = outfiles
    out_gtf = P.snip( out_gtf, ".gz" )
    out_bed = P.snip( out_bed, ".gz" )

    statement = ( "zcat %(infile)s |"
                  " python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=merge-transcripts"
                  "  --log=%(out_gtf)s.log |"
                  " python %(scriptsdir)s/gtf2gff.py"
                  "  --method=promotors"
                  "  --promotor-size=1"
                  "  --genome-file=%(genome_dir)s/%(genome)s"
                  "  --log=%(out_gtf)s.log |"
                  " sed s/promotor/tss/ |"
                  " sed s/exon/gene/ |"
                  " python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=sort --sort-order=gene"
                  "  --log=%(out_gtf)s.log |"
                  " tee %(out_gtf)s |"
                  " python %(scriptsdir)s/gff2bed.py"
                  "  --is-gtf"
                  "  --set-name=gene_id"
                  "  --log=%(out_bed)s.log |"
                  " gzip -f > %(out_bed)s.gz;"
                  " gzip -f %(out_gtf)s" )
    P.run()


@follows( mkdir( "characterize_fantom5_overlap" ) )
@split( classifyMergedLncRNAs, 
        regex( "(.+)/(.+)_final.gtf.gz" ), 
        [ r"characterize_fantom5_overlap/\2_gene_tts.gtf.gz", 
          r"characterize_fantom5_overlap/\2_gene_tts.bed.gz" ] )
def findLncRNAGeneTTS( infile, outfiles ): 
    """
    Returns gff containing the transcription termination site for each gene model
    """
    out_gtf, out_bed = outfiles
    out_gtf = P.snip( out_gtf, ".gz" )
    out_bed = P.snip( out_bed, ".gz" )

    # merge-transcripts returns a single interval spanning the whole gene
    statement = ( "zcat %(infile)s |"    
                  " python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=merge-transcripts"
                  "  --log=%(out_gtf)s.log |"
                  " python %(scriptsdir)s/gtf2gff.py"
                  "  --method=tts"
                  "  --promotor-size=1"
                  "  --genome-file=%(genome_dir)s/%(genome)s"
                  "  --log=%(out_gtf)s.log |"
                  " sed s/promotor/tts/ |"
                  " sed s/exon/gene/ |"
                  " python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=sort --sort-order=gene"
                  "  --log=%(out_gtf)s.log |"
                  " tee %(out_gtf)s |"
                  " python %(scriptsdir)s/gff2bed.py"
                  "  --is-gtf"
                  "  --set-name=gene_id"
                  "  --log=%(out_bed)s.log |"
                  " gzip -f > %(out_bed)s.gz;"
                  " gzip -f %(out_gtf)s" )
    P.run()


@follows( mkdir( "characterize_fantom5_overlap" ) )
@split( classifyMergedLncRNAs, 
        regex( "(.+)/(.+)_final.gtf.gz" ), 
        [ r"characterize_fantom5_overlap/\2_transcript_tss.gtf.gz", 
          r"characterize_fantom5_overlap/\2_transcript_tss.bed.gz" ] )
def findLncRNATranscriptTSS( infile, outfiles ): 
    """
    Returns gff file containing the transcription start site for each transcript.
    Taken from pipeline_annotations.py
    """
    out_gtf, out_bed = outfiles
    out_gtf = P.snip( out_gtf, ".gz" )
    out_bed = P.snip( out_bed, ".gz" )

    # (?) join-exons returns a new interval that spans the whole transcript
    statement = ( "zcat %(infile)s |"    
                  " python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=join-exons"
                  "  --log=%(out_gtf)s.log |"
                  " python %(scriptsdir)s/gtf2gff.py"
                  "  --method=promotors"
                  "  --promotor-size=1"
                  "  --genome-file=%(genome_dir)s/%(genome)s"
                  "  --log=%(out_gtf)s.log |"
                  " sed s/promotor/tss/ |"
                  " python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=sort --sort-order=gene"
                  "  --log=%(out_gtf)s.log |"
                  " tee %(out_gtf)s |"
                  " python %(scriptsdir)s/gff2bed.py"
                  "  --is-gtf"
                  "  --set-name=transcript_id"
                  "  --log=%(out_bed)s.log |"
                  " gzip -f > %(out_bed)s.gz;"
                  " gzip -f %(out_gtf)s" )
    P.run()


@transform( findLncRNAGeneTSS,
            suffix( ".bed.gz" ), 
            add_inputs( os.path.join( PARAMS[ "annotations_dir" ], 
                                      PARAMS_ANNOTATIONS[ "interface_contigs" ] ) ),
            "_1kb.bed.gz" )
def slopLncRNAGeneTSS( infiles, outfile ):
    """
    Create 1kb windows centered on lncRNA gene transcription start sites
    """
    in_bed, in_contigs = infiles
    statement = ( "bedtools slop" 
                   " -b 500"
                   " -i %(in_bed)s"
                   " -g %(in_contigs)s |"
                   " gzip > %(outfile)s" )
    P.run()


@transform( findLncRNAGeneTTS,
            suffix( ".bed.gz" ), 
            add_inputs( os.path.join( PARAMS[ "annotations_dir" ], 
                                      PARAMS_ANNOTATIONS[ "interface_contigs" ] ) ),
            "_1kb.bed.gz" )
def slopLncRNAGeneTTS( infiles, outfile ):
    """
    Create 1kb windows centered on lncRNA gene transcription start sites
    """
    in_bed, in_contigs = infiles
    statement = ( "bedtools slop" 
                   " -b 500"
                   " -i %(in_bed)s"
                   " -g %(in_contigs)s |"
                   " gzip > %(outfile)s" )
    P.run()


@transform( findLncRNATranscriptTSS,
            suffix( ".bed.gz" ), 
            add_inputs( os.path.join( PARAMS[ "annotations_dir" ], 
                                      PARAMS_ANNOTATIONS[ "interface_contigs" ] ) ),
            "_1kb.bed.gz" )
def slopLncRNATranscriptTSS( infiles, outfile ):
    """
    Create 1kb windows centered on lncRNA gene transcription start sites
    """
    in_bed, in_contigs = infiles
    statement = ( "bedtools slop" 
                   " -b 500"
                   " -i %(in_bed)s"
                   " -g %(in_contigs)s |"
                   " gzip > %(outfile)s" )
    P.run()


@transform( [ findLncRNAGeneTSS, findLncRNATranscriptTSS ],
             suffix( ".bed.gz" ),
             add_inputs( os.path.join( PARAMS[ "annotations_dir" ], 
                                       PARAMS_ANNOTATIONS[ "interface_contigs" ] ) ),
             "_4kb.bed.gz" )
def slopLncRNATSS_4kb( infiles, outfile ):
    """
    Create 4kb windows surrounding the TSS of all lncRNA genes|transcripts
    """
    in_bed, in_contigs = infiles 
    statement = ( "bedtools slop" 
                   " -b 2000"
                   " -i %(in_bed)s"
                   " -g %(in_contigs)s |"
                   " gzip > %(outfile)s" )
    P.run()
    

@follows( slopLncRNAGeneTSS,
          slopLncRNAGeneTTS,
          slopLncRNATranscriptTSS,          
          slopLncRNATSS_4kb )
def findLncRNATSS():
    pass


#################################################################################
## subsection: generate bed files of protein coding transcription start sites
#################################################################################                 
@transform( os.path.join( PARAMS["location_transcriptfiles"], "refcoding.gtf.gz" ),
            regex( "(.+)/(.+).gtf.gz" ), 
            r"./characterize_fantom5_overlap/\2_gene_tss.bed.gz" )
def findRefcodingGeneTSS( infile, outfile ): 
    """
    Returns gff containing the transcription start site for each gene model
    Only bed file is dependent outfile, bc used as add_inputs later
    """

    out_bed = outfile
    out_gtf = P.snip( out_bed, "bed.gz" ) + "gtf"
    out_bed = P.snip( out_bed, ".gz" )

    statement = ( "zcat %(infile)s |"
                  " python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=merge-transcripts"
                  "  --log=%(out_gtf)s.log |"
                  " python %(scriptsdir)s/gtf2gff.py"
                  "  --method=promotors"
                  "  --promotor-size=1"
                  "  --genome-file=%(genome_dir)s/%(genome)s"
                  "  --log=%(out_gtf)s.log |"
                  " sed s/promotor/tss/ |"
                  " sed s/exon/gene/ |"
                  " python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=sort --sort-order=gene"
                  "  --log=%(out_gtf)s.log |"
                  " tee %(out_gtf)s |"
                  " python %(scriptsdir)s/gff2bed.py"
                  "  --is-gtf"
                  "  --set-name=gene_id"
                  "  --log=%(out_bed)s.log |"
                  " gzip -f > %(out_bed)s.gz;"
                  " gzip -f %(out_gtf)s" )
    P.run()


@transform( os.path.join( PARAMS["location_transcriptfiles"], "refcoding.gtf.gz" ), 
            regex( "(.+)/(.+).gtf.gz" ), 
            r"./characterize_fantom5_overlap/\2_transcript_tss.bed.gz" )
def findRefcodingTranscriptTSS( infile, outfile ): 
    """
    Returns gff file containing the transcription start site for each transcript.
    Taken from pipeline_annotations.py
    """

    out_bed = outfile
    out_gtf = P.snip( out_bed, "bed.gz" ) + "gtf"
    out_bed = P.snip( out_bed, ".gz" )

    # (?) join-exons returns a new interval that spans the whole transcript
    statement = ( "zcat %(infile)s |"    
                  " python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=join-exons"
                  "  --log=%(out_gtf)s.log |"
                  " python %(scriptsdir)s/gtf2gff.py"
                  "  --method=promotors"
                  "  --promotor-size=1"
                  "  --genome-file=%(genome_dir)s/%(genome)s"
                  "  --log=%(out_gtf)s.log |"
                  " sed s/promotor/tss/ |"
                  " python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=sort --sort-order=gene"
                  "  --log=%(out_gtf)s.log |"
                  " tee %(out_gtf)s |"
                  " python %(scriptsdir)s/gff2bed.py"
                  "  --is-gtf"
                  "  --set-name=transcript_id"
                  "  --log=%(out_bed)s.log |"
                  " gzip -f > %(out_bed)s.gz;"
                  " gzip -f %(out_gtf)s" )
    P.run()


@follows( findRefcodingGeneTSS, 
          findRefcodingTranscriptTSS )
def findRefcodingTSS():
    pass


#################################################################################
## subsection: calculate distance from lncRNA TSS to nearest refcoding TSS
#################################################################################
@transform( [ findLncRNAGeneTSS, findLncRNATranscriptTSS ], 
            regex( "(.+)/(.+)_(gene|transcript)_tss.bed.gz" ),
            add_inputs( findRefcodingGeneTSS ), 
            r"\1/\2_\3_refcoding_distance.tsv.gz" )
def calcLncRNATSSRefcodingDistance( infiles, outfile ):
    """
    Distance between lncRNA Gene/Transcript TSS and nearest protein coding gene
    TSS.
    """
    lnc_tss, refcoding_tss = infiles
    tmpf = P.getTempFilename( "/ifs/scratch" )

    statement = ( "bedtools closest"
                  " -a %(lnc_tss)s"
                  " -b %(refcoding_tss)s"
                  " -s" # enforce strandedness
                  " -d" # report distance of b from a as extra column
                  " > %(tmpf)s"
                  " 2> %(outfile)s.log" )
    P.run()

    P10.postProcessClosestBed( tmpf, outfile, cage=False )
    os.unlink( tmpf )


@transform( [ findLncRNAGeneTSS, findLncRNATranscriptTSS ], 
            regex( "(.+)/(.+)_(gene|transcript)_tss.bed.gz" ),
            add_inputs( findRefcodingGeneTSS ), 
            r"\1/\2_\3_refcoding_Nulldistance.tsv.gz" )
def calcLncRNATSSRefcodingNullDistance( infiles, outfile ):
    """
    Null distribution for distance between lncNRA TSS and nearest robust cage 
    peak calculated as distance to nearest cage peak on oposite strand.
    """
    lnc_tss, refcoding_tss = infiles
    tmpf = P.getTempFilename( "/ifs/scratch" )

    statement = ( "bedtools closest"
                  " -a %(lnc_tss)s"
                  " -b %(refcoding_tss)s"
                  " -S" # enforce opposite strandedness
                  " -d" # report distance of b from a as extra column
                  " > %(tmpf)s"
                  " 2> %(outfile)s.log" )
    P.run()

    P10.postProcessClosestBed( tmpf, outfile, cage=False )
    os.unlink( tmpf )


@collate( [ calcLncRNATSSRefcodingDistance,
            calcLncRNATSSRefcodingNullDistance ], 
          regex( "(.+)/(.+)_refcoding_(distance|Nulldistance).tsv.gz" ), 
          r"\1/lncRNA_TSS_refcoding_distance.tsv.gz" )
def stackLncRNARefcodingTSSDistance( infiles, outfile ):
    tables = " ".join( infiles )
    to_cluster = False
    tmpf = P.getTempFilename(".")
    statement = ( "python %(scriptsdir)s/combine_tables.py"
                  "  --take=10,16" # does nothing in this context!
                  "  --cat=CAT"
                  "  --log=%(outfile)s.log"
                  " %(tables)s |"
                  " cut -f1,5,11,14"
                  " > %(tmpf)s" )
    P.run()

    # split the filename into additional columns
    outf = IOTools.openFile( outfile, "w" )
    outf.write( "interval\tmodel\tlnc_id\trefcoding_id\tdistance\n" )
    header = True
    for line in IOTools.openFile( tmpf ):
        if header: 
            header = False
            continue
        line = line.split()
        # sort filename
        filename = os.path.basename( line[0] )[:-7]
        filename = filename.split("_")
        line_out = [ filename[1], filename[3] ] + line[1:]
        outf.write( "\t".join( line_out ) + "\n" )
    outf.close()
    os.unlink( tmpf )
        

@transform( stackLncRNARefcodingTSSDistance,
            regex( "(.+)/(.+).tsv.gz" ), 
            r"./\2.load" )
def loadStackedLncRNARefcodingTSSDistance( infile, outfile ):
    P.load( infile, outfile )


#################################################################################
## subsection: calculate distance from lncRNA TSS to nearest cage peak
#################################################################################
@transform( [ findLncRNAGeneTSS, findLncRNATranscriptTSS ], 
            regex( "(.+)/(.+)_(gene|transcript)_tss.bed.gz" ),
            add_inputs( liftOverRobustPeaks ), 
            r"\1/\2_\3_robust_cage_distance.tsv.gz" )
def calcLncRNATSSRobustCageDistance( infiles, outfile ):
    """
    Distance between lncRNA TSS and nearest robust cage peak
    """
    lnc_tss, cage_peaks = infiles
    tmpf = P.getTempFilename( "/ifs/scratch" )

    statement = ( "bedtools closest"
                  " -a %(lnc_tss)s"
                  " -b %(cage_peaks)s"
                  " -s" # enforce strandedness
                  " -d" # report distance of b from a as extra column
                  " > %(tmpf)s"
                  " 2> %(outfile)s.log" )
    P.run()

    P10.postProcessClosestBed( tmpf, outfile )
    os.unlink( tmpf )


@transform( [ findLncRNAGeneTSS, findLncRNATranscriptTSS ], 
            regex( "(.+)/(.+)_(gene|transcript)_tss.bed.gz" ),
            add_inputs( liftOverRobustPeaks ), 
            r"\1/\2_\3_robust_cage_Nulldistance.tsv.gz" )
def calcLncRNATSSRobustCageNullDistance( infiles, outfile ):
    """
    Null distribution for distance between lncNRA TSS and nearest robust cage 
    peak calculated as distance to nearest cage peak on oposite strand.
    """
    lnc_tss, cage_peaks = infiles
    tmpf = P.getTempFilename( "/ifs/scratch" )

    statement = ( "bedtools closest"
                  " -a %(lnc_tss)s"
                  " -b %(cage_peaks)s"
                  " -S" # enforce opposite strandedness
                  " -d" # report distance of b from a as extra column
                  " > %(tmpf)s"
                  " 2> %(outfile)s.log" )
    P.run()

    P10.postProcessClosestBed( tmpf, outfile )
    os.unlink( tmpf )


@transform( [ findLncRNAGeneTSS, findLncRNATranscriptTSS ], 
            regex( "(.+)/(.+)_(gene|transcript)_tss.bed.gz" ),
            add_inputs( collateNonCodingCagePeaks_robust ), 
            r"\1/\2_\3_robustnc_cage_distance.tsv.gz" )
def calcLncRNATSSRobustCageDistance_nc( infiles, outfile ):
    """
    Distance between lncRNA TSS and nearest robust cage peak not 
    assigned to an ensembl protein coding gene.
    """
    lnc_tss, cage_peaks = infiles
    tmpf = P.getTempFilename( "/ifs/scratch" )

    statement = ( "bedtools closest"
                  " -a %(lnc_tss)s"
                  " -b %(cage_peaks)s"
                  " -s" # enforce strandedness
                  " -d" # report distance of b from a as extra column
                  " > %(tmpf)s"
                  " 2> %(outfile)s.log" )
    P.run()

    P10.postProcessClosestBed( tmpf, outfile )
    os.unlink( tmpf )


@transform( [ findLncRNAGeneTSS, findLncRNATranscriptTSS ], 
            regex( "(.+)/(.+)_(gene|transcript)_tss.bed.gz" ),
            add_inputs( collateNonCodingCagePeaks_robust ), 
            r"\1/\2_\3_robustnc_cage_Nulldistance.tsv.gz" )
def calcLncRNATSSRobustCageNullDistance_nc( infiles, outfile ):
    """
    Null distribution for distance between lncNRA TSS and nearest robust cage 
    peak calculated as distance to nearest non-coding cage peak on oposite strand.
    """
    lnc_tss, cage_peaks = infiles
    tmpf = P.getTempFilename( "/ifs/scratch" )

    statement = ( "bedtools closest"
                  " -a %(lnc_tss)s"
                  " -b %(cage_peaks)s"
                  " -S" # enforce opposite strandedness
                  " -d" # report distance of b from a as extra column
                  " > %(tmpf)s"
                  " 2> %(outfile)s.log" )
    P.run()

    P10.postProcessClosestBed( tmpf, outfile )
    os.unlink( tmpf )


@transform( [ findLncRNAGeneTSS, findLncRNATranscriptTSS ], 
            regex( "(.+)/(.+)_(gene|transcript)_tss.bed.gz" ),
            add_inputs( liftOverPermissivePeaks ), 
            r"\1/\2_\3_permissive_cage_distance.tsv.gz" )
def calcLncRNATSSPermissiveCageDistance( infiles, outfile ):
    lnc_tss, cage_peaks = infiles
    tmpf = P.getTempFilename( "/ifs/scratch" )

    statement = ( "bedtools closest"
                  " -a %(lnc_tss)s"
                  " -b %(cage_peaks)s"
                  " -s" # enforce strandedness
                  " -d" # report distance of b from a as extra column
                  " > %(tmpf)s"
                  " 2> %(outfile)s.log" )
    P.run()

    P10.postProcessClosestBed( tmpf, outfile )
    os.unlink( tmpf )


@transform( [ findLncRNAGeneTSS, findLncRNATranscriptTSS ], 
            regex( "(.+)/(.+)_(gene|transcript)_tss.bed.gz" ),
            add_inputs( liftOverPermissivePeaks ), 
            r"\1/\2_\3_permissive_cage_Nulldistance.tsv.gz" )
def calcLncRNATSSPermissiveCageNullDistance( infiles, outfile ):
    lnc_tss, cage_peaks = infiles
    tmpf = P.getTempFilename( "/ifs/scratch" )

    statement = ( "bedtools closest"
                  " -a %(lnc_tss)s"
                  " -b %(cage_peaks)s"
                  " -S" # enforce opposite strandedness
                  " -d" # report distance of b from a as extra column
                  " > %(tmpf)s"
                  " 2> %(outfile)s.log" )
    P.run()

    P10.postProcessClosestBed( tmpf, outfile )
    os.unlink( tmpf )

@transform( [ findLncRNAGeneTSS, findLncRNATranscriptTSS ], 
            regex( "(.+)/(.+)_(gene|transcript)_tss.bed.gz" ),
            add_inputs( collateNonCodingCagePeaks_permissive ), 
            r"\1/\2_\3_permissivenc_cage_distance.tsv.gz" )
def calcLncRNATSSPermissiveCageDistance_nc( infiles, outfile ):
    """
    Distance between lncRNA TSS and nearest robust cage peak not 
    assigned to an ensembl protein coding gene.
    """
    lnc_tss, cage_peaks = infiles
    tmpf = P.getTempFilename( "/ifs/scratch" )

    statement = ( "bedtools closest"
                  " -a %(lnc_tss)s"
                  " -b %(cage_peaks)s"
                  " -s" # enforce strandedness
                  " -d" # report distance of b from a as extra column
                  " > %(tmpf)s"
                  " 2> %(outfile)s.log" )
    P.run()

    P10.postProcessClosestBed( tmpf, outfile )
    os.unlink( tmpf )


@transform( [ findLncRNAGeneTSS, findLncRNATranscriptTSS ], 
            regex( "(.+)/(.+)_(gene|transcript)_tss.bed.gz" ),
            add_inputs( collateNonCodingCagePeaks_permissive ), 
            r"\1/\2_\3_permissivenc_cage_Nulldistance.tsv.gz" )
def calcLncRNATSSPermissiveCageNullDistance_nc( infiles, outfile ):
    """
    Null distribution for distance between lncNRA TSS and nearest robust cage 
    peak calculated as distance to nearest non-coding cage peak on oposite strand.
    """
    lnc_tss, cage_peaks = infiles
    tmpf = P.getTempFilename( "/ifs/scratch" )

    statement = ( "bedtools closest"
                  " -a %(lnc_tss)s"
                  " -b %(cage_peaks)s"
                  " -S" # enforce opposite strandedness
                  " -d" # report distance of b from a as extra column
                  " > %(tmpf)s"
                  " 2> %(outfile)s.log" )
    P.run()

    P10.postProcessClosestBed( tmpf, outfile )
    os.unlink( tmpf )


@collate( [ calcLncRNATSSRobustCageDistance, 
            calcLncRNATSSPermissiveCageDistance,
            calcLncRNATSSRobustCageNullDistance,
            calcLncRNATSSPermissiveCageNullDistance ], 
          regex( "(.+)/(.+)_(robust|permissive)_cage_(distance|Nulldistance).tsv.gz" ), 
          r"\1/lncRNA_TSS_cage_distance.tsv.gz" )
def stackLncRNACageDistance( infiles, outfile ):
    tables = " ".join( infiles )
    to_cluster = False
    tmpf = P.getTempFilename(".")
    statement = ( "python %(scriptsdir)s/combine_tables.py"
                  "  --take=10,16" 
                  "  --cat=CAT"
                  "  --log=%(outfile)s.log"
                  " %(tables)s |"
                  " cut -f1,5,11,17"
                  " > %(tmpf)s" )
    P.run()

    # split the filename into additional columns
    outf = IOTools.openFile( outfile, "w" )
    outf.write( "interval\tcage_peak_type\tmodel\tlnc_id\tcage_id\tdistance\n" )
    header = True
    for line in IOTools.openFile( tmpf ):
        if header: 
            header = False
            continue
        line = line.split()
        # sort filename
        filename = os.path.basename( line[0] )[:-7]
        filename = filename.split("_")
        line_out = [ filename[1], filename[2], filename[4] ] + line[1:]
        outf.write( "\t".join( line_out ) + "\n" )
    outf.close()
    os.unlink( tmpf )


@collate( [ calcLncRNATSSRobustCageDistance_nc, 
            calcLncRNATSSPermissiveCageDistance_nc,
            calcLncRNATSSRobustCageNullDistance_nc,
            calcLncRNATSSPermissiveCageNullDistance_nc ], 
          regex( "(.+)/(.+)_(robustnc|permissivenc)_cage_(distance|Nulldistance).tsv.gz" ), 
          r"\1/lncRNA_TSS_noncoding_cage_distance.tsv.gz" )
def stackLncRNACageDistance_nc( infiles, outfile ):
    tables = " ".join( infiles )
    to_cluster = False
    tmpf = P.getTempFilename(".")
    statement = ( "python %(scriptsdir)s/combine_tables.py"
                  "  --take=10,16" 
                  "  --cat=CAT"
                  "  --log=%(outfile)s.log"
                  " %(tables)s |"
                  " cut -f1,5,11,17"
                  " > %(tmpf)s" )
    P.run()

    # split the filename into additional columns
    outf = IOTools.openFile( outfile, "w" )
    outf.write( "interval\tcage_peak_type\tmodel\tlnc_id\tcage_id\tdistance\n" )
    header = True
    for line in IOTools.openFile( tmpf ):
        if header: 
            header = False
            continue
        line = line.split()
        # sort filename
        filename = os.path.basename( line[0] )[:-7]
        filename = filename.split("_")
        line_out = [ filename[1], filename[2], filename[4] ] + line[1:]
        outf.write( "\t".join( line_out ) + "\n" )
    outf.close()
    os.unlink( tmpf )
        

@transform( [ stackLncRNACageDistance, stackLncRNACageDistance_nc ],
            regex( "(.+)/(.+).tsv.gz" ), 
            r"./\2.load" )
def loadStackedLncRNACageDistance( infile, outfile ):
    P.load( infile, outfile )


@collate( [ calcLncRNATSSRobustCageDistance, calcLncRNATSSPermissiveCageDistance ], 
          regex( "(.+)/(.+)_(robust|permissive)_cage_distance.tsv.gz" ), 
          r"\1/\2_cage_distance.tsv.gz" )
def combineLncRNACageDistance( infiles, outfile ):
    """
    
    """
    tables = " ".join( infiles )
    to_cluster = False
    statement = ( "python %(scriptsdir)s/combine_tables.py"
                  "  --columns=4"
                  "  --take=10,16"
                  "  --add-file-prefix"
                  "  --regex-filename='lncRNA_(.+)_cage_distance.tsv.gz'"
                  "  --log=%(outfile)s.log"
                  " %(tables)s |"
                  " gzip > %(outfile)s " )
    P.run()


@collate( [ calcLncRNATSSRobustCageDistance_nc, calcLncRNATSSPermissiveCageDistance_nc ], 
          regex( "(.+)/(.+)_(robustnc|permissivenc)_cage_distance.tsv.gz" ), 
          r"\1/\2_noncoding_cage_distance.tsv.gz" )
def combineLncRNACageDistance_nc( infiles, outfile ):
    """
    
    """
    tables = " ".join( infiles )
    to_cluster = False
    statement = ( "python %(scriptsdir)s/combine_tables.py"
                  "  --columns=4"
                  "  --take=10,16"
                  "  --add-file-prefix"
                  "  --regex-filename='lncRNA_(.+)_cage_distance.tsv.gz'"
                  "  --log=%(outfile)s.log"
                  " %(tables)s |"
                  " gzip > %(outfile)s " )
    P.run()


@transform( [ combineLncRNACageDistance, combineLncRNACageDistance_nc ],
            regex( "(.+)/(.+).tsv.gz" ), 
            r"./\2.load" )
def loadLncRNACageDistance( infile, outfile ):
    P.load( infile, outfile )


@follows( loadLncRNACageDistance, loadStackedLncRNACageDistance )
def calcLncRNACageDistance():
    pass


#################################################################################
## subsection: reset lncRNA TSS based on nearby robust/permissive CAGE peaks
#################################################################################
@collate( loadLncRNACageDistance, 
          regex( ".+/lncRNA_(gene|transcript)_noncoding_cage_distance.load" ), 
          add_inputs( r"./characterize_fantom5_overlap/lncRNA_\1_tss.bed.gz",
                      collateNonCodingCagePeaks_permissive ), 
          r"characterize_fantom5_overlap/lncRNA_\1_permissive_cage_tss.bed.gz" )
def resetLncRNATSS_permissive( infiles, outfile ):
    """
    For all lncRNA gene|transcript tss with a permissive FANTOM CAGE peak
    within 500bp that HAS NOT been associated with a coding ENSEMBL annotation,
    reset tss co-ordinates to FANTOM peak summit co-ordinates
    """
    peak_dist_table, lnc_tss_bed, fantom_bed = infiles[0]

    # fetch dataframe
    peak_dist_table = P.snip( os.path.basename( peak_dist_table ), ".load" )
    statement = ( "SELECT * FROM %s" % peak_dist_table )
    
    peak_dist_df = PU.fetch_DataFrame( statement )

    # establish whether infiles contain genes or transcripts
    if re.search( "transcript", peak_dist_table ):
        assert not re.search( "gene", peak_dist_table)
        reporter = "transcript"
    elif re.search( "gene", peak_dist_table ):
        reporter = "gene"
    else:
        raise Exception( "Unrecognised table: %s" % peak_dist_table)

    # reset TSS
    P10.resetTSS( peak_dist_df, 
                  fantom_bed, 
                  lnc_tss_bed, 
                  outfile,
                  distance = 500, 
                  reporter = reporter, 
                  stringency = "permissivenc" ) # non-coding genes selected


@collate( loadLncRNACageDistance, 
          regex( ".+/lncRNA_(gene|transcript)_noncoding_cage_distance.load" ), 
          add_inputs( r"./characterize_fantom5_overlap/lncRNA_\1_tss.bed.gz",
                      collateNonCodingCagePeaks_permissive ), 
          r"characterize_fantom5_overlap/lncRNA_\1_robust_cage_tss.bed.gz" )
def resetLncRNATSS_robust( infiles, outfile ):
    """
    For all lncRNA gene|transcript tss with a robust FANTOM CAGE peak
    within 500bp that HAS NOT been associated with a coding ENSEMBL annotation,
    reset tss co-ordinates to FANTOM peak summit co-ordinates
    """
    peak_dist_table, lnc_tss_bed, fantom_bed = infiles[0]

    # fetch dataframe
    peak_dist_table = P.snip( os.path.basename( peak_dist_table ), ".load" )
    statement = ( "SELECT * FROM %s" % peak_dist_table )
    peak_dist_df = PU.fetch_DataFrame( statement )

    # establish whether infiles contain genes or transcripts
    if re.search( "transcript", peak_dist_table ):
        assert not re.search( "gene", peak_dist_table)
        reporter = "transcript"
    elif re.search( "gene", peak_dist_table ):
        reporter = "gene"
    else:
        raise Exception( "Unrecognised table: %s" % peak_dist_table)

    # reset TSS
    P10.resetTSS( peak_dist_df, 
                  fantom_bed, 
                  lnc_tss_bed, 
                  outfile,
                  distance = 500, 
                  reporter = reporter, 
                  stringency = "robustnc" )


@transform( [ resetLncRNATSS_robust, resetLncRNATSS_permissive ],
            suffix( ".bed.gz" ), 
             add_inputs( os.path.join( PARAMS[ "annotations_dir" ], 
                                       PARAMS_ANNOTATIONS[ "interface_contigs" ] ) ),
            r"_4kb.bed.gz" )
def slopResetLncRNATSS( infiles, outfile ):
    """
    Create 4kb windows centered on lncRNA gene transcription start sites
    """
    in_bed, in_contigs = infiles
    statement = ( "bedtools slop" 
                   " -b 2000"
                   " -i %(in_bed)s"
                   " -g %(in_contigs)s |"
                   " gzip > %(outfile)s" )
    P.run()


@follows( slopResetLncRNATSS )
def resetLncRNATSS():
    pass


# #################################################################################
# ## subsection: calculate significance of interval overlap using GAT
# #################################################################################
# @collate( [ liftOverRobustPeaks, liftOverPermissivePeaks ],
#           regex( "(.+)/(.+)_cage_peaks_mm10.bed.gz" ),
#           add_inputs( slopLncRNAGeneTSS,
#                       os.path.join( PARAMS[ "annotations_dir" ], 
#                                     PARAMS_ANNOTATIONS[ "interface_contigs_ungapped_bed" ] ), 
#                       os.path.join( PARAMS[ "annotations_dir" ], 
#                                     PARAMS_ANNOTATIONS[ "interface_gc_profile_bed" ] ) ),
#           r"\1/lncRNA_gene_\2_genomic_association.tsv.gz" )
# def runGATForFANTOMPeaks( infiles, outfile ):
#     # again @collate supplies the infiles as a nested list
#     in_cage, in_tss, contig_file, isochore_file = infiles[0]

#     # replace name in TSS intervals
#     if re.search( "permissive", in_cage ):
#         interval_name = "permissive_all"
#     elif re.search( "robust", in_cage ):
#         interval_name = "robust_all"
#     else:
#         raise Exception( "cage peaks must be labelled as either permissive or"
#                          " robust: %s" % in_cage )
#     tmp_tss = P.getTempFile( "." )
#     for line in IOTools.openFile( in_tss ):
#         line = line.split()
#         line[3] = interval_name
#         tmp_tss.write( "\t".join( line ) + "\n" )
#     tmp_tss.close()
#     tmp_tss = tmp_tss.name
          
#     jobOptions = " -l mem_free=10G"
#     statement = ( "gat-run.py"
#                   " --segments=%(in_cage)s"
#                   " --annotations=%(tmp_tss)s"
#                   " --workspace-bed-file=%(contig_file)s"
#                   " --isochore-file=%(isochore_file)s"
#                   " --ignore-segment-tracks"
#                   " --num-samples=1000"
#                   " -v5"
#                   " --log=%(outfile)s.log |"
#                   " gzip > %(outfile)s" )
#     P.run()

#     os.unlink( tmp_tss )


# @collate( [ liftOverRobustPeaks, liftOverPermissivePeaks ],
#           regex( "(.+)/(.+)_cage_peaks_mm10.bed.gz" ),
#           add_inputs( slopLncRNAGeneTSS,
#                       os.path.join( PARAMS[ "annotations_dir" ], 
#                                     PARAMS_ANNOTATIONS[ "interface_contigs_ungapped_bed" ] ), 
#                       os.path.join( PARAMS[ "annotations_dir" ], 
#                                     PARAMS_ANNOTATIONS[ "interface_gc_profile_bed" ] ) ),
#           r"\1/lncRNA_gene_me_\2_genomic_association.tsv.gz" )
# def runGATForFANTOMPeaks_me( infiles, outfile ):
#     # again @collate supplies the infiles as a nested list
#     in_cage, in_tss, contig_file, isochore_file = infiles[0]

#     # remove se loci and rename TSS intervals
#     if re.search( "permissive", in_cage ):
#         interval_name = "permissive_me"
#     elif re.search( "robust", in_cage ):
#         interval_name = "robust_me"
#     else:
#         raise Exception( "cage peaks must be labelled as either permissive or"
#                          " robust: %s" % in_cage )
#     tmp_tss = P.getTempFile( "." )
#     for line in IOTools.openFile( in_tss ):
#         line = line.split()
#         if re.match( "LNCGm", line[3] ):
#             line[3] = interval_name
#             tmp_tss.write( "\t".join( line ) + "\n" )
#         else:
#             continue
#     tmp_tss.close()
#     tmp_tss = tmp_tss.name

#     jobOptions = " -l mem_free=10G"
#     statement = ( "gat-run.py"
#                   " --segments=%(in_cage)s"
#                   " --annotations=%(tmp_tss)s"
#                   " --workspace-bed-file=%(contig_file)s"
#                   " --isochore-file=%(isochore_file)s"
#                   " --ignore-segment-tracks"
#                   " --num-samples=1000"
#                   " -v5"
#                   " --log=%(outfile)s.log |"
#                   " gzip > %(outfile)s" )
#     P.run()
#     os.unlink( tmp_tss )


# @collate( [ liftOverRobustPeaks, liftOverPermissivePeaks ],
#           regex( "(.+)/(.+)_cage_peaks_mm10.bed.gz" ),
#           add_inputs( slopLncRNAGeneTTS,
#                       os.path.join( PARAMS[ "annotations_dir" ], 
#                                     PARAMS_ANNOTATIONS[ "interface_contigs_ungapped_bed" ] ), 
#                       os.path.join( PARAMS[ "annotations_dir" ], 
#                                     PARAMS_ANNOTATIONS[ "interface_gc_profile_bed" ] ) ),
#           r"\1/lncRNA_gene_me_tts_\2_genomic_association.tsv.gz" )
# def runGATForFANTOMPeaks_me_tts( infiles, outfile ):
#     # again @collate supplies the infiles as a nested list
#     in_cage, in_tss, contig_file, isochore_file = infiles[0]

#     # remove se loci and rename TSS intervals
#     if re.search( "permissive", in_cage ):
#         interval_name = "permissive_me_tts"
#     elif re.search( "robust", in_cage ):
#         interval_name = "robust_me_tts"
#     else:
#         raise Exception( "cage peaks must be labelled as either permissive or"
#                          " robust: %s" % in_cage )
#     tmp_tss = P.getTempFile( "." )
#     for line in IOTools.openFile( in_tss ):
#         line = line.split()
#         if re.match( "LNCGm", line[3] ):
#             line[3] = interval_name
#             tmp_tss.write( "\t".join( line ) + "\n" )
#         else:
#             continue
#     tmp_tss.close()
#     tmp_tss = tmp_tss.name

#     jobOptions = " -l mem_free=50G"
#     statement = ( "gat-run.py"
#                   " --segments=%(in_cage)s"
#                   " --annotations=%(tmp_tss)s"
#                   " --workspace-bed-file=%(contig_file)s"
#                   " --isochore-file=%(isochore_file)s"
#                   " --ignore-segment-tracks"
#                   " --num-samples=1000"
#                   " -v5"
#                   " --log=%(outfile)s.log |"
#                   " gzip > %(outfile)s" )
#     P.run()
#     os.unlink( tmp_tss )


# @collate( [ liftOverRobustPeaks, liftOverPermissivePeaks ],
#           regex( "(.+)/(.+)_cage_peaks_mm10.bed.gz" ),
#           add_inputs( os.path.join( PARAMS[ "annotations_dir" ], 
#                           PARAMS_ANNOTATIONS[ "interface_tss_gene_bed" ] ),
#                       os.path.join( PARAMS[ "annotations_dir" ], 
#                                     PARAMS_ANNOTATIONS[ "interface_contigs_ungapped_bed" ] ), 
#                       os.path.join( PARAMS[ "annotations_dir" ], 
#                                     PARAMS_ANNOTATIONS[ "interface_contigs" ] ), 
#                       os.path.join( PARAMS[ "annotations_dir" ], 
#                                     PARAMS_ANNOTATIONS[ "interface_gc_profile_bed" ] ) ),
#           r"\1/refcoding_gene_\2_genomic_association.tsv.gz" )
# def runGATForFANTOMPeaks_refcoding( infiles, outfile ):
#     # again @collate supplies the infiles as a nested list
#     in_cage, in_tss, contig_ungapped, contig_file, isochore_file = infiles[0]

#     # slop the refcoding TSS
#     tmp_tss_1 = P.getTempFilename( "." )
#     tmp_tss_2 = P.getTempFile( "." )
#     statement = ( "bedtools slop" 
#                   " -b 500"
#                   " -i %(in_tss)s"
#                   " -g %(contig_file)s"
#                   " > %(tmp_tss_1)s" )
#     P.run()
#     # replace name in TSS intervals
#     if re.search( "permissive", in_cage ):
#         interval_name = "permissive_refcoding"
#     elif re.search( "robust", in_cage ):
#         interval_name = "robust_refcoding"
#     else:
#         raise Exception( "cage peaks must be labelled as either permissive or"
#                          " robust: %s" % in_cage )
#     for line in IOTools.openFile( tmp_tss_1 ):
#         line = line.split()
#         line[3] = interval_name
#         tmp_tss_2.write( "\t".join( line ) + "\n" )
#     tmp_tss_2.close()
#     tmp_tss_2 = tmp_tss_2.name
          
#     jobOptions = " -l mem_free=100G"
#     statement = ( "gat-run.py"
#                   " --segments=%(in_cage)s"
#                   " --annotations=%(tmp_tss_2)s"
#                   " --workspace-bed-file=%(contig_ungapped)s"
#                   " --isochore-file=%(isochore_file)s"
#                   " --ignore-segment-tracks"
#                   " --num-samples=1000"
#                   " -v5"
#                   " --log=%(outfile)s.log |"
#                   " gzip > %(outfile)s" )
#     P.run()

#     os.unlink( tmp_tss_1 )
#     os.unlink( tmp_tss_2 )

# @merge( [ runGATForFANTOMPeaks, 
#           runGATForFANTOMPeaks_me,
#           runGATForFANTOMPeaks_me_tts,
#           runGATForFANTOMPeaks_refcoding ],
#         "./intervals_fantom_cage/lncRNA_gene_CAGE_GAT_results.tsv.gz" )
# def combineGATFANTOMPeakResults( infiles, outfile ):
#     infiles = " ".join( infiles )
#     statement = ( "python %(scriptsdir)s/combine_tables.py"
#                   "  --cat=analysis_id"
#                   "  --ignore-titles"
#                   "  --log=%(outfile)s.log"
#                   "  %(infiles)s |"
#                   " awk 'NF > 1' |" # Hack, combine tables is outputting the infiles at the start of the file
#                   " gzip > %(outfile)s" )
#     P.run()


# @transform( combineGATFANTOMPeakResults, 
#             regex( "(.+)/(.+).tsv.gz" ), 
#             r"\2.load" )
# def loadGATFANTOMPeakResults( infile, outfile ):
#     P.load( infile, outfile, options="--ignore-column=analysis_id" )

#################################################################################
## subsection: compare follicular expression with Mouse CD19+ Bcells, donor 1
#################################################################################
# extract follicular read counts, RLE normalize using edgeR... return mean
# extract lncRNA ids for those with robust/permissive CAGE support <500bp (lncRNA_TSS_noncoding_cage_distance)
# extract the RLE normalised cage TPM data for Mouse CD19+ Bcells, donor1 load CAGE_id and TPM
# generate dataframe of cage_id, tpm, foll_readcounts, cage_peak_type
@follows( mkdir( "intervals_fantom_cage" ) )
@transform( os.path.join( PARAMS["location_external_datafiles"], 
                          "FANTOM5", 
                          "mm9.cage_peak_tpm.osc.txt.gz" ),
            regex( "(?:.+)/mm9.(.+).osc.txt.gz" ), # ignore group
            r"./intervals_fantom_cage/\1.tsv.gz" )
def extractBcellCAGECounts( infile, outfile ):
    """
    Retrieve the RLE normalized CAGE peak counts for the single B cell sample 
    in the mouse FANTOM5 data... This is $54
    Removed: 01STAT:MAPPED 4993297    02STAT:NORM_FACTOR 1.1936188418439
    """
    # filter commented lines out of the CAGE TPM file
    to_cluster = False
    statement = ( '''
                  zcat %(infile)s |
                  awk '$0 !~ /^#/' |
                  awk '
                  BEGIN {target = "tpm.Mouse%%20CD19%%2b%%20B%%20Cells%%2c%%20donor1.CNhs13531.11856-125A2"};
                  NR==1 { for(i=1;i<=NF;i++){ if($i==target){break} } };
                  { print $1,$i }' |
                  sed 2d |
                  sed 2d |
                  gzip > %(outfile)s
                  ''' )  
    # print statement
    P.run()


@jobs_limit( 1, "RGlobalEnv" )
@transform( summarizeFeatureCounts, 
            regex( "(?:.+)/lncRNA_refcoding_raw_counts.tsv.gz" ),
            r"./intervals_fantom_cage/follicular_raw_counts_rle.tsv.gz" )
def normalizeFollicularRawCounts( infile, outfile ):
    """
    Extract follicular raw read counts and RLE normalize using edgeR
    """
    table = P.snip( os.path.basename( infile ), ".tsv.gz" )
    outfile = P.snip( outfile, ".gz" )

    statement = ( "SELECT"
                  "  gene_id,"
                  "  Bcell_follicular_R11,"
                  "  Bcell_follicular_R12,"
                  "  Bcell_follicular_R13," 
                  "  Bcell_follicular_R14,"
                  "  Bcell_follicular_R15"
                  " FROM %s" % table )
    df = PU.fetch_DataFrame( statement )
    df = df.set_index( "gene_id" )
    weights, rle_df = P10.calcRLENormalisedCountData( df )

    # write various outfiles
    outf_weights = IOTools.openFile( P.snip( outfile, "_rle.tsv" ) 
                                     + "_weights.tsv.gz", "w" )
    outf_weights.write( "\t".join( [str(x) for x in weights] ) )
    outf_weights.close()

    outf_untrans = P.snip( outfile, "_rle.tsv" ) + ".tsv"
    df.to_csv( outf_untrans, sep="\t" )

    rle_df.to_csv( outfile, sep="\t" )

    to_cluster = False
    statement = ( "gzip %(outfile)s;"
                  " gzip %(outf_untrans)s" )
    P.run()


@transform( normalizeFollicularRawCounts, 
            suffix( "_rle.tsv.gz" ),
            "_rle_summary.tsv.gz" )
def averageFollicularRawCounts( infile, outfile ):
    outf = IOTools.openFile( outfile, "w" )
    outf.write( "gene_id\tmean\tmedian\tmax\tmin\n" )

    header = True
    for line in IOTools.openFile( infile ):
        if header:
            header = False
            continue
        line = line.split()
        values = [ float(x) for x in line[1:] ]
        m = P10.robust_mean( values )
        md = P10.robust_median( values )
        mx = max( values )
        mn = min( values )
        line_out = [ line[0], m, md, mx, mn ]
        outf.write( "\t".join( [str(x) for x in line_out] ) + "\n" )
    outf.close()



@merge( [ stackLncRNACageDistance_nc, 
          extractBcellCAGECounts,
          averageFollicularRawCounts ],
        "./intervals_fantom_cage/follicular_cage_comparison_permissive.tsv.gz" )
def compareFollicularVSCAGEExpression_permissive( infiles, outfile ):
    # read in three dataframes
    # join the three dataframes
    # write the output. 
    mapping, cage_counts, follicular_counts = infiles
    outf = P.snip( outfile, ".gz" )
    mapping = P.snip( os.path.basename(mapping), ".tsv.gz" )

    # first read in the cage:gene_id mapping for overlapping intervals
    statement = ( "SELECT lnc_id,cage_id,cage_peak_type"
                  " FROM %(mapping)s"
                  " WHERE model LIKE 'dist%%'" # WARNING won't recognize 'distance'
                  " AND distance<501"
                  " AND cage_peak_type LIKE 'permissivenc'"
                  " AND interval LIKE 'gene'" % locals() )

    df_mapping = PU.fetch_DataFrame( statement )
#    df_mapping.to_csv( "TEMP_permissive", sep="\t" )

    df_mapping.columns = [ "gene_id", "cage_id", "cage_peak_type" ]
    
    # second read in CAGE peak counts
    df_cage_counts = pd.read_table( cage_counts, 
                                    compression = "gzip", 
                                    sep="\s" )
    df_cage_counts.columns= ["cage_id", "cage_coverage"]

    # third, read in the follicular counts
    df_fol_counts = pd.read_table( follicular_counts, 
                                   compression = "gzip", 
                                   sep="\t" )

    # merge df_mapping with df_cage_counts
    df_cage = pd.DataFrame.merge( df_mapping, 
                                  df_cage_counts, 
                                  left_on = "cage_id", 
                                  right_on = "cage_id",
                                  how = "inner" )

    # merge df_cage with df_fol_counts
    df = pd.DataFrame.merge( df_cage, 
                             df_fol_counts, 
                             left_on = "gene_id", 
                             right_on = "gene_id", 
                             how = "inner" )

    df.to_csv( outf, sep = "\t", index=False )
    to_cluster = False
    statement = "gzip -f %(outf)s"
    P.run()


@merge( [ stackLncRNACageDistance_nc, 
          extractBcellCAGECounts,
          averageFollicularRawCounts ],
        "./intervals_fantom_cage/follicular_cage_comparison_robust.tsv.gz" )
def compareFollicularVSCAGEExpression_robust( infiles, outfile ):
    # read in three dataframes
    # join the three dataframes
    # write the output. 
    mapping, cage_counts, follicular_counts = infiles
    outf = P.snip( outfile, ".gz" )
    mapping = P.snip( os.path.basename(mapping), ".tsv.gz" )

    # first read in the cage:gene_id mapping for overlapping intervals
    statement = ( "SELECT lnc_id,cage_id,cage_peak_type"
                  " FROM %(mapping)s"
                  " WHERE model LIKE 'dist%%'" # WARNING won't recognize 'distance'
                  " AND distance<501"
                  " AND cage_peak_type LIKE 'robustnc'"
                  " AND interval LIKE 'gene'" % locals() )

    df_mapping = PU.fetch_DataFrame( statement )
#    df_mapping.to_csv( "TEMP_robust", sep="\t" )


    df_mapping.columns = [ "gene_id", "cage_id", "cage_peak_type" ]
    
    # second read in CAGE peak counts
    df_cage_counts = pd.read_table( cage_counts, 
                                    compression = "gzip", 
                                    sep="\s" )
    df_cage_counts.columns= ["cage_id", "cage_coverage"]

    # third, read in the follicular counts
    df_fol_counts = pd.read_table( follicular_counts, 
                                   compression = "gzip", 
                                   sep="\t" )

    # merge df_mapping with df_cage_counts
    df_cage = pd.DataFrame.merge( df_mapping, 
                                  df_cage_counts, 
                                  left_on = "cage_id", 
                                  right_on = "cage_id",
                                  how = "inner" )

    # merge df_cage with df_fol_counts
    df = pd.DataFrame.merge( df_cage, 
                             df_fol_counts, 
                             left_on = "gene_id", 
                             right_on = "gene_id", 
                             how = "inner" )

    df.to_csv( outf, sep = "\t", index=False )
    to_cluster = False
    statement = "gzip -f %(outf)s"
    P.run()


@transform( [ compareFollicularVSCAGEExpression_permissive, 
              compareFollicularVSCAGEExpression_robust ],
            regex( "(?:.+)/(.+).tsv.gz" ), 
            r"./\1.load" )
def loadFollicularVSCAGEExpression( infile, outfile ):
    P.load( infile, outfile, options="--add-index=gene_id,cage_id" )


# #################################################################################
# ## subsection: summarize lncRNA FANTOM5 comparison
# #################################################################################

# @follows( collateNonCodingCagePeaks, 
#           findLncRNATSS, 
#           calcLncRNACageDistance,
#           resetLncRNATSS )
# #          loadGATFANTOMPeakResults )
# def runLncRNAFANTOM5Comparison():
#     pass



#################################################################################
#################################################################################
#### METASECTION #### Define LncRNAs as Enhancer or Promotor-like ####
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
# Section: Characterize lncRNAs as either eRNA or pRNA
#################################################################################
# normalize bamfiles for merged H3K4me3 H3K4me1 bams
 # calculate me1 me3 coverage across tss intervals using normalised bams
 # output to flatfile for histograms. 
# calculate ratio for expressed lncRNAs
#################################################################################
## subsection: normalise merged filtered ChIP-seq bamfiles
#################################################################################
# NB. When ruffus outputs a list or tuple as the outputs for a task, it doesn't
# seem possible to us either regex() or suffix() to specify a single element of
# that list/tuple as the input to a subsequent task (i.e. when using @transform).
# This applies both to @collate and @files when using a custom generator. 
# Using @collate... as it is most parsimonious solution. 

@follows( mkdir( "characterize_chromatin_state/bamfiles_normalized" ) )
@collate( os.path.join( PARAMS[ "location_chipseq_merged" ], "*.bam" ), 
          regex( r"(.+)/(.+)-(K4|K4me1)-R0_deduped.bwa.bam" ),
          [ r"./characterize_chromatin_state/bamfiles_normalized/\2-K4me3-R0_normalized.bam", 
            r"./characterize_chromatin_state/bamfiles_normalized/\2-K4me1-R0_normalized.bam" ] )
def normalizeK4Bamfiles( infiles, outfiles ):
    """
    Take bamfiles containing pooled samples that passed initial QC steps 
    (output by pipeline_proj010_chipseq.py) and, for each cell type, downsample larger of the 
    two chromatin mark files to match the size of the smaller.
    """
    # specify the two outfiles
    outf_me3, outf_me1 = outfiles

    # find the infiles
    for infile in infiles:
        if re.search( "K4me1", infile ):
            inf_me1 = infile
        else: 
            inf_me3 = infile
    
    E.info( "Files to normalize are as follows: "
            "\n\tinf_me1: %s "
            "\n\tinf_me3: %s "
            "\n\toutf_me1: %s "
            "\n\toutf_me3: %s" % ( inf_me1, inf_me3, outf_me1, outf_me3 ) )
    
    P10.normalizeBamfiles( inf_me1, inf_me3, outf_me1, outf_me3, submit=True )

#################################################################################
## subsection: subset only intergenic lncRNA transcripts
#################################################################################

@follows( mkdir( "characterize_chromatin_state" ) )
@transform( [ slopLncRNATSS_4kb, slopResetLncRNATSS ],
            regex( ".+/(.+)_(gene|gene_robust_cage)_tss_4kb.bed.gz" ),
            add_inputs( classifyMergedLncRNAs ), 
            r"characterize_chromatin_state/\1_\2_intergenic_tss_4kb.bed.gz" )
def subsetIntergenicLncRNATSS( infiles, outfile ):
    """
    Select only lncRNAs that are annotated as intergenic, and outputs
    bedfiles containing the 4kb windows arount TSS.
    ONLY done for original gene models and robust cage adjusted gene models, 
    otherwise there were too many results. 
    """
    lncRNA_tss, lnc_annotation_gtf = infiles

    # return a list of intergenic gene_ids/transcript_ids, depending on infile
    if re.search( "lncRNA_gene", os.path.basename( lncRNA_tss ) ):
        interval = "gene"
    elif re.search( "lncRNA_transcript", os.path.basename( lncRNA_tss ) ):
        interval = "transcript"
    else:
        raise Exception( "Unrecognised infile: %s" % lncRNA_tss )

    intergenic_interval_list = P10.returnSourceIntervals( lnc_annotation_gtf, 
                                                          "intergenic",
                                                          interval )

    # iterate through bedfile and write out only intergenic intervals
    outf = IOTools.openFile( outfile, "w" )
    for bed in Bed.iterator( IOTools.openFile( lncRNA_tss ) ):
        if bed.fields[0] in intergenic_interval_list:
            outf.write( str( bed ) + "\n" )
        else:
            continue
    outf.close()


@transform( subsetIntergenicLncRNATSS,
            suffix( ".bed.gz" ),
            ".gtf.gz" )
def createIntergenicLncRNATSSGTF( infile, outfile ):
    """
    Transform bedfile into GTF file for differential expression testing.
    """
    statement = ( "zcat %(infile)s |"
                  " python %(scriptsdir)s/bed2gff.py"
                  " --as-gtf"
                  " --log=/dev/null |"
                  " gzip > %(outfile)s" )
    P.run()

#################################################################################
## subsection: quantify H3K4me1|H3K4me3 read coverage across lncRNA TSS
#################################################################################
@follows( normalizeK4Bamfiles,
          mkdir( "./characterize_chromatin_state/coverage_gene" ), 
#          mkdir( "./characterize_lncrna_tss/coverage_transcript" ),
          mkdir( "./characterize_chromatin_state/coverage_gene_robust_cage" ) )
#          mkdir( "./characterize_lncrna_tss/coverage_gene_permissive_cage" ),
#          mkdir( "./characterize_lncrna_tss/coverage_transcript_robust_cage" ),
#          mkdir( "./characterize_lncrna_tss/coverage_transcript_permissive_cage" ) )
@product( subsetIntergenicLncRNATSS,
          formatter(".+/lncRNA_(?P<TSSSTATUS>.+)_intergenic_tss_4kb.bed.gz"),
          "./characterize_chromatin_state/bamfiles_normalized/*bam", 
          formatter("(?P<WORKINGDIR>.+)/"
                    "(?P<SUBDIR>.+)/"
                    "(.+)/"
                    "(?P<TISSUE>.+)-(?P<CONDITION>K4me1|K4me3)-R0_normalized.bam" ),
          "{WORKINGDIR[1][0]}/"
          "{SUBDIR[1][0]}/"
          "coverage_{TSSSTATUS[0][0]}/"
          "{TISSUE[1][0]}_{CONDITION[1][0]}_{TSSSTATUS[0][0]}.coverage.bed.gz" )
def calculateK4TSSCoverage( infiles, outfile ):
    tss_bed, K4_bam = infiles
    statement = ( "coverageBed"
                  " -abam %(K4_bam)s"
                  " -b %(tss_bed)s"
                  " | gzip > %(outfile)s" )
    P.run()


@collate( calculateK4TSSCoverage, 
          regex( "(.+)/(.+?)_(.+?)_(.+).coverage.bed.gz" ),
          r"\1/\2_\4.coverage_ratio.tsv.gz" )
def collateK4TSSCoverage( infiles, outfile ):
    """
    Get ChIP read coverage across TSS for each tissue type.
    """
    pseudocount = PARAMS["eRNA_pseudocount"]
    min_cov = PARAMS["eRNA_min_coverage"]
    P.submit( "/ifs/devel/projects/proj010/PipelineProj010",
              "collateK4CoverageBed",
              params= map( str, [ pseudocount, min_cov ] ),
              infiles = infiles,
              outfiles = outfile )


@collate( collateK4TSSCoverage, 
          regex( "(.+)/(.+?)_(.+).coverage_ratio.tsv.gz" ), 
          r"\1/tss_coverage_\3.tsv.gz" )
def combineK4TSSCoverage( infiles, outfile ):
    infiles = " ".join( infiles )
    statement = ( "python %(scriptsdir)s/combine_tables.py"
                  " --columns=1"
                  " --add-file-prefix"
                  " --regex-filename='(.+?)_'"
                  " --log=%(outfile)s.log"
                  " %(infiles)s |"
                  " gzip > %(outfile)s" )
    P.run()


@collate( collateK4TSSCoverage, 
          regex( "(.+)/(.+?)_(.+).coverage_ratio.tsv.gz" ), 
          r"\1/tss_coverage_\3_stacked.tsv.gz" )
def stackK4TSSCoverage( infiles, outfile ):
    tmpf = P.getTempFilename( "." )
    infiles = " ".join( infiles )
    statement = ( "python %(scriptsdir)s/combine_tables.py"
                  " --cat=cell_type"
#                  " --regex-filename='([^/].+?)_'"
                  " --log=%(outfile)s.log"
                  " %(infiles)s "
                  " > %(tmpf)s" )
    P.run()

    # hack... I can't separate the cell-type from the path
    outf = IOTools.openFile( outfile, "w" )
    header = True
    for line in IOTools.openFile( tmpf ):
        if header:
            outf.write( line )
            header = False
            continue
        line = line.split()
        cell_type = os.path.basename( line[0] ).split("_")[0]
        remainder = "\t".join( line[1:] )
        outf.write( cell_type + "\t" + remainder + "\n" )
    outf.close()
    os.unlink( tmpf )


@transform( stackK4TSSCoverage, 
            regex( "(.+)/characterize_chromatin_state/(.+)/(.+).tsv.gz" ),
            r"\1/\3.load" )
def loadK4TSSCoverage_stacked( infile, outfile ):
    P.load( infile, outfile )


@transform( combineK4TSSCoverage, 
            regex( "(.+)/characterize_chromatin_state/(.+)/(.+).tsv.gz" ),
            r"\1/\3.load" )
def loadK4TSSCoverage( infile, outfile ):
    P.load( infile, outfile )


@follows(loadK4TSSCoverage_stacked, 
         loadK4TSSCoverage )
def calculateChIPTSSCoverage():
    pass


############ This part of the pipeline is not currently implemented ########
# Looking at K4 coverage across transcript models did not appreciably alter 
# any of the downstream results. There were too many different options outputs
# for this part of the pipeline, so the transcript level analysis has been 
# commented out. 
# @follows( loadK4TSSCoverage_stacked )
# @transform( stackK4TSSCoverage, 
#             regex( "(.+)/tss_coverage_transcript_(.+).tsv.gz" ), 
#             add_inputs( loadMappedIDs ),
#             r"\1/tss_coverage_transcript_\2_highest.tsv.gz" )
# def findHighestTranscriptK4TSSCoverage( infiles, outfile ):
#     """
#     For each gene find transcript with highest coverage
#     """
#     coverage_table = P.snip( os.path.basename( infiles[0] ), ".tsv.gz" )
#     id_map_table = P.snip( os.path.basename( infiles[1] ), ".load" )
 
#     # select cell types
#     statement = ( "SELECT DISTINCT cell_type FROM %(coverage_table)s" % locals() )
#     cell_types = PU.fetch( statement )
#     cell_types = [ str( x[0] ) for x in cell_types ]  
    
#     # outfile
#     outf = IOTools.openFile( outfile, "w" )
#     outf.write( "gene_id\tcell_type\tlnc_id\tme1_cov\tme3_cov\tratio\tpass\n" )

#     # for each cell type, generate dictionary of tss coverage for each transcript
#     for cell_type in cell_types:
#         gene_dict = collections.defaultdict( list )
#         E.info( "Selecting %s transcript TSS for cell type:"
#                 " %s" % (coverage_table, cell_type) )
#         statement = ( "SELECT a.gene_id," #0
#                              "b.cell_type," #1
#                              "b.lnc_id," #2
#                              "b.me1_cov," #3
#                              "b.me3_cov," #4
#                              "b.ratio," #5
#                              "b.pass" #6
#                       " FROM %(id_map_table)s AS a"
#                       " INNER JOIN %(coverage_table)s AS b"
#                       " ON a.transcript_id=b.lnc_id"
#                       " WHERE b.cell_type='%(cell_type)s'"
#                       " ORDER BY b.lnc_id" % locals() )
#         result = PU.fetch( statement )

#         for entry in result: 
#             gene_id = str( entry[0] )
#             tss = [ str(x) for x in entry ] 
#             gene_dict[ gene_id ].append( tss )
            
#         # iterate through genedict and write out the tss_entry with highest coverage
#         # NB. in case of ties, returns the first tss entry with the highest coverage
#         # (Transcript IDs were named using sorted gtf)
#         E.info( "Writing %s transcript TSS to outfile for cell type:"
#                 " %s" % (coverage_table, cell_type) )
#         for gene_id in gene_dict.keys():
#             # find tss with highest me3 coverage
#             me3_coverage = [ x[3] for x in gene_dict[ gene_id ] ] 
#             me3_highest_loc = me3_coverage.index( max( me3_coverage ) )
#             me3_highest_tss = gene_dict[ gene_id ][ me3_highest_loc ]
#             # find tss with highest me1 coverage
#             me1_coverage = [ x[4] for x in gene_dict[ gene_id ] ]
#             me1_highest_loc = me1_coverage.index( max( me1_coverage ) )
#             me1_highest_tss = gene_dict[ gene_id ][ me1_highest_loc ]

#             # output tss with highest absolute coverage
#             if me3_highest_tss[2] == me1_highest_tss[2]:
#                 outf.write( "\t".join( me3_highest_tss ) + "\n" )
#             elif me3_highest_tss[4] > me1_highest_tss[3]:
#                 outf.write( "\t".join( me3_highest_tss ) + "\n" )
#             elif me1_highest_tss[3] > me3_highest_tss[4]:
#                 outf.write( "\t".join( me1_highest_tss ) + "\n" )
#             else:
#                 # in case of ties, select tss with highest fold change
#                 assert me3_highest_tss[4] == me1_highest_tss[3]
#                 if me3_highest_tss[3] < me1_highest_tss[4]:
#                     outf.write( "\t".join( me3_highest_tss ) + "\n" )
#                 elif me1_highest_tss[4] < me3_highest_tss[3]:
#                     outf.write( "\t".join( me1_highest_tss ) + "\n" )
#                 else:
#                     if me3_highest_tss[2] < me1_highest_tss[2]:
#                         outf.write( "\t".join( me3_highest_tss ) + "\n" )
#                     else:
#                         outf.write( "\t".join( me1_highest_tss ) + "\n" )
#         E.info( "Finished Writing %s transcript TSS to outfile for cell type:"
#                 " %s" % (coverage_table, cell_type) )


# @transform( findHighestTranscriptK4TSSCoverage, 
#             regex( "(.+)/(.+).tsv.gz" ), 
#             r"./\2.load" )
# def loadHighestTranscriptK4TSSCoverage( infile, outfile ):
#     P.load( infile, outfile )         


#################################################################################
## subsection: plot heatmaps of lncRNA TSS chromatin coverage
#################################################################################
# slopLncRNATSS_4kb
### intervals_fantom_cage/lncRNA_gene_tss_4kb.bed.gz
### intervals_fantom_cage/lncRNA_transcript_tss_4kb.bed.gz     

# slopResetLncRNATSS_4kb
### lncRNA_gene_permissive_cage_tss_4kb.bed.gz
### lncRNA_gene_robust_cage_tss_4kb.bed.gz
### lncRNA_transcript_permissive_cage_tss_4kb.bed.gz
### lncRNA_transcript_robust_cage_tss_4kb.bed.gz

# normalizeK4Bamfiles
### characterize_chromatin_state/bamfiles_normalized/*_normalized.bam
@follows( normalizeK4Bamfiles,
          mkdir( "./characterize_chromatin_state/peakshape_no_control_gene" ), 
#          mkdir( "./characterize_chromatin_state/peakshape_no_control_transcript" ),
          mkdir( "./characterize_chromatin_state/peakshape_no_control_gene_robust_cage" ) )
#          mkdir( "./characterize_chromatin_state/peakshape_no_control_gene_permissive_cage" ),
#          mkdir( "./characterize_chromatin_state/peakshape_no_control_transcript_robust_cage" ),
#          mkdir( "./characterize_chromatin_state/peakshape_no_control_transcript_permissive_cage" ) )
@product( subsetIntergenicLncRNATSS,
          formatter(".+/lncRNA_(?P<TSSSTATUS>.+)_intergenic_tss_4kb.bed.gz"),
          "./characterize_chromatin_state/bamfiles_normalized/*bam", 
          formatter("(?P<WORKINGDIR>.+)/"
                    "(?P<SUBDIR>.+)/"
                    "(.+)/"
                    "(?P<TISSUE>.+)-(?P<CONDITION>K4me1|K4me3)-R0_normalized.bam" ),
          "{WORKINGDIR[1][0]}/"
          "{SUBDIR[1][0]}/"
          "peakshape_no_control_{TSSSTATUS[0][0]}/"
          "{TISSUE[1][0]}_{CONDITION[1][0]}_{TSSSTATUS[0][0]}.peakshape.log" )
def calculateK4PeakShape( infiles, outfile ):
    out_stub = P.snip( outfile, ".log" )
    tss_bed, K4_bam = infiles
    
    statement = ( "python %(scriptsdir)s/bam2peakshape.py" 
                  " %(K4_bam)s"
                  " %(tss_bed)s"
                  # "  --method=sort --sort-order=unsorted" # This fails as no --method in bam2peakshape.py
                  "  --sort-order=unsorted" 
                  "  --output-filename-pattern=%(out_stub)s.%%s"
                  "  --log=%(outfile)s"
                  "  --force-output" )
    P.run()


@follows( normalizeK4Bamfiles,
          mkdir( "./characterize_chromatin_state/peakshape_uncentered_no_control_gene" ), 
#          mkdir( "./characterize_chromatin_state/peakshape_uncentered_no_control_transcript" ),
          mkdir( "./characterize_chromatin_state/peakshape_uncentered_no_control_gene_robust_cage" ) )
#          mkdir( "./characterize_chromatin_state/peakshape_uncentered_no_control_gene_permissive_cage" ),
#          mkdir( "./characterize_chromatin_state/peakshape_uncentered_no_control_transcript_robust_cage" ),
#          mkdir( "./characterize_chromatin_state/peakshape_uncentered_no_control_transcript_permissive_cage" ) )
@product( subsetIntergenicLncRNATSS,
          formatter(".+/lncRNA_(?P<TSSSTATUS>.+)_intergenic_tss_4kb.bed.gz"),
          "./characterize_chromatin_state/bamfiles_normalized/*bam", 
          formatter("(?P<WORKINGDIR>.+)/"
                    "(?P<SUBDIR>.+)/"
                    "(.+)/"
                    "(?P<TISSUE>.+)-(?P<CONDITION>K4me1|K4me3)-R0_normalized.bam" ),
          "{WORKINGDIR[1][0]}/"
          "{SUBDIR[1][0]}/"
          "peakshape_uncentered_no_control_{TSSSTATUS[0][0]}/"
          "{TISSUE[1][0]}_{CONDITION[1][0]}_{TSSSTATUS[0][0]}.peakshape.log" )
def calculateK4PeakShape_uncentered( infiles, outfile ):
    out_stub = P.snip( outfile, ".log" )
    tss_bed, K4_bam = infiles
    
    statement = ( "python %(scriptsdir)s/bam2peakshape.py" 
                  " %(K4_bam)s"
                  " %(tss_bed)s"
                  # "  --method=sort --sort-order=unsorted" # This fails as no --method in bam2peakshape.py
                  "  --sort-order=unsorted"
                  "  --centring-method=middle"
                  "  --output-filename-pattern=%(out_stub)s.%%s"
                  "  --log=%(outfile)s"
                  "  --force-output" )
    P.run()


# @transform( calculateK4PeakShape_uncentered,
#             regex( "(.+)/(.+)transcript(.*).peakshape.log" ),
#             add_inputs( r"./tss_coverage_transcript\3_stacked_highest.load" ),
#             r"\1/\2transcript\3_highest.peakshape.log" )
# def subsetTranscriptK4PeakShape_uncentered( infiles, outfile ):
#     # for cell_type of infile[0] retrieve all the lncRNA Transcript IDs in table for infile[1]
#     # iterate through infile[0], output those entries for which there is entry in list from infile[1]
#     in_peakfile = P.snip( infiles[0], ".log" ) + ".matrix_unsorted.gz"
#     transcript_table = P.snip( os.path.basename( infiles[1] ), ".load" )
#     cell_type = os.path.basename( infiles[0] ).split("_")[0]

#     # select transcript_ids from csvdb table pertaining to specified cell_type
#     statement = ( "SELECT lnc_id "
#                   "FROM %(transcript_table)s"
#                   " WHERE cell_type = '%(cell_type)s'" % locals() )
#     transcript_ids = PU.fetch( statement )
#     transcript_ids = [ str(x[0]) for x in transcript_ids ]

#     # iterate through peak file and write out lines that match transcript_ids in list
#     outf = P.snip( outfile, ".log" ) + ".matrix_unsorted.gz"
#     outf = IOTools.openFile( outf, "w" )
#     header = True
#     x = 0
#     for line in IOTools.openFile( in_peakfile ):
#         transcript_id = line.split()[0]
#         if header:
#             outf.write( line )
#             header = False
#         elif transcript_id in transcript_ids:
#             outf.write( line )
#             x += 1
#         else:
#             continue
    
#     outf.close()
#     E.info( "%s intervals written to %s" % ( str( x ), outf ) )
#     P.touch( outfile )

    
# @follows( normalizeK4Bamfiles,
#           mkdir( "./characterize_lncrna_tss/peakshape_K4me1_control_gene" ), 
#           mkdir( "./characterize_lncrna_tss/peakshape_K4me1_control_transcript" ),
#           mkdir( "./characterize_lncrna_tss/peakshape_K4me1_control_gene_robust_cage" ),
#           mkdir( "./characterize_lncrna_tss/peakshape_K4me1_control_gene_permissive_cage" ),
#           mkdir( "./characterize_lncrna_tss/peakshape_K4me1_control_transcript_robust_cage" ),
#           mkdir( "./characterize_lncrna_tss/peakshape_K4me1_control_transcript_permissive_cage" ) )
# @product( subsetIntergenicLncRNATSS,
#           formatter(".+/lncRNA_(?P<TSSSTATUS>.+)_intergenic_tss_4kb.bed.gz"),
#           "./characterize_lncrna_tss/bamfiles_normalized/*bam", 
#           formatter("(?P<WORKINGDIR>.+)/"
#                     "(?P<SUBDIR>.+)/"
#                     "(.+)/"
#                     "(?P<TISSUE>.+)-K4-R0_normalized.bam" ),
#           "{WORKINGDIR[1][0]}/"
#           "{SUBDIR[1][0]}/"
#           "peakshape_K4me1_control_{TSSSTATUS[0][0]}/"
#           "{TISSUE[1][0]}_K4_{TSSSTATUS[0][0]}.peakshape.log" )
# def calcluateK4me3PeakShapeWithK4me1Control( infiles, outfile ):
#     out_stub = P.snip( outfile, ".log" )

#     tss_bed, me3_bam = infiles
#     me1_bam = re.sub( "K4", "K4me1", me3_bam )

#     statement = ( "python %(scriptsdir)s/bam2peakshape.py" 
#                   " %(me3_bam)s"
#                   " %(tss_bed)s"
#                   "  --control-bam-file=%(me1_bam)s"
#                   "  --method=sort --sort-order=peak-height"
#                   "  --output-filename-pattern=%(out_stub)s.%%s"
#                   "  --log=%(outfile)s"
#                   "  --force-output" )
#     P.run()


# @follows( normalizeK4Bamfiles,
#           mkdir( "./characterize_lncrna_tss/peakshape_K4me3_control_gene" ), 
#           mkdir( "./characterize_lncrna_tss/peakshape_K4me3_control_transcript" ),
#           mkdir( "./characterize_lncrna_tss/peakshape_K4me3_control_gene_robust_cage" ),
#           mkdir( "./characterize_lncrna_tss/peakshape_K4me3_control_gene_permissive_cage" ),
#           mkdir( "./characterize_lncrna_tss/peakshape_K4me3_control_transcript_robust_cage" ),
#           mkdir( "./characterize_lncrna_tss/peakshape_K4me3_control_transcript_permissive_cage" ) )
# @product( subsetIntergenicLncRNATSS,
#           formatter(".+/lncRNA_(?P<TSSSTATUS>.+)_intergenic_tss_4kb.bed.gz"),
#           "./characterize_lncrna_tss/bamfiles_normalized/*bam", 
#           formatter("(?P<WORKINGDIR>.+)/"
#                     "(?P<SUBDIR>.+)/"
#                     "(.+)/"
#                     "(?P<TISSUE>.+)-K4me1-R0_normalized.bam" ),
#           "{WORKINGDIR[1][0]}/"
#           "{SUBDIR[1][0]}/"
#           "peakshape_K4me3_control_{TSSSTATUS[0][0]}/"
#           "{TISSUE[1][0]}_K4me1_{TSSSTATUS[0][0]}.peakshape.log" )
# def calcluateK4me1PeakShapeWithK4me3Control( infiles, outfile ):
#     out_stub = P.snip( outfile, ".log" )

#     tss_bed, me1_bam = infiles 
#     me3_bam = re.sub( "K4me1", "K4", me1_bam )

#     statement = ( "python %(scriptsdir)s/bam2peakshape.py" 
#                   " %(me1_bam)s"
#                   " %(tss_bed)s"
#                   "  --control-bam-file=%(me3_bam)s"
#                   "  --method=sort --sort-order=peak-height"
#                   "  --output-filename-pattern=%(out_stub)s.%%s"
#                   "  --log=%(outfile)s"
#                   "  --force-output" )
#     P.run()

# @follows( calculateK4PeakShape,
#           calcluateK4me3PeakShapeWithK4me1Control,
#           calcluateK4me1PeakShapeWithK4me3Control )
# def calculatePeakShapes(): pass


@follows( calculateChIPTSSCoverage )
@collate( [ calculateK4PeakShape, 
            calculateK4PeakShape_uncentered ],
          regex( "(.+)/peakshape(.*)_no_control_gene(.*)/(.+?)_(.+?)_(.+).peakshape.log" ),
          r"\1/peakshape\2_no_control_gene\3/\4_\6.peakshape_highCoverage.png" )
def plotLncRNATSSPeakShape_unsorted_highCoverage( infiles, outfile ):
    """
    """
    ## Specify infiles
    for infile in infiles:
        if re.search( "_K4me3_", infile ):
            in_me3 = re.sub( ".log", ".matrix_unsorted.gz", infile )
        elif re.search( "_K4me1_", infile ):
            in_me1 = re.sub( ".log", ".matrix_unsorted.gz", infile )
        else:
            raise IOError( "Unrecognised infile: %s" % infile )

    ## Find details of sqlite table that specifies lncs with sufficient coverage
    # generate name of coverage table from pathname
    subdir = os.path.abspath( infile ).split( "/" )[-2]
    if re.search( "uncentered", subdir ):
        table_name = "tss_coverage_" + "_".join( subdir.split("_")[4:] )
    else:
        table_name = "tss_coverage_" + "_".join( subdir.split("_")[3:] )

    # generate name of column to filter on
    threshold = os.path.basename( infile ).split("_")[0] + "_pass"
    ratio = os.path.basename( infile ).split("_")[0] + "_ratio"
    
    ## retrieve ordered lncRNAs (plus ratio) that pass coverage threshold
    statement = ( "SELECT lnc_id,%(ratio)s"
                  " FROM %(table_name)s"
                  " WHERE %(threshold)s=1"
                  " ORDER BY %(ratio)s" % locals() )
    dbh = sqlite3.connect( PARAMS["database"] )
    cc = dbh.cursor()
    result = cc.execute( statement ).fetchall()
    lnc_ids = [ str( x[0] ) for x in result ]
    # select a limited number of ratios to make up ylab
    ratios = [ round( float( x[1] ), 3 ) for x in reversed( result ) ]
    ratios_sub = []
    for ratio in enumerate( ratios ):
        if ratio[0] % 100 == 0:
            ratios_sub.append( ratio[1] )
        else:
            ratios_sub.append( "" )
    del ratios_sub[-1]
    ratios_sub.append( round( float( result[1][1] ), 3 ) )
    ratios_sub = [ x for x in reversed( ratios_sub ) ]

    # for the sake of the logfile... 
    tmpf1_name = P.getTempFilename(".")
    tmpf2_name = P.getTempFilename(".")
    tmpf1 = open( tmpf1_name, "wb" )
    tmpf2 = open( tmpf2_name, "wb" )
    pickle.dump( lnc_ids, tmpf1 )
    pickle.dump( ratios_sub, tmpf2 )
    tmpf1.close()
    tmpf2.close()

    ## Submit wrapper for R function
    P.submit( "/ifs/devel/projects/proj010/PipelineProj010",
              "plotMatrixPeakShape",
              params = [tmpf1_name, tmpf2_name ],
              infiles = [ in_me3, in_me1 ],
              outfiles = outfile,
              jobOptions = " -l mem_free=5G" )
    os.unlink( tmpf1_name )
    os.unlink( tmpf2_name )    


# @transform( [ calcluateK4me3PeakShapeWithK4me1Control, 
#               calcluateK4me1PeakShapeWithK4me3Control ],
#             suffix( ".log" ), 
#             "_all.png" )
# def plotLncRNATSSPeakShape_all( infile, outfile ):
#     """
#     Plot heatmaps of peak shape for matrix and control files
#     """
#     in_matrix = re.sub( ".log", ".matrix_peak_height.gz", infile )
#     in_control = re.sub( ".log", ".control_peak_height.gz", infile )
#     P.submit( "/ifs/devel/projects/proj010/PipelineProj010",
#               "plotMatrixAndControlPeakShape",
#               infiles = [ in_matrix, in_control ],
#               outfiles = outfile,
#               jobOptions = " -l mem_free=5G" )


# @follows( loadK4TSSCoverage )    
# @transform( [ calcluateK4me3PeakShapeWithK4me1Control, 
#               calcluateK4me1PeakShapeWithK4me3Control ],
#             suffix( ".log" ), 
#             "_highCoverage.png" )
# def plotLncRNATSSPeakShape_highCoverage( infile, outfile ):
#     """
#     Plot heatmaps of peak shape for matrix and control files
#     that contain only those tss for which there is adequate coverage
#     """
#     ## Find details of sqlite table that specifies lncs with sufficient coverage
#     # generate name of coverage table from pathname
#     subdir = os.path.abspath( infile ).split( "/" )[-2]
#     table_name = "tss_coverage_" + "_".join( subdir.split("_")[3:] )
#     # generate name of column to filter on
#     column = os.path.basename( infile ).split("_")[0] + "_pass"

#     ## Retrieve list of lncRNAs that pass coverage threshold
#     statement = ( "SELECT lnc_id"
#                   " FROM %(table_name)s"
#                   " WHERE %(column)s=1" % locals() )
#     dbh = sqlite3.connect( PARAMS["database"] )
#     cc = dbh.cursor()
#     result = cc.execute( statement ).fetchall()
#     result = [ str( x[0] ) for x in result ]
#     # for the sake of the logfile... 
#     tmpf_name = P.getTempFilename(".")
#     tmpf = open( tmpf_name, "wb" )
#     pickle.dump( result, tmpf )
#     tmpf.close()

#     in_matrix = re.sub( ".log", ".matrix_peak_height.gz", infile )
#     in_control = re.sub( ".log", ".control_peak_height.gz", infile )
#     P.submit( "/ifs/devel/projects/proj010/PipelineProj010",
#               "plotMatrixAndControlPeakShape",
#               params = [tmpf_name, ],
#               infiles = [ in_matrix, in_control ],
#               outfiles = outfile, 
#               jobOptions = " -l mem_free=5G" )
#     os.unlink( tmpf_name )


# @follows( calculateChIPTSSCoverage )
# @collate( subsetTranscriptK4PeakShape_uncentered, 
#           regex( "(.+)/peakshape(.*)_no_control_transcript(.*)/(.+?)_(.+?)_(.+).peakshape.log" ),
#           r"\1/peakshape\2_no_control_transcript\3/\4_\6.peakshape_highCoverage.png" )
# def plotLncRNATranscriptTSSPeakShape_unsorted_highCoverage( infiles, outfile ):
#     ## Specify infiles
#     for infile in infiles:
#         if re.search( "_K4_", infile ):
#             in_me3 = re.sub( ".log", ".matrix_unsorted.gz", infile )
#         elif re.search( "_K4me1_", infile ):
#             in_me1 = re.sub( ".log", ".matrix_unsorted.gz", infile )
#         else:
#             raise IOError( "Unrecognised infile: %s" % infile )

#     # Find details of sqlite table that specifies lncs with sufficient coverage
#     # generate name of coverage table from pathname
#     subdir = os.path.abspath( infile ).split( "/" )[-2]
#     table_name = "tss_coverage_" + "_".join( subdir.split("_")[4:] ) + "_stacked_highest"

#     ## retrieve ordered lncRNAs (plus ratio) that pass coverage threshold
#     cell_type = os.path.basename( in_me3 ).split("_")[0]
#     statement = ( "SELECT lnc_id,ratio"
#                   " FROM %(table_name)s"
#                   " WHERE pass=1 AND cell_type='%(cell_type)s'"
#                   " ORDER BY ratio" % locals() )

#     dbh = sqlite3.connect( PARAMS["database"] )
#     cc = dbh.cursor()
#     result = cc.execute( statement ).fetchall()
#     lnc_ids = [ str( x[0] ) for x in result ]
#     # select a limited number of ratios to make up ylab
#     ratios = [ round( float( x[1] ), 3 ) for x in reversed( result ) ]
#     ratios_sub = []
#     for ratio in enumerate( ratios ):
#         if ratio[0] % 100 == 0:
#             ratios_sub.append( ratio[1] )
#         else:
#             ratios_sub.append( "" )
#     del ratios_sub[-1]
#     ratios_sub.append( round( float( result[1][1] ), 3 ) )
#     ratios_sub = [ x for x in reversed( ratios_sub ) ]

#     # for the sake of the logfile... 
#     tmpf1_name = P.getTempFilename(".")
#     tmpf2_name = P.getTempFilename(".")
#     tmpf1 = open( tmpf1_name, "wb" )
#     tmpf2 = open( tmpf2_name, "wb" )
#     pickle.dump( lnc_ids, tmpf1 )
#     pickle.dump( ratios_sub, tmpf2 )
#     tmpf1.close()
#     tmpf2.close()

#     ## Submit wrapper for R function
#     P.submit( "/ifs/devel/projects/proj010/PipelineProj010",
#               "plotMatrixPeakShape",
#               params = [tmpf1_name, tmpf2_name ],
#               infiles = [ in_me3, in_me1 ],
#               outfiles = outfile,
#               jobOptions = " -l mem_free=5G" )
#     os.unlink( tmpf1_name )
#     os.unlink( tmpf2_name )    


# @follows( plotLncRNATSSPeakShape_all,
#           plotLncRNATSSPeakShape_highCoverage,
#           plotLncRNATSSPeakShape_unsorted_highCoverage,
#           plotLncRNATranscriptTSSPeakShape_unsorted_highCoverage )
# def plotPeakShapes():
#     pass


# #################################################################################
# ## subsction: rank lncRNA loci based on peak height
# #################################################################################
# @transform( [ calcluateK4me3PeakShapeWithK4me1Control,
#               calcluateK4me1PeakShapeWithK4me3Control ],
#             suffix( ".log" ), 
#             ".rank.gz" )
# def rankLncRNAPeaks( infile, outfile ):
#     """
#     Output a file containing lncRNA ID and a rank based on peak height
#     (i.e. position in *.matrix_peak_height.gz)
#     """
#     infile = P.snip( infile, ".log" ) + ".matrix_peak_height.gz"
#     outf = IOTools.openFile( outfile, "w" )
#     outf.write( "lnc_id\trank\n" )
#     rank = 0
#     header = True
#     for line in IOTools.openFile( infile ):
#         if header: 
#             header = False
#             continue
#         rank += 1
#         lnc_id = line.split()[0]
#         outf.write( lnc_id + "\t" + str( rank ) + "\n" )
#     outf.close()


# @collate( rankLncRNAPeaks,
#           regex( "(.+)/peakshape_(.+)/(.+).rank.gz" ),
#           r"\1/peakshape_\2/ranked_peaks_\2.tsv.gz" )
# def collateRankedLncRNAPeaks( infiles, outfile ):
#     infiles = " ".join( infiles )
#     statement = ( "python %(scriptsdir)s/combine_tables.py"
#                   " --columns=1"
#                   " --add-file-prefix"
#                   " --regex-filename='(.+?)_'"
#                   " --log=%(outfile)s.log"
#                   " %(infiles)s |"
#                   " gzip > %(outfile)s" )
#     P.run()


# @transform( collateRankedLncRNAPeaks,
#             regex( "(.+)/(.+)/(.+)/(.+).tsv.gz" ),
#             r"\1/\4.load" )
# def loadRankedLncRNAPeaks( infile, outfile ):
#     P.load( infile, outfile )




#################################################################################
## subsection: intersect lncRNA TSS with significant chromatin peaks
#################################################################################
@follows( mkdir( "./characterize_chromatin_peak_intersection/chipseq_peakfiles" ) )
@transform( glob.glob( os.path.join( PARAMS[ "location_chipseqfiles" ], 
                                     "*.narrowPeak.gz" ) ), 
            regex( "(.+)/(.+).narrowPeak.gz" ), 
            r"./characterize_chromatin_peak_intersection/chipseq_peakfiles/\2.bed.gz" )
def cleanNarrowPeakFiles( infile, outfile ):
    """
    Sort out peak names in narrowPeak files and convert narrowpeak to bed 6
    """
    outf = IOTools.openFile( outfile, "w" )
    for line in IOTools.openFile( infile ):
        line = line.split()
        peak_id = "_".join( line[3].split("_")[-2:] )
        line[3] = peak_id
        # convert to bed6
        line = "\t".join( line[:6] )
        outf.write( line + "\n" )
    outf.close()


@follows( mkdir( "./characterize_chromatin_peak_intersection/lncRNA_gene_peak_distance" ), 
#          mkdir( "./characterize_lncrna_chippeaks/lncRNA_gene_permissive_peak_distance" ), 
          mkdir( "./characterize_chromatin_peak_intersection/lncRNA_gene_robust_peak_distance" ) )
# introducing optional capture groups... (\w)?
#@product( [ findLncRNAGeneTSS, resetLncRNATSS_permissive, resetLncRNATSS_robust ],
@product( [ findLncRNAGeneTSS, resetLncRNATSS_robust ],
          formatter( ".+/lncRNA(.*)_gene(?P<TSS_STATUS>.*?)_(cage_)?tss.bed.gz" ), 
          cleanNarrowPeakFiles,
          formatter("(?P<WORKINGDIR>.+)/"
                    "(?P<SUBDIR>.+)/"
                    "(.+)/"
                    "(?P<TISSUE>.+)-(?P<CONDITION>K4me1|K4me3).bed.gz" ),
          "{WORKINGDIR[1][0]}/"
          "{SUBDIR[1][0]}/"
          "lncRNA_gene{TSS_STATUS[0][0]}_peak_distance/"
          "{TISSUE[1][0]}_{CONDITION[1][0]}_lncRNA_gene{TSS_STATUS[0][0]}_peak_distance.tsv.gz" )
def findNearestK4Peak_geneTSS( infiles, outfile ):
    """
    Finds distance between lncRNA tss and nearest ChIP peak. Uses @product to do
    this for orignal/FANTOM adjusted lncRNA TSS vs all ChIP experiments.
    """
    lnc_tss, chip_peaks = infiles
    statement = ( "bedtools closest"
                  " -a %(lnc_tss)s"
                  " -b %(chip_peaks)s"
                  " -t first" # report the first peak for ties
                  " -d" # report distance of b from a as extra column
                  " | gzip > %(outfile)s"
                  " 2> %(outfile)s.log" )
    P.run()


# @follows( mkdir( "./characterize_lncrna_chippeaks/lncRNA_transcript_highest" ) )
# @split( findLncRNATranscriptTSS, 
#         regex( "(.+)/lncRNA_merged_transcript_tss_bed.gz" ),
#         add_inputs( "./tss_coverage_transcript_stacked_highest.load" ),
#         r"./characterize_lncrna_chippeaks/lncRNA_transcript_highest/*_lncRNA_transcript_tss.bed.gz" )
# def findHighestTSSForEachCellType( infiles, outfiles ):
#     in_bed, in_table = infiles
#     in_table = P.snip( os.path.basename( in_table ), ".load" )
#     out_dir = "./characterize_lncrna_chippeaks/lncRNA_transcript_highest/"
#     out_suffix = "_lncRNA_transcript_tss.bed.gz"
    
#     P10.findHighestTSS( in_bed, in_table, out_dir, out_suffix )


# @follows( mkdir( "./characterize_lncrna_chippeaks/lncRNA_transcript_permissive_highest" ) )
# @split( resetLncRNATSS_permissive,
#         regex( "(.+)/lncRNA_transcript_permissive_cage_tss.bed.gz" ),
#         add_inputs( "./tss_coverage_transcript_permissive_cage_stacked_highest.load" ),
#         r"./characterize_lncrna_chippeaks/lncRNA_transcript_permissive_highest/"
#         "*_lncRNA_transcript_permissive_tss.bed.gz" )
# def findHighestPermissiveTSSForEachCellType( infiles, outfiles ):
#     in_bed, in_table = infiles
#     in_table = P.snip( os.path.basename( in_table ), ".load" )
#     out_dir = "./characterize_lncrna_chippeaks/lncRNA_transcript_permissive_highest/"
#     out_suffix = "_lncRNA_transcript_permissive_tss.bed.gz"
    
#     P10.findHighestTSS( in_bed, in_table, out_dir, out_suffix )


# @follows( mkdir( "./characterize_lncrna_chippeaks/lncRNA_transcript_robust_highest" ) )
# @split( resetLncRNATSS_robust, 
#         regex( "(.+)/lncRNA_transcript_robust_cage_tss.bed.gz" ),
#         add_inputs( "./tss_coverage_transcript_robust_cage_stacked_highest.load" ),
#         r"./characterize_lncrna_chippeaks/lncRNA_transcript_robust_highest/"
#         "*_lncRNA_transcript_robust_tss.bed.gz" )
# def findHighestRobustTSSForEachCellType( infiles, outfiles ):
#     in_bed, in_table = infiles
#     in_table = P.snip( os.path.basename( in_table ), ".load" )
#     out_dir = "./characterize_lncrna_chippeaks/lncRNA_transcript_robust_highest/"
#     out_suffix = "_lncRNA_transcript_robust_tss.bed.gz"
    
#     P10.findHighestTSS( in_bed, in_table, out_dir, out_suffix )


# @follows( mkdir( "./characterize_lncrna_chippeaks/lncRNA_transcript_highest_peak_distance/" ),
#           mkdir( "./characterize_lncrna_chippeaks/lncRNA_transcript_permissive_highest_peak_distance/" ),
#           mkdir( "./characterize_lncrna_chippeaks/lncRNA_transcript_robust_highest_peak_distance/" ) )
# @transform( [ findHighestTSSForEachCellType,
#               findHighestPermissiveTSSForEachCellType,
#               findHighestRobustTSSForEachCellType ], 
#             regex( "(.+)/(.+)/(.+?)_lncRNA_transcript(.*)_tss.bed.gz" ),
#             add_inputs( r"\1/chipseq_peakfiles/\3-K4.bed.gz" ),
#             r"\1/\2_peak_distance/\3_K4_lncRNA_transcript\4_peak_distance.tsv.gz" )
# def findNearestK4me3Peak_transcriptTSS( infiles, outfile ):
#     lnc_tss, chip_peaks = infiles
#     statement = ( "bedtools closest"
#                   " -a %(lnc_tss)s"
#                   " -b %(chip_peaks)s"
#                   " -t first" # report the first peak for ties
#                   " -d" # report distance of b from a as extra column
#                   " | gzip > %(outfile)s"
#                   " 2> %(outfile)s.log" )
#     P.run()


# @follows( findNearestK4me3Peak_transcriptTSS )
# @transform( [ findHighestTSSForEachCellType,
#               findHighestPermissiveTSSForEachCellType,
#               findHighestRobustTSSForEachCellType ], 
#             regex( "(.+)/(.+)/(.+?)_lncRNA_transcript(.*)_tss.bed.gz" ),
#             add_inputs( r"\1/chipseq_peakfiles/\3-K4me1.bed.gz" ),
#             r"\1/\2_peak_distance/\3_K4me1_lncRNA_transcript\4_peak_distance.tsv.gz" )
# def findNearestK4me1Peak_transcriptTSS( infiles, outfile ):
#     lnc_tss, chip_peaks = infiles
#     statement = ( "bedtools closest"
#                   " -a %(lnc_tss)s"
#                   " -b %(chip_peaks)s"
#                   " -t first" # report the first peak for ties
#                   " -d" # report distance of b from a as extra column
#                   " | gzip > %(outfile)s"
#                   " 2> %(outfile)s.log" )
#     P.run()


# @collate( [ findNearestK4Peak_geneTSS, 
#             findNearestK4me3Peak_transcriptTSS, 
#             findNearestK4me1Peak_transcriptTSS ], 
@collate( findNearestK4Peak_geneTSS,
          regex( "(.+)/(.+?)_(.+).tsv.gz" ), 
          r"\1/\3.tsv.gz" )
def stackTSSChIPPeakDistance( infiles, outfile ):
    tmpf = P.getTempFilename( "." )
    infiles = " ".join( infiles )
    headers = [ "cell_type",
                "contig_lnc",
                "start_lnc",
                "end_lnc",
                "lnc_id", 
                "score_lnc",
                "strand_lnc",
                "contig_peak",
                "start_peak",
                "end_peak",
                "peak_id",
                "score_peak",
                "strand_peak",
                "distance",
                "supporting_peak" ]
    statement = ( "python %(scriptsdir)s/combine_tables.py"
                  " --cat=cell_type"
                  " --no-titles"
                  " --log=%(outfile)s.log"
                  " %(infiles)s"
                  " > %(tmpf)s" )
    P.run()

    # hack... I can't separate the cell-type from the path
    outf = IOTools.openFile( outfile, "w" )
    outf.write( "\t".join( headers ) + "\n" )
    for line in IOTools.openFile( tmpf ):
        line = line.split()
        if int( line[-1] ) <= 2000:
            supporting_peak = "1"
        else: 
            supporting_peak = "0"
        cell_type = os.path.basename( line[0] ).split("_")[0]
        remainder = "\t".join( line[1:] )
        outf.write( cell_type + "\t" + remainder + "\t" + supporting_peak + "\n" )
    outf.close()
    os.unlink( tmpf )


@transform( stackTSSChIPPeakDistance, 
            regex( "(.+)/(.+)/(.+).tsv.gz" ),
            r"./\3.load" )
def loadTSSChIPPeakDistance( infile, outfile ):
    P.load( infile, outfile )


#################################################################################
## subsection: discover eRNAs/pRNAs
#################################################################################
@follows( loadTSSChIPPeakDistance,
          mkdir( "characterize_pRNA" ))                                    # pRNA
@split( loadK4TSSCoverage_stacked, 
            regex( "(.+)/tss_coverage_gene(.*?)_(cage_)?stacked.load" ), 
            add_inputs( r"\1/K4me3_lncRNA_gene\2_peak_distance.load" ),       # me3
            r"./characterize_pRNA/*gene\2_pRNA.bed.gz" )                   # pRNA
def findGenepRNAs( infiles, outfiles ):
    coverage_table = P.snip( os.path.basename( infiles[0] ), ".load" )
    peak_table = P.snip( os.path.basename( infiles[1] ), ".load" )
    threshold = 1.0/float( PARAMS["eRNA_fold_change"] )                    # divide
    outdir = "./characterize_pRNA"                                         # pRNA

    # specify outfile suffix
    if coverage_table.endswith( "gene_stacked" ):
        outf_suffix = "gene_pRNA.bed"                                      # pRNA
    elif coverage_table.endswith( "gene_permissive_cage_stacked" ):
        outf_suffix = "gene_permissive_pRNA.bed"                           # pRNA
    elif coverage_table.endswith( "gene_robust_cage_stacked" ):
        outf_suffix = "gene_robust_pRNA.bed"                               # pRNA
    else:
        raise Exception( "Unrecognised coverage table: %s" % coverage_table )
   
    # select cell types 
    statement = ( "SELECT DISTINCT cell_type FROM %(coverage_table)s" % locals() )
    cell_types = PU.fetch( statement )
    cell_types = [ str( x[0] ) for x in cell_types ]  

    # select pRNAs for each cell type
    for cell_type in cell_types:
        outfile = os.path.join( outdir, "_".join( [ cell_type, outf_suffix ] ) )
        statement = ( "SELECT b.contig_lnc,b.start_lnc,b.end_lnc,a.lnc_id,b.score_lnc,b.strand_lnc"
                      " FROM %(coverage_table)s AS a"
                      " INNER JOIN %(peak_table)s AS b"
                      " ON a.lnc_id = b.lnc_id"
                      " WHERE a.ratio < %(threshold)s"                     # less
                      " AND a.pass = 1"
                      " AND a.cell_type = '%(cell_type)s'"
                      " AND b.supporting_peak=1"
                      " AND b.cell_type = '%(cell_type)s'" % locals() )
        # print "\n\n\n" + outfile + "\n" + statement + "\n"

        df = PU.fetch_DataFrame( statement )
        df.to_csv( outfile, sep = "\t", header = False, index = False )

        to_cluster = False
        statement = "gzip -f %(outfile)s"
        P.run()


@follows( loadTSSChIPPeakDistance,
          mkdir( "characterize_eRNA" ))                                  # eRNA
@split( loadK4TSSCoverage_stacked, 
            regex( "(.+)/tss_coverage_gene(.*?)_(cage_)?stacked.load" ), 
            add_inputs( r"\1/K4me1_lncRNA_gene\2_peak_distance.load" ),  # K4me1
            r"./characterize_eRNA/*gene\2_eRNA.bed.gz" )                 # eRNA
def findGeneeRNAs( infiles, outfiles ):
    coverage_table = P.snip( os.path.basename( infiles[0] ), ".load" )
    peak_table = P.snip( os.path.basename( infiles[1] ), ".load" )
    threshold = 1.0*float( PARAMS["eRNA_fold_change"] )                  # multiply
    outdir = "./characterize_eRNA"                                       # eRNA

    # specify outfile suffix
    if coverage_table.endswith( "gene_stacked" ):
        outf_suffix = "gene_eRNA.bed"                                    # eRNA
    elif coverage_table.endswith( "gene_permissive_cage_stacked" ):
        outf_suffix = "gene_permissive_eRNA.bed"                         # eRNA
    elif coverage_table.endswith( "gene_robust_cage_stacked" ):
        outf_suffix = "gene_robust_eRNA.bed"                             # eRNA
    else:
        raise Exception( "Unrecognised coverage table: %s" % coverage_table )
   
    # select cell types 
    statement = ( "SELECT DISTINCT cell_type FROM %(coverage_table)s" % locals() )
    cell_types = PU.fetch( statement )
    cell_types = [ str( x[0] ) for x in cell_types ]  

    # select pRNAs for each cell type
    for cell_type in cell_types:
        outfile = os.path.join( outdir, "_".join( [ cell_type, outf_suffix ] ) )
        statement = ( "SELECT b.contig_lnc,b.start_lnc,b.end_lnc,a.lnc_id,b.score_lnc,b.strand_lnc"
                      " FROM %(coverage_table)s AS a"
                      " INNER JOIN %(peak_table)s AS b"
                      " ON a.lnc_id = b.lnc_id"
                      " WHERE a.ratio > %(threshold)s"                   # greater
                      " AND a.pass = 1"
                      " AND a.cell_type = '%(cell_type)s'"
                      " AND b.supporting_peak=1"
                      " AND b.cell_type = '%(cell_type)s'" % locals() )

        df = PU.fetch_DataFrame( statement )
        df.to_csv( outfile, sep = "\t", header = False, index = False )

        to_cluster = False
        statement = "gzip -f %(outfile)s"
        P.run()


@follows( loadTSSChIPPeakDistance,
          mkdir( "characterize_pRNA" ))                                    # pRNA
@split( loadK4TSSCoverage_stacked, 
            regex( "(.+)/tss_coverage_gene(.*?)_(cage_)?stacked.load" ), 
            add_inputs( r"\1/K4me3_lncRNA_gene\2_peak_distance.load" ),    # me3
            r"./characterize_pRNA/*gene\2_pRNA_ratio.bed.gz" )             # pRNA
def findGenepRNAs_ratio( infiles, outfiles ):
    """
    Exactly the same as above, but outputs the H3K4me1/H3K4me3 ratio as the
    score field in the resulting bedfile. 
    """
    coverage_table = P.snip( os.path.basename( infiles[0] ), ".load" )
    peak_table = P.snip( os.path.basename( infiles[1] ), ".load" )
    threshold = 1.0/float( PARAMS["eRNA_fold_change"] )                    # divide
    outdir = "./characterize_pRNA"                                         # pRNA

    # specify outfile suffix
    if coverage_table.endswith( "gene_stacked" ):
        outf_suffix = "gene_pRNA_ratio.bed"                                # pRNA_ratio
    elif coverage_table.endswith( "gene_permissive_cage_stacked" ):
        outf_suffix = "gene_permissive_pRNA_ratio.bed"                     # pRNA_ratio
    elif coverage_table.endswith( "gene_robust_cage_stacked" ):
        outf_suffix = "gene_robust_pRNA_ratio.bed"                         # pRNA_ratio
    else:
        raise Exception( "Unrecognised coverage table: %s" % coverage_table )
   
    # select cell types 
    statement = ( "SELECT DISTINCT cell_type FROM %(coverage_table)s" % locals() )
    cell_types = PU.fetch( statement )
    cell_types = [ str( x[0] ) for x in cell_types ]  

    # select pRNAs for each cell type
    for cell_type in cell_types:
        outfile = os.path.join( outdir, "_".join( [ cell_type, outf_suffix ] ) )
        statement = ( "SELECT b.contig_lnc,b.start_lnc,b.end_lnc,a.lnc_id,a.ratio,b.strand_lnc"
                      " FROM %(coverage_table)s AS a"
                      " INNER JOIN %(peak_table)s AS b"
                      " ON a.lnc_id = b.lnc_id"
                      " WHERE a.ratio < %(threshold)s"                     # less
                      " AND a.pass = 1"
                      " AND a.cell_type = '%(cell_type)s'"
                      " AND b.supporting_peak=1"
                      " AND b.cell_type = '%(cell_type)s'" % locals() )
        # print "\n\n\n" + outfile + "\n" + statement + "\n"

        df = PU.fetch_DataFrame( statement )
        df.to_csv( outfile, sep = "\t", header = False, index = False )

        to_cluster = False
        statement = "gzip -f %(outfile)s"
        P.run()


@follows( loadTSSChIPPeakDistance,
          mkdir( "characterize_eRNA" ))                                  # eRNA
@split( loadK4TSSCoverage_stacked, 
            regex( "(.+)/tss_coverage_gene(.*?)_(cage_)?stacked.load" ), 
            add_inputs( r"\1/K4me1_lncRNA_gene\2_peak_distance.load" ),  # K4me1
            r"./characterize_eRNA/*gene\2_eRNA_ratio.bed.gz" )           # eRNA
def findGeneeRNAs_ratio( infiles, outfiles ):
    """
    Exactly the same as above, but outputs the H3K4me1/H3K4me3 ratio as the
    score field in the resulting bedfile. 
    """
    coverage_table = P.snip( os.path.basename( infiles[0] ), ".load" )
    peak_table = P.snip( os.path.basename( infiles[1] ), ".load" )
    threshold = 1.0*float( PARAMS["eRNA_fold_change"] )                  # multiply
    outdir = "./characterize_eRNA"                                       # eRNA

    # specify outfile suffix
    if coverage_table.endswith( "gene_stacked" ):
        outf_suffix = "gene_eRNA_ratio.bed"                              # eRNA_ratio
    elif coverage_table.endswith( "gene_permissive_cage_stacked" ):
        outf_suffix = "gene_permissive_eRNA_ratio.bed"                   # eRNA_ratio
    elif coverage_table.endswith( "gene_robust_cage_stacked" ):
        outf_suffix = "gene_robust_eRNA_ratio.bed"                       # eRNA_ratio
    else:
        raise Exception( "Unrecognised coverage table: %s" % coverage_table )
   
    # select cell types 
    statement = ( "SELECT DISTINCT cell_type FROM %(coverage_table)s" % locals() )
    cell_types = PU.fetch( statement )
    cell_types = [ str( x[0] ) for x in cell_types ]  

    # select pRNAs for each cell type
    for cell_type in cell_types:
        outfile = os.path.join( outdir, "_".join( [ cell_type, outf_suffix ] ) )
        statement = ( "SELECT b.contig_lnc,b.start_lnc,b.end_lnc,a.lnc_id,a.ratio,b.strand_lnc"
                      " FROM %(coverage_table)s AS a"
                      " INNER JOIN %(peak_table)s AS b"
                      " ON a.lnc_id = b.lnc_id"
                      " WHERE a.ratio > %(threshold)s"                   # greater
                      " AND a.pass = 1"
                      " AND a.cell_type = '%(cell_type)s'"
                      " AND b.supporting_peak=1"
                      " AND b.cell_type = '%(cell_type)s'" % locals() )

        df = PU.fetch_DataFrame( statement )
        df.to_csv( outfile, sep = "\t", header = False, index = False )

        to_cluster = False
        statement = "gzip -f %(outfile)s"
        P.run()



# @follows( loadTSSChIPPeakDistance,
#           mkdir( "characterize_pRNA" ))                                    # pRNA
# @split( loadHighestTranscriptK4TSSCoverage, 
#             regex( "(.+)/tss_coverage_transcript(.*?)_(cage_)?stacked_highest.load" ), 
#             add_inputs( r"\1/K4_lncRNA_transcript\2_peak_distance.load" ),       # me3
#             r"./characterize_pRNA/*transcript\2_pRNA.bed.gz" )                   # pRNA
# def findTranscriptpRNAs( infiles, outfiles ):
#     coverage_table = P.snip( os.path.basename( infiles[0] ), ".load" )
#     peak_table = P.snip( os.path.basename( infiles[1] ), ".load" )
#     threshold = 1.0/float( PARAMS["eRNA_fold_change"] )                    # divide
#     outdir = "./characterize_pRNA"                                         # pRNA

#     # specify outfile suffix
#     if coverage_table.endswith( "transcript_stacked_highest" ):
#         outf_suffix = "transcript_pRNA.bed"                                      # pRNA
#     elif coverage_table.endswith( "transcript_permissive_cage_stacked_highest" ):
#         outf_suffix = "transcript_permissive_pRNA.bed"                           # pRNA
#     elif coverage_table.endswith( "transcript_robust_cage_stacked_highest" ):
#         outf_suffix = "transcript_robust_pRNA.bed"                               # pRNA
#     else:
#         raise Exception( "Unrecognised coverage table: %s" % coverage_table )
   
#     # select cell types 
#     statement = ( "SELECT DISTINCT cell_type FROM %(coverage_table)s" % locals() )
#     cell_types = PU.fetch( statement )
#     cell_types = [ str( x[0] ) for x in cell_types ]  

#     # select pRNAs for each cell type
#     for cell_type in cell_types:
#         outfile = os.path.join( outdir, "_".join( [ cell_type, outf_suffix ] ) )
#         statement = ( "SELECT b.contig_lnc,b.start_lnc,b.end_lnc,a.gene_id,b.score_lnc,b.strand_lnc"
#                       " FROM %(coverage_table)s AS a"
#                       " INNER JOIN %(peak_table)s AS b"
#                       " ON a.lnc_id = b.lnc_id"
#                       " WHERE a.ratio < %(threshold)s"                     # less
#                       " AND a.pass = 1"
#                       " AND a.cell_type = '%(cell_type)s'"
#                       " AND b.supporting_peak=1"
#                       " AND b.cell_type = '%(cell_type)s'" % locals() )
#         # print "\n\n\n" + outfile + "\n" + statement + "\n"

#         df = PU.fetch_DataFrame( statement )
#         df.to_csv( outfile, sep = "\t", header = False, index = False )

#         to_cluster = False
#         statement = "gzip -f %(outfile)s"
#         P.run()


# @follows( loadTSSChIPPeakDistance,
#           mkdir( "characterize_eRNA" ))                                    # eRNA
# @split( loadHighestTranscriptK4TSSCoverage, 
#             regex( "(.+)/tss_coverage_transcript(.*?)_(cage_)?stacked_highest.load" ), 
#             add_inputs( r"\1/K4me1_lncRNA_transcript\2_peak_distance.load" ),       # me1
#             r"./characterize_eRNA/*transcript\2_eRNA.bed.gz" )                   # eRNA
# def findTranscripteRNAs( infiles, outfiles ):
#     coverage_table = P.snip( os.path.basename( infiles[0] ), ".load" )
#     peak_table = P.snip( os.path.basename( infiles[1] ), ".load" )
#     threshold = 1.0*float( PARAMS["eRNA_fold_change"] )                    # multiply
#     outdir = "./characterize_eRNA"                                         # eRNA

#     # specify outfile suffix
#     if coverage_table.endswith( "transcript_stacked_highest" ):
#         outf_suffix = "transcript_eRNA.bed"                                      # eRNA
#     elif coverage_table.endswith( "transcript_permissive_cage_stacked_highest" ):
#         outf_suffix = "transcript_permissive_eRNA.bed"                           # eRNA
#     elif coverage_table.endswith( "transcript_robust_cage_stacked_highest" ):
#         outf_suffix = "transcript_robust_eRNA.bed"                               # eRNA
#     else:
#         raise Exception( "Unrecognised coverage table: %s" % coverage_table )
   
#     # select cell types 
#     statement = ( "SELECT DISTINCT cell_type FROM %(coverage_table)s" % locals() )
#     cell_types = PU.fetch( statement )
#     cell_types = [ str( x[0] ) for x in cell_types ]  

#     # select pRNAs for each cell type
#     for cell_type in cell_types:
#         outfile = os.path.join( outdir, "_".join( [ cell_type, outf_suffix ] ) )
#         statement = ( "SELECT b.contig_lnc,b.start_lnc,b.end_lnc,a.gene_id,b.score_lnc,b.strand_lnc"
#                       " FROM %(coverage_table)s AS a"
#                       " INNER JOIN %(peak_table)s AS b"
#                       " ON a.lnc_id = b.lnc_id"
#                       " WHERE a.ratio > %(threshold)s"                     # greater
#                       " AND a.pass = 1"
#                       " AND a.cell_type = '%(cell_type)s'"
#                       " AND b.supporting_peak=1"
#                       " AND b.cell_type = '%(cell_type)s'" % locals() )
#         # print "\n\n\n" + outfile + "\n" + statement + "\n"

#         df = PU.fetch_DataFrame( statement )
#         df.to_csv( outfile, sep = "\t", header = False, index = False )

#         to_cluster = False
#         statement = "gzip -f %(outfile)s"
#         P.run()


@collate( [ findGenepRNAs,
            findGeneeRNAs ],
#            findTranscriptpRNAs,
#            findTranscripteRNAs ], 
          regex( "(.+)/(.+?)_(gene|transcript)(.*?)_(eRNA|pRNA).bed.gz" ), 
          r"\1/\3\4_\5.tsv.gz" )
def stackeRNAs( infiles, outfile ):
    outf = IOTools.openFile( outfile, "w" )
    outf.write( "gene_set\tgene_id\n" )
    # iterate through bedfiles write cell_type and gene_id to outfile
    for bedfile in infiles:
        cell_type = os.path.basename( bedfile ).split("_")[0]
        for line in IOTools.openFile( bedfile ):
            gene_id = line.split()[3]
            outf.write( cell_type + "\t" + gene_id + "\n" )
    outf.close()


@transform( stackeRNAs, 
            regex( "(.+)/(.+).tsv.gz" ),
            r"./\2.load" )
def loadStackedeRNAs( infile, outfile ):
    P.load( infile, outfile )


@collate( [ findGenepRNAs,
            findGeneeRNAs ],
#            findTranscriptpRNAs,
#            findTranscripteRNAs ], 
          regex( "(.+)/(.+?)_(gene|transcript)(.*?)_(eRNA|pRNA).bed.gz" ), 
          r"\1/\3\4_\5_collapsed.bed.gz" )
def collapseeRNAs( infiles, outfile ):
    infiles = " ".join( infiles )
    statement = ( "zcat %(infiles)s |"
                  " sort -u |"
                  " gzip > %(outfile)s" )
    P.run()

@follows( mkdir( "characterize_lncrna" ) )
@collate( collapseeRNAs, 
          regex( "(?:.+)/(.*)_collapsed.bed.gz" ),
          add_inputs( classifyMergedLncRNAs ), 
          r"./characterize_lncrna/lncRNA_chromatin_classification.tsv.gz" )
def summarizeCollapsedeRNAs( infiles, outfile ):
    """
    Produce a table that summarizes the chromatin status of all lncRNAs, based
    on all transcript/gene datasets.
    """
    # separate infiles from reference gtf containing all lncRNAs
    lncRNA_gtf = infiles[0][1]
    infile_list = [ x[0] for x in infiles ]

    P10.summarize_eRNAs( lncRNA_gtf, infile_list, outfile )


@transform( summarizeCollapsedeRNAs,
            suffix( ".tsv.gz" ),
            add_inputs( os.path.join( PARAMS[ "location_transcriptfiles" ],
                                      "refcoding.gtf.gz" ) ),
            "_refcoding.tsv.gz" )
def collateCollapsedeRNAAndRefcoding( infiles, outfile ):
    """
    Append refcoding genes to the flat file summarizing eRNAs
    """
    eRNAs, refcoding = infiles

    # iterate through refcoding and pull out all gene_ids
    refcoding_ids = []
    for gtf in GTF.flat_gene_iterator( GTF.iterator( IOTools.openFile( refcoding ) ) ):
        refcoding_ids.append( gtf[0].gene_id )

    # get number of fields in eRNA summary file
    headers = IOTools.openFile( eRNAs ).readline()
    n_fields = len( headers.split()[1:] )

    # write outfile
    outf = IOTools.openFile( outfile, "w" )
    for line in IOTools.openFile( eRNAs ).readlines():
        outf.write( line )
    
    for refcoding_id in refcoding_ids:
        line_out = [ refcoding_id, ] + ["protein_coding"]*n_fields
        outf.write( "\t".join( line_out ) + "\n" )

    outf.close()


@transform( [ summarizeCollapsedeRNAs,
            collateCollapsedeRNAAndRefcoding ],
            regex( "(?:.+)/(.+).tsv.gz" ),
            r"./\1.load" )
def loadCollapsedeRNASummary( infile, outfile ):
    P.load( infile, outfile, options="--add-index=gene_id" )


@follows( collateCollapsedeRNAAndRefcoding,
          stackeRNAs )
def findeRNAsAndpRNAs():
    pass

#################################################################################
## subsection: discover eRNAs/pRNAs based on differential expression
#################################################################################
# do this for the robust_gene_eRNA ONLY
# for each cell stage in Pro, Pre, Immature, Mature, Follicular, Marginal, B1a, Germinal
# get DESeq results from table and write to outfile
# concatenate outfiles. 
# recalculate q values
@follows( mkdir( "characterize_lncrna_deseq" ), 
          mkdir( "characterize_lncrna_deseq/ranked_fold_change" ) )
@split( PARAMS["eRNA_deseq_db"],
        "characterize_lncrna_deseq/ranked_fold_change/*_results.tsv.gz" )
def getDESeqChromatinRatios( infile, outfiles ):
    """
    Get the log2 fold change estimates from the DESeq for each cell type
    """
    outdir = "characterize_lncrna_deseq/ranked_fold_change"
    for cell_type in ["Pro", 
                      "Pre", 
                      "Immature", 
                      "Mature", 
                      "Follicular", 
                      "Marginal", 
                      "B1a", 
                      "Germinal"]:
        table_name = "design" + cell_type + "_vs_lncTSS_feature_counts_deseq_gene_diff"
        outfile = IOTools.openFile( os.path.join( outdir, cell_type + "_results.tsv.gz" ), "w" )

        statement = ("SELECT test_id AS 'gene_id', control_name, treatment_name, l2fold, pvalue"
                     " FROM %(table_name)s" % locals() )
        df = PU.fetch_DataFrame( statement, database = infile )
        df.set_index("gene_id", inplace=True)
        
        # check treatment and control
        assert len(set(df["control_name"])) == 1
        assert df["control_name"][0] == "K4me1"
        assert len(set(df["treatment_name"])) == 1
        assert df["treatment_name"][0] == "K4me3"

        # output only the l2foldchange
        df.to_csv(outfile, columns=["l2fold",], sep ="\t")


@collate( getDESeqChromatinRatios, 
          regex( "(.+)/ranked_fold_change/(.+)_results.tsv.gz" ),
          r"\1/tss_coverage_gene_robust_cage_deseq_l2foldchange.tsv.gz" )
def combineDESeqChromatinRatios( infiles, outfile ):
    infiles = " ".join( infiles )
    statement = ("python %(scriptsdir)s/combine_tables.py"
                 " --columns=1"
                 " --add-file-prefix"
                 " --regex-filename='(.+)_results.tsv.gz'"
                 " --log=%(outfile)s.log"
                 " %(infiles)s |"
                 " gzip > %(outfile)s" )
    P.run()

@transform( combineDESeqChromatinRatios, 
           regex( "(.+)/(.+).tsv.gz" ),
           r"\2.load" )
def loadDESeqChromatinRatios(infile, outfile):
    P.load(infile, outfile, options="--add-index=gene_id")


@follows( mkdir( "characterize_eRNA_deseq" ) )
@split( PARAMS["eRNA_deseq_db"],
        "characterize_eRNA_deseq/*_results.tsv.gz" )
def getDESeqeRNAClassification( infile, outfiles ):
    """
    Select TSS where H3K4me3:H3K4me1 l2 fold-change is < -2
    Filtering on pvalue is done subsequently
    """
    # DESeq results (for all cell types) test me3:me1 ratio
    # Therefore low l2fold change values identify eRNAs. 

    outdir = "characterize_eRNA_deseq"
    direction = "<" # i.e. less than
    fold_change = str( -1 * float(PARAMS["eRNA_deseq_l2_fold_change"]) )
    
    for cell_type in ["Pro", 
                      "Pre", 
                      "Immature", 
                      "Mature", 
                      "Follicular", 
                      "Marginal", 
                      "B1a", 
                      "Germinal"]:
        table_name = "design" + cell_type + "_vs_lncTSS_feature_counts_deseq_gene_diff"
        outfile = IOTools.openFile( os.path.join( outdir, cell_type + "_results.tsv.gz" ), "w" )
        df = P10.getDEIntervals( infile, table_name, direction, fold_change)
        df.to_csv( outfile, sep = "\t", index=False )


@follows( mkdir( "characterize_pRNA_deseq" ), getDESeqeRNAClassification ) # issues with s3 db access
@split( PARAMS["eRNA_deseq_db"],
        "characterize_pRNA_deseq/*_results.tsv.gz" )
def getDESeqpRNAClassification( infile, outfiles ):
    """
    Select TSS where H3K4me3:H3K4me1 l2 fold-change is > 2
    Filtering on pvalue is done subsequently
    """
    # DESeq results (for all cell types) test me3:me1 ratio
    # Therefore low l2fold change values identify eRNAs. 

    outdir = "characterize_pRNA_deseq"
    direction = ">" # i.e. greater than
    fold_change = PARAMS["eRNA_deseq_l2_fold_change"]

    
    for cell_type in ["Pro", 
                      "Pre", 
                      "Immature", 
                      "Mature", 
                      "Follicular", 
                      "Marginal", 
                      "B1a", 
                      "Germinal"]:
        table_name = "design" + cell_type + "_vs_lncTSS_feature_counts_deseq_gene_diff"
        outfile = IOTools.openFile( os.path.join( outdir, cell_type + "_results.tsv.gz" ), "w" )
        df = P10.getDEIntervals( infile, table_name, direction, fold_change)
        df.to_csv( outfile, sep = "\t", index=False )
    

@follows( mkdir("characterize_lncrna_deseq") )
@merge( [getDESeqeRNAClassification, getDESeqpRNAClassification],
        "characterize_lncrna_deseq/lncRNA_TSS_chromatin_state_deseq.tsv.gz" )        
def collateDESeqlncRNAClassification(infiles, outfile):
    """
    Concatenate all tables with significant TSS intervals for each cell type. 
    """
    tmpf = P.getTempFilename("/ifs/scratch")
    outf = P.snip(outfile, ".gz")
    tables = " ".join( infiles )
    statement = ( "python %(scriptsdir)s/combine_tables.py"
                  " --cat=CAT"
                  " --log=%(outfile)s.log"
                  " %(tables)s"
                  " > %(tmpf)s" )
    P.run()


    # Add an adjusted P value to table
    # clear global namespace
    R('''rm(list=ls())''')
    R('''df <- read.table("%(tmpf)s", header=TRUE, sep="\t", stringsAsFactors=FALSE)''' % locals())
    R('''df <- cbind(df, padj = p.adjust(df$pvalue, method="BH"))''')
    R('''write.table(df, file="%(outf)s", quote=FALSE, row.names=FALSE, sep="\t")''' % locals())
    
    to_cluster = False
    statement = "gzip -f %(outf)s"
    P.run()
    
    os.unlink(tmpf)


@transform( collateDESeqlncRNAClassification, 
            regex( "(.+)/(.+).tsv.gz" ),
            r"\2.load" )
def loadDESeqlncRNAClassification( infile, outfile ):
    tempfile = P.getTempFilename("/ifs/scratch")
    tmpf = IOTools.openFile( tempfile, "w" )
    tmpf.write("classification\tcell_type\ttest_id\tcontrol_name\t"
               "treatment_name\tl2fold\tpvalue\tpadj\n")
    for line in IOTools.openFile( infile ):
        if line.startswith("CAT"): continue
        line = line.split()
        chromatin_class, cell_type  = line[0].split("/")
        chromatin_class = chromatin_class.split("_")[1]
        cell_type = cell_type.split("_")[0]
        line_out = [ chromatin_class, cell_type ]
        line_out.extend( line[1:] )
        tmpf.write( "\t".join( line_out ) + "\n" )
    tmpf.close()

    P.load( tempfile, outfile )


@subdivide( resetLncRNATSS_robust,
            regex( "(.+)/lncRNA_gene_robust_cage_tss.bed.gz" ),
            add_inputs( collateDESeqlncRNAClassification ),
            "characterize_eRNA_deseq/*.bed.gz" )
def getDESeqeRNABedfiles( infiles, outfiles ):
    """
    Take BH adjusted pvalues and output a bed file of intervals identified as eRNAs
    for each cell type
    """
    bed_file, deseq_file = infiles
    eRNAs_passed = {"Pro": [], "Pre": [], "Immature": [], "Mature":[], 
                    "Follicular": [], "Marginal": [], "B1a": [], "Germinal": []}
    eRNAs_failed = {"Pro": 0, "Pre": 0, "Immature": 0, "Mature": 0, 
                    "Follicular": 0, "Marginal": 0, "B1a": 0, "Germinal": 0 }

    # Iterate through the padjusted DESeq results 
    # and out put intervals < threhold qvalue to dictionary where key is cell_type
    for line in IOTools.openFile( deseq_file ):
        if line.startswith("CAT"): continue
        line = line.split()
        chromatin_class, cell_type  = line[0].split("/")
        chromatin_class = chromatin_class.split("_")[1]
        cell_type = cell_type.split("_")[0]
        assert cell_type in eRNAs_passed.keys()

        if chromatin_class == "eRNA":
            assert float(line[4]) < 0, "eRNAs should have -ve fold change"
            if float(line[6]) <= float(PARAMS["eRNA_deseq_qvalue"]):
                eRNAs_passed[cell_type].append( line[1] )
            else:
                eRNAs_failed[cell_type] += 1
        else:
            assert chromatin_class == "pRNA", "Check naming in infile"
    E.info( "Intervals that failed qvalue threshold (eRNA):\n %s" % str(eRNAs_failed) )
        
    # load bedfile into dict
    bed_dict = {}
    for bed in Bed.iterator(IOTools.openFile( bed_file )):
        gene_id = bed.fields[0]
        bed_dict[gene_id] = str(bed)

    # iterate through the eRNA dictionary and write bedfile entries to cell type-
    # specific out file
    for cell_type, intervals in eRNAs_passed.iteritems():
        outfile = os.path.join( "characterize_eRNA_deseq", cell_type + "_eRNA.bed.gz" )
        outf = IOTools.openFile( outfile, "w" )
        for interval in intervals:
            outf.write( bed_dict[ interval ] + "\n" )
        outf.close()
        

@subdivide( resetLncRNATSS_robust,
            regex( "(.+)/lncRNA_gene_robust_cage_tss.bed.gz" ),
            add_inputs( collateDESeqlncRNAClassification ),
            "characterize_pRNA_deseq/*.bed.gz" )
def getDESeqpRNABedfiles( infiles, outfiles ):
    """
    Take BH adjusted pvalues and output a bed file of intervals identified as eRNAs
    for each cell type
    """
    bed_file, deseq_file = infiles
    pRNAs_passed = {"Pro": [], "Pre": [], "Immature": [], "Mature":[], 
                    "Follicular": [], "Marginal": [], "B1a": [], "Germinal": []}
    pRNAs_failed = {"Pro": 0, "Pre": 0, "Immature": 0, "Mature": 0, 
                    "Follicular": 0, "Marginal": 0, "B1a": 0, "Germinal": 0 }

    # Iterate through the padjusted DESeq results 
    # and out put intervals < threhold qvalue to dictionary where key is cell_type
    for line in IOTools.openFile( deseq_file ):
        if line.startswith("CAT"): continue
        line = line.split()
        chromatin_class, cell_type  = line[0].split("/")
        chromatin_class = chromatin_class.split("_")[1]
        cell_type = cell_type.split("_")[0]
        assert cell_type in pRNAs_passed.keys()

        if chromatin_class == "pRNA":
            assert float(line[4]) > 0, "pRNAs should have +ve fold change"
            if float(line[6]) <= float(PARAMS["eRNA_deseq_qvalue"]):
                pRNAs_passed[cell_type].append( line[1] )
            else:
                pRNAs_failed[cell_type] += 1
        else:
            assert chromatin_class == "eRNA", "Check naming in infile"
    E.info( "Intervals that failed qvalue threshold (pRNA):\n %s" % str(pRNAs_failed) )
        
    # load bedfile into dict
    bed_dict = {}
    for bed in Bed.iterator(IOTools.openFile( bed_file )):
        gene_id = bed.fields[0]
        bed_dict[gene_id] = str(bed)

    # iterate through the eRNA dictionary and write bedfile entries to cell type-
    # specific out file
    for cell_type, intervals in pRNAs_passed.iteritems():
        outfile = os.path.join( "characterize_pRNA_deseq", cell_type + "_pRNA.bed.gz" )
        outf = IOTools.openFile( outfile, "w" )
        for interval in intervals:
            outf.write( bed_dict[ interval ] + "\n" )
        outf.close()



@collate( [ getDESeqpRNABedfiles, getDESeqeRNABedfiles ],
          regex( "(.+)/(.+?)_(eRNA|pRNA).bed.gz" ), 
          r"\1/gene_robust_\3.tsv.gz" )
def stackeRNAs_deseq( infiles, outfile ):
    outf = IOTools.openFile( outfile, "w" )
    outf.write( "gene_set\tgene_id\n" )
    # iterate through bedfiles write cell_type and gene_id to outfile
    for bedfile in infiles:
        cell_type = os.path.basename( bedfile ).split("_")[0]
        for line in IOTools.openFile( bedfile ):
            gene_id = line.split()[3]
            outf.write( cell_type + "\t" + gene_id + "\n" )
    outf.close()


@collate( [ getDESeqpRNABedfiles, getDESeqeRNABedfiles ],
          regex( "(.+)/(.+?)_(eRNA|pRNA).bed.gz" ), 
          r"\1/gene_robust_\3_collapsed.bed.gz" )
def collapseeRNAs_deseq( infiles, outfile ):
    infiles = " ".join( infiles )
    statement = ( "zcat %(infiles)s |"
                  " sort -u |"
                  " gzip > %(outfile)s" )
    P.run()


@follows( mkdir( "characterize_lncrna_deseq" ) )
@collate( collapseeRNAs_deseq, 
          regex( "(?:.+)/(.*)_collapsed.bed.gz" ),
          add_inputs( classifyMergedLncRNAs ), 
          r"./characterize_lncrna_deseq/lncRNA_chromatin_classification.tsv.gz" )
def summarizeCollapsedeRNAs_deseq( infiles, outfile ):
    """
    Produce a table that summarizes the chromatin status of all lncRNAs, based
    on all transcript/gene datasets.
    """
    # separate infiles from reference gtf containing all lncRNAs
    lncRNA_gtf = infiles[0][1]
    infile_list = [ x[0] for x in infiles ]

    P10.summarize_eRNAs( lncRNA_gtf, infile_list, outfile )


@transform( summarizeCollapsedeRNAs_deseq,
            suffix( ".tsv.gz" ),
            add_inputs( os.path.join( PARAMS[ "location_transcriptfiles" ],
                                      "refcoding.gtf.gz" ) ),
            "_refcoding.tsv.gz" )
def collateCollapsedeRNAAndRefcoding_deseq( infiles, outfile ):
    """
    Append refcoding genes to the flat file summarizing eRNAs
    """
    eRNAs, refcoding = infiles

    # iterate through refcoding and pull out all gene_ids
    refcoding_ids = []
    for gtf in GTF.flat_gene_iterator( GTF.iterator( IOTools.openFile( refcoding ) ) ):
        refcoding_ids.append( gtf[0].gene_id )

    # get number of fields in eRNA summary file
    headers = IOTools.openFile( eRNAs ).readline()
    n_fields = len( headers.split()[1:] )

    # write outfile
    outf = IOTools.openFile( outfile, "w" )
    for line in IOTools.openFile( eRNAs ).readlines():
        outf.write( line )
    
    for refcoding_id in refcoding_ids:
        line_out = [ refcoding_id, ] + ["protein_coding"]*n_fields
        outf.write( "\t".join( line_out ) + "\n" )

    outf.close()


@follows( collateCollapsedeRNAAndRefcoding_deseq,
          stackeRNAs_deseq )
def findDESeqeRNAAndpRNA():
    pass


#################################################################################
## subsection: collate empirical and DESeq chromatin classifications
#################################################################################
@merge( [ collateCollapsedeRNAAndRefcoding,
          collateCollapsedeRNAAndRefcoding_deseq ],
        "lncRNA_refcoding_chromatin_classification.tsv.gz" )
def collateLncRNAChromatinClassifications( infiles, outfile ):
    """
    Combine the chromatin classifications from empirical and DESeq identification
    across robust gene TSS
    """
    tmpf = P.getTempFilename( "/ifs/scratch" )
    tables = " ".join( infiles )
    statement = ( "python %(scriptsdir)s/combine_tables.py"
                  " --columns=1"
                  " --take=4"
                  " --skip-titles"
                  " --header-names=class_empirical,class_deseq"
                  " --log=%(outfile)s.log"
                  " %(tables)s"
                  " > %(tmpf)s" )
    P.run()

    outf = IOTools.openFile( outfile, "w" )
    for line in IOTools.openFile(tmpf):
        if line.startswith("bin"):
            outf.write("gene_id\tclass_empirical\tclass_deseq\tclass_consensus\n")
        else:
            line = line.split()
            if line[1] == line[2]:
                outf.write( "\t".join( [line[0], line[1], line[2], line[1]] ) + "\n" )
            else:
                outf.write( "\t".join( [line[0], line[1], line[2] , "conflicted"] ) + "\n" )
    outf.close()

    os.unlink(tmpf)

@transform( collateLncRNAChromatinClassifications, 
            suffix( ".tsv.gz" ), 
            ".load" )
def loadLncRNAChromatinClassifications( infile, outfile ):
    P.load( infile, outfile, options="--add-index=gene_id" )

#################################################################################
## subsection: plot lncRNA coverage for eRNAs
#################################################################################
@follows( mkdir("characterize_eRNA_plots") )
@subdivide( findGeneeRNAs_ratio, 
            regex( ".+/(follicular)_gene_robust_eRNA_ratio.bed.gz" ),
            add_inputs( r"/ifs/projects/proj010/analysis_chipseq011/pc_merged_bamfiles/merged_deduped_to_keep/\1-K4-R0_deduped.bwa.bam",
                        r"/ifs/projects/proj010/analysis_chipseq011/pc_merged_bamfiles/merged_deduped_to_keep/\1-K4me1-R0_deduped.bwa.bam",
                        r"input_rnaseq_bamfiles_filtered_merged/Bcell-\1-R0.bam" ),
            r"characterize_eRNA_plots/\1_plots.sentinel" )
def ploteRNACoverage(infiles, outfile):
    """
    Take the top five eRNAs - based on ratio - for each cell type and plot
    """
    bedfile, K4me3_bam, K4me1_bam, rna_bam = infiles

    # get the top 5 lncRNA TSS
    lncrna_tss = []
    for line in IOTools.openFile( bedfile ):
        lncrna_tss.append( line.split() )
    lncrna_tss = sorted(lncrna_tss, key=lambda x: x[4])
    lncrna_tss = lncrna_tss[-5:]

    intervals = []
    for i in lncrna_tss:
        gene_id, contig, start, end = i[3], i[0], int(i[1])-5000, int(i[2])+5000
        intervals.append([ gene_id, contig, start, end ])
    
    # create track dictionary
    tracks = { "H3K4me3": {"file": K4me3_bam, "type": "h"},
               "H3K4me1": {"file": K4me1_bam, "type": "h"},
               "RNA-Seq": {"file": rna_bam, "type": "h"} }
    
    # create outfile stub. 
    outf_stub = P.snip( outfile, "_plots.sentinel" )

    P10.GenomePlot( intervals, tracks, outf_stub )


#################################################################################
## subsection: calculate eRNA/pRNA consistency
#################################################################################
@follows(stackeRNAs, mkdir("characterize_eRNA_consistency"))
@split( loadK4TSSCoverage_stacked,
        regex( "(.+)/tss_coverage_gene(.*?)_(cage_)?stacked.load" ),
        add_inputs( r"\1/characterize_eRNA/gene\2_eRNA.tsv.gz" ),
        r"./characterize_eRNA_consistency/*gene\2_eRNA_gsea_summary_table.tsv" )
def calcGeneeRNAConsistency( infiles, outfile ):
    assert os.path.exists( infiles[1] )
    coverage_table, gene_sets = infiles
    coverage_table = P.snip( os.path.basename( coverage_table ), ".load" )

    # specify outfile format
    outdir = "./characterize_eRNA_consistency"
    out_suffix = P.snip( os.path.basename( gene_sets ), ".tsv.gz" ) + "_gsea_summary_table.tsv"

    P10.calculateGeneSetEnrichment( coverage_table, 
                                    gene_sets, 
                                    outdir,
                                    out_suffix )


# @follows(stackeRNAs, mkdir("characterize_eRNA_consistency"))
# @split( loadHighestTranscriptK4TSSCoverage,
#         regex( "(.+)/tss_coverage_transcript(.*?)_(cage_)?stacked_highest.load" ),
#         add_inputs( r"\1/characterize_eRNA/transcript\2_eRNA.tsv.gz" ),
#         r"./characterize_eRNA_consistency/*transcript\2_eRNA_gsea_summary_table.tsv" )
# def calcTranscripteRNAConsistency( infiles, outfile ):
#     assert os.path.exists( infiles[1] )
#     coverage_table, gene_sets = infiles
#     coverage_table = P.snip( os.path.basename( coverage_table ), ".load" )

#     # specify outfile format
#     outdir = "./characterize_eRNA_consistency"
#     out_suffix = P.snip( os.path.basename( gene_sets ), ".tsv.gz" ) + "_gsea_summary_table.tsv"

#     P10.calculateGeneSetEnrichment( coverage_table, 
#                                     gene_sets, 
#                                     outdir,
#                                     out_suffix, 
#                                     transcript = True )


@follows(stackeRNAs, mkdir("characterize_pRNA_consistency"))
@split( loadK4TSSCoverage_stacked,
        regex( "(.+)/tss_coverage_gene(.*?)_(cage_)?stacked.load" ),
        add_inputs( r"\1/characterize_pRNA/gene\2_pRNA.tsv.gz" ),
        r"./characterize_pRNA_consistency/*gene\2_pRNA_gsea_summary_table.tsv" )
def calcGenepRNAConsistency( infiles, outfile ):
    assert os.path.exists( infiles[1] )
    coverage_table, gene_sets = infiles
    coverage_table = P.snip( os.path.basename( coverage_table ), ".load" )

    # specify outfile format
    outdir = "./characterize_pRNA_consistency"
    out_suffix = P.snip( os.path.basename( gene_sets ), ".tsv.gz" ) + "_gsea_summary_table.tsv"

    P10.calculateGeneSetEnrichment( coverage_table, 
                                    gene_sets, 
                                    outdir,
                                    out_suffix )


# @follows(stackeRNAs, mkdir("characterize_pRNA_consistency"))
# @split( loadHighestTranscriptK4TSSCoverage,
#         regex( "(.+)/tss_coverage_transcript(.*?)_(cage_)?stacked_highest.load" ),
#         add_inputs( r"\1/characterize_pRNA/transcript\2_pRNA.tsv.gz" ),
#         r"./characterize_pRNA_consistency/*transcript\2_pRNA_gsea_summary_table.tsv" )
# def calcTranscriptpRNAConsistency( infiles, outfile ):
#     assert os.path.exists( infiles[1] )
#     coverage_table, gene_sets = infiles
#     coverage_table = P.snip( os.path.basename( coverage_table ), ".load" )

#     # specify outfile format
#     outdir = "./characterize_pRNA_consistency"
#     out_suffix = P.snip( os.path.basename( gene_sets ), ".tsv.gz" ) + "_gsea_summary_table.tsv"

#     P10.calculateGeneSetEnrichment( coverage_table, 
#                                     gene_sets, 
#                                     outdir,
#                                     out_suffix, 
#                                     transcript = True )          

@jobs_limit( 1, "RGlobalEnv" )
@collate( [ calcGeneeRNAConsistency,
#            calcTranscripteRNAConsistency,
            calcGenepRNAConsistency ],
#            calcTranscriptpRNAConsistency ],
          regex( "(.+)/(.+?)_(gene|transcript)(.+)_summary_table.tsv" ),
          r"\1/\3\4_matrix.tsv" )
def ploteRNApRNAConsistency( infiles, outfile ):
    """
    Takes GSEA output and creates a table showing the enrichment of genset 
    within the fold change rank of each sample.
    """
    # outfile matrix
    # outf_matrix = P.snip( outfile, ".png" ) + ".tsv"
    outf_matrix = outfile
    # specify the order for the cell_types in resulting outfile. 
    cell_types = [ "pro", 
                   "pre", 
                   "immature", 
                   "mature", 
                   "follicular", 
                   "marginal", 
                   "b1a", 
                   "germinal" ]

    # iterate through infiles and create dataframe
    first = True
    for infile in infiles:
        cell_type = os.path.basename( infile ).split("_")[0]
        # take first infile and create dataframe
        if first:
            first = False
            df = pd.read_table( infile, sep="\t", index_col=0 )
            # select only enrichment scores
            df = pd.DataFrame( { cell_type : df.iloc[:,1] } )
            continue
        # add the remaining infiles as series to existing dataframe
        inf = pd.read_table( infile, sep="\t", index_col=0 )
        inf = inf.iloc[:,1]
        df[ cell_type ] = inf


    # Hack: It turns out that at the selected threshold, it is not possible
    # pRNAs in either marginal or pre samples.
    # These data are therefore missing from the resulting dataframe - add NAs
    row_vals = map(str, df.index.tolist())
    missing = set(cell_types) - set(row_vals)
    for i in missing:
        df.loc[i] = ["NA" for x in df.columns.names]

    # sort dataframe by order specified in cell_types
    df = df.loc[cell_types, cell_types]
    # write dataframe to file. 
    df.to_csv( outf_matrix, sep="\t" )
    
    ## This is currently not implemented... the plots have yet to be made
    # # pickle_dump dataframe and pass to P.submit
    # tmpf_name = P.getTempFilename(".")
    # tmpf = open( tmpf_name, "wb" )
    # pickle.dump( df, tmpf )
    # tmpf.close()

    # # run R plot on cluster
    # P.submit( "/ifs/devel/projects/proj010/PipelineProj010",
    #           "plotConsistencyHeatmap",
    #           infiles = [ tmpf_name, ], 
    #           outfiles = outfile )
    # os.unlink( tmpf_name )


@transform( ploteRNApRNAConsistency, 
            regex( "(?:.+)/gene_robust(.+).tsv" ),
            r"gene_robust\1.load" )
def loadeRNApRNAConsistency( infile, outfile ):

    tmpf = P.getTempFile( "/ifs/scratch" )
    header = True
    for line in IOTools.openFile( infile ):
        line = line.split()
        if header:
            header = False
            line_out = ["cell_type", ] + line
            tmpf.write( "\t".join( line_out ) + "\n" )
        else:
            tmpf.write( "\t".join( line ) + "\n" )
    tmpf.close()

    P.load( tmpf.name, outfile )


#################################################################################
#################################################################################
#################################################################################
# Section: Explore eRNA/pRNA characteristics
#################################################################################
#################################################################################


#################################################################################
#################################################################################
#### METASECTION #### Weighted Gene Co-expression Network Analysis ####
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
# Section: WCGNA on protein coding gene expression
#################################################################################
#################################################################################
@follows()
def loadWGCNALibraries():
    R( '''suppressMessages(library(WGCNA))''' )
    R('''suppressMessages(library(flashClust))''')
    R( '''allowWGCNAThreads(nThreads=8)''' )
#################################################################################
## subsection: filter and clean geneset
#################################################################################
@follows( mkdir( "./expression_wgcna_filter_refcoding" ) )
@transform( [ extractPerSampleFPKMs, 
              logTransformCuffdiffFPKMs ],
            regex( "(.+)/lncRNA_refcoding_cuffdiff_(.*)fpkms.tsv.gz" ),
            r"./expression_wgcna_filter_refcoding/refcoding_cuffdiff_\2wgcna_fpkms.tsv.gz" )
def extractRefcodingFPKMFile( infile, outfile ):
    statement = ( " zcat %(infile)s |"
                  " head -1 |"
                  " gzip"
                  " > %(outfile)s;"
                  " zcat %(infile)s |"
                  " grep ENSMUSG |"
                  " gzip"
                  " >> %(outfile)s" )
    P.run()


@transform( extractRefcodingFPKMFile, 
            suffix( "_fpkms.tsv.gz" ), 
            "_fpkms_filtered_1.tsv.gz" )
def maskOutliersRefcodingFPKMFile( infile, outfile ):
    out_filtered = P.snip( outfile, ".gz" )
    out_masked = re.sub( "filtered", "masked", out_filtered )

    # number of std. dev. from mean at which sample values are masked
    n_std_dev = PARAMS["wgcna_mask_stdev"]

    # open infile as pandas dataframe, with samples as columns and genes as rows
    df = pd.read_table( infile, 
                        compression = "gzip", 
                        index_col = 0, 
                        header = 0, 
                        na_filter = True, 
                        na_values = ["NULL", "NA", "nan", "NaN"] )

    # write outfile containing those genes with masked sample values
    P10.record_outliers( df, out_masked, n_std_dev )

    # mask values that fall greater than n_std_dev from mean
    # axis = 1 means apply across rows
    df = df.apply( P10.mask_outliers, 
                   axis = 1, 
                   raw = True, 
                   std_thresh = n_std_dev )
    
    # write outfiles
    df.to_csv( out_filtered, sep = "\t", na_rep = "NaN" )

    # compress outfiles
    to_cluster = False
    statement = ( "gzip %(out_filtered)s;"
                  " gzip %(out_masked)s" )
    P.run()


@transform( maskOutliersRefcodingFPKMFile, 
            suffix( "_fpkms_filtered_1.tsv.gz" ), 
            "_fpkms_filtered_2.tsv.gz" )
def dropNARefcodingFPKMFile( infile, outfile ):
    out_filtered = P.snip( outfile, ".gz" )
    out_rejected = re.sub( "filtered", "rejected", out_filtered )

    # open infile as pandas dataframe
    df = pd.read_table( infile, 
                        compression = "gzip", 
                        index_col = 0, 
                        header = 0, 
                        na_filter = True, 
                        na_values = ["NULL", "NA", "nan", "NaN"] )

    # minimum number/pptn of samples that must have data
    min_data_limit = P10.set_threshold( PARAMS["wgcna_filter_nan"], 
                                        len( df.columns ) )

    E.info( "Removing genes with too many missing values:\n"
            "\t\tThere are %s samples within dataframe\n"
            "\t\tThe minimum data threshold is set at %s\n"
            "\t\tGenes will therefore be removed if they have data for"
            " fewer than %s samples" % ( str(len(df.columns)), 
                                        PARAMS["wgcna_filter_nan"],
                                        str(min_data_limit) ) )

    # create df of genes that do not meet criteria to outfile_reject
    df_rej = df[ df.apply( P10.drop_nan, 
                       axis = 1, 
                       data_thresh = min_data_limit, 
                       rejected = True ) ]

    # create df of genes that do meet criteria to outfile_filter
    df = df[ df.apply( P10.drop_nan, 
                       axis = 1, 
                       data_thresh = min_data_limit, 
                       rejected = False ) ]

    # write outfiles
    df.to_csv( out_filtered, sep = "\t", na_rep = "NaN" )
    df_rej.to_csv( out_rejected, sep = "\t", na_rep = "NaN" )

    # compress outfiles
    to_cluster = False
    statement = ( "gzip %(out_filtered)s;"
                  " gzip %(out_rejected)s" )
    P.run()
    

@follows( mkdir( "expression_wgcna_filtered_fpkm" ) )
@transform( dropNARefcodingFPKMFile, 
            regex( "(?:.+)/(.+)_fpkms_filtered_2.tsv.gz" ), 
            r"./expression_wgcna_filtered_fpkm/\1_high_fpkm.tsv.gz" )
def dropLowFPKMRefcodingFPKMFile( infile, outfile ):
    out_filtered = P.snip( outfile, ".gz" )
    out_rejected = out_filtered + "_rejected"

    # open infile as pandas dataframe
    df = pd.read_table( infile, 
                        compression = "gzip", 
                        index_col = 0, 
                        header = 0, 
                        na_filter = True, 
                        na_values = ["NULL", "NA", "nan", "NaN"] )
    
    # specify i) the minimum fpkm acceptable, ii) the minimumn number of samples 
    # that must have an fpkm value greater than specified minimum.
    
    min_fpkm = PARAMS["wgcna_filter_min_fpkm"]
    min_n_samples = P10.set_threshold( PARAMS["wgcna_filter_min_samples_fpkm"], 
                                       len( df.columns ) )
    E.info( "Removing genes based on an FPKM threshold:\n"
            "\t\tFPKM threshold set at %s\n"
            "\t\tSample threshold set at %s\n"
            "\t\tThere are %s samples in the dataframe\n"
            "\t\tTherefore genes will be removed if they have fewer than %s\n"
            "\t\tsamples with an FPKM above %s" % ( min_fpkm, 
                                                    PARAMS["wgcna_filter_min_samples_fpkm"], 
                                                    str(len(df.columns)), 
                                                    min_n_samples, 
                                                    min_fpkm ) )

    # create df of genes that do not meet criteria: outfile_reject
    df_rej = df[ df.apply( P10.drop_fpkm, 
                           axis = 1, 
                           min_fpkm = min_fpkm, 
                           fpkm_thresh = min_n_samples,                           
                           rejected = True ) ]

    # create df of genes that do meet criteria: outfile_filter
    df = df[ df.apply( P10.drop_fpkm, 
                       axis = 1, 
                       min_fpkm = min_fpkm, 
                       fpkm_thresh = min_n_samples,                           
                       rejected = False ) ]

    # write outfiles
    df.to_csv( out_filtered, sep = "\t", na_rep = "NaN" )
    df_rej.to_csv( out_rejected, sep = "\t", na_rep = "NaN" )

    # compress outfiles
    to_cluster = False
    statement = ( "gzip %(out_filtered)s;"
                  " gzip %(out_rejected)s" )
    P.run()


@follows( mkdir( "expression_wgcna_filtered_high_cv" ) )
@transform( dropLowFPKMRefcodingFPKMFile, 
            regex( "(?:.+)/(.+)_high_fpkm.tsv.gz" ), 
            r"./expression_wgcna_filtered_high_cv/\1_high_cv.tsv.gz" )
def dropLowCVRefcodingFPKMFile( infile, outfile ):
    out_filtered = P.snip( outfile, ".gz" )
    out_rejected = out_filtered + "_rejected"

    # open infile as pandas dataframe,
    df = pd.read_table( infile, 
                        compression = "gzip", 
                        index_col = 0, 
                        header = 0, 
                        na_filter = True, 
                        na_values = ["NULL", "NA", "nan", "NaN"] )

    # specify the minimum coefficient of variation that is acceptable for a gene
    # to be retained (pptn).
    cv_thresh = PARAMS["wgcna_filter_min_cv"]

    E.info( "Removing genes based on their coefficient of variation:\n"
            "\t\tRemoving genes with a cv of less than %s percent"
            % str( float( PARAMS["wgcna_filter_min_cv"] ) * 100 ) )

    # create df of genes that do not meet criteria: outfile_reject
    df_rej = df[ df.apply( P10.drop_cv,
                           axis = 1, 
                           cv_thresh = cv_thresh,
                           rejected = True ) ]
    
    # create df of genes that do not meet criteria: outfile_filter
    df = df[ df.apply( P10.drop_cv, 
                       axis = 1, 
                       cv_thresh = cv_thresh,
                       rejected = False ) ]

    # write outfiles
    df.to_csv( out_filtered, sep = "\t", na_rep = "NaN" )
    df_rej.to_csv( out_rejected, sep = "\t", na_rep = "NaN" )

    # compress outfiles
    to_cluster = False
    statement = ( "gzip %(out_filtered)s;"
                  " gzip %(out_rejected)s" )
    P.run()

    ##Extra
    # write the cv of the rejected samples
    out_r_cv = re.sub( "rejected", "rejected_CV", out_rejected )
    out_f_cv = re.sub( ".tsv", "filtered_CV.tsv", out_filtered )
    series_r_cv = df_rej.apply( P10.calc_cv, axis = 1 )
    series_f_cv = df.apply( P10.calc_cv, axis = 1 )
    series_r_cv.to_csv( out_r_cv, sep = "\t" )
    series_f_cv.to_csv( out_f_cv, sep = "\t" )


# @transform( dropLowFPKMRefcodingFPKMFile, 
#             suffix( "_fpkms_filtered_3.tsv.gz" ), 
#             "_fpkms_filtered_logged.tsv.gz" )
# def log2Transform( infile, outfile ):
#     outf = P.snip( outfile, ".gz" )

#     # open infile as pandas dataframe,
#     df = pd.read_table( infile, 
#                         compression = "gzip", 
#                         index_col = 0, 
#                         header = 0, 
#                         na_filter = True, 
#                         na_values = ["NULL", "NA", "nan", "NaN"] )

#     # log2 transform
#     df = df.apply( lambda x: np.log2( x + 0.0001 ), axis = 1 )

#     # write to outfile
#     df.to_csv( outf, sep = "\t", na_rep = "NaN" )

#     # compress outfiles
#     to_cluster = False
#     statement = ( "gzip %(outf)s;" )
                  
#     P.run()

@follows( dropLowCVRefcodingFPKMFile )
def filterRefcodingFPKMs_wgcna():
    pass

#################################################################################
## subsection: select rlog transformed data for wgcna analysis
#################################################################################
@follows( mkdir( "expression_wgcna_filtered_rlog" ), dropLowFPKMRefcodingFPKMFile )
@transform( rlogTransformData, 
            regex( "(.+)/lncRNA_(.+)_rlog_trans_counts.tsv.gz" ),
            add_inputs( "expression_wgcna_filtered_fpkm/refcoding_cuffdiff_wgcna_high_fpkm.tsv.gz" ),
            r"expression_wgcna_filtered_rlog/\2_wgcna_rlog.tsv.gz" )
def findRLogTransformedDataForWGCNA( infiles, outfile ):
    """
    Extract Rlog transformed count data for genes that passed the FPKM filtering
    thresholds for WGCNA.
    """
    rlog_file, fpkm_file = infiles

    # retrieve gene_ids for genes that passed fpkm filtering
    good_gene_ids = []
    for line in IOTools.openFile(fpkm_file):
        if line.startswith("gene_id"): continue
        good_gene_ids.append( line.split()[0] )

    outf = IOTools.openFile(outfile, "w")
    for line in IOTools.openFile(rlog_file):
        if line.startswith("gene_id"):
            outf.write( line )
        elif line.split()[0] in good_gene_ids:
            outf.write( line )
        else:
            continue
    outf.close()            


#################################################################################
## subsection: perform WGCNA on consistently highly expressed genes
#################################################################################
#  takes the filtered geneset following dropLowFPKMRefcodingFPKMFile
@jobs_limit( 1, "wgcna" )
#@follows( loadWGCNALibraries )
@transform( [ dropLowFPKMRefcodingFPKMFile, 
              dropLowCVRefcodingFPKMFile, 
              findRLogTransformedDataForWGCNA ],
            suffix( ".tsv.gz" ), 
            "_gsg.csv.gz" )
def findGoodSamplesGenes( infile, outfile ):
    out_rejected = P.snip( outfile, ".csv.gz" ) + "_rejected.csv.gz" 
    tmp_file = P.getTempFilename( "." )

    P10.wgcnaRunGoodSampleGenes( infile, tmp_file, out_rejected )

    # because write table misses a column name for rownames
    to_cluster = False
    statement = ( "zcat %(tmp_file)s |"
                  " sed '1 s/^/gene_id,/' |"
                  " gzip"
                  " > %(outfile)s" )
    P.run()
    os.unlink( tmp_file )


@jobs_limit( 1, "wgcna" )
#@follows( loadWGCNALibraries )
@split( findGoodSamplesGenes, 
            regex( "(.+)_gsg.csv.gz" ), 
            [ r"\1_sample_clustering.png", 
              r"\1_sample_clustering_sampleTree.rds",
              r"\1_sample_clustering_dat.rds"] )
def clusterSamples( infile, outfiles ):
    outf_png, outf_tree, outf_dat = outfiles
    P10.wgcnaClusterSamples( infile, outf_png, outf_tree, outf_dat )


@jobs_limit( 1, "wgcna" )
@collate( clusterSamples,
          regex( "(.+)_sample_clustering_(.+).rds" ), 
          r"\1_DF.rds" )
def removeOutlierSamples( infiles, outfile ):
    # assign infiles based on their expected filenames
    for infile in infiles:
        if re.search( "sampleTree.rds", infile ):
            in_tree = infile
        elif re.search( "dat.rds", infile ):
            in_dat = infile
        else: 
            continue

    # Set parameters for removing outlying samples
    cut_height = 15000
    min_size = 10

    # the log(x+1) transformed data and rlog transformed data
    # have different cut parameters
    if re.search( "_l_", os.path.basename( in_dat ) ):
        cut_height = False
    elif re.search( "rlog", os.path.basename( in_dat ) ):
        cut_height = False

    # pass parameters to function to remove outliers
    P10.wgcnaRemoveOutlierSamples( in_tree, 
                                   in_dat, 
                                   outfile, 
                                   cut_height = cut_height,
                                   min_size = min_size )
    
@jobs_limit( 1, "wgcna" )
@split( removeOutlierSamples,
        regex( "(.+)_DF.rds" ),
            [ r"\1_scaleIndependence.png",
              r"\1_meanConnectivity.png" ] )
def pickSoftThreshold( infile, outfiles ):
    out_scaleInd, out_meanCon = outfiles

    # set a vector containing potential soft thresholds
    # if PARAMS[""]:
    #     powers = PARAMS[""]
    # else:
    #     powers = range(1, 12)
    #     powers.extend( range(12, 25, 2) )
    powers = range(1, 12)
    powers.extend( range(12, 25, 2) )

    # default threshold for selecting power
    threshold = 0.9

    # create diagnostic plots for soft threshold selection
    P10.wgcnaPlotSoftThreshold( infile, 
                                powers, 
                                threshold, 
                                out_scaleInd, 
                                out_meanCon )



#@follows( loadWGCNALibraries )
@jobs_limit( 1, "wgcna" )
@follows( pickSoftThreshold )
@split( removeOutlierSamples, 
        regex( "(.+)_DF.rds" ),
        [ r"\1_ADJ_matrix.rds",
          r"\1_TOM_matrix.rds", 
          r"\1_scaledCON_hist.png" ] )
#          r"\1_NC.rds" ] )
def calcWeightedGeneCoexpressionNetwork( infile, outfiles ):
    out_adj, out_tom, out_scaledCon = outfiles # , out_nc = outfiles

    # set a power for creating scale free network
    inf = os.path.basename( infile )
    if re.search( "cuffdiff_wgcna_high_fpkm", inf ):
        power = PARAMS["wgcna_power_fpkm"]
    elif re.search( "cuffdiff_l_wgcna_high_fpkm", inf ):
        power = PARAMS["wgcna_power_l_fpkm"]
    elif re.search( "cuffdiff_wgcna_high_cv", inf ):
        power = PARAMS["wgcna_power_high_cv"]
    elif re.search( "cuffdiff_l_wgcna_high_cv", inf ):
        power = PARAMS["wgcna_power_l_high_cv"]
    elif re.search( "wgcna_rlog", inf ):
        power = PARAMS["wgcna_power_rlog"]
    else:
        raise IOError( "Unrecognised infile %s" % os.path.abspath( infile ) )

    # calculate the required matrices and diagnostics
    P10.wgcnaCalcCoexpressionNetwork( infile, 
                                      out_adj, 
                                      out_tom, 
                                      out_scaledCon, 
                                      power )


#@follows( loadWGCNALibraries )
@jobs_limit( 1, "wgcna" )
@split( calcWeightedGeneCoexpressionNetwork,
            regex( "(.+)_TOM_matrix.rds" ),
            [ r"\1_TOMdist_geneTree.rds", 
              r"\1_TOMdist_geneTree.png" ] )
#              r"\1_TOMdist_plot.png" ] )
def clusterTOM( infile, outfiles ):
    out_geneTree, out_dendro = outfiles # , out_TOMplot = outfiles

    P10.wgcnaClusterTOM( infile, out_geneTree, out_dendro )
    # Option to do this with P.Submit is disabled.
    #  params = [ infile, out_geneTree, out_dendro ]
    #  log = out_geneTree + ".log"


    #  P.submit( "/ifs/devel/projects/proj010/PipelineProj010",
    #            "wgcnaClusterTOM",
    #            params,
    #            logfile = log, 
    #            jobOptions = "-l mem_free=10G" )


# 1) cutreeDynamicTree with recursive module assignment dissabled
@jobs_limit( 1, "wgcna" )
#@follows( loadWGCNALibraries )
@split( ( calcWeightedGeneCoexpressionNetwork, clusterTOM ),
          regex( "(.+)_TOMdist_geneTree.rds" ), 
          add_inputs( r"\1_TOM_matrix.rds" ), 
          [ r"\1_dynamic_false_moduleColours.rds", 
            r"\1_dynamic_false_dendroAndColours.png" ] )
def treeCut_depthF( infiles, outfiles ):
    """
    Run cutreeDynamic with method = tree and deepSplit = FALSE
    """
    in_tree, in_tom = infiles
    out_modules, out_tree_png = outfiles

    # specify other parameters
    method = "tree"
    deep_split = False
    min_cluster_size = PARAMS["wgcna_treecut_depthf_mincluster"]

    P10.wgcnaCutreeDynamic( in_tree, 
                            in_tom, 
                            method, 
                            out_modules, 
                            out_tree_png, 
                            deep_split, 
                            min_cluster_size )

@jobs_limit( 1, "wgcna" )
@split( ( treeCut_depthF, removeOutlierSamples ), 
          regex( "(.+)_dynamic_false_moduleColours.rds" ), 
          add_inputs( r"\1_DF.rds" ),
          [ r"\1_dynamic_false_EG1.rds", r"\1_dynamic_false_EG1_cluster.png" ] )
def calcEigengenes_treeCut_depthF( infiles, outfiles ):
    """
    Run moduleEigengenes to obtain list of characteristics for the first principal 
    component of each module. 
    """
    in_modules, in_dat = infiles
    out_eigen, out_eigenplot = outfiles

    # set parameters
    merge_dist = PARAMS["wgcna_treecut_depthf_mergedist"]

    P10.wgcnaModuleEigengenes( in_dat,
                               in_modules, 
                               out_eigen, 
                               out_eigenplot, 
                               merge_dist )

@jobs_limit( 1, "wgcna" )
@follows( calcEigengenes_treeCut_depthF )
@split( ( treeCut_depthF, removeOutlierSamples ), 
        regex( "(.+)_dynamic_false_moduleColours.rds" ), 
        add_inputs( r"\1_DF.rds" ), 
        [ r"\1_dynamic_false_mergedModuleColours.rds", 
          r"\1_dynamic_false_merged_EG1.rds" ] )
def merge_treeCut_depthF( infiles, outfiles ):
    """
    Merge similar modules based on the pearson correlation of their eigengenes
    """
    in_modules, in_dat = infiles
    out_modules, out_eigen = outfiles

    # set parameters
    merge_dist = PARAMS["wgcna_treecut_depthf_mergedist"]

    P10.wgcnaMergeCloseModules( in_dat, 
                                in_modules, 
                                out_modules, 
                                out_eigen, 
                                merge_dist )
    
@jobs_limit( 1, "wgcna" )
@collate( [ clusterTOM, merge_treeCut_depthF ], 
          regex( "(.+)_TOMdist_geneTree.rds" ), 
          add_inputs( r"\1_dynamic_false_mergedModuleColours.rds" ), 
          r"\1_dynamic_false_mergedDendroAndColours.png" )
def plot_dendro_merge_treeCut_depthF( infiles, outfile ):
    """
    Plot a dendrogram for the merged modules
    """
    in_geneTree, in_modules = infiles[0]

    P10.wgcnaPlotDendroAndColors( in_geneTree, 
                                  in_modules, 
                                  outfile, 
                                  groupLabels = '"Dynamic\nFalse"' )


# 2) cutreeDynamicTree with recursive module assignment enabled
###@follows( plot_dendro_merge_treeCut_depthF )
#@follows( loadWGCNALibraries )
@jobs_limit( 1, "wgcna" )
@split( ( calcWeightedGeneCoexpressionNetwork, clusterTOM ),
          regex( "(.+)_TOMdist_geneTree.rds" ), 
          add_inputs( r"\1_TOM_matrix.rds" ), 
          [ r"\1_dynamic_true_moduleColours.rds", 
            r"\1_dynamic_true_dendroAndColours.png" ] )
def treeCut_depthT( infiles, outfiles ):
    """
    Run cutreeDynamic with method = tree and deepSplit = TRUE
    """
    # when using collate and add_inputs the infiles are a nested list
    in_tree, in_tom = infiles
    out_modules, out_tree_png = outfiles

    # specify other parameters
    method = "tree"
    deep_split = True
    min_cluster_size = PARAMS["wgcna_treecut_deptht_mincluster"]

    P10.wgcnaCutreeDynamic( in_tree, 
                            in_tom, 
                            method, 
                            out_modules, 
                            out_tree_png, 
                            deep_split, 
                            min_cluster_size )

@jobs_limit( 1, "wgcna" )
@split( ( treeCut_depthT, removeOutlierSamples ), 
          regex( "(.+)_dynamic_true_moduleColours.rds" ), 
          add_inputs( r"\1_DF.rds" ),
          [ r"\1_dynamic_true_EG1.rds", r"\1_dynamic_true_EG1_cluster.png" ] )
def calcEigengenes_treeCut_depthT( infiles, outfiles ):
    """
    Run moduleEigengenes to obtain list of characteristics for the first principal 
    component of each module. 
    """
    in_modules, in_dat = infiles
    out_eigen, out_eigenplot = outfiles

    # set parameters
    merge_dist = PARAMS["wgcna_treecut_deptht_mergedist"]

    P10.wgcnaModuleEigengenes( in_dat,
                               in_modules, 
                               out_eigen, 
                               out_eigenplot, 
                               merge_dist )

@jobs_limit( 1, "wgcna" )
@follows( calcEigengenes_treeCut_depthT )
@split( ( treeCut_depthT, removeOutlierSamples ), 
        regex( "(.+)_dynamic_true_moduleColours.rds" ), 
        add_inputs( r"\1_DF.rds" ), 
        [ r"\1_dynamic_true_mergedModuleColours.rds", 
          r"\1_dynamic_true_merged_EG1.rds" ] )
def merge_treeCut_depthT( infiles, outfiles ):
    """
    Merge similar modules based on the pearson correlation of their eigengenes
    """
    in_modules, in_dat = infiles
    out_modules, out_eigen = outfiles

    # set parameters
    merge_dist = PARAMS["wgcna_treecut_deptht_mergedist"]

    P10.wgcnaMergeCloseModules( in_dat, 
                                in_modules, 
                                out_modules, 
                                out_eigen, 
                                merge_dist )
    
@jobs_limit( 1, "wgcna" )
@collate( [ clusterTOM, merge_treeCut_depthT ], 
          regex( "(.+)_TOMdist_geneTree.rds" ), 
          add_inputs( r"\1_dynamic_true_mergedModuleColours.rds" ), 
          r"\1_dynamic_true_mergedDendroAndColours.png" )
def plot_dendro_merge_treeCut_depthT( infiles, outfile ):
    """
    Plot a dendrogram for the merged modules
    """
    in_geneTree, in_modules = infiles[0]

    P10.wgcnaPlotDendroAndColors( in_geneTree, 
                                  in_modules, 
                                  outfile, 
                                  groupLabels = '"Dynamic\nTrue"' )


# 3) cutreeHybrid with deepSplit set Low
###@follows( plot_dendro_merge_treeCut_depthF )
#@follows( loadWGCNALibraries )
@jobs_limit( 1, "wgcna" )
@split( ( calcWeightedGeneCoexpressionNetwork, clusterTOM ),
          regex( "(.+)_TOMdist_geneTree.rds" ), 
          add_inputs( r"\1_TOM_matrix.rds" ), 
          [ r"\1_hybrid_low_moduleColours.rds", 
            r"\1_hybrid_low_dendroAndColours.png" ] )
def hybridCut_depthLow( infiles, outfiles ):
    """
    Run cutreeDynamic with method = hybrid and deepSplit set low
    """
    # when using collate and add_inputs the infiles are a nested list
    in_tree, in_tom = infiles
    out_modules, out_tree_png = outfiles

    # specify other parameters
    method = "hybrid"
    deep_split = PARAMS["wgcna_hybridcut_depthlow_deepsplit"]
    min_cluster_size = PARAMS["wgcna_hybridcut_depthlow_mincluster"]

    P10.wgcnaCutreeDynamic( in_tree, 
                            in_tom, 
                            method, 
                            out_modules, 
                            out_tree_png, 
                            deep_split, 
                            min_cluster_size )

@jobs_limit( 1, "wgcna" )
@split( ( hybridCut_depthLow, removeOutlierSamples ), 
          regex( "(.+)_hybrid_low_moduleColours.rds" ), 
          add_inputs( r"\1_DF.rds" ),
          [ r"\1_hybrid_low_EG1.rds", r"\1_hybrid_low_EG1_cluster.png" ] )
def calcEigengenes_hybridCut_depthLow( infiles, outfiles ):
    """
    Run moduleEigengenes to obtain list of characteristics for the first principal 
    component of each module. 
    """
    in_modules, in_dat = infiles
    out_eigen, out_eigenplot = outfiles

    # set parameters
    merge_dist = PARAMS["wgcna_hybridcut_depthlow_mergedist"]

    P10.wgcnaModuleEigengenes( in_dat,
                               in_modules, 
                               out_eigen, 
                               out_eigenplot, 
                               merge_dist )


@jobs_limit( 1, "wgcna" )
@follows( calcEigengenes_hybridCut_depthLow )
@split( ( hybridCut_depthLow, removeOutlierSamples ), 
        regex( "(.+)_hybrid_low_moduleColours.rds" ), 
        add_inputs( r"\1_DF.rds" ), 
        [ r"\1_hybrid_low_mergedModuleColours.rds", 
          r"\1_hybrid_low_merged_EG1.rds" ] )
def merge_hybridCut_depthLow( infiles, outfiles ):
    """
    Merge similar modules based on the pearson correlation of their eigengenes
    """
    in_modules, in_dat = infiles
    out_modules, out_eigen = outfiles

    # set parameters
    merge_dist = PARAMS["wgcna_hybridcut_depthlow_mergedist"]

    P10.wgcnaMergeCloseModules( in_dat, 
                                in_modules, 
                                out_modules, 
                                out_eigen, 
                                merge_dist )
    

@jobs_limit( 1, "wgcna" )
@collate( [ clusterTOM, merge_hybridCut_depthLow ], 
          regex( "(.+)_TOMdist_geneTree.rds" ), 
          add_inputs( r"\1_hybrid_low_mergedModuleColours.rds" ), 
          r"\1_hybrid_low_mergedDendroAndColours.png" )
def plot_dendro_merge_hybridCut_depthLow( infiles, outfile ):
    """
    Plot a dendrogram for the merged modules
    """
    in_geneTree, in_modules = infiles[0]

    P10.wgcnaPlotDendroAndColors( in_geneTree, 
                                  in_modules, 
                                  outfile, 
                                  groupLabels = '"Hybrid\nLow"' )


# 4) cutreeHybrid with deepSplit set High
###@follows( plot_dendro_merge_hybridCut_depthLow )
#@follows( loadWGCNALibraries )
@jobs_limit( 1, "wgcna" )
@split( ( calcWeightedGeneCoexpressionNetwork, clusterTOM ),
          regex( "(.+)_TOMdist_geneTree.rds" ), 
          add_inputs( r"\1_TOM_matrix.rds" ), 
          [ r"\1_hybrid_high_moduleColours.rds", 
            r"\1_hybrid_high_dendroAndColours.png" ] )
def hybridCut_depthHigh( infiles, outfiles ):
    """
    Run cutreeDynamic with method = hybrid and deepSplit set high
    """
    # when using collate and add_inputs the infiles are a nested list
    in_tree, in_tom = infiles
    out_modules, out_tree_png = outfiles

    # specify other parameters
    method = "hybrid"
    deep_split = PARAMS["wgcna_hybridcut_depthhigh_deepsplit"]
    min_cluster_size = PARAMS["wgcna_hybridcut_depthhigh_mincluster"]

    P10.wgcnaCutreeDynamic( in_tree, 
                            in_tom, 
                            method, 
                            out_modules, 
                            out_tree_png, 
                            deep_split, 
                            min_cluster_size )


@jobs_limit( 1, "wgcna" )
@split( ( hybridCut_depthHigh, removeOutlierSamples ), 
          regex( "(.+)_hybrid_high_moduleColours.rds" ), 
          add_inputs( r"\1_DF.rds" ),
          [ r"\1_hybrid_high_EG1.rds", r"\1_hybrid_high_EG1_cluster.png" ] )
def calcEigengenes_hybridCut_depthHigh( infiles, outfiles ):
    """
    Run moduleEigengenes to obtain list of characteristics for the first principal 
    component of each module. 
    """
    in_modules, in_dat = infiles
    out_eigen, out_eigenplot = outfiles

    # set parameters
    merge_dist = PARAMS["wgcna_hybridcut_depthhigh_mergedist"]

    P10.wgcnaModuleEigengenes( in_dat,
                               in_modules, 
                               out_eigen, 
                               out_eigenplot, 
                               merge_dist )


@jobs_limit( 1, "wgcna" )
@follows( calcEigengenes_hybridCut_depthHigh )
@split( ( hybridCut_depthHigh, removeOutlierSamples ), 
        regex( "(.+)_hybrid_high_moduleColours.rds" ), 
        add_inputs( r"\1_DF.rds" ), 
        [ r"\1_hybrid_high_mergedModuleColours.rds", 
          r"\1_hybrid_high_merged_EG1.rds" ] )
def merge_hybridCut_depthHigh( infiles, outfiles ):
    """
    Merge similar modules based on the pearson correlation of their eigengenes
    """
    in_modules, in_dat = infiles
    out_modules, out_eigen = outfiles

    # set parameters
    merge_dist = PARAMS["wgcna_hybridcut_depthhigh_mergedist"]

    P10.wgcnaMergeCloseModules( in_dat, 
                                in_modules, 
                                out_modules, 
                                out_eigen, 
                                merge_dist )
    

@jobs_limit( 1, "wgcna" )
@collate( [ clusterTOM, merge_hybridCut_depthHigh ], 
          regex( "(.+)_TOMdist_geneTree.rds" ), 
          add_inputs( r"\1_hybrid_high_mergedModuleColours.rds" ), 
          r"\1_hybrid_high_mergedDendroAndColours.png" )
def plot_dendro_merge_hybridCut_depthHigh( infiles, outfile ):
    """
    Plot a dendrogram for the merged modules
    """
    in_geneTree, in_modules = infiles[0]

    P10.wgcnaPlotDendroAndColors( in_geneTree, 
                                  in_modules, 
                                  outfile, 
                                  groupLabels = '"Hybrid\nHigh"' )


#@follows( loadWGCNALibraries )
@jobs_limit( 1, "wgcna" )
@collate( [ clusterTOM, 
            merge_treeCut_depthF, 
            merge_treeCut_depthT, 
            merge_hybridCut_depthLow,
            merge_hybridCut_depthHigh ], 
          regex( "(.+)_TOMdist_geneTree.rds" ), 
          add_inputs( r"\1_dynamic_false_mergedModuleColours.rds", 
                      r"\1_dynamic_true_mergedModuleColours.rds", 
                      r"\1_hybrid_low_mergedModuleColours.rds", 
                      r"\1_hybrid_high_mergedModuleColours.rds" ),
          r"\1_treeCutComparison.png" )
def compareMergedModules( infiles, outfile):
    in_geneTree = infiles[0][0]
    in_modules = infiles[0][1:]

    group_labels = 'c("Dynamic_false","Dynamic_true","Hybrid_low","Hybrid_high")'

    P10.wgcnaPlotDendroAndColors( in_geneTree, 
                                  in_modules, 
                                  outfile, 
                                  groupLabels = group_labels, 
                                  main = '"Comparison of TreeCut approaches"' )


@follows( plot_dendro_merge_hybridCut_depthHigh,
          plot_dendro_merge_hybridCut_depthLow,
          plot_dendro_merge_treeCut_depthT,
          plot_dendro_merge_treeCut_depthF )
def runWGCNA():
    pass



#################################################################################
## subsection: output module information for dynamic_false
#################################################################################

#@follows( loadWGCNALibraries )
@jobs_limit( 1, "wgcna" )
@collate( [ merge_treeCut_depthF, findGoodSamplesGenes ],
            regex( r"(.+)/(.+)_dynamic_false_mergedModuleColours.rds" ), 
            add_inputs( r"\1/\2_gsg.csv.gz" ),
            r"\1/\2_module_assignment.tsv.gz" )
def extractModuleAssignment( infiles, outfile ):
    in_modules, in_dat = infiles[0]
    tmpf = P.getTempFilename( "." )
    
    P10.wgcnaExtractModuleAssignment( in_modules, in_dat, tmpf )

    to_cluster = False
    statement = ( "cat %(tmpf)s |"
                  " sed s/TRUE/1/g |"
                  " sed s/FALSE/0/g |"
                  " gzip > %(outfile)s" )
    P.run()
    os.unlink( tmpf )

    
@transform( extractModuleAssignment, 
            regex( "(.+)/(.+).tsv.gz" ), 
            r"./\2.load" )
def loadModuleAssignment( infile, outfile ):
    P.load( infile, outfile )


@transform(loadModuleAssignment,
           suffix("rlog_module_assignment.load"),
           "rlog_module_assignment_flattened.load")
def flattenModuleAssignment(infile, outfile):
    """
    """
    tmpf = P.getTempFilename(".")
    table = P.snip(os.path.basename(infile), ".load")
    out_tab = P.snip(os.path.basename(outfile), ".load")

    statement = "SELECT * FROM %s" % table
    df = PU.fetch_DataFrame(statement)
    df.set_index("gene_id", inplace=True)

    outf = IOTools.openFile(tmpf, "w")
    outf.write("gene_id\tmodule\n")
    # returns tuple with index as [0] and row val series as [1]
    for row in df.iterrows():
        gene_id = row[0]
        module = row[1][row[1] == 1]
        assert len(module.index) == 1, "Gene assigned to multiple modules..."
        module = module.index[0]
        outf.write(gene_id + "\t" + module + "\n")
    outf.close()

    P.load(tmpf, outfile, options="--table=%s" % out_tab)


@jobs_limit( 1, "wgcna" )
@follows( loadModuleAssignment,
          loadSummaryRlogTransformedCountData, # median rlog transformed counts
          loadReadCounts, # median FPKMs
          mkdir( "expression_wgcna_filtered_fpkm/refcoding_cuffdiff_l_wgcna_high_fpkm_module_heatmaps_zscore" ),
          mkdir( "expression_wgcna_filtered_fpkm/refcoding_cuffdiff_wgcna_high_fpkm_module_heatmaps_zscore" ),
          mkdir( "expression_wgcna_filtered_high_cv/refcoding_cuffdiff_l_wgcna_high_cv_module_heatmaps_zscore" ),
          mkdir( "expression_wgcna_filtered_high_cv/refcoding_cuffdiff_wgcna_high_cv_module_heatmaps_zscore" ),
          mkdir( "expression_wgcna_filtered_rlog/refcoding_wgcna_rlog_module_heatmaps_zscore" ) )
@split( extractModuleAssignment, 
        regex("(.+)/(.+)_module_assignment.tsv.gz"),
        add_inputs( extractMedianRlogTransformedCountData,
                    extractMedianFPKMs,
                    extractMedianLogTransformedFPKMs ), 
        r"\1/\2_module_heatmaps_zscore/clusteredHeatmapModule*.png" )
def plotModuleHeatmaps_zscore( infiles, outfiles ):
    """
    Retrieve module assignment and expression value, 
    plot heatmap for each module
    """
#    P.submit( "/ifs/devel/projects/proj010/PipelineProj010",
#              "plotZscoreHeatmaps",
#              infiles=infiles,
#              outfiles=outfiles )
    in_modules = P.snip(os.path.basename( infiles[0] ), ".tsv.gz" )


    infile = os.path.basename( infiles[0] )
    if re.match( "refcoding_cuffdiff_l_wgcna_high_fpkm", infile ):
        infile_2 = infiles[3]
        assert infile_2.endswith("_median_l_fpkms.tsv.gz")
    elif re.match( "refcoding_cuffdiff_wgcna_high_fpkm", infile ):
        infile_2 = infiles[2]
        assert infile_2.endswith("_median_fpkms.tsv.gz")
    elif re.match( "refcoding_cuffdiff_l_wgcna_high_cv", infile ):
        infile_2 = infiles[3]
        assert infile_2.endswith("_median_l_fpkms.tsv.gz")
    elif re.match( "refcoding_cuffdiff_wgcna_high_cv", infile ):
        infile_2 = infiles[2]
        assert infile_2.endswith("_median_fpkms.tsv.gz")
    elif re.match( "refcoding_wgcna_rlog", infile ):
        infile_2 = infiles[1]
        assert infile_2.endswith("rlog_trans_median_counts.tsv.gz")
    else:
        raise IOError( "Unrecognised infile %s" % os.path.abspath( infile ) )
        
    in_expression = P.snip(os.path.basename( infile_2 ), ".tsv.gz" )
    out_stub = P.snip(infiles[0], "_module_assignment.tsv.gz") + "_module_heatmaps_zscore"

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
        if infile_2.endswith("_median_l_fpkms.tsv.gz"):
            statement = ( "SELECT"
                          " a.gene_id,"
                          " a.Bcell_pro AS pro,"
                          " a.Bcell_pre AS pre,"
                          " a.Bcell_immature AS immature," 
                          " a.Bcell_mature AS mature,"
                          " a.Bcell_follicular AS follicular," 
                          " a.Bcell_marginal AS marginal,"
                          " a.Bcell_b1a AS b1a,"
                          " a.Bcell_germinal AS germinal"
                          " FROM %(in_expression)s AS a"
                          " INNER JOIN %(in_modules)s AS b"
                          " ON a.gene_id = b.gene_id"
                          " WHERE b.%(module)s = 1" % locals() )
        else:
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
        # df = pandas2ri.py2ri( df )
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


# @follows( loadWGCNALibraries )
# @jobs_limit( 1, "wgcna" )
# @follows( mkdir( "./expression_wgcna_filtered_fpkm/module_heatmaps" ), 
#           mkdir( "./expression_wgcna_filtered_high_cv/module_heatmaps" ) )
# @split( [clusterSamples, findGoodSamplesGenes, merge_treeCut_depthF ], 
#           regex( r"(.+)/(.+)_sample_clustering_sampleTree.rds" ), 
#           add_inputs( r"\1/\2_gsg.csv.gz", 
#                       r"\1/\2_dynamic_false_mergedModuleColours.rds" ),
#           r"\1/module_heatmaps/clusteredHeatmapModule*" )
# def plotModuleHeatmaps( infiles, outfiles ):
#     in_tree, in_dat, in_modules = infiles
#     out_dir = "./expression_wgcna_filtered_dynamicTrue_modules/module_heatmaps"
    
#     P10.wgcnaPlotModuleSpecificHeatmaps( in_dat, in_modules, in_tree, out_dir )


@jobs_limit( 1, "wgcna" )
@follows( mkdir( "expression_wgcna_filtered_fpkm/refcoding_cuffdiff_l_wgcna_high_fpkm_module_eigenplots" ),
          mkdir( "expression_wgcna_filtered_fpkm/refcoding_cuffdiff_wgcna_high_fpkm_module_eigenplots" ),
          mkdir( "expression_wgcna_filtered_high_cv/refcoding_cuffdiff_l_wgcna_high_cv_module_eigenplots" ),
          mkdir( "expression_wgcna_filtered_high_cv/refcoding_cuffdiff_wgcna_high_cv_module_eigenplots" ),
          mkdir( "expression_wgcna_filtered_rlog/refcoding_wgcna_rlog_module_eigenplots" ) )
@split( [merge_treeCut_depthF, findGoodSamplesGenes ],
        regex( r"(.+)/(.+)_dynamic_false_merged_EG1.rds" ), 
        add_inputs( r"\1/\2_gsg.csv.gz" ),
        r"\1/\2_module_eigenplots/moduleEigenGene*" )
def plotModuleEigengenes( infiles, outfiles ):
    in_eig, in_dat = infiles
    out_dir = P.snip( in_dat, "_gsg.csv.gz" ) + "_module_eigenplots"

    # based on removeOutlierSamples, the transformed data have no samples to 
    # remove, whilst the untransformed do...
    if re.search( "cuffdiff_l_wgcna", os.path.basename( in_eig ) ):
        to_remove = []
    elif re.search( "wgcna_rlog", os.path.basename( in_eig ) ):
        to_remove = []
    elif re.search( "cuffdiff_wgcna", os.path.basename( in_eig ) ):
        to_remove = PARAMS[ "wgcna_to_remove" ].split(",")
    else:
        raise IOError( "Unrecognised infile %s" % os.path.abspath( in_eig ) )

    P10.wgcnaPlotModuleEigengenes( in_eig, in_dat, out_dir, to_remove )


@jobs_limit( 1, "wgcna" )
@collate( [ findGoodSamplesGenes, merge_treeCut_depthF ], 
          regex( "(.+)/(.+)_gsg.csv.gz" ), 
          add_inputs( r"\1/\2_dynamic_false_merged_EG1.rds" ), 
          r"\1/\2_eigengenes.tsv" )
def extractModuleEigengenes( infiles, outfile ):
    in_csv, in_eigen = infiles[0]

    in_eig = os.path.basename( in_eigen )
    if re.search( "cuffdiff_l_wgcna", os.path.basename( in_eig ) ):
        to_remove = []
    elif re.search( "wgcna_rlog", os.path.basename( in_eig ) ):
        to_remove = []
    elif re.search( "cuffdiff_wgcna", os.path.basename( in_eig ) ):
        to_remove = PARAMS[ "wgcna_to_remove" ].split(",")
    else:
        raise IOError( "Unrecognised infile %s" % os.path.abspath( in_eig ) )

    P10.extractEigengenes_R( in_csv, in_eigen, to_remove, outfile )



# @transform( extractModuleEigengenes, 
#             regex( "(.+)/(.+).tsv" ), 
#             r"./\2.load" )
# def loadModuleEigengenes( infile, outfile ):
#     P.load( infile, outfile )


@transform( extractModuleEigengenes, suffix( ".tsv" ), "_stacked.tsv" )
def stackModuleEigengenes( infile, outfile ):
    """
    Use pandas to melt table of module eigengenes
    """    
    df = pd.read_table( infile, header = 0 )
    # melt, requires id_vars to be a column not the index
    df = pd.melt( df, 
                  id_vars = "module_id", 
                  var_name = "sample_id", 
                  value_name = "expression" )
    df = df.set_index( "module_id" )

    # Split the sample_id column, using Series.str.spit()
    # NB. sample_ids are in the form tissue.condition.replicate
    # This changes to tissue_condition_replicate in csvdb
    split_col = df["sample_id"].str.split(".")
    df["cell_type"] = split_col.str[1]
    df["replicate"] = split_col.str[2]

    # write table to outfile
    df.to_csv( outfile, sep="\t" )


@transform( stackModuleEigengenes, regex( "(?:.+)/(.+).tsv" ), r"./\1.load" )
def loadStackedModuleEigengenes( infile, outfile ):
    P.load( infile, outfile )


# NOT CURRENTLY IMPLEMENTED
# # correlate module eigengenes with RIN numbers
# @transform( loadModuleEigengenes, 
#             regex( "(.+)/(.+).load" ),
#             add_inputs( os.path.join( PARAMS[ "location_external_datafiles" ],
#                                       "rna_qc_rin.csv" )  ),
#             r"./expression_wgcna_filtered_dynamicTrue_modules/\2_RINcorrelations.tsv" )
# def calculateEigengeneRINCorrelation( infiles, outfile ):
#     in_eig, in_rin = infiles
#     out_dir = "./expression_wgcna_filtered_dynamicTrue_modules" 

#     # open infiles as dataframe, series
#     eig_table = P.snip( os.path.basename( in_eig ), ".load" ) 
#     statement = ( "SELECT * FROM %(eig_table)s;" % locals() )
#     eig_df = PU.fetch_DataFrame( statement )
#     eig_df = eig_df.set_index( "module_id" )

#     rin_df = pd.read_table( in_rin, 
#                             index_col = 0, 
#                             header = 0, 
#                             na_values = "NA", 
#                             sep = ",")
#     rin_df = rin_df.T
#     rin_vals = rin_df.loc["RIN"]

#     P10.calcEigenRINCorrelation( eig_df, rin_vals, out_dir, outfile )
    

# @transform( calculateEigengeneRINCorrelation, 
#             regex( "(.+)/(.+).tsv" ), 
#             r"./\2.load" )
# def loadEigenRINCorrelation( infile, outfile ):
#     P.load( infile, outfile )

#################################################################################
## subsection: GO for module eigengenes
#################################################################################
# This section only takes forward:
# i) rlog transformed data
# ii) dynamic tree cut results (with recursive clustering set to false)
@follows(mkdir("expression_wgcna_filtered_rlog/refcoding_wgcna_rlog_module_GO_enrichment"))
@split(loadModuleAssignment,
       regex(".*(refcoding_wgcna_rlog_module)_assignment.load"),
       r"expression_wgcna_filtered_rlog/\1_GO_enrichment/\1_*tsv")
def splitModulesForGOEnrichment(infile, outfiles):
    """
    Run GOseq with correction for gene length
    """
    out_dir = "expression_wgcna_filtered_rlog/refcoding_wgcna_rlog_module_GO_enrichment"
    table = P.snip(os.path.basename(infile), ".load")

    statement = "SELECT * FROM %(table)s" % locals()
    df = PU.fetch_DataFrame(statement)
    df.set_index("gene_id", inplace=True)

    for module in df.columns.tolist():
        E.info("Running GOseq for module: %s" % module)
        genes = pd.DataFrame(df[module])
        genes.columns = [module,]
        outf = os.path.join(out_dir, table + "_" + module + "_GO_enrichment.tsv")
        genes.to_csv(IOTools.openFile(outf, "w"), sep="\t")


# The following error is given when running getModuleGOEnrichment on grey module...
# Error in if (min(fv) < lower_bound) fv = fv - min(fv) + lower_bound : 
#   missing value where TRUE/FALSE needed


@jobs_limit( 1, "RGlobalEnv")
@transform(splitModulesForGOEnrichment,
           regex("(.+)/refcoding_wgcna_rlog_module_assignment_b(.+).tsv"),
           add_inputs(extractLocusCoordinates),
           r"\1/refcoding_wgcna_rlog_module_assignment_b\2.tsv.gz")
def getModuleGOEnrichment_b(infiles, outfile):
    """
    Run GOSeq for individual modules in series 
    """
    infile, gene_sizes = infiles
    gene_sizes = P.snip(os.path.basename(gene_sizes), ".tsv.gz")
    statement = ("SELECT gene_id,Location FROM %s" % gene_sizes)
    df = PU.fetch_DataFrame(statement)
    df.set_index("gene_id", inplace=True)
    P10.runGOSeq(infile, df, outfile) #, submit=True)


@jobs_limit( 1, "RGlobalEnv")
@follows(getModuleGOEnrichment_b)
@transform(splitModulesForGOEnrichment,
           regex("(.+)/refcoding_wgcna_rlog_module_assignment_dark(.+).tsv"),
           add_inputs(extractLocusCoordinates),
           r"\1/refcoding_wgcna_rlog_module_assignment_dark\2.tsv.gz")
def getModuleGOEnrichment_dark(infiles, outfile):
    """
    Run GOSeq for individual modules in series 
    """
    infile, gene_sizes = infiles
    gene_sizes = P.snip(os.path.basename(gene_sizes), ".tsv.gz")
    statement = ("SELECT gene_id,Location FROM %s" % gene_sizes)
    df = PU.fetch_DataFrame(statement)
    df.set_index("gene_id", inplace=True)
    P10.runGOSeq(infile, df, outfile) #, submit=True)


@jobs_limit( 1, "RGlobalEnv")
@follows(getModuleGOEnrichment_dark)
@transform(splitModulesForGOEnrichment,
           regex("(.+)/refcoding_wgcna_rlog_module_assignment_grey60(.+).tsv"),
           add_inputs(extractLocusCoordinates),
           r"\1/refcoding_wgcna_rlog_module_assignment_grey60\2.tsv.gz")
def getModuleGOEnrichment_grey(infiles, outfile):
    """
    Run GOSeq for individual modules in series 
    """
    infile, gene_sizes = infiles
    gene_sizes = P.snip(os.path.basename(gene_sizes), ".tsv.gz")
    statement = ("SELECT gene_id,Location FROM %s" % gene_sizes)
    df = PU.fetch_DataFrame(statement)
    df.set_index("gene_id", inplace=True)
    P10.runGOSeq(infile, df, outfile) #, submit=True)

@jobs_limit( 1, "RGlobalEnv")
@follows(getModuleGOEnrichment_grey)
@transform(splitModulesForGOEnrichment,
           regex("(.+)/refcoding_wgcna_rlog_module_assignment_f(.+).tsv"),
           add_inputs(extractLocusCoordinates),
           r"\1/refcoding_wgcna_rlog_module_assignment_f\2.tsv.gz")
def getModuleGOEnrichment_f(infiles, outfile):
    """
    Run GOSeq for individual modules in series 
    """
    infile, gene_sizes = infiles
    gene_sizes = P.snip(os.path.basename(gene_sizes), ".tsv.gz")
    statement = ("SELECT gene_id,Location FROM %s" % gene_sizes)
    df = PU.fetch_DataFrame(statement)
    df.set_index("gene_id", inplace=True)
    P10.runGOSeq(infile, df, outfile) #, submit=True)


@jobs_limit( 1, "RGlobalEnv")
@follows(getModuleGOEnrichment_f)
@transform(splitModulesForGOEnrichment,
           regex("(.+)/refcoding_wgcna_rlog_module_assignment_l(.+).tsv"),
           add_inputs(extractLocusCoordinates),
           r"\1/refcoding_wgcna_rlog_module_assignment_l\2.tsv.gz")
def getModuleGOEnrichment_l(infiles, outfile):
    """
    Run GOSeq for individual modules in series 
    """
    infile, gene_sizes = infiles
    gene_sizes = P.snip(os.path.basename(gene_sizes), ".tsv.gz")
    statement = ("SELECT gene_id,Location FROM %s" % gene_sizes)
    df = PU.fetch_DataFrame(statement)
    df.set_index("gene_id", inplace=True)
    P10.runGOSeq(infile, df, outfile) #, submit=True)


@jobs_limit( 1, "RGlobalEnv")
@follows(getModuleGOEnrichment_l)
@transform(splitModulesForGOEnrichment,
           regex("(.+)/refcoding_wgcna_rlog_module_assignment_m(.+).tsv"),
           add_inputs(extractLocusCoordinates),
           r"\1/refcoding_wgcna_rlog_module_assignment_m\2.tsv.gz")
def getModuleGOEnrichment_m(infiles, outfile):
    """
    Run GOSeq for individual modules in series 
    """
    infile, gene_sizes = infiles
    gene_sizes = P.snip(os.path.basename(gene_sizes), ".tsv.gz")
    statement = ("SELECT gene_id,Location FROM %s" % gene_sizes)
    df = PU.fetch_DataFrame(statement)
    df.set_index("gene_id", inplace=True)
    P10.runGOSeq(infile, df, outfile) #, submit=True)


@jobs_limit( 1, "RGlobalEnv")
@follows(getModuleGOEnrichment_m)
@transform(splitModulesForGOEnrichment,
           regex("(.+)/refcoding_wgcna_rlog_module_assignment_p(.+).tsv"),
           add_inputs(extractLocusCoordinates),
           r"\1/refcoding_wgcna_rlog_module_assignment_p\2.tsv.gz")
def getModuleGOEnrichment_p(infiles, outfile):
    """
    Run GOSeq for individual modules in series 
    """
    infile, gene_sizes = infiles
    gene_sizes = P.snip(os.path.basename(gene_sizes), ".tsv.gz")
    statement = ("SELECT gene_id,Location FROM %s" % gene_sizes)
    df = PU.fetch_DataFrame(statement)
    df.set_index("gene_id", inplace=True)
    P10.runGOSeq(infile, df, outfile) #, submit=True)


@jobs_limit( 1, "RGlobalEnv")
@follows(getModuleGOEnrichment_p)
@transform(splitModulesForGOEnrichment,
           regex("(.+)/refcoding_wgcna_rlog_module_assignment_t(.+).tsv"),
           add_inputs(extractLocusCoordinates),
           r"\1/refcoding_wgcna_rlog_module_assignment_t\2.tsv.gz")
def getModuleGOEnrichment_t(infiles, outfile):
    """
    Run GOSeq for individual modules in series 
    """
    infile, gene_sizes = infiles
    gene_sizes = P.snip(os.path.basename(gene_sizes), ".tsv.gz")
    statement = ("SELECT gene_id,Location FROM %s" % gene_sizes)
    df = PU.fetch_DataFrame(statement)
    df.set_index("gene_id", inplace=True)
    P10.runGOSeq(infile, df, outfile) #, submit=True)


# @follows(mkdir("expression_wgcna_filtered_rlog/refcoding_wgcna_rlog_module_GO_enrichment"))
# @split(loadModuleAssignment,
#        regex(".*(refcoding_wgcna_rlog_module)_assignment.load"),
#        r"expression_wgcna_filtered_rlog/\1_GO_enrichment/\1_*tsv.gz")
# def getModuleGOEnrichment(infile, outfiles):
#     """
#     Run GOseq with correction for gene length
#     """
#     goseq = importr("goseq")
#     grdevices = importr("grDevices")

#     out_dir = "expression_wgcna_filtered_rlog/refcoding_wgcna_rlog_module_GO_enrichment"
#     table = P.snip(os.path.basename(infile), ".load")

#     statement = "SELECT * FROM %(table)s" % locals()
#     df = PU.fetch_DataFrame(statement)
#     df.set_index("gene_id", inplace=True)

#     background = map(str, df.index.tolist())
#     for module in df.columns.tolist():
#         E.info("Running GOseq for module: %s" % module)
#         outf = os.path.join(out_dir, table + "_" + module + "_GO_enrichment.tsv.gz")
#         out_plot = P.snip(outf, ".tsv.gz") + ".png"
#         grdevices.png(file=out_plot)

#         genes = df[module]
#         genes = robjects.IntVector(genes)
#         genes.names = background

#         pwf = goseq.nullp(genes, "mm10", "ensGene")
#         GO_wall = goseq.goseq(pwf, "mm10", "ensGene")

#         df_out = pandas2ri.ri2pandas(GO_wall)
#         df_out["module"] = len(df_out.index)*[module,]
#         df_out.to_csv(IOTools.openFile(outf, "w"), sep="\t")
#         grdevices.dev_off()


@collate([ getModuleGOEnrichment_b, 
           getModuleGOEnrichment_dark, 
           getModuleGOEnrichment_grey, 
           getModuleGOEnrichment_f, 
           getModuleGOEnrichment_l, 
           getModuleGOEnrichment_m, 
           getModuleGOEnrichment_p, 
           getModuleGOEnrichment_t ],
         regex("(.+)/(refcoding_wgcna_rlog_module)_.+_GO_enrichment.tsv.gz"),
         r"\1/\2_GO_enrichment_unadj.tsv.gz")
def collateModuleGOEnrichment(infiles, outfile):
    """
    """
    file_names = " ".join(infiles)
    statement = ("python %(scriptsdir)s/combine_tables.py"
                 "  --cat=Module_filename"
                 "  --log=%(outfile)s.log"
                 " %(file_names)s |"
                 " gzip > %(outfile)s")
    P.run()


@transform(collateModuleGOEnrichment, suffix("_unadj.tsv.gz"), ".tsv.gz")
def padjustModuleGOEnrichment(infile, outfile):
    """
    Run BH correction for GO enrichment pvalues. 
    This runs separate tests for over- & under-representation.
    """
    stats = importr("stats")
    df = pd.read_table(infile, compression="gzip", sep="\t")
#    df.drop("Unnamed: 1", inplace=True, axis=1)
    df.drop("Module_filename", inplace=True, axis=1)

    over_p = df["over_represented_pvalue"]
    over_p = robjects.FloatVector(over_p)
    over_p_adj = stats.p_adjust(over_p, method="BH")
    # numpy2ri.activate is implicit due to imported modules...
    # over_p_adj = rpyn.ri2numpy(over_p_adj)
    df["over_represented_padj"] = list(over_p_adj)

    under_p = df["under_represented_pvalue"]
    under_p = robjects.FloatVector(under_p)
    under_p_adj = stats.p_adjust(under_p, method="BH")
    # under_p_adj = rpyn.rpy2numpy(under_p_adj)
    df["under_represented_padj"] = list(under_p_adj)

    df.to_csv(IOTools.openFile(outfile, "w"),
              sep="\t",
              index=False)


@transform(padjustModuleGOEnrichment, regex(".+/(.+).tsv.gz"), r"\1.load")
def loadModuleGOEnrichment(infile, outfile):
    P.load(infile, outfile)


#################################################################################
## subsection: correlation between lncRNAs and module eigengenes
#################################################################################

@transform( extractModuleEigengenes, 
            regex( "(.+)_rlog_eigengenes.tsv" ), 
            add_inputs( extractPerSampleFPKMs ), 
            r"\1_rlog_eigen_lncRNA_pearsonCor.tsv.gz" )
def calculateLncRNAEigengeneCorrelation( infiles, outfile ):
    eig_table, fpkm_table = infiles
    outfile_rej = P.snip( outfile, ".tsv.gz" ) + "_rejected.tsv.gz"
    min_n = 20

    # retrieve lncRNA expression data
    # fpkm_table = P.snip( fpkm_table, ".load" )
    # regex = "LNC"
    # statement = ( "SELECT * FROM %(fpkm_table)s WHERE gene_id LIKE '%%LNC%%'" % locals() )
    # lnc_df = PU.fetch_DataFrame( statement )
    # lnc_df = lnc_df.set_index( "gene_id" )

    # hack... I can't subset a dataframe based on re match to index
    # the eigengene_ids are "." separated, whilst the refcoding_ids are "-" separatedx
    tmpf = P.getTempFilename( "." )
    to_cluster = False
    statement = ( "zcat %(fpkm_table)s |"
                  " awk 'NR==1 || $1 ~ /^LNC/' |"
                  " sed '1 s/Bcell-/Bcell\./g' |"
                  " sed '1 s/-R/\.R/g' "
                  "> %(tmpf)s" )
    P.run()

    lnc_df = pd.read_table( tmpf, 
                            index_col = 0, 
                            header = 0, 
                            na_filter = True, 
                            na_values = [ "NULL", "NA", "nan", "NaN", "N/A" ] )
    
    # retrieve eigengenes
    eig_df = pd.read_table( eig_table, 
                            index_col = 0,
                            header = 0,
                            na_filter = True,
                            na_values = [ "NULL", "NA", "nan", "NaN", "N/A" ] )

    E.info( "Eigen_df sample_ids: %s" % eig_df.columns )
    E.info( "LncRNA_df sample_ids: %s" % lnc_df.columns )

    # pass both dataframes to function that returns nested dict of pearson_r vals
    P10.calculateCorrelations_df( eig_df, lnc_df, min_n, outfile, outfile_rej )
    os.unlink( tmpf )

@transform( calculateLncRNAEigengeneCorrelation,
            regex( "(.+)/(.+).tsv.gz"),
            r"./\2.load" )
def loadLncRNAEigengeneCorrelation( infile, outfile ):
    P.load( infile, outfile )

#################################################################################
## subsection: find highest correlated wgcna module for each lncRNA
#################################################################################
@follows( loadLncRNAEigengeneCorrelation )
@transform( calculateLncRNAEigengeneCorrelation, 
            suffix( "_pearsonCor.tsv.gz" ), 
            "_highest_pearsonCor.tsv.gz" )
def findLncRNAHighestEigengeneCorrelation( infile, outfile ):
    """
    Create boolean table specifying highest correlating wgcna module for each
    lncRNA, provided highest correlation is above minimum threshold (0.8).
    """
    table = P.snip( os.path.basename( infile ), ".tsv.gz" )
    statement = ( "SELECT * FROM %s" % table )

    # retrieve dataframe of lncRNA vs module eigengene pearson correlations
    df = PU.fetch_DataFrame( statement )
    df.set_index( "gene_id", inplace = True )

    # specify the minium correlation threshold at which to associate a
    # lncRNA with a module
    #     min_cor = PARAMS["wgcna_min_cor"]
    min_cor = 0

    df_bool = df.apply( P10.returnHighestValue, axis = 1, min_thresh = min_cor )

    df_bool.to_csv( IOTools.openFile(outfile, "w"), sep = "\t" )



@transform( findLncRNAHighestEigengeneCorrelation, 
            regex( "(?:.+)/(.+).tsv.gz" ), 
            r"\1.load" )
def loadLncRNAHighestEigengeneCorrelation( infile, outfile ):
    P.load( infile, outfile, options="--add-index=gene_id" )


@transform(loadLncRNAEigengeneCorrelation,
           suffix(".load"),
           "_max_cor.load")
def loadLncRNAMaxEigengeneCorrelation(infile, outfile):
    """
    Iterate over the lncRNA eigengene correlations, load the 
    module and corr value for the module with the highest abs corr
    for each lncRNA. 
    """
    tmpf = P.getTempFile(".")
    table = P.snip(os.path.basename(infile), ".load")
    out_tab = P.snip(os.path.basename(outfile), ".load")

    statement = "SELECT * FROM %s" % table
    df = PU.fetch_DataFrame(statement)
    df.set_index("gene_id", inplace=True)

    tmpf.write("gene_id\tWGCNA_Module\tPearson_Correlation\n")
    for i in df.iterrows():
        gene_id = i[0]
        module = abs(i[1].astype("float")).idxmax()[2:]
        corr = i[1][abs(i[1]).idxmax()]
        tmpf.write("\t".join([gene_id, module, str(corr)]) + "\n")
    tmpf.close()
    tmpf = tmpf.name

    P.load(tmpf, outfile, options="--table=%s" % out_tab)
    os.unlink(tmpf)


@transform(loadLncRNAHighestEigengeneCorrelation, suffix(".load"), "_flattened.load")
def flattenLncRNAHighestEigengeneCorrelation(infile, outfile):
    """
    """
    tmpf = P.getTempFile(".")
    table = P.snip(os.path.basename(infile), ".load")
    out_tab = P.snip(os.path.basename(outfile), ".load")

    statement = "SELECT * FROM %s" % table
    df = PU.fetch_DataFrame(statement)
    df.set_index("gene_id", inplace=True)

    outf = IOTools.openFile(tmpf, "w")
    outf.write("lnc_id\tmodule\n")
    for row in df.iterrows():
        lnc_id = row[0]
        module = row[1][row[1] == "True"]
        assert len(module.index) < 2, "lncRNA assigned to multiple modules"
        if len(module.index) == 0:
            module = "None"
        else:
            module = module.index[0]
        outf.write(lnc_id + "\t" + module + "\n")
    outf.close()

    P.load(tmpf, outfile, options="--table=%s" % out_tab)
    os.unlink(tmpf)


@transform(loadNearestProteinCodingGenes, 
           regex(".*/lncRNA_nearest_refcoding_gene.load"),
           add_inputs( summarizelncRNARefcodingBiotypes,
                       flattenModuleAssignment,
                       flattenLncRNAHighestEigengeneCorrelation ),
           r"expression_wgcna_filtered_rlog/lncRNA_refcoding_module_overlap.tsv.gz")
def calculateLncRNARefcodingModuleOverlap(infiles, outfile):
    """
    """
    # get table names
    infiles = [P.snip(os.path.basename(x), ".load") for x in infiles]
    nn_genes, lnc_biotypes, pc_modules, lnc_modules = infiles
    
    # get lncRNA module assignment
    statement = "SELECT * FROM %s" % lnc_modules
    df_lnc = PU.fetch_DataFrame(statement)
    df_lnc.set_index("lnc_id", inplace=True)

    # get protein coding module assignment
    statement = "SELECT * FROM %s" % pc_modules
    df_pc = PU.fetch_DataFrame(statement)
    df_pc.set_index("gene_id", inplace=True)

    # get table of proximal protein coding genes
    statement = ("SELECT a.gene_id, a.id3, a.id5"
                 " FROM %(nn_genes)s AS a"
                 " INNER JOIN %(lnc_biotypes)s AS b"
                 " ON a.gene_id = b.gene_id"
                 " WHERE b.biotype = 'intergenic'" % locals())
    df_nn = PU.fetch_DataFrame(statement)
    df_nn.set_index("gene_id", inplace=True)

    # set up out table
    df_out = {"lnc_id": [], 
              "lnc_module": [], 
              "id3_module": [], 
              "id5_module": [], 
              "module_overlap": []}

    # set up outfile for catching lncRNAs without proximal coding loci
    out_fail = P.snip(outfile, ".tsv.gz") + "_failed.tsv.gz"
    outf = IOTools.openFile(out_fail, "w")
    outf.write("lncID\tID5\tID3\n")

    for lncRNA in df_nn.iterrows():
        # get lnc ID, 3'ID,  5'ID
        lnc_id = lncRNA[0]
        id3 = lncRNA[1].loc["id3"]
        id5 = lncRNA[1].loc["id5"]
        # some lncRNAs have no proximal genes...
        if id3 == None and id5 == None:
            E.warn("LncRNA %s has no proximal coding genes" % lnc_id)
            outf.write("\t".join([lnc_id, str(id5), str(id3)]) + "\n")

        # extract module affiliation for lncRNA 
        lnc_module = df_lnc.loc[lnc_id][0]
        # strip ME prefix!
        assert lnc_module.startswith("ME"), "No prefix..."
        lnc_module = lnc_module[2:]

        # extract 5' 3' module affiliation
        try:
            id3_module = df_pc.loc[id3][0]
        except KeyError:
            # proximal 3' gene not in wgcna analysis (low expression)
            id3_module = "None"
        except ValueError:
            # lncRNA has no recorded proximal 5' gene 
            assert id3 == None
            E.warn("LncRNA %s has no proximal 3' gene" % lnc_id)
            outf.write("\t".join([lnc_id, str(id5), str(id3)]) + "\n")
            id3_module = "None"
        try:
            id5_module = df_pc.loc[id5][0]
        except KeyError:
            # proximal 5' gene not in wgcna analysis (low expression)
            id5_module = "None"
        except ValueError:
            # lncRNA has no recorded proximal 5' gene 
            assert id5 == None
            outf.write("\t".join([lnc_id, str(id5), str(id3)]) + "\n")
            E.warn("LncRNA %s has no proximal 5' gene" % lnc_id)
            id_module = "None"

        # add module affiliation to df_out
        df_out["lnc_id"].append(lnc_id)
        df_out["lnc_module"].append(lnc_module)
        df_out["id3_module"].append(id3_module)
        df_out["id5_module"].append(id5_module)
        if lnc_module != "None":
            if lnc_module in [id3_module, id5_module]:
                df_out["module_overlap"].append(1)
            else:
                df_out["module_overlap"].append(0)
        else:
            df_out["module_overlap"].append(0)
    
    df_out = pd.DataFrame(df_out)
    df_out.to_csv(IOTools.openFile(outfile, "w"), sep="\t", index=False)
    outf.close()


@transform(calculateLncRNARefcodingModuleOverlap,
           regex(".+/(.+).tsv.gz"),
           r"\1.load")
def loadLncRNARefcodingModuleOverlap(infile, outfile):
    P.load(infile, outfile)


###############################################################################
## subsection: Find module assoc. lncRNAs proximal to module assoc. GO genes
###############################################################################

@transform(loadModuleGOEnrichment,
           regex("(.+)_module_GO_enrichment.load"),
           add_inputs(loadLncRNAEigengeneCorrelation, 
                      loadNearestProteinCodingGenes,
                      loadLncRNAChromatinClassifications),
           r"expression_wgcna_filtered_rlog/\1_lncRNA_with_proximal_GO_genes.tsv.gz")
def findProximalModuleGOGenes(infiles, outfile):
    """
    Select enriched GO terms for wgcna module, retrieve gene_ids for GO terms.
    Select all lncRNAs correlated (>0.8) with wgcna module, get proximal 
    protein coding genes. Intersect proximal genes with GO genes. 
    """

    def _getGOID(gene_id, id_table, term_table):
        """
        Receive ensembl gene_id and return GO ID and GO terms
        """
        ids = id_table.loc[gene_id]
        # if multiple entries in table, will return series, else string
        # if type(ids) == str:
        #     # query table mapping go ids to terms. 
        #     terms = term_table.loc[ids]
        #     go_ids = ids
        #     go_terms = terms
        # elif type(ids) == pd.core.series.Series:
        if type(ids) == pd.core.series.Series:
            # there is only one go term for gene_id
            ids = list(set(ids))
        elif type(ids) == pd.core.frame.DataFrame:
            # one gene_id maps to severl GO categories.
            ids = ids["go_id"]
        else: 
            raise IOError ("Unrecognised type returned by .loc %s" % ids.__class__)

        terms = []
        for i in ids:
            terms.append(str(term_table.loc[i][0]))
            go_ids = ",".join(ids)
            go_terms = ",".join(terms)

        return go_ids, go_terms

    infiles = [P.snip(os.path.basename(x), ".load") for x in infiles]
    GO_cat, eigen_cor, nn_genes, chr_class = infiles

    biomaRt = importr("biomaRt")

    # get modules
    statement = "SELECT module FROM %(GO_cat)s GROUP BY module" % locals()
    modules = PU.fetch(statement)
    modules = [str(x[0]) for x in modules]

    # set up a table that maps GO ids to GO terms
    statement = ("SELECT category, term FROM %(GO_cat)s" % locals())
    df_terms = PU.fetch_DataFrame(statement)
    df_terms.drop_duplicates(inplace=True)
    df_terms.set_index("category", inplace=True)

    # set up list to take details of lncRNAs with proximal genes in GO set
    out_list = []
    out_list_headers = ["lnc_id", 
                        "lnc_chromatin_classification", 
                        "pr_id", 
                        "module", 
                        "go_ids",
                        "go_terms"]

    for module in modules:
        # fetch enriched (p<0.05) GO terms for module
        statement = ("SELECT category from %(GO_cat)s"
                     " WHERE module = '%(module)s'"
                     " AND over_represented_padj < 0.05" % locals())
        go_terms = PU.fetch(statement)
        go_terms = robjects.StrVector([str(x[0]) for x in go_terms])
        E.info("Fetched %i GO terms for module %s" % (len(go_terms), module))

        if len(go_terms) == 0:
            E.warn("There are no enriched GO categories for module %s" % module)
            continue

        # fetch GO associated genes from biomart
        ensembl = biomaRt.useMart("ensembl", dataset="mmusculus_gene_ensembl")

        attributes = robjects.StrVector(["ensembl_gene_id", "go_id"])
        gene_data = biomaRt.getBM(attributes=attributes,
                                  filters="go_id",
                                  values=go_terms,
                                  mart=ensembl)
        gene_data = pandas2ri.ri2pandas(gene_data)

        E.info("There are %i GO genes for module %s" %(len(gene_data.index), module))

        module_go_genes = list(set(map(str, gene_data["ensembl_gene_id"].tolist())))

        # fetch lncRNAs that are highly correlated with module
        statement = ("SELECT a.gene_id,b.class_empirical AS Classification,c.id3,c.id5"
                     " FROM %(eigen_cor)s AS a"
                     " INNER JOIN %(chr_class)s AS b"
                     " ON a.gene_id = b.gene_id"
                     " INNER JOIN %(nn_genes)s As c"
                     " ON a.gene_id = c.gene_id"
                     " WHERE ABS(a.ME%(module)s) > 0.8" % locals())
        df_lncRNA = PU.fetch_DataFrame(statement)

        if len(df_lncRNA.index) == 0:
            E.warn("There are no highly correlated lncRNAs for module %s" % module)
            continue

        # set gene_data table for accessing go_id
        gene_data.set_index("ensembl_gene_id", inplace=True)

        success = 0
        # extract lncRNAs where nearest neighbour gene is in GO set
        for row in df_lncRNA.iterrows():
            lnc = map(str, row[1].tolist())
            lnc_id, cr_class, id3, id5 = lnc
            if id3 in module_go_genes:
                go_ids, go_terms = _getGOID(id3, gene_data, df_terms)
                out_list.append([lnc_id, cr_class, id3, module, go_ids, go_terms])
                success += 1
            if id5 in module_go_genes:
                go_ids, go_terms = _getGOID(id5, gene_data, df_terms)
                out_list.append([lnc_id, cr_class, id5, module, go_ids, go_terms])
                success += 1
            else:
                continue
        E.info("Finished searching for proximal genes in module %s,"
               " there were %i overlaps" % (module, success))

    # contstruct output dataframe
    df_out = pd.DataFrame(out_list, columns = out_list_headers)
    df_out.to_csv(IOTools.openFile(outfile, "w"), index=False, sep="\t")


@transform(findProximalModuleGOGenes, regex(".+/(.+).tsv.gz"), r"\1.load")
def loadProximalModuleGOGenes(infile, outfile):
    P.load(infile, outfile, options="--add-index=lnc_id")



# # ## Not Yet Implemented ##
# # The idea is to find out whether highly correlated lncRNAs are significantly more
# # likely to appear in the same module as their nearest neighbour protein coding gene. 

# # def findeRNAnnGeneModuleOverlap( infile, outfile ):
# #     """
# #     """
# # # get boolean table... iterate through modules... select lncRNAs correlated with each module...
# # # retrieve nearest neighbour genes for each lncRNA in module. 
# # # select lncRNAs that are multi-exon, single-exon, eRNA, pRNA.
# # # check whether nearest neighbour gene is in module or not. 
# # # for each module return number of correlated lncRNA with nearest neighbour in module and without. 

# # # lncRNA_classification (From loadLncRNAGenomicPositionSummary)
# # # lncRNA_chromatin_classification (From loadCollapsedeRNASummary)
# # # lncRNA_nearest_refcoding_gene (From findNearestProteinCodingGenes)
# # # refcoding_wgcna_high_fpkm_eigen_lncRNA_highest_pearsonCor
# # # get table with lncRNA_status, nn_protein_gene, modules all

# # # select first a table full of lncRNAs their classification, nearest protein coding gene, module assignment. 
# # """
# # SELECT 
# # A.classification AS genomic_position,
# # B.gene AS chromatin_classification,
# # C.closest_id AS nearest_protein_coding,
# # D.gene_id,D.MEbrown
# # FROM refcoding_wgcna_high_fpkm_eigen_lncRNA_highest_pearsonCor AS D
# # INNER JOIN lncRNA_classification AS A ON D.gene_id = A.gene_id
# # INNER JOIN lncRNA_chromatin_classification AS B ON D.gene_id = B.gene_id
# # INNER JOIN lncRNA_nearest_refcoding_gene AS C ON D.gene_id = C.gene_id
# # """

# # # select second a table with protein coding gene module assignment 
# # """
# # SELECT * FROM refcoding_wgcna_high_fpkm_module_assignment
# # """

# # # Take two dataframes
# # # subset dataframe 1 to find lncRNAs if interest (intergenic, eRNA, pRNA, me)
# # # iterate through modules (columns) select refcoding gene for lncRNA in module...
# # # check whether refcoding gene in module in second dataframe. 

# #################################################################################
# ## subsection: plot dendrograms of wgcna modules and correlated lncRNAs
# #################################################################################

# @follows( mkdir( "./expression_wgcna_filtered_fpkm/module_dendrograms" ), 
#           mkdir( "./expression_wgcna_filtered_high_cv/module_dendrograms" ),
#           loadLncRNAEigengeneCorrelation, 
#           loadRefcodingCommonNames )
# @collate( [ calculateLncRNAEigengeneCorrelation, extractModuleAssignment ],
#           regex( "(.+)/(.+)(_eigen_lncRNA_pearsonCor.tsv.gz|_module_assignment.tsv.gz)" ),
#           add_inputs( loadPerSampleFPKMs ), 
#           r"\1/module_dendrograms/\2.sentinel" )
# @jobs_limit( 1 )
# def plotModuleDendrograms( infiles, outfile ):
#     """
#     Iterate through modules extract module genes, and most highly correlated
#     lncRNAs plot PCAs of relatedness between genes.
#     """
#     # specify infiles... .load files weren't used to preserve different analyses
#     # collate pairs the add_input with each of the infiles in a pairwise manner
#     in_correlations, in_modules = [ x[0] for x in infiles ]
#     in_fpkms = infiles[0][1]
#     in_correlations = P.snip( os.path.basename( in_correlations ), ".tsv.gz" )
#     in_modules = P.snip( os.path.basename( in_modules ), ".tsv.gz" )
#     in_fpkms = P.snip( os.path.basename( in_fpkms ), ".load" )
#     # check files are assigned correctly
#     assert re.search( "pearsonCor", in_correlations ), "Incorrect file assigned %s" % in_correlations
#     assert re.search( "module_assignment", in_modules ), "Incorrect file assigned %s" % in_modules

#     # get the number of lncRNA to plot on each PCA (pass as string)
#     n_pca = PARAMS["wgcna_n_pca"]

#     # generate dataframe of fpkm values for all genes, removing outlier samples
#     statement = ( "SELECT a.*, b.gene_name"
#                   " FROM %(in_fpkms)s as a"
#                   " INNER JOIN lncRNA_refcoding_gene_names as b"
#                   " ON a.gene_id = b.gene_id" % locals() )
#     df_fpkm = PU.fetch_DataFrame( statement )
#     df_fpkm.set_index( "gene_id", inplace=True )
#     if PARAMS["wgcna_to_remove"]:
#         for sample_id in PARAMS["wgcna_to_remove"].split(","):
#             # csvdb has sample_ids separated by underscores
#             sample_id = re.sub( "\.", "_", sample_id )
#             df_fpkm.drop( sample_id, axis=1, inplace=True )

#     # create a list of modules from in_modules
#     statement = ( "SELECT * FROM %(in_modules)s LIMIT 5" % locals() )
#     df = PU.fetch_DataFrame( statement )
#     modules = list( df.columns.values )
#     modules.remove( "gene_id" )

#     # iterate through modules, plotting PCA for each one
#     # WARNING: Table of module assignments has no ME prefix
#     for module in modules:
#         # get list of top 10 most highly correlated lncRNAs
#         statement = ( "SELECT gene_id,ME%(module)s"
#                       " FROM %(in_correlations)s"
#                       " ORDER BY ABS(ME%(module)s) DESC"
#                       " LIMIT %(n_pca)s" % locals() )
#         # PU_fetch returns nested list of unicode
#         result = PU.fetch( statement )
#         lnc_genes = [ str(x[0]) for x in result ]

#         # create a list of +vely and -vely correlated lncRNA_ids
#         pos_cor_lnc = []
#         neg_cor_lnc = []
#         for lnc in result:
#             if float(lnc[1]) <= 0:
#                 neg_cor_lnc.append( str(lnc[0]) )
#             else:
#                 pos_cor_lnc.append( str(lnc[0]) )

#         # get list of all gene_ids in module
#         statement = ( "SELECT gene_id"
#                       " FROM %(in_modules)s"
#                       " WHERE %(module)s = 1" % locals() )
#         module_genes = [ str(x[0]) for x in PU.fetch( statement ) ]

#         # merge module gene_ids and highly correlated lncRNAs
#         all_gene_ids = lnc_genes + module_genes

#         # subset the expression dataframe based on merged gene list
#         df_sub = df_fpkm.ix[ all_gene_ids ]
#         # drop the index based on gene_ids and replace with gene names
#         df_sub.reset_index( drop=True, inplace=True )
#         df_sub.set_index( "gene_name", inplace=True )    

#         # for now, write the subset dataframes to a flatfile
#         outf_df = P.snip( outfile, ".sentinel" ) + "_" + module + ".tsv"
#         df_sub.to_csv( outf_df, sep="\t", na_rep = "NULL" )

#         # statement = ( "sed -i 's/NULL/NA/g' %s" % outf_df )
#         # to_cluster = False
#         # P.run()

#         P10.wgcnaPlotModuleDendrograms( df_sub, pos_cor_lnc, neg_cor_lnc, outf_df, module )
#         E.info( "Successfully plotted module: %s" % module )
#     P.touch( outfile )


# #################################################################################
# ## subsection: run GO analysis for modules
# #################################################################################
# @follows( mkdir( "./expression_wgcna_filtered_fpkm/transfac_analysis" ),
#           mkdir( "expression_wgcna_filtered_high_cv/transfac_analysis" ) )
# @split( extractModuleAssignment, 
#         regex( "(.+)/(.+)_(.+)/(.+)_module_assignment.tsv.gz" ), 
#         r"\1/\2_\3/transfac_analysis/*.foreground.tsv" )        
# def splitWGCNAModules_transfac( infile, outfiles ):
#     """
#     Split module table into separate gene lists for transfacmatch pipeline.
#     """
#     module_table = P.snip( os.path.basename( infile ), ".tsv.gz" )

#     # fetch module names from table 
#     statement = ("SELECT * FROM %(module_table)s LIMIT 5" % locals() )
#     df = PU.fetch_DataFrame( statement )
#     df.drop("gene_id", 1, inplace=True)
#     module_ids = list( df.columns.values )   

#     # generate background gene list
#     statement = ( "SELECT gene_id FROM %(module_table)s" % locals() )
#     bg_genes = [ str(x[0]) for x in PU.fetch( statement ) ]

#     # iterate through modules and write out a gene_list for each. 
#     for module in module_ids:
#         statement = ( "SELECT gene_id"
#                       " FROM %(module_table)s"
#                       " WHERE %(module)s = 1" % locals() )
#         module_genes = [ str(x[0]) for x in PU.fetch( statement ) ]

#         # write to outfile
#         out_name = module + ".foreground.tsv"
#         out_dir = os.path.dirname( infile )
#         outf = os.path.join( out_dir, "transfac_analysis", out_name )
#         outf = IOTools.openFile( outf, "w" )
#         outf.write( "gene_id\n%s\n" % ( "\n".join( sorted( module_genes ) ) ) )
#         outf.close()

#         # write all genes in WGCNA analysis out to background file
#         outf_bg = re.sub( "foreground", "background", outf.name )
#         outf_bg = IOTools.openFile( outf_bg, "w" )
#         outf_bg.write( "gene_id\n%s\n" % ( "\n".join( sorted( bg_genes ) ) ) )
#         outf_bg.close()

# #################################################################################
# ## subsection: run GO analysis for modules
# #################################################################################
# # loadModuleAssignment regex(wgcna_high_fpkm_module_assignment) # gives the gene_lists
# # add_inputs dropLowFPKMRefcodingFPKMFile # gives the background list
# # PARAMSANNOTATIONS["interface_go"] # gives the go assignments.  
# @follows( mkdir( "./expression_wgcna_filtered_fpkm/go_analysis" ), 
#           mkdir( "expression_wgcna_filtered_high_cv/go_analysis" ),
#           loadModuleAssignment )
# @split( extractModuleAssignment, 
#         regex( "(.+)/(.+)_(.+)/(.+)_module_assignment.tsv.gz" ), 
#         r"\1/\2_\3/go_analysis/*_foreground_genes.tsv" )        
# def splitWGCNAModules( infile, outfiles ):
#     """
#     Split module table into separate gene lists... allows modules to be run in
#     parallel
#     """
#     module_table = P.snip( os.path.basename( infile ), ".tsv.gz" )

#     # fetch module names from table 
#     statement = ("SELECT * FROM %(module_table)s LIMIT 5" % locals() )
#     df = PU.fetch_DataFrame( statement )
#     df.drop("gene_id", 1, inplace=True)
#     module_ids = list( df.columns.values )
    
#     # iterate through modules and write out a gene_list for each. 
#     for module in module_ids:
#         statement = ( "SELECT gene_id"
#                       " FROM %(module_table)s"
#                       " WHERE %(module)s = 1" % locals() )
#         module_genes = [ str(x[0]) for x in PU.fetch( statement ) ]

#         # write to outfile
#         out_name = module + "_foreground_genes.tsv"
#         out_dir = os.path.dirname( infile )
#         outf = os.path.join( out_dir, "go_analysis", out_name )
#         outf = IOTools.openFile( outf, "w" )
#         outf.write( "gene_id\n%s\n" % ( "\n".join( sorted( module_genes ) ) ) )
#         outf.close()


# @transform( splitWGCNAModules, 
#             regex( "(.+)/(.+)_(.+)/go_analysis/(.+)_foreground_genes.tsv" ), 
#             add_inputs( r"\1/\2_\3/refcoding_wgcna_high_\3.tsv.gz" ),            
#             # pick one... cell_location is output last
#             r"\1/\2_\3/go_analysis/\4.cell_location.results" )
# def runGOForWGCNAModules( infiles, outfile ):
#     """
#     Run GO.py for each wgcna module, using all genes in cluster as background
#     """
#     in_module, in_background = infiles
#     out_dir = os.path.dirname( outfile )
#     module = P.snip( os.path.basename( outfile ), ".cell_location.results" )

#     # generate background gene_list
#     tmp_bg = P.getTempFilename( "." )
#     E.info( "Writing Background to %s" % tmp_bg )
#     to_cluster = False
#     statement = ( " zcat %(in_background)s |"
#                   " cut -f1"
#                   " >> %(tmp_bg)s" % locals() )
#     print statement
#     P.run()

#     # get gene ontology information from pipeline annotations
#     go_assignments = os.path.join( PARAMS["annotations_dir"], 
#                                    PARAMS_ANNOTATIONS["interface_go"] )
#     go_ontology = os.path.join( PARAMS["annotations_dir"],
#                                 PARAMS_ANNOTATIONS["interface_go_obo"] )

#     # run GO
#     E.info( "Running GO enrichment analysis for module: %s" % module )
#     to_cluster = True
#     statement = ( "python %(scriptsdir)s/runGO.py"
#                   " --filename-input=%(go_assignments)s"
#                   " --genes-tsv-file=%(in_module)s"
#                   " --background-tsv-file=%(tmp_bg)s"
#                   " --method=sample --sample-size=1000"
#                   " --fdr"
#                   " --log=%(in_module)s.log"
#                   " --filename-ontology=%(go_ontology)s"
#                   " --output-filename-pattern='%(out_dir)s/%(module)s.%%(go)s.%%(section)s'" )
#     P.run()

#     E.info( "Completed GO enrichment analysis for %s" % outfile )

#     os.unlink( tmp_bg )


# # @split( extractModuleAssignment, 
# #         regex( "(.+)/(.+)_(.+)/(.+)_module_assignment.tsv.gz" ), 
# #         add_inputs( r"\1/\2_\3/refcoding_wgcna_high_\3.tsv.gz" ),
# #         r"\1/\2_\3/go_analysis/*.results" )
# # def runGOForWGCNAModules( infiles, outfiles ):
# #     in_modules, in_background = infiles
# #     out_dir = os.path.join( os.path.dirname( in_background ), "go_analysis" )

# #     # generate background gene_list
# #     tmp_bg = P.getTempFilename( "." )
# #     E.info( "Writing Background to %s" % tmp_bg )
# #     statement = ( "zcat %(in_background)s |"
# #                   " sed 1d |"
# #                   " cut -f1"
# #                   " > %(tmp_bg)s" )
# #     P.run()

# #     # doesn't work as runGO.py doesn't recognise %(set)s
# #     # statement = ( "python %(scriptsdir)s/runGO.py"
# #     #               " --filename-input=%(go_assignments)s"
# #     #               " --genes-tsv-file=%(in_modules)s"
# #     #               " --background-tsv-file=%(tmp_bg)s"
# #     #               " --method=sample --sample-size=10000"
# #     #               " --fdr"
# #     #               " --filename-ontology=%(go_ontology)s"
# #     #               " --output-filename-pattern='go_analysis/%%(set)s.%%(go)s.%%(section)s'" )
# #     # P.run()

# #     module_table = P.snip( os.path.basename( in_modules ), ".tsv.gz" )
# #     # get modules...
# #     statement = ("SELECT * FROM %(module_table)s LIMIT 5" % locals() )
# #     df = PU.fetch_DataFrame( statement )
# #     df.drop("gene_id", 1, inplace=True)
# #     module_ids = list( df.columns.values )

# #     # iterate through modules and run go analysis for each in turn
# #     for module in module_ids:
# #         E.info( "Running GO enrichment analysis for module: %s" % module )
# #         # get tempfile
# #         tmpf = P.getTempFile( "/ifs/scratch" )
# #         statement = ( "SELECT gene_id FROM %(module_table)s"
# #                       " WHERE %(module)s = 1" % locals() )
# #         module_genes = [ str(x[0]) for x in PU.fetch( statement ) ]
# #         for gene in module_genes:
# #             tmpf.write( gene + "\n" )
# #         tmpf.close()
# #         tmpf = tmpf.name

# #         statement = ( "python %(scriptsdir)s/runGO.py"
# #                       " --filename-input=%(go_assignments)s"
# #                       " --genes-tsv-file=%(tmpf)s"
# #                       # " --background-tsv-file=%(tmp_bg)s"
# #                       " --method=sample --sample-size=10000"
# #                       " --fdr"
# #                       " --filename-ontology=%(go_ontology)s"
# #                       " --output-filename-pattern='%(out_dir)s/%(module)s.%%(go)s.%%(section)s'" )
# #         P.run()
# #         os.unlink( tmpf )
# #         E.info( "Completed GO enrichment analysis for module: %s" % module )

# #     # os.unlink( tmp_bg )

# #################################################################################
# ## subsection: run GAT for TFBS overlap with module gene promoters 
# #################################################################################
# # Extract refcoding regions (1kb upstr & 1kb dstr of TSS)
# # get module assignment and output module as 4th column of bed file
# # run GAT for each TF ChIP dataset using the bedfile of modules. 
# @follows( mkdir( "characterize_wgcna_module_TFBS_overlap_promoters" ) )
# @transform( os.path.join( PARAMS["annotations_dir"], 
#                           PARAMS_ANNOTATIONS["interface_tss_gene_bed"] ),
#             regex( "(.+)/(.+).bed.gz" ), 
#             add_inputs( os.path.join( PARAMS[ "annotations_dir" ], 
#                                       PARAMS_ANNOTATIONS[ "interface_contigs" ] ) ),
#             r"characterize_wgcna_module_TFBS_overlap_promoters/refcoding_tss.bed.gz" )
# def slopRefcodingTSS( infiles, outfile ):
#     """
#     Create interval of regions around refcoding TSS... 1kb upstr, 1kb dstr
#     """
#     tss_file, contig_file = infiles
#     statement = ( "bedtools slop"
#                   " -l 1000"
#                   " -r 1000"
#                   " -s"
#                   " -i %(tss_file)s"
#                   " -g %(contig_file)s |"
#                   " gzip > %(outfile)s" )
#     P.run()


# @transform( loadModuleAssignment, 
#             regex( "(.*)/refcoding_wgcna_high_fpkm_module_assignment.load" ),
#             add_inputs( slopRefcodingTSS ), 
#             r"characterize_wgcna_module_TFBS_overlap_promoters/module_promotors.bed.gz" )
# def assignPromotersToModules( infiles, outfile ):
#     """
#     Get table of module assignment (high fpkm only) and refcoding tss, output 
#     a bed file with module assignment as 4th column
#     """
#     module_table, refcoding_tss_bed = infiles

#     P10.assignIntervalsToModules( module_table, refcoding_tss_bed, outfile )


# @transform( os.path.join(PARAMS["location_tfbs_files"], "*.bed.gz"),
#             regex( "(?:.+)/(.+).mm10.bed.gz" ),
#             add_inputs( assignPromotersToModules, 
#                         os.path.join(PARAMS["annotations_dir"], 
#                                      PARAMS_ANNOTATIONS["interface_contigs_ungapped_bed"]),
#                         os.path.join(PARAMS["annotations_dir"], 
#                                      PARAMS_ANNOTATIONS["interface_gc_profile_bed"]) ),
#             r"characterize_wgcna_module_TFBS_overlap_promoters/\1/\1.genomic_association.tsv.gz" )
# def runGATForTFBSInPromoters( infiles, outfile ):
#     """
#     Run GAT for the different TF ChIP sets against the promotor regions for the 
#     different wgcna modules
#     """
#     in_tf_chip, in_promotors, contig_file, isochore_file = infiles

#     # out directories are not made explicitly because there are too many
#     out_dir = os.path.dirname( outfile )
#     if not os.path.exists( out_dir ):
#         os.makedirs( out_dir )

#     jobOptions = "-l mem_free=10G"
#     statement = ( "gat-run.py"
#                   " --segments=%(in_tf_chip)s"
#                   " --annotations=%(in_promotors)s"
#                   " --workspace-bed-file=%(contig_file)s"
#                   " --isochore-file=%(isochore_file)s"
#                   " --ignore-segment-tracks"
#                   " --truncate-workspace-to-annotations"
#                   " --num-samples=10000"
#                   " -v5"
#                   " --log=%(outfile)s.log |"
#                   " gzip > %(outfile)s" )
#     P.run()


# #################################################################################
# ## subsection: run GAT for TFBS overlap with module gene bodies 
# #################################################################################

#################################################################################
#### METASECTION #### Search for miRNA recognition elements ####
#################################################################################
#################################################################################
# Section: Generate TargetScan Input Files
#################################################################################
@follows(mkdir("characterize_miRNA"))
@transform(classifyMergedLncRNAs, 
           regex(".+/(.+)_final.gtf.gz"),
           r"characterize_miRNA/\1.fasta.gz")
def getLncRNAFasta(infile, outfile):
    """
    Get Fasta file of collapsed lncRNA transcripts
    """
    genome = os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fasta")
    statement = ("zcat %(infile)s |"
                 " python %(scriptsdir)s/gtf2gtf.py"
                 "  --method=merge-exons"
                 "  --log=%(outfile)s.log |"
                 " python %(scriptsdir)s/gtf2gtf.py"
                 "  --method=set-transcript-to-gene"
                 "  --log=%(outfile)s.log |"
                 " python %(scriptsdir)s/gff2fasta.py"
                 "  --genome-file=%(genome)s"
                 "  --log=%(outfile)s.log"
                 "  --is-gtf |"
                 " gzip > %(outfile)s")
    P.run()



@merge([getLncRNAFasta, 
        os.path.join(PARAMS["location_mirna"], "ensembl_3utr.fasta.gz")],
       "characterize_miRNA/lncRNA_refcoding_utr.fasta")
def createTargetScanFasta(infiles, outfile):
    """
    Merge the 3' UTRs from biomart with collapsed lncRNA fasta files, format 
    output as required by targetscan.
    Output columns: transcript_id 10090 sequence, where 10090 is mm id number. 
    """
    lnc_fasta, pc_fasta = infiles

    outf = IOTools.openFile(outfile, "w")

    for line in FastaIterator.FastaIterator(IOTools.openFile(lnc_fasta)):
        gene_id = line.title.split()[0]
        outf.write(gene_id + "\t10090\t" + line.sequence + "\n")

    for line in FastaIterator.FastaIterator(IOTools.openFile(pc_fasta)):
        outf.write(line.title + "\t10090\t" + line.sequence + "\n")
    
    outf.close()


@follows(mkdir("characterize_miRNA"))
@transform(os.path.join(PARAMS["location_mirna"], "miR_Family_Info.txt"),
           regex(".+/(.+)_Family_Info.txt"),
           r"characterize_miRNA/\1_table.txt")
def createmiRNATable(infile, outfile):
    """
    Parse the table of miRNAs downloaded from TargetScan, output only mm miRNAs
    (10090), output only columns required for TS input file.
    """
    statement = ("cat %(infile)s |"
                 " awk '$3 == 10090' |"
                 " cut -f1,2,3"
                 " > %(outfile)s")
    P.run()


@merge([createTargetScanFasta, createmiRNATable],
       "characterize_miRNA/lncRNA_refcoding_predictedTargets.txt")
def runTargetScan(infiles, outfile):
    """
    Run TargetScan using only mouse miRNAs, against collapsed lncRNA transcripts
    and refcoding 3'UTRs
    """
    UTR_file, miRNA_file = infiles
    executable = "/ifs/home/jethro/apps/targetscan/targetscan_60.pl"
    statement = ("perl %(executable)s"
                 " %(miRNA_file)s"
                 " %(UTR_file)s"
                 " %(outfile)s")
    P.run()


@split(runTargetScan, 
       regex("(.+)/(lncRNA)_(refcoding)_(predictedTargets).txt"), 
       [r"\1/\2_\4_summary.p", r"\1/\3_\4_summary.p"])
def summarizeTargetScan(infile, outfiles):
    """
    Parse targetscan output and collate the miRNA_family_ID mapping to each
    locus.
    """
    out_lnc = collections.defaultdict(list)
    out_pc = collections.defaultdict(list)

    for line in IOTools.openFile(infile):
        if line.startswith("a_Gene_ID"): 
            continue
        elif line.startswith("ENSMUSG"):
            gene_id, miRNA_id = line.split()[:2]
            gene_id = gene_id.split("|")[0]
            out_pc[gene_id].append(miRNA_id)
        elif line.startswith("LNCG"):
            gene_id, miRNA_id = line.split()[:2]
            out_lnc[gene_id].append(miRNA_id)
        else:
            raise ValueError("Unrecognized gene_id: %s" % line.split()[0]) 

    pickle.dump(out_lnc, open(outfiles[0], "wb"))
    pickle.dump(out_pc, open(outfiles[1], "wb"))


@follows(mkdir("characterize_miRNA/lncRNA"))
@subdivide(summarizeTargetScan,
           regex("(.+)/lncRNA_predictedTargets_summary.p"),
           r"\1/lncRNA/*.tsv")
def splitTargetScanLncRNAs(infile, outfiles):
    """
    Write out flatfile for each lncRNA for embarassingly parallel processing. 
    """
    inf_dict = pickle.load(open(infile, "rb"))
    for gene in inf_dict.iterkeys(): 
        outf = os.path.join("characterize_miRNA/lncRNA", gene + ".tsv")
        outf = IOTools.openFile(outf, "w")
        outf.write("\t".join(inf_dict[gene]))
        outf.close()


@follows(mkdir("characterize_miRNA/lncRNA_refcoding"))
@transform("characterize_miRNA/lncRNA/*.tsv", 
           regex("(.+)/(.+).tsv"), 
           add_inputs("characterize_miRNA/refcoding_predictedTargets_summary.p"),
           r"\1/lncRNA_refcoding/\2_refcoding.tsv")
def calculateRefcodingLncRNAmiRNAOverlap(infiles, outfile):
    """
    For each lncRNA, write out overlap with each refcoding gene. 
    """
    P10.calcLncRNARefcodingOverlap(infiles, outfile, submit=True)


@merge("characterize_miRNA/lncRNA_refcoding/*_refcoding.tsv",
       "characterize_miRNA/lncRNA_refcoding_miRNA_overlap.tsv.gz")
def collateRefcodingLncRNAmiRNAOverlap(infiles, outfile):
    """
    Join lncRNA refcoding miRNA overlap files into a single table
    """
    tables = " ".join(infiles)
    job_options = "-l mem_free=20G"
    statement = ("python %(scriptsdir)s/combine_tables.py"
                 " --columns=1"
                 " --log=%(outfile)s.log"
                 " %(tables)s |"
                 " gzip > %(outfile)s")
    P.run()


@transform(collateRefcodingLncRNAmiRNAOverlap,
           suffix(".tsv.gz"),
           add_inputs(flattenModuleAssignment),
           "_plus_modules.tsv.gz")
def selectmiRNAOverlapForWGCNAGenes(infiles, outfile):
    """
    Subset the miRNA results table so that it only contains refcoding genes
    included in wgcna analysis.
    """
    miRNA_table, module_table = infiles
    module_table = P.snip(os.path.basename(module_table), ".load")
    statement = "SELECT * FROM %(module_table)s" % locals()
    df_modules = PU.fetch_DataFrame(statement)
    df_modules.set_index("gene_id", inplace=True)

    df = pd.read_table(miRNA_table, header=0, index_col=0, compression="gzip")

    # subset df so it only contains genes that are in wgcna analysis
    df_out = df_modules.merge(df, how="inner", left_index=True, right_index=True)
    df_out.to_csv(IOTools.openFile(outfile,"w" ), sep="\t", index_label="gene_id")


@transform(selectmiRNAOverlapForWGCNAGenes,
           regex("(.+)/(.+).tsv.gz"),
           add_inputs(loadLncRNAChromatinClassifications,
                      flattenLncRNAHighestEigengeneCorrelation),
           r"\1/pRNA_miRNA_enrichment.tsv.gz")
def calculatepRNAmiRNAEnrichment(infiles, outfile):
    """
    Calculate pRNA enrichment
    """
    overlap_tab, lnc_class, lnc_module = infiles
    lnc_class = P.snip(os.path.basename(lnc_class), ".load")
    lnc_module = P.snip(os.path.basename(lnc_module), ".load")

    true_overlap, null_dist, null_mean, df = P10.getLncRNA_refcoding_miRNA_overlap(lnc_module,
                                                                                      lnc_class,
                                                                                      overlap_tab,
                                                                                      lnc_type="pRNA",
                                                                                      permutations=1000)
    # print outfiles
    out_true = P.snip(outfile, ".tsv.gz") + "_true_dist.tsv.gz"
    out_true = IOTools.openFile(out_true, "w")
    out_true.write("\t".join(true_overlap))
    out_true.close()

    out_nd = P.snip(outfile, ".tsv.gz") + "_null_dist.tsv.gz"
    out_nd = IOTools.openFile(out_nd, "w")
    out_nd.write("\t".join(null_dist))
    out_nd.close()

    out_nm = P.snip(outfile, ".tsv.gz") + "_null_mean.tsv.gz"
    out_nm = IOTools.openFile(out_nm, "w")
    out_nm.write("\t".join(null_mean))
    out_nm.close()
    
    df.to_csv(IOTools.openFile(outfile, "w"), sep="\t", index=True)


@transform(selectmiRNAOverlapForWGCNAGenes,
           regex("(.+)/(.+).tsv.gz"),
           add_inputs(loadLncRNAChromatinClassifications,
                      flattenLncRNAHighestEigengeneCorrelation),
           r"\1/eRNA_miRNA_enrichment.tsv.gz")
def calculateeRNAmiRNAEnrichment(infiles, outfile):
    """
    Calculate eRNA enrichment
    """
    overlap_tab, lnc_class, lnc_module = infiles
    lnc_class = P.snip(os.path.basename(lnc_class), ".load")
    lnc_module = P.snip(os.path.basename(lnc_module), ".load")

    true_overlap, null_dist, null_mean, df = P10.getLncRNA_refcoding_miRNA_overlap(lnc_module,
                                                                                      lnc_class,
                                                                                      overlap_tab,
                                                                                      lnc_type="eRNA",
                                                                                      permutations=1000)
    # print outfiles
    out_true = P.snip(outfile, ".tsv.gz") + "_true_dist.tsv.gz"
    out_true = IOTools.openFile(out_true, "w")
    out_true.write("\t".join(true_overlap))
    out_true.close()

    out_nd = P.snip(outfile, ".tsv.gz") + "_null_dist.tsv.gz"
    out_nd = IOTools.openFile(out_nd, "w")
    out_nd.write("\t".join(null_dist))
    out_nd.close()

    out_nm = P.snip(outfile, ".tsv.gz") + "_null_mean.tsv.gz"
    out_nm = IOTools.openFile(out_nm, "w")
    out_nm.write("\t".join(null_mean))
    out_nm.close()
    
    df.to_csv(IOTools.openFile(outfile, "w"), sep="\t", index=True)

@follows(calculatepRNAmiRNAEnrichment, calculateeRNAmiRNAEnrichment)
def calculatemiRNAEnrichment():
    pass

###############################################################################
# Code that took too long to run
###############################################################################
# @follows(mkdir("characterize_miRNA/lncRNA_refcoding"))
# @product(splitTargetScanLncRNAs,
#          formatter("(.+).tsv"),
#          splitTargetScanRefcoding,
#          formatter("(.+).tsv"),
#          "{path[0][0]}/lncRNA_refcoding/{basename[0][0]}_vs_{basename[1][0]}.tsv")
# def calculateRefcodingLncRNAmiRNAOverlap(infiles, outfile):
#     """
#     For each lncRNA refcoding combination, calculate the percent of total miRNA
#     sites that are shared between the two loci.
#     """
#     lncRNA, pc_gene = infiles
#     lnc_miRNAs = IOTools.openFile(lncRNA).readline().split()
#     pc_miRNAs = IOTools.openFile(pc_gene).readline().split()


#     total = len(lnc_miRNAs) + len(pc_miRNAs)
#     # number of lncRNA miRNA sites also found in pc_gene
#     for miRNA in lnc_miRNAs:
#         if miRNA in pc_miRNAs:
#             shared += 1
#         else:
#             continue
#     # number of pc miRNA sites also found in lncRNA
#     for miRNA in pc_miRNAs:
#         if miRNA in lnc_miRNAs:
#             shared += 1
#         else:
#             continue
#     pptn = float(shared)/float(total)
#     outf = IOTools.openFile(outfile, "w")
#     outf.write(pptn)
#     outf.close()


# @merge(calculateRefcodingLncRNAmiRNAOverlap, 
#        "characterize_miRNA/lncRNA_refcoding_miRNA_overlap.tsv.gz")
# def collateRefcodingLncRNAmiRNAOverlap(infiles, outfile):
#     """
#     Create single dataframe containing lncRNA refcoding miRNA overlap
#     """
#     lnc_set = set()
#     pc_set = set()
#     for infile in infiles:
#         infile = P.snip(os.path.basename(infile), ".tsv")
#         lncRNA, refcoding = infile.split("_vs_")
#         lnc_set.add(lncRNA)
#         pc_set.add(refcoding)
#     lncRNAs = list(lnc_set)
#     pc_genes = list(pc_set)

#     # create boolean df
#     df = pd.DataFrame(np.zeros((len(pc_genes), len(lncRNAs))))
#     df.columns = lncRNAs
#     df.index = pc_genes

#     # populate df based on the proportion of miRNA sites that overlap
#     for infile in infiles:
#         pptn = float(IOTools.openFile(infile).readline())
#         infile = P.snip(os.path.basename(infile), ".tsv")
#         lncRNA, refcoding = infile.split("_vs_")
#         df.loc[refcoding, lncRNA] = pptn

#     # output pickled dataframe as backup. 


# @collate(summarizeTargetScan,
#          regex("(.+)/(.+)_predictedTargets_summary.p"),
#          add_inputs(flattenModuleAssignment),
#          r"\1/lncRNA_refcoding_miRNA_overlap.tsv.gz")
# def calculateLncRNARefcodingmiRNAOverlap(infiles, outfile):
#     """
#     """
#     # Hack - I can never get @collate output to work as input
#     module_table = infiles[0][1]
#     infiles = [x[0] for x in infiles]
#     for inf in infiles:
#         if re.search("lncRNA", os.path.basename(inf)):
#             lnc_dict =inf
#         elif re.search("refcoding", os.path.basename(inf)):
#             pc_dict = inf
#         else:
#             raise IOError("Unrecognized infile %s" % inf)

#     lnc_dict = pickle.load(open(lnc_dict, "rb"))
#     pc_dict = pickle.load(open(pc_dict, "rb"))

#     # create boolean df
#     df = pd.DataFrame(np.zeros((len(pc_dict.keys()), len(lnc_dict.keys()))))
#     df.columns = lnc_dict.keys()
#     df.index = pc_dict.keys()

#     # populate df based on the proportion of miRNA sites that overlap
#     for lncRNA in lnc_dict.keys():
#         for pc_gene in pc_dict.keys():
#             shared = 0
#             # get lists of miRNA recognition sites
#             lnc_miRNAs = lnc_dict[lncRNA]
#             pc_miRNAs = pc_dict[pc_gene]
#             # total number of miRNA sites across both loci
#             total = len(lnc_miRNAs) + len(pc_miRNAs)
#             # number of lncRNA miRNA sites also found in pc_gene
#             for miRNA in lnc_miRNAs:
#                 if miRNA in pc_miRNAs:
#                     shared += 1
#                 else:
#                     continue
#             # number of pc miRNA sites also found in lncRNA
#             for miRNA in pc_miRNAs:
#                 if miRNA in lnc_miRNAs:
#                     shared += 1
#                 else:
#                     continue
#             # the pptn of total miRNA sites shared between lncRNA and pc gene
#             pptn = float(shared)/float(total)
#             df.loc[pc_gene, lncRNA] = pptn

#     # output pickled overlap table
#     backup = P.snip(outfile, ".tsv.gz") + "_backup.p"
#     pickle.dump(df, open(backup, "wb"))

#     # df = pickle.load(open("characterize_miRNA/lncRNA_refcoding_miRNA_overlap_backup.p", "rb"))

#     # protein coding genes with sufficient coverage to be in wgcna analysis
#     module_table = P.snip(os.path.basename(module_table), ".load")
#     statement = "SELECT * FROM %(module_table)s" % locals()
#     df_modules = PU.fetch_DataFrame(statement)
#     df_modules.set_index("gene_id", inplace=True)

#     # subset df so it only contains genes that are in wgcna analysis
#     df_out = df_modules.merge(df, how="inner", left_index=True, right_index=True)
    
#     df_out.to_csv(IOTools.openFile(outfile,"w" ), sep="\t", index_label="gene_id")


# @transform(calculateLncRNARefcodingmiRNAOverlap,
#            suffix(".tsv.gz"),
#            add_inputs(flattenLncRNAHighestEigengeneCorrelation,
#                       loadLncRNAChromatinClassifications),
#            "_null.p")
# def calculateNullLncRNARefcodingmiRNAOverlap(infiles, outfile):
#     """
#     """
#     overlap_tab, lnc_module, lnc_class = infiles
#     lnc_module = P.snip(os.path.basename(lnc_module), ".load")
#     lnc_class = P.snip(os.path.basename(lnc_class), ".load")

#     # get pRNA gene_ids and most highly correlated module...
#     statement = ("SELECT a.*, b.class_empirical AS classification"
#                  " FROM %(lnc_module)s AS a"
#                  " INNER JOIN %(lnc_class)s AS b"
#                  " ON a.lnc_id = b.gene_id"
#                  " WHERE classification = 'pRNA'" % locals())
#     df_lnc = PU.fetch_DataFrame(statement)
#     lncs = map(str, df_lnc["lnc_id"].tolist())
#     # lncs = lncs + ["module",]

#     modules = map(str, [x[2:] for x in df_lnc["module"].tolist()])
#     lnc_modules = zip(lncs, modules)

#     # get boolean table of lncRNA/pc gene miRNA overlap
#     df = pd.read_table(overlap_tab, index_col = 0, compression="gzip")
#     # subset pRNAs
#     df = df.transpose()
#     df = df[df.index.isin(lncs)]
#     df = df.transpose()

#     # iterate through the pRNAs and calculate pptn of module genes that share
#     # a miRNA with lncRNA.
#     true_overlap = []
#     for lnc, module in lnc_modules:
#         df_tmp = df[df["module"] == module]
#         overlap = df_tmp[lnc].tolist()
#         pptn_overlap = float(sum(overlap))/len(overlap)
#         true_overlap.append(pptn_overlap)
    
#     # calculate null distribution for module genes/lncRNAs that share miRNA
#     null_dist = []
#     # extract module assignment 
#     modules = map(str, df["module"].tolist())
#     df_null = df.copy()
#     for i in range(0, 1000):
#         # permute protein coding module assignment
#         null_dist = np.random.permutation(df.index)
#         print null_dist

#         print type(null_dist)

#         df_null.iloc[null_dist,:]
#         df_null["module"] = modules
        
#         # iterate through pRNAs and calc pptn of module genes with shared miRNA
#         null_overlap = []
#         for lnc, module in lnc_modules:
#             df_tmp = df_null[df_null["module"] == module]
#             overlap = df_tmp[lnc].tolist()
#             pptn_overlap = float(sum(overlap))/len(overlap)
#             null_overlap.append(pptn_overlap)
#         null_dist.append(null_overlap)

#     out_dict = {"true_overlap": true_overlap, "null_dist": null_dist}
#     pickle.dump(out_dict, open(outfile, "wb"))        


#################################################################################
#### METASECTION #### Consider Transcription Factor Knock-Outs ####
#################################################################################
# Analyse transcriptome data for TF knock-out studies
#################################################################################
#################################################################################
# Section: PAX5 KO differentially expressed genes
#################################################################################
@follows( mkdir( "expression_pax5_knock_out" ) )
@jobs_limit( 1 )
@transform( os.path.join( PARAMS["location_pax5_files"], "design*.tsv" ), 
            regex( "(.+)/design(.+).tsv" ), 
            add_inputs( os.path.join( PARAMS["location_pax5_files"],
                                      "lncRNA.feature_counts.tsv.gz" ),
                          mapTranscriptToGeneID ), 
            r"expression_pax5_knock_out/\2_count_table.tsv.gz" )
def filterPax5CountTable( infiles, outfile ):
    """
    Removes antisense and antisense_downstream transcripts (data are not stranded).
    Runs the standard filtering of count tables employed in the differential 
    expression pipeline (loadTagData, filterTagData).
    Outputs a single count data table ready for differential expression 
    testing. 
    """
    P10.filterCountTable( infiles[1], infiles[0], infiles[2], outfile )


@jobs_limit( 1 )
@transform( filterPax5CountTable, 
            regex( "(.+)/(.+)_count_table.tsv.gz" ),
            add_inputs( os.path.join( PARAMS["location_pax5_files"], 
                                      r"design\2.tsv" ) ),
            r"\1/\2_pax5_ko_results.tsv" )
def runDESeq2OnPax5Knockouts( infiles, outfile ):
    """
    Runs R script for barebones DESeq2 analysis
    """
    count_table = os.path.abspath( infiles[0] )
    design_file = os.path.abspath( infiles[1] )
    outfile_stub = P.snip( os.path.abspath( outfile ), "_results.tsv" )
    statement = ( "Rscript /ifs/devel/projects/proj010/pax5_deseq2.R"
                  " %(count_table)s"
                  " %(design_file)s"
                  " %(outfile_stub)s" % locals() )
    P.run()


@merge( runDESeq2OnPax5Knockouts, 
        "pax5_knockout_deseq2_results.load" )
def loadDESeq2Pax5Results( infiles, outfile ):
    """
    Expects two infiles (one for mature one for pro), calculates BH fdr across
    both datasets. 
    """
    padj_results = P10.padjDESeq2Results( infiles )
    P.load( padj_results, outfile)
        

#################################################################################
## subsection: PAX5 TF lncRNA intersection
#################################################################################
@follows( mkdir( "intersect_lncrna_pax5" ) )
@merge( [ os.path.join(PARAMS["location_pax5_files"], 
                         "pro-pax5-R1_VS_none_peaks.narrowPeak"),
          os.path.join(PARAMS["location_pax5_files"], 
                       "pro-pax5-R2_VS_none_peaks.narrowPeak") ],
        os.path.join("intersect_lncrna_pax5", "pro_pax5.bed.gz" ) )
def extractPAX5ProTFBS( infiles, outfile ):
    """
    Get significant PAX5 pro ChIP peaks and output as bed3 format
    These samples had no input...
    """
    pro_pax5_R1, pro_pax5_R2 = infiles

    pro_pax5_R1_peaks = P.getTempFilename("/ifs/scratch")
    pro_pax5_R2_peaks = P.getTempFilename("/ifs/scratch")

    pro_pax5_R1_intersect = P.getTempFilename("/ifs/scratch")
    pro_pax5_R2_intersect = P.getTempFilename("/ifs/scratch")

    pro_pax5_R0 = P.getTempFilename("/ifs/scratch")

    # filter macs2 narrowpeak file to retrieve values with -log10(p) > 10
    # get consensus peaks from each run
    # merge consensus peaks
    statement = ( "cat %(pro_pax5_R1)s | awk '$8 > 10' > %(pro_pax5_R1_peaks)s;"
                  " cat %(pro_pax5_R2)s | awk '$8 > 10' > %(pro_pax5_R2_peaks)s;"
                  " bedtools intersect"
                  "  -a %(pro_pax5_R1_peaks)s"
                  "  -b %(pro_pax5_R2_peaks)s"
                  "  -wa "
                  "  > %(pro_pax5_R1_intersect)s;"
                  " checkpoint;"
                  " bedtools intersect"
                  "  -a %(pro_pax5_R2_peaks)s"
                  "  -b %(pro_pax5_R1_peaks)s" 
                  "  -wa "
                  "  > %(pro_pax5_R2_intersect)s;"
                  " checkpoint;"
                  " cat %(pro_pax5_R1_intersect)s > %(pro_pax5_R0)s;"
                  " cat %(pro_pax5_R2_intersect)s >> %(pro_pax5_R0)s;"
                  " bedtools sort -i %(pro_pax5_R0)s | bedtools merge -i stdin |"
                  " gzip > %(outfile)s" )
    # print statement % locals()
    P.run()

    os.unlink( pro_pax5_R1_peaks )
    os.unlink( pro_pax5_R2_peaks )
    os.unlink( pro_pax5_R1_intersect )
    os.unlink( pro_pax5_R2_intersect )
    os.unlink( pro_pax5_R0 )


@follows( mkdir( "intersect_lncrna_pax5" ) )
@merge( [ os.path.join(PARAMS["location_pax5_files"], 
                         "mat-pax5-R1_VS_mat-input-R1_peaks.narrowPeak.gz"),
          os.path.join(PARAMS["location_pax5_files"], 
                       "mat-pax5-R2_VS_mat-input-R1_peaks.narrowPeak.gz") ],
        os.path.join("intersect_lncrna_pax5", "mat_pax5.bed.gz" ) )
def extractPAX5MatureTFBS( infiles, outfile ):
    """
    Get significant PAX5 pro ChIP peaks and output as bed3 format
    These samples had a single input, one was paired end and one single end.
    """ 
    mat_pax5_R1, mat_pax5_R2 = infiles

    mat_pax5_R1_peaks = P.getTempFilename("/ifs/scratch")
    mat_pax5_R2_peaks = P.getTempFilename("/ifs/scratch")

    mat_pax5_R1_intersect = P.getTempFilename("/ifs/scratch")
    mat_pax5_R2_intersect = P.getTempFilename("/ifs/scratch")

    mat_pax5_R0 = P.getTempFilename("/ifs/scratch")

    # filter macs2 narrowpeak file to retrieve values with -log10(p) > 10
    # get consensus peaks from each run
    # merge consensus peaks
    statement = ( "zcat %(mat_pax5_R1)s | awk '$8 > 10' > %(mat_pax5_R1_peaks)s;"
                  " zcat %(mat_pax5_R2)s | awk '$8 > 10' > %(mat_pax5_R2_peaks)s;"
                  " bedtools intersect"
                  "  -a %(mat_pax5_R1_peaks)s"
                  "  -b %(mat_pax5_R2_peaks)s"
                  "  -wa "
                  "  > %(mat_pax5_R1_intersect)s;"
                  " checkpoint;"
                  " bedtools intersect"
                  "  -a %(mat_pax5_R2_peaks)s"
                  "  -b %(mat_pax5_R1_peaks)s"
                  "  -wa "
                  "  > %(mat_pax5_R2_intersect)s;"
                  " checkpoint;"
                  " cat %(mat_pax5_R1_intersect)s > %(mat_pax5_R0)s;"
                  " cat %(mat_pax5_R2_intersect)s >> %(mat_pax5_R0)s;"
                  " bedtools sort -i %(mat_pax5_R0)s | bedtools merge -i stdin |"
                  " gzip > %(outfile)s" )
    P.run()

    os.unlink( mat_pax5_R1_peaks )
    os.unlink( mat_pax5_R2_peaks )
    os.unlink( mat_pax5_R1_intersect )
    os.unlink( mat_pax5_R2_intersect )
    os.unlink( mat_pax5_R0 )

#################################################################################
## subsection: PAX5 TF Published Intervals
#################################################################################

@follows(mkdir("intersect_lncrna_pax5"))
@split("/ifs/projects/proj010/external_data/PAX5/GSE38046/inline-supplementary-material-6.xls",
       ["intersect_lncrna_pax5/mat_pax5_pub.mm9.bed.gz",
        "intersect_lncrna_pax5/pro_pax5_pub.mm9.bed.gz"])
def extractPublishedPAX5Peaks(infile, outfiles):
    """
    """
    out_mat, out_pro = outfiles
    
    # extract mature B cell peaks
    df_pro = pd.io.excel.read_excel(io=infile, 
                                   sheetname="Pax5 peaks_Pro-B cells",
                                   header=None,
                                   skiprows=[0, 1, 2, 3])
    df_pro = df_pro.ix[:,0:2]
    outf_pro = IOTools.openFile(out_pro, "w")
    df_pro.to_csv(outf_pro, header = False, index = False, sep = "\t")

    df_mat = pd.io.excel.read_excel(io=infile,
                                    sheetname="Pax5 peaks_Mat B Cells",
                                    header=None,
                                    skiprows=[0, 1, 2, 3])
    df_mat = df_mat.ix[:,0:2]
    outf_mat = IOTools.openFile(out_mat, "w")
    df_mat.to_csv(outf_mat, header = False, index = False, sep = "\t")



@transform(extractPublishedPAX5Peaks, 
           suffix(".mm9.bed.gz"), 
           add_inputs( "/ifs/mirror/ucsc/mm9/liftOver/mm9ToMm10.over.chain.gz" ),
           ".bed.gz")
def liftOverPublishedPAX5Peaks(infiles, outfile):
    """
    LiftOver published PAX5 TFBS from mm9 to mm10
    """

    infile, chainfile = infiles
    outfile = P.snip(outfile, ".gz")
    unlifted = re.sub("pub", "unlifted", outfile)

    statement = ("liftOver"
                 " %(infile)s"
                 " %(chainfile)s"
                 " %(outfile)s"
                 " %(unlifted)s;"
                 " gzip -f %(outfile)s;"
                 " gzip -f %(unlifted)s")
    P.run()


@transform(liftOverPublishedPAX5Peaks, 
           suffix(".bed.gz"),
           add_inputs(os.path.join(PARAMS["location_transcriptfiles"],
                                   "refcoding.gtf.gz"),
                      os.path.join(PARAMS["annotations_dir"],
                                   PARAMS_ANNOTATIONS["interface_contigs"])),
           "_noncoding.bed.gz")
def getNonCodingPublishedPAX5Peaks(infiles, outfile):
    """
    Use bedtools intersect to remove the TFBS that intersect a protein coding 
    gene locus.
    """
    tfbs_bed, refcoding_gtf, genome_file = infiles
    tmpf = P.getTempFilename(".")

    # create intervals from gtf gene models of protein coding genes
    statement = ("zcat %(refcoding_gtf)s |"
                  " python %(scriptsdir)s/gtf2gtf.py"
                  "  --method=merge-transcripts"
                  "  --log=%(outfile)s.log |"
                  " bedtools slop "
                  "  -i stdin"
                  "  -l 1000"
                  "  -r 0" # necessary to state -r if -l given
                  "  -g %(genome_file)s |"
                  " python %(scriptsdir)s/gff2bed.py"
                  "  --is-gtf"
                  "  --log=%(outfile)s.log"
                  " > %(tmpf)s" )
    P.run()

    # remove TFBS intervals that intersect protein coding gene intervals                 
    statement = ("bedtools intersect"
                 " -a %(tfbs_bed)s"
                 " -b %(tmpf)s"
                 " -v |"
                 " gzip > %(outfile)s")
    P.run()

    os.unlink(tmpf)


@merge(getNonCodingPublishedPAX5Peaks,
       "intersect_lncrna_pax5/pax5_pub_noncoding.bed.gz")
def mergeNonCodingPublishedPAX5Peaks(infiles, outfile):
    """
    Use Bedtools to merge the two bed files, creating a consensus set of
    intervals that represents all published PAX5 tfbs in pro and mat B cells
    that doesn't intersect a coding gene.
    """
    # cat bed files into a single tmpfile, sort tmpfile, merge tmpfile. 
    infiles = " ".join(infiles)
    to_cluster=False
    tmpf1 = P.getTempFilename("/ifs/scratch")
    tmpf2 = P.getTempFilename("/ifs/scratch")
    statement = ("zcat %(infiles)s > %(tmpf1)s;"
                 " bedtools sort -i %(tmpf1)s > %(tmpf2)s;"
                 " bedtools merge -i %(tmpf2)s | gzip > %(outfile)s")
    P.run()


@merge(liftOverPublishedPAX5Peaks,
       "intersect_lncrna_pax5/pax5_pub.bed.gz")
def mergePublishedPAX5Peaks(infiles, outfile):
    """
    Use Bedtools to merge the two bed files, creating a consensus set of
    intervals that represents all published PAX5 tfbs in pro and mat B cells
    """
    # cat bed files into a single tmpfile, sort tmpfile, merge tmpfile. 
    infiles = " ".join(infiles)
    to_cluster=False
    tmpf1 = P.getTempFilename("/ifs/scratch")
    tmpf2 = P.getTempFilename("/ifs/scratch")
    statement = ("zcat %(infiles)s > %(tmpf1)s;"
                 " bedtools sort -i %(tmpf1)s > %(tmpf2)s;"
                 " bedtools merge -i %(tmpf2)s | gzip > %(outfile)s")
    P.run()


#################################################################################
## subsection: PAX5 TF lncRNA intersection
#################################################################################

@follows( mkdir( "intersect_lncrna_pax5" ) )
@transform( [ extractPAX5ProTFBS, 
              extractPAX5MatureTFBS, 
              liftOverPublishedPAX5Peaks,
              getNonCodingPublishedPAX5Peaks ],
            regex("(.+)/(.+).bed.gz"),
            add_inputs( classifyMergedLncRNAs, # genes_gtf
                        os.path.join(PARAMS["annotations_dir"], 
                                     PARAMS_ANNOTATIONS["interface_contigs"]), # genome_file_ungapped
                        os.path.join(PARAMS["annotations_dir"], 
                                     PARAMS_ANNOTATIONS["interface_contigs_ungapped_bed"]), # genome_file
                        os.path.join(PARAMS["annotations_dir"], 
                                     PARAMS_ANNOTATIONS["interface_gc_profile_bed"]) ), # isochore_file
            r"intersect_lncrna_pax5/\2_coverage.tsv.gz" )
def getlncRNAPAX5Overlap( infiles, outfile ):
    """
    Calculate overlap statistics for lncRNAs and PAX5 binding sites
    """
    tf_bed, genes_gtf, genome_file, genome_file_ungapped, isochore_file = infiles
    outf_stub = P.snip( os.path.abspath(outfile), "_coverage.tsv.gz" )
    P10.compareOverlaps( genes_gtf, 
                         tf_bed, 
                         genome_file, 
                         genome_file_ungapped, 
                         isochore_file, 
                         outf_stub )


@collate( getlncRNAPAX5Overlap,
          regex( "(.+)/(.+)_pax5(.*)_coverage.tsv.gz" ),
          r"\1/pax5\3_lncRNA_coverage.tsv.gz" )
def collatelncRNAPAX5Overlap( infiles, outfile ):
    """
    Combine boolean lists of PAX5 coverage of lncRNAs for Pro and Mature B cells. 
    """
    infiles = " ".join( infiles )
    statement = ("python %(scriptsdir)s/combine_tables.py"
                 " --columns=1"
                 " --log=%(outfile)s.log"
                 " %(infiles)s |"
                 " gzip > %(outfile)s" )
    P.run()


@transform( collatelncRNAPAX5Overlap,
            regex( "(.+)/(.+).tsv.gz" ),
            r"\2.load" )
def loadlncRNAPAX5Overlap( infile, outfile ):
    P.load( infile, outfile, options="--add-index=gene_id" )


#################################################################################
## subsection: generate count data for leukemia B cell data
#################################################################################
@follows( mkdir( "./expression_pax5_leukemia_cells" ) )
@transform(combineRefcodingAndLncRNA,
           regex("(.*)lncRNA_refcoding.gtf.gz"),
            r"expression_pax5_leukemia_cells/lncRNA_refcoding.gtf.gz" )
def removeAntisenseLncRNAs(infile, outfile):
    """
    On the assumption that data are not stranded, antisense and antisense_dstr
    lncRNAs are removed from the gtf prior to running featureCounts.
    """
    outf = IOTools.openFile(outfile, "w")
    for gtf in GTF.iterator(IOTools.openFile(infile)):
        if gtf.source in ["antisense", "antisense_downstream"]:
            continue
        else:
            outf.write( str(gtf) + "\n" )
    outf.close()


@transform( os.path.join( PARAMS["location_bamfiles_leukemia"], "*.bam"), 
            regex( "(.+)/(.+).bam" ), 
            add_inputs( removeAntisenseLncRNAs ), 
            r"expression_pax5_leukemia_cells/\2_vs_lncRNA_refcoding.tsv.gz" )
def buildFeatureCounts_leukemia( infiles, outfile ):
    """
    Runs PipelineRnaseq.runFeatureCounts() 
    Library is unstranded (i.e. featureCounts strand=0)
    Discard multimapping reads (as suggested by SA at CSAMA) is default
    Minimum mapping quality (-Q) set to 10
    """
    bamfile, annotations = infiles
    PipelineRnaseq.runFeatureCounts(
        annotations,
        bamfile,
        outfile,
        nthreads=4,
        strand=0,
        options='-Q10' )


@collate( buildFeatureCounts_leukemia, 
          regex( "(.+)/(.+)_vs_(.+).tsv.gz" ),
          r"\1/\3_raw_counts.tsv.gz" )
def summarizeFeatureCounts_leukemia( infiles, outfile ):
    """
    Collate count data into a single file to be read as matrix of counts
    (Lifted out of pipeline_rnaseqdiffexpression.py, with minor changes)
    """
    infiles = " ".join( infiles )
    statement = ( "python %(scriptsdir)s/combine_tables.py"
                  " --columns=1"
                  " --take=7"
                  " --use-file-prefix"
                  " --regex-filename='(.+)_vs.+.tsv.gz'"
                  " --log=%(outfile)s.log"
                  " %(infiles)s"
                  " | sed 's/Geneid/gene_id/'"
                  " | gzip > %(outfile)s" )
    P.run()

@merge(buildFeatureCounts_leukemia, "expression_pax5_leukemia_cells/design.tsv")
def generateLeukemiaDesignFile(infiles, outfile):
    """
    Create design file based on the featureCounts outfile names
    """
    outf = IOTools.openFile(outfile, "w")
    outf.write("track\tinclude\tgroup\tpair\n")
    
    for inf in infiles:
        track = P.snip(os.path.basename(inf), ".star_vs_lncRNA_refcoding.tsv.gz")
        include = "1"
        group = track.split("-")[1]
        pair = "0"
        outf.write("\t".join([track, include, group, pair]) + "\n")

    outf.close()


@jobs_limit( 1, "RGlobalEnv" )
@transform( summarizeFeatureCounts_leukemia, 
            regex( "(.+)/(.+)_raw_counts.tsv.gz" ),
            add_inputs( generateLeukemiaDesignFile ),
            r"\1/\2_BAll_results.tsv" )
def runDESeq2OnLeukemiaCells( infiles, outfile ):
    """
    Runs R script for barebones DESeq2 analysis
    """
    count_table = os.path.abspath( infiles[0] )
    design_file = os.path.abspath( infiles[1] )
    outfile_stub = P.snip( os.path.abspath( outfile ), "_results.tsv" )
    statement = ( "Rscript /ifs/devel/projects/proj010/all_deseq2.R"
                  " %(count_table)s"
                  " %(design_file)s"
                  " %(outfile_stub)s" % locals() )
    P.run()


@transform( runDESeq2OnLeukemiaCells, 
            regex( ".+/(.+)_(.+)_BAll_results.tsv" ), 
            r"pax5_leukemia_deseq2_results.load" )
def loadDESeq2LeukemiaResults( infile, outfile ):
    """
    Loads the differential expression results for ebf1 KO in Mature B cells
    """
    tmpf = P.getTempFilename( "/ifs/scratch" )

    # Add column for cell type and load
    df = pd.read_table( infile, sep = ",", index_col=0 )
    df.index.name = "gene_id"
    df.to_csv( tmpf, sep="\t" )

    P.load( tmpf, outfile, options="--add-index=gene_id" )

    os.unlink( tmpf )


@follows(mkdir("expression_pax5_leukemia_cells"))
@follows(mkdir("expression_pax5_leukemia_cells/merged_bam"))
@collate( os.path.join( PARAMS["location_bamfiles_leukemia"], "*.bam"), 
            regex( "(.+)/Bcell-(.+)-R(.).star.bam" ), 
            r"expression_pax5_leukemia_cells/merged_bam/Bcell-\2-R0.bam" )
def mergeLeukemiaBamfiles(infiles, outfile):
    """
    Merge Bamfiles. 
    """
    infiles = " ".join(infiles)
    P10.mergeBam(infiles, outfile)


@collate(mergeLeukemiaBamfiles, regex("(.+)/(.+).bam"), r"\1/normalized.sentinel")
def normalizeLeukemiaBamfiles(infiles, outfile):
    for infile in infiles:
        if re.search("-dox-", infile):
            in_dox = infile
        elif re.search("-all-", infile):
            in_all = infile
        else:
            raise IOError("Unrecognized infile: %s" % infile)
    
    out_dox = P.snip(in_dox, ".bam") + "_normalized.bam"
    out_all = P.snip(in_all, ".bam") + "_normalized.bam"

    # count reads in bamfile
    dox_count, all_count = 0, 0
    for read in pysam.Samfile(in_dox, "rb"):
        dox_count += 1
    for read in pysam.Samfile(in_all, "rb"):
        all_count += 1
    E.info("There are %i dox reads" % dox_count)
    E.info("There are %i all reads" % all_count)

    # normalize larger bamfile
    if dox_count > all_count:
        E.info("Randomly downsampling dox bamfile")
        P10.normalize(in_dox, dox_count, out_dox, all_count)
        os.symlink(os.path.abspath(in_all), out_all)
        pysam.index(out_all)
    else:
        E.info("Randomly downsampling all bamfile")
        P10.normalize(in_all, all_count, out_all, dox_count)
        os.symlink(in_dox, out_dox)
        pysam.index(out_dox)


#################################################################################
## subsection: find wgcna module enrichment in B-ALL DE results
#################################################################################

@follows(mkdir("expression_pax5_leukemia_cells_wgcna_module_enrichment"))
@transform(flattenModuleAssignment,
           regex("(.+).load"),
           add_inputs(loadDESeq2LeukemiaResults),
           r"expression_pax5_leukemia_cells_wgcna_module_enrichment/module_gsea_enrichment.tsv")
def getWGCNAModuleEnrichment_BAll(infiles, outfile):
    """
    Check to see which WGCNA modules are enriched in DE results comparing 
    untreated and treated B-ALL gene expression.
    """
    modules, de_results = [P.snip(os.path.basename(x), ".load") for x in infiles]

    # fetch module assignment
    statement = ("SELECT * FROM %(modules)s" % locals())
    df_module = PU.fetch_DataFrame(statement)
    
    # fetch ABSOLUTE fold-change 
    statement = ("SELECT b.gene_id, ABS(b.log2FoldChange) AS l2f"
                 " FROM refcoding_wgcna_rlog_module_assignment_flattened AS a"
                 " LEFT JOIN pax5_leukemia_deseq2_results AS b"
                 " ON a.gene_id = b.gene_id")
    df_fc = PU.fetch_DataFrame(statement)
    df_fc.set_index("gene_id", inplace=True)
    df_fc.dropna(inplace=True)

    # push to R env. 
    df_fc = pandas2ri.py2ri(df_fc)
    df_module = pandas2ri.py2ri(df_module)

    # run GSEA enrichment using piano
    piano = importr("piano")
    myGSC = piano.loadGSC(df_module)
    gsaRes = piano.runGSA(df_fc, gsc=myGSC)

    # output GSEA results
    df_out = pandas2ri.ri2pandas(piano.GSAsummaryTable(gsaRes))
    
    print df_out.__class__
    df_out.to_csv(outfile, sep="\t", index=False)


@transform(loadlncRNAPAX5Overlap,
           regex("(.+)_pub_noncoding_lncRNA_coverage.load"),
           add_inputs(getWGCNAModuleEnrichment_BAll,
                      flattenLncRNAHighestEigengeneCorrelation),
           r"expression_pax5_leukemia_cells_wgcna_module_enrichment/enriched_module_tfbs_overlap.tsv")
def getEnrichedModuleTFBSOverlap_lncRNA(infiles, outfile):
    """
    
    """
    pro_overlap, gsea_res, lnc_modules = infiles
    pro_overlap = P.snip(os.path.basename(pro_overlap), ".load")
    lnc_modules = P.snip(os.path.basename(lnc_modules), ".load")

    # fetch the lncRNAs with pax5 binding in pro B cells
    statement = ("SELECT a.gene_id,"
                 " a.pro_pax5_pub_noncoding AS pro_pax5,"
                 " b.module"
                 " FROM %(pro_overlap)s AS a"
                 " INNER JOIN %(lnc_modules)s AS b"
                 " ON a.gene_id = b.lnc_id" % locals())
    df = PU.fetch_DataFrame(statement)
    df.set_index("gene_id", inplace=True)

    # find the wgcna modules enriched in PRE vs B-ALL comparison
    enriched_modules = []
    for line in IOTools.openFile(gsea_res):
        if line.startswith("Name"): continue
        line = line.split()
        if float(line[4]) < 0.05:
            enriched_modules.append(line[0])
    enriched_modules = ["ME" + x for x in enriched_modules]
    E.info("Enriched modules: %s" % "\t".join(enriched_modules))

    df["enriched_module"] = df["module"].isin(enriched_modules)
    df["pax5_overlap"] = df["pro_pax5"] == "TRUE"
    df.drop("module", axis=1, inplace=True)
    df.drop("pro_pax5", axis=1, inplace=True)
    df = df.astype(int)
    df.to_csv(outfile, sep="\t")

    # output fishers exact test results
    out_fisher = P.snip(outfile, ".tsv") + "_fishers.rds"
    df_r = pandas2ri.ri2py(df) # com.convert_to_r_dataframe(df)
    # df_sum = pandas2ri.ri2pandas(R["table"](df_r))
    # print df_sum
    ft = R["fisher.test"](R["table"](df_r))
    R["saveRDS"](ft, file=out_fisher)
    
    
#################################################################################
# Section: EBF1 KO differentially expressed genes
#################################################################################
@follows( mkdir( "expression_ebf1_knock_out" ) )
@jobs_limit( 1 )
@transform( os.path.join( PARAMS["location_ebf1_files"], "design*.tsv" ), 
            regex( "(.+)/design(.+).tsv" ), 
            add_inputs( os.path.join( PARAMS["location_ebf1_files"],
                                      "lncRNA.feature_counts.tsv.gz" ),
                          mapTranscriptToGeneID ), 
            r"expression_ebf1_knock_out/\2_count_table.tsv.gz" )
def filterEbf1CountTable( infiles, outfile ):
    """
    Doesn't remove any lncRNAs (data are fr-firststrand (dUTP)).
    Runs the standard filtering of count tables employed in the differential 
    expression pipeline (loadTagData, filterTagData).
    Outputs a single count data table ready for differential expression 
    testing. 
    """
    design_table, count_table, biotype_table = infiles
    P10.filterCountTable( count_table, 
                          design_table, 
                          biotype_table, 
                          outfile, 
                          biotypes_to_remove = [])


@jobs_limit( 1 )
@transform( filterEbf1CountTable, 
            regex( "(.+)/(.+)_count_table.tsv.gz" ),
            add_inputs( os.path.join( PARAMS["location_ebf1_files"], 
                                      r"design\2.tsv" ) ),
            r"\1/\2_ebf1_ko_results.tsv" )
def runDESeq2OnEbf1Knockouts( infiles, outfile ):
    """
    Runs R script for barebones DESeq2 analysis using script written for Pax5
    There are no replicates for the knockout, DESeq2 handles this automatically
    by automatically calculating dispersion across all samples rather than within
    condition. This is more conservative than a standard DESeq2 run.
    """
    count_table = os.path.abspath( infiles[0] )
    design_file = os.path.abspath( infiles[1] )
    outfile_stub = P.snip( os.path.abspath( outfile ), "_results.tsv" )
    statement = ( "Rscript /ifs/devel/projects/proj010/pax5_deseq2.R"
                  " %(count_table)s"
                  " %(design_file)s"
                  " %(outfile_stub)s" % locals() )
    P.run()


@transform( runDESeq2OnEbf1Knockouts, 
            regex( ".+/(.+)_(.+)_ko_results.tsv" ), 
            r"\2_knockout_deseq2_results.load" )
def loadDESeq2Ebf1Results( infile, outfile ):
    """
    Loads the differential expression results for ebf1 KO in Mature B cells
    """
    tmpf = P.getTempFilename( "/ifs/scratch" )

    # Add column for cell type and load
    df = pd.read_table( infile, sep = ",", index_col=0 )
    df["Bcell"] = "Mature"
    df.index.name = "gene_id"
    df.to_csv( tmpf, sep="\t" )

    P.load( tmpf, outfile, options="--add-index=gene_id" )

    os.unlink( tmpf )


#################################################################################
## subsection: EBF1 TF lncRNA intersection
#################################################################################
@follows( mkdir( "intersect_lncrna_ebf1" ) )
@transform( os.path.join(PARAMS["location_ebf1_files"], 
                         "spl-ebf1-R1_VS_spl-input-R1_peaks.narrowPeak.gz"),
            regex( "(.+)/(.+)-(.+)-R1_VS_spl-input-R1_peaks.narrowPeak.gz" ),
            r"intersect_lncrna_ebf1/\2_\3.bed.gz" )
def extractEBF1TFBS( infile, outfile ):
    """
    Get significant ChIP peaks and output as bed3 format
    """
    # filter macs2 narrowpeak file to retrieve values with -log10(p) > 10
    statement = ( "zcat %(infile)s |"
                  " awk '$8 > 10' |"
                  " cut -f1,2,3 |"
                  " gzip > %(outfile)s" )
    P.run()


@follows( mkdir( "intersect_lncrna_ebf1" ) )
@transform( extractEBF1TFBS,
            regex("(.+)/(.+).bed.gz"),
            add_inputs( classifyMergedLncRNAs, # genes_gtf
                        os.path.join(PARAMS["annotations_dir"], 
                                     PARAMS_ANNOTATIONS["interface_contigs"]), # genome_file_ungapped
                        os.path.join(PARAMS["annotations_dir"], 
                                     PARAMS_ANNOTATIONS["interface_contigs_ungapped_bed"]), # genome_file
                        os.path.join(PARAMS["annotations_dir"], 
                                     PARAMS_ANNOTATIONS["interface_gc_profile_bed"]) ), # isochore_file
            r"intersect_lncrna_ebf1/\2_coverage.tsv.gz" )
def getlncRNAEBF1Overlap( infiles, outfile ):
    """
    Calculate overlap statistics for lncRNAs and EBF1 binding sites
    """
    tf_bed, genes_gtf, genome_file, genome_file_ungapped, isochore_file = infiles
    outf_stub = P.snip( os.path.abspath(outfile), "_coverage.tsv.gz" )
    P10.compareOverlaps( genes_gtf, 
                         tf_bed, 
                         genome_file, 
                         genome_file_ungapped, 
                         isochore_file, 
                         outf_stub )


@transform( getlncRNAEBF1Overlap,
            regex( "(.+)/spl_(.+)_(.+).tsv.gz" ),
            r"\2_lncRNA_\3.load" )
def loadlncRNAEBF1Overlap( infile, outfile ):
    P.load( infile, outfile, options="--add-index=gene_id" )


@follows( loadlncRNAPAX5Overlap, loadlncRNAEBF1Overlap )
def runTFBSOverlap():
    pass


#################################################################################
#################################################################################
#### METASECTION #### Generate BigWigs for plotting  ####
#################################################################################
#################################################################################

@follows(mkdir("output_bigwigfiles"), mkdir("output_bigwigfiles/chromatin_files"))
@transform( "characterize_chromatin_state/bamfiles_normalized/*bam",
            regex( ".+/(.+)-(K4me1-R0|K4me3-R0)_normalized.bam" ),
            r"output_bigwigfiles/chromatin_files/\1-\2.bw")
def createChromatinBigWigs(infile, outfile):
    """
    Convert the merged chromatin bamfiles into bigwigs. These do not need to be 
    normalised as they are already size adjusted based on their corresponding pairs
    """
    bamfile = infile
    contig_file = os.path.join(PARAMS["annotations_dir"], "contigs.tsv")
    tmpf = P.getTempFilename("/ifs/scratch")
    statement = (" bedtools bamtobed"
                 "  -i %(bamfile)s"
                 "  -split |"
                 " bedtools genomecov"
                 "  -i stdin"
                 "  -g %(contig_file)s"
                 "  -bg"
                 "  > %(tmpf)s"
                 "  2> %(outfile)s.log;"
                 " bedGraphToBigWig"
                 "  %(tmpf)s" # in.bedGraph
                 "  %(contig_file)s" # chrom.sizes
                 "  %(outfile)s" # out.bw
                 "  2>> %(outfile)s.log")
    P.run()


@follows(mkdir("output_bigwigfiles"), mkdir("output_bigwigfiles/pooled_counts"))
@transform(os.path.join(PARAMS["location_bamfiles_filtered_merged"], "*.bam"), 
           regex("(.+)/(.+).bam"),
           add_inputs(combineRefcodingAndLncRNA),
           r"output_bigwigfiles/pooled_counts/\2_vs_lncRNA_refcoding.tsv.gz")
def buildFeatureCounts_pooled(infiles, outfile):
    """
    Runs PipelineRnaseq.runFeatureCounts() 
    Assuming library is fr-firststrand (i.e. featureCounts strand=2)
    Discard multimapping reads (as suggested by SA at CSAMA)
    Minimum mapping quality (-Q) set to 10
    """
    bamfile, annotations = infiles
    PipelineRnaseq.runFeatureCounts(
        annotations,
        bamfile,
        outfile,
        nthreads=4,
        strand=2,
        options='-Q10' )


@collate( buildFeatureCounts_pooled, 
          regex( "(.+)/(.+)_vs_(.+).tsv.gz" ),
          r"\1/\3_raw_counts.tsv.gz" )
def summarizeFeatureCounts_pooled( infiles, outfile ):
    """
    Collate count data into a single file to be read as matrix of counts
    (Lifted out of pipeline_rnaseqdiffexpression.py, with minor changes)
    """
    infiles = " ".join( infiles )
    statement = ( "python %(scriptsdir)s/combine_tables.py"
                  " --columns=1"
                  " --take=7"
                  " --use-file-prefix"
                  " --regex-filename='(.+)_vs.+.tsv.gz'"
                  " --log=%(outfile)s.log"
                  " %(infiles)s"
                  " | sed 's/Geneid/gene_id/'"
                  " | gzip > %(outfile)s" )
    P.run()


@transform(summarizeFeatureCounts_pooled,
           suffix("_raw_counts.tsv.gz"),
           "_normalization_factors.tsv")
def calcFeatureCountsNormalizationFactors_pooled(infile, outfile):
    """
    Use DESeq to calculate normalization factors based on featureCounts data
    """
    # get count dataframe
    df = pd.read_table(infile, header=0, index_col=0, compression="gzip")
    headers = list(df.columns.values)
    df = pandas2ri.py2ri(df) # com.convert_to_r_dataframe(df)

    # calculate robust normalization factors using DESeq
    DESeq = importr("DESeq")
    size_factors = DESeq.estimateSizeFactorsForMatrix(df)

    # output size factors
    df_out = pd.DataFrame({"sample_id": headers,
                           "normalization": size_factors.tolist()})
    df_out.to_csv(outfile, sep="\t", index=False)


@follows(mkdir("output_bigwigfiles/rnaseq_files"))
@transform(os.path.join(PARAMS["location_bamfiles_filtered_merged"], "*.bam"), 
           regex("(.+)/(.+).bam"),
           add_inputs(calcFeatureCountsNormalizationFactors_pooled),
           r"output_bigwigfiles/rnaseq_files/\2.bw")
def createRNAseqBigWigs(infiles, outfile):
    """
    Use bedtools bamtobed to convert bamfile to bedfile of split intervals,
    then bedtools genomeCov to calculate normalized bedgraph files, then
    kent tool's bedGraphToBigWig
    NB. Normalization is done on the split intervals... whereas estimation of
    the normalization factor was done on the paired read counts.
    """
    bamfile, size_factors = infiles
    sample_id = P.snip(os.path.basename(bamfile), ".bam")
    contig_file = os.path.join(PARAMS["annotations_dir"], "contigs.tsv")

    # get normalization factor note that this is 1/deseq size_factor
    df = pd.read_table(size_factors, header=0)
    df.set_index("sample_id", inplace=True)
    factor = df.loc[sample_id, "normalization"]
    factor = str(1.0/factor)

    # generate bigwig (bedGraphToBigWig can't be piped)
    tmpf = P.getTempFilename("/ifs/scratch")
    statement = (" bedtools bamtobed"
                 "  -i %(bamfile)s"
                 "  -split |"
                 " bedtools genomecov"
                 "  -i stdin"
                 "  -g %(contig_file)s"
                 "  -bg"
                 "  -scale %(factor)s"
                 "  > %(tmpf)s"
                 "  2> %(outfile)s.log;"
                 " bedGraphToBigWig"
                 "  %(tmpf)s" # in.bedGraph
                 "  %(contig_file)s" # chrom.sizes
                 "  %(outfile)s" # out.bw
                 "  2>> %(outfile)s.log")
    P.run()


@follows(createChromatinBigWigs, createRNAseqBigWigs)
def createBigWigs():
    pass
#################################################################################
@transform(os.path.join(PARAMS["location_bamfiles_filtered_merged"], "*.bam"), 
           suffix("R0.bam"),
#           add_inputs(calcFeatureCountsNormalizationFactors_pooled),
           add_inputs("output_bigwigfiles/pooled_counts/lncRNA_refcoding_normalization_factors.tsv"),
           "_normalized.bamx")
def normalize_pooled_rnaseq_bamfiles(infiles, outfile):
    """
    Using the normalization factors calculated above, downsample the bamfiles. 
    """
    bamfile, normalization_file = infiles
    P10.normalizeBams(bamfile, normalization_file, outfile, submit=False)


#################################################################################
#### METASECTION #### Creating Output Tables ####
#################################################################################
#################################################################################
# Section: A bed12 of lncRNAs...
#################################################################################
@follows(mkdir("output_tables"))
@transform(classifyMergedLncRNAs,
           regex("(.+)/(.+).gtf.gz"),
           r"output_tables/\2_bed.gz")
def createLncRNABed12(infile, outfile):
    """
    Create table that approximates lncRNA bed12 file, except that there's an extra
    name field containing transcript ID. 
    """
    P10.gtfToBed12(infile, outfile)


@transform(createLncRNABed12,
           regex(".+/(.+).gz"),
           r"\1.load")
def loadLncRNABed12(infile, outfile):
    """
    """
    header_names=("contig,start,end,gene_id,transcript_id,score,strand,"
                  "thickStart,thickEnd,colourRGB,blockCount,blockSizes,blockStarts")
    P.load(infile, outfile, options="--header-names=%s" % header_names)





#################################################################################
# Section: UCSC Track Hub Files
#################################################################################
# convert ChIP-seq bamfiles to bw
@follows(mkdir("output_ucsc_track_hub"),
         mkdir("output_ucsc_track_hub/chipseq_bigwigs"))
@transform(os.path.join(PARAMS["location_chipseq_merged"], "*.bwa.bam"),
           regex(".+/(.+)K4me1-R0_deduped.bwa.bam"),
           r"output_ucsc_track_hub/chipseq_bigwigs/\1H3K4me1.bw")
def convertH3K4me1BamfilesToBigWigs(infile, outfile):
    """
    Use bam2wiggle.py to convert (merged, deduped) chipseq bamfiles to bigwigs
    """
    # job_options = "-l mem_free=30G"
    # statement = ("python %(scriptsdir)s/bam2wiggle.py"
    #              " --output-format=bigwig"
    #              " --output-filename-pattern=%(outfile)s"
    #              " %(infile)s")
    # P.run()
    P10.bamToBigWig(infile, outfile, job_options="-l mem_free=35G", submit=True)


@follows(mkdir("output_ucsc_track_hub"),
         mkdir("output_ucsc_track_hub/chipseq_bigwigs"))
@transform(os.path.join(PARAMS["location_chipseq_merged"], "*.bwa.bam"),
           regex(".+/(.+)input-R0_deduped.bwa.bam"),
           r"output_ucsc_track_hub/chipseq_bigwigs/\1input.bw")
def convertInputBamfilesToBigWigs(infile, outfile):
    """
    Use bam2wiggle.py to convert (merged, deduped) chipseq bamfiles to bigwigs
    """
    P10.bamToBigWig(infile, outfile, job_options="-l mem_free=35G", submit=True)


@follows(mkdir("output_ucsc_track_hub"),
         mkdir("output_ucsc_track_hub/chipseq_bigwigs"))
@transform(os.path.join(PARAMS["location_chipseq_merged"], "*.bwa.bam"),
           regex(".+/(.+)K4-R0_deduped.bwa.bam"),
           r"output_ucsc_track_hub/chipseq_bigwigs/\1H3K4me3.bw")
def convertH3K4me3BamfilesToBigWigs(infile, outfile):
    """
    Use bam2wiggle.py to convert (merged, deduped) chipseq bamfiles to bigwigs
    """
    P10.bamToBigWig(infile, outfile, job_options="-l mem_free=35G", submit=True)


@follows(mkdir("output_ucsc_track_hub"),
         mkdir("output_ucsc_track_hub/chipseq_peak_intervals"))
@transform(cleanNarrowPeakFiles,
           regex(".+/(.+)-(K4me1|K4me3).bed.gz"),
           add_inputs(os.path.join(PARAMS["annotations_dir"], "contigs.tsv")),
           r"output_ucsc_track_hub/chipseq_peak_intervals/\1-H3\2.bb")
def convertChIPSeqPeaksToBigBed(infiles, outfile):
    """
    Use bedToBigBed to convert peak intervals to bigbed format. 
    """
    bedfile, contigfile = infiles 
    tmpf = P.getTempFilename("/ifs/scratch")
    job_options = "-l mem_free=5G"
    statement = ("zcat %(bedfile)s | sort -k1,1 -k2,2n > %(tmpf)s;"
                 " bedToBigBed %(tmpf)s %(contigfile)s %(outfile)s")
    P.run()
    os.unlink(tmpf)


@follows(mkdir("output_ucsc_track_hub"),
         mkdir("output_ucsc_track_hub/chipseq_peak_intervals"))
@collate(cleanNarrowPeakFiles,
         regex(".+/(.+)-(.+).bed.gz"),
           add_inputs(os.path.join(PARAMS["annotations_dir"], "contigs.tsv")),
         r"output_ucsc_track_hub/chipseq_peak_intervals/H3\2.bb")
def convertMergedChIPSeqPeaksToBigBed(infiles, outfile):
    """
    Merge the bed files for each chromatin mark. Convert merged file to bigBed.
    """
    contig_file = infiles[0][1]
    infiles = " ".join([x[0] for x in infiles])

    tmpf = P.getTempFilename("/ifs/scratch")
    statement = ("zcat %(infiles)s |"
                 " sort -k1,1 -k2,2n |"
                 " bedtools merge -i stdin |"
                 " sort -k1,1 -k2,2n > %(tmpf)s;"
                 "bedToBigBed %(tmpf)s %(contig_file)s %(outfile)s")
    P.run()
    os.unlink(tmpf)


@follows(mkdir("output_ucsc_track_hub"),
         mkdir("output_ucsc_track_hub/rnaseq_bigwigs"))
@transform([os.path.join(PARAMS["location_bamfiles_filtered"], "*.bam"),
            os.path.join(PARAMS["location_bamfiles_filtered_merged"], "*.bam")],
           regex(".+/Bcell-(.+)-(.+).bam"),
           r"output_ucsc_track_hub/rnaseq_bigwigs/\1-\2.bw")
def convertRNASeqBamfilesToBigWigs(infile, outfile):
    """
    Use bam2wiggle.py to convert (merged, deduped) chipseq bamfiles to bigwigs
    """
    P10.bamToBigWig(infile, outfile, job_options="-l mem_free=35G", submit=True)


@follows(mkdir("output_ucsc_track_hub"),
         mkdir("output_ucsc_track_hub/lncRNA_bigbed2"))
@split(classifyMergedLncRNAs,
       regex(".+/lncRNA_final.gtf.gz"),
       add_inputs(os.path.join(PARAMS["annotations_dir"], "contigs.tsv")),
       ["output_ucsc_track_hub/lncRNA_bigbed2/lncRNA_gene_models.bb",
        "output_ucsc_track_hub/lncRNA_bigbed2/lncRNA_transcript_models.bb"])
def convertlncRNAToBigBed(infiles, outfiles):
    gtffile, contigs = infiles
    out_gene, out_transcript = outfiles

    # create an .es file for indexing the bigbed by gene/transcript name
    out_es = P.snip(out_gene, "_gene_models.bb") + ".es"
    with IOTools.openFile(out_es, "w") as outes:
        outes.write(
            "table hg18KGchr7\n"
            "\"UCSC Genes for chr7 with color plus GeneSymbol and SwissProtID\"\n"
            "(\n"
            "string  chrom;\"Reference sequence chromosome or scaffold\"\n"
            "uint    chromStart;\"Start position of feature on chromosome\"\n"
            "uint    chromEnd;\"End position of feature on chromosome\"\n"
            "string  name;\"Name of gene\"\n"
            "uint    score;\"Score\"\n"
            "char[1] strand;\"+ or - for strand\"\n"
            "uint    thickStart;\"Coding region start\"\n"
            "uint    thickEnd;\"Coding region end\"\n"
            "uint  reserved;\"Green on + strand, Red on - strand\"\n"
            "int blockCount;     \"Number of blocks\"\n"
            "int[blockCount] blockSizes; \"Comma separated list of block sizes\"\n"
            "int[blockCount] chromStarts; \"Start positions relative to chromStart\"\n"
            ")")

    # merge transcripts to get consensus gene models
    tmpf = "output_ucsc_track_hub/lncRNA_bigbed2/lncRNA_gene.gtf" #P.getTempFilename("/ifs/scratch")
    to_cluster = False
    statement = ("zcat %(gtffile)s |"
                 " python %(scriptsdir)s/gtf2gtf.py"
                 "  --method=merge-exons"
                 "  --log=/dev/null |"
                 " python %(scriptsdir)s/gtf2gtf.py"
                 "  --method=set-transcript-to-gene"
                 "  --log=/dev/null"
                 " > %(tmpf)s")
    P.run()

    # create bed12 
    tmp_gene = "output_ucsc_track_hub/lncRNA_bigbed2/lncRNA_gene.bed" #P.getTempFilename("/ifs/scratch")
    PipelineLncRNA.gtfToBed12(tmpf, tmp_gene, "gene")

    # convert to bigBed
    tmp_gene_1 = "output_ucsc_track_hub/lncRNA_bigbed2/lncRNA_gene_sorted.bed" #P.getTempFilename("/ifs/scratch")
    statement = ("cat %(tmp_gene)s |"
                 " sort -k1,1 -k2,2n"
                 " > %(tmp_gene_1)s;"
                 " /ifs/home/jethro/Downloads/bedToBigBed -type=bed12 -as=%(out_es)s -extraIndex=name %(tmp_gene_1)s %(contigs)s %(out_gene)s")
                 #" /ifs/home/jethro/Downloads/bedToBigBed %(tmp_gene_1)s %(contigs)s %(out_gene)s")
    P.run()

    tmp_transcript = "output_ucsc_track_hub/lncRNA_bigbed2/lncRNA_transcript.bed" #P.getTempFilename("/ifs/scratch")
    PipelineLncRNA.gtfToBed12(gtffile, tmp_transcript, "transcript")
    tmp_transcript_1 = "output_ucsc_track_hub/lncRNA_bigbed2/lncRNA_transcript_sorted.bed" #P.getTempFilename("/ifs/scratch")
    statement = ("cat %(tmp_transcript)s |"
                 " sort -k1,1 -k2,2n"
                 " > %(tmp_transcript_1)s;"
                 " /ifs/home/jethro/Downloads/bedToBigBed -type=bed12 -as=%(out_es)s -extraIndex=name %(tmp_transcript_1)s %(contigs)s %(out_transcript)s")
                 #" bedToBigBed %(tmp_transcript_1)s %(contigs)s %(out_transcript)s")
    P.run()


@follows(convertH3K4me1BamfilesToBigWigs,
         convertInputBamfilesToBigWigs,
         convertH3K4me3BamfilesToBigWigs,
         convertChIPSeqPeaksToBigBed,
         convertRNASeqBamfilesToBigWigs,
         convertlncRNAToBigBed)
def createUCSCTrackFiles():
    pass


@merge([convertH3K4me1BamfilesToBigWigs,
        convertInputBamfilesToBigWigs,
        convertH3K4me3BamfilesToBigWigs,
        convertChIPSeqPeaksToBigBed,
        convertMergedChIPSeqPeaksToBigBed,
        convertRNASeqBamfilesToBigWigs,
        convertlncRNAToBigBed],
       "/ifs/projects/proj010/web/UCSC_Track_Tub/hub.txt")
def createUCSCTrackHub(infiles, outfile):
    """
    Iterate through all infiles and create dictionary that contains 
    {filetype: filename,}. 
    Specify trackhub prefix as 'Bcell_lncRNA'
    Specify project name as 'LncRNAs in B cell biology'
    Leave project_id as blank, as it's not used
    Specify hub_id as 'ucsc_track_hub'
    """

    # create track hub dict
    export_files = collections.OrderedDict()

    # set up different composite track lists.
    for i in ["lncrna_loci", "rnaseq_pro", "rnaseq_pre", "rnaseq_immature",
              "rnaseq_mature", "rnaseq_follicular", "rnaseq_marginal", "rnaseq_b1a",
              "rnaseq_germinal", "chip_pro", "chip_pre", "chip_immature", "chip_mature",
              "chip_follicular", "chip_marginal", "chip_b1a", "chip_germinal",
              "H3K4me1_loci", "H3K4me3_loci"]:
        export_files[i] = []

    # iterate through the output directories, put the respective files into the 
    # respective composite track list. 
    for inf in infiles:
        if inf.endswith(".log"):
            continue
        dirname = os.path.split(os.path.dirname(inf))[1]
        # lncRNA gene/transcript models
        if dirname == "lncRNA_bigbed":
            export_files["lncrna_loci"].append(inf)
        # chipseq peak intervals
        elif dirname == "chipseq_peak_intervals":
            file_name = os.path.basename(inf)
            if re.search("H3K4me1", inf):
                export_files["H3K4me1_loci"].append(inf)
            elif re.search("H3K4me3", inf):
                export_files["H3K4me3_loci"].append(inf)
            else:
                raise ValueError("Unrecognised file %s" % inf)
        # rnaseq alignment bigwigs
        elif re.search("rnaseq", dirname):
            cell_type = os.path.basename(inf).split("-")[0]
            comp_track_list = "rnaseq_" + cell_type
            export_files[comp_track_list].append(inf)
        # chipseq alignment bigwigs
        elif dirname == "chipseq_bigwigs":
            # ignore inputs
            if re.search("input", os.path.basename(inf)):
                continue
            cell_type = os.path.basename(inf).split("-")[0]
            comp_track_list = "chip_" + cell_type
            export_files[comp_track_list].append(inf)
        else:
            raise ValueError("Unrecognised output directory %s" % dirname)

    # sort the data in each composite track list
    for track in export_files.keys():
        export_files[track] = sorted(export_files[track])

    P10.publish_tracks(export_files,
                       PARAMS,
                       prefix="Bcell-",
                       project_name="LncRNAs in B cells",
                       hub_id = "UCSC_Track_Hub")


    # export_files = collections.defaultdict(list)
    # for infile in infiles:
    #     if infile.endswith("bw"):
    #         export_files["bigWig"].append(infile)
    #     elif infile.endswith("bb"):
    #         export_files["bigBed"].append(infile)
    #     elif infile.endswith("bb12"):
    #         export_files["bigBed12"].append(infile)
    #     else:
    #         raise IOError("Unrecognized file type %s" % os.path.basename(infile))

    # PipelinePublishing.publish_tracks(export_files,
    #                                   PARAMS,
    #                                   prefix="Bcell-",
    #                                   project_name="LncRNAs in B cell biology",
    #                                   hub_id = "UCSC_track_hub")

# #################################################################################
# #################################################################################
# #################################################################################
# ## primary targets
# #################################################################################
# @follows( )
# def full(): pass

# @follows( loadFPKMLncRNAProteinCodingPearsonCorrelations,
#           loadFPKMMeanLncRNAProteinCodingPearsonCorrelations )
# #          loadCountLncRNAProteinCodingPearsonCorrelations,
# #          loadResidualLncRNAProteinCodingPearsonCorrelations,
# #          loadFittedLncRNAProteinCodingPearsonCorrelations )
# def calculateCorrelations(): pass

# #################################################################################
# #################################################################################
# #################################################################################
# ## primary targets
# #################################################################################
@follows( filterMELncRNA,
          filterSELncRNA,
          filterLncRNA,)
def runLncRNAFiltering():
    pass

#@follows( runLncRNALocationStats, 
          



@follows( mkdir( "report" ) )
def build_report():
    '''build report from scratch.'''

    E.info( "starting report build process from scratch" )
    P.run_report( clean = True )

@follows( mkdir( "report" ) )
def update_report():
    '''update report.'''

    E.info( "updating report" )
    P.run_report( clean = False )

@follows( update_report )
def publish_report():
    '''publish report.'''

    # files to be made available on publishing. 
    export_files = {
        # wgcna module heatmaps
        os.path.join(PARAMS["report_prefix_dir"], "wgcna", "high_fpkm", "module_heatmaps"):
        glob.glob("./expression_wgcna_filtered_high_fpkm/module_heatmaps_zscore/*"),
        os.path.join(PARAMS["report_prefix_dir"], "wgcna", "high_cv", "module_heatmaps"):
        glob.glob("./expression_wgcna_filtered_high_cv/module_heatmaps_zscore/*"),
        # wgcna dendrograms of lncRNA position within modules
        os.path.join(PARAMS["report_prefix_dir"], "wgcna", "high_fpkm", "module_dendrograms"):
        glob.glob("./expression_wgcna_filtered_high_fpkm/module_dendrograms/*"),
        os.path.join(PARAMS["report_prefix_dir"], "wgcna", "high_cv", "module_heatmaps"):
        glob.glob("./expression_wgcna_filtered_high_cv/module_heatmaps_zscore/*"),      
    }

    E.info( "publishing report" )
    P.publish_report()

if __name__== "__main__":
    # P.checkFiles( ("genome.fasta", "genome.idx" ) )
    sys.exit( P.main(sys.argv) )
