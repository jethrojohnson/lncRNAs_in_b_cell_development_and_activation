#################################################################################
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
==================================================
PipelineProj010QC.py - custom QC tasks for proj010
==================================================

:Author: Jethro Johnson
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This module file contains generic QC functions for use in the script 
pipeline_proj010_lncrna.py

"""

import gzip, os, re
import pysam
import random

import CGATPipelines.Pipeline as P
import CGAT.Experiment as E
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools

#################################################################################
# section: general functions
#################################################################################

def filterBam( inf, outf, params ):
    """
    This is set up to be run with P.submit; it takes the infile and outfile as 
    elements of a list passed to params.
    The infile is filtered using pysam. Any read that are flagged as unmapped 
    (0x4), or as secondary (0x100) are removed. The latter does not apply to BWA
    mapped reads when the remove_non_unique option is set when running the 
    mapping pipeline.
    """
    infile, outfile = params
    samfile = pysam.Samfile( infile, "rb" )
    mappedreads = pysam.Samfile( outfile, "wb", template=samfile )
    for read in samfile.fetch():
        if read.is_unmapped:
            continue
        elif read.is_secondary:
            continue
        else:
            mappedreads.write( read )
    mappedreads.close()
    samfile.close()

    
def combineTables( file_names, table_name, headers = False, transpose = False ):

    if headers: 
        if transpose: 
            statement = '''
                        python %(scriptsdir)s/combine_tables.py
                            --columns=1
                            --skip-titles
                            --header-names=%(headers)s
                            --missing-value=0
                            --sort-keys=numeric
                            --log=%(table_name)s.log
                            %(file_names)s |
                        python %(scriptsdir)s/table2table.py
                            --transpose
                        > %(table_name)s
                        '''
        else: 
            statement = '''
                        python %(scriptsdir)s/combine_tables.py
                            --columns=1
                            --skip-titles
                            --header-names=%(headers)s
                            --missing-value=0
                            --sort-keys=numeric
                            --log=%(table_name)s.log
                            --stdout=%(table_name)s
                            %(file_names)s
                        '''
    else:
        if transpose: 
            statement = '''
                        python %(scriptsdir)s/combine_tables.py
                            --cat=CAT
                            --missing-value=0
                            --sort-keys=numeric
                            --log=%(table_name)s.log
                            %(file_names)s |
                        python %(scriptsdir)s/table2table.py
                            --transpose
                        > %(table_name)s
                        '''
        else: 
            statement = '''
                        python %(scriptsdir)s/combine_tables.py
                            --cat=CAT
                            --missing-value=0
                            --sort-keys=numeric
                            --log=%(table_name)s.log
                            --stdout=%(table_name)s
                            %(file_names)s
                        '''
    P.run()
    return table_name

## The following doesn't work becayse the headers mapped to the first statement
## are not passed to the second statement.
"""
    if headers:
        combine_table_parameters = ''' 
                                   --columns=1
                                   --skip-titles
                                   --header-names= %s
                                   ''' % headers

    else: 
        combine_table_parameters = '''
                                   --cat=CAT
  
                                 '''

    if transpose:
        statement = '''
             python %(scriptsdir)s/combine_tables.py
                 %(combine_table_parameters)s
                 --missing-value=0
                 --log=%(table_name)s.log
                 %(file_names)s |
             python %(scriptsdir)s/table2table.py
                 --transpose
             > %(table_name)s
            '''
        P.run()

    else:
        statement = '''
             python /ifs/devel/jethro/cgat/scripts/combine_tables.py
             %(combine_table_parameters)s
             --missing-value=0
             --log=%(table_name)s.log
             --stdout=%(table_name)s
             %(file_names)s
            '''
        P.run()

    return table_name
"""


## It would be nice to turn extractTables into a class, but I don't know how to
## pass arguments to the __call__ method. 
#class extractTables( object ): 
#"""
#Base class for extracting tables from text files, will take a regular expression
#to identify the start of the table, then subsequently returns all lines up to 
#the first empty line. 
#"""
#statement = None
#
#def __init__( self, *args, **kwargs ):
#    pass
#
#def __call__( ### how to pas arguments to a call function?? )


# the following is based on function loadPicardMetrics in PipelineMappingQC.py
def extractTables( infile, outfile, table_id ):
    infile_lines = IOTools.openFile( infile, "r" ).readlines()
    outf = IOTools.openFile( outfile, "w" )
    # create regular expression object that denotes start of histogram
    regex = re.compile( table_id ) 

    # then search through file and find start of histogram
    # (enumerate() returns a dict where key is an integer, and value is a line)
    for line_number, line in enumerate( infile_lines ): 
        # (.search() searches for string and returns MatchObject instance)
        if regex.search( line ):
            histogram_lines = infile_lines[line_number+1:]
            break

    # then search through file and find end of histogram
    for line_number, line in enumerate(histogram_lines):
        if not line.strip():
            histogram_lines = histogram_lines[:line_number]
            break

    # put in a check to make sure there are right number of fields
#    field_number = len( histogram_lines[0][:-1].split("\t") )

    # then write histogram to outfile
#    for line in histogram_lines:
#        fields = len( line[:-1].split("\t") )
#        if fields == field_number:
#            outf.write( line )
#        else:
#            raise ValueError( "file has %s fields, not %s " % ( fields, field_number ) )

    for line in histogram_lines: 
        outf.write( line )

    outf.close()
    return outfile

def concatTables( infiles, outfile, headers, keep_table = False ): 
    if keep_table: 
        outf = IOTools.openFile( outfile, "w" )
    else:
        outf = P.getTempFile()
    outf.write( headers )
    for inf in infiles:
        lines = IOTools.openFile( inf ).readlines()
        line = lines[0]
        outf.write( line )
    outf.close()
    return outf.name

#################################################################################
# section: calculating coverage over unannotated regions of the genome
#################################################################################

def createUnannotatedRegions( annotation_file, 
                              contig_file, 
                              outfile, 
                              increase, 
                              filter_contigs ): 
    to_cluster = True
    statement = ''' 
        zcat %(annotation_file)s |
        python %(scriptsdir)s/gtf2gtf.py
            --method=merge-transcripts
            --log=%(outfile)s.log |
        python %(scriptsdir)s/gff2bed.py |
        slopBed 
            -b %(increase)s 
            -i stdin 
            -g %(contig_file)s |
        complementBed
            -i stdin
            -g %(contig_file)s |
        awk '$1 !~ /%(filter_contigs)s/' |
        gzip > %(outfile)s
        '''
    P.run()


def createGenomeWindows( genome, outfile, windows, filter_contigs ):
    to_cluster = True
    filter_contigs = "'" + filter_contigs + "'"
    statement = ( "python %(scriptsdir)s/index2bed.py"
                  "  --genome=%(genome)s"
                  "  --fixed-width-windows=%(windows)s"
                  "  --remove-regex=%(filter_contigs)s"
                  "  --log=%(outfile)s.log |"
                  " gzip > %(outfile)s" )
    P.run()


def createFeatureWindows( feature, genome_windows, outfile, repeat_regions ):
    to_cluster = True
    statement = '''
        bedtools intersect
            -wa
            -a %(genome_windows)s
            -b %(feature)s
            2> %(outfile)s.log |
        bedtools intersect
            -v
            -a stdin
            -b %(repeat_regions)s
            2>> %(outfile)s.log |
        gzip > %(outfile)s
        '''
    P.run()


def calculateCoverage( bamfile, bedfile, outfile ):
    to_cluster = True
    statement = '''
        coverageBed
            -abam %(bamfile)s
            -b %(bedfile)s |
        cut -f4 |
        python /ifs/devel/jethro/cgat/scripts/data2histogram.py
            --bin-size=1
            --no-titles
            --log=%(outfile)s.log |
        gzip > %(outfile)s
        '''
    P.run()


#################################################################################
# section: running Picard tools
#################################################################################

def getPicardOptions():
    return "-pe dedicated 3 -R y -l mem_free=1.4G -l picard=1"

def runPicardMarkDuplicates( bamfile, 
                             outfile, 
                             retain_bam = False, 
                             remove_duplicates=True ):
    to_cluster = True
    job_options = getPicardOptions()
    if retain_bam:
        if remove_duplicates:
            outfile = ( P.snip( outfile, "_duplicate.metrics" ) )
            statement = '''
                MarkDuplicates
                    INPUT=%(bamfile)s
                    OUTPUT=%(outfile)s_deduped.bam
                    METRICS_FILE=%(outfile)s_duplicate.metrics
                    REMOVE_DUPLICATES=true
                    ASSUME_SORTED=true
                    VALIDATION_STRINGENCY=SILENT; 
                checkpoint;
                samtools index %(outfile)s_deduped.bam
                '''
        else:
            outfile = ( P.snip( outfile, "_duplicate.metrics" ) )
            statement = '''
                MarkDuplicates
                    INPUT=%(bamfile)s
                    OUTPUT=%(outfile)s_marked.bam
                    METRICS_FILE=%(outfile)s_duplicate.metrics
                    ASSUME_SORTED=true
                    VALIDATION_STRINGENCY=SILENT; 
                checkpoint;
                samtools index %(outfile)s_deduped.bam
                '''
    else: 
        statement = '''
            MarkDuplicates
                INPUT=%(bamfile)s
                OUTPUT=/dev/null
                METRICS_FILE=%(outfile)s
                ASSUME_SORTED=true
                VALIDATION_STRINGENCY=SILENT
            '''
    P.run()

def runFastqToSam( infile, bamfile, filename, quality_format="Standard" ):
    to_cluster = True
    sample = P.snip( infile, ".fastq.1.gz" )
    job_options = getPicardOptions()
    statement = ( "FastqToSam"
                  " FASTQ=%(infile)s"
                  " FASTQ2=%(filename)s.fastq.2.gz"
                  " OUTPUT=/dev/stdout"
                  " QUALITY_FORMAT=%(quality_format)s"
                  " SAMPLE_NAME=%(filename)s | "
                  "SortSam"
                  " INPUT=/dev/stdin" 
                  " SORT_ORDER=coordinate"
                  " OUTPUT=%(bamfile)s" )
    P.run()
    return bamfile

def runEstimateLibraryComplexity( bamfile, outfile ):
    to_cluster = True
    job_options = getPicardOptions()
    statement = ("EstimateLibraryComplexity" 
                     " I=%(bamfile)s"
                     " O=%(outfile)s" )
    P.run()

def extractPicardStats( infile, 
                        outfile, 
                        table_id, 
                        track, 
                        track_to_replace ):
    output = extractTables( infile, outfile, table_id )
    statement = ( "sed -i 's/%(track_to_replace)s/%(track)s/' %(output)s " )
    P.run()

def combinePicardStats( file_names, outfile, table_name ):
    combined_tables = combineTables( file_names, table_name )
    statement = ( "sed -i 's|./qc_unmapped_stats/||' %(table_name)s;"
                  " sed -i 's|.complexity_stats||' %(table_name)s" )
    P.run()
    return combined_tables

def combinePicardHistograms( file_names, table_name ): 
    statement = ( "python %(scriptsdir)s/combine_tables.py"
                  " --columns=1"
                  " --missing-value=0"
                  " --sort-keys=numeric"
                  " --log=%(table_name)s.log"
                  " --stdout=%(table_name)s"
                  " %(file_names)s" )
    P.run()
    return table_name

#################################################################################
# section: running CHANCE
#################################################################################

def combineCHANCEStats( infiles, outfile  ):
    outfile = open( outfile, "w" )

    headers = [ "sample_id", ]
    for line in IOTools.openFile( infiles[0], "r" ).readlines():
        if line.strip():
            line = line.strip()
            headers.append( line.split(",")[0] )
        else: break
    outfile.write( "\t".join( headers ) + "\n" )

    for infile in infiles:
        inf = P.snip( os.path.basename( infile ), ".chance" )
        infile_lines = IOTools.openFile( infile, "r" ).readlines()

        results = [ inf, ]
        for line in infile_lines:    
            if line.strip():
                line = line.strip()
                results.append( line.split(",")[1] )
            else: break
        outfile.write( "\t".join( results ) + "\n" )
        
#################################################################################
# section: running counting reads across gene models
#################################################################################

def calculateIntronExonReadRatio( bamfile, gene_models, exons, outfile, feature="genes" ):

    tmpf_1 = P.getTempFilename("/scratch") # will contain exon read counts. 
#    P.info( "tmpf_1 is %s" % tmpf_1 )
    tmpf_2 = P.getTempFilename("/scratch") # will contain intron read counts 
#    P.info( "tmpf_2 is %s" % tmpf_2 )
    tmpf_3 = P.getTempFilename(".") # will contain both
#    P.info( "tmpf_3 is %s" % tmpf_3 )

    to_cluster = True
    statement = (  # write number of reads across exons and exon length to tmpf_1
        "zcat %(gene_models)s |"
        " python %(scriptsdir)s/gtf2table.py"
        " --reporter=%(feature)s"
        " --bam-file=%(bamfile)s"
        " --counter=read-counts"
        " --column-prefix='count_'"
        " --counter=length"
        " --column-prefix='length_'"
        " --log=%(outfile)s.log |"
        " cut -f 1,3,14 |"
        " sed 1d |"
        " sort -k 1,1"
        " > %(tmpf_1)s;" 

        # convert introns to exons (ripped from pipeline_mapping)
        " zcat %(gene_models)s |" 
        " python %(scriptsdir)s/gtf2gtf.py"
        " --method=exons2introns"
        " --intron-min-length=100"
        " --intron-border=10"
        " --log=%(outfile)s.log |"
        " python %(scriptsdir)s/gff2gff.py"
        " --crop=%(exons)s"
        " --log=%(outfile)s.log |"
        " python %(scriptsdir)s/gtf2gtf.py"
        " --method=set-transcript-to-gene"
        " --log=%(outfile)s.log |"
        " sed 's/\\tintron\\t/\\texon\\t/g' |"
        
        # write number of reads across introns and intron length to tmpf_2
        " python %(scriptsdir)s/gtf2table.py" 
        " --reporter=%(feature)s"
        " --bam-file=%(bamfile)s"
        " --counter=read-counts"
        " --column-prefix='count_'"
        " --counter=length"
        " --column-prefix='length_'"
        " --log=%(outfile)s.log |"
        " cut -f 1,3,14 |"
        " sed 1d |"
        " sort -k 1,1"
        " > %(tmpf_2)s;" 

        # join tmpf_1 and tmpf_2 together on the basis of gene_id
        " join -j 1 %(tmpf_1)s %(tmpf_2)s"
        " > %(tmpf_3)s"
        " 2>> %(outfile)s.log" )
    P.run()

    outf = gzip.open( outfile, "wb" )
    outf.write("gene_id\tie_ratio\n" )

    inf = open( tmpf_3, "r" ).readlines()
    for line in inf:
        line = line.split()
        gene_id = line[0]
        exon_coverage = (float(line[1])/float(line[2]))
        intron_coverage = (float(line[3])/float(line[4]))
        if exon_coverage == 0 and intron_coverage == 0:
            outf.write( gene_id + "\tNA\n" )
        elif exon_coverage == 0 and intron_coverage != 0:
            outf.write( gene_id + "\tInf\n" )
        else:
            intron_exon_ratio = (intron_coverage/exon_coverage)
            outf.write( gene_id + "\t" + str(intron_exon_ratio) + "\n" )
    outf.close()

    os.unlink(tmpf_1)
    os.unlink(tmpf_2)
    os.unlink(tmpf_3)

