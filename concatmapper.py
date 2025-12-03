#!/usr/bin/env python

"""
Copyright 2025  Daryl Gohl

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from scipy.interpolate import interp1d
import numpy as np
import pysam

__version__ = "1.0"

def get_args(x):
    x.add_argument("-q", "--fastq_file",
                   type = str,
                   default = None,
                   metavar = '',
                   help="Input path for fastq sequencing reads file [required].")
    x.add_argument("-r", "--fasta_file",
                   type = str,
                   default = None,
                   metavar = '',
                   help="Input path for fasta reference file [required].")
    x.add_argument("-o", "--output_dir",
                   type = str,
                   default = '',
                   metavar = '',
                   help = "Output directory for sam file and plot (default: same folder as input fastq file)")
    x.add_argument("-n", "--output_file",
                   type = str,
                   default = '',
                   metavar = '',
                   help = "Output name for sam file and plot (default: same name as input fastq file)")
    x.add_argument("-m", "--min_length",
                   type = int,
                   default = 100,
                   metavar = '',
                   help="Minimum mapped read length to plot (default: 100)")
    x.add_argument("-u", "--unsorted",
                   type = str,
                   default = False,
                   metavar = '',
                   help="Plot from unsorted sam file (default: False - plot sorted sam file)")                   
    x.add_argument("-l", "--line_spacing",
                   type = float,
                   default = 0.02,
                   metavar = '',
                   help = "Radial spacing of each read on plot (default: 0.2)")
    x.add_argument("-w", "--line_width",
                   type = float,
                   default = 0.75,
                   metavar = '',
                   help = "Line width of each read on plot (default: 0.75)")
    x.add_argument("-c", "--circle_size",
                   type = float,
                   default = 0.45,
                   metavar = '',
                   help = "Size of central circle (default: 0.45)")
    x.add_argument("-s", "--fig_size",
                   type = float,
                   default = 10,
                   metavar = '',
                   help = "Size of figure (default: 10)")
    x.add_argument("-x", "--clip",
                  type = str,
                  default = False,
                  metavar = '',
                  help = "Plot clipped portion of reads (default: False)")
    x.add_argument("-f", "--figure_format",
                   type = str,
                   default = "pdf",
                   metavar = '',
                   help = "Format of saved figure, supported formats: eps, jpeg, jpg, pdf, pgf, png, ps, raw, rgba, svg, svgz, tif, tiff.")
#####################
    args = x.parse_args()
    return args

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('True', 'yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('False', 'no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def main():
    # Parses command line arguments
    argparser = argparse.ArgumentParser(description = "Map and plot reads against a circular reference (v" + __version__ +")\n" + \
                                                  "by Daryl Gohl\n" + \
                                                  "This program takes in a fasta reference file and a fastq sequencing file and concatenates the reference sequence (two copies of the reference sequence repeated in tandem) and then uses minimap2 to map reads to the concatenated reference, generating a .sam file that can be used for subsequent plotting.",
                                    add_help = True,
                                    epilog ='')
    args = get_args(argparser)
    args = vars(args)

    #Possible input args
    fastq_filename = args['fastq_file']
    fasta_reference = args['fasta_file']
    out_folder = args['output_dir']
    out_file = args['output_file']
    incrementor = args['line_spacing']
    width = args['line_width']
    Min_length = int(args['min_length'])
    Fig_size = args['fig_size']
    start_coord = args['circle_size']
    Fig_format = args['figure_format']
    Include_clipped = str2bool(args['clip'])
    Sorted_Unsorted = str2bool(args['unsorted'])

    #Set output directory
    if out_folder == '':
      out_dir = os.path.dirname(fastq_filename)
    else:
      out_dir = out_folder

    os.chdir(out_dir)

    #Set output file name
    if out_file == '':
      out_file = os.path.basename(fastq_filename)[:-6] + ".sam"
    else:
      out_file = out_file  + ".sam"

    #Read in reference and generate concatenated reference file
    refid = SeqIO.read(fasta_reference, "fasta").id
    refseq = SeqIO.read(fasta_reference, "fasta").seq
    #print(refid)
    #print(refseq)
    refid_new = str(refid) + "_CONCAT"
    concat_seq = str(refseq)+str(refseq)
    reference_len = len(str(refseq))

    record = SeqRecord(
        Seq(concat_seq),
        id=refid_new,
        name=refid_new,
        description="concatenated reference file for ConcatMap",
    )

    with open(refid_new, "w") as output_handle:
        SeqIO.write(record, output_handle, "fasta")

    #Map with minimap2 and generate sam file
    execute = "minimap2 --sam-hit-only -a " + refid_new + " " + fastq_filename + " > " + out_file
    os.system(execute)

    filename = out_file
    ref_length = reference_len

    Plot_clipped = str2bool(Include_clipped)

    start_coord_orig = start_coord

    #SAM File
    #unsorted_samfile = pysam.AlignmentFile(filename, "r")
    out_name = out_file[:-4] + "_sorted.sam"
    pysam.samtools.sort("-o", out_name, filename, catch_stdout=False)
    if Sorted_Unsorted != False:
    	samfile = pysam.AlignmentFile(filename, "r")
    else:
    	samfile = pysam.AlignmentFile(out_name, "r")

    #reads
    reference_start = []
    reference_end = []
    clipped_start = []
    clipped_end = []
    clipped_start_len = []
    clipped_end_len = []
    for r in samfile.fetch(until_eof=True):
        x = r.query_alignment_start #position of start of alignment in read
        y = r.query_alignment_end #position of end of alignment in read
        z = r.query_alignment_length #length of alignment in read
        q = r.infer_read_length() #length of read including soft-clipped bases
        x1 = r.reference_start #position of start of alignment in reference
        y1 = r.reference_end #position of end of alignment in reference
        z1 = r.reference_length #length alignment to reference
        qs0 = x1-x #start position of clipped bases relative to the reference upstream
        qs1 = x #end position of clipped bases relative to the reference upstream
        cs_len = x #length of sequence clipped off the upstream end
        ce_len = q-y #length of sequence clipped off the downstream end
        qe0 = y #start position of clipped bases relative to the reference downstream
        qe1 = y1+(q-y) #end position of clipped bases relative to the reference downstream
        if z1 > Min_length:
            reference_start.append(x1)
            reference_end.append(y1)
            clipped_start.append(qs0) #clipped sequence start relative to reference
            clipped_end.append(qe1) #clipped sequence end relative to reference
            clipped_start_len.append(cs_len) #length of upstream clipped sequence
            clipped_end_len.append(ce_len) #length of downstream clipped sequence
            #print str(r.infer_read_length())
            #pr = "clip start,end = " + str(qs0) + "," + str(qe1) + " mapped start, end = " + str(x1) + "," + str(y1)
            #print pr


    ###MAPPED READS
    #Dedup reference
    break_span = []
    ref_start_collapsed = []
    ref_end_collapsed = []
    for i, item in enumerate(reference_start):
        if item < ref_length and reference_end[i] > ref_length:
            break_span.append(True)
        else:
            break_span.append(False)
        if item > ref_length:
            ref_start_collapsed.append(item-ref_length)
        else:
            ref_start_collapsed.append(item)
        if reference_end[i] > ref_length:
            ref_end_collapsed.append(reference_end[i]-ref_length)
        else:
            ref_end_collapsed.append(reference_end[i])

    #Convert to polar coordinated (0-360)
    deg_min = 0
    deg_max = 360
    deg_per_base = 360.0/ref_length
    deg_start = []
    deg_end = []
    for i, item in enumerate(ref_start_collapsed):
        deg_end.append(360-item*deg_per_base) #Note: Subtracting from 360 makes degrees 0 to 360 go in more intuitive clockwise orientation rather than default counterclockwise.
        deg_start.append(360-ref_end_collapsed[i]*deg_per_base) #Note: Subtracting from 360 makes degrees 0 to 360 go in more intuitive clockwise orientation rather than default counterclockwise.


    ###FULL READS (with clipped relative to reference)
    #First, clean up <0 or >ref_length cases
    #If <0, add ref_length
    #If >ref_length, subtract ref_length
    #This will turn these cases into break spanning reads which will be handled normally below
    clipped_start_clean = []
    clipped_end_clean = []
    mod = []
    for i, item in enumerate(clipped_start):
      if item < 0:
        clipped_start_clean.append(item + ref_length)
        clipped_end_clean.append(clipped_end[i] + ref_length)
        mod.append(True)
      elif clipped_end[i] > 2*ref_length:
        clipped_start_clean.append(item - ref_length)
        clipped_end_clean.append(clipped_end[i] - ref_length)
        mod.append(False)
      else:
        clipped_start_clean.append(item)
        clipped_end_clean.append(clipped_end[i])
        mod.append(False)

    clipped_break_span = []
    clipped_start_collapsed = []
    clipped_end_collapsed = []
    for i, item in enumerate(clipped_start_clean):
        if item < ref_length and clipped_end_clean[i] > ref_length:
            clipped_break_span.append(True)
        else:
            clipped_break_span.append(False)
        if item > ref_length:
            clipped_start_collapsed.append(item-ref_length)
        else:
            clipped_start_collapsed.append(item)
        if clipped_end[i] > ref_length:
            clipped_end_collapsed.append(clipped_end_clean[i]-ref_length)
        else:
            clipped_end_collapsed.append(clipped_end_clean[i])

    #Convert to polar coordinated (0-360)
    deg_min = 0
    deg_max = 360
    deg_per_base = 360.0/ref_length
    clipped_deg_start = []
    clipped_deg_end = []
    for i, item in enumerate(clipped_start_collapsed):
        clipped_deg_end.append(360-item*deg_per_base)
        clipped_deg_start.append(360-clipped_end_collapsed[i]*deg_per_base)

    #Plotting data
    with plt.style.context('ggplot'):
        fig = plt.figure(figsize=(Fig_size,Fig_size))
        #my_dpi = 300
        #fig = plt.figure(figsize=(2400/my_dpi, 2400/my_dpi), dpi=my_dpi)
        ax = fig.add_subplot(111, projection="polar")
        ax.grid(False)
        ax.set_rticks([])
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.set_theta_zero_location('N')
        ax.set_facecolor('white')
        ax.axis('off')

        # Connect two points with a curve
        for curve in [[[0, 360], [(start_coord-2*incrementor), (start_coord-2*incrementor)]]]:
            curve[0] = np.deg2rad(curve[0])
            x = np.linspace( curve[0][0], curve[0][1], 500)
            y = interp1d( curve[0], curve[1])( x)
            ax.plot(x, y, linewidth = 5)
        for curve in [[[0, 0], [0.0, 0.0]]]:
            curve[0] = np.deg2rad(curve[0])
            x = np.linspace( curve[0][0], curve[0][1], 500)
            y = interp1d( curve[0], curve[1])( x)
            ax.plot(x, y)

        #Draw clipped lines
        if Plot_clipped != False:
          for i, item in enumerate(clipped_break_span):
              if mod[i] == True:
                item = break_span[i]
              if item == False:
                  for curve in [[[clipped_deg_start[i], clipped_deg_end[i]], [start_coord, start_coord]]]:
                      curve[0] = np.deg2rad(curve[0])
                      x = np.linspace( curve[0][0], curve[0][1], 500)
                      y = interp1d( curve[0], curve[1])( x)
                      ax.plot(x, y, color='red', linewidth = width)
              else:
                  for curve in [[[360, clipped_deg_start[i]], [start_coord, start_coord]]]:
                      curve[0] = np.deg2rad(curve[0])
                      x = np.linspace( curve[0][0], curve[0][1], 500)
                      y = interp1d( curve[0], curve[1])( x)
                      ax.plot(x, y, color='red', linewidth = width)
                  for curve in [[[clipped_deg_end[i], 0], [start_coord, start_coord]]]:
                      curve[0] = np.deg2rad(curve[0])
                      x = np.linspace( curve[0][0], curve[0][1], 500)
                      y = interp1d( curve[0], curve[1])( x)
                      ax.plot(x, y, color='red', linewidth = width)
              start_coord = start_coord + incrementor

        #Mapped reads
        start_coord = start_coord_orig
        for i, item in enumerate(break_span):
            if item == False:
                for curve in [[[deg_start[i], deg_end[i]], [start_coord, start_coord]]]:
                    curve[0] = np.deg2rad(curve[0])
                    x = np.linspace( curve[0][0], curve[0][1], 500)
                    y = interp1d( curve[0], curve[1])( x)
                    ax.plot(x, y, color='grey', linewidth = width)
            else:
                for curve in [[[360, deg_start[i]], [start_coord, start_coord]]]:
                    curve[0] = np.deg2rad(curve[0])
                    x = np.linspace( curve[0][0], curve[0][1], 500)
                    y = interp1d( curve[0], curve[1])( x)
                    ax.plot(x, y, color='grey', linewidth = width)
                for curve in [[[deg_end[i], 0], [start_coord, start_coord]]]:
                    curve[0] = np.deg2rad(curve[0])
                    x = np.linspace( curve[0][0], curve[0][1], 500)
                    y = interp1d( curve[0], curve[1])( x)
                    ax.plot(x, y, color='grey', linewidth = width)
            start_coord = start_coord + incrementor

    #Set figure output file name
    file = out_file[:-4] + "." + Fig_format
    fig_name = file

    plt.savefig(fig_name,bbox_inches='tight')

if __name__ == "__main__":
    main()
