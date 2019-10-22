#!/usr/bin/env python

import argparse
## Parse command line arguments
parser = argparse.ArgumentParser(description='Count allele frequencies for single cells.')
parser.add_argument("-v", "--verbose", help="set verbose mode", action="store_true")
parser.add_argument("--snps", help="SNPs file, tab separated file with no header and 4 columns: chr, pos (0-indexed), reference base, alternate base",required=True);
parser.add_argument("--barcodes", help="barcodes file containing only barcodes of interest from the CB tag in the bam file; one barcode per line",required=True);
parser.add_argument("--output-format", help="output format for resulting tables; mm for Matrix Market, hdf for HDF", default="mm", choices=['mm'], dest="output_format"); # TODO: add hdf
parser.add_argument("--bamfile", help="sorted and indexed bamfile from single cell experiment")
parser.add_argument("--max-depth", help="max_depth argument for fileup", type=int, default=99999999, dest="maxdepth")
parser.add_argument("--output-prefix", help="prefix of output files", default="", dest="prefix")
parser.add_argument("--min-base-quality", help="minimum base quality", default=0, dest="minbbasequal")
parser.add_argument("--min-mapping-quality", help="minimum mapping quality", default=0, dest="minmappingquality")
parser.add_argument("--cell-barcode-tag", help="cellular barcode tag", dest="cellbarcodetag", default="CB")
parser.add_argument("--adjust-capq-threshold", help="adjust capq threshold", dest="adjustcapqthresbold", default=0)
args = parser.parse_args()

## This takes a while, delays help if at the top
import pandas as pd
import pysam
import sys
import scipy
import scipy.io
import scipy.sparse
import numpy
import os.path

## Read positions
if args.verbose:
    sys.stdout.write("reading snps file... ")
    sys.stdout.flush();
if not os.path.exists(args.snps):
    print("Error: snps file does not exist");
    sys.exit(1);
with open(args.snps, 'r') as f:
    postable = pd.read_table(f, header=None, sep=' ', names=['chr','pos','ref','alt'])
if args.verbose:
    sys.stdout.write("done\n")
    sys.stdout.flush();
postable["posindex"] = postable["chr"].map(str) + ':' +  postable["pos"].map(str)
    
if args.verbose:
    sys.stdout.write("building positions index... ")
    sys.stdout.flush();
pos_index = {}
for index, row in postable.iterrows():
    pos_index[row['posindex']] = index;
if args.verbose:
    sys.stdout.write("done\n")
    sys.stdout.flush();

## Read barcodes
if args.verbose:
    sys.stdout.write("reading barcodes file... ");
    sys.stdout.flush();
if not os.path.exists(args.barcodes):
    print("Error: barcodes file does not exist");
    sys.exit(1);
with open(args.barcodes,'r') as f:
    barcodetable = pd.read_table(f, header=None, sep=' ', names=['barcode'])
if args.verbose:
    sys.stdout.write("reading barcodes file... ");
    sys.stdout.write("done\n")

if args.verbose:
    sys.stdout.write("building barcodes index... ");
    sys.stdout.flush()
barcode_index = {}
for index, row in barcodetable.iterrows():
    barcode_index[row['barcode']] = index;
if args.verbose:
    sys.stdout.write("done\n");


## Generate dok arrays to store the data
# cells x positions
ncells = barcodetable.shape[0];
npositions = postable.shape[0];
# dok arrays allow fast access
refarray = scipy.sparse.dok_matrix((ncells,npositions),'uint32')
altarray = scipy.sparse.dok_matrix((ncells,npositions),'uint32')
covarray = scipy.sparse.dok_matrix((ncells,npositions),'uint32')


# start pileup
if args.verbose:
    sys.stdout.write("perfoming pileup...");
    sys.stdout.flush();

# open the bam file (must be sorted and indexed)
samfile = pysam.AlignmentFile(args.bamfile,"rb")

# do pileups
i=0;
for pileupcolumn in samfile.pileup(max_depth=args.maxdepth, min_base_quality=args.minbasequal, min_mapping_quality=args.minmappingquality, dest_capq_threshold=args.adjustcapqthreshold):
    i += 1;
    if ( args.verbose and i % 100000 == 0):
        sys.stdout.write('.');
        sys.stdout.flush();
    posindex = str(pileupcolumn.reference_name) +':' + str(pileupcolumn.reference_pos + 1 ); # +1 for correct coordinates
    if posindex in pos_index:
        posi = pos_index[posindex];
        refallele = postable.iloc[posi,]['ref']
        altallele = postable.iloc[posi,]['alt']
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                try:
                    cbtagvalue  = pileupread.alignment.get_tag(args.cellbarcodetag)
                    if cbtagvalue in barcode_index:
                        barcodei = barcode_index[cbtagvalue];
                        readbase = pileupread.alignment.query_sequence[pileupread.query_position];
                        covarray[barcodei,posi] += 1;
                        if (refallele == readbase):
                            refarray[barcodei,posi] += 1;
                        elif (altallele == readbase):
                            altarray[barcodei,posi] += 1;
                except:
                    pass

if args.verbose:
    sys.stdout.write("done\n");

## close the sam file
samfile.close()

# Save coverage as triplets
if args.verbose:
    sys.stdout.write('saving matrices...');
    sys.stdout.flush()

if args.output_format == 'mm':
    # cells x positions
    # ncells = barcodetable.shape[0];
    # npositions = postable.shape[0];
    #postable["posindex"]
    #barcodetable["barcode"]
    
    # todo save the names
    
    scipy.io.mmwrite(target=args.prefix + "covmat.mtx", a=covarray, field='integer');
    scipy.io.mmwrite(target=args.prefix + "refmat.mtx", a=refarray, field='integer');
    scipy.io.mmwrite(target=args.prefix + "altmat.mtx", a=altarray, field='integer');
else:
    sys.stdout.write('unsupported output format');
    sys.stdout.flush()
    
if args.verbose:
    sys.stdout.write('done\n');
    sys.stdout.flush()

# terminate
sys.exit(0)
