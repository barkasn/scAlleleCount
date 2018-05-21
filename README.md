# scAlleleCount

## Installation
scAlleleCount requires python2 with the pandas, pysam, scipy and numpy modules. These are common modules that your administrator has probably already installed. If you don't have them installed the fastest way to get a suitable environment is to install anaconda.

On ubuntu you can install anaconda using the following commands. Please note that the link is up to date as of May 2018, go to https://anaconda.org/ for an updated link.
```{sh}
cd 
wget https://repo.anaconda.com/archive/Anaconda2-5.1.0-Linux-x86_64.sh
chmod u+x Anaconda2-5.1.0-Linux-x86_64.sh
./Anaconda2-5.1.0-Linux-x86_64.sh
```

Install pysam (the only requirement that doesn't come with the default anaconda installation)
```{sh}
conda install -c bioconda pysam
```

## Running scAlleleCount.py
You can call scAlleleCount.py from the command line or using the getFastCellAlleleCount() R function. In both cases you need an indexed and sorted bamfile (see bamtools for information on preparing this) and two tables, a table of snps to interrogate and a set of barcodes (CB) tags corresponding to cells to count.

### Calling scAlleleCount.py from the command line
To call scAlleleCount.py from the command line you need to prepare a SNPS and a BARCODES file for input. The --help parameter provides more information on input arguments. The output will be save with the provided PREFIX (which can be an absolute or relative path) appended with the suffixes: covmat.mtx, refmat.mtx, altmat.mtx, for the coverage, reference allele and alternate allele matrices respectively. \*.mtx files can be read in R with the Matrix::readMM() function.

```{sh}
user@server:~$ ./scAlleleCount.py --help
usage: scAlleleCount.py [-h] [-v] --snps SNPS --barcodes BARCODES
                        [--output-format {mm}] [--max-depth MAXDEPTH]
                        [--output-prefix PREFIX]
                        bamfile

Count allele frequencies for single cells.

positional arguments:
  bamfile               sorted and indexed bamfile from single cell experiment

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         set verbose mode
  --snps SNPS           SNPs file, tab separated file with no header and 4
                        columns: chr, pos (0-indexed), reference base,
                        alternate base
  --barcodes BARCODES   barcodes file containing only barcodes of interest
                        from the CB tag in the bam file; one barcode per line
  --output-format {mm}  output format for resulting tables; mm for Matrix
                        Market, hdf for HDF
  --max-depth MAXDEPTH  max_depth argument for fileup
  --output-prefix PREFIX
                        prefix of output files
```

Example BARCODES file contents:

```{sh}
AAACATACACCACA-2
AAACATACACCGAT-2
AAACATACACCTAG-2
AAACATACACGGGA-4
AAACATACAGTCTG-2
AAACATACCACACA-4
AAACATACCACTCC-4
AAACATACCCGTAA-4
AAACATACCGACAT-1
AAACATACCTCCAC-1
```

Example SNPS file contents:
```{sh}
10 93603 C T
10 93816 C T
10 93945 G A
10 94026 G A
10 94083 C T
10 94545 C T
10 94870 A G
10 94872 A C
10 95074 G A
10 95096 A G
```

### Using getFastCellAlleleCount() 
getFastCellAlleleCount() is a wrapper around scAlleleCount.py. Provided with a barcodes character vector and a dataframe of SNPs, it will generate the intermediate files in a temporary directory, run scAlelleCount.py, read back the output and return it.

Example Usage
```{r}
bamFile <- "mybam.chr10.bam"
snps <- read.table('chr10snps.txt',stringsAsFactors=F)
cellBarcodes <- read.table('barcodes.txt',stringsAsFactors=F)$V1

> head(snps)
   V1    V2 V3 V4
 1 10 93603  C  T
 2 10 93816  C  T
 3 10 93945  G  A
 4 10 94026  G  A
 5 10 94083  C  T
 6 10 94545  C  T
 
 > head(cellBarcodes)
 [1] "AAACATACACCACA-2" "AAACATACACCGAT-2" "AAACATACACCTAG-2" "AAACATACACGGGA-4"
 [5] "AAACATACAGTCTG-2" "AAACATACCACACA-4"


x <- getFastCellAlleleCount(snps, bamFile, cellBarcodes)

> str(x,1)
 List of 3
  $ refmat:Formal class 'dgTMatrix' [package "Matrix"] with 6 slots
  $ altmat:Formal class 'dgTMatrix' [package "Matrix"] with 6 slots
  $ covmat:Formal class 'dgTMatrix' [package "Matrix"] with 6 slots

```
