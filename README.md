# CopywriteR 

Current methods for detection of copy number aberrations (CNA) from targeted
sequencing data use the depth of coverage of captured exons only. Accurate CNA
determination is complicated by uneven genomic distribution and non-uniform
capture efficiency of targeted exons. Here we present a new tool, CopywriteR, which
eludes these problems by exploiting ‘off-target’ sequence reads. CopywriteR allows
for extracting uniformly distributed copy number information, can be used
without reference and can be applied to sequencing data obtained from various
techniques including chromatin immunoprecipitation and target enrichment on
small-size gene panels. CopywriteR outperforms existing methods and constitutes a
widely applicable alternative to available tools.

## Requirements:

CopywriteR was developed for UNIX based systems (including OSX) and requires the following command line tools:

- Samtools (http://samtools.sourceforge.net/) - To test Samtools: `$ samtools`
- Bedtools (http://bedtools.readthedocs.org/) - To test Bedtools: `$ bedtools --version`
- MACS 1.4 (http://liulab.dfci.harvard.edu/MACS/). To test MACS: `$ macs14 --version`

The following R-packages are required to run CopywriteR:

- Rsamtools
- CGHcall
- snowfall
- IRanges
- matrixStats
- data.table
- gtools

Some of these can be installed from bioconductor.org:

    > source("http://bioconductor.org/biocLite.R")
    > biocLite(c('Rsamtools', 'CGHcall', 'snowfall', 'IRanges'))

The remaining R-packages are available through CRAN:

    > install.packages(c('matrixStats', 'data.table', 'gtools'))

## Installation R-package:

After installing the required tools as described above you can download the pre-compiled CopywriteR R-package.
The package can be installed from the command line using the following command:

    $ R CMD INSTALL CopywriteR*.tar.gz

## CopywriteR usage:

Load the CopywriteR package in R using:

    > library("CopywriteR")

CopywriteR contains three main functions:

preCopywriteR will generate bin, mapability, GC-content, blacklist and capture regions .bed files for any specified bin size and for available reference genomes.
For every combination of reference genome and bin size, a separate set of such files needs to be created.
preCopywriteR takes pre-assembled 1kb bin mapability, GC-content and blacklist .bed files as input.
Available reference genomes are hg19, mm9 and mm10.
These can be downloaded from the release section.

    preCopywriteR(blackGCMapaFolder, outputFolder, binSize, reference)

CopywriteR will generate separate tables with compensated read counts and normalized compensated read counts (after correction for GC-content, mapability and removal of blacklisted regions).

    CopywriteR(bamFolder, destinationFolder, referenceFolder, whichControl, ncpu, captureRegionsBedFile)

plotCNA performs segmentation, calling and plotting of copy number profiles using the CGHcall package.

    plotCNA(destinationFolder, set.nchrom)

For more details see R-package manual.
Alternatively, one of the following commands can be used to show help files for the corresponding function:

    > ?preCopywriteR
    > ?CopywriteR
    > ?plotCNA

## Troubleshooting

There are a number of requirements for your CopywriteR analysis to run successfully. These are discussed below.

### Chromosome names

CopywriteR by default assumes that the chromosome names in .bam files are "1", "2", ... "X", and "Y".
These chromosome names are incorporated in the bin, mapability, GC-content, blacklist and capture regions .bed files by preCopywriteR.
CopywriteR can also be applied to .bam files with different chromosome names.
In this case, the supporting .bed files need to have the same chromosome notation.
In case of non-matching chromosome notations, an error message will be displayed and the analysis will be aborted.
Non-matching chromosome names between .bam and supporting .bed files can be matched by changing the supporting .bed files, for instance in UNIX using awk as follows:

    > awk 'BEGIN { FS = OFS = "\t"; } { print "chr"$1, $2, $3; }' bins.bed > ../NEWFOLDER/bins.bed
    > awk 'BEGIN { FS = OFS = "\t"; } { print "chr"$1, $2, $3; }' blacklist.bed > ../NEWFOLDER/blacklist.bed
    > awk 'BEGIN { FS = OFS = "\t"; } { print "chr"$1, $2, $3, $4, $5; }' GC_content.bed > ../NEWFOLDER/GC_content.bed
    > awk 'BEGIN { FS = OFS = "\t"; } { print "chr"$1, $2, $3, $4; }' mapability.bed > ../NEWFOLDER/mapability.bed

Alternatively, chromosome names in .bam files can be adjusted using bedtools as follows (UNIX only):

    $ samtools view -H in.bam | awk 'BEGIN { FS = OFS = "\t"; } {if ($1 == "@SQ") { gsub("SN:chr", "SN:", $2); print $1, $2, $3; } else print; }' | samtools reheader - in.bam > out.bam

Please note that gsub in awk works similar to the gsub command in R.
Also, changing the chromosome names in the .bam header is sufficient as the chromosome names in the body of the file are in fact references to the chromosome names in the header.

### Number of chromosomes

CGHcall fails to run on the Y-chromosome when it has too few data points. If this occurs, we recommend a re-run of plotCNA with set.nchrom set as the total amount of chromosomes minus 1 (i.e., 23 for the human genome). CGHcall will then ignore the Y-chromosome.

## Contact

We have tried to make the CopywriteR code readable and its use as easy as possible. If any questions arise regarding the package, or if you want to report any bugs, please do not hesitate and contact:

- [Thomas Kuilman](mailto:t.kuilman@nki.nl)
- [Oscar Krijgsman](mailto:o.krijgsman@nki.nl)

Thomas and Oscar are working in the laboratory of Prof. Dr. Daniel S. Peeper.

- [Lab website](http://research.nki.nl/peeperlab/)


## Reported bugs

- None

## Changes and additions we are currently working on

- [ ] Address all NOTES when building package
- [ ] Solve problem with 0-based ranges
- [ ] Remove blacklisted regions by excluding regions and not excluding bins
- [ ] Check problem with not using all bams as samples (RB)
- [ ] Include helper file folder in log-file
- [ ] Change CopywriteR to allow custom bins?
- [ ] Compile into bioConductor package
- [x] Solve issue with data RB (implemented in source code)
- [x] Remove potential to overwrite existing folders and files (implemented in source code)
- [x] Change handling of NAs in plotCNA function (implemented in source code)
- [x] Make names consistent and apply Rlint to code (implemented in source code)
- [x] Use BiocParallel for parallel computing (implemented in source code)
- [x] Check why sex chromosomes are not plotted (solved in source code)
- [x] Add track line to .igv output file (#track viewLimits=-3:3 graphType=heatmap color=255,0,0) (implemented in source code)
- [x] Check log.txt output file (implemented in source code)
- [x] Check that all chromosomes have same length (implemented in source code)
- [x] Remove possibility for error when allowing more cpus than samples (implemented in source code)
- [x] Check for platform-compatibility when handling sample.files variable (implemented in source code)
- [x] Check potential plotting error when running on custom bin files (implemented in source code)
- [x] Remove dependency for MACS 1.4 and use chipseq package instead (implemented in source code)
- [x] Change -Inf in log2-table to small value for visualization in IGV (implemented in source code)
- [x] Provide the amount of sequence reads used by CopywriteR (implemented in source code)
- [x] Check rm() function if no capture region file is specified (implemented in source code)
- [x] Provide the option to keep intermediary files; remove by default (implemented in source code)
- [x] Check writing permissions on bam.folder (implemented in source code)
- [x] Include check for presence .bai files (implemented in source code)
- [x] Use file.paths for all file paths (implemented in source code)
- [x] Change output format for raw data from .txt to .igv (implemented in source code)
- [x] Remove bug in determining the total chromosome number with missing sex chromosomes (implemented in source code)
- [x] Remove dependency for CGHcall and use only DNAcopy for segmentation (implemented in source code)
- [x] Provide plotting function for DNAcopy output (implemented in source code)
- [x] Make keeping intermediate .bam and MACS files optional (implemented in source code)
- [x] Change input structure to allow more flexibility (implemented in source code)
- [x] Clean up code (implemented in source code)
- [x] Remove dependency for bedtools (replaced by GenomicRanges in source code)
- [x] Remove dependency for Samtools (replaced by Rsamtools in source code)
- [x] Include a check for integrity of BAM files (implemented in source code)
- [x] Change output format to igv2.0 (implemented in source code)
- [x] Remove bug in determining chromosome name prefixes (fixed in source code)
- [x] Debug error when using 1kb bin size
- [x] Support alternative chromosome names (i.e., "chr1" instead of "1")
- [x] Extract binSize from bins.bed file
- [x] Support for relative path names
- [x] Remove requirement for a trailing `/` in folder path names
- [x] Add 1kb mapability and GC-content files for mouse mm9 and mm10 genomes
- [x] Make captureRegionsBedFile optional
- [x] Allow processing of single-end sequences
- [x] Increase speed for generating bins in `preCopywriteR`