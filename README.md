# ENCODER

Current methods for the detection of copy number aberrations (CNAs) from whole-exome sequencing (WES) data use the depth of coverage of captured exons only.
Accurate CNA determination is complicated by the uneven distribution of exons throughout the genome and by the non-uniform sequence capture.
Therefore, we have developed ENCODER (ENhanced COpy number Detection from Exome Reads), which eludes these problems by exploiting the ‘off-target’ sequence reads.
ENCODER allows the extraction of uniformly distributed copy number information, and outperforms methods based on exonic sequence read counts, particularly on samples of low quality.

## Requirements:

ENCODER was developed for UNIX based systems (including OSX) and requires the following command line tools:

- Samtools (http://samtools.sourceforge.net/) - To test Samtools: `$ samtools`
- Bedtools (http://bedtools.readthedocs.org/) - To test Bedtools: `$ bedtools --version`
- MACS 1.4 (http://liulab.dfci.harvard.edu/MACS/). To test MACS: `$ macs14 --version`

The following R-packages are required to run ENCODER:

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

After installing the required tools as described above you can download the pre-compiled ENCODER R-package and annotation files.
The package can be installed from the command line using the following command:

    $ R CMD INSTALL ENCODER*.tar.gz

## ENCODER usage:

Load the ENCODER package in R using:

    > library("ENCODER")

ENCODER contains three main functions:

preENCODER will generate mapability and GC-content files for a specified bin size and reference genome.
For every combination of reference genome and bin size, a new set of mapability and GC-content files needs to be created.
preENCODER takes pre-assembled 1kb bin mapability and GC-content files as input.
These can be downloaded from the release section.
Available reference genomes are hg19, mm9 and mm10.

    preENCODER(MAPA_GC_location, outputFolder, binSize, reference)

ENCODER will generate separate tables with compensated read counts and normalized compensated read counts (after correction for GC-content, mapability and removal of blacklisted regions).

    ENCODER(bamFolder, destinationFolder, referenceFolder, whichControl, ncpu, captureRegionsBedFile)

plotCNA performs segmentation, calling and plotting of copy number profiles using the CGHcall package.

    plotCNA(destinationFolder)

For more details see R-package manual.
Alternatively, one of the following commands can be used to show help files for the corresponding function:

    > ?preENCODER
    > ?ENCODER
    > ?plotCNA

## Troubleshooting

There are a number of requirements for your ENCODER analysis to run successfully. These are discussed below.

### Chromosome names

ENCODER assumes by default that the chromosome names in .bam files are "1" through "20", "X", and "Y".
However, the bin, mappability, GC-content, blacklist and capture regions .bed files can be adjusted to match the notation in the .bam files, and ENCODER will also work.
ENCODER checks whether the same set of chromosome names have been used throughout all .bam files that are to be analyzed.
If they don't match, an error will be displayed and ENCODER aborts the analysis.
A similar error will be displayed if the chromosome names do not match those of the bin, mappability, GC-content, blacklist and capture regions .bed files.
Chromosome names in .bam files can be changed to the above-mentioned chromosome notation using bedtools as follows (UNIX only):

    $ samtools view -H in.bam | awk 'BEGIN { FS = OFS = "\t"; } {if ($1 == "@SQ") { gsub("SN:chr", "SN:", $2); print $1, $2, $3; } else print; }' | samtools reheader - in.bam > out.bam

Please note that gsub in awk works similar to the gsub command in R.
Also, changing the chromosome names in the .bam header is sufficient as the chromosome names in the body of the file are in fact references to the chromosome names in the header.

## Contact

We have tried to make the ENCODER code readable and its use as easy as possible. If any questions arise regarding the package, or if you want to report any bugs, please do not hesitate and contact:

- [Thomas Kuilman](mailto:t.kuilman@nki.nl)
- [Oscar Krijgsman](mailto:o.krijgsman@nki.nl)

Thomas and Oscar are working in the laboratory of Prof. Dr. Daniel S. Peeper.

- [Lab website](http://research.nki.nl/peeperlab/)


## Reported bugs

- None

## Changes and additions we are currently working on

- [ ] Make keeping intermediate .bam and MACS files optional
- [ ] Check warning message "Setting LC_CTYPE failed, using C" with some installations of R 
- [ ] Implement different input structure to indicate which bam files should be used as references
- [ ] Change from MACS 1.4 to other ChIP seq tool available in R (chipseq from bioconductor?)
- [ ] Change from Samtools to Rsamtools
- [ ] Remove all unix specific functions and commands
- [ ] Clean up code
- [ ] Compile into bioConductor package
- [ ] Make names consistent and apply Rlint to code
- [ ] Support alternative chromosome names (i.e., "chr1" instead of "1"; implemented in source code)
- [ ] Extract binSize from bins.bed file (implemented in source code)
- [ ] Support for relative path names (implemented in source code)
- [ ] Remove requirement for a trailing `/` in folder path names (implemented in source code)
- [x] Add 1kb mapability and GC-content files for mouse mm9 and mm10 genomes
- [x] Make captureRegionsBedFile optional
- [x] Allow processing of single-end sequences
- [x] Increase speed for generating bins in `preENCODER`