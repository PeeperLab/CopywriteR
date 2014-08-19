ENCODER
=======

Current methods for the detection of copy number aberrations (CNAs) from whole-exome sequencing (WES) data use the depth of coverage of captured exons only. Accurate CNA determination is complicated by the uneven distribution of exons throughout the genome and by the non-uniform sequence capture. Therefore, we have developed ENCODER (ENhanced COpy number Detection from Exome Reads), which eludes these problems by exploiting the ‘off-target’ sequence reads. ENCODER allows the extraction of uniformly distributed copy number information, and outperforms methods based on exonic sequence read counts, particularly on samples of low quality.

# Requirements:

ENCODER was developed for UNIX based systems (including OSX) and requires the following tools to be installed on your system: 

- Samtools (http://samtools.sourceforge.net/). To test Samtools: `$ samtools`

- Bedtools (http://bedtools.readthedocs.org/). To test Bedtools `$ bedtools --version`

- MACS 1.4 (http://liulab.dfci.harvard.edu/MACS/). To test MACS `$ macs14 --version`

- Multiple R-packages available from bioconductor.org. 
 Executing the following code in R will install or update the required packages: 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`> source("http://bioconductor.org/biocLite.R")` 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`> biocLite(c('Rsamtools', 'CGHcall', 'snowfall', 'IRanges'))` 

- Additional R-packages are available through CRAN.
 Executing the following code in R will install or update the remaining required packages: 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`> install.packages(c('matrixStats', 'data.table', 'gtools'))` 

# Installation R-package:

After installing the required tools as described above you can download the pre-compiled ENCODER R-package and annotation files. The package can be installed from the command line using the following command:

&nbsp;&nbsp;&nbsp;`$ R CMD INSTALL ENCODER.tar.gz` 


# ENCODER usage:

Start R and load the ENCODER package using: 

&nbsp;&nbsp;&nbsp;`> library("ENCODER")`

ENCODER contains three main functions: 

- preENCODER  
The function `preENCODER(MAPA_GC_location, outputFolder, binSize, reference)` will generate reference files for mapability, GC-content and bin files. 
This function should be run for each combination of reference (e.g. hg19, mm10) and binSize. 

- ENCODER  
The function `ENCODER(bamFolder, destinationFolder, referenceFolder, whichControl, captureRegionsBedFile, ncpu)` will generate separate tables with compensated read counts and normalized compensated read counts (after correction for GC-content, mapability and removal of blacklisted regions).

- CNAprofiles  
The function `CNAprofile(destinationFolder)` performs segmentation, calling and plotting of copy number profiles using the CGHcall package.

For more details see R-package. `> ?preENCODER`, `> ?ENCODER`, and `> ?CNAprofile`  in R will show help files and descriptions for each of the functions. 


# Contact

We have tried to make the ENCODER code readable and its use as easy as possible. If any questions arise regarding the package, or if you want to report any bugs, please do not hesitate and contact:

- Thomas Kuilman (t.kuilman@nki.nl)
- Oscar Krijgsman (o.krijgsman@nki.nl)

Thomas and Oscar are working in the laboratory of Prof. Dr. Daniel S. Peeper.

- Lab website (http://research.nki.nl/peeperlab/)

# Reported bugs

Major bugs

- 140814 - None...


Minor bugs

- 140814 - Destination folder needs complete path, `./` alone does not work. 
- 140814 - Currently all folder paths in `> ENCODER()` need a trailing `/`. 

# Changes and additions we are currently working on

- [ ] Extract binSize from bins.bed file
- [ ] Different input structure to indicate which bam files should be used as controls
- [ ] Change from MACS 1.4 to other ChIP seq tool available in R (chipseq from bioconductor?)
- [ ] Change from Samtools to Rsamtools
- [ ] Remove all unix specific functions and commands
- [ ] Clean up code
- [ ] Compile into bioConductor package
- [x] Make captureRegionsBedFile optional
- [x] Allow processing of single-end sequences
- [x] Increase speed for generating bins in `preENCODER`






