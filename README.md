ENCODER
=======

ENCODER: off-target sequence reads for DNA copy number detection by whole-exome sequencing

Current methods for the detection of copy number aberrations (CNA) from whole-exome sequencing (WES) data use the depth of coverage of captured exons only. Accurate CNA determination is complicated by the uneven distribution of exons throughout the genome and non-uniform sequence capture. We developed ENCODER (ENhanced COpy number Detection from Exome Reads), which eludes these problems by exploiting the ‘off-target’ sequence reads. ENCODER allows the extraction of uniformly distributed copy number information. In addition, ENCODER outperforms methods based on exonic sequence read counts, particularly on samples of low quality.


# Requirements:

ENCODER was developed for UNIX based systems (including MAC OSX) and requires the following tools to be installed on your system: 

- Samtools (http://samtools.sourceforge.net/). To test Samtools: `$ samtools`

- Bedtools (http://bedtools.readthedocs.org/). To test Bedtools `$ bedtools --version`

- MACS 1.4 (http://liulab.dfci.harvard.edu/MACS/). To test MACS `$ macs14 --version`

- Multiple R-packages available from bioconductor.org. 
 Executing the following code in R will install or update the required packages: 

`> source("http://bioconductor.org/biocLite.R")` 

`> biocLite(c('Rsamtools', 'CGHcall', 'snowfall', 'doParallel', 'IRanges'))` 

`> update.packages(repos=biocinstallRepos(), ask=FALSE)` 


# Installation R-package:

After installing the required tools as described above you can download the pre-compiled ENCODER R-package and annotation files. 

The pre-compiled ENCODER R-package can be installed from the command line with the following command: `$ R CMD INSTALL ENCODER.tar.gz`. 



# ENCODER usage:

- preENCODER
The function `preENCODER()` 

- ENCODER
The function `ENCODER()`

- CNAprofiles
The function `CNAprofile()`


# Contact

Developers:
- Thomas Kuilman (t.kuilman@nki.nl)
- Oscar Krijgsman (o.krijgsman@nki.nl)
Peeper group

Division of Molecular Oncology (H5)

The Netherlands Cancer Institute

Plesmanlaan 121

1066 CX Amsterdam

The Netherlands








