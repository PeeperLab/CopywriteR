############################################################
######  Pipeline for copy number alteration analysis  ######
############################################################

## Dependencies:
## Samtools
## MACS1.4.2
## Bedtools
## Rsamtools
## SeqCNV.R
## snowfall
## CGHcall
## doParallel (and its dependencies)
## IRanges

## Use sorted bed-files: sort -k1,1 -k2,2n in.bed > in.sorted.bed
## Optional: captureRegionsBedFile
## Eindig met Stanford file (ID, chr, start, stop)
## optional: extended output
## Std output
## Does not work with 1 sample?

## To do list ENCODER R-package
## Change from doParallel to snowFall

input <- inputStructure("./", "./", c(1,2,3,4,5,6,1,2,3,4,5,6), "hg19", "/srv/nfs4/medoid-bulk/NKI/t.kuilman/Resources/Agilent_SureSelect_All_exon_50M_hg19_edit.bed", 12)

input <- inputStructure("/srv/nfs4/paranoid-home/NKI/t.kuilman/PDX_cellines_pipeline/BAMfiles/", "/srv/nfs4/paranoid-home/NKI/t.kuilman/PDX_cellines_pipeline/", c(1,2,3,4,5,6,1,2,3,4,5,6), "hg19", "/srv/nfs4/medoid-bulk/NKI/t.kuilman/Resources/Agilent_SureSelect_All_exon_50M_hg19_edit.bed", 12)

inputStructure <- function(bamfolder, destinationfolder, whichControl, reference, captureRegionsBedFile, ncpu) {
	if(!(reference == "hg19" || reference == "mm9" || reference == "mm10")) {
		stop("Unknown reference")
	}
	system(paste0("mkdir ", destinationfolder, "/CNAprofiles/"))
	bam_list <- list.files(path = bamfolder, pattern = ".bam$")
	sink(file = paste0(destinationfolder, "CNAprofiles/log.txt"), type = c("output", "message"))
	for(i in 1:length(bam_list)) {
		cat("The control for sample", bam_list[i], "will be", bam_list[whichControl[i]], "\n")
	}
	cat("The reference for this analysis is", reference, "\n")
	cat("The capture region file is", captureRegionsBedFile, "\n")
	cat("This analysis will be run on", ncpu, "cpus", "\n\n\n")
	sink()
	for(i in 1:length(bam_list)) {
		cat("The control for sample", bam_list[i], "will be", bam_list[whichControl[i]], "\n")
	}
	cat("The reference for this analysis is", reference, "\n")
	cat("The capture region file is", captureRegionsBedFile, "\n")
	cat("This analysis will be run on", ncpu, "cpus", "\n\n\n")
	return(list(bamfolder = bamfolder, destinationfolder = destinationfolder, whichControl = whichControl, reference = reference, captureRegionsBedFile = captureRegionsBedFile, ncpu = ncpu))
}




# ENCODER <- function(inputFile, outputFolder, reference, binSize, bedfile, ncpu){

ENCODER <- function(inputStructure) {
	BINSIZE = 20000
	
	start_time <- Sys.time()
	sink(file = paste0(inputStructure$destinationfolder, "CNAprofiles/log.txt"), append = TRUE, type = c("output", "message"))
	options(width = 150)
	system(paste0("mkdir ", inputStructure$destinationfolder, "/CNAprofiles/BamBaiMacsFiles/"))

	# Load library for parallel computing
	library(snowfall)
	
	# Set the reference-dependent files
	seqCNVFile <- "/srv/nfs4/medoid-bulk/NKI/t.kuilman/R_scripts_others/gcmappa/SeqCNV/SeqCNV.R"
	captureRegionsBedFile <- inputStructure$captureRegionsBedFile
	if(inputStructure$reference == "hg19") {
		windowBedFile <- "/srv/nfs4/medoid-bulk/NKI/t.kuilman/Resources/Genome_bin_0_20000_avelds.bed"
		gcContentBedFile <- "/srv/nfs4/medoid-bulk/NKI/t.kuilman/R_scripts_others/gcmappa/gccontent-20000-hg19.txt"
		mappabilityBedFile <- "/srv/nfs4/medoid-bulk/NKI/t.kuilman/R_scripts_others/gcmappa/mapa51-20K-hg19.bed"
		blacklistBedFile <- "/srv/nfs4/medoid-bulk/NKI/t.kuilman/R_scripts_others/gcmappa/hg19-blacklist-nochr.bed"
		nchrom = 24
	}
	if(inputStructure$reference == "mm9") {
		windowBedFile <- "/srv/nfs4/medoid-bulk/NKI/t.kuilman/Resources/Genome_bin_0_20000_avelds_mm9.bed"
		gcContentBedFile <- "/srv/nfs4/medoid-bulk/NKI/t.kuilman/R_scripts_others/gcmappa/gccontent-20000-mm9.txt"
		mappabilityBedFile <- "/srv/nfs4/medoid-bulk/NKI/t.kuilman/R_scripts_others/gcmappa/mapa51-20K-mm9.bed"
		blacklistBedFile <- "/srv/nfs4/medoid-bulk/NKI/t.kuilman/R_scripts_others/gcmappa/mm9-blacklist-nochr.bed"
		nchrom = 21
	}
	if(inputStructure$reference == "mm10") {
		windowBedFile <- "/srv/nfs4/medoid-bulk/NKI/t.kuilman/Resources/Genome_bin_0_20000_avelds_mm10.bed"
		gcContentBedFile <- "/srv/nfs4/medoid-bulk/NKI/t.kuilman/R_scripts_others/gcmappa/gccontent-20000-mm10.txt"
		mappabilityBedFile <- "/srv/nfs4/medoid-bulk/NKI/t.kuilman/R_scripts_others/gcmappa/mapa51-20K-mm10.bed"
		blacklistBedFile <- "/srv/nfs4/medoid-bulk/NKI/t.kuilman/R_scripts_others/gcmappa/mm10-blacklist-nochr.bed"
		nchrom = 21
	}
	
	# Create list of .bam files
	setwd(inputStructure$bamfolder)
	bam_list <- list.files(path = inputStructure$bamfolder, pattern = ".bam$")
	cat(bam_list, "\n", sep = "\n")

	# Index .bam files
	ibam <- function(bam_list) {
		system(paste("samtools index", bam_list))
		paste("samtools index", bam_list)
	}
	sfInit(parallel=TRUE, cpus = inputStructure$ncpu)
	toLog <- sfLapply(bam_list, ibam)
	sfStop()
	cat(unlist(toLog), "\n", sep = "\n")

	# Remove anomalous reads and reads with Phred < 37
	properreads <- function(bam_list, inputStructure) {	
		system(paste0("samtools view -b -f 2 -q 37 ", bam_list, " > ", inputStructure$destinationfolder, "CNAprofiles/BamBaiMacsFiles/", gsub(".bam$", "_properreads.bam", bam_list)))
		paste0("samtools view -b -f 2 -q 37 ", bam_list, " > ", inputStructure$destinationfolder, "CNAprofiles/BamBaiMacsFiles/", gsub(".bam$", "_properreads.bam", bam_list))
	}
	sfInit(parallel=TRUE, cpus = inputStructure$ncpu)
	toLog <- sfLapply(bam_list, properreads, inputStructure)
	sfStop()
	cat(unlist(toLog), "\n", sep = "\n")

	################
	statistics <- matrix(nrow = length(bam_list), ncol = 7, dimnames = list(paste(bam_list), c("Total", "Total properpair", "Unmapable", "Mitochondrion", "All chromosomes", "Rest", "In peaks")))
	i <- c(1:length(bam_list))
	stats <- function(i, bam_list) {
		as.numeric(system(paste0("samtools view -c ", bam_list[i]), intern = TRUE))
	}
	sfInit(parallel=TRUE, cpus = inputStructure$ncpu)
	res <- sfSapply(i, stats, bam_list)
	sfStop()
	statistics[,1] <- res
	################
	
	# Create new .bam list
	setwd(paste0(inputStructure$destinationfolder, "CNAprofiles/BamBaiMacsFiles/"))
	bam_list <- list.files(path = paste0(inputStructure$destinationfolder, "CNAprofiles/BamBaiMacsFiles/"), pattern = "_properreads.bam$")
	cat(bam_list, "\n", sep = "\n")

	# Index _properreads.bam files
	iproperreads <- function(bam_list) {
		system(paste("samtools index", bam_list))
		paste("samtools index", bam_list)
	}
	sfInit(parallel=TRUE, cpus = inputStructure$ncpu)
	toLog <- sfLapply(bam_list, iproperreads)
	sfStop()
	cat(unlist(toLog), "\n", sep = "\n")

	################
	for (i in 1:length(bam_list)) {
		table <- as.matrix(system(paste0("samtools idxstats ", bam_list[i]), intern = TRUE), nrow = 26, ncol = 4)
		table <- strsplit(table, "\t")
		table <- do.call(rbind, table)
		MIT <- as.integer(table[grep("M", table[,1]),3])
		CHR <- sum(as.integer(table[,3])) - MIT
		UNMAP <- sum(as.integer(table[,4]))
		statistics[i,2] <- CHR + MIT + UNMAP
		statistics[i,3] <- UNMAP
		statistics[i,4] <- MIT
		statistics[i,5] <- CHR
	}
	################
	
	# Create list with numbers of controls
	controlNumbers <- unique(inputStructure$whichControl)
	
	# Call peaks in .bam file of control sample
	i <- c(1:length(controlNumbers))	
	macs14 <- function(i, controlNumbers, bam_list) {
		system(paste0("macs14 -t " , bam_list[controlNumbers[i]], " -n MACS", controlNumbers[i], " -g hs --nolambda"))
		paste0("macs14 -t " , bam_list[controlNumbers[i]], " -n MACS", controlNumbers[i], " -g hs --nolambda")
	}
	sfInit(parallel=TRUE, cpus = min(length(controlNumbers), inputStructure$ncpu))
	toLog <- sfLapply(i, macs14, controlNumbers, bam_list)
	sfStop()
	cat(unlist(toLog), "\n", sep = "\n")
	
	# Remove peak-regions from .bam files
	i <- c(1:length(bam_list))
	peakrm <- function(i, bam_list, inputStructure) {
		system(paste0("bedtools intersect -abam ", bam_list[i], " -b MACS", inputStructure$whichControl[i], "_peaks.bed -v > ", gsub(".bam$", "_peakrm.bam", bam_list[i])))
		paste0("bedtools intersect -abam ", bam_list[i], " -b MACS", inputStructure$whichControl[i], "_peaks.bed -v > ", gsub(".bam$", "_peakrm.bam", bam_list[i]))
	}
	sfInit(parallel=TRUE, cpus = inputStructure$ncpu)
	toLog <- sfSapply(i, peakrm, bam_list, inputStructure)
	sfStop()
	cat(unlist(toLog), "\n", sep = "\n")

	# Create new .bam list
	bam_list <- list.files(path = paste0(inputStructure$destinationfolder, "CNAprofiles/BamBaiMacsFiles/"), pattern = "_peakrm.bam$")
	cat(bam_list, "\n", sep = "\n")

	# Index _peakrm.bam files
	i <- c(1:length(bam_list))
	ipeakrm <- function(i, bam_list) {
		system(paste("samtools index", bam_list[i]))
		paste("samtools index", bam_list[i])
	}
	sfInit(parallel=TRUE, cpus = inputStructure$ncpu)
	toLog <- sfSapply(i, ipeakrm, bam_list)
	sfStop()
	cat(unlist(toLog), "\n", sep = "\n")

	################
	for (i in 1:length(bam_list)) {
		table <- as.matrix(system(paste0("samtools idxstats ", bam_list[i]), intern = TRUE), nrow = 26, ncol = 4)
		table <- strsplit(table, "\t")
		table <- do.call(rbind, table)
		MIT <- as.integer(table[grep("M", table[,1]),3])
		CHR <- sum(as.integer(table[,3])) - MIT
		statistics[i,6] <- CHR
		statistics[i,7] <- as.integer(statistics[i,5]) - CHR
	}
	################
	print(statistics)
	cat("\n\n")
	
	# Create read_count matrix
	library(Rsamtools)
	
	bed <- read.table(file = windowBedFile, sep = "\t") ######
	read_count <- matrix(data = 0, ncol = 4, nrow = nrow(bed))
	read_count[,1] <- paste(bed[,1], paste(bed[,2], bed[,3], sep = "-"), sep = ":")
	read_count[,2] <- gsub("chr", "", paste(bed[,1]))
	read_count[,3] <- paste(bed[,2])
	read_count[,4] <- paste(bed[,3])

	# Calculate the number of reads per bin
	i <- c(1:length(bam_list))
	scanbam <- function(i, bam_list, bed) {
		library(Rsamtools)
		which <- RangedData(space = bed[,1], IRanges(bed[,2], bed[,3]))
		param <- ScanBamParam(which = which, what = c("pos"))
		bamreads <- scanBam(file = paste(bam_list[i]), param = param)
		readmap <- matrix(data = 0, ncol = 1, nrow = nrow(bed))
		for (j in 1:length(bamreads)) {
			readmap[j] <- length(unlist(bamreads[j]))
		}
		return(list(readmap, paste0("Rsamtools finished calculating reads per bin in sample ", i, " out of ", length(bam_list), "; number of bins = ", length(bamreads))))
	}
	sfInit(parallel=TRUE, cpus = length(bam_list))
	res <- sfSapply(i, scanbam, bam_list, bed)
	sfStop()
	for(i in seq(1,2*length(bam_list),2)) {
		read_count <- cbind(read_count, res[[i]])
	}
	cat(unlist(res[2,]), "\n", sep = "\n")
	sink()
	
	# Compensate for removal of reads in peak regions
	read_count <- cbind(read_count, read_count[,-c(1, 2, 3, 4)], read_count[,-c(1, 2, 3, 4)])
	for(controlNumber in controlNumbers) {	

		# Calculate the cumulative length of peak regions per bin
		intersection <- system(paste0("bedtools intersect -a ", windowBedFile, " -b ", inputStructure$destinationfolder,
			"CNAprofiles/BamBaiMacsFiles/MACS", controlNumber, "_peaks.bed -wao"), intern = TRUE)
		intersection <- strsplit(intersection, "\t")
		intersection <- do.call(rbind, intersection)
		intersection <- intersection[,c(1,2,3,9)]
		vec <- vector()
		for (i in 2:nrow(intersection)) {
			if(intersection[i,1] == intersection[i-1,1] && intersection[i,2] == intersection[i-1,2] && intersection[i,3] == intersection[i-1,3] ) {
				intersection[i,4] <- as.integer(intersection[i,4]) + as.integer(intersection[i-1,4])
				vec <- append(vec, i-1)
			}
		}
		intersection <- intersection[-vec,]

		# Check correspondence bins, calculate fraction of remaining bin length (after peak region removal), and calculate compensated read numbers
		controlFor <- grep(controlNumber, inputStructure$whichControl)
		for (i in 1:nrow(read_count)) {
			if(read_count[i,2] == intersection[i,1] && read_count[i,3] == intersection[i,2] && read_count[i,4] == intersection[i,3] ) {
				fraction_of_bin <- (BINSIZE-as.numeric(intersection[i,4])) / BINSIZE
				read_count[i,(4 + 2 * length(bam_list) + controlFor)] <- fraction_of_bin
				if(fraction_of_bin != 0) {
					read_count[i,(4 + controlFor)] <- as.numeric(read_count[i,(4 + length(bam_list) + controlFor)]) / fraction_of_bin
				}
				else {
					read_count[i,(4 + controlFor)] <- 0
				}
			}
			else {
				stop("Bins do not correspond")
			}
		}
	}
	colnames(read_count) <- c("BinID", "Chromosome", "StartPos", "StopPos", sub(pattern = "$", ".compensated", paste(bam_list)),
		paste(bam_list), sub(pattern = "$", ".fractionOfBin", paste(bam_list)))

	# Create output file
	write.table(read_count, file = paste0(inputStructure$destinationfolder, "CNAprofiles/", "read_counts_compensated.txt"), row.names = FALSE, col.names = TRUE, sep = "\t")

	# Create histograms of fraction_of_bin (fraction of length in bins (after peak region removal)
	system(paste0("mkdir ", inputStructure$destinationfolder, "CNAprofiles/qc/"))
	for(i in 1:length(bam_list)) {
		pdf(paste0(inputStructure$destinationfolder, "CNAprofiles/qc/fraction_of_bin_", i, ".pdf"), width=7, height=7)
		plot(ecdf(as.numeric(read_count[,4+(2*length(bam_list))+i])), verticals = TRUE, ylab = "Fraction of bins",
			xlab = "Remaining fraction of bin after peak removal", main = "Cumulative distribution of remaining bin fraction")
		dev.off()
	}

	###########################
	###  Script by a.velds  ###
	###########################
	
	require(iterators)
	require(foreach)
	library(doParallel)
	registerDoParallel(cores = inputStructure$ncpu)
	
	source(seqCNVFile)
	read_count <- read_count[,1:(4 + length(bam_list))]
	for(i in 5:ncol(read_count)) {
		write.table(read_count[,c(2,3,4,i)], colnames(read_count)[i], quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
	}
	f <- .findCovFiles("bam.compensated")
	
	gc <- read.delim(gcContentBedFile)[,c(1,2,3,5)]
	colnames(gc) <- c("chr","start","end","gc")
	mapa <- read.delim(mappabilityBedFile, header=FALSE, col.names=c("chr","start","end","mapa"))
	black <- read.delim(blacklistBedFile, header=F, col.names=c("chr", "start", "end"))
	
	data <- .loadCovData(f, gc=gc, mapa=mapa, black=black, excludechr="MT", datacol=4)
	#fix the names
	sampnames <- paste(sub("$","", f), paste0(round(colSums(data$cov) / 1e6,1),"M"), sep="_")
	colnames(data$cov) <- sampnames
	usepoints <- !(data$anno$chr %in% c("X","Y","MT", "chrX", "chrY", "chrM"))
	
	system.time(ratios <- foreach(x=iter(data$cov, by="col"), .combine="cbind") %dopar% {
		 .tng(data.frame(count=x[,1], gc=data$anno$gc, mapa=data$anno$mapa), use=usepoints, correctmapa=T, plot=paste0(inputStructure$destinationfolder, "CNAprofiles/qc/", colnames(x), ".png"))
	})
	
	colnames(ratios) <- sampnames
	rd <- list(ratios=ratios[!data$anno$black & !is.na(data$anno$mapa) & data$anno$mapa > .2, ], anno=data$anno[!data$anno$black & !is.na(data$anno$mapa) & data$anno$mapa > .2, ])
	
	sink(file = paste0(inputStructure$destinationfolder, "CNAprofiles/log.txt"), append = TRUE, type = c("output", "message"))
	rd2 <- rd
	for(i in 1:length(bam_list)) {
		rd2$ratios[,i] <- rd$ratios[,i] - rd$ratios[,inputStructure$whichControl[i]]
		cat("Relative log2-values are calculated for sample", colnames(rd2$ratios)[i], "with control", colnames(rd2$ratios)[inputStructure$whichControl[i]], "\n")
	}
	colnames(rd2$ratios) <- gsub("$", ".rel", colnames(rd2$ratios))
	cat("\n\n")

	###########################
	###  Script by a.velds  ###
	###########################
	
	# Create table with corrected log2 values and write to file
	read_count <- matrix(nrow = nrow(rd$ratios))
	read_count <- cbind(read_count, rd$anno[1:3])
	read_count[,1] <- paste(rd$anno[,1], paste(rd$anno[,2], rd$anno[,3], sep = "-"), sep = ":")
	read_count <- cbind(read_count, rd$ratios)
	read_count <- cbind(read_count, rd2$ratios)
	colnames(read_count) <- c("BinID", "Chromosome", "StartPos", "StopPos", colnames(read_count[5:length(colnames(read_count))]))
	
	read_count <- read_count[-which(rowSums(is.na(read_count[,-c(1:4)])) > 0),]
	read_count[,2] <- gsub("X", nchrom - 1, read_count[,2])
	read_count[,2] <- gsub("Y", nchrom, read_count[,2])#########################################
	read_count[,2] <- as.integer(read_count[,2])
	read_count[read_count == -Inf] <- -.Machine$integer.max/2
	read_count[read_count == Inf] <- .Machine$integer.max/2
	
	write.table(read_count, paste0(inputStructure$destinationfolder, "CNAprofiles/log2ratio_compensated_corrected.txt"), sep="\t", row.names=F, quote=F)

	# Calculate statistics of overlap of peaks and capture regions
	for(controlNumber in controlNumbers) {
		intersection <- system(paste0("bedtools intersect -a ", captureRegionsBedFile, " -b ", inputStructure$destinationfolder, "CNAprofiles/BamBaiMacsFiles/MACS", controlNumber, "_peaks.bed -wao"), intern = TRUE)
		intersection <- strsplit(intersection, "\t")
		intersection <- do.call(rbind, intersection)
		intersection <- intersection[,c(1,2,3,9)]
		vec <- vector()
		for (i in 2:nrow(intersection)) {
			if(intersection[i,1] == intersection[i-1,1] && intersection[i,2] == intersection[i-1,2] && intersection[i,3] == intersection[i-1,3] ) {
				intersection[i,4] <- as.integer(intersection[i,4]) + as.integer(intersection[i-1,4])
				vec <- append(vec, i-1)
			}
		}
		intersection <- intersection[-vec,]
		
		cat("Number of exons covered by peaks in sample ", bam_list[controlNumber], ": ", sum(intersection[,4] != "0"), "\n")
	
		intersection <- system(paste0("bedtools intersect -a ", inputStructure$destinationfolder, "CNAprofiles/BamBaiMacsFiles/MACS", controlNumber, "_peaks.bed -b ", captureRegionsBedFile, " -wao"), intern = TRUE)
		intersection <- strsplit(intersection, "\t")
		intersection <- do.call(rbind, intersection)
		intersection <- intersection[,c(1,2,3,9)]
		vec <- vector()
		for (i in 2:nrow(intersection)) {
			if(intersection[i,1] == intersection[i-1,1] && intersection[i,2] == intersection[i-1,2] && intersection[i,3] == intersection[i-1,3] ) {
				intersection[i,4] <- as.integer(intersection[i,4]) + as.integer(intersection[i-1,4])
				vec <- append(vec, i-1)
			}
		}
		intersection <- intersection[-vec,]
		
		cat("Number of peaks covered by exons in sample", bam_list[controlNumber], ": ", sum(intersection[,4] != "0"), "\n")
		cat("Total number of exons in sample", bam_list[controlNumber], ": ", nrow(read.table(captureRegionsBedFile, sep = "\t")), "\n")
		cat("Total number of peaks", bam_list[controlNumber], ":", nrow(read.table(paste0(inputStructure$destinationfolder, "CNAprofiles/BamBaiMacsFiles/MACS", controlNumber, "_peaks.bed"), sep = "\t")), "\n\n")
		
	}
	sink()
	cat("Total calculation time: ", Sys.time() - start_time, "\n\n")
	save(inputStructure, file = paste0(inputStructure$destinationfolder, "CNAprofiles/input.Rdata"))

	quit(save = "yes")
}
