ENCODER <- function(bamFolder, destinationFolder, referenceFolder, whichControl, ncpu, captureRegionsBedFile) {
	
	start_time <- Sys.time()

	## Make folder path absolute
	bamFolder <- tools::file_path_as_absolute(bamFolder)
	destinationFolder <- tools::file_path_as_absolute(destinationFolder)
	referenceFolder <- tools::file_path_as_absolute(referenceFolder)
	
	## Make folder path independent of trailing /
	bamFolder <- paste(unlist(strsplit(gsub("/$", "", bamFolder), "/")), "", sep = "/", collapse = "")
	destinationFolder <- paste(unlist(strsplit(gsub("/$", "", destinationFolder), "/")), "", sep = "/", collapse = "")
	referenceFolder <- paste(unlist(strsplit(gsub("/$", "", referenceFolder), "/")), "", sep = "/", collapse = "")
	
	if(missing(captureRegionsBedFile)) {
		captureRegionsBedFile <- "not specified"
	}
	
	##############################################################
	## Generate inputStructure to run ENCODER and check folders ##
	##############################################################
	
	# Check all folders
	if(!file.exists(bamFolder)){
		stop("The bamFolder could not be found. Please change your bamFolder path.")
	}
	
	if(!file.exists(destinationFolder)){
		stop("The destinationFolder could not be found. Please change your destinationFolder path.")
	}
	
	if(!file.exists(referenceFolder)){
		stop("The referenceFolder could not be found. Please change your referenceFolder path or run `preENCODER` to generate the required folder with GC-content and mapability files for your desired bin size.")
	}
	
	if(!file.exists(captureRegionsBedFile) & captureRegionsBedFile != "not specified") {
		stop("The captureRegionsBedFile could not be found. Please change your captureRegionsBedFile path.")
	}
	
	windowBedFile <- paste0(referenceFolder, "bins.bed")
	blacklistBedFile <- paste0(referenceFolder, "blacklist.bed")
	gcContentBedFile <- paste0(referenceFolder, "GC_content.bed")
	mappabilityBedFile <- paste0(referenceFolder, "mapability.bed")

	bed <- read.table(file = windowBedFile, as.is = TRUE, sep = "\t")
	nchrom <- length(unique(bed$V1))
	binSize <- bed$V3[1]
	
	# List all input files and write to log
	dir.create(paste0(destinationFolder, "/CNAprofiles/"))
	
	# Create data.frame with bam-files and corresponding reference and index of reference
	bam_list <- list.files(path = bamFolder, pattern = ".bam$")
	control_list <- whichControl
	
	sink(file = paste0(destinationFolder, "CNAprofiles/log.txt"), type = c("output", "message"))
	for(i in 1:length(bam_list)) {
		cat("The control for sample", bam_list[i], "will be", bam_list[control_list[i]], "\n") #####
	}
	cat("The binSize for this analysis is", binSize, "\n")
	cat("The capture region file is", captureRegionsBedFile, "\n")
	cat("This analysis will be run on", ncpu, "cpus", "\n\n\n")
	sink()
	for(i in 1:length(bam_list)) {
		cat("The control for sample", bam_list[i], "will be", bam_list[control_list[i]], "\n")
	}
	cat("The capture region file is", captureRegionsBedFile, "\n")
	cat("This analysis will be run on", ncpu, "cpus", "\n\n\n")
	
						# Test for compatibilty chromosome names
						prefixes <- vector(mode = "character")
						for(bam in bam_list) {
							header <- scanBamHeader(paste0(bamFolder, bam))
							prefixes <- append(prefixes, names(header[[1]]$targets)[1])
						}
						if(!all(prefixes == prefixes[1])) {
							stop("The bam files have different chromosome names. Please adjust the bam-files accordingly.")
						} else {
							prefix <- prefixes[1]
						}
						## Add match checks for bin, mappability, GC-content, blacklist and capture regions .bed files
					
	sink(file = paste0(destinationFolder, "CNAprofiles/log.txt"), append = TRUE, type = c("output", "message"))
	options(width = 150)
	
	dir.create(paste0(destinationFolder, "CNAprofiles/BamBaiMacsFiles/"))
		
	###############################
	## Run the ENCODER algortihm ##
	###############################

	# Create list of .bam files
	setwd(bamFolder)
	cat(bam_list, "\n", sep = "\n")

	# Index .bam files
	ibam <- function(bam_list) {
		system(paste("samtools index", bam_list))
		paste("samtools index", bam_list)
	}
	sfInit(parallel=TRUE, cpus = ncpu)
	toLog <- sfLapply(bam_list, ibam)
	sfStop()
	cat(unlist(toLog), "\n", sep = "\n")
	
	# Check whether BAMs are paired-end
	numberpairedendreads <- function(bam_list) {
		system(paste0("samtools view -f 1 -c ", bam_list), intern = TRUE)
	}
	sfInit(parallel=TRUE, cpus = ncpu)
	pairedEnd <- sfLapply(bam_list, numberpairedendreads)
	sfStop()
	pairedEnd <- ifelse(unlist(pairedEnd) > 0, TRUE, FALSE)
	for(i in 1:length(bam_list)) {
		cat("Paired-end sequencing for sample ", bam_list[i], ": ", unlist(pairedEnd)[i], "\n", sep = "")
	}
	cat("\n\n")

	# Remove anomalous reads and reads with Phred < 37
	i <- c(1:length(bam_list))
	properreads <- function(i, bam_list, destinationFolder, pairedEnd) {
		if(pairedEnd[i]) {
			system(paste0("samtools view -b -f 2 -q 37 ", bam_list[i], " > ", destinationFolder, "CNAprofiles/BamBaiMacsFiles/", gsub(".bam$", "_properreads.bam", bam_list[i])))
			paste0("samtools view -b -f 2 -q 37 ", bam_list[i], " > ", destinationFolder, "CNAprofiles/BamBaiMacsFiles/", gsub(".bam$", "_properreads.bam", bam_list[i]))
		}
		else {
			system(paste0("samtools view -b -q 37 ", bam_list[i], " > ", destinationFolder, "CNAprofiles/BamBaiMacsFiles/", gsub(".bam$", "_properreads.bam", bam_list[i])))
			paste0("samtools view -b -q 37 ", bam_list[i], " > ", destinationFolder, "CNAprofiles/BamBaiMacsFiles/", gsub(".bam$", "_properreads.bam", bam_list[i]))
		}
	}
	sfInit(parallel=TRUE, cpus = ncpu)
	toLog <- sfLapply(i, properreads, bam_list, destinationFolder, pairedEnd)
	sfStop()
	cat(unlist(toLog), "\n", sep = "\n")

	################
	statistics <- matrix(nrow = length(bam_list), ncol = 7, dimnames = list(paste(bam_list), c("Total", "Total properpair", "Unmapable", "Mitochondrion", "All chromosomes", "Rest", "In peaks")))
	i <- c(1:length(bam_list))
	stats <- function(i, bam_list) {
		as.numeric(system(paste0("samtools view -c ", bam_list[i]), intern = TRUE))
	}
	sfInit(parallel=TRUE, cpus = ncpu)
	res <- sfSapply(i, stats, bam_list)
	sfStop()
	statistics[,1] <- res
	################
	
	# Create new .bam list
	setwd(paste0(destinationFolder, "CNAprofiles/BamBaiMacsFiles/"))
	bam_list <- gsub(".bam$", "_properreads.bam", bam_list)
	cat(bam_list, "\n", sep = "\n")

	# Index _properreads.bam files
	iproperreads <- function(bam_list) {
		system(paste("samtools index", bam_list))
		paste("samtools index", bam_list)
	}
	sfInit(parallel=TRUE, cpus = ncpu)
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
	controlNumbers <- unique(control_list) #####
	
	# Call peaks in .bam file of control sample
	macs14 <- function(controlNumbers, bam_list) {
		system(paste0("macs14 -t " , bam_list[controlNumbers], " -n MACS", controlNumbers, " -g hs --nolambda"))
		paste0("macs14 -t " , bam_list[controlNumbers], " -n MACS", controlNumbers, " -g hs --nolambda")
	}
	sfInit(parallel=TRUE, cpus = min(length(controlNumbers), ncpu))
	toLog <- sfLapply(controlNumbers , macs14, bam_list)
	sfStop()
	cat(unlist(toLog), "\n", sep = "\n")
	
	# Remove peak-regions from .bam files
	i <- 1:length(bam_list)
	peakrm <- function(i, bam_list, control_list) {
		system(paste0("bedtools intersect -abam ", bam_list[i], " -b MACS", control_list[i], "_peaks.bed -v > ", gsub(".bam$", "_peakrm.bam", bam_list[i])))
		paste0("bedtools intersect -abam ", bam_list[i], " -b MACS", control_list[i], "_peaks.bed -v > ", gsub(".bam$", "_peakrm.bam", bam_list[i]))
	}
	sfInit(parallel=TRUE, cpus = ncpu)
	toLog <- sfSapply(i, peakrm, bam_list, control_list)
	sfStop()
	cat(unlist(toLog), "\n", sep = "\n")

	# Create new .bam list
	bam_list <- gsub(".bam$", "_peakrm.bam", bam_list)
	cat(bam_list, "\n", sep = "\n")

	# Index _peakrm.bam files
	ipeakrm <- function(bam_list) {
		system(paste("samtools index", bam_list))
		paste("samtools index", bam_list)
	}
	sfInit(parallel=TRUE, cpus = ncpu)
	toLog <- sfSapply(bam_list, ipeakrm)
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
	
	bed <- read.table(file = windowBedFile, sep = "\t")
	read_count <- matrix(data = 0, ncol = 4, nrow = nrow(bed))
	read_count[,1] <- paste(bed[,1], paste(bed[,2], bed[,3], sep = "-"), sep = ":")
	read_count[,2] <- gsub("chr", "", paste(bed[,1]))
	read_count[,3] <- paste(bed[,2])
	read_count[,4] <- paste(bed[,3])

	# Calculate the number of reads per bin
	scanbam <- function(bam_list, bed) {
		library(Rsamtools)
		which <- RangedData(space = bed[,1], IRanges(bed[,2], bed[,3]))
		param <- ScanBamParam(which = which, what = c("pos"))
		bamreads <- scanBam(file = paste(bam_list), param = param)
		readmap <- matrix(data = 0, ncol = 1, nrow = nrow(bed))
		for (j in 1:length(bamreads)) {
			readmap[j] <- length(unlist(bamreads[j]))
		}
		return(list(readmap, paste0("Rsamtools finished calculating reads per bin in sample ", bam_list, "; number of bins = ", length(bamreads))))
	}
	sfInit(parallel=TRUE, cpus = ncpu)
	res <- sfSapply(bam_list, scanbam, bed)
	sfStop()
	for(i in seq(1,2*length(bam_list),2)) {
		read_count <- cbind(read_count, res[[i]])
	}
	cat(unlist(res[2,]), "\n", sep = "\n")
	sink()
	
	# Compensate for removal of reads in peak regions
	read_count <- cbind(read_count, read_count[,-c(1, 2, 3, 4)], read_count[,-c(1, 2, 3, 4)])
	
	for(controlNumber in controlNumbers) {

		# Calculate overlap of peaks with bins using bedtools
		intersection <- system(paste0("bedtools intersect -a ", windowBedFile, " -b ", destinationFolder,
			"CNAprofiles/BamBaiMacsFiles/MACS", controlNumber, "_peaks.bed -wao"), intern = TRUE)
		intersection <- strsplit(intersection, "\t")
		intersection <- as.data.frame(do.call(rbind, intersection))
		intersection <- intersection[,c(1,2,3,9)]
		colnames(intersection) <- c("Chromosome", "StartPos", "StopPos", "Overlap")
		intersection$Overlap <- as.numeric(as.character(intersection$Overlap))
		
		## Calculate the cumulative overlap using data.table
		intersection <- data.table(intersection)
		intersection <- intersection[,lapply(.SD, sum), by = c("Chromosome", "StartPos", "StopPos")]
		intersection <- as.data.frame(intersection)

		# Check correspondence bins, calculate fraction of remaining bin length (after peak region removal), and calculate compensated read numbers
		controlFor <- grep(controlNumber, whichControl)
		if(all(read_count[,2] == intersection[,1] & read_count[,3] == intersection[,2] & read_count[,4] == intersection[,3])) {
			fraction_of_bin <- (binSize-as.numeric(intersection[,4])) / binSize
			read_count[,(4 + 2 * length(bam_list) + controlFor)] <- fraction_of_bin
			for (i in 1:nrow(read_count)) {
				if(fraction_of_bin[i] != 0) {
					read_count[i,(4 + controlFor)] <- as.numeric(read_count[i,(4 + length(bam_list) + controlFor)]) / fraction_of_bin[i]
				}
				else {
					read_count[i,(4 + controlFor)] <- 0
				}
			}
		}
		else {
			stop("Bins do not correspond")
		}
	}
	
	colnames(read_count) <- c("BinID", "Chromosome", "StartPos", "StopPos", sub(pattern = "$", ".compensated", paste(bam_list)),
		paste(bam_list), sub(pattern = "$", ".fractionOfBin", paste(bam_list)))

	# Create output file
	write.table(read_count, file = paste0(destinationFolder, "CNAprofiles/", "read_counts_compensated.txt"), row.names = FALSE, col.names = TRUE, sep = "\t")

	# Create histograms of fraction_of_bin (fraction of length in bins (after peak region removal)
	dir.create(paste0(destinationFolder, "CNAprofiles/qc/"))
	for(i in 1:length(bam_list)) {
		pdf(paste0(destinationFolder, "CNAprofiles/qc/fraction_of_bin_", i, ".pdf"), width=7, height=7)
		plot(ecdf(as.numeric(read_count[,4+(2*length(bam_list))+i])), verticals = TRUE, ylab = "Fraction of bins",
			xlab = "Remaining fraction of bin after peak removal", main = "Cumulative distribution of remaining bin fraction")
		dev.off()
	}

	#############################################
	## Normalize for GC-content and mapability ##
	#############################################

	read_count <- read_count[,1:(4 + length(bam_list))]
	for(i in 5:ncol(read_count)) {
		write.table(read_count[,c(2,3,4,i)], colnames(read_count)[i], quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
	}
	f <- list.files(pattern = "bam.compensated")
	
	gc <- read.delim(gcContentBedFile)[,c(1,2,3,5)]
	colnames(gc) <- c("chr","start","end","gc")
	mapa <- read.delim(mappabilityBedFile, header=FALSE, col.names=c("chr","start","end","mapa"))
	black <- read.delim(blacklistBedFile, header=F, col.names=c("chr", "start", "end"))
	
	data <- .loadCovData(f, gc=gc, mapa=mapa, black=black, excludechr="MT", datacol=4)
	#fix the names
	sampnames <- paste(sub("$","", f), paste0(round(colSums(data$cov) / 1e6,1),"M"), sep="_")
	colnames(data$cov) <- sampnames
	usepoints <- !(data$anno$chr %in% c("X","Y","MT", "chrX", "chrY", "chrM"))
	
	# Perform normalization (in .tng helper function)
	i <- c(1:ncol(data$cov))
	normalizeRC <- function(i, data, .tng, usepoints, destinationFolder) {	
		.tng(data.frame(count = data$cov[,i], gc = data$anno$gc, mapa = data$anno$mapa), use = usepoints, correctmapa = TRUE, plot = paste0(destinationFolder, "CNAprofiles/qc/", colnames(data$cov)[i], ".png"))
	}
	sfInit(parallel=TRUE, cpus = ncpu)
	ratios <- sfLapply(i, normalizeRC, data, .tng, usepoints, destinationFolder)
	sfStop()
	ratios <- matrix(unlist(ratios), ncol = length(bam_list))
	
	colnames(ratios) <- sampnames
	
	rd <- list(ratios=ratios[!data$anno$black & !is.na(data$anno$mapa) & data$anno$mapa > .2, ], anno=data$anno[!data$anno$black & !is.na(data$anno$mapa) & data$anno$mapa > .2, ])
	
	sink(file = paste0(destinationFolder, "CNAprofiles/log.txt"), append = TRUE, type = c("output", "message"))
	rd2 <- rd
	for(i in 1:length(bam_list)) {
		rd2$ratios[,i] <- rd$ratios[,i] - rd$ratios[,whichControl[i]]
		cat("Relative log2-values are calculated for sample", colnames(rd2$ratios)[i], "with control", colnames(rd2$ratios)[whichControl[i]], "\n")
	}
	colnames(rd2$ratios) <- gsub("$", ".rel", colnames(rd2$ratios))
	cat("\n\n")

	###################
	## Create output ##
	###################
	
	# Create table with corrected log2 values and write to file
	read_count <- matrix(nrow = nrow(rd$ratios))
	read_count <- cbind(read_count, rd$anno[1:3])
	read_count[,1] <- paste(rd$anno[,1], paste(rd$anno[,2], rd$anno[,3], sep = "-"), sep = ":")
	read_count <- cbind(read_count, rd$ratios, rd2$ratios)
	colnames(read_count) <- c("BinID", "Chromosome", "StartPos", "StopPos", colnames(read_count[5:length(colnames(read_count))]))
	
	read_count <- read_count[-which(rowSums(is.na(read_count[,-c(1:4)])) > 0),]
	read_count[,2] <- gsub("X", nchrom - 1, read_count[,2])
	read_count[,2] <- gsub("Y", nchrom, read_count[,2])
	read_count[,2] <- as.integer(read_count[,2])
	read_count[read_count == -Inf] <- -.Machine$integer.max/2
	read_count[read_count == Inf] <- .Machine$integer.max/2
	
	write.table(read_count, paste0(destinationFolder, "CNAprofiles/log2ratio_compensated_corrected.txt"), sep="\t", row.names=F, quote=F)

	##############################################################################
	## Calculate overlap with captureRegionsBedFile for quality control purpose ##
	##############################################################################
	
	if(captureRegionsBedFile != "not specified") {
		for(controlNumber in controlNumbers) {
			intersection <- system(paste0("bedtools intersect -a ", captureRegionsBedFile, " -b ", destinationFolder, "CNAprofiles/BamBaiMacsFiles/MACS", controlNumber, "_peaks.bed -wao"), intern = TRUE)
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
		
			intersection <- system(paste0("bedtools intersect -a ", destinationFolder, "CNAprofiles/BamBaiMacsFiles/MACS", controlNumber, "_peaks.bed -b ", captureRegionsBedFile, " -wao"), intern = TRUE)
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
			cat("Total number of peaks", bam_list[controlNumber], ":", nrow(read.table(paste0(destinationFolder, "CNAprofiles/BamBaiMacsFiles/MACS", controlNumber, "_peaks.bed"), sep = "\t")), "\n\n")	
		}
	}
	
	sink()
	cat("Total calculation time: ", Sys.time() - start_time, "\n\n")
	
	inputStructure <- list(destinationFolder = destinationFolder, ncpu = ncpu, nchrom = nchrom)
	save(inputStructure, file = paste0(destinationFolder, "CNAprofiles/input.Rdata"))

}
