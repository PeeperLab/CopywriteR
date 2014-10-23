plotCNA <- function(destination.folder, set.nchrom = "determined.from.reference") {


#   read.counts <- read.counts[-which(rowSums(is.na(read.counts[, -c(1:4)])) > 0),]
#   read.counts[, 2] <- gsub(prefixes[1], "", read.counts[, 2])
#   read.counts[, 2] <- gsub("X", nchrom - 1, read.counts[, 2])
#   read.counts[, 2] <- gsub("Y", nchrom, read.counts[, 2])
#   read.counts[, 2] <- as.integer(read.counts[, 2])
#   read.counts[read.counts == -Inf] <- -.Machine$integer.max/2
#   read.counts[read.counts == Inf] <- .Machine$integer.max/2
# 



	## Make folder path absolute
	destination.folder <- tools::file_path_as_absolute(destination.folder)

	## Make folder path independent of trailing /
	destination.folder <- paste(unlist(strsplit(gsub("/$", "", destination.folder), "/")), "", sep = "/", collapse = "")
	
	# Check all folders
	if(!file.exists(destination.folder)){
		stop("The destination folder could not be found. Please change your destination.folder path.")
	}

	load(paste0(destination.folder, "input.Rdata"), .GlobalEnv)
	
	if(set.nchrom != "determined.from.reference") {
		inputStructure$nchrom <- set.nchrom
	}

	# Read data
	read_count <- read.table(file = paste0(inputStructure$destination.folder, "log2ratio_compensated_corrected.txt"), sep = "\t", header = TRUE, check.names = FALSE)

	# Run CGHcall
	raw <- make_cghRaw(read_count)
	prep <- preprocess(raw, maxmiss = 30, nchrom = inputStructure$nchrom)
	nor <-  normalize(prep, method = "median", smoothOutliers = TRUE)
	seg <-  segmentData(nor, method = "DNAcopy", nperm = 2000, undo.splits = "sdundo", min.width = 5, undo.SD = 1, clen = 25, relSDlong = 5)
	segnorm <- postsegnormalize(seg, inter = c(-1,1))
	listcalls <- CGHcall(segnorm, nclass = 5, robustsig = "yes", cellularity = 1, ncpus = inputStructure$ncpu)
	calls <- ExpandCGHcall(listcalls, segnorm, divide = 5, memeff = FALSE)
	save(calls, file = paste0(inputStructure$destination.folder, "calls.Rdata"))
	
	# Make whole-genome plots
	system(paste0("mkdir ", inputStructure$destination.folder, "Call_plots"))
	colnames_read_count <- gsub(" ", ".", colnames(read_count)[5:ncol(read_count)])
	for (i in 1:(ncol(read_count) - 4)) {
		setwd(paste0(inputStructure$destination.folder, "Call_plots/"))
		system(paste0("mkdir ",i, "_", colnames_read_count[i]))
		setwd(paste0(i, "_", colnames_read_count[i]))
		png(filename=paste0(i, "_freqonly_", colnames_read_count[i],".png"), width=2*480, height=480)
			plot(nor[,i], dotres = 1)
		dev.off()	
		pdf(file=paste0(i, "_freqonly_", colnames_read_count[i],".pdf"), width=2*7, height=7)
			plot(nor[,i], dotres = 1)
		dev.off()
		png(filename=paste0(i, "_", colnames_read_count[i],".png"), width=2*480, height=480)
			plot(calls[,i], dotres=1)
		dev.off()
		pdf(file=paste0(i, "_", colnames_read_count[i],".pdf"), width=2*7, height=7)
			plot(calls[,i], dotres=1)
		dev.off()
		png(filename=paste0(i, "_segmented_", colnames_read_count[i],".png"), width=2*480, height=480)
			plot(seg[,i], dotres=1)
		dev.off()
		pdf(file=paste0(i, "_segmented_", colnames_read_count[i],".pdf"), width=2*7, height=7)
			plot(seg[,i], dotres=1)
		dev.off()
			
		# Make per-genome plots
		for (j in 1:inputStructure$nchrom){
			print(j)
			png(filename=paste0("chromo_", j, "_", colnames_read_count[i],"_freqonly.png"), width=2*480, height=480)
				plot(nor[chromosomes(nor) == j,i], dotres=1)
			dev.off()
			pdf(file=paste0("chromo_", j, "_", colnames_read_count[i],"_freqonly.pdf"), width=2*7, height=7)
				plot(nor[chromosomes(nor) == j,i], dotres=1)
			dev.off()
			png(filename=paste0("chromo_", j, "_", colnames_read_count[i],"_seg.png"), width=2*480, height=480)
				plot(seg[chromosomes(seg) == j,i], dotres=1)
			dev.off()
			pdf(file=paste0("chromo_", j, "_", colnames_read_count[i],"_seg.pdf"), width=2*7, height=7)
				plot(seg[chromosomes(seg) == j,i], dotres=1)
			dev.off()
			png(filename=paste0("chromo_", j, "_", colnames_read_count[i],".png"), width=2*480, height=480)
				plot(calls[chromosomes(calls)==j,i], dotres=1)
			dev.off()
			pdf(file=paste0("chromo_", j, "_", colnames_read_count[i],".pdf"), width=2*7, height=7)
				plot(calls[chromosomes(calls)==j,i], dotres=1)
			dev.off()
		}
	}
}