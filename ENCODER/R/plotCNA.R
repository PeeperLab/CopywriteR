plotCNA <- function(destinationFolder, set.nchrom = "determined.from.reference") {
	load(paste0(destinationFolder, "CNAprofiles/input.Rdata"), .GlobalEnv)
	
	if(set.nchrom != "determined.from.reference" & is.integer(set.nchrom)) {
		inputStructure$nchrom <- set.nchrom
	}

	# Read data
	read_count <- read.table(file = paste0(inputStructure$destinationFolder, "CNAprofiles/log2ratio_compensated_corrected.txt"), sep = "\t", header = TRUE, check.names = FALSE)

	# Run CGHcall
	raw <- make_cghRaw(read_count)
	prep <- preprocess(raw, maxmiss = 30, nchrom = inputStructure$nchrom)
	nor <-  normalize(prep, method = "median", smoothOutliers = TRUE)
	seg <-  segmentData(nor, method = "DNAcopy", nperm = 2000, undo.splits = "sdundo", min.width = 5, undo.SD = 1, clen = 25, relSDlong = 5)
	segnorm <- postsegnormalize(seg, inter = c(-1,1))
	listcalls <- CGHcall(segnorm, nclass = 5, robustsig = "yes", cellularity = 1, ncpus = inputStructure$ncpu)
	calls <- ExpandCGHcall(listcalls, segnorm, divide = 5, memeff = FALSE)
	save(calls, file = paste0(inputStructure$destinationFolder, "CNAprofiles/calls.Rdata"))
	
	# Make whole-genome plots
	system(paste0("mkdir ", inputStructure$destinationFolder, "CNAprofiles/Call_plots"))
	colnames_read_count <- gsub(" ", ".", colnames(read_count)[5:ncol(read_count)])
	for (i in 1:(ncol(read_count) - 4)) {
		setwd(paste0(inputStructure$destinationFolder, "CNAprofiles/Call_plots/"))
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