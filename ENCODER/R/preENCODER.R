### Generate annotation files (GC-content, mappability and bed files) for desired binSize

preENCODER<-function(blackGCMapaFolder, outputFolder, binSize, reference){
	
	## Make folder path independent of trailing /
	blackGCMapaFolder <- paste(unlist(strsplit(gsub("/$", "", blackGCMapaFolder), "/")), "", sep = "/", collapse = "")
	outputFolder <- paste(unlist(strsplit(gsub("/$", "", outputFolder), "/")), "", sep = "/", collapse = "")
	
	## Checks
	# Check whether folders exist.
	if (!(reference == "hg19" || reference == "mm10" || reference == "mm9")){
		stop("The reference is not recognised. Please provide a suitable reference.")
	}

	if(file.exists(blackGCMapaFolder) == FALSE) {
		stop("The annotation folder could not be found. Please change your blackGCMapaFolder path.")
	} else {
		cat("Reference folder", blackGCMapaFolder, "detected", "\n")
	}

	if(file.exists(paste(outputFolder)) == FALSE) {
		stop("The outputFolder could not be found. Please change your outputPath.")
	} else {
		cat("Reference folder", outputFolder, "detected", "\n")
	}
		
	## Check if binSize is factor of 1kb
	if(!.is.wholenumber(binSize/1000)){
		stop("Please provide a binSize which is a multiple of 1000.")
	}
	
	## Generate files with desired binSize (MAPA, GC, bed, blacklist)
	# Load blacklist, CG-content and mapability files
	load(paste0(blackGCMapaFolder, "mapability.rda"))
	load(paste0(blackGCMapaFolder, "GCcontent.rda"))
	bed_file<-read.table(paste0(blackGCMapaFolder, "blacklist.bed"), as.is = TRUE, sep="\t")

	# Create bins with desired bin size
	MERGEBINNUMBER <- binSize/1000
	
	newBin <- NULL
	options(warn = -1)
	options(scipen = 999)
	for(chr in unique(mapa$chromosome)) {
		col2 <- colMins(matrix(mapa$start[mapa$chromosome == chr], nrow = MERGEBINNUMBER))
		col3 <- colMaxs(matrix(mapa$end[mapa$chromosome == chr], nrow = MERGEBINNUMBER))
		tmp <- cbind(chr, col2, col3)
		tmp <- tmp[1:(nrow(tmp)-1),]
		newBin <- rbind(newBin, tmp)
	}
	options(scipen = 0)
	options(warn = 0)
	
	cat("Generated", binSize, "bp bins for all chromosomes", "\n")
	
	# Create mapabillity file with desired bin size
	MERGEBINNUMBER <- binSize/1000

	newMapa <- NULL
	options(warn = -1)
	options(scipen = 999)
	for(chr in unique(mapa$chromosome)) {
		col2 <- colMins(matrix(mapa$start[mapa$chromosome == chr], nrow = MERGEBINNUMBER))
		col3 <- colMaxs(matrix(mapa$end[mapa$chromosome == chr], nrow = MERGEBINNUMBER))
		col4 <- colMeans(matrix(mapa$mapability[mapa$chromosome == chr], nrow = MERGEBINNUMBER))
		tmp <- cbind(chr, col2, col3, col4)
		tmp <- tmp[1:(nrow(tmp)-1),]
		newMapa <- rbind(newMapa, tmp)
	}
	options(scipen = 0)
	options(warn = 0)
	
	cat("Generated mapability file for bin size of", binSize,"bp", "\n")	
	
	# Create GC-content file with desired bin size
	newGC <- NULL
	options(warn = -1)
	options(scipen = 999)
	for(chr in unique(mapa$chromosome)) {
		col2 <- colMins(matrix(GC$start[GC$chromosome == chr], nrow = MERGEBINNUMBER))
		col3 <- colMaxs(matrix(GC$end[GC$chromosome == chr], nrow = MERGEBINNUMBER))
		col4 <- colMeans(matrix(GC$ATcontent[GC$chromosome == chr], nrow = MERGEBINNUMBER))
		col5 <- colMeans(matrix(GC$GCcontent[GC$chromosome == chr], nrow = MERGEBINNUMBER))
		tmp <- cbind(chr, col2, col3, col4, col5)
		tmp <- tmp[1:(nrow(tmp)-1),]
		newGC <- rbind(newGC, tmp)
	}
	options(scipen = 0)
	options(warn = 0)

	cat("Generated GC-content file for bin size of", binSize,"bp", "\n")		
	
	# Create folder for output files
	file_name <- paste0(reference, "_", binSize/1000, "kb/")
	dir.create(paste0(outputFolder, file_name))
	if(file.exists(paste0(outputFolder, file_name)) == FALSE) {
		stop("No output folder created, please check argument outputFolder and the corresponding folder permissions.")
	}

	## Write files to folder
	# Blacklist
	write.table(bed_file, file = paste0(outputFolder, file_name, "blacklist.bed"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
	# GC-content
	write.table(newGC, file = paste0(outputFolder, file_name, "GC_content.bed"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
	# Mapability
	write.table(newMapa, file = paste0(outputFolder, file_name, "mapability.bed"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")	
	# bed file with bins
	write.table(newBin, file = paste0(outputFolder, file_name, "bins.bed"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
	
}

