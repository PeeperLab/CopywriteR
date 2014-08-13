### Generate annotation files (GC-content, mappability and bed files) for desired binSize

preENCODER<-function(MAPA_GC_location, outputFolder, binSize, reference){
	
	## Check MAPA_GC_location for 1kb files
	## For hg19 (human)
	if (reference=="hg19"){
		if(file.exists(paste(MAPA_GC_location, "/hg19_1kb/", sep=""))){
			cat("Reference folder hg19_1kb detected", "\n")
		}
	}
	## For mouse (mm10)
	if (reference=="mm10"){
		if(file.exists(paste(MAPA_GC_location, "/mm10_1kb/", sep=""))){
			cat("Reference folder mm10_1kb detected", "\n")
		}
	}
	if (!(reference == "hg19" || reference == "mm10")){
		stop("The reference is not recognised. Please provide a suitable reference (mm10 or hg19).")
	}
	
	## Check if binSize is factor of 1kb
	ifelse(!.is.wholenumber(binSize/10000), ){
		stop("Please provide a binSize which exceed a factor of 1000.")
	}
	
	## Generate files with desired binSize (MAPA, GC, bed, blacklist)
	
	

	mapa<-load(paste(MAPA_GC_location, "/hg19_1kb/mapa_1kb.rda", sep=""))
	
	
}





## To run
preENCODER("/Users/o.krijgsman/Documents/PostDoc/Projects/ENCODER_package/Annotation",	"/Users/o.krijgsman/Documents/PostDoc/Projects/ENCODER_package/testruns/", 20000, "hg19")
