### Generate annotation files (GC-content, mappability and bed files) for desired binSize

preENCODER<-function(MAPA_GC_location, outputFolder, binSize, reference){
	
	## Check MAPA_GC_location for 1kb files
	if (reference=="hg19"){
		
	}
	
	if (reference=="mm10"){
		
	}
	if (!(reference == "hg19" || reference == "mm10"){
		stop("The reference is not recognised. Please provide a suitable reference (mm10 or hg19).")
	}
	
	
	## Check binSize
	ifelse(is.interger(binSize/10000), )
	
	## Generate files with desired binSize (MAPA, GC, bed, blacklist)
	
}





## To run
preENCODER("/Users/o.krijgsman/Documents/PostDoc/Projects/ENCODER_package/Annotation",	"/Users/o.krijgsman/Documents/PostDoc/Projects/ENCODER_package/testruns/", 20000, "hg19")
