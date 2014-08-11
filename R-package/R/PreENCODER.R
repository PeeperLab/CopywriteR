### Generate annotation files (GC-content, mappability and bed files) for desired binSize
##

## To run
preENCODER("/Users/o.krijgsman/Documents/PostDoc/Projects/ENCODER_package/Annotation",	"/Users/o.krijgsman/Documents/PostDoc/Projects/ENCODER_package/testruns/", 20000, "hg19")




preENCODER<-function(MAPA_GC_location, outputFolder, binSize, reference){
	
	## Check MAPA_GC_location for 1kb files
	if (reference=="hg19"){
		
	}
	
	if (reference=="mm10"){
	
	}
	
	
	## Check binSize
	ifelse(is.interger(binSize/2), )
	
	## Generate files with desired binSize (MAPA, GC, bed, blacklist)
	
}