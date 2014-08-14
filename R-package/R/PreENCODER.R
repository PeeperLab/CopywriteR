### Generate annotation files (GC-content, mappability and bed files) for desired binSize

preENCODER<-function(MAPA_GC_location, outputFolder, binSize, reference){
	
	
	####### Checks
	## Check whether folders exist.
	if(file.exists(paste(outputFolder))==FALSE){
	stop("The outputFolder could not be found. Please change your outputPath.")}
	
	if(file.exists(paste(MAPA_GC_location))==FALSE){
	stop("The annotation folder could not be found. Please change your MAPA_GC_location path.")}
		
	## Check MAPA_GC_location for 1kb files
	# For hg19 (human)
	if (reference=="hg19"){
		if(file.exists(paste(MAPA_GC_location, "/hg19_1kb/", sep=""))){
			cat("Reference folder hg19_1kb detected", "\n")
		}
	}
	# For mouse (mm10)
	if (reference=="mm10"){
		if(file.exists(paste(MAPA_GC_location, "/mm10_1kb/", sep=""))){
			cat("Reference folder mm10_1kb detected", "\n")
		}
	}
	if (!(reference == "hg19" || reference == "mm10")){
		stop("The reference is not recognised. Please provide a suitable reference (mm10 or hg19).")
	}
	
	## Check if binSize is factor of 1kb
	ifelse(!.is.wholenumber(binSize/1000), ){
		stop("Please provide a binSize which is a multiplication of 1000.")
	}
	
	####### 
	## Generate files with desired binSize (MAPA, GC, bed, blacklist)
	
	# Load blacklist, CG-content and mapability files
	load(paste(MAPA_GC_location, "/hg19_1kb/mapa_1kb.rda", sep=""))
	load(paste(MAPA_GC_location, "/hg19_1kb/GCcontent_1kb.rda", sep=""))
	bed_file<-read.table(paste(MAPA_GC_location, "/hg19_1kb/hg19-blacklist-nochr.bed", sep=""), sep="\t")
	Chr_length<-read.table(file="/Users/o.krijgsman/Desktop/Chr_length.txt", sep="\t", header=F)

	## Create output files

	###### Create bins with desired bin size
	numChr<-length(unique(GC2$Chromosome))

	Chr_length<-matrix(data=0, ncol=2, nrow=numChr)
	for (i in 1:numChr){
		Chr_length[i,1]<-unique(GC2$Chromosome)[i]
		Chr_length[i,2]<- GC2$End[max(which(GC2$Chromosome==unique(GC2$Chromosome)[i]))]
	}

	L_bin<-floor(Chr_length[,2]/binSize)
	bed_file<-matrix(data=0, ncol=3, nrow=sum(L_bin))

	for (i in 1:numChr){
		for (j in 1:L_bin[i]){
			# Chromosome 1
			if (i==1){
				bed_file[j,1]<-paste(Chr_length[i,1])
				bed_file[j,2]<-(j*binSize)-(binSize-1)
				bed_file[j,3]<-(j*binSize)
			}
			## Chromosomes 2 : Last
			if (i>1){
				bed_file[(sum(L_bin[1:(i-1)]))+j,1]<-paste(Chr_length[i,1])
				bed_file[(sum(L_bin[1:(i-1)]))+j,2]<-(j*binSize)-(binSize-1)
				bed_file[(sum(L_bin[1:(i-1)]))+j,3]<-(j*binSize)
			}
		}
	}
	cat("Generated", binSize, "bp bins for all", numChr,"chromosomes", "\n")
	
	
	###### Create mapabillity file with desired bin size



	
	cat("Generated mapabillity file for binSize", binSize,"of bp", "\n")	
	
	
	
	
	###### Create GC-content file with desired bin size




	cat("Generated GC-content file for binSize", binSize,"of bp", "\n")		
	
	
	
	# Create folder for output files
	file_name<-paste(reference, "_", binSize/1000, "kb", sep="")
	dir.create(paste(outputFolder, file_name, sep=""))
	if(file.exists(paste(outputFolder, file_name, sep=""))==FALSE){
		stop("No output folder created, please check `outputFolder` and folder permissions.")
	}

	## Write files to folder
	# Blacklist
	write.table(bed_file, file=paste(outputFolder, file_name,"/blacklist.bed", sep="\t"), quote=F, row.names=F)
	# GC-content
	write.table(?????, file=paste(outputFolder, file_name,"/GC_content.txt", sep="\t"), quote=F, row.names=F)
	# Mapability
	write.table(?????, file=paste(outputFolder, file_name,"/mapability.bed", sep="\t"), quote=F, row.names=F)	
	# bed file with bins
	write.table(bed_file, file=paste(outputFolder, file_name,"/bins.bed", sep="\t"), quote=F, row.names=F, col.names=F)

	
}

