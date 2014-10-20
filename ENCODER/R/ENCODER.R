ENCODER <- function(bam.folder, destination.folder, reference.folder,
                    which.control, ncpu, capture.regions.file) {
  
  start.time <- Sys.time()
  
  ## Check for presence capture.regions.file and make path absolute if present
  if (missing(capture.regions.file)) {
    capture.regions.file <- "not specified"
  } else {
    capture.regions.file <- tools::file_path_as_absolute(capture.regions.file)  
  }
  
  # ## Check for presence keep.intermediairy.files and set to FALSE by default
  # if (missing(keep.intermediairy.files)) {
  #   keep.intermediairy.files <- FALSE
  # }

  ## Make folder paths absolute
  bam.folder <- tools::file_path_as_absolute(bam.folder)
  destination.folder <- tools::file_path_as_absolute(destination.folder)
  reference.folder <- tools::file_path_as_absolute(reference.folder)
  
  ## Make folder path independent of trailing /
  bam.folder <- paste0(gsub("/$", "", bam.folder), "/")
  destination.folder <- paste0(gsub("/$", "", destination.folder), "/")
  reference.folder <- paste0(gsub("/$", "", reference.folder), "/")
    
  ##############################################################
  ## Generate inputStructure to run ENCODER and check folders ##
  ##############################################################
  
  ## Check the existence of folders and files
  if (!file.exists(bam.folder)) {
    stop(paste("The bam.folder could not be found. Please change your",
               "bam.folder path."))
  }
  
  if (!file.exists(destination.folder)) {
    stop(paste("The destination.folder could not be found. Please change your",
               "destination.folder path."))
  }
  
  if (!file.exists(reference.folder)) {
    stop(paste("The reference.folder could not be found. Please change your",
               "reference.folder path or run `preENCODER` to generate the",
               "required folder with GC-content and mapability files for your",
               "desired bin size."))
  }
  
  if (!file.exists(capture.regions.file) &
      capture.regions.file != "not specified") {
    stop(paste("The capture.regions.file could not be found. Please change",
               "your capture.regions.file path."))
  }

  ## Create paths to helper files
  bin.file <- paste0(reference.folder, "bins.bed")
  blacklist.file <- paste0(reference.folder, "blacklist.bed")
  gc.content.file <- paste0(reference.folder, "GC_content.bed")
  mappability.file <- paste0(reference.folder, "mapability.bed")

  ## Retrieve number of chromosomes and bin size from bin.file helper file
  bin.bed <- read.table(file = bin.file, as.is = TRUE, sep = "\t")
  colnames(bin.bed) <- c("chr", "start", "end")
  chrom <- unique(bin.bed$chr)
  nchrom <- length(chrom)
  bin.size <- bin.bed$end[1]
  
  ## List all input files and write to log
  dir.create(paste0(destination.folder, "CNAprofiles/"))
  
  ## Create lists with bam-files and index of corresponding reference
  bam.list <- list.files(path = bam.folder, pattern = ".bam$")
  control.list <- which.control

  ## Provide output for log file
  sink(file = paste0(destination.folder, "CNAprofiles/log.txt"),
       type = c("output", "message"))
       
  for(i in 1:length(bam.list)) {
    cat("The control for sample", bam.list[i], "will be",
        bam.list[control.list[i]], "\n")
  }
  cat("The bin.size for this analysis is", bin.size, "\n")
  cat("The capture region file is", capture.regions.file, "\n")
  cat("This analysis will be run on", ncpu, "cpus", "\n\n\n")
  sink()
  for(i in 1:length(bam.list)) {
    cat("The control for sample", bam.list[i], "will be",
        bam.list[control.list[i]], "\n")
  }
  cat("The capture region file is", capture.regions.file, "\n")
  cat("This analysis will be run on", ncpu, "cpus", "\n\n\n")
  
  ## Test for compatibilty chromosome names
  prefixes <- vector(mode = "character")

  tryCatch({
    for(bam in bam.list) {
      header <- scanBamHeader(paste0(bam.folder, bam))
      chr.names <- names(header[[1]]$targets)[1]
      prefixes <- append(prefixes, gsub("[[:digit:]]|X|Y", "", chr.names)[1])
    }
  }, error = function(e) {
    cat("ERROR: The BAM file header of file", bam,
        "is corrupted or truncated.\n")
    cat("ERROR: Please rebuild this BAM file or exclude it from analysis.\n")
    stop("Stopping execution of the remaining part of the script...")    
  })

  if (!all(prefixes == prefixes[1])) {
    stop(paste("The bam files have different chromosome names. Please adjust",
               "the .bam files such that they contain the same chromosome",
               "notation."))
  } else {
    prefixes <- prefixes[1]
    chr.names <- unlist(strsplit(readLines(con <- file(bin.file), n = 1),
                                 "\t"))[1]
    prefixes[2] <- gsub("[[:digit:]]|X|Y", "", chr.names)
    close(con)
    chr.names <- unlist(strsplit(readLines(con <- file(blacklist.file), n = 1),
                                 "\t"))[1]
    prefixes[3] <- gsub("[[:digit:]]|X|Y", "", chr.names)
    close(con)
    chr.names <- unlist(strsplit(readLines(con <- file(gc.content.file), n = 1),
                                 "\t"))[1]
    prefixes[4] <- gsub("[[:digit:]]|X|Y", "", chr.names)
    close(con)
    chr.names <- unlist(strsplit(readLines(con <- file(mappability.file),
                                           n = 1), "\t"))[1]
    prefixes[5] <- gsub("[[:digit:]]|X|Y", "", chr.names)
    close(con)
    if (!all(prefixes == prefixes[1])) {
      stop(paste("The bam files and supporting .bed files have different",
                 "chromosome names. Please adjust the input files such that",
                 "they contain the same chromosome notation."))
    }
  }
          
  sink(file = paste0(destination.folder, "CNAprofiles/log.txt"), append = TRUE,
       type = c("output", "message"))
  options(width = 150)
  
  dir.create(paste0(destination.folder, "CNAprofiles/BamBaiMacsFiles/"))
    
  ###############################
  ## Run the ENCODER algortihm ##
  ###############################

  ## Create list of .bam files
  setwd(bam.folder)
  cat(bam.list, "\n", sep = "\n")

  ## Index .bam files
  IndexBam <- function(bam.list) {
    indexBam(bam.list)
    paste0("indexBam(", bam.list, ")")
  }
  sfInit(parallel=TRUE, cpus = ncpu)
  sfLibrary(Rsamtools)
  to.log <- sfSapply(bam.list, IndexBam)
  sfStop()
  cat(to.log, "\n", sep = "\n")
  
  ## Check whether BAMs are paired-end
  NumberPairedEndReads <- function(bam.list) {
    flag <- scanBamFlag(isPaired = TRUE)
    param <- ScanBamParam(flag = flag)
    countBam(bam.list, param = param)
  }
  sfInit(parallel=TRUE, cpus = ncpu)
  sfLibrary(Rsamtools)
  is.paired.end <- sfLapply(bam.list, NumberPairedEndReads)
  sfStop()
  is.paired.end <- Reduce(function(x, y) {merge(x, y, all = TRUE)},
                          is.paired.end)
  is.paired.end <- ifelse(is.paired.end$records > 0, TRUE, FALSE)
  for(i in 1:length(bam.list)) {
    cat("Paired-end sequencing for sample ", bam.list[i], ": ",
        is.paired.end[i], "\n", sep = "")
  }
  cat("\n\n")

  ## Remove anomalous reads and reads with Phred < 37
  i <- c(1:length(bam.list))
  ProperReads <- function(i, bam.list, destination.folder, is.paired.end) {
    if (is.paired.end[i]) {
      flag <- scanBamFlag(isProperPair = TRUE)
      param <- ScanBamParam(flag = flag, what = "mapq")
      filter <- FilterRules(list(isHighQual = function(x) x$mapq >= 37))
      filterBam(bam.list[i], paste0(destination.folder,
                                    "CNAprofiles/BamBaiMacsFiles/",
                                    gsub(".bam$", "_properreads.bam",
                                         bam.list[i])), filter = filter,
                                    indexDestination = TRUE, param = param)
      paste0("filterBam(", bam.list[i], ", ", destination.folder,
             "CNAprofiles/BamBaiMacsFiles/", gsub(".bam$", "_properreads.bam",
                                                  bam.list[i]),
             ", filter = filter, indexDestination = TRUE, param = param)")
    } else {
      param <- ScanBamParam(what = "mapq")
      filter <- FilterRules(list(isHighQual = function(x) x$mapq >= 37))
      filterBam(bam.list[i], paste0(destination.folder,
                                    "CNAprofiles/BamBaiMacsFiles/",
                                    gsub(".bam$", "_properreads.bam",
                                         bam.list[i])), filter = filter,
                                    indexDestination = TRUE, param = param)
      paste0("filterBam(", bam.list[i], ", ", destination.folder,
             "CNAprofiles/BamBaiMacsFiles/", gsub(".bam$", "_properreads.bam",
                                                  bam.list[i]),
             ", filter = filter, indexDestination = TRUE, param = param)")
    }
  }
  sfInit(parallel=TRUE, cpus = ncpu)
  sfLibrary(Rsamtools)
  to.log <- sfSapply(i, ProperReads, bam.list, destination.folder,
                     is.paired.end)
  sfStop()
  cat(to.log, "\n", sep = "\n")

  ################
  statistics <- matrix(nrow = length(bam.list), ncol = 6,
                       dimnames = list(bam.list, c("Total", "Total properreads",
                                                   "Unmapable / Mitochondrial",
                                                   "On chromosomes",
                                                   "Not in peaks", "In peaks")))
  Stats <- function(bam.list) {
    countBam(bam.list)
  }
  sfInit(parallel=TRUE, cpus = ncpu)
  sfLibrary(Rsamtools)
  res <- sfSapply(bam.list, Stats)
  sfStop()
  statistics[,1] <- unlist(data.frame(t(res))$records)
  ################
  
  ## Create new .bam list
  setwd(paste0(destination.folder, "CNAprofiles/BamBaiMacsFiles/"))
  bam.list <- gsub(".bam$", "_properreads.bam", bam.list)
  cat(bam.list, "\n", sep = "\n")

  ################
  Stats <- function(bam.list, chrom) {
    all.reads <- countBam(bam.list)$records
    which <- GRanges(seqnames = chrom,
                     ranges = IRanges(start = rep(1, length(chrom)),
                     end = rep(536870912, length(chrom)))) # 536870912 is max
    what <- c("pos")
    param <- ScanBamParam(which = which, what = what)
    chrom.reads <- countBam(file = bam.list, param = param)
    chrom.reads <- sum(chrom.reads$records)
    c(all.reads = all.reads, chrom.reads = chrom.reads)
  }
  sfInit(parallel=TRUE, cpus = ncpu)
  sfLibrary(Rsamtools)
  res <- data.frame(t(sfSapply(bam.list, Stats, chrom)))
  sfStop()
  statistics[, 2] <- res$all.reads
  statistics[, 3] <- res$all.reads - res$chrom.reads
  statistics[, 4] <- res$chrom.reads
  ################
  
  ## Create list with numbers of controls
  control.numbers <- unique(control.list) #####
  
  ## Call peaks in .bam file of control sample
  Macs14 <- function(control.numbers, bam.list) {
    system(paste0("macs14 -t " , bam.list[control.numbers], " -n MACS",
                  control.numbers, " -g hs --nolambda"))
    paste0("macs14 -t " , bam.list[control.numbers], " -n MACS",
           control.numbers, " -g hs --nolambda")
  }
  sfInit(parallel=TRUE, cpus = min(length(control.numbers), ncpu))
  to.log <- sfSapply(control.numbers , Macs14, bam.list)
  sfStop()
  cat(to.log, "\n", sep = "\n")
  
  ## Alternative for bedtools
#   outside.peaks.grange <- list()
#   all.grange <- GRanges(seqnames = chrom, ranges = IRanges(start = rep(1, length(chrom)), end = rep(536870912, length(chrom))))
#   for(control in 1:1) {
#     bed <- read.table(file = paste0("MACS", control, "_peaks.bed"), as.is = TRUE, sep = "\t")
#     colnames(bed) <- c("chr", "start", "end", "peak.name")
#     peaks.grange <- with(bed, GRanges(Rle(chr), IRanges(start, end)))
#     outside.peaks.grange[control] <- setdiff(all.grange, peaks.grange)
#   }
  # i <- 1:length(bam.list)
  # Intersect <- function(i, bam.list, control.list) {
    # bed <- read.table(file = , as.is = TRUE, sep = "\t")
  # }
  
  ## Remove peak-regions from .bam files
  i <- 1:length(bam.list)
  PeakRm <- function(i, bam.list, control.list) {
    system(paste0("bedtools intersect -abam ", bam.list[i], " -b MACS", control.list[i], "_peaks.bed -v > ", gsub(".bam$", "_peakrm.bam", bam.list[i])))
    paste0("bedtools intersect -abam ", bam.list[i], " -b MACS", control.list[i], "_peaks.bed -v > ", gsub(".bam$", "_peakrm.bam", bam.list[i]))
  }
  sfInit(parallel=TRUE, cpus = ncpu)
  to.log <- sfSapply(i, PeakRm, bam.list, control.list)
  sfStop()
  cat(to.log, "\n", sep = "\n")

  ## Remove _properreads.bam(.bai) files
  # if (!keep.intermediairy.files) {
  #   file.remove(bam.list)
  #  file.remove(gsub(".bam$", ".bam.bai", bam.list))
  #  file.remove(list.files(pattern = "MACS"))
  # }

  ## Create new .bam list
  bam.list <- gsub(".bam$", "_peakrm.bam", bam.list)
  cat(bam.list, "\n", sep = "\n")
  
  ## Index _peakrm.bam files
  IndexPeakRm <- function(bam.list) {
    indexBam(bam.list)
    paste0("indexBam(", bam.list, ")")
  }
  sfInit(parallel=TRUE, cpus = ncpu)
  sfLibrary(Rsamtools)
  to.log <- sfSapply(bam.list, IndexPeakRm)
  sfStop()
  cat(to.log, "\n", sep = "\n")

  ################
  Stats <- function(bam.list, chrom) {
    which <- GRanges(seqnames = chrom, ranges = IRanges(start = rep(1, length(chrom)), end = rep(536870912, length(chrom)))) ## 536870912 is the maximum number for 'end'
    what <- c("pos")
    param <- ScanBamParam(which = which, what = what)
    chrom.reads <- countBam(file = bam.list, param = param)
    c(chrom.reads = sum(chrom.reads$records))
  }
  sfInit(parallel=TRUE, cpus = ncpu)
  sfLibrary(Rsamtools)
  res <- data.frame(t(data.frame(sfSapply(bam.list, Stats, chrom, simplify = FALSE))))
  sfStop()
  statistics[, 5] <- res$chrom.reads
  statistics[, 6] <- statistics[, 4] - res$chrom.reads
  ################
  
  print(statistics)
  cat("\n\n")
  
  ## Create read.count matrix
  bed <- read.table(file = bin.file, sep = "\t")
  read.count <- matrix(nrow = nrow(bed))
  read.count[,1] <- paste(bed[,1], paste(bed[,2], bed[,3], sep = "-"), sep = ":")
  read.count <- cbind(read.count, bed[, 1:3])

  ## Calculate the number of reads per bin
  CalculateDepthOfCoverage <- function(bam.list, bed) {
    library(Rsamtools)
    which <- RangedData(space = bed[,1], IRanges(bed[,2], bed[,3]))
    param <- ScanBamParam(which = which, what = c("pos"))
    bamreads <- scanBam(file = paste(bam.list), param = param)
    readmap <- matrix(data = 0, ncol = 1, nrow = nrow(bed))
    for (j in 1:length(bamreads)) {
      readmap[j] <- length(unlist(bamreads[j]))
    }
    return(list(readmap, paste0("Rsamtools finished calculating reads per bin in sample ", bam.list, "; number of bins = ", length(bamreads))))
  }
  sfInit(parallel=TRUE, cpus = ncpu)
  res <- sfSapply(bam.list, CalculateDepthOfCoverage, bed)
  sfStop()
  for(i in seq(1,2*length(bam.list),2)) {
    read.count <- cbind(read.count, res[[i]])
  }
  cat(unlist(res[2,]), "\n", sep = "\n")
  sink()
  
  ## Compensate for removal of reads in peak regions
  read.count <- cbind(read.count, read.count[,-c(1, 2, 3, 4)], read.count[,-c(1, 2, 3, 4)])
  
  for(control.number in control.numbers) {

    # Calculate overlap of peaks with bins using bedtools
    intersection <- system(paste0("bedtools intersect -a ", bin.file, " -b ", destination.folder,
      "CNAprofiles/BamBaiMacsFiles/MACS", control.number, "_peaks.bed -wao"), intern = TRUE)
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
    control.for <- grep(control.number, which.control)
    if (all(read.count[,2] == intersection[,1] & read.count[,3] == intersection[,2] & read.count[,4] == intersection[,3])) {
      fraction_of_bin <- (bin.size-as.numeric(intersection[,4])) / bin.size
      read.count[, (4 + 2 * length(bam.list) + control.for)] <- fraction_of_bin
      for(control in control.for) {
        read.count[, (4 + control)] <- ifelse(fraction_of_bin != 0
          , as.numeric(read.count[, (4 + length(bam.list) + control)]) / fraction_of_bin
          , 0)
      }
    }
    else {
      stop("Bins do not correspond")
    }
  }
  
  colnames(read.count) <- c("BinID", "Chromosome", "StartPos", "StopPos", sub(pattern = "$", ".compensated", paste(bam.list)),
    paste(bam.list), sub(pattern = "$", ".fractionOfBin", paste(bam.list)))

  ## Create output file
  output <- read.count
  output[, 2] <- gsub(prefixes[1], "", output[, 2])
  output[,2] <- gsub("X", nchrom - 1, output[,2])
  output[,2] <- gsub("Y", nchrom, output[,2])
  write.table(output, file = paste0(destination.folder, "CNAprofiles/", "read.counts_compensated.txt"), row.names = FALSE, col.names = TRUE, sep = "\t")

  ## Create histograms of fraction_of_bin (fraction of length in bins (after peak region removal)
  dir.create(paste0(destination.folder, "CNAprofiles/qc/"))
  for(i in 1:length(bam.list)) {
    pdf(paste0(destination.folder, "CNAprofiles/qc/fraction_of_bin_", i, ".pdf"), width=7, height=7)
    plot(ecdf(as.numeric(read.count[,4+(2*length(bam.list))+i])), verticals = TRUE, ylab = "Fraction of bins",
      xlab = "Remaining fraction of bin after peak removal", main = "Cumulative distribution of remaining bin fraction")
    dev.off()
  }

  #############################################
  ## Normalize for GC-content and mapability ##
  #############################################

  ## Read files
  read.count <- read.count[,1:(4 + length(bam.list))]
  for(i in 5:ncol(read.count)) {
    write.table(read.count[,c(2,3,4,i)], colnames(read.count)[i], quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  }
  f <- list.files(pattern = "bam.compensated")
  
  gc <- read.delim(gc.content.file)[,c(1,2,3,5)]
  colnames(gc) <- c("chr","start","end","gc")
  mapa <- read.delim(mappability.file, header=FALSE, col.names=c("chr","start","end","mapa"))
  black <- read.delim(blacklist.file, header=F, col.names=c("chr", "start", "end"))
  
  data <- .loadCovData(f, gc=gc, mapa=mapa, black=black, excludechr="MT", datacol=4)
  sampnames <- paste(sub("$","", f), paste0(round(colSums(data$cov) / 1e6,1),"M"), sep="_")
  colnames(data$cov) <- sampnames
  usepoints <- !(data$anno$chr %in% c("X","Y","MT", "chrX", "chrY", "chrM"))
  
  ## Perform normalization (in .tng helper function)
  tryCatch({
    i <- c(1:ncol(data$cov))
    NormalizeDOC <- function(i, data, .tng, usepoints, destination.folder) {  
      .tng(data.frame(count = data$cov[,i], gc = data$anno$gc, mapa = data$anno$mapa), use = usepoints, correctmapa = TRUE, plot = paste0(destination.folder, "CNAprofiles/qc/", colnames(data$cov)[i], ".png"))
    }
    sfInit(parallel=TRUE, cpus = ncpu)
    ratios <- sfLapply(i, NormalizeDOC, data, .tng, usepoints, destination.folder)
    sfStop()
    ratios <- matrix(unlist(ratios), ncol = length(bam.list))
  }, error = function(e) {
    cat("ERROR: The GC-content and mappability normalization did not work due to a failure to calculate loesses.\n")
    cat("ERROR: This can generally be solved by using larger bin sizes.\n")
    stop("Stopping execution of the remaining part of the script...")    
  })
  
  colnames(ratios) <- sampnames
  
  rd <- list(ratios=ratios[!data$anno$black & !is.na(data$anno$mapa) & data$anno$mapa > .2, ], anno=data$anno[!data$anno$black & !is.na(data$anno$mapa) & data$anno$mapa > .2, ])
  
  sink(file = paste0(destination.folder, "CNAprofiles/log.txt"), append = TRUE, type = c("output", "message"))
  rd2 <- rd
  for(i in 1:length(bam.list)) {
    rd2$ratios[,i] <- rd$ratios[,i] - rd$ratios[,which.control[i]]
    cat("Relative log2-values are calculated for sample", colnames(rd2$ratios)[i], "with control", colnames(rd2$ratios)[which.control[i]], "\n")
  }
  colnames(rd2$ratios) <- gsub("$", ".rel", colnames(rd2$ratios))
  cat("\n\n")

  ###################
  ## Create output ##
  ###################
  
  ## Create table with corrected log2 values and write to file
  read.count <- matrix(nrow = nrow(rd$ratios))
  read.count <- cbind(read.count, rd$anno[1:3])
  read.count[,1] <- paste(rd$anno[,1], paste(rd$anno[,2], rd$anno[,3], sep = "-"), sep = ":")
  read.count <- cbind(read.count, rd$ratios, rd2$ratios)
  colnames(read.count) <- c("BinID", "Chromosome", "StartPos", "StopPos", colnames(read.count[5:length(colnames(read.count))]))
  
  read.count <- read.count[-which(rowSums(is.na(read.count[,-c(1:4)])) > 0),]
  read.count[, 2] <- gsub(prefixes[1], "", read.count[, 2])
  read.count[, 2] <- gsub("X", nchrom - 1, read.count[, 2])
  read.count[, 2] <- gsub("Y", nchrom, read.count[, 2])
  read.count[, 2] <- as.integer(read.count[, 2])
  read.count[read.count == -Inf] <- -.Machine$integer.max/2
  read.count[read.count == Inf] <- .Machine$integer.max/2
  
  write.table(read.count, paste0(destination.folder, "CNAprofiles/log2ratio_compensated_corrected.txt"), sep="\t", row.names = FALSE, quote = FALSE)

  ##############################################################################
  ## Calculate overlap with capture.regions.file for quality control purpose ##
  ##############################################################################
  
  if (capture.regions.file != "not specified") {
    for(control.number in control.numbers) {
      intersection <- system(paste0("bedtools intersect -a ", capture.regions.file, " -b ", destination.folder, "CNAprofiles/BamBaiMacsFiles/MACS", control.number, "_peaks.bed -wao"), intern = TRUE)
      intersection <- strsplit(intersection, "\t")
      intersection <- do.call(rbind, intersection)
      intersection <- intersection[,c(1,2,3,9)]
      vec <- vector()
      for (i in 2:nrow(intersection)) {
        if (intersection[i,1] == intersection[i-1,1] && intersection[i,2] == intersection[i-1,2] && intersection[i,3] == intersection[i-1,3] ) {
          intersection[i,4] <- as.integer(intersection[i,4]) + as.integer(intersection[i-1,4])
          vec <- append(vec, i-1)
        }
      }
      intersection <- intersection[-vec,]
      
      cat("Number of exons covered by peaks in sample ", bam.list[control.number], ": ", sum(intersection[,4] != "0"), "\n")
    
      intersection <- system(paste0("bedtools intersect -a ", destination.folder, "CNAprofiles/BamBaiMacsFiles/MACS", control.number, "_peaks.bed -b ", capture.regions.file, " -wao"), intern = TRUE)
      intersection <- strsplit(intersection, "\t")
      intersection <- do.call(rbind, intersection)
      intersection <- intersection[,c(1,2,3,9)]
      vec <- vector()
      for (i in 2:nrow(intersection)) {
        if (intersection[i,1] == intersection[i-1,1] && intersection[i,2] == intersection[i-1,2] && intersection[i,3] == intersection[i-1,3] ) {
          intersection[i,4] <- as.integer(intersection[i,4]) + as.integer(intersection[i-1,4])
          vec <- append(vec, i-1)
        }
      }
      intersection <- intersection[-vec,]
      
      cat("Number of peaks covered by exons in sample", bam.list[control.number], ": ", sum(intersection[,4] != "0"), "\n")
      cat("Total number of exons in sample", bam.list[control.number], ": ", nrow(read.table(capture.regions.file, sep = "\t")), "\n")
      cat("Total number of peaks", bam.list[control.number], ":", nrow(read.table(paste0(destination.folder, "CNAprofiles/BamBaiMacsFiles/MACS", control.number, "_peaks.bed"), sep = "\t")), "\n\n")  
    }
  }
  
  sink()
  ## Remove BamBaiMacsFiles folder
  # if (!keep.intermediairy.files) {
  #   unlink(paste0(destination.folder, "CNAprofiles/BamBaiMacsFiles/"))
  # }
  cat("Total calculation time: ", Sys.time() - start.time, "\n\n")
  
  inputStructure <- list(destination.folder = destination.folder, ncpu = ncpu, nchrom = nchrom)
  save(inputStructure, file = paste0(destination.folder, "CNAprofiles/input.Rdata"))

}
