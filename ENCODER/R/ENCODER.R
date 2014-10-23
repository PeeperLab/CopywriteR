ENCODER <- function(sample.control, destination.folder,
                    reference.folder, ncpu, capture.regions.file) {
  
  ##########################
  ## Check and initialize ##
  ##########################

  start.time <- Sys.time()
  
  ## Make capture.regions.file path absolute if present
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
  sample.control <- apply(sample.control, c(1, 2), function(x) {
    if (!is.na(x)) {
      tools::file_path_as_absolute(x)
    } else {
      x <- NA	
    }
  })
  sample.control <- data.frame(sample.control, stringsAsFactors = FALSE)
  colnames(sample.control) <- c("samples", "controls")
  destination.folder <- tools::file_path_as_absolute(destination.folder)
  reference.folder <- tools::file_path_as_absolute(reference.folder)
  
  ## Make folder path independent of trailing /
  destination.folder <- paste0(gsub("/$", "", destination.folder), "/")
  reference.folder <- paste0(gsub("/$", "", reference.folder), "/")

  ## Check the existence of folders and files
  invisible(apply(sample.control, c(1, 2), function(x) {
    if (!is.na(x) & !file.exists(x)) {
      stop("The file ", x, " could not be found. Please change the path to ",
           "this file.")
    }
  }))
  
  if (!file.exists(destination.folder)) {
    stop("The destination.folder could not be found. Please change your ",
         "destination.folder path.")
  }
  
  if (!file.exists(reference.folder)) {
    stop("The reference.folder could not be found. Please change your ",
         "reference.folder path or run `preENCODER` to generate the required ",
         "folder with GC-content and mapability files for your desired bin ",
         "size.")
  }
  
  if (!file.exists(capture.regions.file) &
      capture.regions.file != "not specified") {
    stop("The capture.regions.file could not be found. Please change your ",
         "capture.regions.file path.")
  }

  ## Create lists with BAM files and index of corresponding control
  sample.paths <- unlist(sample.control)
  sample.paths <- unique(sample.paths[!is.na(sample.paths)])
  sample.files <- unname(sapply(sample.paths, function(x) {
    x <- unlist(strsplit(x, "/"))[length(unlist(strsplit(x, "/")))]
  }))
  single.index <- sample.control$samples[is.na(sample.control$controls)]
  single.index <- match(single.index, sample.paths)
  dual.index <- sample.control[!is.na(sample.control$controls), ]
  dual.index[, ] <- apply(dual.index, c(1, 2), function(x) {
    match(x, sample.paths)
  })

  ## Create paths to helper files
  bin.file <- paste0(reference.folder, "bins.bed")
  blacklist.file <- paste0(reference.folder, "blacklist.bed")
  gc.content.file <- paste0(reference.folder, "GC_content.bed")
  mapability.file <- paste0(reference.folder, "mapability.bed")

  ## Retrieve number of chromosomes and bin size from bin.bed helper file
  bin.bed <- read.table(file = bin.file, as.is = TRUE, sep = "\t")
  colnames(bin.bed) <- c("Chromosome", "Start", "End")
  chrom <- unique(bin.bed$Chromosome)
  nchrom <- length(chrom)
  bin.size <- bin.bed$End[1]
  
  ## Create folders
  destination.folder <- paste0(destination.folder, "CNAprofiles/")
  tryCatch({
    if (!file.exists(paste0(destination.folder, "BamBaiMacsFiles/"))) {
      dir.create(paste0(destination.folder, "BamBaiMacsFiles/"),
                 recursive = TRUE)
    }
  }, warning = function(e) {
    cat("ERROR: You do not have write permissions in the destination folder.\n")
    stop("Stopping execution of the remaining part of the script...")
  })
  
  ## Provide output for log file
  sink(file = paste0(destination.folder, "log.txt"),
       type = c("output", "message"))
  options(width = 150)

  cat("The following samples will be analyzed single-channel:\n")
  cat(sample.files[single.index], sep = "\n")
  cat("\nThe following samples will be analyzed dual-channel:\n")
  cat(paste0("sample: ", sample.files[dual.index$samples],
             ";\tmatching control: ", sample.files[dual.index$controls]),
      sep = "\n")
  cat("The bin size for this analysis is", bin.size, "\n")
  cat("The capture region file is", capture.regions.file, "\n")
  cat("This analysis will be run on", ncpu, "cpus", "\n\n\n")
  sink()
  
  cat("The following samples will be analyzed single-channel:\n")
  cat(sample.files[single.index], sep = "\n")
  cat("\nThe following samples will be analyzed dual-channel:\n")
  cat(paste0("sample: ", sample.files[dual.index$samples],
             ";\tmatching control: ", sample.files[dual.index$controls]),
      sep = "\n")
  cat("The bin size for this analysis is", bin.size, "\n")
  cat("The capture region file is", capture.regions.file, "\n")
  cat("This analysis will be run on", ncpu, "cpus", "\n\n\n")
  
  ## Test for compatibilty chromosome names
  prefixes <- vector(mode = "character")

  tryCatch({
    for(samp in sample.paths) {
      header <- scanBamHeader(samp)
      chr.names <- names(header[[1]]$targets)[1]
      prefixes <- append(prefixes, gsub("[[:digit:]]|X|Y", "", chr.names)[1])
    }
  }, error = function(e) {
    cat("ERROR: The BAM file header of file", samp, "is corrupted or",
        "truncated.\n")
    cat("ERROR: Please rebuild this BAM file or exclude it from analysis.\n")
    stop("Stopping execution of the remaining part of the script...")    
  })

  if (!all(prefixes == prefixes[1])) {
    stop("The bam files have different chromosome names. Please adjust the ",
         ".bam files such that they contain the same chromosome notation.")
  } else {
    prefixes <- prefixes[1]
    
    con <- file(bin.file)
    chr.names <- unlist(strsplit(readLines(con, n = 1), "\t"))[1]
    prefixes[2] <- gsub("[[:digit:]]|X|Y", "", chr.names)
    close(con)
    
    con <- file(blacklist.file)
    chr.names <- unlist(strsplit(readLines(con, n = 1), "\t"))[1]
    prefixes[3] <- gsub("[[:digit:]]|X|Y", "", chr.names)
    close(con)
    
    con <- file(gc.content.file)
    chr.names <- unlist(strsplit(readLines(con, n = 1), "\t"))[1]
    prefixes[4] <- gsub("[[:digit:]]|X|Y", "", chr.names)
    close(con)
    
    con <- file(mapability.file)
    chr.names <- unlist(strsplit(readLines(con, n = 1), "\t"))[1]
    prefixes[5] <- gsub("[[:digit:]]|X|Y", "", chr.names)
    close(con)
    
    if (!all(prefixes == prefixes[1])) {
      stop("The bam files and supporting .bed files have different chromosome ",
           "names. Please adjust the input files such that they contain the ",
           "same chromosome notation.")
    }
  }  ## UP TO HERE

  ########################################################
  ## Calculate depth of coverage using off-target reads ##
  ########################################################

  sink(file = paste0(destination.folder, "log.txt"), append = TRUE,
       type = c("output", "message"))
    
  ## Create list of .bam files
  cat(sample.files, "\n", sep = "\n")

  ## Index .bam files
  IndexBam <- function(sample.paths) {
    indexBam(sample.paths)
    paste0("indexBam(\"", sample.paths, "\")")
  }
  sfInit(parallel=TRUE, cpus = ncpu)
  sfLibrary(Rsamtools)
  to.log <- sfSapply(sample.paths, IndexBam)
  sfStop()
  cat(to.log, "\n", sep = "\n")
  
  ## Check whether BAMs are paired-end
  NumberPairedEndReads <- function(sample.paths) {
    flag <- scanBamFlag(isPaired = TRUE)
    param <- ScanBamParam(flag = flag)
    countBam(sample.paths, param = param)
  }
  sfInit(parallel=TRUE, cpus = ncpu)
  sfLibrary(Rsamtools)
  is.paired.end <- sfLapply(sample.paths, NumberPairedEndReads)
  sfStop()
  is.paired.end <- Reduce(function(x, y) {
                            merge(x, y, all = TRUE)
                          }, is.paired.end)
  is.paired.end <- ifelse(is.paired.end$records > 0, TRUE, FALSE)
  for(i in 1:length(sample.files)) {
    cat("Paired-end sequencing for sample ", sample.files[i], ": ",
        is.paired.end[i], "\n", sep = "")
  }
  cat("\n\n")

  ## Remove anomalous reads and reads with Phred < 37
  i <- c(1:length(sample.paths))
  ProperReads <- function(i, sample.paths, destination.folder, sample.files,
                          is.paired.end) {
    if (is.paired.end[i]) {
      flag <- scanBamFlag(isProperPair = TRUE)
      param <- ScanBamParam(flag = flag, what = "mapq")
      filter <- FilterRules(list(isHighQual = function(x) {
                                   x$mapq >= 37
                                 }))
      filterBam(sample.paths[i],
                paste0(destination.folder, "BamBaiMacsFiles/",
                       gsub(".bam$", "_properreads.bam", sample.files[i])),
                filter = filter, indexDestination = TRUE, param = param)
      paste0("filterBam(\"", sample.paths[i], "\", \"", destination.folder,
             "BamBaiMacsFiles/", gsub(".bam$", "_properreads.bam",
                                      sample.files[i]), "\", filter = filter, ",
             "indexDestination = TRUE, param = param)")
    } else {
      param <- ScanBamParam(what = "mapq")
      filter <- FilterRules(list(isHighQual = function(x) {
                                   x$mapq >= 37
                                 }))
      filterBam(sample.paths[i],
                paste0(destination.folder, "BamBaiMacsFiles/",
                       gsub(".bam$", "_properreads.bam", sample.files[i])),
                filter = filter, indexDestination = TRUE, param = param)
      paste0("filterBam(\"", sample.files[i], "\", \"", destination.folder,
             "BamBaiMacsFiles/", gsub(".bam$", "_properreads.bam",
                                      sample.files[i]), "\", filter = filter, ",
             "indexDestination = TRUE, param = param)")
    }
  }
  sfInit(parallel=TRUE, cpus = ncpu)
  sfLibrary(Rsamtools)
  to.log <- sfSapply(i, ProperReads, sample.paths, destination.folder,
                     sample.files, is.paired.end)
  sfStop()
  cat(to.log, "\n", sep = "\n")

  ## Read count statistics
  Stats <- function(sample.paths) {
    countBam(sample.paths)
  }
  sfInit(parallel=TRUE, cpus = ncpu)
  sfLibrary(Rsamtools)
  res <- sfSapply(sample.paths, Stats)
  sfStop()
  statistics <- data.frame(t(res))[, "records", drop = FALSE]
  rownames(statistics) <- sample.files
  statistics <- within(statistics, {
    total <- records
    rm(records)
  })
  
  ## Create new .bam list
  setwd(paste0(destination.folder, "BamBaiMacsFiles/"))
  sample.files <- gsub(".bam$", "_properreads.bam", sample.files)
  cat(sample.files, "\n", sep = "\n")

  ## Read count statistics
  Stats <- function(sample.files) {
    countBam(sample.files)$records
  }
  sfInit(parallel=TRUE, cpus = ncpu)
  sfLibrary(Rsamtools)
  res <- sfSapply(sample.files, Stats)
  sfStop()
  statistics <- within(statistics, {
    total.properreads <- res
  })
  
  ## Create list with numbers of controls
  control.numbers <- unique(dual.index$controls)
  
  ## Call peaks in .bam file of control sample
  Macs14 <- function(control.numbers, sample.files) {
    system(paste0("macs14 -t " , sample.files[control.numbers], " -n MACS",
                  control.numbers, " -g hs --nolambda"))
    paste0("macs14 -t " , sample.files[control.numbers], " -n MACS",
           control.numbers, " -g hs --nolambda")
  }
  sfInit(parallel=TRUE, cpus = min(length(control.numbers), ncpu))
  to.log <- sfSapply(control.numbers , Macs14, sample.files)
  sfStop()
  cat(to.log, "\n", sep = "\n")
  
  ## Alternative for bedtools
  # Create GRanges file for bins
  bin.grange <- with(bin.bed,
                     GRanges(seqnames = Chromosome,
                             ranges = IRanges(start = Start, end = End)))
  
  i <- c(1:length(sample.files))
  CalculateDepthOfCoverage <- function(i, sample.files, control.list,
                                       bin.grange, bin.size) {
    # Create GRanges object of peak files
    bed <- read.table(file = paste0("MACS", control.list[i], "_peaks.bed"),
                      as.is = TRUE, sep = "\t")
    colnames(bed) <- c("Chromosome", "Start", "End")
    peak.grange <- with(bed, GRanges(seqnames = Chromosome,
                        ranges = IRanges(start = Start, end = End)))
    
    # Calculate setdiff without reducing ranges
    outside.peak.grange <- split(bin.grange,
                                 rep_len(c(1, 2),
                                         length.out = length(bin.grange)))
    outside.peak.grange <- lapply(outside.peak.grange, function(x) {
      setdiff(x, peak.grange)
    })
    outside.peak.grange <- c(outside.peak.grange[[1]], outside.peak.grange[[2]])
    outside.peak.grange <- outside.peak.grange[order(outside.peak.grange)]
    
    # Fix the 1-based coordinate system
    ranges(outside.peak.grange)@width[1] <-
      (ranges(outside.peak.grange)@width[1] + 1L)
    ranges(outside.peak.grange)@start[1] <- 0L
    
    # countBam on remainder of bins
    param <- ScanBamParam(which = outside.peak.grange, what = c("pos"))
    counts <- countBam(sample.files[i], param = param)
    
    # Fix MT as levels in factor seqnames (remainder from setdiff operation)
    counts$space <- as.factor(as.character(counts$space))
    
    # Aggregate counts in bins
    counts <- data.table(counts)
    counts <- within(counts, {
      aggregate.factor <- ceiling(end / bin.size)
      range.length <- end - start
    })
    counts <- counts[, list(start = min(start),
                            end = max(end),
                            records = sum(records),
                            range.length = sum(range.length)),
                     by = c("space", "aggregate.factor")]
    counts <- data.frame(counts)
    
    # Replace bins by real bins & calculate compensated read counts
    if (all.equal(counts$space, as.factor(seqnames(bin.grange)))) {
      counts <- within(counts, {
        Chromosome <- as.factor(seqnames(bin.grange))
        Start <- ranges(bin.grange)@start
        End <- ranges(bin.grange)@start + ranges(bin.grange)@width - 1L
        Feature <- paste0(Chromosome, ":", paste0(Start, "-", End))
        assign(paste0("read.counts.", sample.files[i]), records)
        assign(paste0("read.counts.compensated.", sample.files[i]),
               records / (range.length / bin.size))
        assign(paste0("fraction.of.bin.", sample.files[i]),
               range.length / bin.size)
        rm(space, start, end, records, aggregate.factor, range.length)
      })
    } else {
      stop("Chromosome names do not match between the counts and bin.grange ",
           "variables. Stopping execution of the remaining part of the ",
           "script...")
    }

    # Return
    return(list(counts[, paste0("read.counts.", sample.files[i]), drop = FALSE],
                counts[, paste0("read.counts.compensated.", sample.files[i]),
                       drop = FALSE],
                counts[, paste0("fraction.of.bin.", sample.files[i]),
                       drop = FALSE],
                paste0("Rsamtools finished calculating reads per bin in ",
                       "sample ", sample.files[i], "; number of bins = ",
                       nrow(counts))))
  }
  sfInit(parallel=TRUE, cpus = ncpu)
  sfLibrary(Rsamtools)
  sfLibrary(data.table)
  res <- sfSapply(i, CalculateDepthOfCoverage, sample.files,
                  dual.index$controls, bin.grange, bin.size)
  sfStop()
  read.counts <- data.frame(Chromosome = as.factor(seqnames(bin.grange)), 
                            Start = ranges(bin.grange)@start,
                            End = ranges(bin.grange)@start +
                              ranges(bin.grange)@width - 1L,
                            Feature = paste0(Chromosome, ":",
                                              paste0(Start, "-", End)))
  read.counts <- cbind(read.counts[, ], Reduce(cbind, res[1, ]),
                           Reduce(cbind, res[2, ]), Reduce(cbind, res[3, ]))
  cat(unlist(res[4, ]), "\n", sep = "\n")
  sink()
   
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ## Remove peak-regions from .bam files
  i <- 1:length(sample.files)
  PeakRm <- function(i, sample.files, control.list) {
    system(paste0("bedtools intersect -abam ", sample.files[i], " -b MACS",
                  control.list[i], "_peaks.bed -v > ",
                  gsub(".bam$", "_peakrm.bam", sample.files[i])))
    paste0("bedtools intersect -abam ", sample.files[i], " -b MACS",
                  control.list[i], "_peaks.bed -v > ",
                  gsub(".bam$", "_peakrm.bam", sample.files[i]))
  }
  sfInit(parallel=TRUE, cpus = ncpu)
  to.log <- sfSapply(i, PeakRm, sample.files, dual.index$controls)
  sfStop()
  cat(to.log, "\n", sep = "\n")

  ## Remove _properreads.bam(.bai) files
  # if (!keep.intermediairy.files) {
  #   file.remove(sample.files)
  #  file.remove(gsub(".bam$", ".bam.bai", sample.files))
  #  file.remove(list.files(pattern = "MACS"))
  # }

  ## Create new .bam list
  sample.files <- gsub(".bam$", "_peakrm.bam", sample.files)
  cat(sample.files, "\n", sep = "\n")
  
  ## Index _peakrm.bam files
  IndexPeakRm <- function(sample.files) {
    indexBam(sample.files)
    paste0("indexBam(", sample.files, ")")
  }
  sfInit(parallel=TRUE, cpus = ncpu)
  sfLibrary(Rsamtools)
  to.log <- sfSapply(sample.files, IndexPeakRm)
  sfStop()
  cat(to.log, "\n", sep = "\n")

  ## Read count statistics
  Stats <- function(sample.files) {
    countBam(sample.files)$records
  }
  sfInit(parallel=TRUE, cpus = ncpu)
  sfLibrary(Rsamtools)
  res <- sfSapply(sample.files, Stats)
  sfStop()
  statistics <- within(statistics, {
    not.in.peaks <- res
    in.peaks <- total.properreads - not.in.peaks
  })
  print(statistics)
  cat("\n\n")
  
  ## Create read.counts matrix
  read.counts <- matrix(nrow = nrow(bin.bed))
  read.counts[,1] <- paste(bin.bed[,1],
                          paste(bin.bed[,2], bin.bed[,3], sep = "-"),
                          sep = ":")
  read.counts <- cbind(read.counts, bin.bed[, 1:3])

  ## Calculate the number of reads per bin
  CalculateDepthOfCoverage <- function(sample.files, bin.bed) {
    library(Rsamtools)
    which <- RangedData(space = bin.bed[,1], IRanges(bin.bed[,2], bin.bed[,3]))
    param <- ScanBamParam(which = which, what = c("pos"))
    bamreads <- scanBam(file = paste(sample.files), param = param)
    readmap <- matrix(data = 0, ncol = 1, nrow = nrow(bin.bed))
    for (j in 1:length(bamreads)) {
      readmap[j] <- length(unlist(bamreads[j]))
    }
    return(list(readmap, paste0("Rsamtools finished calculating reads per bin ",
                                "in sample ", sample.files, "; number of bins = ",
                                length(bamreads))))
  }
  sfInit(parallel=TRUE, cpus = ncpu)
  res <- sfSapply(sample.files, CalculateDepthOfCoverage, bin.bed)
  sfStop()
  for(i in seq(1, 2 * length(sample.files), 2)) {
    read.counts <- cbind(read.counts, res[[i]])
  }
  cat(unlist(res[2, ]), "\n", sep = "\n")
  sink()
  
  ## Compensate for removal of reads in peak regions
  read.counts <- cbind(read.counts, read.counts[, -c(1, 2, 3, 4)],
                      read.counts[, -c(1, 2, 3, 4)])
  
  for(control.number in control.numbers) {

    # Calculate overlap of peaks with bins using bedtools
    intersection <- system(paste0("bedtools intersect -a ", bin.file, " -b ",
                                  destination.folder,
                                  "BamBaiMacsFiles/MACS",
                                  control.number, "_peaks.bed -wao"),
                           intern = TRUE)
    intersection <- strsplit(intersection, "\t")
    intersection <- as.data.frame(do.call(rbind, intersection))
    intersection <- intersection[, c(1, 2, 3, 9)]
    colnames(intersection) <- c("Chromosome", "StartPos", "StopPos", "Overlap")
    intersection$Overlap <- as.numeric(as.character(intersection$Overlap))
    
    ## Calculate the cumulative overlap using data.table
    intersection <- data.table(intersection)
    intersection <- intersection[, lapply(.SD, sum),
                                 by = c("Chromosome", "StartPos", "StopPos")]
    intersection <- as.data.frame(intersection)

    # Check correspondence bins, calculate fraction of remaining bin length
    # (after peak region removal), and calculate compensated read numbers
    control.for <- grep(control.number, which.control)
    if (all(read.counts[, 2] == intersection[, 1] &
            read.counts[, 3] == intersection[, 2] &
            read.counts[, 4] == intersection[, 3])) { ## Make into one check?
      fraction.of.bin <- (bin.size - as.numeric(intersection[, 4])) / bin.size
      read.counts[, (4 + 2 * length(sample.files) + control.for)] <- fraction.of.bin
      for(control in control.for) {
        read.counts[, (4 + control)] <- ifelse(fraction.of.bin != 0,
                                              as.numeric(read.counts[, (4 + length(sample.files) + control)]) / fraction.of.bin,
                                              0)
      }
    } else {
      stop("Bins do not correspond")
    }
  }
  
  colnames(read.counts) <- c("BinID", "Chromosome", "StartPos", "StopPos",
                            sub(pattern = "$", ".compensated", paste(sample.files)),
                            paste(sample.files), sub(pattern = "$",
                                                 ".fractionOfBin",
                                                 paste(sample.files)))

  ## Create output file
  output <- read.counts
  output[, 2] <- gsub(prefixes[1], "", output[, 2])
  output[, 2] <- gsub("X", nchrom - 1, output[, 2])
  output[, 2] <- gsub("Y", nchrom, output[, 2])
  write.table(output, file = paste0(destination.folder,
                                    "read_counts_compensated.txt"),
                                    row.names = FALSE, col.names = TRUE,
                                    sep = "\t")

  ## Create histograms of fraction.of.bin
  ## (fraction of length in bins (after peak region removal)
  dir.create(paste0(destination.folder, "qc/"))
  for(i in 1:length(sample.files)) {
    pdf(paste0(destination.folder, "qc/fraction.of.bin_", i,
               ".pdf"), width=7, height=7)
    plot(ecdf(as.numeric(read.counts[,4+(2*length(sample.files))+i])),
         verticals = TRUE, ylab = "Fraction of bins", 
         xlab = "Remaining fraction of bin after peak removal",
         main = "Cumulative distribution of remaining bin fraction")
    dev.off()
  }

  #############################################
  ## Normalize for GC-content and mapability ##
  #############################################

  ## Read files
  read.counts <- read.counts[,1:(4 + length(sample.files))]
  for(i in 5:ncol(read.counts)) {
    write.table(read.counts[,c(2,3,4,i)], colnames(read.counts)[i], quote = FALSE,
                sep = "\t", row.names = FALSE, col.names = FALSE)
  }
  f <- list.files(pattern = "bam.compensated")
  
  gc <- read.delim(gc.content.file)[, c(1, 2, 3, 5)]
  colnames(gc) <- c("chr", "start", "end", "gc")
  mapa <- read.delim(mapability.file, header = FALSE,
                     col.names = c("chr", "start", "end", "mapa"))
  black <- read.delim(blacklist.file, header = FALSE,
                      col.names = c("chr", "start", "end"))
  
  data <- .loadCovData(f, gc = gc, mapa = mapa, black = black,
                       excludechr = "MT", datacol = 4)
  sampnames <- paste(sub("$", "", f), paste0(round(colSums(data$cov) / 1e6, 1),
                                             "M"), sep = "_")
  colnames(data$cov) <- sampnames
  usepoints <- !(data$anno$chr %in% c("X","Y","MT", "chrX", "chrY", "chrM"))
  
  ## Perform normalization (in .tng helper function)
  tryCatch({
    i <- c(1:ncol(data$cov))
    NormalizeDOC <- function(i, data, .tng, usepoints, destination.folder) {  
      .tng(data.frame(count = data$cov[,i], gc = data$anno$gc,
                      mapa = data$anno$mapa), use = usepoints,
           correctmapa = TRUE, plot = paste0(destination.folder,
                                             "qc/",
                                             colnames(data$cov)[i], ".png"))
    }
    sfInit(parallel=TRUE, cpus = ncpu)
    ratios <- sfLapply(i, NormalizeDOC, data, .tng, usepoints,
                       destination.folder)
    sfStop()
    ratios <- matrix(unlist(ratios), ncol = length(sample.files))
  }, error = function(e) {
    cat("ERROR: The GC-content and mapability normalization did not work due",
        "to a failure to calculate loesses.\n")
    cat("ERROR: This can generally be solved by using larger bin sizes.\n")
    stop("Stopping execution of the remaining part of the script...")    
  })
  
  colnames(ratios) <- sampnames
  
  rd <- list(ratios = ratios[!data$anno$black & !is.na(data$anno$mapa) &
                             data$anno$mapa > .2, ],
             anno = data$anno[!data$anno$black & !is.na(data$anno$mapa) &
                              data$anno$mapa > .2, ])
  
  sink(file = paste0(destination.folder, "log.txt"), append = TRUE,
       type = c("output", "message"))
  rd2 <- rd
  for(i in 1:length(sample.files)) {
    rd2$ratios[,i] <- rd$ratios[,i] - rd$ratios[,which.control[i]]
    cat("Relative log2-values are calculated for sample",
        colnames(rd2$ratios)[i], "with control",
        colnames(rd2$ratios)[which.control[i]], "\n")
  }
  colnames(rd2$ratios) <- gsub("$", ".rel", colnames(rd2$ratios))
  cat("\n\n")

  ###################
  ## Create output ##
  ###################
  
  ## Create table with corrected log2 values and write to file
  read.counts <- matrix(nrow = nrow(rd$ratios))
  read.counts <- cbind(read.counts, rd$anno[1:3])
  read.counts[,1] <- paste(rd$anno[,1], paste(rd$anno[,2], rd$anno[,3],
                                             sep = "-"), sep = ":")
  read.counts <- cbind(read.counts, rd$ratios, rd2$ratios)
  colnames(read.counts) <- c("BinID", "Chromosome", "StartPos", "StopPos",
                            colnames(read.counts[5:nrow(read.counts)]))
  
  read.counts <- read.counts[-which(rowSums(is.na(read.counts[, -c(1:4)])) > 0),]
  read.counts[, 2] <- gsub(prefixes[1], "", read.counts[, 2])
  read.counts[, 2] <- gsub("X", nchrom - 1, read.counts[, 2])
  read.counts[, 2] <- gsub("Y", nchrom, read.counts[, 2])
  read.counts[, 2] <- as.integer(read.counts[, 2])
  read.counts[read.counts == -Inf] <- -.Machine$integer.max/2
  read.counts[read.counts == Inf] <- .Machine$integer.max/2
  
  write.table(read.counts, paste0(destination.folder,
                                 "log2ratio_compensated_corrected.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  ##############################################################################
  ## Calculate overlap with capture.regions.file for quality control purpose ##
  ##############################################################################
  
  if (capture.regions.file != "not specified") {
    for(control.number in control.numbers) {
      intersection <- system(paste0("bedtools intersect -a ",
                                    capture.regions.file, " -b ",
                                    destination.folder,
                                    "BamBaiMacsFiles/MACS",
                                    control.number, "_peaks.bed -wao"),
                             intern = TRUE)
      intersection <- strsplit(intersection, "\t")
      intersection <- do.call(rbind, intersection)
      intersection <- intersection[,c(1,2,3,9)]
      vec <- vector()
      for (i in 2:nrow(intersection)) {
        if (intersection[i,1] == intersection[i-1,1] && intersection[i,2] == intersection[i-1,2] && intersection[i,3] == intersection[i-1,3] ) { ## Better check!
          intersection[i,4] <- as.integer(intersection[i,4]) + as.integer(intersection[i-1,4])
          vec <- append(vec, i-1)
        }
      }
      intersection <- intersection[-vec,]
      
      cat("Number of exons covered by peaks in sample ",
          sample.files[control.number], ": ", sum(intersection[,4] != "0"), "\n")
    
      intersection <- system(paste0("bedtools intersect -a ",
                                    destination.folder,
                                    "BamBaiMacsFiles/MACS",
                                    control.number, "_peaks.bed -b ",
                                    capture.regions.file, " -wao"),
                             intern = TRUE)
      intersection <- strsplit(intersection, "\t")
      intersection <- do.call(rbind, intersection)
      intersection <- intersection[,c(1,2,3,9)]
      vec <- vector()
      for (i in 2:nrow(intersection)) {
        if (intersection[i, 1] == intersection[i - 1, 1] && intersection[i, 2] == intersection[i - 1, 2] && intersection[i, 3] == intersection[i - 1, 3] ) { ## Better check!
          intersection[i,4] <- as.integer(intersection[i,4]) + as.integer(intersection[i-1,4])
          vec <- append(vec, i-1)
        }
      }
      intersection <- intersection[-vec,]
      
      cat("Number of peaks covered by exons in sample",
          sample.files[control.number], ": ", sum(intersection[,4] != "0"), "\n")
      cat("Total number of exons in sample", sample.files[control.number], ": ",
          nrow(read.table(capture.regions.file, sep = "\t")), "\n")
      cat("Total number of peaks", sample.files[control.number], ":",
          nrow(read.table(paste0(destination.folder,
                                 "BamBaiMacsFiles/MACS",
                                 control.number, "_peaks.bed"), sep = "\t")),
          "\n\n")  
    }
  }
  
  sink()
  ## Remove BamBaiMacsFiles folder
  # if (!keep.intermediairy.files) {
  #   unlink(paste0(destination.folder, "BamBaiMacsFiles/"))
  # }
  cat("Total calculation time: ", Sys.time() - start.time, "\n\n")
  
  inputStructure <- list(destination.folder = destination.folder, ncpu = ncpu,
                         nchrom = nchrom)
  save(inputStructure, file = paste0(destination.folder, "input.Rdata"))

}
