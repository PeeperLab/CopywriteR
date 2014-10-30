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
  sample.control <- apply(sample.control, c(1, 2), tools::file_path_as_absolute)
  sample.control <- data.frame(sample.control, stringsAsFactors = FALSE)
  colnames(sample.control) <- c("samples", "controls")
  destination.folder <- tools::file_path_as_absolute(destination.folder)
  reference.folder <- tools::file_path_as_absolute(reference.folder)
  
  ## Add trailing / to folder paths
  destination.folder <- paste0(destination.folder, "/")
  reference.folder <- paste0(reference.folder, "/")

  ## Check the existence of folders and files
  invisible(apply(sample.control, c(1, 2), function(x) {
    if (!file.exists(x)) {
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
  control.indices <- match(sample.control$controls, sample.control$samples)

  ## Create paths to helper files
  bin.file <- paste0(reference.folder, "bins.bed")
  blacklist.file <- paste0(reference.folder, "blacklist.bed")
  gc.content.file <- paste0(reference.folder, "GC_content.bed")
  mapability.file <- paste0(reference.folder, "mapability.bed")

  ## Retrieve number of chromosomes and bin size from bin.bed helper file
  bin.bed <- read.table(file = bin.file, as.is = TRUE, sep = "\t")
  colnames(bin.bed) <- c("Chromosome", "Start", "End")
  nchrom <- length(unique(bin.bed$Chromosome))
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

  cat("The following samples will be analyzed:\n")
  cat(paste0("sample: ", sample.files,
             ";\tmatching control: ", sample.files[control.indices]),
      sep = "\n")
  cat("The bin size for this analysis is", bin.size, "\n")
  cat("The capture region file is", capture.regions.file, "\n")
  cat("This analysis will be run on", ncpu, "cpus", "\n\n\n")
  sink()
  
  cat("The following samples will be analyzed:\n")
  cat(paste0("sample: ", sample.files,
             ";\tmatching control: ", sample.files[control.indices]),
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
  }
  
  ## Garbage collection
  rm(samp, header, con, chr.names, bin.file, reference.folder)

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
  cat("\n\n")
  
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
  cat("\n\n")

  ## Garbage collection
  rm(Stats, is.paired.end, sample.paths, ProperReads,
     NumberPairedEndReads, IndexBam)
  
  ## Create list with numbers of controls
  control.uniq.indices <- unique(control.indices)
  
  ## Call peaks in .bam file of control sample
  Macs14 <- function(control.uniq.indices, sample.files) {
    system(paste0("macs14 -t " , sample.files[control.uniq.indices], " -n MACS",
                  control.uniq.indices, " -g hs --nolambda"))
    paste0("macs14 -t " , sample.files[control.uniq.indices], " -n MACS",
           control.uniq.indices, " -g hs --nolambda")
  }
  sfInit(parallel=TRUE, cpus = min(length(control.uniq.indices), ncpu))
  to.log <- sfSapply(control.uniq.indices , Macs14, sample.files)
  sfStop()
  cat(to.log, "\n", sep = "\n")
  
  ## Read count statistics
  Stats <- function(sample.files, bin.bed) {
    all.reads <- countBam(sample.files)$records
    which <- with(bin.bed,
                  reduce(GRanges(seqnames = Chromosome,
                                 ranges = IRanges(start = Start, end = End))))
    what <- c("pos")
    param <- ScanBamParam(which = which, what = what)
    chrom.reads <- countBam(file = sample.files, param = param)
    chrom.reads <- sum(chrom.reads$records)
    c(all.reads = all.reads, chrom.reads = chrom.reads)
  }
  sfInit(parallel=TRUE, cpus = ncpu)
  sfLibrary(Rsamtools)
  res <- data.frame(t(sfSapply(sample.files, Stats, bin.bed)))
  sfStop()
  statistics <- within(statistics, {
    total.properreads<- res$all.reads
    unmapable.or.mitochondrial <- res$all.reads - res$chrom.reads
    on.chromosomes <- res$chrom.reads
  })
  print(statistics)
  cat("\n\n")
  
  ## Alternative for bedtools
  # Create GRanges file for bins
  bin.grange <- with(bin.bed,
                     GRanges(seqnames = Chromosome,
                             ranges = IRanges(start = Start, end = End)))
  
  i <- c(1:length(sample.files))
  CalculateDepthOfCoverage <- function(i, sample.files, control.indices,
                                       bin.grange, bin.size) {
    # Create GRanges object of peak files
    bed <- read.table(file = paste0("MACS", control.indices[i], "_peaks.bed"),
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
    sample.files[i] <- gsub("_properreads.bam$", ".bam", sample.files[i])
    
    if (all.equal(counts$space, as(seqnames(bin.grange), "factor"))) {
      counts <- within(counts, {
        Chromosome <- as(seqnames(bin.grange), "factor")
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
    return(list(counts[, paste0("read.counts.compensated.", sample.files[i]),
                       drop = FALSE],
                counts[, paste0("read.counts.", sample.files[i]), drop = FALSE],
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
                  control.indices, bin.grange, bin.size)
  sfStop()
  read.counts <- data.frame(Chromosome = as(seqnames(bin.grange), "factor"),
                            Start = ranges(bin.grange)@start,
                            End = ranges(bin.grange)@start +
                              ranges(bin.grange)@width - 1L)
  read.counts <- within(read.counts, {
    Feature = paste0(Chromosome, ":", paste0(Start, "-", End))
  })
  read.counts <- cbind(read.counts[, ], Reduce(cbind, res[1, ]),
                       Reduce(cbind, res[2, ]), Reduce(cbind, res[3, ]))
  cat(unlist(res[4, ]), "\n", sep = "\n")
   
  write.table(read.counts, file = paste0(destination.folder,
                                    "read_counts.txt"),
                                    row.names = FALSE, col.names = TRUE,
                                    sep = "\t")

  ## Create histograms of fraction.of.bin
  ## (fraction of length in bins (after peak region removal)
  sample.files <- gsub("_properreads.bam$", ".bam", sample.files)
  
  dir.create(paste0(destination.folder, "qc/"))
  for(i in 1:length(sample.files)) {
    pdf(paste0(destination.folder, "qc/fraction.of.bin.", sample.files[i],
               ".pdf"), width=7, height=7)
      plot(ecdf(as.numeric(read.counts[,4+(2*length(sample.files))+i])),
           verticals = TRUE, ylab = "Fraction of bins", 
           xlab = "Remaining fraction of bin after peak removal",
           main = "Cumulative distribution of remaining bin fraction")
    dev.off()
  }
  
  ## Garbage collection
  rm(statistics, bin.bed, bin.grange, to.log, res, Macs14, control.indices,
     CalculateDepthOfCoverage, bin.size)

  #############################################
  ## Normalize for GC-content and mapability ##
  #############################################

  ## Read files
  read.counts <- read.counts[,1:(4 + length(sample.files))]
  for(i in 5:ncol(read.counts)) {
    write.table(read.counts[, c(1, 2, 3, i)], colnames(read.counts)[i],
                quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  }
  f <- list.files(pattern = "read.counts.compensated")
  
  gc <- read.delim(gc.content.file)[, c(1, 2, 3, 5)]
  colnames(gc) <- c("chr", "start", "end", "gc")
  mapa <- read.delim(mapability.file, header = FALSE,
                     col.names = c("chr", "start", "end", "mapa"))
  black <- read.delim(blacklist.file, header = FALSE,
                      col.names = c("chr", "start", "end"))
  
  data <- .loadCovData(f, gc = gc, mapa = mapa, black = black,
                       excludechr = "MT", datacol = 4)
  colnames(data$cov) <- f
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
    log2.read.counts <- matrix(unlist(ratios), ncol = length(sample.files))
  }, error = function(e) {
    cat("ERROR: The GC-content and mapability normalization did not work due",
        "to a failure to calculate loesses.\n")
    cat("ERROR: This can generally be solved by using larger bin sizes.\n")
    stop("Stopping execution of the remaining part of the script...")    
  })
  
  colnames(log2.read.counts) <- paste0("log2.", gsub(".compensated", "", f))
  
  ###################
  ## Create output ##
  ###################
  selection <- !data$anno$black & !is.na(data$anno$mapa) & data$anno$mapa > .2
  log2.read.counts <- data.frame(data$anno[selection, c("chr", "start", "end")],
                                 log2.read.counts[selection, ])
  
  ## Create table with corrected log2 values and write to file
  log2.read.counts <- cbind(with(log2.read.counts, cbind(
    Chromosome = chr,
    Start = start,
    End = end,
    Feature = paste0(chr, ":", paste0(start, "-", end))
  )), log2.read.counts[, 4:ncol(log2.read.counts)])
  
  write.table(log2.read.counts, paste0(destination.folder,
                                       "log2_read_counts.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  ## Garbage collection
  rm(f, gc, mapa, black, data, usepoints, selection, blacklist.file,
     mapability.file, gc.content.file, read.counts, ratios, NormalizeDOC,
     log2.read.counts, i)

  #############################################################################
  ## Calculate overlap with capture.regions.file for quality control purpose ##
  #############################################################################
  
  if (capture.regions.file != "not specified") {
    captured.bed <- read.table(capture.regions.file, as.is = TRUE, sep = "\t")
    colnames(captured.bed) <- c("Chromosome", "Start", "End")
    captured.grange <- with(captured.bed,
                            GRanges(seqnames = Chromosome,
                                    ranges = IRanges(start = Start, end = End)))
    for (control.index in control.uniq.indices) {
      peak.bed <- read.table(file = paste0("MACS", control.index, "_peaks.bed"),
                             as.is = TRUE, sep = "\t")
      colnames(peak.bed) <- c("Chromosome", "Start", "End")
      peak.grange <- with(peak.bed, GRanges(seqnames = Chromosome,
                        ranges = IRanges(start = Start, end = End)))

      overlap <- findOverlaps(captured.grange, peak.grange)
      
      cat("Number of capture regions covered by peaks in sample",
          length(unique(queryHits(overlap))), "\n")
      cat("Number of peaks covered by capture regions in sample",
          length(unique(subjectHits(overlap))), "\n")          
      cat("Total number of capture regions in sample",
          sample.files[control.index], ": ", length(captured.grange), "\n")
      cat("Total number of peaks", sample.files[control.index], ":",
          length(peak.grange), "\n\n")
    }
  }
  
  sink()
  ## Remove BamBaiMacsFiles folder
  # if (!keep.intermediairy.files) {
  #   unlink(paste0(destination.folder, "BamBaiMacsFiles/"))
  # }
  
  ## Garbage collection
  rm(overlap, captured.bed, peak.bed, capture.regions.file, peak.grange,
     sample.files, control.uniq.indices, control.index, captured.grange)
  
  cat("Total calculation time: ", Sys.time() - start.time, "\n\n")
  
  inputStructure <- list(sample.control = sample.control,
                         ncpu = ncpu,
                         nchrom = nchrom,
                         prefix = prefixes[1])
  save(inputStructure, file = paste0(destination.folder, "input.Rdata"))

}