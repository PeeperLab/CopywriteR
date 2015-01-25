ENCODER <- function(sample.control, destination.folder, reference.folder, ncpu,
                    capture.regions.file, keep.intermediairy.files = FALSE) {
  
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
  
  ## Make folder paths absolute
  sample.control <- apply(sample.control, c(1, 2), tools::file_path_as_absolute)
  sample.control <- data.frame(sample.control, stringsAsFactors = FALSE)
  colnames(sample.control) <- c("samples", "controls")
  destination.folder <- tools::file_path_as_absolute(destination.folder)
  reference.folder <- tools::file_path_as_absolute(reference.folder)
  
  ## Check the existence of folders and files
  invisible(apply(sample.control, c(1, 2), function(x) {
    if (!file.exists(x)) {
      stop("The file ", x, " could not be found.\nPlease change the path to ",
           "this file.")
    }
  }))
  
  if (!file.exists(destination.folder)) {
    stop("The destination.folder could not be found.\nPlease change your ",
         "destination.folder path.")
  }
  
  if (!file.exists(reference.folder)) {
    stop("The reference.folder could not be found.\nPlease change your ",
         "reference.folder path or run `preENCODER` to generate the required ",
         "folder with GC-content and mapability files for your desired bin ",
         "size.")
  }
  
  if (!file.exists(capture.regions.file) &
      capture.regions.file != "not specified") {
    stop("The capture.regions.file could not be found.\nPlease change your ",
         "capture.regions.file path.")
  }

  ## Check for write permissions in the destination folder
  if (file.access(destination.folder, 2) == -1) {
    stop("You do not have write permission in the destination folder.")
  }

  ## Create lists with BAM files and index of corresponding control
  sample.paths <- unlist(sample.control)
  sample.paths <- unique(sample.paths[!is.na(sample.paths)])
  sample.files <- unname(sapply(sample.paths, function(x) {
    x <- unlist(strsplit(x, "/"))[length(unlist(strsplit(x, "/")))]
  }))
  control.indices <- match(sample.control$controls, sample.control$samples)

  ## Create paths to helper files
  bin.file <- file.path(reference.folder, "bins.bed")
  blacklist.file <- file.path(reference.folder, "blacklist.bed")
  gc.content.file <- file.path(reference.folder, "GC_content.bed")
  mapability.file <- file.path(reference.folder, "mapability.bed")

  ## Retrieve number of chromosomes and bin size from bin.bed helper file
  bin.bed <- read.table(file = bin.file, as.is = TRUE, sep = "\t")
  colnames(bin.bed) <- c("Chromosome", "Start", "End")
  nchrom <- length(grep("[0-9]", unique(bin.bed$Chromosome)))
  bin.size <- bin.bed$End[1]
  
  ## Create folders
  destination.folder <- file.path(destination.folder, "CNAprofiles")
  tryCatch({
    if (!file.exists(file.path(destination.folder, "BamBaiMacsFiles"))) {
      dir.create(file.path(destination.folder, "BamBaiMacsFiles"),
                 recursive = TRUE)
    }
  }, warning = function(e) {
    stop("You do not have write permissions in the destination folder.\n",
         "Stopping execution of the remaining part of the script...")
  })
  
  ## Provide output for log file
  sink(file = file.path(destination.folder, "log.txt"),
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
      prefixes <- append(prefixes, gsub("[[:digit:]]|X|Y|M|T", "", chr.names)[1])
    }
  }, error = function(e) {
    stop("The BAM file header of file", samp, "is corrupted or truncated.\n",
         "Please rebuild this BAM file or exclude it from analysis.\n",
         "Stopping execution of the remaining part of the script...")    
  })

  if (!all(prefixes == prefixes[1])) {
    stop("The bam files have different chromosome names.\n",
         "Please adjust the .bam files such that they contain the same ",
         "chromosome notation.")
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
           "names.\n",
           "Please adjust the input files such that they contain the same ",
           "chromosome notation.")
    }
  }
  
  ## Garbage collection
  rm(chr.names, con, header, reference.folder, samp)

  ########################################################
  ## Calculate depth of coverage using off-target reads ##
  ########################################################

  sink(file = file.path(destination.folder, "log.txt"), append = TRUE,
       type = c("output", "message"))
    
  ## Create list of .bam files
  cat(sample.files, "\n", sep = "\n")

  ## Index .bam files
  if (!any(list.files(pattern = ".bam$") ==
           gsub(".bai$", "", list.files(pattern = ".bai$")))) {
    if (file.access(".", 2) == -1) {
      stop("The .bam files are not indexed and you do not have write ",
           "permission in (one of) the folder(s) where the .bam files ",
           "are located.")
    }
    IndexBam <- function(sample.paths) {
      indexBam(sample.paths)
      paste0("indexBam(\"", sample.paths, "\")")
    }
    sfInit(parallel=TRUE, cpus = ncpu)
    sfLibrary(Rsamtools)
    to.log <- sfSapply(sample.paths, IndexBam)
    sfStop()
    cat(to.log, "\n", sep = "\n")
    
    ## Garbage collection
    rm(IndexBam)
  }
  
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
                file.path(destination.folder, "BamBaiMacsFiles",
                          gsub(".bam$", "_properreads.bam", sample.files[i])),
                filter = filter, indexDestination = TRUE, param = param)
      paste0("filterBam(\"", sample.paths[i], "\", \"",
             file.path(destination.folder, "BamBaiMacsFiles",
                       gsub(".bam$", "_properreads.bam", sample.files[i])),
             "\", filter = filter, ", "indexDestination = TRUE, param = param)")
    } else {
      param <- ScanBamParam(what = "mapq")
      filter <- FilterRules(list(isHighQual = function(x) {
                                   x$mapq >= 37
                                 }))
      filterBam(sample.paths[i],
                file.path(destination.folder, "BamBaiMacsFiles/",
                          gsub(".bam$", "_properreads.bam", sample.files[i])),
                filter = filter, indexDestination = TRUE, param = param)
      paste0("filterBam(\"", sample.paths[i], "\", \"",
             file.path(destination.folder, "BamBaiMacsFiles",
                       gsub(".bam$", "_properreads.bam", sample.files[i])),
             "\", filter = filter, ", "indexDestination = TRUE, param = param)")
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
  setwd(file.path(destination.folder, "BamBaiMacsFiles"))
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

  ## Garbage collection
  rm(is.paired.end, NumberPairedEndReads, ProperReads, sample.paths, Stats)
  
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
 
    # countBam on remainder of bins
    param <- ScanBamParam(which = reduce(outside.peak.grange), what = c("pos"))
    counts.ENCODER <- countBam(sample.files[i], param = param)
    counts.ENCODER <- sum(counts.ENCODER$records)
    
    # Fix MT as levels in factor seqnames (remainder from setdiff operation)
    counts$space <- as.factor(as.character(counts$space))
    
    # Replace bins by real bins & calculate compensated read counts
    sample.files[i] <- gsub("_properreads.bam$", ".bam", sample.files[i])
    
    # Aggregate counts and total length per bin
    counts.grange <- GRanges(seqnames = counts$space,
                             ranges = IRanges(start = counts$start,
                                              end = counts$end),
                             records = counts$records)
    overlaps <- findOverlaps(counts.grange, bin.grange, minoverlap = 2L)
    index <- subjectHits(overlaps)
    records <- counts.grange[queryHits(overlaps)]@elementMetadata@listData$records
    lengths <- ranges(pintersect(counts.grange[queryHits(overlaps)],
                                 bin.grange[subjectHits(overlaps)]))@width
    aggregate.data.table <- data.table(index, records, lengths)
    aggregate.data.table <- aggregate.data.table[, list(records = sum(records),
                                                         lengths = sum(lengths)),
                                                 by = c("index")]
    aggregate.data.table <- aggregate.data.table[match(1:length(bin.grange),
                                                       aggregate.data.table$index), ]
    aggregate.data.table$records[which(is.na(aggregate.data.table$index))] <- 0
    aggregate.data.table$lengths[which(is.na(aggregate.data.table$index))] <- 0

    # Replace bins by real bins & calculate compensated read counts
    counts <- within(aggregate.data.table, {
      Chromosome <- as(seqnames(bin.grange), "factor")
      Start <- ranges(bin.grange)@start
      End <- ranges(bin.grange)@start + ranges(bin.grange)@width - 1L
      Feature <- paste0(Chromosome, ":", paste0(Start, "-", End))
      assign(paste0("read.counts.", sample.files[i]), records)
      assign(paste0("read.counts.compensated.", sample.files[i]),
             records / (lengths / (bin.size + 1)))
      assign(paste0("fraction.of.bin.", sample.files[i]),
             lengths / (bin.size + 1))
      rm(index, lengths, records)
    })
    counts <- data.frame(counts, check.names = FALSE)
  
    # Return
    return(list(counts[, paste0("read.counts.compensated.", sample.files[i]),
                       drop = FALSE],
                counts[, paste0("read.counts.", sample.files[i]), drop = FALSE],
                counts[, paste0("fraction.of.bin.", sample.files[i]),
                       drop = FALSE],
                paste0("Rsamtools finished calculating reads per bin in ",
                       "sample ", sample.files[i], "; number of bins = ",
                       nrow(counts)),
                counts.ENCODER))
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
                       
  ## Remove potential NAs introduced by peaks spanning entire bins
  read.counts[, 5:ncol(read.counts)] <-
    apply(read.counts[, 5:ncol(read.counts), drop = FALSE],
          c(1, 2), function(x) {
    if (is.na(x)) {
      x <- 0
    } else {
      x <- x
    }
  })

  cat(unlist(res[4, ]), "\n", sep = "\n")
  
  statistics <- within(statistics, {
    off.target <- unlist(res[5, ])
    on.target <- on.chromosomes - off.target
    rm(on.chromosomes)
  })
  print(statistics)
  cat("\n\n")

  ## Add counts.ENCODER
   
  write.table(read.counts, file = file.path(destination.folder,
                                    "read_counts.txt"),
                                    row.names = FALSE, col.names = TRUE,
                                    sep = "\t")

  ## Create histograms of fraction.of.bin
  ## (fraction of length in bins (after peak region removal)
  sample.files <- gsub("_properreads.bam$", ".bam", sample.files)
  
  dir.create(file.path(destination.folder, "qc"))
  for(i in 1:length(sample.files)) {
    pdf(file.path(destination.folder, "qc", paste0("fraction.of.bin.",
                                                sample.files[i], ".pdf")),
        width=7, height=7)
      plot(ecdf(as.numeric(read.counts[,4+(2*length(sample.files))+i])),
           verticals = TRUE, ylab = "Fraction of bins", 
           xlab = "Remaining fraction of bin after peak removal",
           main = "Cumulative distribution of remaining bin fraction")
    dev.off()
  }
  
  ## Garbage collection
  rm(bin.bed, bin.grange, CalculateDepthOfCoverage, control.indices, Macs14,
     res, statistics, to.log)

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
           correctmapa = TRUE, plot = file.path(destination.folder, "qc",
                                             paste0(colnames(data$cov)[i],
                                                    ".png")))
    }
    sfInit(parallel=TRUE, cpus = ncpu)
    ratios <- sfLapply(i, NormalizeDOC, data, .tng, usepoints,
                       destination.folder)
    sfStop()
    log2.read.counts <- matrix(unlist(ratios), ncol = length(sample.files))
  }, error = function(e) {
    stop("The GC-content and mapability normalization did not work due to a ",
        "failure to calculate loesses.\n",
        "This can generally be solved by using larger bin sizes.\n",
        "Stopping execution of the remaining part of the script...")
  })
  
  colnames(log2.read.counts) <- paste0("log2.", gsub(".compensated", "", f))
  
  ###################
  ## Create output ##
  ###################
  selection <- !data$anno$black & !is.na(data$anno$mapa) & data$anno$mapa > .2
  log2.read.counts <- data.frame(data$anno[selection, c("chr", "start", "end")],
                                 log2.read.counts[selection, ],
                                 check.names = FALSE)
  
  ## Create table with corrected log2 values and write to file
  log2.read.counts <- cbind(with(log2.read.counts, cbind(
    Chromosome = chr,
    Start = start,
    End = end,
    Feature = paste0(chr, ":", paste0(start, "-", end))
  )), log2.read.counts[, 4:ncol(log2.read.counts)])
  
  write.table(log2.read.counts, file.path(destination.folder,
                                       "log2_read_counts.igv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  ## Garbage collection
  rm(black, blacklist.file, data, f, gc, gc.content.file, i, log2.read.counts,
     mapa, mapability.file, NormalizeDOC, ratios, read.counts, selection,
     usepoints)

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
    
    ## Garbage collection
    rm(captured.bed, captured.grange, capture.regions.file, overlap, peak.bed,
       peak.grange)
  }
  
  sink()
  
  ## Remove BamBaiMacsFiles folder
  if (!keep.intermediairy.files) {
    unlink(file.path(destination.folder, "BamBaiMacsFiles"), recursive = TRUE)
  }
  
  ## Garbage collection
  rm(control.index, control.uniq.indices, sample.files)
  
  cat("Total calculation time: ", Sys.time() - start.time, "\n\n")
  
  inputStructure <- list(sample.control = sample.control,
                         ncpu = ncpu,
                         nchrom = nchrom,
                         prefix = prefixes[1],
                         bin.file = bin.file,
                         bin.size = bin.size)
  save(inputStructure, file = file.path(destination.folder, "input.Rdata"))

}