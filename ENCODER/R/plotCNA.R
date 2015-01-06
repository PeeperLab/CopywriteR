plotCNA <- function(destination.folder, smoothed = TRUE, sample.plot, y.min,
                    y.max, ...) {
  
  #############Change 20 kb bins
  
  ## Make folder path absolute
  destination.folder <- tools::file_path_as_absolute(destination.folder)

  ## Add trailing / to folder paths
  destination.folder <- paste0(destination.folder, "/CNAprofiles/")

  ## Check the existence of folders and files
  if (!file.exists(destination.folder)) {
    stop("The destination folder could not be found. Please change your ",
         "destination.folder path.")
  }

  load(paste0(destination.folder, "input.Rdata"), .GlobalEnv)
  
  ## Set variables
  nchrom <- inputStructure$nchrom
  prefix <- inputStructure$prefix
  ncpu <- inputStructure$ncpu
  bin.file <- inputStructure$bin.file
  bin.size <- inputStructure$bin.size
  
  ## Set sample.plot
  if (missing(sample.plot)) {
    sample.plot <- apply(inputStructure$sample.control, c(1, 2), function(x) {
      x <- paste0("log2.read.counts.", gsub(".*/", "", x))
    })
    all.samples <- unique(as.vector(sample.plot))
    sample.plot <- rbind(data.frame(sample.plot, stringsAsFactors = FALSE),
                         data.frame(samples = all.samples,
                                    controls = rep(NA, length(all.samples)),
                                    stringsAsFactors = FALSE))
  } else {
    sample.plot[, ] <- apply(sample.plot, c(1, 2), function(x) {
      if (!is.na(x)) {
        x <- paste0("log2.read.counts.",
                    gsub("_properreads", "", gsub(".*/", "", x)))
      } else {
        x <- as.character(x)
      }
    })
  }
  
  ## Read data
  log2.read.counts <- read.table(file = paste0(destination.folder,
                                               "log2_read_counts.txt"),
                                 sep = "\t", header = TRUE, check.names = FALSE)

  ## Remove prefix and convert chromosome names to integers
  log2.read.counts <-
    log2.read.counts[-which(rowSums(is.na(log2.read.counts[, -c(1:4)])) > 0), ]
  log2.read.counts$Chromosome <- gsub(prefix, "", log2.read.counts$Chromosome)
  log2.read.counts$Chromosome <-
    gsub("X", nchrom - 1, log2.read.counts$Chromosome)
  log2.read.counts$Chromosome <- gsub("Y", nchrom, log2.read.counts$Chromosome)
  log2.read.counts$Chromosome <- as.integer(log2.read.counts$Chromosome)
  
  ## Fix behaviour of DNAcopy with 'outlier' values
  log2.read.counts[, 5:ncol(log2.read.counts)] <-
    apply(log2.read.counts[, 5:ncol(log2.read.counts), drop = FALSE],
          c(1, 2), function(x) {
    if (x < -5) {
      x <- -5
    } else if (x > 10) {
      x <- 10
    } else {
      x <- x
    }
  })
  
  ## Create table with values to be plotted
  if (all(na.omit(unlist(sample.plot)) %in% colnames(log2.read.counts))) {
    plotting.values <- within(log2.read.counts, {
      # Loop through sample.plot calculating the relevant absolute and relative
      # log2.read.counts values. Loop used for readability. Loop in reverse
      # order to get the ordering in data.frame correct
      for (i in nrow(sample.plot):1) {
        if (!is.na(sample.plot$controls[i])) {
          assign(paste0(sample.plot$samples[i], ".vs.",
                        sample.plot$controls[i]),
                 log2.read.counts[, sample.plot$samples[i]] -
                 log2.read.counts[, sample.plot$controls[i]])
        } else {
          assign(paste0(sample.plot$samples[i], ".vs.none"),
                 log2.read.counts[, sample.plot$samples[i]])
        }
      }
    rm(list = c("i", colnames(log2.read.counts[5:ncol(log2.read.counts)])))
    })
  } else {
    stop("One of the samples in sample.plot refers to a BAM file that has not ",
         "been processed in ENCODER. Please make sure that you have provided ",
         "the correct input files or re-run ENCODER accordingly.")
  }
  
  ## Apply DNAcopy
  CNA.object <- CNA(plotting.values[, 5:ncol(plotting.values)],
                    plotting.values$Chromosome,
                    rowMeans(plotting.values[, c("Start", "End")]),
                    data.type = "logratio",
                    sampleid =
                      colnames(plotting.values)[5:ncol(plotting.values)])
  if (smoothed) {
    CNA.object <- smooth.CNA(CNA.object)  
  }
  segment.CNA.object <- segment(CNA.object, verbose = 1, ...)
  save(segment.CNA.object, file = paste0(destination.folder, "segment.Rdata"))
  
  ## Calculate the chromosome lengths from the bin.bed file
  bin.bed <- read.table(file = bin.file, as.is = TRUE, sep = "\t")
  colnames(bin.bed) <- c("Chromosome", "Start", "End")
  bin.bed$Chromosome <- gsub(prefix, "", bin.bed$Chromosome)
  bin.bed$Chromosome <- gsub("X", nchrom - 1, bin.bed$Chromosome)
  bin.bed$Chromosome <- gsub("Y", nchrom, bin.bed$Chromosome)
  bin.bed$Chromosome <- as.integer(bin.bed$Chromosome)
  bin.bed <- with(bin.bed,
                  reduce(GRanges(seqnames = Chromosome,
                                 ranges = IRanges(start = Start, end = End))))
  bin.bed <- bin.bed[order(as.integer(seqlevels(bin.bed)))]
  chrom.lengths <- data.frame(Chromosome = 1:length(seqlevels(bin.bed)),
                              Length = ranges(bin.bed)@width - 1)
  chrom.lengths <- within(chrom.lengths, {
    CumSum <- c(0, cumsum(Length)[1:(nrow(chrom.lengths) - 1)])
  })
  # save(chrom.lengths, file = paste0(destination.folder, "chrom_lengths.Rdata"))

  ## Create plots
  segment.CNA.object$output <- within(segment.CNA.object$output, {
    start.position.chrom <-
      chrom.lengths$CumSum[match(as.integer(chrom), chrom.lengths$Chromosome)]
    loc.start <- loc.start + start.position.chrom
    loc.end <- loc.end + start.position.chrom
    rm(start.position.chrom)
  })
  segment.CNA.object$data <- within(segment.CNA.object$data, {
    start.position.chrom <-
      chrom.lengths$CumSum[match(as.integer(chrom), chrom.lengths$Chromosome)]
    maploc <- maploc + start.position.chrom
    rm(start.position.chrom)
  })

  # Get sample names
  samples <- colnames(segment.CNA.object$data)
  samples <- samples[3:length(samples)]

  # Add trailing / to folder paths
  plot.folder <- paste0(destination.folder, "plots/")
  dir.create(plot.folder)
  setwd(plot.folder)
  
  # Set boundaries plots
  if (missing(y.min)) {
    y.min <- -3
  }
  if (missing(y.max)) {
    y.max <- 2
  }

  # Loop through samples using lapply
  invisible(lapply(samples, function(x) {
    # Select sample
    select.sample <- x
    current.sample <- segment.CNA.object
    current.sample$data <- current.sample$data[, c("chrom", "maploc",
                                                   select.sample)]
    current.sample$output <-
      current.sample$output[current.sample$output$ID == select.sample, ]
  
    # Create and set new directory
    dir.create(paste0(plot.folder, select.sample))
    setwd(paste0(plot.folder, select.sample))
  
    # Loop through chromosomes using lapply
    invisible(lapply(c(as.list(1:nchrom), list(1:nchrom)), function(x) {
      # Select chromosome
      select.chrom <- x
      current.sample$data <-
        current.sample$data[current.sample$data$chrom %in% select.chrom, ]
      current.sample$output <-
        current.sample$output[current.sample$output$chrom %in% select.chrom, ]
  
      # Set variables
      genome.position.min <-
        chrom.lengths$CumSum[match(min(select.chrom), chrom.lengths$Chromosome)]
      genome.position.max <-
        chrom.lengths$CumSum[match(max(select.chrom), chrom.lengths$Chromosome)] + 
        chrom.lengths$Length[match(max(select.chrom), chrom.lengths$Chromosome)]
    
      # Plot data
      if (length(x) > 1) {
        name.chrom <- "all.chrom"
      } else {
        name.chrom <- paste0("chrom.", x)
      }
      pdf(file = paste0(name.chrom, ".pdf"), width = 14, height = 7)
        plot(current.sample$data[, "maploc"],
             current.sample$data[, select.sample], ylim = c(y.min, y.max),
             pch = ".", cex = 1,
             xlim = c(genome.position.min, genome.position.max), xaxs = "i",
             xlab = "", ylab = "log2 value", xaxt = "n", main = select.sample)
      
        points(current.sample$data$maploc[current.sample$data[, select.sample] > y.max],
               rep.int(y.max, sum(current.sample$data[, select.sample] > y.max)),
               cex = 0.25, col = "green", pch = 2)
        points(current.sample$data$maploc[current.sample$data[, select.sample] < y.min],
               rep.int(y.min, sum(current.sample$data[, select.sample] < y.min)),
               cex = 0.25, col = "red", pch = 6)
            
        invisible(apply(current.sample$output, 1, function(x) {
          segments(x0 = as.numeric(x["loc.start"]),
                   y0 = as.numeric(x["seg.mean"]),
                   x1 = as.numeric(x["loc.end"]),
                   col = "red", pch = ".", cex = 0.2)
        }))
        par(xpd = TRUE)
        text(x = 0.96 * (genome.position.max - genome.position.min) +
             genome.position.min,
             y = 0.075 * (y.max - y.min) + y.max,
             labels = paste0("mad = ", round(madDiff(current.sample$data[, select.sample]), 3)))
        text(x = 0.035 * (genome.position.max - genome.position.min) +
             genome.position.min,
             y = 0.075 * (y.max - y.min) + y.max,
             labels = paste0(bin.size / 1000, " kb bins"))
        par(xpd = FALSE)
        ticks <- (chrom.lengths$CumSum + chrom.lengths$Length / 2)[select.chrom]
        axis(1, at = ticks, labels = select.chrom)
        if (length(select.chrom) > 1) {
          abline(v = chrom.lengths$CumSum[2:nrow(chrom.lengths)],
                 lty = "dotted")
        }
      dev.off()
    }))
  }))
}
