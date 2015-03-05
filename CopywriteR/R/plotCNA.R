plotCNA <- function(destination.folder, smoothed = TRUE, sample.plot, y.min,
                    y.max, ...) {

    start.time <- Sys.time()

    ## Make destination folder path absolute
    destination.folder <- tools::file_path_as_absolute(destination.folder)
    destination.folder <- file.path(destination.folder, "CNAprofiles")

    ## Provide output to log
    flog.appender(appender.file(file.path(destination.folder,
                                          "CopywriteR.log")))
    args.ellipsis <- unlist(list(...))
    flog.info(paste0("plotCNA was run using the following commands:", "\n\n",
                     "plotCNA(destination.folder = \"",
                     dirname(destination.folder),
                     "\", smoothed = ", smoothed,
                     if (!missing(sample.plot)) {", sample.plot = sample.plot"},
                     if (!missing(y.min)) {paste0(", y.min = ", y.min)},
                     if (!missing(y.max)) {paste0(", y.max = ", y.max)},
                     if (!missing(...)) {paste0(", ", paste(names(args.ellipsis),
                                                "=", args.ellipsis,
                                                collapse = ", "))},
                     ")"))

    ## Check the existence of folders and files
    if (!file.exists(destination.folder)) {
        stop(.wrap("The destination folder could not be found. Please change",
                   "the path specified in", sQuote(destination.folder)))
    }

    # Pre-define inputStructure variable to avoid it from raising NOTES during
    # R CMD CHECK (variable is loaded from file below)
    inputStructure <- NULL

    # Load inputStructure variable
    load(file.path(destination.folder, "input.Rdata"))

    ## Set variables
    chromosomes <- inputStructure$chromosomes
    nautosomes <- length(grep("[0-9]", chromosomes))
    prefix <- inputStructure$prefix
    bin.size <- inputStructure$bin.size

    ## Set sample.plot
    if (missing(sample.plot)) {
        sample.plot <- apply(inputStructure$sample.control, c(1, 2),
                             function(x) {
            x <- paste0("log2.", basename(x))
        })
        all.samples <- unique(as.vector(sample.plot))
        sample.plot <- rbind(data.frame(sample.plot, stringsAsFactors = FALSE),
                             data.frame(samples = all.samples,
                                        controls = rep(NA, length(all.samples)), 
                                        stringsAsFactors = FALSE))
    } else {
        colnames(sample.plot) <- c("samples", "controls")
        sample.plot[, ] <- apply(sample.plot, c(1, 2), function(x) {
            if (!is.na(x)) {
                x <- paste0("log2.read.counts.",
                            gsub("_properreads", "", basename(x)))
            } else {
                x <- as.character(x)
            }
        })
    }

    ## Read data
    log2.read.counts <- read.table(file = file.path(destination.folder,
                                                    "log2_read_counts.igv"),
                                   sep = "\t", header = TRUE,
                                   check.names = FALSE)

    ## Remove prefix and convert chromosome names to integers
    log2.read.counts$Chromosome <- gsub(prefix, "", log2.read.counts$Chromosome)
    chromosomes <- gsub(prefix, "", chromosomes)
    
    log2.read.counts$Chromosome <- gsub("X", nautosomes + 1,
                                        log2.read.counts$Chromosome)
    chromosomes <- gsub("X", nautosomes + 1, chromosomes)
    
    log2.read.counts$Chromosome <- gsub("Y", nautosomes + 2,
                                        log2.read.counts$Chromosome)
    chromosomes <- gsub("Y", nautosomes + 2, chromosomes)
    
    log2.read.counts$Chromosome <- as.integer(log2.read.counts$Chromosome)
    chromosomes <- sort(as.integer(chromosomes))

    ## Fix behaviour of DNAcopy with 'outlier' values
    log2.read.counts[, 5:ncol(log2.read.counts)] <-
        apply(log2.read.counts[, 5:ncol(log2.read.counts), drop = FALSE],
              c(1, 2), function(x) {
        if (is.na(x)) {
            x <- NA
        } else if (x < -5) {
            x <- -5
        } else if (x > 10) {
            x <- 10
        } else {
            x <- x
        }
    })

    ## Create table with values to be plotted
    if (all(na.omit(unlist(sample.plot)) %in% colnames(log2.read.counts))) {
        plotting.values <-
            log2.read.counts[, c("Chromosome", "Start", "End", "Feature")]
        for (i in seq_len(nrow(sample.plot))) {
            if (!is.na(sample.plot$controls[i])) {
                plotting.values[, paste0(sample.plot$samples[i], ".vs.",
                                         sample.plot$controls[i])] <- 
                    log2.read.counts[, sample.plot$samples[i]] -
                    log2.read.counts[, sample.plot$controls[i]]
            } else {
                plotting.values[, paste0(sample.plot$samples[i], ".vs.none")] <-
                    log2.read.counts[, sample.plot$samples[i]]
            }
        }
    } else {
        stop(.wrap("One of the samples in", sQuote(sample.plot), "refers to a",
                   "BAM file that has not been processed in CopywriteR. Please",
                   "make sure that you have provided the correct input files",
                   "or re-run CopywriteR accordingly."))
    }

    ## Apply DNAcopy
    CNA.object <- CNA(plotting.values[, 5:ncol(plotting.values), drop = FALSE],
                      plotting.values$Chromosome,
                      rowMeans(plotting.values[, c("Start", "End")]),
                      data.type = "logratio", 
                      sampleid = colnames(plotting.values)[5:ncol(plotting.values)])
    if (smoothed) {
        CNA.object <- smooth.CNA(CNA.object)
    }
    segment.CNA.object <- segment(CNA.object, verbose = 1, ...)
    save(segment.CNA.object, file = file.path(destination.folder,
                                              "segment.Rdata"))

    ## Calculate the chromosome lengths from the bin.bed file
    chrom.lengths <-
        scanBamHeader(inputStructure$sample.control$samples[1])[[1]]$targets
    names(chrom.lengths) <- gsub(prefix, "", names(chrom.lengths))    
    chrom.lengths <- data.frame(Chromosome = names(chrom.lengths),
                                Length = chrom.lengths, row.names = NULL)
    suppressWarnings(chrom.lengths <- chrom.lengths[chrom.lengths$Chromosome %in% c("X", "Y") |
                                                    !is.na(as.integer(chrom.lengths$Chromosome)), ])
    chrom.lengths <- chrom.lengths[mixedorder(chrom.lengths$Chromosome), ]
    chrom.lengths$Chromosome <- gsub("X", nautosomes + 1,
                                     chrom.lengths$Chromosome)
    chrom.lengths$Chromosome <- gsub("Y", nautosomes + 2,
                                     chrom.lengths$Chromosome)
    chrom.lengths[, "CumSum"] <- c(0, cumsum(as.numeric(chrom.lengths$Length))[seq_len(nrow(chrom.lengths) - 1)])
    
    ## Create plots
    segment.CNA.object$output[, "start.position.chrom"] <-
        chrom.lengths$CumSum[match(as.integer(segment.CNA.object$output$chrom),
                                   chrom.lengths$Chromosome)]
    segment.CNA.object$output$loc.start <- segment.CNA.object$output$loc.start + segment.CNA.object$output$start.position.chrom
    segment.CNA.object$output$loc.end <- segment.CNA.object$output$loc.end + segment.CNA.object$output$start.position.chrom
    segment.CNA.object$output$start.position.chrom <- NULL
    
    segment.CNA.object$data[, "start.position.chrom"] <-
        chrom.lengths$CumSum[match(as.integer(segment.CNA.object$data$chrom),
                                   chrom.lengths$Chromosome)]
    segment.CNA.object$data$maploc <- segment.CNA.object$data$maploc + segment.CNA.object$data$start.position.chrom
    segment.CNA.object$data$start.position.chrom <- NULL

    # Get sample names
    samples <- colnames(segment.CNA.object$data)
    samples <- samples[3:length(samples)]

    # Create plots folder
    plot.folder <- file.path(destination.folder, "plots")
    dir.create(plot.folder)
    setwd(plot.folder)

    # Set boundaries plots
    if (missing(y.min)) {
        y.min <- -2
    }
    if (missing(y.max)) {
        y.max <- 5
    }

    # Loop through samples using lapply
    invisible(lapply(samples, function(x) {
        # Select sample
        select.sample <- x
        current.sample <- segment.CNA.object
        current.sample$data <-
            current.sample$data[, c("chrom", "maploc", select.sample)]
        current.sample$output <-
            current.sample$output[current.sample$output$ID == select.sample, ]

        # Create and set new directory
        dir.create(file.path(plot.folder, select.sample))
        setwd(file.path(plot.folder, select.sample))

        # Loop through chromosomes using lapply
        invisible(lapply(c(as.list(chromosomes), list(chromosomes)),
                         function(x) {
            # Select chromosome
            select.chrom <- x
            current.sample$data <- current.sample$data[current.sample$data$chrom %in% select.chrom, ]
            current.sample$data <- current.sample$data[!is.na(current.sample$data[, select.sample]), ]
            current.sample$output <- current.sample$output[current.sample$output$chrom %in% select.chrom, ]

            # Set variables
            genome.position.min <-
                chrom.lengths$CumSum[match(min(select.chrom),
                                           chrom.lengths$Chromosome)]
            genome.position.max <-
                chrom.lengths$CumSum[match(max(select.chrom),
                                           chrom.lengths$Chromosome)] +
                                     chrom.lengths$Length[match(max(select.chrom), chrom.lengths$Chromosome)]

            # Plot data
            if (length(x) > 1) {
                name.chrom <- "all_chrom"
            } else {
                name.chrom <- paste0("chrom_", x)
            }
            pdf(file = paste0(name.chrom, ".pdf"), width = 14, height = 7)
            plot(current.sample$data[, "maploc"],
                 current.sample$data[, select.sample], ylim = c(y.min, y.max),
                 pch = ".", cex = 1, xlim = c(genome.position.min, 
                                              genome.position.max),
                 xaxs = "i", xlab = "", ylab = "log2 value", xaxt = "n",
                 main = select.sample)

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
            text(x = 0.96 * (genome.position.max - genome.position.min) + genome.position.min,
                 y = 0.075 * (y.max - y.min) + y.max,
                 labels = paste0("mad = ", round(madDiff(current.sample$data[, select.sample]), 3)))
            text(x = 0.035 * (genome.position.max - genome.position.min) + genome.position.min,
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
    flog.info(paste("Total calculation time of plotCNA was",
                    round(difftime(Sys.time(), start.time, units = "mins"), 2),
                    "minutes"))
    cat("Total calculation time of CopywriteR was: ",
        Sys.time() - start.time, "\n\n")
}
