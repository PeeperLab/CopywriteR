CopywriteR <- function(sample.control, destination.folder, reference.folder,
                       bp.param, capture.regions.file,
                       keep.intermediary.files = FALSE) {

    ##########################
    ## Check and initialise ##
    ##########################

    start.time <- Sys.time()

    ## Restore work directory upon exit
    wd.orig <- getwd()
    on.exit(setwd(wd.orig))
    
    ## Make capture.regions.file path absolute if present
    if (missing(capture.regions.file)) {
        capture.regions.file <- "not specified"
    } else {
        capture.regions.file <- tools::file_path_as_absolute(capture.regions.file)
    }

    ## Make folder paths absolute
    sample.control <- apply(sample.control, c(1, 2),
                            tools::file_path_as_absolute)
    sample.control <- data.frame(sample.control, stringsAsFactors = FALSE)
    colnames(sample.control) <- c("samples", "controls")
    destination.folder <- tools::file_path_as_absolute(destination.folder)
    reference.folder <- tools::file_path_as_absolute(reference.folder)

    ## Check the existence of folders and files
    invisible(apply(sample.control, c(1, 2), function(x) {
        if (!file.exists(x)) {
            stop(.wrap("The file", sQuote(x), "could not be found. Please",
                       "change the path to this file."))
        }
    }))

    if (!file.exists(destination.folder)) {
        stop(.wrap("The destination folder could not be found. Please change",
                   "the path specified in", sQuote(destination.folder)))
    }

    if (!file.exists(reference.folder)) {
        stop(.wrap("The reference folder could not be found. Please change",
                   "the path specified in", sQuote(reference.folder), "or run",
                   "'preCopywriteR' to generate the required folder with",
                   "GC-content and mappability files for your desired bin",
                   "size."))
    }

    if (!file.exists(capture.regions.file) & capture.regions.file != "not specified") {
        stop(.wrap("The capture regions file could not be found. Please change",
                   "the path specified in", sQuote(capture.regions.file)))
    }

    ## Check for write permissions in the destination folder
    if (file.access(destination.folder, 2) == -1) {
        stop(.wrap("You do not have write permission in the destination",
                   "folder."))
    }

    ## Create lists with BAM files and index of corresponding control
    sample.paths <- unlist(sample.control)
    sample.paths <- unique(sample.paths[!is.na(sample.paths)])
    sample.files <- basename(sample.paths)
    control.indices <- match(sample.control$controls, sample.paths)
    sample.indices <- match(sample.control$samples, sample.paths)

    ## Load helper files
    blacklist.grange <- NULL
    GC.mappa.grange <- NULL
    load(file.path(reference.folder, "blacklist.rda"))
    load(file.path(reference.folder, "GC_mappability.rda"))

    ## Calculate the maximal number of CPUs to be used
    ncpu <- bpworkers(bp.param)

    ## Retrieve number of chromosomes and bin size from GC.mappa.grange object
    chromosomes <- seqlevels(GC.mappa.grange)
    bin.size <- width(GC.mappa.grange)[1]

    ## Create folders
    destination.folder <- file.path(destination.folder, "CNAprofiles")
    tryCatch({
        if (!file.exists(file.path(destination.folder))) {
            dir.create(file.path(destination.folder),
                       recursive = TRUE)
        } else {
            stop(.wrap("The folder",
                       sQuote(file.path(destination.folder, "CNAprofiles")),
                       "already exists. Please remove it, or (in case you",
                       "still need it), rename it to prevent files from being",
                       "overwritten."))
        }
    }, warning = function(e) {
        stop(.wrap("You do not have write permissions in the destination",
                   "folder. Stopping execution of the remaining part of the",
                   "script..."))
    })
    dir.create(file.path(destination.folder, "BamBaiPeaksFiles"),
               recursive = TRUE)

    ## Provide output to log
    flog.appender(appender.file(file.path(destination.folder,
                                          "CopywriteR.log")))
    flog.info(paste("Running CopywriteR version",
                    as(packageVersion("CopywriteR"), "character"), "..."))
    flog.info(paste0("CopywriteR was run using the following commands:", "\n\n",
                     "CopywriteR(sample.control = sample.control, ",
                     "destination.folder = \"", dirname(destination.folder),
                     "\", reference.folder = \"", reference.folder,
                     "\", BPPARAM = bp.param, capture.regions.file = \"",
                     capture.regions.file, "\", keep.intermediary.files = ",
                     keep.intermediary.files, ")"))
    flog.info("The value of bp.param was:", getClass(bp.param), capture = TRUE)
    flog.info("The value of sample.control was:", sample.control,
              capture = TRUE)
    flog.info(paste("The bin size for this analysis is", bin.size))
    flog.info(paste("The capture region file is", capture.regions.file))
    flog.info(paste("This analysis will be run on", ncpu, "cpus"))

    cat(.wrap("The following samples will be analyzed:"), "\n")
    cat(paste("sample:", sample.files[sample.indices], ";", "\t", "matching",
               "control:", sample.files[control.indices]), sep = "\n")
    cat(.wrap("The bin size for this analysis is", bin.size), "\n")
    cat(.wrap("The capture region file is", capture.regions.file), "\n")
    cat(.wrap("This analysis will be run on", ncpu, "cpus"), "\n")

    ## Test for compatibility chromosome names
    prefixes <- vector(mode = "character")
    chr.names <- NULL
    chr.lengths <- NULL
    chr.sort.mode <- NULL

    tryCatch({
        for (samp in sample.paths) {
            header <- scanBamHeader(samp)
            chr.sort.mode <- c(chr.sort.mode, header[[1]]$text$'@HD'[2])
            current.chr.names <- names(header[[1]]$targets)
            chr.names <- c(chr.names, current.chr.names)
            chr.lengths <- c(chr.lengths, header[[1]]$targets)
            prefixes <- append(prefixes, gsub("[[:digit:]]|X|Y|M|T",
                                              "", chr.names[1])[1])
        }
    }, error = function(e) {
        stop(.wrap("The BAM file header of file", sQuote(samp), "is corrupted",
                   "or truncated. Please rebuild this BAM file or exclude it",
                   "from analysis. Stopping execution of the remaining part of",
                   "the script..."))
    })

    if (!all(chr.sort.mode == "SO:coordinate")) {
        stop(.wrap("The following .bam files are unsorted:"), "\n",
             paste(sample.paths[which(!chr.sort.mode == "SO:coordinate")], "\n",
             "Please sort these .bam files based on coordinates",
             collapse = "\n"))
    }
    if (!all(prefixes == prefixes[1])) {
        stop(.wrap("The bam files have different chromosome name prefixes.",
                   "Please adjust the .bam files such that they contain the",
                   "same chromosome notation."))
    } else if (!length(unique(chr.names)) * length(sample.paths) == length(chr.names)) {
        stop(.wrap("The bam files have been mapped to different reference",
                   "genomes. Please run only .bam files mapped to the same",
                   "reference genome together."))
    } else if (!length(unique(chr.lengths)) * length(sample.paths) == length(chr.lengths)) {
        stop(.wrap("The bam files have different chromosome names. Please",
                   "adjust the .bam files such that they contain the same",
                   "chromosome notation."))
    } else {
        prefixes <- prefixes[1]

        chr.names <- seqlevels(GC.mappa.grange)
        prefixes[2] <- gsub("[[:digit:]]|X|Y", "", chr.names)[1]

        chr.names <- seqlevels(blacklist.grange)
        prefixes[3] <- gsub("[[:digit:]]|X|Y", "", chr.names)[1]

        if (!all(prefixes == prefixes[1])) {
            stop(.wrap("The bam files and supporting .bed files have different",
                       "chromosome names. Please adjust the input files such",
                       "that they contain the same chromosome notation."))
        }
    }

    ########################################################
    ## Calculate depth of coverage using off-target reads ##
    ########################################################

    ## Create list of .bam files
    flog.info("CopywriteR will analyze the following (unique) samples:",
              sample.files, capture = TRUE)

    ## Index .bam files
    if (!all(file.exists(gsub("$", ".bai", sample.paths)))) {
        if (file.access(".", 2) == -1) {
            stop(.wrap("The .bam files are not indexed and you do not have",
                       "write permission in (one of) the folder(s) where the",
                       ".bam files are located."))
        }
        IndexBam <- function(sample.paths) {
            indexBam(sample.paths)
            paste0("indexBam(\"", sample.paths, "\")")
        }
        to.log <- bplapply(sample.paths, IndexBam, BPPARAM = bp.param)
        lapply(to.log, flog.info)
    }

    ## Check whether BAMs are paired-end
    NumberPairedEndReads <- function(sample.paths) {
    
        bam <- open(BamFile(sample.paths, yieldSize = 1))
        close(bam)
        what <- c("flag")
        param <- ScanBamParam(what = what)
        bam <- readGAlignments(bam, param = param)
        intToBits(mcols(bam)$flag)[1] == 01
    }
    is.paired.end <- bplapply(sample.paths, NumberPairedEndReads,
                              BPPARAM = bp.param)
    is.paired.end <- unlist(is.paired.end)
    for (i in seq_along(sample.files)) {
        flog.info(paste0("Paired-end sequencing for sample ", sample.files[i],
                         ": ", is.paired.end[i]))
    }

    ## Remove anomalous reads and reads with Phred < 37
    i <- c(seq_along(sample.paths))
    ProperReads <- function(i, sample.paths, destination.folder, sample.files,
                            is.paired.end) {

        if (is.paired.end[i]) {
            flag <- scanBamFlag(isProperPair = TRUE)
            param <- ScanBamParam(flag = flag, what = "mapq")
            filter <- FilterRules(list(isHighQual = function(x) {
                x$mapq >= 37
            }))
            filterBam(sample.paths[i], file.path(destination.folder,
                                                 "BamBaiPeaksFiles",
                                                 gsub(".bam$",
                                                      "_properreads.bam",
                                                      sample.files[i])),
                      filter = filter, indexDestination = TRUE, param = param)
            paste0("filterBam(\"", sample.paths[i], "\", \"",
                   file.path(destination.folder, "BamBaiPeaksFiles",
                             gsub(".bam$", "_properreads.bam",
                                  sample.files[i])),
                   "\", filter = filter, indexDestination = TRUE, ",
                   "param = param)")
        } else {
            param <- ScanBamParam(what = "mapq")
            filter <- FilterRules(list(isHighQual = function(x) {
                x$mapq >= 37
            }))
            filterBam(sample.paths[i], file.path(destination.folder,
                                                 "BamBaiPeaksFiles",
                                                 gsub(".bam$",
                                                      "_properreads.bam",
                                                      sample.files[i])),
                      filter = filter, indexDestination = TRUE, param = param)
            paste0("filterBam(\"", sample.paths[i], "\", \"",
                   file.path(destination.folder, "BamBaiPeaksFiles",
                             gsub(".bam$", "_properreads.bam",
                                  sample.files[i])),
                   "\", filter = filter, indexDestination = TRUE, ",
                   "param = param)")
        }
    }
    to.log <- bplapply(i, ProperReads, sample.paths, destination.folder,
                       sample.files, is.paired.end, BPPARAM = bp.param)
    lapply(to.log, flog.info)

    ## Read count statistics
    Stats.1 <- function(sample.paths) {
        countBam(sample.paths)
    }
    res <- bplapply(sample.paths, Stats.1, BPPARAM = bp.param)
    res <- Reduce(function(x,y) {rbind(x,y)}, res)
    statistics <- res[, "records", drop = FALSE]
    rownames(statistics) <- sample.files
    statistics[, "total"] <- statistics$records
    statistics[, "records"] <- NULL

    ## Create new .bam list
    setwd(file.path(destination.folder, "BamBaiPeaksFiles"))
    sample.files <- gsub(".bam$", "_properreads.bam", sample.files)
    flog.info(paste("CopywriteR has finished filtering low-quality and",
                    "anomalous reads in the following (unique) samples:"),
                    sample.files, capture = TRUE)

    ## Read count statistics
    Stats.2 <- function(sample.files) {
        countBam(sample.files)$records
    }
    res <- bplapply(sample.files, Stats.2, BPPARAM = bp.param)
    res <- Reduce(function(x,y) {rbind(x,y)}, res)
    statistics[, "total.properreads"] <- res

    ## Create list with numbers of controls
    control.uniq.indices <- unique(control.indices)

    ## Call peaks in .bam file of control sample
    DetectPeaks <- function(control.uniq.indices, sample.files, prefix,
                            used.chromosomes, .peakCutoff, destination.folder) {

        ## j represents the minimal peak width
        j <- 100

        ## resolution is the resolution at which peaks are determined
        resolution <- 20000

        ## Initialize
        merged.bed <- NULL
        chromosomes <- scanBamHeader(sample.files[control.uniq.indices])[[1]][["targets"]]
        chromosomes <- chromosomes[names(chromosomes) %in% used.chromosomes]
        cov.all <- coverage(BamFile(sample.files[control.uniq.indices]))
        
        ## Obtain read length
        bam <- open(BamFile(sample.files[control.uniq.indices], yieldSize = 1))
        close(bam)
        bam <- readGAlignments(bam)
        read.length <- qwidth(bam)[1]
        
        for (selection in seq_along(chromosomes)) {

            ## Read bam file per chromosome and calculate coverage
            cov.chr <- cov.all[[names(chromosomes)[selection]]]
            cov.chr <- as.vector(cov.chr) # Makes subsetting faster

            ## Calculate peak detection and extension cutoffs per bin
            peak.detection.cutoff <- vector(length = length(cov.chr)%/%resolution)
            cov.chr.subsets <- NULL
            for (i in seq_len(length(cov.chr)%/%resolution)) {
                cov.chr.subsets[[i]] <- cov.chr[((i - 1) * resolution + 1):(i * resolution)]
                peak.detection.cutoff[i] <- ceiling(.peakCutoff(cov.chr.subsets[[i]]))
            }

            ## Fill in missing values (zeroes and NAs) in peak.detection.cutoffs
            peak.detection.cutoff[which(peak.detection.cutoff == 0)] <- NA
            no.values.cutoffs <- which(is.na(peak.detection.cutoff))
            no.values.replacements <- vector(length = length(no.values.cutoffs))
            for (i in seq_along(no.values.cutoffs)) {
                index <- no.values.cutoffs[i]
                while (is.na(peak.detection.cutoff[index]) & index > 1) {
                    index = index - 1
                }
                lower.value <- peak.detection.cutoff[index]
                index <- no.values.cutoffs[i]
                while (is.na(peak.detection.cutoff[index]) & index < length(peak.detection.cutoff)) {
                    index = index + 1
                }
                upper.value <- peak.detection.cutoff[index]
                no.values.replacements[i] <- ceiling(mean(c(lower.value,
                                                            upper.value),
                                                          na.rm = TRUE))
            }
            peak.detection.cutoff[no.values.cutoffs] <- no.values.replacements

            ## Create islands and concatenate; select islands that are bigger than certain width (j)
            # Use mapply for simultaneously iterating lists and vectors
            # If no reads are present anywhere on a chromosome,
            # peak.detection.cutoff is NaN for all bins -> test for this
            if (!is.na(peak.detection.cutoff[1])) {
            
                peak.ranges <- mapply(function(x, z) {
                    x <- slice(x, lower = peak.detection.cutoff[z])
                    shift(ranges(x), resolution * (z - 1))
                }, cov.chr.subsets, seq_along(peak.detection.cutoff))
                peak.ranges <- Reduce(function(x, y) {
                    c(x, y)
                }, peak.ranges)
                peak.ranges <- reduce(peak.ranges)
                peak.ranges <- peak.ranges[width(peak.ranges) > j, ]
                peak.ranges <- peak.ranges[end(peak.ranges) < chromosomes[selection]%/%resolution * resolution, ]
                
                ## Create RleViews object and calculate peakSummary
                peaks.ranges.rleviews <- Views(Rle(cov.chr), peak.ranges)
                peaks <- peakSummary(peaks.ranges.rleviews)
                
            } else {
                peaks <- data.frame()
            }

            if (nrow(peaks) > 0) {
                test <- data.frame(seqnames = names(chromosomes)[selection],
                                   start = start(peaks), end = end(peaks))

                ## Reiterate over peaks to check for for large differences in peak detection cutoffs
                retest.peak.ranges <- apply(test, 1, function(x) {
                    left.lower.boundary <- max(0, (as.integer(x["start"]) - (resolution + 1)))
                    left.higher.boundary <- max(0, (as.integer(x["start"]) - 1))
                    right.lower.boundary <- min(chromosomes[selection],
                                                (as.integer(x["end"]) + 1))
                    right.higher.boundary <- min(chromosomes[selection],
                                                 (as.integer(x["end"]) + (resolution + 1)))
                    left.peakCutoff <- ceiling(.peakCutoff(cov.chr[left.lower.boundary:left.higher.boundary]))
                    right.peakCutoff <- ceiling(.peakCutoff(cov.chr[right.lower.boundary:right.higher.boundary]))
                    max.peakCutoff <- max(left.peakCutoff, right.peakCutoff)
                    tmp <- slice(cov.chr[as.integer(x["start"]):as.integer(x["end"])],
                                 lower = max.peakCutoff)
                    shift(ranges(tmp), as.integer(x["start"]) - 1)
                })

                retest.peak.ranges <- Reduce(function(x, y) {
                    c(x, y)
                }, retest.peak.ranges)
                retest.peak.ranges <- reduce(retest.peak.ranges)

                ## Select all with width of more than 100 and within bin regions
                retest.peak.ranges <- retest.peak.ranges[width(retest.peak.ranges) > j, ]
                retest.peak.ranges <- retest.peak.ranges[end(retest.peak.ranges) < chromosomes[selection]%/%resolution * resolution, ]

                ## Create RleViews object and calculate peakSummary
                retest.peaks.ranges.rleviews <- Views(S4Vectors::Rle(cov.chr),
                                                      retest.peak.ranges)
                retest.peaks <- peakSummary(retest.peaks.ranges.rleviews)

                test <- cbind(start(retest.peaks), end(retest.peaks))
                colnames(test) <- c("start", "end")

                if (nrow(test) > 0) {
                    ## Calculate extension cutoff for every peak
                    lower.cutoff.peaks <- apply(test, 1, function(x) {
                        left.lower.boundary <- max(0, (as.integer(x["start"]) - (resolution + 1)))
                        left.higher.boundary <- max(0, (as.integer(x["start"]) - 1))
                        right.lower.boundary <- min(chromosomes[selection],
                                                    (as.integer(x["end"]) + 1))
                        right.higher.boundary <- min(chromosomes[selection],
                                                     (as.integer(x["end"]) + (resolution + 1)))
                        left.peakCutoff <- floor(.peakCutoff(cov.chr[left.lower.boundary:left.higher.boundary],
                                                 fdr.cutoff = 0.1))
                        right.peakCutoff <- floor(.peakCutoff(cov.chr[right.lower.boundary:right.higher.boundary],
                                                  fdr.cutoff = 0.1))
                        min(left.peakCutoff, right.peakCutoff)
                    })
                    lower.cutoff.peaks <- unlist(lower.cutoff.peaks)

                    if (nrow(test) > 0) {
                        for (i in seq_len(nrow(test))) {
                            index <- test[i, "start"]
                            while (cov.chr[index] > lower.cutoff.peaks[i]) {
                                index = index - 1
                            }
                            test[i, "start"] <- index - read.length
                            index <- test[i, "end"]
                            while (cov.chr[index] > lower.cutoff.peaks[i]) {
                                index = index + 1
                            }
                            test[i, "end"] <- index + read.length
                        }
                        test <- as(data.frame(seqnames = names(chromosomes)[selection],
                                              test), "GRanges")
                        test <- test[order(test)]
                        test <- reduce(test)
                        test <- as(test, "data.frame")
                        test[, "width"] <- NULL
                        test[, "strand"] <- NULL
                        merged.bed <- rbind(merged.bed, test)
                    }
                }
            }
        }

        ## Write data
        write.table(merged.bed, file = file.path(destination.folder,
                                                 "BamBaiPeaksFiles",
                                                 paste0("peaks",
                                                        control.uniq.indices,
                                                        ".bed")), sep = "\t", 
                    row.names = FALSE, col.names = FALSE, quote = FALSE)

        paste0("Analysis of peaks in sample ",
               sample.files[control.uniq.indices],
               " is done; output file: peaks", control.uniq.indices, ".bed")
    }
    # current.bp.param <- bp.param
    # bp.param$workers <- ifelse(bp.param$workers < length(control.uniq.indices),
    #                           bp.param$workers, length(control.uniq.indices))
    to.log <- bplapply(control.uniq.indices, DetectPeaks, sample.files,
                       prefixes[1], chromosomes, .peakCutoff,
                       destination.folder, BPPARAM = bp.param)
    # bp.param <- current.bp.param
    lapply(to.log, flog.info)

    ## Read count statistics
    Stats.3 <- function(sample.files, GC.mappa.grange) {
        all.reads <- countBam(sample.files)$records
        which <- reduce(GC.mappa.grange)
        what <- c("pos")
        param <- ScanBamParam(which = which, what = what)
        chrom.reads <- countBam(file = sample.files, param = param)
        chrom.reads <- sum(chrom.reads$records)
        c(all.reads = all.reads, chrom.reads = chrom.reads)
    }
    res <- bplapply(sample.files, Stats.3, GC.mappa.grange, BPPARAM = bp.param)
    res <- res <- data.frame(do.call(Map, c(rbind, res)))
    statistics[, "unmappable.or.mitochondrial"] <- res$all.reads - res$chrom.reads
    statistics[, "on.chromosomes"] <- res$chrom.reads

    ## Calculate coverage
    i <- c(seq_along(sample.indices))
    CalculateDepthOfCoverage <- function(i, sample.files, control.indices,
                                         sample.indices, GC.mappa.grange,
                                         bin.size) {

        # Create GRanges object of peak files
        if (file.info(paste0("peaks", control.indices[i], ".bed"))$size != 0) {
            bed <- read.table(file = paste0("peaks", control.indices[i], ".bed"),
                              as.is = TRUE, sep = "\t")
            colnames(bed) <- c("Seqnames", "Start", "End")
            peak.grange <- makeGRangesFromDataFrame(bed)

            # Calculate setdiff without reducing ranges
            outside.peak.grange <- split(GC.mappa.grange, rep_len(c(1, 2),
                                         length.out = length(GC.mappa.grange)))
            outside.peak.grange <- lapply(outside.peak.grange, function(x) {
                setdiff(x, peak.grange)
            })
            outside.peak.grange <- c(outside.peak.grange[[1]],
                                     outside.peak.grange[[2]])
            outside.peak.grange <- outside.peak.grange[order(outside.peak.grange)]
        } else {
            outside.peak.grange <- GC.mappa.grange
        }

        # countBam on remainder of bins
        param <- ScanBamParam(which = outside.peak.grange, what = c("pos"))
        counts <- countBam(sample.files[sample.indices[i]], param = param)

        # countBam on remainder of bins
        param <- ScanBamParam(which = reduce(outside.peak.grange),
                              what = c("pos"))
        counts.CopywriteR <- countBam(sample.files[sample.indices[i]],
                                      param = param)
        counts.CopywriteR <- sum(counts.CopywriteR$records)

        # Fix MT as levels in factor seqnames (remainder from setdiff operation)
        counts$space <- as.factor(as.character(counts$space))
        colnames(counts)[1] <- "seqnames"

        # Replace bins by real bins & calculate compensated read counts
        sample.files[sample.indices[i]] <- gsub("_properreads.bam$", ".bam",
                                                sample.files[sample.indices[i]])

        # Aggregate counts and total length per bin
        counts.grange <- makeGRangesFromDataFrame(counts[, c("seqnames",
                                                             "start", "end",
                                                             "records")],
                                                  keep.extra.columns = TRUE)
        overlaps <- findOverlaps(counts.grange, GC.mappa.grange, minoverlap = 1L)
        index <- subjectHits(overlaps)
        records <- mcols(counts.grange[queryHits(overlaps)])$records
        lengths <- width(pintersect(counts.grange[queryHits(overlaps)],
                                    GC.mappa.grange[subjectHits(overlaps)]))
        aggregate.data.table <- data.table(index, records, lengths)
        aggregate.data.table <- aggregate.data.table[, list(records = sum(records),
                                                            lengths = sum(lengths)),
                                                     by = c("index")]
        aggregate.data.table <- aggregate.data.table[match(seq_along(GC.mappa.grange),
                                                           aggregate.data.table$index), ]
        aggregate.data.table$records[which(is.na(aggregate.data.table$index))] <- 0
        aggregate.data.table$lengths[which(is.na(aggregate.data.table$index))] <- 0

        # Replace bins by real bins & calculate compensated read counts
        counts <- aggregate.data.table
        counts[, "Chromosome"] <- as(seqnames(GC.mappa.grange), "factor")
        counts[, "Start"] <- start(GC.mappa.grange)
        counts[, "End"] <- end(GC.mappa.grange)
        counts[, "Feature"] <- paste0(counts$Chromosome, ":",
                                      paste0(counts$Start, "-", counts$End))
        counts[, paste0("read.counts.", sample.files[sample.indices[i]])] <-
            counts$records
        counts[, paste0("read.counts.compensated.",
                        sample.files[sample.indices[i]])] <- 
            counts$records / (counts$lengths / bin.size)
        counts[, paste0("fraction.of.bin.", sample.files[sample.indices[i]])] <-
            counts$lengths / bin.size
        counts$index <- NULL
        counts$lengths <- NULL
        counts$records <- NULL
        counts <- data.frame(counts, check.names = FALSE)

        # Return
        return(list(counts[, paste0("read.counts.compensated.",
                                    sample.files[sample.indices[i]]),
                                    drop = FALSE],
                    counts[, paste0("read.counts.",
                                    sample.files[sample.indices[i]]), 
                                    drop = FALSE],
                    counts[, paste0("fraction.of.bin.",
                                    sample.files[sample.indices[i]]),
                                    drop = FALSE],
                    paste0("Rsamtools finished calculating reads per bin in ", 
                           "sample ", sample.files[sample.indices[i]],
                           "; number of bins = ", nrow(counts)),
                    counts.CopywriteR))
    }
    res <- bplapply(i, CalculateDepthOfCoverage, sample.files, control.indices,
                    sample.indices, GC.mappa.grange, bin.size,
                    BPPARAM = bp.param)
    read.counts <- data.frame(Chromosome = as(seqnames(GC.mappa.grange),
                                              "factor"),
                              Start = start(GC.mappa.grange),
                              End = end(GC.mappa.grange))
    read.counts[, "Feature"] <- paste0(read.counts$Chromosome, ":",
                                       paste0(read.counts$Start, "-",
                                              read.counts$End))
    ## ‘Map’ applies a function to the corresponding elements of given vectors.
    res <- do.call(Map, c(cbind, res))
    read.counts <- cbind(read.counts[, ], Reduce(cbind, res[1:3]))

    ## Remove potential NAs introduced by peaks spanning entire bins
    read.counts[is.na(read.counts)] <- 0
    
    sapply(res[[4]], flog.info)

    ## Read count statistics
    statistics[, "off.target"] <- unlist(res[5])
    statistics[, "on.target"] <- statistics$on.chromosomes - statistics$off.target
    statistics$on.chromosomes <- NULL
    flog.info(paste("The following number of sequence reads were found in the",
                    "samples analyzed by CopywriteR:"), statistics,
              capture = TRUE)

    write.table(read.counts[mixedorder(read.counts$Chromosome), ],
                file = file.path(destination.folder, "read_counts.txt"),
                row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")

    ## Create histograms of fraction.of.bin
    ## (fraction of length in bins (after peak region removal)
    sample.files <- gsub("_properreads.bam$", ".bam", sample.files)

    dir.create(file.path(destination.folder, "qc"))
    for (i in seq_along(sample.indices)) {
        pdf(file.path(destination.folder, "qc",
                      paste0("fraction.of.bin.",
                             sample.files[sample.indices[i]], ".pdf")),
            width = 7, height = 7)
        plot(ecdf(as.numeric(read.counts[, 4 + (2 * length(sample.indices)) + i])),
             verticals = TRUE, ylab = "Fraction of bins",
             xlab = "Remaining fraction of bin after peak removal", 
             main = "Cumulative distribution of remaining bin fraction")
        dev.off()
    }

    ##############################################
    ## Normalize for GC-content and mappability ##
    ##############################################

    ## Create data input file for correction function
    read.counts <- read.counts[, seq_len(4 + length(sample.indices))]
    data <- list(cov = read.counts[, 5:ncol(read.counts), drop = FALSE],
                 anno = read.counts[, c("Chromosome", "Start",
                                        "End", "Feature")])
    data$anno[, "mappa"] <- GC.mappa.grange$mappability
    data$anno[, "gc"] <- GC.mappa.grange$GCcontent
    data$anno[, "black"] <- overlapsAny(GC.mappa.grange, blacklist.grange)

    usepoints <- !(data$anno$Chromosome %in% c("X", "Y", "MT",
                                               "chrX", "chrY", "chrM"))

    ## Perform GC-content and mappability corrections (in .tng helper function)
    tryCatch({
        i <- c(seq_len(ncol(data$cov)))
        NormalizeDOC <- function(i, data, .tng, usepoints, destination.folder) {
            .tng(data.frame(count = data$cov[, i], gc = data$anno$gc,
                            mappa = data$anno$mappa),
                 use = usepoints & data$cov[, i] != 0, correctmappa = TRUE,
                 plot = file.path(destination.folder, "qc",
                                  paste0(colnames(data$cov)[i], ".png")))
        }
        ratios <- bplapply(i, NormalizeDOC, data, .tng, usepoints,
                           destination.folder, BPPARAM = bp.param)
        log2.read.counts <- matrix(unlist(ratios), ncol = length(sample.indices))
    }, error = function(e) {
        stop(.wrap("The GC-content and mappability normalization did not work",
                   "due to a failure to calculate loesses. This can generally",
                   "be solved by using larger bin sizes. Stopping execution of",
                   "the remaining part of the script..."))
    })

    colnames(log2.read.counts) <- paste0("log2.", sample.files[sample.indices])

    ###################
    ## Create output ##
    ###################

    ## Create table with corrected log2 values and write to file
    selection <- !data$anno$black & !is.na(data$anno$mappa) & data$anno$mappa > 0.2
    log2.read.counts <- data.frame(data$anno[selection, c("Chromosome", "Start",
                                                          "End", "Feature")],
                                   log2.read.counts[selection, , drop = FALSE],
                                   check.names = FALSE)
    log2.read.counts$Chromosome <- as.vector(log2.read.counts$Chromosome)

    ## Replace -/+Inf values to -/+large values for compatibility with IGV browser
    log2.read.counts[log2.read.counts == -Inf] <- -.Machine$integer.max / 2
    log2.read.counts[log2.read.counts == Inf] <- .Machine$integer.max / 2

    ## Create tracking line for viewing options in IGV
    igv.track.line <- "#track viewLimits=-3:3 graphType=heatmap color=255,0,0"
    
    ## Write output
    write.table(igv.track.line, file = file.path(destination.folder,
                                                 "log2_read_counts.igv"),
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    ## Suppress warning that you are 'appending column names to file'
    suppressWarnings(write.table(log2.read.counts[mixedorder(log2.read.counts$Chromosome), ],
                                 file = file.path(destination.folder,
                                                  "log2_read_counts.igv"),
                                 append = TRUE, sep = "\t",
                                 row.names = FALSE, quote = FALSE))

    #############################################################################
    ## Calculate overlap with capture.regions.file for quality control purpose ##
    #############################################################################

    ## Calculate overlap with capture regions
    if (capture.regions.file != "not specified") {
        captured.bed <- read.table(capture.regions.file, as.is = TRUE,
                                   sep = "\t")
        colnames(captured.bed) <- c("seqnames", "Start", "End")
        captured.grange <- makeGRangesFromDataFrame(captured.bed)
				seqlevels(captured.grange) <- gsub("^", unique(prefixes),
																	         gsub(unique(prefixes), "",
																				        seqlevels(captured.grange)))
        
        for (control.index in control.uniq.indices) {
            peak.bed <- read.table(file = paste0("peaks", control.index,
                                                 ".bed"),
                                   as.is = TRUE, sep = "\t")
            colnames(peak.bed) <- c("seqnames", "Start", "End")
            peak.grange <- makeGRangesFromDataFrame(peak.bed)

            overlap <- findOverlaps(captured.grange, peak.grange)

            flog.info(paste0("Number of capture regions covered by peaks in ",
                             "sample ", sample.files[control.index], ": ",
                             length(unique(queryHits(overlap)))))
            flog.info(paste0("Number of peaks covered by capture regions in ",
                             "sample ", sample.files[control.index], ": ",
                             length(unique(subjectHits(overlap)))))
            flog.info(paste0("Total number of capture regions in sample ",
                             sample.files[control.index], ": ",
                             length(captured.grange)))
            flog.info(paste0("Total number of peaks in sample ",
                             sample.files[control.index], ": ",
                             length(peak.grange)))
        }
    }

    ## Remove BamBaiPeaksFiles folder
    if (!keep.intermediary.files) {
        flog.info(paste("The 'keep.intermediary.files variable' was set to",
                        "FALSE; removing temporary files ..."))
        unlink(file.path(destination.folder, "BamBaiPeaksFiles"),
               recursive = TRUE)
    }

    flog.info(paste("Total calculation time of CopywriteR was",
                    round(difftime(Sys.time(), start.time, units = "hours"), 2),
                    "hours"))
    cat("Total calculation time of CopywriteR was: ",
        Sys.time() - start.time, "\n\n")

    inputStructure <- list(sample.control = sample.control,
                           chromosomes = chromosomes, prefix = prefixes[1],
                           bin.size = bin.size)
    save(inputStructure, file = file.path(destination.folder, "input.Rdata"))

}
