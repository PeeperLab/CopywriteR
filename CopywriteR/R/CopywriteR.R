CopywriteR <- function(sample.control, destination.folder, reference.folder,
                       ncpu, bpparam, capture.regions.file,
                       keep.intermediairy.files = FALSE) {

    ##########################
    ## Check and initialise ##
    ##########################

    start.time <- Sys.time()
    
    ## Set old options back on exit
    old.options <- options()
    on.exit(options(old.options))

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
            stop("The file ", x, " could not be found.\nPlease change the ",
                 "path to this file.")
        }
    }))

    if (!file.exists(destination.folder)) {
        stop("The destination.folder could not be found.\nPlease change your ",
             "destination.folder path.")
    }

    if (!file.exists(reference.folder)) {
        stop("The reference.folder could not be found.\nPlease change your ",
             "reference.folder path or run `preCopywriteR` to generate the ",
             "required folder with GC-content and mapability files for your ",
             "desired bin size.")
    }

    if (!file.exists(capture.regions.file) & capture.regions.file != "not specified") {
        stop("The capture.regions.file could not be found.\nPlease change ",
             "your capture.regions.file path.")
    }

    ## Check for write permissions in the destination folder
    if (file.access(destination.folder, 2) == -1) {
        stop("You do not have write permission in the destination folder.")
    }

    ## Create lists with BAM files and index of corresponding control
    sample.paths <- unlist(sample.control)
    sample.paths <- unique(sample.paths[!is.na(sample.paths)])
    sample.files <- basename(sample.paths)
    control.indices <- match(sample.control$controls, sample.control$samples)

    ## Create paths to helper files
    bin.file <- file.path(reference.folder, "bins.bed")
    blacklist.file <- file.path(reference.folder, "blacklist.bed")
    gc.content.file <- file.path(reference.folder, "GC_content.bed")
    mapability.file <- file.path(reference.folder, "mapability.bed")

    ## Cap the number of cpus to be used to the number of samples
    ncpu <- if (ncpu < nrow(sample.control)) ncpu else nrow(sample.control)

    ## Retrieve number of chromosomes and bin size from bin.bed helper file
    bin.bed <- read.table(file = bin.file, as.is = TRUE, sep = "\t")
    colnames(bin.bed) <- c("Chromosome", "Start", "End")
    chromosomes <- unique(bin.bed$Chromosome)
    bin.size <- bin.bed$End[1]

    ## Create folders
    destination.folder <- file.path(destination.folder, "CNAprofiles")
    tryCatch({
				if (!file.exists(file.path(destination.folder))) {
						dir.create(file.path(destination.folder),
											 recursive = TRUE)
				} else {
						stop("The folder ", file.path(destination.folder, "CNAprofiles"),
								 " already exists. Please remove it, or (in case you still need ",
								 "it), rename it to prevent files from being overwritten.")
				}
    }, warning = function(e) {
        stop("You do not have write permissions in the destination folder.\n",
             "Stopping execution of the remaining part of the script...")
    })
		dir.create(file.path(destination.folder, "BamBaiPeaksFiles"),
							 recursive = TRUE)

    ## Provide output for log file
    sink(file = file.path(destination.folder, "log.txt"),
         type = c("output", "message"))
    options(width = 150)

    cat("The following samples will be analyzed:\n")
    cat(paste0("sample: ", sample.files, ";\tmatching control: ",
               sample.files[control.indices]), sep = "\n")
    cat("The bin size for this analysis is", bin.size, "\n")
    cat("The capture region file is", capture.regions.file, "\n")
    cat("This analysis will be run on", ncpu, "cpus", "\n\n\n")
    sink()

    cat("The following samples will be analyzed:\n")
    cat(paste0("sample: ", sample.files, ";\tmatching control: ",
               sample.files[control.indices]), sep = "\n")
    cat("The bin size for this analysis is", bin.size, "\n")
    cat("The capture region file is", capture.regions.file, "\n")
    cat("This analysis will be run on", ncpu, "cpus", "\n\n\n")

    ## Test for compatibilty chromosome names
    prefixes <- vector(mode = "character")
    chr.names <- NULL
    chr.lengths <- NULL

    tryCatch({
        for (samp in sample.paths) {
            header <- scanBamHeader(samp)
            current.chr.names <- names(header[[1]]$targets)
            chr.names <- c(chr.names, current.chr.names)
            chr.lengths <- c(chr.lengths, header[[1]]$targets)
            prefixes <- append(prefixes, gsub("[[:digit:]]|X|Y|M|T",
                                              "", chr.names[1])[1])
        }
    }, error = function(e) {
        stop("The BAM file header of file", samp, "is corrupted or ",
             "truncated.\nPlease rebuild this BAM file or exclude it from ",
             "analysis.\nStopping execution of the remaining part of the ",
             "script...")
    })

    if (!all(prefixes == prefixes[1])) {
        stop("The bam files have different chromosome name prefixes.\nPlease ",
             "adjust the .bam files such that they contain the same ",
             "chromosome notation.")
    } else if (!length(unique(chr.names)) * length(sample.paths) == length(chr.names)) {
        stop("The bam files have been mapped to different reference genomes.\n",
             "Please run only .bam files mapped to the same reference genome ", 
             "together.")
    } else if (!length(unique(chr.lengths)) * length(sample.paths) == length(chr.lengths)) {
        stop("The bam files have different chromosome names.\nPlease adjust ",
             "the .bam files such that they contain the same chromosome ",
             "notation.")
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
            stop("The bam files and supporting .bed files have different ",
                 "chromosome names.\nPlease adjust the input files such that ",
                 "they contain the same chromosome notation.")
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
    if (!any(list.files(pattern = ".bam$") == gsub(".bai$", "", list.files(pattern = ".bai$")))) {
        if (file.access(".", 2) == -1) {
            stop("The .bam files are not indexed and you do not have write ",
                 "permission in (one of) the folder(s) where the .bam files ",
                 "are located.")
        }
        IndexBam <- function(sample.paths) {
            library(Rsamtools)
            indexBam(sample.paths)
            paste0("indexBam(\"", sample.paths, "\")")
        }
        to.log <- bplapply(sample.paths, IndexBam, BPPARAM = bpparam)
        cat(unlist(to.log), "\n", sep = "\n")

        ## Garbage collection
        rm(IndexBam)
    }

    ## Check whether BAMs are paired-end
    NumberPairedEndReads <- function(sample.paths) {
        library(Rsamtools)
        flag <- scanBamFlag(isPaired = TRUE)
        param <- ScanBamParam(flag = flag)
        countBam(sample.paths, param = param)
    }
    is.paired.end <- bplapply(sample.paths, NumberPairedEndReads,
                              BPPARAM = bpparam)
    is.paired.end <- Reduce(function(x, y) {
        merge(x, y, all = TRUE)
    }, is.paired.end)
    is.paired.end <- ifelse(is.paired.end$records > 0, TRUE, FALSE)
    for (i in 1:length(sample.files)) {
        cat("Paired-end sequencing for sample ", sample.files[i], ": ",
            is.paired.end[i], "\n", sep = "")
    }
    cat("\n\n")

    ## Remove anomalous reads and reads with Phred < 37
    i <- c(1:length(sample.paths))
    ProperReads <- function(i, sample.paths, destination.folder, sample.files,
                            is.paired.end) {
        library(Rsamtools)
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
                             gsub(".bam$", "_properreads.bam", sample.files[i])),
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
                             gsub(".bam$", "_properreads.bam", sample.files[i])),
                             "\", filter = filter, indexDestination = TRUE, ",
                             "param = param)")
        }
    }
    to.log <- bplapply(i, ProperReads, sample.paths, destination.folder,
                       sample.files, is.paired.end, BPPARAM = bpparam)
    cat(unlist(to.log), "\n", sep = "\n")

    ## Read count statistics
    Stats <- function(sample.paths) {
        library(Rsamtools)
        countBam(sample.paths)
    }
    res <- bplapply(sample.paths, Stats, BPPARAM = bpparam)
    res <- Reduce(function(x,y) {rbind(x,y)}, res)
    statistics <- res[, "records", drop = FALSE]
    rownames(statistics) <- sample.files
    statistics <- within(statistics, {
        total <- records
        rm(records)
    })

    ## Create new .bam list
    setwd(file.path(destination.folder, "BamBaiPeaksFiles"))
    sample.files <- gsub(".bam$", "_properreads.bam", sample.files)
    cat(sample.files, "\n", sep = "\n")

    ## Read count statistics
    Stats <- function(sample.files) {
        library(Rsamtools)
        countBam(sample.files)$records
    }
    res <- bplapply(sample.files, Stats, BPPARAM = bpparam)
    res <- Reduce(function(x,y) {rbind(x,y)}, res)
    statistics <- within(statistics, {
        total.properreads <- res
    })

    ## Garbage collection
    rm(is.paired.end, NumberPairedEndReads, ProperReads, sample.paths, Stats)

    ## Create list with numbers of controls
    control.uniq.indices <- unique(control.indices)

    ## Call peaks in .bam file of control sample
    DetectPeaks <- function(control.uniq.indices, sample.files, prefix,
                            used.chromosomes, .peakCutoff, destination.folder) {

        library(Rsamtools)
        library(chipseq)
        library(GenomicRanges)
        library(GenomicAlignments)

        ## j represents the minimal peak width
        j <- 100

        ## resolution is the resolution at which peaks are determined
        resolution <- 20000

        ## Initialize
        merged.bed <- NULL
        chromosomes <- scanBamHeader(sample.files[control.uniq.indices])[[1]][["targets"]]
        chromosomes <- chromosomes[names(chromosomes) %in% used.chromosomes]

        for (selection in 1:length(chromosomes)) {

            ## Read bam file per chromosome and calculate coverage
            which <- GRanges(names(chromosomes)[selection],
                             IRanges(1, unname(chromosomes)[selection]))
            what <- c("rname", "pos", "strand", "qwidth")
            param <- ScanBamParam(which = which, what = what)
            bam <- readGAlignments(sample.files[control.uniq.indices],
                                   param = param)
            cov.chr <- coverage(bam)
            cov.chr <- cov.chr@listData[names(chromosomes)[selection]][[1]]
            cov.chr <- as.vector(cov.chr) # Makes subsetting faster

            ## Calculate peak detection and extension cutoffs per bin
            peak.detection.cutoff <- vector(length = length(cov.chr)%/%resolution)
            cov.chr.subsets <- NULL
            for (i in 1:(length(cov.chr)%/%resolution)) {
                cov.chr.subsets[[i]] <- cov.chr[(((i - 1) * resolution)):(i * resolution)]
                peak.detection.cutoff[i] <- ceiling(.peakCutoff(cov.chr.subsets[[i]]))
            }

            ## Fill in missing values (zeroes and NAs) in peak.detection.cutoffs
            peak.detection.cutoff[which(peak.detection.cutoff == 0)] <- NA
            no.values.cutoffs <- which(is.na(peak.detection.cutoff))
            no.values.replacements <- vector(length = length(no.values.cutoffs))
            for (i in 1:length(no.values.cutoffs)) {
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
										print(x)
										x <- slice(x, lower = peak.detection.cutoff[z])
										shift(x@ranges, resolution * (z - 1))
								}, cov.chr.subsets, 1:length(peak.detection.cutoff))
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
                    shift(tmp@ranges, as.integer(x["start"]) - 1)
                })

                retest.peak.ranges <- Reduce(function(x, y) {
                    c(x, y)
                }, retest.peak.ranges)
                retest.peak.ranges <- reduce(retest.peak.ranges)

                ## Select all with width of more than 100 and within bin regions
                retest.peak.ranges <- retest.peak.ranges[width(retest.peak.ranges) > j, ]
                retest.peak.ranges <- retest.peak.ranges[end(retest.peak.ranges) < chromosomes[selection]%/%resolution * resolution, ]

                ## Create RleViews object and calculate peakSummary
                retest.peaks.ranges.rleviews <- Views(Rle(cov.chr),
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

                    read.length <- qwidth(bam)[1]
                    if (nrow(test) > 0) {
                        for (i in 1:nrow(test)) {
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
                        test <- within(test, {
                            rm(width, strand)
                        })
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
    to.log <- bplapply(control.uniq.indices, DetectPeaks, sample.files,
                       prefixes[1], unique(bin.bed$Chromosome), .peakCutoff,
                       destination.folder, BPPARAM = bpparam)
    cat(unlist(to.log), "\n", sep = "\n")

    ## Read count statistics
    Stats <- function(sample.files, bin.bed) {
        library(Rsamtools)
        all.reads <- countBam(sample.files)$records
        which <- with(bin.bed, reduce(GRanges(seqnames = Chromosome,
                                              ranges = IRanges(start = Start,
                                                               end = End))))
        what <- c("pos")
        param <- ScanBamParam(which = which, what = what)
        chrom.reads <- countBam(file = sample.files, param = param)
        chrom.reads <- sum(chrom.reads$records)
        c(all.reads = all.reads, chrom.reads = chrom.reads)
    }
    res <- bplapply(sample.files, Stats, bin.bed, BPPARAM = bpparam)
    res <- data.frame(Reduce(function(x,y) {rbind(x,y)}, res))
    statistics <- within(statistics, {
        total.properreads <- res$all.reads
        unmapable.or.mitochondrial <- res$all.reads - res$chrom.reads
        on.chromosomes <- res$chrom.reads
    })

    ## Alternative for bedtools
    # Create GRanges file for bins
    bin.grange <- with(bin.bed, GRanges(seqnames = Chromosome,
                                        ranges = IRanges(start = Start,
                                                         end = End)))
    bin.grange <- bin.grange[order(bin.grange)]
                                                  

    i <- c(1:length(sample.files))
    CalculateDepthOfCoverage <- function(i, sample.files, control.indices,
                                         bin.grange, bin.size) {

        library(Rsamtools)
        library(data.table)
        
        # Create GRanges object of peak files
        if (file.info(paste0("peaks", control.indices[i], ".bed"))$size != 0) {
            bed <- read.table(file = paste0("peaks", control.indices[i], ".bed"),
                              as.is = TRUE, sep = "\t")
            colnames(bed) <- c("Chromosome", "Start", "End")
            peak.grange <- with(bed, GRanges(seqnames = Chromosome,
                                             ranges = IRanges(start = Start,
                                                              end = End)))

            # Calculate setdiff without reducing ranges
            outside.peak.grange <- split(bin.grange, rep_len(c(1, 2),
                                         length.out = length(bin.grange)))
            outside.peak.grange <- lapply(outside.peak.grange, function(x) {
                setdiff(x, peak.grange)
            })
            outside.peak.grange <- c(outside.peak.grange[[1]],
                                     outside.peak.grange[[2]])
            outside.peak.grange <- outside.peak.grange[order(outside.peak.grange)]
            
            # Fix the 1-based coordinate system
            ranges(outside.peak.grange)@width[1] <- (ranges(outside.peak.grange)@width[1] + 1L)
            ranges(outside.peak.grange)@start[1] <- 0L
        } else {
            outside.peak.grange <- bin.grange
        }

        # countBam on remainder of bins
        param <- ScanBamParam(which = outside.peak.grange, what = c("pos"))
        counts <- countBam(sample.files[i], param = param)

        # countBam on remainder of bins
        param <- ScanBamParam(which = reduce(outside.peak.grange),
                              what = c("pos"))
        counts.CopywriteR <- countBam(sample.files[i], param = param)
        counts.CopywriteR <- sum(counts.CopywriteR$records)

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
                   records/(lengths/(bin.size + 1)))
            assign(paste0("fraction.of.bin.", sample.files[i]),
                   lengths/(bin.size + 1))
            rm(index, lengths, records)
        })
        counts <- data.frame(counts, check.names = FALSE)

        # Return
        return(list(counts[, paste0("read.counts.compensated.",
                                    sample.files[i]), drop = FALSE],
                    counts[, paste0("read.counts.", sample.files[i]), 
                                    drop = FALSE],
                    counts[, paste0("fraction.of.bin.", sample.files[i]),
                                    drop = FALSE],
                    paste0("Rsamtools finished calculating reads per bin in ", 
                           "sample ", sample.files[i], "; number of bins = ",
                           nrow(counts)),
                    counts.CopywriteR))
    }
    res <- bplapply(i, CalculateDepthOfCoverage, sample.files, control.indices,
                    bin.grange, bin.size, BPPARAM = bpparam)
    read.counts <- data.frame(Chromosome = as(seqnames(bin.grange), "factor"),
                              Start = ranges(bin.grange)@start,
                              End = ranges(bin.grange)@start + ranges(bin.grange)@width - 1L)
    read.counts <- within(read.counts, {
        Feature = paste0(Chromosome, ":", paste0(Start, "-", End))
    })
    ## ‘Map’ applies a function to the corresponding elements of given vectors.
    res <- do.call(Map, c(cbind, res))
    read.counts <- cbind(read.counts[, ], Reduce(cbind, res[1:3]))

    ## Remove potential NAs introduced by peaks spanning entire bins
    res[is.na(res)] <- 0
    
    cat(unlist(res[4]), "\n", sep = "\n")

    ## Read count statistics
    statistics <- within(statistics, {
        off.target <- unlist(res[5])
        on.target <- on.chromosomes - off.target
        rm(on.chromosomes)
    })
    print(statistics)
    cat("\n\n")

    write.table(read.counts[mixedorder(read.counts$Chromosome), ],
                file = file.path(destination.folder, "read_counts.txt"),
                row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")

    ## Create histograms of fraction.of.bin
    ## (fraction of length in bins (after peak region removal)
    sample.files <- gsub("_properreads.bam$", ".bam", sample.files)

    dir.create(file.path(destination.folder, "qc"))
    for (i in 1:length(sample.files)) {
        pdf(file.path(destination.folder, "qc",
                      paste0("fraction.of.bin.", sample.files[i], ".pdf")),
            width = 7, height = 7)
        plot(ecdf(as.numeric(read.counts[, 4 + (2 * length(sample.files)) + i])),
             verticals = TRUE, ylab = "Fraction of bins",
             xlab = "Remaining fraction of bin after peak removal", 
             main = "Cumulative distribution of remaining bin fraction")
        dev.off()
    }

    ## Garbage collection
    rm(bin.bed, bin.grange, CalculateDepthOfCoverage, control.indices, res,
       statistics, to.log)

    #############################################
    ## Normalize for GC-content and mapability ##
    #############################################

    ## Read files
    read.counts <- read.counts[, 1:(4 + length(sample.files))]
    for (i in 5:ncol(read.counts)) {
        write.table(read.counts[, c(1, 2, 3, i)], colnames(read.counts)[i],
                    quote = FALSE, sep = "\t", row.names = FALSE,
                    col.names = FALSE)
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
    usepoints <- !(data$anno$chr %in% c("X", "Y", "MT", "chrX", "chrY", "chrM"))

    ## Perform normalization (in .tng helper function)
    tryCatch({
        i <- c(1:ncol(data$cov))
        NormalizeDOC <- function(i, data, .tng, usepoints, destination.folder) {
            .tng(data.frame(count = data$cov[, i], gc = data$anno$gc,
                            mapa = data$anno$mapa),
                 use = usepoints & data$cov[, i] != 0, correctmapa = TRUE,
                 plot = file.path(destination.folder, "qc",
                                  paste0(colnames(data$cov)[i], ".png")))
        }
        ratios <- bplapply(i, NormalizeDOC, data, .tng, usepoints,
                           destination.folder, BPPARAM = bpparam)
        log2.read.counts <- matrix(unlist(ratios), ncol = length(sample.files))
    }, error = function(e) {
        stop("The GC-content and mapability normalization did not work due to ",
             "a failure to calculate loesses.\nThis can generally be solved ",
             "by using larger bin sizes.\nStopping execution of the remaining ",
             "part of the script...")
    })

    colnames(log2.read.counts) <- paste0("log2.", gsub(".compensated", "", f))

    ###################
    ## Create output ##
    ###################
    selection <- !data$anno$black & !is.na(data$anno$mapa) & data$anno$mapa > 0.2
    log2.read.counts <- data.frame(data$anno[selection, c("chr", "start", "end")],
                                   log2.read.counts[selection, , drop = FALSE],
                                   check.names = FALSE)

    ## Create table with corrected log2 values and write to file
    log2.read.counts <- cbind(with(log2.read.counts,
                                   cbind(Chromosome = as.vector(chr),
                                         Start = start,
                                         End = end,
                                         Feature = paste0(chr, ":",
                                                          paste0(start, "-",
                                                                 end)))),
                              log2.read.counts[, 4:ncol(log2.read.counts), drop = FALSE])

    ## Replace -/+Inf values to -/+large values for compatibility with IGV browser
    log2.read.counts[log2.read.counts == -Inf] <- -.Machine$integer.max/2
    log2.read.counts[log2.read.counts == Inf] <- .Machine$integer.max/2

    write.table(log2.read.counts[mixedorder(log2.read.counts$Chromosome), ],
                file.path(destination.folder, "log2_read_counts.igv"),
                sep = "\t", row.names = FALSE, quote = FALSE)

    ## Garbage collection
    rm(black, blacklist.file, data, f, gc, gc.content.file, i, log2.read.counts,
       mapa, mapability.file, NormalizeDOC, ratios, read.counts, selection, 
       usepoints)

    #############################################################################
    ## Calculate overlap with capture.regions.file for quality control purpose ##
    #############################################################################

    ## Calculate overlap with capture regions
    if (capture.regions.file != "not specified") {
        captured.bed <- read.table(capture.regions.file, as.is = TRUE,
                                   sep = "\t")
        colnames(captured.bed) <- c("Chromosome", "Start", "End")
        captured.grange <- with(captured.bed,
                                GRanges(seqnames = Chromosome,
                                        ranges = IRanges(start = Start,
                                                         end = End)))
        for (control.index in control.uniq.indices) {
            peak.bed <- read.table(file = paste0("peaks", control.index,
                                                 ".bed"),
                                   as.is = TRUE, sep = "\t")
            colnames(peak.bed) <- c("Chromosome", "Start", "End")
            peak.grange <- with(peak.bed,
                                GRanges(seqnames = Chromosome,
                                        ranges = IRanges(start = Start,
                                                         end = End)))

            overlap <- findOverlaps(captured.grange, peak.grange)

            cat("Number of capture regions covered by peaks in sample",
                length(unique(queryHits(overlap))), "\n")
            cat("Number of peaks covered by capture regions in sample",
                length(unique(subjectHits(overlap))), "\n")
            cat("Total number of capture regions in sample",
                sample.files[control.index], ": ", length(captured.grange),
                "\n")
            cat("Total number of peaks", sample.files[control.index], ":",
                length(peak.grange), "\n\n")
        }

        ## Garbage collection
        rm(captured.bed, captured.grange, capture.regions.file, overlap,
           peak.bed, peak.grange)
    }

    sink()

    ## Remove BamBaiPeaksFiles folder
    if (!keep.intermediairy.files) {
        unlink(file.path(destination.folder, "BamBaiPeaksFiles"),
               recursive = TRUE)
    }

    ## Garbage collection
    rm(control.index, control.uniq.indices, sample.files)

    cat("Total calculation time: ", Sys.time() - start.time, "\n\n")

    inputStructure <- list(sample.control = sample.control, ncpu = ncpu,
                           chromosomes = chromosomes, prefix = prefixes[1],
                           bin.file = bin.file, bin.size = bin.size)
    save(inputStructure, file = file.path(destination.folder, "input.Rdata"))

}