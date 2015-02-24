preCopywriteR <- function(output.folder, bin.size, ref.genome) {

    ## Make folder path absolute
    output.folder <- tools::file_path_as_absolute(output.folder)

    ## Checks
    # Check whether ref.genome is supported
    if (!(ref.genome == "hg19" || ref.genome == "mm10" || ref.genome == "mm9")) {
        stop("The ref.genome is not recognised. Please provide a suitable",
             "ref.genome.")
    }

    # Check the existence of the output folder
    if (file.exists(output.folder) == FALSE) {
        stop("The output.folder could not be found. Please change the output",
             "folder path.")
    } else {
        cat("The output folder", output.folder, "has been detected", "\n")
    }

    # Check if bin.size is factor of 1kb
    if (!.is.wholenumber(bin.size/1000)) {
        stop("Please provide a bin size which is a multiple of 1000.")
    }

    ## Generate files with desired bin.size (mappa, GC, bed, blacklist)
    # Pre-define mappa and GC variables to avoid them raising NOTES during
    # R CMD CHECK (variables are loaded from file below)
    mappa <- NULL
    GC <- NULL
    
    # Load GC.grange, GC.mappa.grange and blacklist.grange variables
    black.GC.mappa.folder <- getPathHelperFiles(ref.genome)
    load(file.path(black.GC.mappa.folder, "GC_mappability.rda"))
    load(file.path(black.GC.mappa.folder, "blacklist.rda"))

    # Create bins with desired bin size
    MERGEBINNUMBER <- bin.size/1000

    # Create bins
    custom.bin <- data.frame()
    
    old.options <- options()
    on.exit(options(old.options))
    options(warn = -1, scipen = 999)
    
    for (chr in seqnames(GC.mappa.grange)@values) {
        selection <- as(seqnames(GC.mappa.grange) == chr, "vector")
        start.bin <- colMins(matrix(start(GC.mappa.grange)[selection],
                                    nrow = MERGEBINNUMBER))
        end.bin <- colMaxs(matrix(end(GC.mappa.grange)[selection],
                                  nrow = MERGEBINNUMBER))
        ATcontent.bin <- colMeans(matrix(GC.mappa.grange$ATcontent[selection],
                                         nrow = MERGEBINNUMBER))
        GCcontent.bin <- colMeans(matrix(GC.mappa.grange$GCcontent[selection],
                                         nrow = MERGEBINNUMBER))
        mappability.bin <- colMeans(matrix(GC.mappa.grange$mappability[selection],
                                           nrow = MERGEBINNUMBER))
        chr.bin <- cbind(seqnames = chr, start = start.bin, end = end.bin,
                     ATcontent = ATcontent.bin, GCcontent = GCcontent.bin,
                     mappability = mappability.bin)
        chr.bin <- chr.bin[1:(nrow(chr.bin) - 1), ]
        custom.bin <- rbind(custom.bin, chr.bin)
    }
    custom.bin[, 2:ncol(custom.bin)] <-
        apply(custom.bin[, 2:ncol(custom.bin)], c(1, 2), as.numeric)
    custom.bin[, 4:ncol(custom.bin)] <-
        apply(custom.bin[, 4:ncol(custom.bin)], c(1, 2), round, 3)
    GC.mappa.grange <- makeGRangesFromDataFrame(custom.bin,
                                                keep.extra.columns = TRUE)
    
    cat("Generated", bin.size, "bp bins for all chromosomes with matching",
        "GC-content and mappability data", "\n")

    ## Create folder for output files
    file.name <- paste0(ref.genome, "_", bin.size/1000, "kb")
    dir.create(file.path(output.folder, file.name))
    if (file.exists(file.path(output.folder, file.name)) == FALSE) {
        stop("The output folder could not be created, please check argument",
             "output.folder and the corresponding folder permissions.")
    }

    ## Write files to folder
    # Blacklist
    save(blacklist.grange, file = file.path(output.folder, file.name,
                                            "blacklist.rda"), compress = "xz")
    # GC-content and mappability
    save(GC.mappa.grange,
         file = file.path(output.folder, file.name, "GC_mappability.rda"),
         compress = "xz")

}