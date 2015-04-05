preCopywriteR <- function(output.folder, bin.size, ref.genome, prefix = "") {

    ## Make folder path absolute
    output.folder <- tools::file_path_as_absolute(output.folder)

    ## Checks

    # Check the existence of the output folder
    if (file.exists(output.folder) == FALSE) {
        stop(.wrap("The output folder could not be found. Please change the",
                   "path specified in", sQuote(output.folder)))
    } else {
        cat(.wrap("The output folder", sQuote(output.folder), "has been",
                  "detected"), "\n")
    }

    # Check if bin.size is factor of 1kb
    if (!.is.wholenumber(bin.size/1000)) {
        stop(.wrap("Please provide a bin size which is a multiple of 1000."))
    }

    ## Generate files with desired bin.size (mappa, GC, bed, blacklist)
    # Pre-define GC.mappa.grange and blacklist.grange variables to avoid them
    # from raising NOTES during R CMD CHECK (variables are loaded from file
    # below)
    GC.mappa.grange <- NULL
    blacklist.grange <- NULL
    
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
    
    for (chr in seqlevels(GC.mappa.grange)) {
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
        chr.bin <- cbind(seqnames = paste0(prefix, chr), start = start.bin,
                         end = end.bin, ATcontent = ATcontent.bin,
                         GCcontent = GCcontent.bin,
                         mappability = mappability.bin)
        chr.bin <- chr.bin[seq_len(nrow(chr.bin) - 1), ]
        custom.bin <- rbind(custom.bin, chr.bin)
    }
    custom.bin[, 2:ncol(custom.bin)] <-
        apply(custom.bin[, 2:ncol(custom.bin)], c(1, 2), as.numeric)
    custom.bin[, 4:ncol(custom.bin)] <-
        apply(custom.bin[, 4:ncol(custom.bin)], c(1, 2), round, 3)
    GC.mappa.grange <- makeGRangesFromDataFrame(custom.bin,
                                                keep.extra.columns = TRUE)
    
    cat(.wrap("Generated GC-content and mappability data at", bin.size, "bp",
              "resolution..."), "\n")
    
    # Generate blacklist object
    blacklist.grange <- as(blacklist.grange, "data.frame")
    blacklist.grange$seqnames <- paste0(prefix, blacklist.grange$seqnames)
    blacklist.grange <- makeGRangesFromDataFrame(blacklist.grange)
    cat(.wrap("Generated blacklist file..."), "\n")

    ## Create folder for output files
    file.name <- paste0(ref.genome, "_", bin.size/1000, "kb",
                        if (prefix != "") {paste0("_", prefix)})
    dir.create(file.path(output.folder, file.name))
    if (file.exists(file.path(output.folder, file.name)) == FALSE) {
        stop(.wrap("The output folder could not be created; please check that",
                   "you have permission to write in", sQuote(output.folder)))
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