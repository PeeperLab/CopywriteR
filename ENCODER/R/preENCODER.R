preENCODER <- function(black.GC.mapa.folder, output.folder, bin.size,
                       reference) {

  ## Make folder path absolute
  black.GC.mapa.folder <- tools::file_path_as_absolute(black.GC.mapa.folder)
  output.folder <- tools::file_path_as_absolute(output.folder)

  ## Checks
  # Check whether folders exist.
  if (!(reference == "hg19" || reference == "mm10" || reference == "mm9")) {
    stop("The reference is not recognised. Please provide a suitable",
         "reference.")
  }

  if (file.exists(black.GC.mapa.folder) == FALSE) {
    stop("The annotation folder could not be found. Please change your",
         "black.GC.mapa.folder path.")
  } else {
    cat("Reference folder", black.GC.mapa.folder, "detected", "\n")
  }

  if (file.exists(output.folder) == FALSE) {
    stop("The output.folder could not be found. Please change the output",
         "folder path.")
  } else {
    cat("Reference folder", output.folder, "detected", "\n")
  }

  # Check if bin.size is factor of 1kb
  if (!.is.wholenumber(bin.size/1000)) {
    stop("Please provide a bin.size which is a multiple of 1000.")
  }

  ## Generate files with desired bin.size (MAPA, GC, bed, blacklist)
  # Load blacklist, CG-content and mapability files
  load(file.path(black.GC.mapa.folder, "mapability.rda"))
  load(file.path(black.GC.mapa.folder, "GCcontent.rda"))
  bed_file <- read.table(file.path(black.GC.mapa.folder, "blacklist.bed"),
                         as.is = TRUE, sep = "\t")

  # Create bins with desired bin size
  MERGEBINNUMBER <- bin.size/1000

  newBin <- NULL
  options(warn = -1)
  options(scipen = 999)
  for (chr in unique(mapa$chromosome)) {
    col2 <- colMins(matrix(mapa$start[mapa$chromosome == chr],
                           nrow = MERGEBINNUMBER))
    col3 <- colMaxs(matrix(mapa$end[mapa$chromosome == chr],
                           nrow = MERGEBINNUMBER))
    tmp <- cbind(chr, col2, col3)
    tmp <- tmp[1:(nrow(tmp) - 1), ]
    newBin <- rbind(newBin, tmp)
  }
  options(scipen = 0)
  options(warn = 0)

  cat("Generated", bin.size, "bp bins for all chromosomes", "\n")

  # Create mapabillity file with desired bin size
  MERGEBINNUMBER <- bin.size/1000

  newMapa <- NULL
  options(warn = -1)
  options(scipen = 999)
  for (chr in unique(mapa$chromosome)) {
    col2 <- colMins(matrix(mapa$start[mapa$chromosome == chr],
                           nrow = MERGEBINNUMBER))
    col3 <- colMaxs(matrix(mapa$end[mapa$chromosome == chr],
                           nrow = MERGEBINNUMBER))
    col4 <- colMeans(matrix(mapa$mapability[mapa$chromosome == chr],
                            nrow = MERGEBINNUMBER))
    tmp <- cbind(chr, col2, col3, col4)
    tmp <- tmp[1:(nrow(tmp) - 1), ]
    newMapa <- rbind(newMapa, tmp)
  }
  options(scipen = 0)
  options(warn = 0)

  cat("Generated mapability file for bin size of", bin.size, "bp", "\n")

  # Create GC-content file with desired bin size
  newGC <- NULL
  options(warn = -1)
  options(scipen = 999)
  for (chr in unique(mapa$chromosome)) {
    col2 <- colMins(matrix(GC$start[GC$chromosome == chr],
                           nrow = MERGEBINNUMBER))
    col3 <- colMaxs(matrix(GC$end[GC$chromosome == chr], nrow = MERGEBINNUMBER))
    col4 <- colMeans(matrix(GC$ATcontent[GC$chromosome == chr],
                            nrow = MERGEBINNUMBER))
    col5 <- colMeans(matrix(GC$GCcontent[GC$chromosome == chr],
                            nrow = MERGEBINNUMBER))
    tmp <- cbind(chr, col2, col3, col4, col5)
    tmp <- tmp[1:(nrow(tmp) - 1), ]
    newGC <- rbind(newGC, tmp)
  }
  options(scipen = 0)
  options(warn = 0)

  cat("Generated GC-content file for bin size of", bin.size, "bp", "\n")

  ## Create folder for output files
  file_name <- paste0(reference, "_", bin.size/1000, "kb/")
  dir.create(file.path(output.folder, file_name))
  if (file.exists(file.path(output.folder, file_name)) == FALSE) {
    stop("No output folder created, please check argument output.folder and",
         "the corresponding folder permissions.")
  }

  ## Write files to folder
  # Blacklist
  write.table(bed_file, file = file.path(output.folder, file_name,
                                         "blacklist.bed"), quote = FALSE,
              row.names = FALSE, col.names = FALSE, sep = "\t")
  # GC-content
  write.table(newGC, file = file.path(output.folder, file_name,
                                      "GC_content.bed"), quote = FALSE,
              row.names = FALSE, col.names = FALSE, sep = "\t")
  # Mapability
  write.table(newMapa, file = file.path(output.folder, file_name,
                                        "mapability.bed"), quote = FALSE,
              row.names = FALSE, col.names = FALSE, sep = "\t")
  # bed file with bins
  write.table(newBin, file = file.path(output.folder, file_name, "bins.bed"),
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

}
