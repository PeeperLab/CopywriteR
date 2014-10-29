plotCNA <- function(destination.folder, set.nchrom, sample.plot) {

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
  
  ## Set nchrom
  if (missing(set.nchrom)) {
    nchrom <- inputStructure$nchrom
  } else {
    nchrom <- set.nchrom
  }
  
  prefix <- inputStructure$prefix
  ncpu <- inputStructure$ncpu
  
  
  ## Set sample.plot
  if (missing(sample.plot)) {
    all.samples <- unique(unlist(inputStructure$sample.control))
    sample.plot <- rbind(inputStructure$sample.control,
                         data.frame(samples = all.samples,
                                    controls = rep(NA, length(all.samples))))
  }
  sample.plot[, ] <- apply(sample.plot, c(1, 2), function(x) {
    if (!is.na(x)) {
      x <- paste0("log2.read.counts.",
                  unlist(strsplit(x, "/"))[length(unlist(strsplit(x, "/")))])
    } else {
      x <- x
    }
  })

  ## Read data
  log2.read.counts <- read.table(file = paste0(destination.folder,
                                               "log2_read_counts.txt"),
                                 sep = "\t", header = TRUE, check.names = FALSE)

  ## Remove prefix and convert chromosome names to integers
  log2.read.counts <-
    log2.read.counts[-which(rowSums(is.na(log2.read.counts[, -c(1:4)])) > 0), ]
  log2.read.counts[, 2] <- gsub(prefix, "", log2.read.counts[, 2])
  log2.read.counts[, 2] <- gsub("X", nchrom - 1, log2.read.counts[, 2])
  log2.read.counts[, 2] <- gsub("Y", nchrom, log2.read.counts[, 2])
  log2.read.counts[, 2] <- as.integer(log2.read.counts[, 2])
  
  ## Fix behaviour of DNAcopy with very low values
  log2.read.counts[, 5:ncol(log2.read.counts)] <-
    apply(log2.read.counts[, 5:ncol(log2.read.counts)], c(1, 2), function(x) {
    if (x < -100) {
      x <- -100
    } else if (x > 100) {
      x <- 100
    } else {
      x <- x
    }
  })
  
  ## Create table with values to be plotted
  if (all(vector(sample.plot) %in% colnames(log2.read.counts))) {
  } else {
    stop("One of the samples in sample.plot refers to a BAM file that has not ",
         "been processed in ENCODER. Please make sure that you have provided ",
         "the correct input files or re-run ENCODER accordingly.")
  }
  
  ## Apply DNAcopy
  CNA.object <- CNA(log2.read.counts[, 5:ncol(log2.read.counts)],
                    log2.read.counts$Chromosome,
                    rowMeans(log2.read.counts[, c("StartPos", "StopPos")]),
                    data.type = "logratio",
                    sampleid =
                      colnames(log2.read.counts)[5:ncol(log2.read.counts)])
 
  smoothed.CNA.object <- smooth.CNA(CNA.object)
  segment.CNA.object <- segment(smoothed.CNA.object, verbose = 1)
  
  


  
  









  colnames(log2.read.counts) <- paste0("log2.", gsub(".compensated", "", sampnames))





  ## Convert names in sample.plot to the corresponding names in log2.read.counts 

  ## Test whether all sample files have been processed before
  if (all(vector(sample.plot) %in% colnames(log2.read.counts))) {
  } else {
    stop("One of the samples in sample.plot refers to a BAM file that has not ",
         "been processed in ENCODER. Please make sure that you have provided ",
         "the correct input files or re-run ENCODER accordingly.")
  }

  # Run CGHcall
  raw <- make_cghRaw(read_count)
  prep <- preprocess(raw, maxmiss = 30, nchrom = inputStructure$nchrom)
  nor <-  normalize(prep, method = "median", smoothOutliers = TRUE)
  seg <-  segmentData(nor, method = "DNAcopy", nperm = 2000, undo.splits = "sdundo", min.width = 5, undo.SD = 1, clen = 25, relSDlong = 5)
  segnorm <- postsegnormalize(seg, inter = c(-1,1))
  listcalls <- CGHcall(segnorm, nclass = 5, robustsig = "yes", cellularity = 1, ncpus = inputStructure$ncpu)
  calls <- ExpandCGHcall(listcalls, segnorm, divide = 5, memeff = FALSE)
  save(calls, file = paste0(inputStructure$destination.folder, "calls.Rdata"))
  
  # Make whole-genome plots
  system(paste0("mkdir ", inputStructure$destination.folder, "Call_plots"))
  colnames_read_count <- gsub(" ", ".", colnames(read_count)[5:ncol(read_count)])
  for (i in 1:(ncol(read_count) - 4)) {
    setwd(paste0(inputStructure$destination.folder, "Call_plots/"))
    system(paste0("mkdir ",i, "_", colnames_read_count[i]))
    setwd(paste0(i, "_", colnames_read_count[i]))
    png(filename=paste0(i, "_freqonly_", colnames_read_count[i],".png"), width=2*480, height=480)
      plot(nor[,i], dotres = 1)
    dev.off()  
    pdf(file=paste0(i, "_freqonly_", colnames_read_count[i],".pdf"), width=2*7, height=7)
      plot(nor[,i], dotres = 1)
    dev.off()
    png(filename=paste0(i, "_", colnames_read_count[i],".png"), width=2*480, height=480)
      plot(calls[,i], dotres=1)
    dev.off()
    pdf(file=paste0(i, "_", colnames_read_count[i],".pdf"), width=2*7, height=7)
      plot(calls[,i], dotres=1)
    dev.off()
    png(filename=paste0(i, "_segmented_", colnames_read_count[i],".png"), width=2*480, height=480)
      plot(seg[,i], dotres=1)
    dev.off()
    pdf(file=paste0(i, "_segmented_", colnames_read_count[i],".pdf"), width=2*7, height=7)
      plot(seg[,i], dotres=1)
    dev.off()
      
    # Make per-genome plots
    for (j in 1:inputStructure$nchrom){
      print(j)
      png(filename=paste0("chromo_", j, "_", colnames_read_count[i],"_freqonly.png"), width=2*480, height=480)
        plot(nor[chromosomes(nor) == j,i], dotres=1)
      dev.off()
      pdf(file=paste0("chromo_", j, "_", colnames_read_count[i],"_freqonly.pdf"), width=2*7, height=7)
        plot(nor[chromosomes(nor) == j,i], dotres=1)
      dev.off()
      png(filename=paste0("chromo_", j, "_", colnames_read_count[i],"_seg.png"), width=2*480, height=480)
        plot(seg[chromosomes(seg) == j,i], dotres=1)
      dev.off()
      pdf(file=paste0("chromo_", j, "_", colnames_read_count[i],"_seg.pdf"), width=2*7, height=7)
        plot(seg[chromosomes(seg) == j,i], dotres=1)
      dev.off()
      png(filename=paste0("chromo_", j, "_", colnames_read_count[i],".png"), width=2*480, height=480)
        plot(calls[chromosomes(calls)==j,i], dotres=1)
      dev.off()
      pdf(file=paste0("chromo_", j, "_", colnames_read_count[i],".pdf"), width=2*7, height=7)
        plot(calls[chromosomes(calls)==j,i], dotres=1)
      dev.off()
    }
  }
}