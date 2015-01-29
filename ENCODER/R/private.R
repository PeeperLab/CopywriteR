.is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol

.loadCovData <- function(files, gc = NULL, mapa = NULL, black = NULL, excludechr = NULL, datacol = 5) {
  d <- lapply(files, read.delim, header = F)
  cnames <- .removeCommonFix(files)
  names(d) <- cnames

  #checks
  if (length(unique(sapply(d, nrow))) != 1) {
    stop("Files have variable number of lines and are not compatible")
  }

  #naturally sort data 
  d <- lapply(d, function(x) {
    x[, 1] <- factor(x[, 1], levels = mixedsort(levels(x[, 1])))
    x[order(x[, 1], x[, 2]), ]
  })

  #check row order in data files
  for (i in 2:length(files)) {
    if (!isTRUE(all.equal(d[[1]][, 1:3], d[[i]][, 1:3], check.attributes = FALSE))) 
      stop(paste("Datafiles 1 and", i, "are not in the same order"))
  }

  cov <- sapply(d, function(x) x[!(x[, 1] %in% excludechr), datacol])
  anno <- data.frame(droplevels(d[[1]][!(d[[1]][, 1] %in% excludechr), 1:3]))
  colnames(anno) <- c("chr", "start", "end")
  annoid <- paste(anno$chr, anno$start, anno$end, sep = ":")
  if (!is.null(gc)) {
    anno$gc <- gc[match(annoid, paste(gc$chr, gc$start, gc$end, sep = ":")), "gc"]
  }
  if (!is.null(mapa)) {
    anno$mapa <- mapa[match(annoid, paste(mapa$chr, mapa$start, mapa$end, sep = ":")), "mapa"]
  }
  if (!is.null(black)) {
    anno$black <- rep(F, nrow(anno))
    for (c in levels(anno$chr)) {
      ar <- IRanges(anno$start[anno$chr == c], anno$end[anno$chr == c])
      br <- IRanges(black$start[black$chr == c], black$end[black$chr == c])
      anno$black[anno$chr == c] <- ar %over% br
    }
  }
  list(cov = cov, anno = anno)
}

.peakCutoff.dynamic <- function(cov, fdr.cutoff = 0.0001, k = 2:150) {
	length.y <- length(cov)
	y <- tabulate(cov)
	y <- append(y, length.y - sum(y), 0)
	names(y) <- 0:(length(y) - 1)
	z <- sort(y, decreasing = TRUE)
	lambda <- as.integer(names(z[1])) + 1 # Add one as otherwise 0 is possible outcome
	second.largest <- as.integer(names(z[2])) + 1
	if (lambda == 1 & second.largest != 2) { # If second largest value after 0 is not 1, use other lambda
		lambda <- second.largest + 1
	}
	n <- exp(log(y[1] + 1) - dpois(1, lambda, log = TRUE)) # 1 Added to y[1]; otherwise if y[1] = 0 outcome is 0
	exp.fd <- n * ppois(k - 1, lambda, lower.tail = FALSE)
	obs.d <- integer(length(k))
	for (i in seq_along(k)) {
		obs.d[i] <- sum(y[as.integer(names(y)) >= k[i]])
	}
	FDR <- ifelse(obs.d == 0, 0, exp.fd/obs.d)
	fdr.ok <- which(FDR < fdr.cutoff)
	if (length(fdr.ok) < 1) 
		stop("No cutoff with low enough FDR found")
	fdr.chosen <- fdr.ok[1]
	k[fdr.chosen - 1] + (FDR[fdr.chosen - 1] - fdr.cutoff)/(FDR[fdr.chosen - 1] - FDR[fdr.chosen])
}

.removeCommonFix <- function(names, distance = 1) {
  l <- strsplit(names, "")

  #clip prefix
  pclip <- 1
  while (length(unique(sapply(l, "[", pclip))) <= distance) {
    pclip <- pclip + 1
  }

  #reverse strings for end clip pos
  l <- lapply(l, rev)
  eclip <- 1
  while (length(unique(sapply(l, "[", eclip))) <= distance) {
    eclip <- eclip + 1
  }

  sapply(names, function(x) substr(x, pclip, nchar(x) - eclip), USE.NAMES = F)

}

.tng <- function(df, use, correctmapa = TRUE, plot = NULL, verbose = T) {
  #tests
  if (!is.logical(use) && length(use) == nrow(df)) 
    stop("use should be logicval vector with same size as df")
  #df colums?
  
  if (!is.null(plot)) {
    if (!is.logical(plot)) {
      if (verbose) 
        cat("Plotting to file", plot, "\n")
      png(plot, width = 700, height = 1400)
      par(mfrow = c(2, 1))
      on.exit(dev.off())
      plot <- TRUE
    } else if (plot) {
      par(mfrow = c(2, 1))
    }
  }

  #exclude contains the points to exclude in the 
  #fitting (usually sex chromosomes and blacklisted regions)
  #gc fits also excludes the low mappability data

  #correct gc using double lowess
  gcuse <- (use & !is.na(df$mapa) & df$mapa > 0.8 & !is.na(df$gc) & df$gc > 0)
  rough <- loess(count ~ gc, data = df, subset = gcuse, span = 0.03)
  i <- seq(0, 1, by = 0.001)
  final <- loess(predict(rough, i) ~ i, span = 0.3)
  normv <- predict(final, df$gc)
  df$countgcloess <- df$count/(normv/median(normv, na.rm = T))

  if (plot) {
    plot(count ~ gc, data = df, subset = gcuse, ylim = quantile(df$count[gcuse], c(1e-04, 0.999)), xlim = c(0, 1), pch = ".")
    points(count ~ gc, data = df, subset = !gcuse, col = rgb(1, 0, 0, 0.3), pch = ".")
    lines(i, predict(rough, i), col = "green")
    points(df$gc, normv, col = "red", pch = ".")
  }

  #correct mapa using linear function that intercepts zero
  #if(correctmapa) {
  #mapause <- (use & !is.na(df$mapa))
  #lm(countgcloess~0+mapa, data=df, subset=mapause) ->fll
  #if(verbose) print(summary(fll))

  #if (plot) {
  #  plot(countgcloess ~ mapa, data=df, subset=mapause, ylim=quantile(df$countgcloess, c(0.0001, .999), na.rm=T), pch=".")
  #  points(countgcloess ~ mapa, data=df, subset=!mapause, col=rgb(1,0,0,.3), pch=".")
  #  abline(0, fll$coef, col=2)
  #}

  #correct mapa using double lowess -> paired end sequencing
  if (correctmapa) {
    mapause <- (use & !is.na(df$mapa))
    rough <- loess(countgcloess ~ mapa, data = df, subset = mapause, span = 0.03)
    i <- seq(0, 1, by = 0.001)
    final <- loess(predict(rough, i) ~ i, span = 0.3)
    normv <- predict(final, df$mapa)
    df$countgcmapaloess <- df$countgcloess/(normv/median(normv, na.rm = T))

    if (plot) {
      plot(countgcloess ~ mapa, data = df, subset = mapause, ylim = quantile(df$countgcloess[mapause], c(1e-04, 0.999), na.rm = T), xlim = c(0, 
        1), pch = ".")
      points(countgcloess ~ mapa, data = df, subset = !mapause, col = rgb(1, 0, 0, 0.3), pch = ".")
      lines(i, predict(rough, i), col = "green")
      points(df$mapa, normv, subset = (!is.na(df$mapa) && !is.na(normv)), col = "red", pch = ".")
    }

    return(log2(df$countgcmapaloess/median(df$countgcmapaloess[use], na.rm = T)))
  } else {
    #corerct agains median value (exluding sex chr)
    log2(df$countgcloess/median(df$countgcloess[use], na.rm = T))
  }
}