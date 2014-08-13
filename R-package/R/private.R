
.is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol


.findCovFiles <- function(pattern,...) {
	list.files(pattern=pattern, ...)
}

.loadCovData <- function(files, gc=NULL, mapa=NULL, black=NULL, excludechr=NULL, datacol=5) {
	d <- lapply(files, read.delim, header=F)
	cnames <- .removeCommonFix(files)
	names(d) <- cnames

	#checks
	if(length(unique(sapply(d,nrow))) != 1) {
		stop("Files have variable number of lines and are not compatible")
	}

	#naturally sort data 
	require(gtools)
	d <- lapply(d, function(x) {
		x[,1] <- factor(x[,1], levels=mixedsort(levels(x[,1])))
		x[order(x[,1], x[,2]),]
		})

	#check row order in data files
	for (i in 2:length(files)) {
		if(!isTRUE(all.equal(d[[1]][,1:3], d[[i]][,1:3], check.attributes=FALSE))) stop(paste("Datafiles 1 and", i,"are not in the same order"))
	}

	cov <- sapply(d, function(x) x[!(x[,1] %in% excludechr),datacol])
	anno <- data.frame(droplevels(d[[1]][!(d[[1]][,1] %in% excludechr),1:3]))
	colnames(anno) <- c("chr","start","end")
	annoid <- paste(anno$chr, anno$start, anno$end, sep=":")
	if(!is.null(gc)) {
		anno$gc <- gc[match(annoid, paste(gc$chr, gc$start, gc$end, sep=":")),"gc"]
	}
	if(!is.null(mapa)) {
		anno$mapa <- mapa[match(annoid, paste(mapa$chr, mapa$start, mapa$end, sep=":")),"mapa"]
	}
	if(!is.null(black)) {
		suppressPackageStartupMessages(require(IRanges))
		anno$black <- rep(F, nrow(anno))
		for (c in levels(anno$chr)) {
			ar <- IRanges(anno$start[anno$chr == c], anno$end[anno$chr == c])
			br <- IRanges(black$start[black$chr == c], black$end[black$chr == c])
			anno$black[anno$chr == c] <- ar %over% br
		}
	}
	list(cov=cov, anno=anno)
}

.normalizeSampleGC <- function(x, gc, maxpoints=15000, plot=F) {
	x <- 2^x
	valid <- is.finite(x) & !is.na(gc) #& x > 0
	use <- which(valid)
	if(length(use) > maxpoints) use <- sample(use, maxpoints)

	df <- data.frame(x=x[use],gc=gc[use])
	fit <- loess(x~gc, data=df)
	normv <- rep(NA, length(x))
	normv[valid] <- predict(fit, data.frame(gc=gc[valid]))
	if(plot) {
		plot(x~gc, data=df)
		lines(normv, col=2)
	}
	log2(x/(normv/median(normv, na.rm=T)))
}

.normalizeGC <- function(ratios) {
	m <- ratios$ratios
	for (i in 1:ncol(m))
	m[,i] <- .normalizeSampleGC(m[,i],ratios$anno$X5_pct_gc)
	list(ratios=m, anno=ratios$anno)
}

.tng <- function(df, use, correctmapa=TRUE,  plot=NULL, verbose=T) {
	#tests
	if(!is.logical(use) && length(use) ==nrow(df))
		stop("use should be logicval vector with same size as df")
	#df colums?

	if(!is.null(plot)) {
		if(!is.logical(plot)) {
			if(verbose) cat("Plotting to file", plot,"\n")
			png(plot, width=700, height=1400)
			par(mfrow=c(2,1))
			on.exit(dev.off())
			plot <- TRUE
		} else if(plot) {
			par(mfrow=c(2,1))
		}
	}

	#exclude contains the points to exclude in the 
	#fitting (usually sex chromosomes and blacklisted regions)
	# gc fits  also excludes the low mappability data

	#correct gc using double lowess
	gcuse <- (use & !is.na(df$mapa) & df$mapa > .8 & !is.na(df$gc) & df$gc > 0)
	rough <- loess(count ~ gc, data=df, subset=gcuse, span = 0.03)
	i <- seq(0, 1, by = 0.001)
	final <- loess(predict(rough, i) ~ i, span = 0.3)
	normv <- predict(final, df$gc)
	df$countgcloess <- df$count/(normv/median(normv, na.rm=T))

	if(plot) {
		plot(count ~ gc, data=df, subset=gcuse, ylim=quantile(df$count[gcuse], c(0.0001, .999)), xlim=c(0,1), pch=".")
		points(count ~ gc, data=df, subset=!gcuse, col=rgb(1,0,0,.3), pch=".")
		lines(i, predict(rough, i), col="green")
		points(df$gc, normv, col="red", pch=".")
	}

	#correct mapa using linear function that intercepts zero
	#if(correctmapa) {
	#mapause <- (use & !is.na(df$mapa))
	#lm(countgcloess~0+mapa, data=df, subset=mapause) ->fll
	#if(verbose) print(summary(fll))

	#if (plot) {
	#	plot(countgcloess ~ mapa, data=df, subset=mapause, ylim=quantile(df$countgcloess, c(0.0001, .999), na.rm=T), pch=".")
	#	points(countgcloess ~ mapa, data=df, subset=!mapause, col=rgb(1,0,0,.3), pch=".")
	#	abline(0, fll$coef, col=2)
	#}

	#correct mapa using double lowess -> paired end sequencing
	if(correctmapa) {
		mapause <- (use & !is.na(df$mapa))
		rough <- loess(countgcloess ~ mapa, data=df, subset=mapause, span = 0.03)
		i <- seq(0, 1, by = 0.001)
		final <- loess(predict(rough, i) ~ i, span = 0.3)
		normv <- predict(final, df$mapa)
		df$countgcmapaloess <- df$countgcloess/(normv/median(normv, na.rm=T))
	
		if (plot) {
			plot(countgcloess ~ mapa, data=df, subset=mapause, ylim=quantile(df$countgcloess[mapause], c(0.0001, .999), na.rm=T), xlim=c(0,1), pch=".")
			points(countgcloess ~ mapa, data=df, subset=!mapause, col=rgb(1,0,0,.3), pch=".")
			lines(i, predict(rough, i), col="green")
			points(df$mapa, normv, subset = (!is.na(df$mapa) && !is.na(normv)), col="red", pch=".")
		}
	
		return(log2(df$countgcmapaloess / median(df$countgcmapaloess[use], na.rm=T)))
	} else {
		#corerct agains median value (exluding sex chr)
		log2(df$countgcloess / median(df$countgcloess[use], na.rm=T))
	}
}


.removeCommonFix <- function(names, distance=1) {
	l <- strsplit(names,"")

	#clip prefix
	pclip <- 1
	while( length(unique(sapply(l, "[", pclip))) <= distance) {
		pclip <- pclip + 1
	}

	#reverse strings for end clip pos
	l <- lapply(l, rev)
	eclip <- 1
	while( length(unique(sapply(l, "[", eclip))) <= distance) {
		eclip <- eclip + 1
	}

	sapply(names, function(x) substr(x, pclip, nchar(x) - eclip), USE.NAMES=F)

}

