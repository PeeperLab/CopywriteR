.is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
    abs(x - round(x)) < tol
}

.peakCutoff <- function(cov, fdr.cutoff = 1e-04, k = 2:150) {
    length.y <- length(cov)
    y <- tabulate(cov)
    y <- append(y, length.y - sum(y), 0)
    names(y) <- 0:(length(y) - 1)
    z <- sort(y, decreasing = TRUE)
    # Add 1 as otherwise 0 is possible outcome
    lambda <- as.integer(names(z[1])) + 1
    second.largest <- as.integer(names(z[2])) + 1
    # If second largest value after 0 is not 1, use other lambda
    if (lambda == 1 & second.largest != 2) {
        lambda <- second.largest + 1
    }
    # Add 1 to y[1]; otherwise if y[1] = 0 outcome is 0
    n <- exp(log(y[1] + 1) - dpois(1, lambda, log = TRUE))
    exp.fd <- n * ppois(k - 1, lambda, lower.tail = FALSE)
    obs.d <- rep(0, length(k))
    names.y <- as.integer(names(y))
    for (i in seq_along(k)) {
        tmp <- sum(y[names.y >= k[i]])
        if (tmp != 0) {
            obs.d[i] <- tmp
        } else {
            break
        }
    }
    FDR <- ifelse(obs.d == 0, 0, exp.fd/obs.d)
    fdr.ok <- which(FDR < fdr.cutoff)
    # Return 0 is there are no cutoffs low enough
    if (length(fdr.ok) < 1) {
        return(0)
    }
    fdr.chosen <- fdr.ok[1]
    tmp <- k[fdr.chosen - 1] + (FDR[fdr.chosen - 1] - fdr.cutoff)/(FDR[fdr.chosen - 1] - FDR[fdr.chosen])
    ## Always give output 0 if tmp is numeric(0)
    if (length(tmp) > 0) {
        tmp
    } else {
        0
    }
}

.tng <- function(df, use, correctmappa = TRUE, plot = NULL, verbose = TRUE) {
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
    gcuse <- (use & !is.na(df$mappa) & df$mappa > 0.8 & !is.na(df$gc) & df$gc > 0)
    rough <- loess(count ~ gc, data = df, subset = gcuse, span = 0.03)
    i <- seq(0, 1, by = 0.001)
    final <- loess(predict(rough, i) ~ i, span = 0.3)
    normv <- predict(final, df$gc)
    df$countgcloess <- df$count/(normv/median(normv, na.rm = TRUE))

    if (plot) {
        plot(count ~ gc, data = df, subset = gcuse,
             ylim = quantile(df$count[gcuse], c(1e-04, 0.999)), xlim = c(0, 1),
             pch = ".")
        points(count ~ gc, data = df, subset = !gcuse, col = rgb(1, 0, 0, 0.3),
               pch = ".")
        lines(i, predict(rough, i), col = "green")
        points(df$gc, normv, col = "red", pch = ".")
    }

    #correct mappa using linear function that intercepts zero
    #if(correctmappa) {
    #mappause <- (use & !is.na(df$mappa))
    #lm(countgcloess~0+mappa, data=df, subset=mappause) ->fll
    #if(verbose) print(summary(fll))

    #if (plot) {
    #  plot(countgcloess ~ mappa, data=df, subset=mappause,
    #       ylim=quantile(df$countgcloess, c(0.0001, .999), na.rm=T), pch=".")
    #  points(countgcloess ~ mappa, data=df, subset=!mappause,
    #         col=rgb(1,0,0,.3), pch=".")
    #  abline(0, fll$coef, col=2)
    #}

    #correct mappa using double lowess -> paired end sequencing
    if (correctmappa) {
        mappause <- (use & !is.na(df$mappa))
        rough <- loess(countgcloess ~ mappa, data = df, subset = mappause,
                       span = 0.03)
        i <- seq(0, 1, by = 0.001)
        final <- loess(predict(rough, i) ~ i, span = 0.3)
        normv <- predict(final, df$mappa)
        df$countgcmappaloess <- df$countgcloess/(normv/median(normv,
                                                             na.rm = TRUE))

        if (plot) {
            plot(countgcloess ~ mappa, data = df, subset = mappause,
                 ylim = quantile(df$countgcloess[mappause], c(1e-04, 0.999),
                                 na.rm = TRUE), xlim = c(0, 1), pch = ".")
            points(countgcloess ~ mappa, data = df, subset = !mappause,
                   col = rgb(1, 0, 0, 0.3), pch = ".")
            lines(i, predict(rough, i), col = "green")
            subset.points <- !is.na(df$mappa) && !is.na(normv)
            points(df$mappa[subset.points], normv[subset.points], col = "red",
                   pch = ".")
        }

        return(log2(df$countgcmappaloess/median(df$countgcmappaloess[use],
                                               na.rm = TRUE)))
    } else {
        #corerct agains median value (exluding sex chr)
        log2(df$countgcloess/median(df$countgcloess[use], na.rm = TRUE))
    }
}

.wrap <- function(...) {
    file.sep <- .Platform$file.sep
    splitted <- paste(unlist(strsplit(paste(...), split = file.sep)),
                                      collapse = paste0(file.sep, " "))
    splitted.pasted <- paste(strwrap(paste(splitted),
                                     exdent = 2), collapse = "\n")
    gsub(paste0(file.sep, " "), file.sep, splitted.pasted)
}
