"plot.SPSbalan" <-
function (x, ...) 
{
    if (x$youtype == "contin" && x$bins > 1) {
        form <- as.formula(paste(x$trtm, "~", x$xvar, "|", x$qbin))
        bwplot(form, data = x$df3)
    }
}

"plot.SPSloess" <-
function (x, tcol = "blue", ucol = "red", dcol = "green3", ...) 
{
    # First plot...
    plot(x$losub0[, 1], x$losub0[, 2], type = "p", ann = FALSE, 
        col = tcol, pch = 21)
    lines(x$losub0[, 1], x$losub0[, 4], type = "l", col = tcol, 
        lty = "solid")
    points(x$losub1[, 1], x$losub1[, 2], type = "p", col = ucol, 
        pch = 24)
    lines(x$losub1[, 1], x$losub1[, 4], type = "l", col = ucol, 
        lty = "dashed")
    title(main = paste("Outcomes and PS Loess Fit, Span =", x$span), 
        ylab = paste("Observed and Smoothed", x$yvar), xlab = "Estimated Propensity Score")
    opar <- par(ask = dev.interactive(orNone = TRUE))
    # Second plot...
    plot(x$losub0[, 1], x$losub0[, 4], type = "l", ann = FALSE, 
        col = tcol, lty = "solid")
    lines(x$losub1[, 1], x$losub1[, 4], type = "l", col = ucol, 
        lty = "dashed")
    title(main = paste("PS Loess Fit, Span =", x$span), ylab = paste("Smoothed", 
        x$yvar), xlab = "Estimated Propensity Score")
    # Third plot...
    plot(x$logrid[, 1], 0.1 * x$logrid[, 10], type = "p", ann = FALSE, 
        pch = 21)
    lines(x$logrid[, 1], x$logrid[, 11], type = "l", col = dcol, 
        lty = "solid")
    title(main = "PS Probability Density", xlab = "Estimated Propensity Score")
    par(opar)
    NULL
}

"plot.SPSoutco" <-
function (x, ...) 
{
    # First plot...
    barplot(t(as.matrix(x$binfreq[, 2:3])), beside = TRUE, names.arg = as.vector(x$binfreq[, 
        1]), main = "Frequencies by PS Bin")
    if (x$youtype == "contin") {
        opar <- par(ask = dev.interactive(orNone = TRUE))
        # Second plot...
        barplot(t(as.matrix(x$binmean[, 2:3])), beside = TRUE, 
            names.arg = as.vector(x$binmean[, 1]), main = "Outcome Means by PS Bin")
        # Third plot...
        plot(x$pbindif, x$pbinsde, ann = FALSE, type = "n", ylim = c(0, 
            max(x$pbinsde, na.rm = TRUE)))
        symbols(x$pbindif, x$pbinsde, circles = x$pbinsiz, inches = 0.25, 
            add = TRUE)
        par(lty = 1)
        abline(v = 0)
        par(lty = 2)
        abline(v = x$awbdif)
        par(lty = 3)
        abline(v = x$wwbdif)
        title(main = paste("Supervised Propensity Scoring: BINS =", 
            x$bins), xlab = "Within Bin Treatment Difference", 
            ylab = "Difference Standard Deviation", sub = paste("Number of Informative PS Bins =", 
                length(na.omit(x$pbinsde))))
        par(opar)
    }
}

"plot.SPSsmoot" <-
function (x, tcol = "blue", ucol = "red", dcol = "green3", ...) 
{
    # First plot...
    plot(x$spsub0[, 1], x$spsub0[, 2], type = "p", ann = FALSE, 
        col = tcol, pch = 21)
    lines(x$spsub0[, 1], x$spsub0[, 4], type = "l", col = tcol, 
        lty = "solid")
    points(x$spsub1[, 1], x$spsub1[, 2], type = "p", col = ucol, 
        pch = 24)
    lines(x$spsub1[, 1], x$spsub1[, 4], type = "l", col = ucol, 
        lty = "dashed")
    title(main = paste("Outcomes and PS Smoothing Splines, df =", 
        x$df), ylab = paste("Observed and Smoothed", x$yvar), 
        xlab = "Estimated Propensity Score")
    opar <- par(ask = dev.interactive(orNone = TRUE))
    # Second plot...
    plot(x$spsub0[, 1], x$spsub0[, 4], type = "l", ann = FALSE, 
        col = tcol, lty = "solid")
    lines(x$spsub1[, 1], x$spsub1[, 4], type = "l", col = ucol, 
        lty = "dashed")
    title(main = paste("PS Smoothing Splines, df =", x$df), ylab = paste("Smoothed", 
        x$yvar), xlab = "Estimated Propensity Score")
    # Third plot...
    plot(x$ssgrid[, 1], 0.1 * x$ssgrid[, 10], type = "p", ann = FALSE, 
        pch = 21)
    lines(x$ssgrid[, 1], x$ssgrid[, 11], type = "l", col = dcol, 
        lty = "solid")
    title(main = "PS Probability Density", xlab = "Estimated Propensity Score")
    par(opar)
    NULL
}

"plot.UPSaltdd" <- 
function(x, breaks="Sturges", ...)
{
    opar <- par(no.readonly = TRUE, ask = dev.interactive(orNone = TRUE))
    on.exit(par(opar))
    if( x$NNobj=="NA" ) par(mfrow = c(2, 1)) else par(mfrow = c(2, 2))
    nx <- length(x$altdd[,1])
    if( x$NNobj!="NA") {
        xmin <- min(x$alxmin, x$nnlxmin)
        xmax <- max(x$alxmax, x$nnlxmax)
        }
    else {
        xmin <- x$alxmin
        xmax <- x$alxmax
        }
    plot(x$altdcdf, x$qq, ann = FALSE, type = "l", xlim = c(xmin,xmax))
    par(lty=1)
    title(main=paste("Artificial LTD CDF"),
        ylab = "Probability Less Than",
        xlab =paste("From ", nx, "Random, Informative Clusters within", x$reps, "reps of", x$clus))
    hist(x$altdcdf, breaks=breaks, xlim=c(xmin,xmax),
        main=paste("Artificial LTD Histogram"), ylab = "Probability Density",
        xlab = paste("Artificial Effects of", x$trtm, "on", x$yvar, "using Random Clusters.")
        )
    if( x$NNobj!="NA" ) {
        nx <- length(x$nnltdd[,1])
        plot(x$nnltdcdf, x$nq, ann = FALSE, type = "l", xlim = c(xmin,xmax))
        par(lty=1)
        title(main=paste("Nearest Neighbor LTD CDF"),
            ylab = "Probability Less Than",
            xlab =paste("From ", nx, "Informative Clusters within", x$clus))
        hist(x$nnltdcdf, breaks=breaks, xlim=c(xmin,xmax),
            main=paste("Nearest Neighbor LTD Histogram"), ylab = "Probability Density",
            xlab = paste("Local Effects of", x$trtm, "on", x$yvar, "using Relevant Clusters.")
            )
        }
    }

"plot.UPShclus" <-
function (x, ...) 
{
    if (x$method == "diana") {
        plot(x$upshcl, main = "Unsupervised Divisive Hierarchy", 
            sub = paste("Divisive Coefficient = ", round(x$upshcl$dc, 
                digits = 2)))
    }
    else if (x$method == "agnes") {
        plot(x$upshcl, main = "Unsupervised Agglomerative Hierarchy", 
            sub = paste("Agglomerative Coefficient = ", round(x$upshcl$ac, 
                digits = 2)))
    }
    else {
        plclust(x$upshcl, main = "Unsupervised Clustering Hierarchy")
    }
}

"plot.UPSivadj" <-
function (x, ...) 
{
    if (x$youtype == "contin") {
        inchsz <- 2 * max(x$pbinsiz)^2/x$symsiz
        plot(x$pbinpsp, x$pbinout, ann = FALSE, type = "n")
        symbols(x$pbinpsp, x$pbinout, circles = x$pbinsiz, inches = inchsz, 
            add = TRUE)
        par(lty = 1)
        if (x$actclust > 1) 
            abline(x$ivfit)
        title(main = paste("Unsupervised IV Clusters =", x$actclust), 
            xlab = "Within-Cluster Treatment Percentage (PS)", 
            ylab = "Observed LATE", sub = "Symbol Area proportional to Cluster Size")
    }
}

"plot.UPSnnltd" <-
function (x, pballs = TRUE, nnplot = "snob", nnalpha = 1.4, ...) 
{
    if (x$youtype == "contin") {
        UPSpars <- get("UPSaccum.pars")
        if (nnplot != "snob" && nnplot != "dens" && nnplot != "cdf" && nnplot != "seq")
            nnplot <- "all"
        opar <- par(no.readonly = TRUE, ask = dev.interactive(orNone = TRUE))
        on.exit(par(opar))
        if (nnplot == "all")
            par(mfrow=c(2,2))
        else
            par(mfrow=c(1,1))
        if (nnplot == "all" || nnplot == "dens" || nnplot == "cdf" || nnplot == "seq") {
            nm <- complete.cases(x$pbinsde)
            nn <- as.vector(x$pbinsiz[nm]^2)
            xx <- as.vector(x$pbindif[nm])
            ww <- as.vector(x$pbinsde[nm]^-2)
            ww <- round( sum(nn)* ww / sum(ww) )
            ltdfit <- ssden(~xx, weight = ww, alpha = nnalpha, maxiter = 30)
            xp <- seq(as.numeric(UPSpars[9]),as.numeric(UPSpars[10]), len=101)
        }
        if (nnplot == "all" || nnplot == "seq" || nnplot == "snob") {
            inchsz <- 2 * max(x$pbinsiz)^2/x$symsiz
            plot(x$pbindif, x$pbinsde, ann = FALSE, type = "n",
                ylim = c(0, as.numeric(UPSpars[8])),
                xlim = c(as.numeric(UPSpars[9]), as.numeric(UPSpars[10])))
            symbols(x$pbindif, x$pbinsde, circles = x$pbinsiz, inches = inchsz, 
                add = TRUE)
            par(lty = 2, col = "red")
            abline(v = x$awbdif)
            par(lty = 3, col = "green3")
            abline(v = x$wwbdif)
            par(lty = 1, col = "black")
            abline(v = 0)
            if (pballs) {
                symbols(x$awbdif, x$awbsde, circles = x$awbsde, inches = 0.25, 
                    bg = "red", add = TRUE)
                symbols(x$wwbdif, x$wwbsde, circles = x$wwbsde, inches = 0.25, 
                    bg = "green3", add = TRUE)
            }
            title(main = paste("Unsupervised NN/LTD Snow Balls for",
                length(na.omit(x$pbinsde)), "of", x$numclust, "Clusters"),
                xlab = "Local Treatment Difference (LTD)", ylab = "LTD Standard Error")
            if (nnplot == "seq") {
                cat("\nPress the Enter key to view the DENS nnplot...")
                scan()
            }
        }
        if (nnplot == "all" || nnplot == "seq" || nnplot == "dens") {
            plot(xp, dssden(ltdfit, xp), ann=FALSE, type="l", lty=1, col="blue")
            par(lty=2, col="red")
            abline(v=x$awbdif)
            par(lty=3, col="green3")
            abline(v=x$wwbdif)
            par(lty=1, col="black")
            abline(v=0)           
            title(main=paste("NN/LTD Weighted Density for", length(na.omit(x$pbinsde)),
               "Informative Clusters"), xlab = "Local Treatment Difference (LTD)",
                ylab = "gss Probability Density")
            if (nnplot == "seq") {
                cat("\nPress the Enter key to view the CDF nnplot...")
                scan()
            }
        }
        if (nnplot == "all" || nnplot == "seq" || nnplot == "cdf") {
            qq <- seq(0,1,len=51)
            cdf <- qssden(ltdfit, qq)
            plot(cdf, qq, ann=FALSE, type="l", lty=1, col="blue")
            par(lty=2, col="red")
            abline(v=x$awbdif)
            par(lty=3, col="green3")
            abline(v=x$wwbdif)
            par(lty=1, col="black")
            abline(v=0)            
            title(main=paste("NN/LTD Weighted CDF for", length(na.omit(x$pbinsde)),
               "Informative Clusters"), xlab = "Local Treatment Difference (LTD)",
                ylab = "gss Probability Less Than")
        }
    }
}

"print.SPSbalan" <-
function (x, ...) 
{
    cat("\nSPSbalan: Checking for Within-Bin BALANCE on X Covariates\n")
    cat("Data Frame input:", x$dframe, "\n")
    cat("Treatment factor:", x$trtm, "\n")
    cat("Number of PS Bins:", x$bins, "\n")
    cat("Covariate X variable:", x$xvar, "\n")
    cat("\nTest for Overall Covariate Differences between Treatment Groups:\n")
    print(x$form)
    if (x$youtype == "contin") {
        print(summary(x$aovdiff))
        if (x$bins > 1) {
            cat("\nTest for Covariate Differences between Treatments within Bins:\n")
            print(x$form2)
            print(summary(x$bindiff))
        }
    }
    else {
        print(x$tab)
        print(chisq.test(x$tab))
        cat("\nTest for Covariate Differences between Treatments within Bins")
        cat("\n    Cumulative ChiSquare = ", round(x$cumchi, 
            digits = 2))
        cat("\n    Cumulative df = ", x$cumdf)
        cat("\n    Cumulative p.value = ", round(1 - pchisq(x$cumchi, 
            x$cumdf), digits = 4), "\n\n")
    }
}

"print.SPSloess" <-
function (x, ...) 
{
    cat("\n\nSPSloess Object: Supervised PS LOESS Smoothing Analysis\n")
    cat("Data Frame input:", x$dframe, "\n")
    cat("Treatment factor:", x$trtm, "\n")
    cat("Propensity Score variable:", x$pscr, "\n")
    cat("Outcome variable:", x$yvar, "\n")
    cat("\nOverall Average Outcome Means by Treatment")
    print(x$rawmean)
    cat("\nOverall Treatment Frequencies")
    print(x$rawfreq)
    cat("\nRaw Average Outcome Difference =", x$rawmean[2, ] - 
        x$rawmean[1, ])
    cat("\nStandard Deviation of Difference =", sqrt(x$rawvars[2, 
        ] + x$rawvars[1, ]))
    cat("\nTest of Raw Treatment Difference:\n")
    print(summary(x$aovdiff))
    cat("\n\nLOESS Smoothing Span =", x$span)
    cat("\nMean Spline Weighted Difference =", x$lotdif)
    cat("\nStandard Error of Spline Difference =", x$lotsde, 
        "\n\n")
}

"print.SPSlogit" <-
function (x, ...) 
{
    cat("\nSPSlogit Object: Supervised Estimation of Propensity Scores\n")
    cat("Data Frame  input:", x$dfname, "\n")
    cat("Data Frame output:", x$dfoutnam, "\n")
    cat("Treatment Factor:", x$trtm, "\n")
    cat("Logistic Regression Formula:\n")
    print(x$form)
    cat("\n    Variable containing fitted PScores =", x$pfit)
    cat("\n    Variable containing PS Ranks =", x$prnk)
    cat("\n    PScore Bin Number Factor =", x$qbin)
    cat("\n    Number of PS Bins =", x$bins, "\n")
    print(summary(x$glmobj))
}

"print.SPSoutco" <-
function (x, ...) 
{
    cat("\nSPSoutco Object: Within Bin Outcome LTD\n")
    cat("Data Frame input:", x$dframe, "\n")
    cat("Outcome Y variable:", x$yvar, "\n")
    cat("Treatment difference:", x$PStdif, "\n")
    cat("Number of PS Bins:", x$bins, "\n")
    cat("Number of Informative PS Bins =", length(na.omit(x$pbinsde)), 
        "\n")
    cat("\n    Overall Raw Average Treatment Difference =", x$ratdif)
    cat("\n    Standard Deviation of this Raw Difference =", 
        x$ratsde)
    cat("\n    Average Within Bin Treatment Difference =", x$awbdif)
    cat("\n    Standard Deviation of Average Within Bin Difference =", 
        x$awbsde)
    cat("\n    Inverse Variance Weighted Difference =", x$wwbdif)
    cat("\n    Standard Deviation of this Weighted Difference =", 
        x$wwbsde)
    cat("\n\nTest for Raw / Unadjusted Treatment Difference:\n")
    print(x$form)
    if (x$youtype == "contin") {
        print(summary(x$aovdiff))
        if (x$bins > 1) {
            cat("\nTest for Treatment Difference within Bins:\n")
            print(x$form2)
            print(summary(x$bindiff))
        }
    }
    else {
        print(x$tab)
        print(chisq.test(x$tab))
        cat("\nTest for Treatment Difference within Bins")
        cat("\n    Cumulative ChiSquare = ", round(x$cumchi, 
            digits = 2))
        cat("\n    Cumulative df = ", x$cumdf)
        cat("\n    Cumulative p.value = ", round(1 - pchisq(x$cumchi, 
            x$cumdf), digits = 4), "\n\n")
    }
}

"print.SPSsmoot" <-
function (x, ...) 
{
    cat("\n\nSPSsmoot Object: Supervised PS Smoothing Spline Analysis\n")
    cat("Data Frame input:", x$dframe, "\n")
    cat("Treatment factor:", x$trtm, "\n")
    cat("Propensity Score variable:", x$pscr, "\n")
    cat("Outcome variable:", x$yvar, "\n")
    cat("\nOverall Average Outcome Means by Treatment")
    print(x$rawmean)
    cat("\nOverall Treatment Frequencies")
    print(x$rawfreq)
    cat("\nRaw Average Outcome Difference =", x$rawmean[2, ] - 
        x$rawmean[1, ])
    cat("\nStandard Deviation of Difference =", sqrt(x$rawvars[2, 
        ] + x$rawvars[1, ]))
    cat("\nTest of Raw Treatment Difference:\n")
    print(summary(x$aovdiff))
    cat("\n\nSmoothing Spline Degrees-of-Freedom =", x$df)
    cat("\nMean Spline Weighted Difference =", x$sptdif)
    cat("\nStandard Error of Spline Difference =", x$sptsde, 
        "\n\n")
}

"print.UPSaltdd" <-
function(x, ...)
{
    cat("\nUPSaltdd Object: Artificial Distribution of LTDs for random clusters...\n")
    cat("Data Frame:", x$dframe, "\n")
    cat("Outcome Variable:", x$yvar, "\n")
    cat("Treatment Factor:", x$trtm, "\n")
    cat("Scedasticity assumption:", x$scedas, "\n")
    cat("Number of Replications:", x$reps, "\n")
    cat("Number of Clusters per Replication:", x$clus, "\n")
    cat("Total Number of Informative Clusters =", length(x$altdd[,1]), "\n" )
    cat("\n    Mean Artificial Treatment Difference =", mean(x$altdd[,1]))
    cat("\n    Number of Smoothed Sample Quantiles  =", length(x$altdcdf))
    cat("\n    Mean of Smoothed Sample Quantiles    =", mean(x$altdcdf))
    cat("\n    Std. Deviation of Sample Quantiles   =", var(x$altdcdf)^0.5, "\n")
    if( x$NNobj!="NA" ) {
        cat("\nUPSnnltd Object:", x$NNobj, "\n")
        cat("Number of Informative Clusters =", length(x$nnltdd[,1]), "\n" )
        cat("\n    Mean of Observed LTD Distribution    =", mean(x$nnltdd[,1]))
        cat("\n    Number of Smoothed Sample Quantiles  =", length(x$nnltdcdf))
        cat("\n    Mean of Smoothed Sample Quantiles    =", mean(x$nnltdcdf))
        cat("\n    Std. Deviation of Sample Quantiles   =", var(x$nnltdcdf)^0.5, "\n")
        }
    }

"print.UPShclus" <-
function (x, ...) 
{
    cat("\n\nUPShclus object: Unsupervised Hierarchical Clustering\n")
    cat("\nData Frame input:", x$dframe)
    cat("\nClustering algorithm used:", x$method)
    cat("\nCovariate X variables:")
    print(x$xvars, quote = FALSE)
}

"print.UPSivadj" <-
function (x, ...) 
{
    cat("\nUPSivadj Object: Clustering Instrumental Variable (IV) Adjustment\n")
    cat("Hierarchical Clustering object:", x$hiclus, "\n")
    cat("Data Frame input:", x$dframe, "\n")
    cat("Treatment variable:", x$trtm, "\n")
    cat("Outcome variable:", x$yvar, "\n")
    cat("Maximum Number of Clusters:", x$numclust, "\n")
    if (x$actclust != x$numclust) 
        cat("Actual Number of Clusters =", x$actclust, "\n")
    cat("\n    Overall Raw Average Outcome by Treatment\n")
    print(x$rawmean)
    cat("\n    Overall Treatment Frequency\n")
    print(x$rawfreq)
    cat("\n    Predicted Outcome at Treament Percentage Zero =", 
        x$ivtzero)
    cat("\n    Standard Error at Treament Percentage Zero =", 
        x$ivtzsde)
    cat("\n    Predicted Outcome at Treament Percentage 100 = ", 
        x$ivt100p)
    cat("\n    Standard Error at Treament Percentage 100 =", 
        x$ivt1pse)
    cat("\n    Predicted Overall Treatment Outcome Difference =", 
        x$ivtdiff)
    cat("\n    Standard Error of Overall Outcome Difference =", 
        x$ivtdsde)
    cat("\n\nWithin Cluster Average Outcomes\n")
    print(x$binmean)
    cat("\nWithin Cluster Treatment Frequencies\n")
    print(x$binfreq)
    if (x$youtype == "contin" && x$actclust > 1) {
        cat("\nWithin Cluster Treatment Percentages\n")
        print(x$pbinpsp)
        cat("\nCluster Symbol Radii (Root Total Sizes)\n")
        print(x$pbinsiz)
    }
}

"print.UPSnnltd" <-
function (x, ...) 
{
    cat("\nUPSnnltd Object: Nearest Neighbor Local Treatment Differences\n")
    cat("Hierarchical Clustering object:", x$hiclus, "\n")
    cat("Data Frame input:", x$dframe, "\n")
    cat("Outcome variable:", x$yvar, "\n")
    cat("Treatment difference:", x$PStdif, "\n")
    cat("Scedasticity assumption:", x$scedas, "\n")
    if (length(x$sig2) == 1) 
        cat("Homoscedastic Sigma:", sqrt(x$sig2), "\n")
    cat("Maximum Number of Clusters:", x$numclust, "\n")
    if (x$actclust != x$numclust) 
        cat("Actual Number of Clusters =", x$actclust, "\n")
    cat("Number of Informative Clusters =", length(na.omit(x$pbinsde)), 
        "\n")
    cat("\n    Overall Raw Average Treatment Difference =", x$ratdif)
    cat("\n    Standard Deviation of this Raw Difference =", 
        x$ratsde)
    cat("\n    Average Within Cluster Treatment Difference =", 
        x$awbdif)
    cat("\n    Standard Deviation of Average Within Cluster Difference =", 
        x$awbsde)
    cat("\n    Inverse Variance Weighted Difference =", x$wwbdif)
    cat("\n    Standard Deviation of this Weighted Difference =", 
        x$wwbsde)
    cat("\n\nTest for Raw / Unadjusted Treatment Difference:\n")
    print(x$form)
    if (x$youtype == "contin") {
        print(x$aovdiff)
        if (x$actclust > 1) {
            cat("\nTest for Treatment Difference within Clusters:\n")
            print(x$form2)
            print(x$bindiff)
        }
    }
    else {
        print(x$tab)
        print(chisq.test(x$tab))
        cat("\nTest for Treatment Difference within Clusters")
        cat("\n    Cumulative ChiSquare = ", round(x$cumchi, 
            digits = 2))
        cat("\n    Cumulative df = ", x$cumdf)
        cat("\n    Cumulative p.value = ", round(1 - pchisq(x$cumchi, 
            x$cumdf), digits = 4), "\n\n")
    }
}

"SPSbalan" <-
function (dframe, trtm, qbin, xvar, faclev = 3) 
{
    if (missing(dframe) || !inherits(dframe, "data.frame")) 
        stop("First argument to SPSbalan must be an existing Data Frame.")
    if (missing(trtm)) 
        stop("Second argument to SPSbalan must name the Treatment factor.")
    trtm <- deparse(substitute(trtm))
    if (!is.element(trtm, dimnames(dframe)[[2]])) 
        stop("Treatment factor must be an existing Data Frame variable.")
    if (missing(qbin)) 
        stop("Third argument to SPSbalan must be the PS Bin Number variable.")
    qbin <- deparse(substitute(qbin))
    if (!is.element(qbin, dimnames(dframe)[[2]])) 
        stop("PS Bin number factor must be an existing Data Frame variable.")
    if (missing(xvar)) 
        stop("Fourth argument to SPSbalan must name a X predictor variable.")
    xvar <- deparse(substitute(xvar))
    if (!is.element(xvar, dimnames(dframe)[[2]])) 
        stop("Predictor X variable must be an existing Data Frame variable.")
    bins <- length(table(dframe[, qbin]))
    SPSolist <- list(dframe = deparse(substitute(dframe)), trtm = trtm, 
        qbin = qbin, bins = bins, xvar = xvar)
    form <- as.formula(paste(xvar, "~", trtm))
    SPSolist <- c(SPSolist, list(form = form))
    if (length(table(dframe[, xvar])) > faclev) {
        youtype <- "contin"
        aovdiff <- aov(form, dframe, na.action = na.omit)
        SPSolist <- c(SPSolist, list(youtype = youtype, faclev = faclev, 
            aovdiff = aovdiff))
        if (bins > 1) {
            form <- as.formula(paste(xvar, "~", qbin, "+", trtm, 
                "%in%", qbin))
            aovdiff <- aov(form, dframe, na.action = na.omit)
            df3 <- as.data.frame(cbind(dframe[, xvar], as.numeric(as.character(dframe[, 
                trtm])), dframe[, qbin]))
            df3[, 2] <- as.factor(df3[, 2])
            df3[, 3] <- as.factor(df3[, 3])
            names(df3) <- c(xvar, trtm, qbin)
            SPSolist <- c(SPSolist, list(form2 = form, bindiff = aovdiff, 
                df3 = df3))
        }
    }
    else {
        youtype <- "factor"
        df3 <- na.omit(as.data.frame(cbind(as.factor(dframe[, 
            xvar]), dframe[, trtm], dframe[, qbin])))
        names(df3) <- c(xvar, trtm, qbin)
        tab <- table(df3[, 1], df3[, 2])
        SPSolist <- c(SPSolist, list(youtype = youtype, faclev = faclev, 
            factab = tab))
        tab <- table(df3[, 1], df3[, 2], df3[, 3])
        cumchi <- 0
        cumdf <- 0
        for (i in 1:bins) {
            ht <- chisq.test(tab[, , i])
            if (ht$statistic < Inf) {
                cumchi <- cumchi + ht$statistic
                cumdf <- cumdf + ht$parameters
            }
        }
        SPSolist <- c(SPSolist, list(tab = tab, cumchi = cumchi, 
            cumdf = cumdf))
    }
    class(SPSolist) <- "SPSbalan"
    SPSolist
}

"SPSloess" <-
function (dframe, trtm, pscr, yvar, faclev = 3, deg = 2, span = 0.75, 
    fam = "symmetric") 
{
    if (missing(dframe) || !inherits(dframe, "data.frame")) 
        stop("First argument to SPSloess must be an existing Data Frame.")
    if (missing(trtm)) 
        stop("Second argument to SPSloess must name the Treatment factor.")
    trtm <- deparse(substitute(trtm))
    if (!is.element(trtm, dimnames(dframe)[[2]])) 
        stop("Treatment factor must be an existing Data Frame variable.")
    if (missing(pscr)) 
        stop("Third argument to SPSlowes must name the PScore variable.")
    pscr <- deparse(substitute(pscr))
    if (!is.element(pscr, dimnames(dframe)[[2]])) 
        stop("Propensity Score must be an existing Data Frame variable.")
    if (missing(yvar)) 
        stop("Fourth argument to SPSloess must name the Target Outcome.")
    yvar <- deparse(substitute(yvar))
    if (!is.element(yvar, dimnames(dframe)[[2]])) 
        stop("Target Outcome must be an existing Data Frame variable.")
    if (length(table(dframe[, yvar])) <= faclev) 
        stop("SPSloess is intended for use only with continuous Outcomes.")
    SPSolist <- list(dframe = deparse(substitute(dframe)), trtm = trtm, 
        pscr = pscr, yvar = yvar)
    PSmean <- as.matrix(tapply(dframe[, yvar], dframe[, trtm], 
        na.rm = TRUE, mean))
    PStrtm <- paste(trtm, "=", dimnames(PSmean)[[1]])
    PSvars <- as.matrix(tapply(dframe[, yvar], dframe[, trtm], 
        na.rm = TRUE, var))
    PSfreq <- as.matrix(table(dframe[, trtm]))
    PSvars <- PSvars/PSfreq
    form <- as.formula(paste(yvar, "~", trtm))
    aovdiff <- aov(form, dframe, na.action = na.omit)
    SPSolist <- c(SPSolist, list(rawmean = PSmean, rawvars = PSvars, 
        rawfreq = PSfreq, form = form, aovdiff = aovdiff))
    m <- order(dframe[, pscr])
    lofit <- as.data.frame(cbind(dframe[m, pscr], dframe[m, yvar], 
        as.numeric(as.character(dframe[m, trtm])), 0))
    lofit <- na.omit(lofit)
    names(lofit) <- c("PS", "YVAR", "TRTM", "FIT")
    logrid <- as.data.frame(cbind(seq(0.005, 0.995, length = 100), 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    names(logrid) <- c("PS", "F0", "S0", "C0", "F1", "S1", "C1", 
        "DIF", "SED", "HST", "DEN")
    losub0 <- subset.data.frame(lofit, TRTM == 0)
    losub1 <- subset.data.frame(lofit, TRTM == 1)
    loobj0 <- loess(YVAR ~ PS, losub0, family = fam, degree = deg, 
        span = span)
    losub0$FIT <- fitted.values(loobj0)
    fit0 <- predict(loobj0, logrid$PS, se = TRUE)
    logrid$F0 <- fit0$fit
    logrid$S0 <- fit0$se.fit
    logrid$C0 <- hist(losub0$PS, breaks = seq(0, 1, length = 101), 
        plot = FALSE)$counts
    loobj1 <- loess(YVAR ~ PS, losub1, family = fam, degree = deg, 
        span = span)
    losub1$FIT <- fitted.values(loobj1)
    fit1 <- predict(loobj1, logrid$PS, se = TRUE)
    logrid$F1 <- fit1$fit
    logrid$S1 <- fit1$se.fit
    logrid$C1 <- hist(losub1$PS, breaks = seq(0, 1, length = 101), 
        plot = FALSE)$counts
    logrid$DIF <- logrid$F1 - logrid$F0
    logrid$SED <- sqrt(logrid$S0^2 + logrid$S1^2)
    logrid$HST <- hist(lofit$PS, breaks = seq(0, 1, length = 101), 
        plot = FALSE)$counts
    logrid$DEN <- density(lofit$PS, bw = "nrd0", kernel = "gaussian", 
        n = 100, from = 0.005, to = 0.995)$y
    grid <- na.omit(logrid)
    if (sum(grid$DEN) > 0) 
        grid$DEN <- grid$DEN/sum(grid$DEN)
    PSmean <- sum(grid$DIF * grid$DEN)
    PSster <- sum(grid$SED * grid$DEN)
    SPSolist <- c(SPSolist, list(logrid = logrid, losub0 = losub0, 
        losub1 = losub1, span = span, lotdif = PSmean, lotsde = PSster))
    class(SPSolist) <- "SPSloess"
    SPSolist
}

"SPSlogit" <-
function (dframe, form, pfit, prnk, qbin, bins = 5, appn = "") 
{
    if (missing(dframe) || !inherits(dframe, "data.frame")) 
        stop("First argument to SPSlogit must be an existing Data Frame.")
    dfname <- deparse(substitute(dframe))
    if (missing(form) || class(form) != "formula") 
        stop("Second argument to SPSlogit must be a valid formula.")
    trtm <- deparse(form[[2]])
    if (!is.element(trtm, dimnames(dframe)[[2]])) 
        stop("Response variable in the SPSlogit formula must be an existing treatment factor.")
    dframe[, trtm] <- as.factor(dframe[, trtm])
    if (missing(pfit)) 
        stop("Third argument to SPSlogit must name the PScore variable.")
    if (missing(prnk)) 
        stop("Fourth argument to SPSlogit must name the PS Rank variable.")
    if (missing(qbin)) 
        stop("Fifth argument to SPSlogit must name the Bin Number variable.")
    glmobj <- glm(form, family = binomial(), data = dframe)
    df3 <- as.data.frame(fitted.values(glmobj))
    pfit <- deparse(substitute(pfit))
    names(df3) <- pfit
    prnk <- deparse(substitute(prnk))
    df3[, prnk] <- rank(df3[, pfit], na.last = TRUE)
    qbin <- deparse(substitute(qbin))
    df3[, qbin] <- factor(1 + floor((bins * df3[, prnk])/(1 + 
        length(df3[, prnk]))))
    dframe <- merge(dframe, df3, by.x = "row.names", by.y = "row.names", 
        all.x = TRUE)
    if (appn == "") 
        dfoutnam <- dfname
    else dfoutnam <- appn
    assign(dfoutnam, dframe, env = .GlobalEnv)
    SPSolist <- list(dfname = dfname, dfoutnam = dfoutnam, trtm = trtm, 
        form = form, pfit = pfit, prnk = prnk, qbin = qbin, bins = bins, 
        glmobj = glmobj)
    class(SPSolist) <- "SPSlogit"
    SPSolist
}

"SPSnbins" <-
function (dframe, prnk, qbin, bins = 8) 
{
    if (missing(dframe) || !inherits(dframe, "data.frame")) 
        stop("First argument to SPSnbins must be an existing Data Frame.")
    prnk <- deparse(substitute(prnk))
    if (missing(prnk)) 
        stop("Second argument to SPSnbins must name the patient rank variable.")
    if (!is.element(prnk, dimnames(dframe)[[2]])) 
        stop("SPSnbins patient ranks must be stored in an existing variable.")
    qbin <- deparse(substitute(qbin))
    if (missing(qbin)) 
        stop("Third argument to SPSnbins must name the bin number variable.")
    subjects <- length(dframe[, prnk])
    if (bins < 2) 
        bins <- 2
    if (bins > floor(subjects/2)) 
        bins <- floor(subjects/2)
    dframe[, qbin] <- factor(1 + floor((bins * dframe[, prnk])/(1 + 
        subjects)))
    dframe
}

"SPSoutco" <-
function (dframe, trtm, qbin, yvar, faclev = 3) 
{
    if (missing(dframe) || !inherits(dframe, "data.frame")) 
        stop("First argument to SPSoutco must be an existing Data Frame.")
    if (missing(trtm)) 
        stop("Second argument to SPSoutco must name the Treatment factor.")
    trtm <- deparse(substitute(trtm))
    if (!is.element(trtm, dimnames(dframe)[[2]])) 
        stop("Treatment factor must be an existing Data Frame variable.")
    if (missing(qbin)) 
        stop("Third argument to SPSoutco must be the PS Bin Number variable.")
    qbin <- deparse(substitute(qbin))
    if (!is.element(qbin, dimnames(dframe)[[2]])) 
        stop("PS Bin number must be an existing Data Frame variable.")
    if (missing(yvar)) 
        stop("Fourth argument to SPSoutco must name the Outcome variable.")
    yvar <- deparse(substitute(yvar))
    if (!is.element(yvar, dimnames(dframe)[[2]])) 
        stop("Target Outcome or Covariate must be an existing Data Frame variable.")
    bins <- length(table(dframe[, qbin]))
    SPSolist <- list(dframe = deparse(substitute(dframe)), trtm = trtm, 
        yvar = yvar, bins = bins)
    PSmean <- as.matrix(tapply(dframe[, yvar], dframe[, trtm], 
        na.rm = TRUE, mean))
    PStrtm <- paste(trtm, "=", dimnames(PSmean)[[1]])
    PStdif <- paste(PStrtm[2], "minus", PStrtm[1])
    PSvars <- as.matrix(tapply(dframe[, yvar], dframe[, trtm], 
        na.rm = TRUE, var))
    PSfreq <- as.matrix(table(dframe[, trtm]))
    PSvars <- PSvars/PSfreq
    SPSolist <- c(SPSolist, list(PStdif = PStdif, rawmean = PSmean, 
        rawvars = PSvars, rawfreq = PSfreq))
    RATdif <- sum(PSmean[2, ] - PSmean[1, ])
    RATsde <- sum(sqrt(PSvars[2, ] + PSvars[1, ]))
    SPSolist <- c(SPSolist, list(ratdif = RATdif, ratsde = RATsde))
    PSmean <- tapply(dframe[, yvar], list(dframe[, qbin], dframe[, 
        trtm]), na.rm = TRUE, mean)
    PSmean <- cbind(matrix(1:bins, bins, 1), PSmean)
    dimnames(PSmean) <- list(1:bins, c("BIN", PStrtm[1], PStrtm[2]))
    PSmean <- as.data.frame(PSmean)
    PSmean[, "BIN"] <- as.factor(PSmean[, "BIN"])
    PSvars <- as.matrix(tapply(dframe[, yvar], list(dframe[, 
        qbin], dframe[, trtm]), na.rm = TRUE, var))
    PSfreq <- as.matrix(table(dframe[, qbin], dframe[, trtm]))
    PSvars <- PSvars/PSfreq
    PSfreq <- cbind(matrix(1:bins, bins, 1), PSfreq)
    dimnames(PSfreq) <- list(1:bins, c("BIN", PStrtm[1], PStrtm[2]))
    PSfreq <- as.data.frame(PSfreq)
    PSfreq[, "BIN"] <- as.factor(PSfreq[, "BIN"])
    SPSolist <- c(SPSolist, list(binmean = PSmean, binvars = PSvars, 
        binfreq = PSfreq))
    awbdif <- sum(na.omit(((PSmean[, 3] - PSmean[, 2]) * (PSfreq[, 
        3] + PSfreq[, 2]))/sum(PSfreq[, 3] + PSfreq[, 2])))
    awbsde <- sqrt(sum(na.omit((PSvars[, 2] + PSvars[, 1]) * 
        (PSfreq[, 3] + PSfreq[, 2])^2))/(sum(PSfreq[, 3] + PSfreq[, 
        2]))^2)
    wwbdif <- sum(na.omit((PSmean[, 3] - PSmean[, 2])/(PSvars[, 
        2] + PSvars[, 1])/sum(1/na.omit(PSvars[, 2] + PSvars[, 
        1]))))
    wwbsde <- sqrt(1/sum(1/na.omit(PSvars[, 2] + PSvars[, 1])))
    form <- as.formula(paste(yvar, "~", trtm))
    SPSolist <- c(SPSolist, list(awbdif = awbdif, awbsde = awbsde, 
        wwbdif = wwbdif, wwbsde = wwbsde, form = form))
    if (length(table(dframe[, yvar])) > faclev) {
        youtype <- "contin"
        aovdiff <- aov(form, dframe, na.action = na.omit)
        SPSolist <- c(SPSolist, list(youtype = youtype, faclev = faclev, 
            aovdiff = aovdiff))
        if (bins > 1) {
            form <- as.formula(paste(yvar, "~", qbin, "+", trtm, 
                "%in%", qbin))
            aovdiff <- aov(form, dframe, na.action = na.omit)
            SPSolist <- c(SPSolist, list(form2 = form, bindiff = aovdiff))
        }
        pbindif <- PSmean[, 3] - PSmean[, 2]
        pbinsde <- sqrt(PSvars[, 2] + PSvars[, 1])
        pbinsiz <- sqrt(PSfreq[, 3] + PSfreq[, 2])
        SPSolist <- c(SPSolist, list(pbindif = pbindif, pbinsde = pbinsde, 
            pbinsiz = pbinsiz))
    }
    else {
        youtype <- "factor"
        df3 <- as.data.frame(cbind(dframe[, yvar], dframe[, trtm]))
        df3[, 3] <- dframe[, qbin]
        df3 <- na.omit(df3)
        names(df3) <- c(yvar, trtm, qbin)
        df3[, 1] <- as.factor(df3[, 1])
        tab <- table(df3[, 1], df3[, 2])
        SPSolist <- c(SPSolist, list(youtype = youtype, faclev = faclev, 
            factab = tab))
        tab <- table(df3[, 1], df3[, 2], df3[, 3])
        cumchi <- 0
        cumdf <- 0
        for (i in 1:bins) {
            ht <- chisq.test(tab[, , i])
            if (!is.na(ht$statistic) && ht$statistic < Inf) {
                cumchi <- cumchi + ht$statistic
                cumdf <- cumdf + ht$parameters
            }
        }
        SPSolist <- c(SPSolist, list(tab = tab, cumchi = cumchi, 
            cumdf = cumdf))
    }
    class(SPSolist) <- "SPSoutco"
    SPSolist
}

"SPSsmoot" <-
function (dframe, trtm, pscr, yvar, faclev = 3, df = 5, spar = NULL, 
    cv = FALSE, penalty = 1) 
{
    if (missing(dframe) || !inherits(dframe, "data.frame")) 
        stop("First argument to SPSsmoot must be a Data Frame name.")
    if (missing(trtm)) 
        stop("Second argument to SPSsmoot must name the Treatment factor.")
    trtm <- deparse(substitute(trtm))
    if (!is.element(trtm, dimnames(dframe)[[2]])) 
        stop("Treatment factor must be an existing Data Frame variable.")
    if (missing(pscr)) 
        stop("Third argument to SPSsmoot must be fitted Propensity Scores.")
    pscr <- deparse(substitute(pscr))
    if (!is.element(pscr, dimnames(dframe)[[2]])) 
        stop("Propensity Score must be an existing Data Frame variable.")
    if (missing(yvar)) 
        stop("Fourth argument to SPSsmoot must name the Target Outcome.")
    yvar <- deparse(substitute(yvar))
    if (!is.element(yvar, dimnames(dframe)[[2]])) 
        stop("Target Outcome must be an existing Data Frame variable.")
    if (length(table(dframe[, yvar])) <= faclev) 
        stop("SPSsmoot is intended for use only with continuous Outcomes.")
    SPSolist <- list(dframe = deparse(substitute(dframe)), trtm = trtm, 
        pscr = pscr, yvar = yvar)
    PSmean <- as.matrix(tapply(dframe[, yvar], dframe[, trtm], 
        na.rm = TRUE, mean))
    PSvars <- as.matrix(tapply(dframe[, yvar], dframe[, trtm], 
        na.rm = TRUE, var))
    PSfreq <- as.matrix(table(dframe[, trtm]))
    PSvars <- PSvars/PSfreq
    form <- as.formula(paste(yvar, "~", trtm))
    aovdiff <- aov(form, dframe, na.action = na.omit)
    SPSolist <- c(SPSolist, list(rawmean = PSmean, rawvars = PSvars, 
        rawfreq = PSfreq, form = form, aovdiff = aovdiff))
    ssfit <- as.data.frame(cbind(dframe[, pscr], dframe[, yvar], 
        as.numeric(as.character(dframe[, trtm])), 0, 0))
    ssfit <- na.omit(ssfit)
    names(ssfit) <- c("PS", "YVAR", "TRTM", "FIT", "SEP")
    ssgrid <- as.data.frame(cbind(seq(0.005, 0.995, length = 100), 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    names(ssgrid) <- c("PS", "F0", "S0", "C0", "F1", "S1", "C1", 
        "DIF", "SED", "HST", "DEN")
    spsub0 <- subset.data.frame(ssfit, TRTM == 0)
    spsub1 <- subset.data.frame(ssfit, TRTM == 1)
    ssobj0 <- smooth.spline(spsub0$PS, spsub0$YVAR, df = df, 
        spar = spar, cv = cv, penalty = penalty)
    res <- (ssobj0$yin - ssobj0$y)/(1 - ssobj0$lev)
    sigma <- sqrt(var(res))
    spsub0 <- as.data.frame(cbind(ssobj0$x, ssobj0$yin, 0, ssobj0$y, 
        sigma * sqrt(ssobj0$lev)))
    names(spsub0) <- c("PS", "YAVG", "TRTM", "FIT", "SEP")
    ssgrid$F0 <- predict(ssobj0, ssgrid$PS)$y
    ssobj0 <- smooth.spline(spsub0$PS, spsub0$SEP, df = df, spar = spar, 
        cv = cv, penalty = penalty)
    ssgrid$S0 <- predict(ssobj0, ssgrid$PS)$y
    ssobj1 <- smooth.spline(spsub1$PS, spsub1$YVAR, df = df, 
        spar = spar, cv = cv, penalty = penalty)
    res <- (ssobj1$yin - ssobj1$y)/(1 - ssobj1$lev)
    sigma <- sqrt(var(res))
    spsub1 <- as.data.frame(cbind(ssobj1$x, ssobj1$yin, 1, ssobj1$y, 
        sigma * sqrt(ssobj1$lev)))
    names(spsub1) <- c("PS", "YAVG", "TRTM", "FIT", "SEP")
    ssgrid$F1 <- predict(ssobj1, ssgrid$PS)$y
    ssobj1 <- smooth.spline(spsub1$PS, spsub1$SEP, df = df, spar = spar, 
        cv = cv, penalty = penalty)
    ssgrid$S1 <- predict(ssobj1, ssgrid$PS)$y
    ssgrid$DIF <- ssgrid$F1 - ssgrid$F0
    ssgrid$SED <- sqrt(ssgrid$S0^2 + ssgrid$S1^2)
    ssgrid$C0 <- hist(spsub0$PS, breaks = seq(0, 1, length = 101), 
        plot = FALSE, probability = FALSE)$counts
    ssgrid$C1 <- hist(spsub1$PS, breaks = seq(0, 1, length = 101), 
        plot = FALSE, probability = FALSE)$counts
    ssgrid$HST <- hist(ssfit$PS, breaks = seq(0, 1, length = 101), 
        plot = FALSE, probability = TRUE)$counts
    ssgrid$DEN <- density(ssfit$PS, bw = "nrd0", kernel = "gaussian", 
        n = 100, from = 0.005, to = 0.995)$y
    grid <- na.omit(ssgrid)
    if (sum(grid$DEN) > 0) 
        grid$DEN <- grid$DEN/sum(grid$DEN)
    PSmean <- sum(grid$DIF * grid$DEN)
    PSster <- sum(grid$SED * grid$DEN)
    SPSolist <- c(SPSolist, list(ssgrid = ssgrid, spsub0 = spsub0, 
        spsub1 = spsub1, df = df, sptdif = PSmean, sptsde = PSster))
    class(SPSolist) <- "SPSsmoot"
    SPSolist
}

"UPlinint" <-
function(q, xmin, n, x, w)
{
    # linear interpolation within a CDF...
    
    if( q < 0 || q > 1 )
        stop("Desired quantiles must be between 0 and 1, inclusive.")
    cdf <- xmin
    while( TRUE ) {
        if( q == 0 ) break
        if( q == 1 ) {
            cdf <- x[n]
            break
            }
        j <- 0
        qlo <- 0
        qhi <- w[1] + 0.0000001
        while( q > qhi ) {
             j <- j+1
             if( j > (n-1) ) {
                 j <- j-1
                 break
                 }
             qlo <- qhi
             qhi <- qhi + w[(j+1)]
             }
        if( j== 0 )
            cdf <- xmin + (q-qlo)*(x[1]-xmin)/(qhi-qlo)
        else
            cdf <- x[j] + (q-qlo)*(x[(j+1)]-x[j])/(qhi-qlo)
        break
        }    
    cdf
}

"UPSaccum" <-
function (hiclus, dframe, trtm, yvar, faclev = 3, scedas = "hete", 
    accobj = "UPSframe") 
{
    if (missing(hiclus) || (!inherits(hiclus$upshcl, "diana") && 
        !inherits(hiclus$upshcl, "agnes") && !inherits(hiclus$upshcl, 
        "hclust"))) 
        stop("First argument to UPSaccum must be a diana, agnes or hclust object.")
    hiclus <- deparse(substitute(hiclus))
    if (missing(dframe) || !inherits(dframe, "data.frame")) 
        stop("Second argument to UPSaccum must be an existing Data Frame.")
    if (missing(trtm)) 
        stop("Third argument to UPSaccum must name the Treatment factor.")
    trtm <- deparse(substitute(trtm))
    if (!is.element(trtm, dimnames(dframe)[[2]])) 
        stop("Treatment factor must be an existing Data Frame variable.")
    if (length(table(dframe[, trtm])) != 2) 
        stop("Treatment factor must assume exactly two different levels.")
    if (missing(yvar)) 
        stop("Fourth argument to UPSaccum must name the Target variable.")
    yvar <- deparse(substitute(yvar))
    if (!is.element(yvar, dimnames(dframe)[[2]])) 
        stop("Target Outcome or Covariate must be an existing Data Frame variable.")
    dframe <- deparse(substitute(dframe))
    if (scedas != "homo") 
        scedas <- "hete"
    nnymax <- 0
    nnxmin <- 0
    nnxmax <- 0
    UPSaccum.pars <<- cbind(hiclus, dframe, trtm, yvar, faclev, 
        scedas, accobj, nnymax, nnxmin, nnxmax)
    accdf <- as.data.frame(cbind("NONE", hiclus, dframe, trtm, 
        yvar, 1, 0, 0))
    names(accdf) <- c("NNIV", "hicl", "dfrm", "trtm", "yvar", 
        "bins", "tdif", "tdse")
    assign(accobj, accdf, env = .GlobalEnv)
}

"UPSaltdd" <-
function(dframe, trtm, yvar, faclev=3, scedas="homo", NNobj=NA, clus=50, reps=10, seed=12345)
{
    # Compute the Artificial LTD Distribution for random patient clusterings  ...i.e. include
    # the least or less relevant comparisons as well as those neurtal, more or most relevant.
    
    if(missing(dframe)||!inherits(dframe,"data.frame"))
        stop("First argument to UPSaltdd must be an existing Data Frame.")
    if(missing(trtm))
        stop("Second argument to UPSaltdd must name the Treatment Factor.")
    trtm <- deparse(substitute(trtm))
    if(!is.element(trtm,dimnames(dframe)[[2]]))
        stop("Treatment Factor must be present in the Data Frame.")
    if(length(table(dframe[,trtm]))!=2)
        stop("Treatment Factor must assume exactly two different levels.")
    if(missing(yvar))
        stop("Third argument to UPSaltdd must name a continuous Outcome Variable.")
    yvar <- deparse(substitute(yvar))
    if(!is.element(yvar,dimnames(dframe)[[2]]))
        stop("Outcome Variable must be present in the Data Frame.")
    if(scedas!="hete") scedas <- "homo" # variances either homoscedastic or heteroscedastic
    NNobjnam <- deparse(substitute(NNobj))
    if(NNobjnam!="NA"&&!inherits(NNobj,"UPSnnltd"))
        stop("The NNobj argument to UPSaltdd must either be NA or an existing UPSnnltd object.")
    if(NNobjnam!="NA")
       clus <- NNobj$numclust
    if(clus < 2)
       clus <- 2
    ytvars <- c(yvar,trtm)
    ytdata <- dframe[,ytvars]
    if(length(table(ytdata[,1])) <= faclev)
        stop("To be considered continuous, Outcome Variable must assume more than faclev values.")
    pats <- length(ytdata[,1])
    altdg <- "Agroups"
    ALmean <- as.matrix(tapply(ytdata[,1], ytdata[,2], na.rm = TRUE, mean))
    ALtrtm <- paste(trtm, "=", dimnames(ALmean)[[1]])
    # Start forming UPSaltdd output list...
    dframe <- deparse(substitute(dframe))

    ALolist <- list(dframe=dframe, trtm=trtm, yvar=yvar, faclev=faclev, scedas=scedas,
                    NNobj=NNobjnam, clus=clus, reps=reps, pats=pats, seed=seed)
    set.seed(seed)   # Set seed for Monte Carlo pseudo random sequence...

    for(i in 1:reps) {

       xrand <- as.data.frame(cbind(rnorm(pats),rnorm(pats)))
       crand <- kmeans(xrand, clus)
       crand <- as.data.frame(crand$cluster)

       ALmean <- as.matrix(tapply(ytdata[,1], list(crand[,1], ytdata[,2]), na.rm = TRUE, mean))
       ALmean <- cbind(matrix(1:clus, clus, 1), ALmean)
       dimnames(ALmean) <- list(1:clus, c("ALC", ALtrtm[1], ALtrtm[2]))
       ALmean <- as.data.frame(ALmean)
       ALmean[,"ALC"] <- as.factor(ALmean[,"ALC"])
       ALvars <- as.matrix(tapply(ytdata[,1], list(crand[,1], ytdata[,2]), na.rm = TRUE, var))
       ALfreq <- as.matrix(table(crand[,1], ytdata[,2]))
       ALvars <- ALvars / ALfreq  # redefined below when scedas=="homo"
       ALvars[ALvars==Inf] <- NA  

       altdif <- na.omit( ALmean[,3] - ALmean[,2] )
       if( scedas=="homo" ) {
           ALvars <- 1 / as.matrix(table(crand[,1], ytdata[,2]))
           ALvars[ALvars==Inf] <- NA
           }
       altdwt <- ALvars[,2] + ALvars[,1]
       altdwt <- na.omit( 1 / altdwt )
       altdwt <- altdwt / sum(altdwt)
       if( i == 1 ) {
           altdacc <- as.matrix(na.omit(cbind(altdif,altdwt)))
           }
       else {
           altdacc <- rbind( altdacc, as.matrix(na.omit(cbind(altdif,altdwt))))
           }
       }
    alxmin <- min(altdacc[,1])
    alxmax <- max(altdacc[,1])
    alymax <- max(altdacc[,2])
    ALolist <- c(ALolist, list(altdd=altdacc, alxmin=alxmin, alxmax=alxmax, alymax=alymax))
    
    xx <- as.vector(altdacc[,1])  
    yy <- as.vector(altdacc[,2])       
    oo <- order(xx)
    xx <- xx[oo]
    yy <- yy[oo]
    nx <- length(xx)
    # Locate and eliminate duplicate xx values...
    j <- 0
    if( nx > 1 ) {
        for(i in seq(nx, 2, by=-1) ) {
            if( xx[(i-1)]==xx[i] ) {
                yy[(i-1)] <- yy[(i-1)]+yy[i]
                xx[i] <- alxmax + 9999
                j <- j+1
                }
            }
        }
    if( j > 0 ) {
        oo <- order(xx)
        xx <- xx[oo]
        yy <- yy[oo]
        nx <- nx - j
        }
    yy <- yy/sum(yy[1:nx])
    # Locate mid-points between jumps...
    if( nx > 2 ) {
        for(i in 2:(nx-1)) {
            xx[i] <- (xx[i]+xx[(i+1)])/2
            }
        }
    qq <- seq(0,1,len=min(1001,1+round(nx*4)))
    cdf <- qq
    for(i in 1:length(qq)) {
        cdf[i] <- UPlinint(qq[i], alxmin, nx, xx, yy)
        }

    ALolist <- c(ALolist, list(altdcdf=cdf, qq=qq))
    
    if(NNobjnam!="NA") {
       nnltdif <- NNobj$pbindif
       nnltdwt <- NNobj$pbinsde
       nnltdsm <- min(mean(nnltdwt, na.rm=TRUE), min(nnltdwt[nnltdwt!=0.0]))
       nnltdwt[nnltdwt==0.0] <- nnltdsm
       nnltdwt <- 1 / nnltdwt^2
       nnltdacc <- as.matrix(na.omit(cbind(nnltdif,nnltdwt)))
       nnltdacc[,2] <- nnltdacc[,2] / sum(nnltdacc[,2])
       nnlxmin <- min(nnltdacc[,1])
       nnlxmax <- max(nnltdacc[,1])
       nnlymax <- max(nnltdacc[,2])
       ALolist <- c(ALolist, list(nnltdd=nnltdacc, nnlxmin=nnlxmin, nnlxmax=nnlxmax, nnlymax=nnlymax))
    
       xx <- as.vector(nnltdacc[,1])  
       yy <- as.vector(nnltdacc[,2])       
       oo <- order(xx)
       xx <- xx[oo]
       yy <- yy[oo]
       nx <- length(xx)
       # Locate and eliminate duplicate xx values...
       j <- 0
       if( nx > 1 ) {
          for(i in seq(nx, 2, by=-1) ) {
               if( xx[(i-1)]==xx[i] ) {
                   yy[(i-1)] <- yy[(i-1)]+yy[i]
                   xx[i] <- alxmax + 9999
                   j <- j+1
                   }
               }
           }
       if( j > 0 ) {
           oo <- order(xx)
           xx <- xx[oo]
           yy <- yy[oo]
           nx <- nx - j
           }
       yy <- yy/sum(yy[1:nx])
       # Locate mid-points between jumps...
       if( nx > 2 ) {
           for(i in 2:(nx-1)) {
               xx[i] <- (xx[i]+xx[(i+1)])/2
               }
           }
       qq <- seq(0,1,len=min(501,1+round(nx*4)))
       cdf <- qq
       for(i in 1:length(qq)) {
           cdf[i] <- UPlinint(qq[i], alxmin, nx, xx, yy)
           }

       ALolist <- c(ALolist, list(nnltdcdf=cdf, nq=qq))
       }
    
    class(ALolist) <- "UPSaltdd"
    ALolist
}

"UPSgraph" <-
function (nncol = "red", nwcol = "green3", ivcol = "blue", ...) 
{
    UPSpars <- get("UPSaccum.pars")
    hiclus <- UPSpars[1]
    dframe <- UPSpars[2]
    trtm <- UPSpars[3]
    yvar <- UPSpars[4]
    faclev <- as.numeric(UPSpars[5])
    scedas <- UPSpars[6]
    UPSdf <- get(UPSpars[7])
    if (UPSdf$NNIV[1] == "NONE") 
        UPSdf <- UPSdf[2:length(UPSdf$NNIV), ]
    UPSdf$bins <- as.numeric(as.character(UPSdf$bins))
    UPSdf$tdif <- as.numeric(as.character(UPSdf$tdif))
    UPSdf$tdse <- as.numeric(as.character(UPSdf$tdse))
    m <- order(UPSdf$hicl, UPSdf$dfrm, UPSdf$trtm, UPSdf$yvar, 
        UPSdf$NNIV, UPSdf$bins)
    UPSdf <- UPSdf[m, ]
    assign(UPSpars[7], UPSdf, env = .GlobalEnv)
    plot(UPSdf$bins, UPSdf$tdif, ylim = c(min(0, min(UPSdf$tdif - 
        2 * UPSdf$tdse)), max(0, max(UPSdf$tdif + 2 * UPSdf$tdse))), 
        log = "x", ann = FALSE, type = "n")
    symbols(UPSdf$bins, UPSdf$tdif, boxplots = cbind(0, 0, 2 * 
        UPSdf$tdse, 2 * UPSdf$tdse, 0), inches = FALSE, add = TRUE)
    x <- UPSdf[UPSdf$NNIV == "NW", ]
    points(x$bins, x$tdif, pch = 21, col = nwcol)
    x <- UPSdf[UPSdf$NNIV == "IV", ]
    points(x$bins, x$tdif, pch = 21, col = ivcol)
    x <- UPSdf[UPSdf$NNIV == "NN", ]
    points(x$bins, x$tdif, pch = 21, col = nncol)
    abline(h = UPSdf$tdif[1], lty = 2, col = nncol)
    title(main = "Unsupervised LTD Distributiun Sensitivity", 
        xlab = "Number of Clusters", ylab = "Mean LTD +/-2 Sigma LTD")
}

"UPShclus" <-
function (dframe, xvars, method = "diana") 
{
    if (missing(dframe) || !inherits(dframe, "data.frame")) 
        stop("First argument to UPShclus must be a Data Frame name.")
    if (missing(xvars)) 
        stop("Second argument to UPShclus must be a list of X variables.")
    xpc <- prcomp(dframe[, xvars], retx = TRUE, center = TRUE, 
        scale. = TRUE)
    if (method == "diana") {
        upshcl <- diana(dist(xpc$x/xpc$sdev), metric = "euclidean", 
            stand = TRUE, keep.diss = FALSE, keep.data = FALSE)
    }
    else if (method == "agnes") {
        upshcl <- agnes(dist(xpc$x/xpc$sdev), diss = FALSE, metric = "euclidean", 
            method = "complete", stand = TRUE, keep.diss = FALSE, 
            keep.data = FALSE)
    }
    else {
        upshcl <- hclust(dist(xpc$x/xpc$sdev), "complete")
    }
    dframe <- deparse(substitute(dframe))
    HCLolist <- list(dframe = dframe, xvars = xvars, method = method, 
        upshcl = upshcl)
    class(HCLolist) <- "UPShclus"
    HCLolist
}

"UPSivadj" <-
function (numclust) 
{
    if (missing(numclust)) 
        stop("The argument to UPSivadj must specify a Number of Clusters.")
    if (numclust < 3) 
        numclust <- 3
    UPSpars <- get("UPSaccum.pars")
    hiclus <- get(UPSpars[1])
    dframe <- get(UPSpars[2])
    trtm <- UPSpars[3]
    yvar <- UPSpars[4]
    faclev <- as.numeric(UPSpars[5])
    scedas <- UPSpars[6]
    UPSdf <- get(UPSpars[7])
    hclbin <- "HclusBin"
    hbins <- as.data.frame(cutree(hiclus$upshcl, k = numclust))
    names(hbins) <- hclbin
    dfivadj <- merge(dframe, hbins, by.x = "row.names", by.y = "row.names", 
        all.x = TRUE)
    dfivadj[, hclbin] <- as.factor(dfivadj[, hclbin])
    bins <- length(table(dfivadj[, hclbin]))
    IVolist <- list(hiclus = UPSpars[1], dframe = UPSpars[2], 
        trtm = trtm, yvar = yvar, numclust = numclust, actclust = bins, 
        scedas = scedas)
    PSmean <- as.matrix(tapply(dfivadj[, yvar], dfivadj[, trtm], 
        na.rm = TRUE, mean))
    PStrtm <- paste(trtm, "=", dimnames(PSmean)[[1]])
    PSfreq <- as.matrix(table(dfivadj[, trtm]))
    IVolist <- c(IVolist, list(ivhbindf = hbins, rawmean = PSmean, 
        rawfreq = PSfreq))
    PSmean <- tapply(dfivadj[, yvar], dfivadj[, hclbin], na.rm = TRUE, 
        mean)
    PSmean <- cbind(matrix(1:bins, bins, 1), PSmean)
    dimnames(PSmean) <- list(1:bins, c("BIN", "OUT"))
    PSmean <- as.data.frame(PSmean)
    PSmean[, "BIN"] <- as.factor(PSmean[, "BIN"])
    PSfreq <- as.matrix(table(dfivadj[, hclbin], dfivadj[, trtm]))
    PSfreq <- cbind(matrix(1:bins, bins, 1), PSfreq)
    dimnames(PSfreq) <- list(1:bins, c("BIN", PStrtm[1], PStrtm[2]))
    PSfreq <- as.data.frame(PSfreq)
    PSfreq[, "BIN"] <- as.factor(PSfreq[, "BIN"])
    IVolist <- c(IVolist, list(binmean = PSmean, binfreq = PSfreq))
    if (length(table(dfivadj[, yvar])) > faclev) {
        youtype <- "contin"
        IVolist <- c(IVolist, list(youtype = youtype, faclev = faclev))
        pbinout <- PSmean[, 2]
        pbinpsp <- PSfreq[, 3] + PSfreq[, 2]
        pbinpsp <- 100 * PSfreq[, 3]/pbinpsp
        pbinsiz <- sqrt(PSfreq[, 3] + PSfreq[, 2])
        symsiz <- length(dframe[, yvar])
        IVolist <- c(IVolist, list(pbinout = pbinout, pbinpsp = pbinpsp, 
            pbinsiz = pbinsiz, symsiz = symsiz))
        if (bins > 1) {
            ivfit <- lm(pbinout ~ pbinpsp, weight = pbinsiz^2)
            ivsum <- summary(ivfit)
            ivtzero <- ivsum$coefficients[[1]]
            ivtzsde <- ivsum$coefficients[[3]]
            ivtdiff <- 100 * ivsum$coefficients[[2]]
            ivtdsde <- 100 * ivsum$coefficients[[4]]
            ivsum <- predict(ivfit, data.frame(pbinpsp = 100), 
                se.fit = TRUE)
            ivt100p <- ivsum$fit
            ivt1pse <- ivsum$se.fit
            IVolist <- c(IVolist, list(ivfit = ivfit, ivtzero = ivtzero, 
                ivtzsde = ivtzsde, ivtdiff = ivtdiff, ivtdsde = ivtdsde, 
                ivt100p = ivt100p, ivt1pse = ivt1pse))
        }
    }
    else {
        youtype <- "factor"
        IVolist <- c(IVolist, list(youtype = youtype, faclev = faclev))
    }
    accnew <- as.data.frame(cbind("IV", UPSpars[1], UPSpars[2], 
        trtm, yvar, bins, ivtdiff, ivtdsde))
    names(accnew) <- c("NNIV", "hicl", "dfrm", "trtm", "yvar", 
        "bins", "tdif", "tdse")
    UPSdf <- as.data.frame(rbind(UPSdf, accnew))
    assign(UPSpars[7], UPSdf, env = .GlobalEnv)
    class(IVolist) <- "UPSivadj"
    IVolist
}

"UPSnnltd" <-
function (numclust) 
{
    if (missing(numclust)) 
        stop("The argument to UPSnnltd must specify a Number of Clusters.")
    if (numclust < 1) 
        numclust <- 1
    UPSpars <- get("UPSaccum.pars")
    hiclus <- get(UPSpars[1])
    dframe <- get(UPSpars[2])
    trtm <- UPSpars[3]
    yvar <- UPSpars[4]
    faclev <- as.numeric(UPSpars[5])
    scedas <- UPSpars[6]
    UPSdf <- get(UPSpars[7])
    hclbin <- "HclusBin"
    hbins <- as.data.frame(cutree(hiclus$upshcl, k = numclust))
    names(hbins) <- hclbin
    dfnnltd <- merge(dframe, hbins, by.x = "row.names", by.y = "row.names", 
        all.x = TRUE)
    dfnnltd[, hclbin] <- as.factor(dfnnltd[, hclbin])
    bins <- length(table(dfnnltd[, hclbin]))
    NNolist <- list(hiclus = UPSpars[1], dframe = UPSpars[2], 
        trtm = trtm, yvar = yvar, numclust = numclust, actclust = bins, 
        scedas = scedas)
    PSmean <- as.matrix(tapply(dfnnltd[, yvar], dfnnltd[, trtm], 
        na.rm = TRUE, mean))
    PStrtm <- paste(trtm, "=", dimnames(PSmean)[[1]])
    PStdif <- paste(PStrtm[2], "minus", PStrtm[1])
    PSvars <- as.matrix(tapply(dfnnltd[, yvar], dfnnltd[, trtm], 
        na.rm = TRUE, var))
    PSfreq <- as.matrix(table(dfnnltd[, trtm]))
    PSvars <- PSvars/PSfreq
    NNolist <- c(NNolist, list(PStdif = PStdif, nnhbindf = hbins, 
        rawmean = PSmean, rawvars = PSvars, rawfreq = PSfreq))
    RATdif <- sum(PSmean[2, ] - PSmean[1, ])
    RATsde <- sum(sqrt(PSvars[2, ] + PSvars[1, ]))
    NNolist <- c(NNolist, list(ratdif = RATdif, ratsde = RATsde))
    PSmean <- tapply(dfnnltd[, yvar], list(dfnnltd[, hclbin], 
        dfnnltd[, trtm]), na.rm = TRUE, mean)
    PSmean <- cbind(matrix(1:bins, bins, 1), PSmean)
    dimnames(PSmean) <- list(1:bins, c("BIN", PStrtm[1], PStrtm[2]))
    PSmean <- as.data.frame(PSmean)
    PSmean[, "BIN"] <- as.factor(PSmean[, "BIN"])
    PSvars <- as.matrix(tapply(dfnnltd[, yvar], list(dfnnltd[, 
        hclbin], dfnnltd[, trtm]), na.rm = TRUE, var))
    PSfreq <- as.matrix(table(dfnnltd[, hclbin], dfnnltd[, trtm]))
    PSvars <- PSvars/PSfreq
    PSvars[PSvars == Inf] <- NA
    PSfreq <- cbind(matrix(1:bins, bins, 1), PSfreq)
    dimnames(PSfreq) <- list(1:bins, c("BIN", PStrtm[1], PStrtm[2]))
    PSfreq <- as.data.frame(PSfreq)
    PSfreq[, "BIN"] <- as.factor(PSfreq[, "BIN"])
    NNolist <- c(NNolist, list(binmean = PSmean, binvars = PSvars, 
        binfreq = PSfreq))
    if (length(table(dfnnltd[, yvar])) > faclev) {
        youtype <- "contin"
        form <- as.formula(paste(yvar, "~", trtm))
        aovdiff <- aov(form, dfnnltd, na.action = na.omit)
        NNolist <- c(NNolist, list(youtype = youtype, faclev = faclev, 
            aovdiff = invisible(summary(aovdiff))))
        sig2 <- 0
        if (bins > 1) {
            form <- as.formula(paste(yvar, "~", hclbin, "+", 
                trtm, "%in%", hclbin))
            aovdiff <- aov(form, dfnnltd, na.action = na.omit)
            raov <- as.matrix(aovdiff$residuals)
            sig2 <- apply(raov, 2, na.rm = TRUE, var) * (length(raov) - 
                1)/aovdiff$df.residual
            NNolist <- c(NNolist, list(form2 = form, bindiff = invisible(summary(aovdiff)), 
                sig2 = sig2))
        }
        pbindif <- PSmean[, 3] - PSmean[, 2]
        nnxmin <- min(pbindif, na.rm = TRUE)
        nnxmax <- max(pbindif, na.rm = TRUE)
        if (scedas == "homo" && sig2 > 0) {
            PSvars <- sig2/as.matrix(table(dfnnltd[, hclbin], 
                dfnnltd[, trtm]))
            PSvars[PSvars == Inf] <- NA
        }
        pbinsde <- sqrt(PSvars[, 2] + PSvars[, 1])
        nnymax <- max(pbinsde, na.rm = TRUE)
        pbinsiz <- sqrt(PSfreq[, 3] + PSfreq[, 2])
        symsiz <- length(dframe[, yvar])
        NNolist <- c(NNolist, list(pbindif = pbindif, pbinsde = pbinsde, 
            pbinsiz = pbinsiz, symsiz = symsiz))
    }
    else {
        youtype <- "factor"
        df3 <- as.data.frame(cbind(dfnnltd[, yvar], dfnnltd[, 
            trtm]))
        df3[, 3] <- dfnnltd[, hclbin]
        df3 <- na.omit(df3)
        names(df3) <- c(yvar, trtm, hclbin)
        df3[, 1] <- as.factor(df3[, 1])
        tab <- table(df3[, 1], df3[, 2])
        NNolist <- c(NNolist, list(youtype = youtype, faclev = faclev, 
            factab = tab))
        tab <- table(df3[, 1], df3[, 2], df3[, 3])
        cumchi <- 0
        cumdf <- 0
        for (i in 1:bins) {
            ht <- chisq.test(tab[, , i])
            if (!is.na(ht$statistic) && ht$statistic < Inf) {
                cumchi <- cumchi + ht$statistic
                cumdf <- cumdf + ht$parameters
            }
        }
        nnymax <- 0
        nnxmin <- 0
        nnxmax <- 0
        NNolist <- c(NNolist, list(cumchi = cumchi, cumdf = cumdf))
    }
    awbdif <- sum(na.omit(((PSmean[, 3] - PSmean[, 2]) * (PSfreq[, 
        3] + PSfreq[, 2]))/sum(PSfreq[, 3] + PSfreq[, 2])))
    awbsde <- sqrt(sum(na.omit((PSvars[, 2] + PSvars[, 1]) * 
        (PSfreq[, 3] + PSfreq[, 2])^2))/(sum(PSfreq[, 3] + PSfreq[, 
        2]))^2)
    wwbdif <- sum(na.omit((PSmean[, 3] - PSmean[, 2])/(PSvars[, 
        2] + PSvars[, 1])/sum(1/na.omit(PSvars[, 2] + PSvars[, 
        1]))))
    wwbsde <- sqrt(1/sum(1/na.omit(PSvars[, 2] + PSvars[, 1])))
    form <- as.formula(paste(yvar, "~", trtm))
    NNolist <- c(NNolist, list(awbdif = awbdif, awbsde = awbsde, 
        wwbdif = wwbdif, wwbsde = wwbsde, form = form))
    accnew <- as.data.frame(cbind("NN", UPSpars[1], UPSpars[2], 
        trtm, yvar, bins, awbdif, awbsde))
    names(accnew) <- c("NNIV", "hicl", "dfrm", "trtm", "yvar", 
        "bins", "tdif", "tdse")
    UPSdf <- as.data.frame(rbind(UPSdf, accnew))
    accnew <- as.data.frame(cbind("NW", UPSpars[1], UPSpars[2], 
        trtm, yvar, bins, wwbdif, wwbsde))
    names(accnew) <- c("NNIV", "hicl", "dfrm", "trtm", "yvar", 
        "bins", "tdif", "tdse")
    UPSdf <- as.data.frame(rbind(UPSdf, accnew))
    assign(UPSpars[7], UPSdf, env = .GlobalEnv)
    if (as.numeric(UPSpars[8]) < nnymax) 
        UPSaccum.pars[8] <<- nnymax
    if (as.numeric(UPSpars[9]) > nnxmin) 
        UPSaccum.pars[9] <<- nnxmin
    if (as.numeric(UPSpars[10]) < nnxmax) 
        UPSaccum.pars[10] <<- nnxmax
    class(NNolist) <- "UPSnnltd"
    NNolist
}
