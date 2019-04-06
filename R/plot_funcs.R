#' Trace plots of MCMC samples
#'
#' Creates plots of iterations vs. sampled values for each chain. Plots are laid
#' out in a grid with columns containing plots grouped per MCMC method and rows
#' containing the plots grouped by parameter.
#'
#' @param samples Matrix containing MCMC samples. Each column contains a
#'  seperate chain of a unique combination of a MCMC method specified in
#'  \code{mcmcNames} and parameter specified in \code{parNames}. Columns
#'  must be named according to the convention "\code{mcmcName} \code{parName}"
#'  (e.g. "RWMH A").
#' @param mcmcNames Vector or list of strings or epxressions of MCMC method
#'  names. Must be list if expressions are used.
#' @param parNames Vector or list of strings or epxressions of parameter names.
#'  Must be list if expressions are used.
#' @param iters Indexing vector used to select the iterations to plot.
#' @param cols Vector of elements that can be passed to \code{\link{plot}} as
#'  \code{col} (color) values applied per parameter. Elements are repeated to
#'  match length of parNames.
#'
#' @export
PlotTraces <- function(samples, mcmcNames, parNames, iters=-(1:1000),
                       cols = (1:8)){
  numMcmc <- length(mcmcNames)
  numPars <- length(parNames)
  par(mfcol = c(numPars, numMcmc), mar = c(4, 4, 2.5, 2))
  cols <- rep(cols, ceiling(length(parNames) / length(cols)))
  for (mcmcName in mcmcNames) {
    c <- 0
    for (parName in parNames) {
      c <- c + 1
      colName <- paste(mcmcName, parName, sep = " ")
      if (colName %in% colnames(samples)) {
        stats::plot.ts(
          samples[iters, colName],
          col = cols[c],
          axes = TRUE,
          ann = FALSE
        )
        if (paste(parName) == paste(parNames[1][[1]])) {
          title(main = mcmcName, font.main = 1, cex.main = 1.2, line = 1.2)
        }
        title(xlab = "Draws", ylab = parName, line = 2.5)
      }else{
        warning(sprintf("%s is not included in the sample\n", colName))
      }
    }
  }
}

#' Histogram plots of MCMC samples
#'
#' Creates histogram per MCMC chain. Plots are laid out in a grid with columns
#' containing plots grouped per MCMC method and rows containing the plots
#' grouped by parameter.
#'
#' @inheritParams PlotTraces
#' @param burn Integer specifying the num of iterations to discard of each
#'  chain.
#' @param bins Integer specifying the num of bins to create for the histogram.
#'
#' @export
PlotHists <- function(samples, mcmcNames, parNames, burn=1000, bins=50){
  numMcmc <- length(mcmcNames)
  numPars <- length(parNames)
  par(mfcol = c(numPars, numMcmc), mar = c(4, 4, 2.5, 2))
  for (mcmcName in mcmcNames) {
    for (parName in parNames) {
      colName <- paste(mcmcName, parName, sep = " ")
      if (colName %in% colnames(samples)) {
        samples_ <- samples[-(1:burn), colName]
        hist(
          samples[-(1:burn), colName],
          axes = FALSE,
          col = "slategray1",
          breaks = seq(min(samples_), max(samples_), length.out = bins + 1),
          xlab = "",
          ylab = "",
          main = ""
        )
        title(xlab = parName, main = "", line = 2.2)
        title(ylab = "Frequency", main = "", line = 2)
        if ((paste(parName) == paste(parNames[1][[1]]))) {
          title(main = mcmcName, cex.main = 1.2, font.main = 1, line = 1.2)
        }
        axis(side = 1, line = -0.0, tck = -0.03, mgp = c(3, .5, 0))
        axis(side = 2, line = -0.35, tck = -0.03, mgp = c(3, .5, 0))
      }else{
        warning(sprintf("%s is not included the samples", colName))
      }
    }
  }
}

#' Scatter plots of joint posterior distributions
#'
#' Creates scatter plots of sample values of one parameter vs another.
#'
#' @param samples Matrix containing MCMC samples. Each column contains a
#'  seperate chain of a parameter specified in \code{paramPairs} and is named
#'  accordingly (e.g. given a \code{paramPair} c("A", "B") \code{samples} is
#'  expected to contain a column named "A" and a column named "B" ).
#' @param paramPairs A list of \code{parNames} (see e.g.
#'  \code{\link{PlotTraces}} for a description of a parNames object).
#' @param lbParams Named numeric vector representing the lower plot limits of
#'  the parameters. Names correspond to the \code{parName}s specified in
#'  \code{paramPairs}.
#' @param ubParams Named numeric vector representing the upper plot limits of
#'  the parameters. Names correspond to the \code{parName}s specified in
#'  \code{paramPairs}.
#' @inheritParams PlotHists
#'
#' @export
PlotJointDists <- function(samples, paramPairs, lbParams, ubParams, burn=1000){
  numPairs <- length(paramPairs)
  par(mfrow = c(ceiling(numPairs / 2), 2), mar = c(4, 4, 2, 2))
  for (pair in paramPairs) {
    xParam <- paste(pair[1])
    yParam <- paste(pair[2])
    if ((xParam %in% colnames(samples)) & (yParam %in% colnames(samples))) {
      plot(
        samples[-(1:burn), xParam],
        samples[-(1:burn), yParam],
        cex = c(0.1, 0.1),
        ylab = "",
        xlab = "",
        xlim = c(lbParams[xParam], ubParams[xParam]),
        ylim = c(lbParams[yParam], ubParams[yParam])
      )
      title(xlab = pair[1], ylab = pair[2], cex.lab = 1.2, main = "",
            line = 2.5)
    }else{
      warning(sprintf("%s is not included in the sample\n", paramPairs))
    }
  }
}

#' Autocorrelation function plots of MCMC samples
#'
#' Plots the autocorrelation function (ACF) (see \code{\link[stats]{acf}}) per
#' MCMC chain. Plots are laid out in a grid with columns containing plots
#' grouped per MCMC method and rows containing the plots grouped by parameter.
#'
#' @inheritParams PlotHists
#' @param lags Integer specifying the num of lags upto which to compute the ACF.
#'
#' @export
PlotACFs <- function(samples, mcmcNames, parNames, burn=1000, lags=50){
  numMcmc <- length(mcmcNames)
  numPars <- length(parNames)
  par(mfcol = c(numPars, numMcmc), mar = c(4, 4, 2.5, 2))
  for (mcmcName in mcmcNames) {
    for (parName in parNames) {
      colName <- paste(mcmcName, parName, sep = " ")
      if (colName %in% colnames(samples)) {
        unused <- stats::acf(
          samples[-(1:burn), colName],
          lwd = 2.5,
          ylab = "",
          axes = FALSE,
          ylim = c(-0.05, 1),
          main = "",
          lag.max = lags,
          ci = 0,
          col = "deepskyblue2",
          xlab = ""
        )
        title(main = mcmcName, font.main = 1, cex.main = 1.2, line = 1.2)
        title(xlab = "Lags", line = 2)
        title(ylab = bquote(rho(.(parName[[1]]))), cex.lab = 1.2, line = 2)
        axis(side = 1, line = -.25, tck = -.03, mgp = c(3, .5, 0))
        axis(side = 2, line = -.1, tck = -.03, mgp = c(3, .5, 0))
      }else{
        warning(sprintf("%s is not included the samples", colName))
      }

    }
  }
}

#' Plot highest posterior density region over time
#'
#' Plots the (100 * \code{alpha})% highest posterior density (HPD) boundaries of
#' a statstic for which posterior samples are provided over a certain time
#' frame.
#'
#' @param statDraws Matrix of posterior samples of shape \code{numDraws} by
#'  \code{numObs}, so each row represents one posterior sample of a time series
#'  of the statistic.
#' @param observations Numeric vector of length \code{numObs} to plot as
#'  observational proxy for the statistic.
#' @param dates Vector of strings of length \code{numObs} representing the
#' dates at which the statsitic is observed.
#' @param startDate String representing the starting date for the plot. In
#'  combination with \code{endDate} it allows the user to specify a specific
#'  date range contained in \code{dates} for display purposes. If NULL the
#'  first date in dates is used. If not NULL it must be an element of dates.
#' @param endDate Same as \code{startDate}, but for the end of the date range.
#' @param alpha Numeric specifying the posterior probability mass contained in
#'  the HPD region.
#' @param ylab Y axis label passed to \code{\link{plot}}.
#' @param statStr String to describe the statistic.
#' @param obsStr String to describe the observation series.
#' @param statsML Numeric vector of length \code{numObs} (Optional). Time series
#'  of the statsitic's value under maximum likelihood (ML) parameters.
#' @param zoomDate String specifying a date for which to display a small
#'  histogram of the posterior sample of the statistic in the top right corner
#'  of the plot (Optional). If NULL, no such histogram is displayed.
#' @param xy1 Location specifiers for drawing the top line from the
#'  \code{zoomDate} inside the main plot to the histogram (Optional). If NULL,
#'  no such line is displayed.
#' @param xy2 Same as \code{xy1}, but for the bottom line.
#' @param modeCol Color for the mode. Passed to \code{col} argument for
#'  \code{\link{plot}}.
#' @param MLCol Color for the \code{statsML}. Passed to \code{col} argument
#'  for \code{\link{plot}}.
#' @param fillCol Color to fill the HPD region. Passed to \code{col} argument
#'  for \code{\link{polygon}}.
#' @param borderCol Color for the borders of the HPD region. Passed to
#'  \code{col} argument \code{\link{plot}}.
#'
#' @export
PlotHPDOverTime <- function(statDraws, observations, dates, startDate=NULL,
                            endDate=NULL,
                            alpha=0.99, ylab=NULL, statStr="", obsStr="Obs",
                            statsML=NULL, zoomDate=NULL, xy1=NULL, xy2=NULL,
                            modeCol="red", MLCol="blue", fillCol="pink1",
                            borderCol="red2"){
  mcmcStats <- coda::mcmc(statDraws)
  hpd <- coda::HPDinterval(mcmcStats, prob = alpha)
  statsUb <- hpd[, 2]
  statsLb <- hpd[, 1]
  numObs <- dim(statDraws)[2]
  statsMode <- numeric(numObs)
  for (i in 1:numObs) {
    statsMode[i] <- Mode(mcmcStats[, i])
  }

  startIdx <- which(dates == startDate)
  endIdx <- which(dates == endDate)
  idxs <- (startIdx:endIdx)

  par(mfrow = c(1, 1), mar = c(5, 5, 4, 4))
  plot(
    x = dates[idxs],
    observations[idxs],
    xaxt = "n",
    type = "n",
    xlab = "Years",
    ylab = ylab,
    lab = c(10, 5, 7),
    ylim = c(0, 100)
  )
  axis.Date(1, at = seq(min(dates[idxs]), max(dates[idxs]), by = "2 mon"),
            format = "%Y-%m")
  polygon(
    c(dates[idxs], rev(dates[idxs])),
    c(statsUb[idxs], rev(statsLb[idxs])),
    col = fillCol,
    border = NA
  )
  lines(dates[idxs], abs(observations[idxs]))
  lines(dates[idxs], statsMode[idxs], lwd = 1.5, col = modeCol)
  lines(dates[idxs], statsUb[idxs], col = borderCol, lty = 2)
  lines(dates[idxs], statsLb[idxs], col = borderCol, lty = 2)
  if (!is.null(statsML)) {
    lines(dates[idxs], statsML[idxs], lty = 3, col = MLCol, lwd = 1.5)
  }

  leg <- c(obsStr, paste("Mode", statStr))
  if (!is.null(statsML)) leg <- c(leg, paste("ML estimate", statStr))
  leg <- c(leg, paste(sprintf("%i", alpha * 100), "% HPD region", sep = ""))
  if (!is.null(statsML)) {
    legend(
      "topleft",
      legend = leg,
      bty = "n",
      col = c("black", modeCol, MLCol, NA),
      lty = c(1, 1, 2, 3),
      lwd = c(1, 1.5, 1.5, 1),
      pch = c(NA, NA, NA, 22),
      pt.bg = c(NA, NA, NA, NA),
      fill = c(0, 0, 0, fillCol),
      border = c(NA, NA, NA, borderCol)
    )
  }else{
    legend(
      "topleft",
      legend = leg,
      bty = "n",
      col = c("black", modeCol, NA),
      lty = c(1, 1, 3),
      lwd = c(1, 1.5, 1),
      pch = c(NA, NA, 22),
      pt.bg = c(NA, NA, NA),
      fill = c(0, 0, fillCol),
      border = c(NA, NA, borderCol)
    )
  }

  if (!is.null(zoomDate)) {
    if (!is.null(xy1)) {
      lines(xy1, lty = 2)
    }
    if (!is.null(xy2)) {
      lines(xy2, lty = 2)
    }

    u <- par("usr")
    v <- c(grconvertX(u[1:2], "user", "ndc"),
           grconvertY(u[3:4], "user", "ndc"))
    v <- c((v[1] + v[2]) / 1.5, v[2], (v[3] + v[4]) / 1.8, v[4])
    par(fig = v, new = TRUE, mar = c(0, 0, 0, 0))
    zoomIdx <- which(dates[idxs] == as.Date(zoomDate))
    x <- statDraws[, idxs][, zoomIdx]
    hist(x, xlab = "", ylab = "", axes = FALSE, prob = TRUE,
         mar = c(1.6, 1.6, .1, .1), breaks = 50, main = "")
    title_ = bquote(.(as.character(dates[idxs][zoomIdx])))
    title(main = title_, font.main = 1, line = -1)
    mtext(text = statStr, side = 1, line = 1.5)
    mtext(text = "Probability",  side = 2, line = 1.5)
    axis(side = 1, line = -0.1, tck = -0.03, mgp = c(3, .5, 0))
    axis(side = 2, line = -0.1, tck = -0.03, mgp = c(3, .5, 0))

    Density <- stats::density(x, adjust = 2)
    lines(stats::density(x), col = "green3")
    leg <- c(
      "Fitted density",
      paste(sprintf("%i", alpha * 100), "% HPD region", sep = "")
    )
    legend("right", legend = leg, bty = "n", col = c("green3", "red2"),
           lty = c(1, 2), lwd = c(1, 2), cex = 1)
    lines(c(statsLb[idxs][zoomIdx], statsLb[idxs][zoomIdx]),
          c(0, 100 * Density$y[15]), col = "red2", lwd = 2)
    lines(c(statsUb[idxs][zoomIdx], statsUb[idxs][zoomIdx]),
          c(100 * Density$y[15], 0), col = "red2", lwd = 2)
  }
}
