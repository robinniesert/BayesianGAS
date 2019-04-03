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
        plot.ts(
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

PlotHistograms <- function(samples, mcmcNames, parNames, burn=1000, bins=50){
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
        if (paste(parName) == paste(expression(upsilon))) {
          title(xlab = parName, cex.lab = 1.2, main = "", line = 2.2)
        }else{
          title(xlab = parName, main = "", line = 2.2)
        }
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

PlotACFs <- function(samples, mcmcNames, parNames, burn=1000, lags=50){
  numMcmc <- length(mcmcNames)
  numPars <- length(parNames)
  par(mfcol = c(numPars, numMcmc), mar = c(4, 4, 2.5, 2))
  for (mcmcName in mcmcNames) {
    for (parName in parNames) {
      colName <- paste(mcmcName, parName, sep = " ")
      if (colName %in% colnames(samples)) {
        acf <- acf(
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

PlotHPDOverTime <- function(statDraws, observations, times, startDate, endDate,
                            alpha=0.99, ylab=NULL, statStr="", obsStr="Obs",
                            statsML=NULL, zoomDate=NULL, xy1=NULL, xy2=NULL,
                            modeCol="red", MLCol="blue", fillCol="pink1",
                            borderCol="red2"){
  mcmcStats <- coda::mcmc(statDraws)
  hpd <- coda::HPDinterval(mcmcStats, prob = alpha)
  statsUb <- hpd[, 2]
  statsLb <- hpd[, 1]
  statsMode <- numeric(numObs)
  for (i in 1:numObs) {
    statsMode[i] <- Mode(mcmcStats[, i])
  }

  startIdx <- which(times == startDate)
  endIdx <- which(times == endDate)
  idxs <- (startIdx:endIdx)

  par(mfrow = c(1, 1), mar = c(5, 5, 4, 4))
  plot(
    x = times[idxs],
    observations[idxs],
    xaxt = "n",
    type = "n",
    xlab = "Years",
    ylab = ylab,
    lab = c(10, 5, 7),
    ylim = c(0, 100)
  )
  axis.Date(1, at = seq(min(times[idxs]), max(times[idxs]), by = "2 mon"),
            format = "%Y-%m")
  polygon(
    c(times[idxs], rev(times[idxs])),
    c(statsUb[idxs], rev(statsLb[idxs])),
    col = fillCol,
    border = NA
  )
  lines(times[idxs], abs(observations[idxs]))
  lines(times[idxs], statsMode[idxs], lwd = 1.5, col = modeCol)
  lines(times[idxs], statsUb[idxs], col = borderCol, lty = 2)
  lines(times[idxs], statsLb[idxs], col = borderCol, lty = 2)
  if (!is.null(statsML)) {
    lines(times[idxs], statsML[idxs], lty = 3, col = MLCol, lwd = 1.5)
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
    if (!is.null(xy1) & !is.null(xy2)) {
      lines(xy1, lty = 2)
      lines(xy2, lty = 2)
    }

    u <- par("usr")
    v <- c(grconvertX(u[1:2], "user", "ndc"),
           grconvertY(u[3:4], "user", "ndc"))
    v <- c((v[1] + v[2]) / 1.5, v[2], (v[3] + v[4]) / 1.8, v[4])
    par(fig = v, new = TRUE, mar = c(0, 0, 0, 0))
    zoomIdx <- which(times[idxs] == as.Date(zoomDate))
    x <- statDraws[, idxs][, zoomIdx]
    hist(x, xlab = "", ylab = "", axes = FALSE, prob = TRUE,
         mar = c(1.6, 1.6, .1, .1), breaks = 50, main = "")
    title_ = bquote(.(as.character(times[idxs][zoomIdx])))
    title(main = title_, font.main = 1, line = -1)
    mtext(text = statStr, side = 1, line = 1.5)
    mtext(text = "Probability",  side = 2, line = 1.5)
    axis(side = 1, line = -0.1, tck = -0.03, mgp = c(3, .5, 0))
    axis(side = 2, line = -0.1, tck = -0.03, mgp = c(3, .5, 0))

    Density <- density(x, adjust = 2)
    lines(density(x), col = "green3")
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
