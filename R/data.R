#' Prices of the S&P 500 index
#'
#' A dataset containing 10 years worth of daily prices for the S&P 500 index -
#' from 1997-04-24 till 2017-04-21, recorded in reverse chronological order.
#'
#' @docType data
#' @keywords datasets
#' @name spData
#' @usage data(SP500)
#' @format A data frame called \code{spData} with 5032 rows and 7 columns:
#'  \itemize{
#'    \item Date. trading date
#'    \item Open. opening price for the day
#'    \item High. highest price recorded for the day
#'    \item Low. lowest price recorded for the day
#'    \item Close. closing price for the day adjusted for splits
#'    \item Volume. trading volume for the day
#'    \item Adj.Close. closing price for the day adjusted for both dividends and
#'     splits
#'  }
#' @source \url{https://finance.yahoo.com/quote/\%5EGSPC/history?p=\%5EGSPC}
NULL
