.supportedModels <- c(
  "BetaGenTEGARCH",
  "BetaTEGARCH",
  "DPMP",
  "DPMP1-I",
  "DPMP1-H",
  "DPMP1-Inv",
  "DPMP2-I",
  "DPMP2-H",
  "DPMP2-Inv",
  "DPMP3-I",
  "DPMP3-H",
  "DPMP3-Inv"
)

.GetModelStr <- function(model){
  if (is.character(model)) {
    modelStr <- model
  }else{
    if (grepl("Rcpp", class(model)[1])) {
      modelStr <- strsplit(class(model)[1], "_")[[1]][2]
    }else{
      modelStr <- class(model)[1]
    }
  }
  return(modelStr)
}

.CheckModelValidity <- function(model){
  modelStr <- .GetModelStr(model)
  return(modelStr %in% .supportedModels)
}

#' @describeIn CreateModel Used exclusively with sring argument.
CreateModelFromStr <- function(model, params=NA){
  if (!is.character(model)) {
    stop("Pass string for model argument.")
  }else{
    if (!any(is.na(params))) {
      model <- new(GASModel, model, params)
    }else{
      model <- new(GASModel, model)
    }
    return(model)
  }
}

#' R wrapper to create GASModel objects.
#'
#' Generates GASModel objects from a string, or wraps \code{GASModel$SetParams}
#' function if the passed object is of class \code{GASModel}.
#'
#' @param model String or GASModel object.
#' @param params Vector or list of initial or new parameters. If NA and the
#' \code{model} argument is a string the model's default parameters are used.
#' If NA and the \code{model} argument is a GASModel object the object's current
#' parameters are maintained.
#'
#' @export
CreateModel <- function(model, params=NA){
  modelValid <- .CheckModelValidity(model)
  if (modelValid) {
    if (is.character(model)) {
      model <- CreateModelFromStr(model, params)
    }else{
      if (!any(is.na(params))) {
        model$SetParams(params)
      }
    }
    return(model)
  }else{
    stop("Specified model is invalid.")
  }
}

#' Computes the mode.
#'
#' Returns the mode of a sample \code{x} using 0.01 width bins.
#'
#' @param x Vector representing a sample.
#'
#' @examples Mode(c(1, 1.8,2, 1.996, 1.991, 2.0005, 2.3, 2.1, 2.5)) # 2
#'
#' @export
Mode <- function(x) {
  ux <- stats::na.omit(unique(x))
  ux <- round(ux, 2)
  ux[which.max(tabulate(match(round(x, 2), ux)))]
}
