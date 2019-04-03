supportedModels <- c("BetaGenTEGARCH", "BetaTEGRACH")

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
  return(modelStr %in% supportedModels)
}

CreateModelFromStr <- function(modelStr, params=NA){
  if (!is.character(modelStr)) {
    stop("Pass string for model argument.")
  }else{
    if (!any(is.na(params))) {
      model <- new(GASModel, modelStr, params)
    }else{
      model <- new(GASModel, modelStr)
    }
    return(model)
  }
}

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

Mode <- function(x) {
  ux <- na.omit(unique(x))
  ux <- round(ux, 2)
  ux[which.max(tabulate(match(round(x, 2), ux)))]
}
