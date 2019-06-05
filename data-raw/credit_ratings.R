rawData <- read.csv("data-raw/CreditRatings.csv");

# Set the rating classes for each observation
igGrades <- c("BBB-", "BBB", "BBB+", "A-", "A", "A+", "AA-", "AA", "AA+",
               "AAA")
sigGrades <- c("BB+", "BB", "BB-", "B-", "B", "B+", "CCC-", "CCC", "CCC+",
                "CC", "C")
defGrades <- c("D", "SD")
igIdx <- which(rawData$splticrm %in% igGrades)
sigIdx <- which(rawData$splticrm %in% sigGrades)
defIdx <- which(rawData$splticrm %in% defGrades)
ratingClass <- vector(mode = "character", length = dim(rawData)[1])
ratingClass[igIdx] <- "IG"
ratingClass[sigIdx] <- "SIG"
ratingClass[defIdx] <- "D"
rawData$ratingclass <- as.factor(ratingClass);

# Set the transition type if any, per observation (takes quite some time)
rawData$transition <- vector(mode = "character", length = dim(rawData)[1])
companyTics <- unique(rawData[ratingClass != 0, "tic"])
pb <- txtProgressBar(min = 0, max = length(companyTics), style = 3)
i <- 0
for (tic in companyTics) {
  compIdx <- which(rawData$tic == tic)
  if (length(unique(rawData[compIdx, "ratingclass"])) > 1) {
    currentClass <- rawData[compIdx, "ratingclass"][1]
    for (j in 2:length(compIdx)) {
      class <- rawData[compIdx, "ratingclass"][j]
      if (class != currentClass) {
        if (currentClass == "IG" && class == "SIG") {
          rawData[compIdx, "transition"][j] <- "IGtoSIG"
        } else if (currentClass == "IG" && class == "D") {
          rawData[compIdx, "transition"][j] <- "IGtoD"
        } else if (currentClass == "SIG" && class == "IG") {
          rawData[compIdx, "transition"][j] <- "SIGtoIG"
        } else if (currentClass == "SIG" && class == "D") {
          rawData[compIdx, "transition"][j] <- "SIGtoD"
        }
      }
      currentClass <- class
    }
  }
  i <- i + 1
  setTxtProgressBar(pb, i)
}
rawData$transition <- as.factor(rawData$transition)
save(rawData, file = "data-raw/tmp1CreditRatings.RData")
#load("data-raw/tmp1CreditRatings.RData")

library("lubridate")
rawData$datadate <- as.Date(rawData$datadate,  "%d/%m/%Y")
# Order the dataframe chronlogically
rawData <- rawData[order(rawData$datadate),]
# Assign a numerical value to each month in the dataset
dates <- floor_date(rawData$datadate, unit = "months")
rawData$nummonth <- cumsum(!duplicated(dates))
# Select the rows at which rating class transtions occur
transitionData <- rawData[rawData[, "transition"] != "",]
# Check if any company transitions more than once in a given month. This
# validates the assumption of an orderly counting process at the monthly freq.
for (i in unique(transitionData$nummonth)) {
  if (any(duplicated(transitionData[transitionData$nummonth == i, "tic"]))) {
    print(i, transitionData[transitionData$nummonth == i, "tic"])
  }
}

# Create dummy variables for the transition types
transitionData$IGtoSIG <- numeric(dim(transitionData)[1])
transitionData$IGtoD <- numeric(dim(transitionData)[1])
transitionData$SIGtoIG <- numeric(dim(transitionData)[1])
transitionData$SIGtoD <- numeric(dim(transitionData)[1])
transitionData$IGtoSIG[transitionData$transition == "IGtoSIG"] <- 1
transitionData$IGtoD[transitionData$transition == "IGtoD"] <- 1
transitionData$SIGtoIG[transitionData$transition == "SIGtoIG"] <- 1
transitionData$SIGtoD[transitionData$transition == "SIGtoD"] <- 1

# Reset numerical month values to where the first transition event occured
transitionData$nummonth <- transitionData$nummonth - 86 # start 1986/01/31

# Determine when transition events occured
eventBoolIdx = !duplicated(transitionData$nummonth)
events <- cumsum(eventBoolIdx)
numEvents <- max(events)
# Directly compute the time difference in months (\tau_{t} - \tau_{t-1}).
# We'll never need the level values of \tau for computations.
diffTau <- c(1, diff(transitionData$nummonth[eventBoolIdx]))

# Determine number of companies that transiton and how many companies could have
# potentially transitioned for each possible rating transition at each event.
companies <- unique(transitionData$tic)
companyID <- vector(mode = "integer", length = dim(transitionData)[1])
numCompanies <- length(companies)
# Vector to keep track of the state of all companies. Value is 1 if the company
# is in the IG class, 2 if in the SIG class and 3 implies default.
companyStates <- numeric(numCompanies)
# Set company IDs and initialize the companyStates vector (derriving state
# backwards from the first observed transitions).
for (i in 1:numCompanies) {
  ix <- which(transitionData$tic == companies[i])
  companyID[ix] <- i
  if (transitionData[as.vector(ix)[1], "IGtoSIG"] == 1) {
    companyStates[i] <- 1
  } else if ((transitionData[as.vector(ix)[1], "IGtoD"] == 1)) {
    companyStates[i] <- 1
  } else if ((transitionData[as.vector(ix)[1], "SIGtoIG"] == 1)) {
    companyStates[i] <- 2
  } else if ((transitionData[as.vector(ix)[1], "SIGtoD"] == 1)) {
    companyStates[i] <- 2
  }
}
transitionData$companyID <- companyID

# Iterate over all events setting transitions and potential transitions
transitions <- matrix(0, numEvents, 4)
possibleTranstions <- matrix(0, numEvents, 4)
for (i in 1:numEvents) {
  possibleTranstions[i, 1] <- sum(companyStates == 1)
  possibleTranstions[i, 2] <- sum(companyStates == 1)
  possibleTranstions[i, 3] <- sum(companyStates == 2)
  possibleTranstions[i, 4] <- sum(companyStates == 2)
  ix <- as.vector(events == i)
  companyStates[transitionData[ix, "companyID"]] <- (
    transitionData[ix, "IGtoSIG"] * 2 +
      transitionData[ix, "IGtoD"] * 3 +
      transitionData[ix, "SIGtoIG"] +
      transitionData[ix, "SIGtoD"] * 3
  )
  transitions[i, ] <- colSums(
    transitionData[ix, c("IGtoSIG", "IGtoD", "SIGtoIG", "SIGtoD")])
}

# Save data
save(transitions, possibleTranstions, transitionData,
     file = "data/CreditRatings.RData")
