library(markovchain)

N <- 7; p <- 1/2; q <- 1 - p;
stateNames <- as.character( 0:N )
## Be careful here, because states numbered from 0,
## but R indexes from 1
transMatrix <- matrix(0, N+1,N+1)
transMatrix[1,1] <- q
transMatrix[1,2] <- p
transMatrix[N+1, N  ] <- q
transMatrix[N+1, N+1] <- p
for (row in 2:N) {
    transMatrix[row, row-1] <- ((row-1)/N)*q
    transMatrix[row,row] <- ((N-(row-1))/N)*q + ((row-1)/N)*p
    transMatrix[row,row+1] <- ((N-(row-1))/N)*p
}

startState <- "0"

altEhernfest <- new("markovchain", transitionMatrix=transMatrix,
                 states=stateNames, name="AltEhernfest")


pathLength <- c(199, 399, 599, 999)
nTrials <- c(200, 400, 800, 1600)

eps <- 0.01

results <- matrix(0, length(pathLength), length(nTrials))

for (pL in 1:length(pathLength)) {
for  (nT in 1:length(nTrials)) {

outlierThreshold <- floor( eps * (pathLength[pL] + 1) ) + 1

nEpsOutlier<- 0

for (i in 1:nTrials[nT]) {
            
    path <- rmarkovchain(n=pathLength[pL],
                         object = altEhernfest,
                         t0 = startState)
    fullPath <- c(startState, path)
    labelFullPath <- as.numeric(fullPath)
    labelStartState <- as.numeric(startState)
    if (!identical(sort(labelFullPath)[outlierThreshold], labelStartState)) {
                nEpsOutlier <- nEpsOutlier + 1
    }
}

probEpsOutlier <- nEpsOutlier/nTrials[nT]

results[pL, nT] <- probEpsOutlier
}
}

results

## NAME: altEhernfestOutliers.R
## USAGE: within R, at interactive prompt
##        source("altEhernfestOutliers.R")
## REQUIRED ARGUMENTS: none
## OPTIONS: none
## DESCRIPTION: An R script using the markovchain library to 
##              simulate an alternative Ehrenfest urn model.            
## DIAGNOSTICS: some diagnostics implicit from markovchain
## CONFIGURATION AND ENVIRONMENT:  library(markovchain)
## DEPENDENCIES:  library(markovchain)
## INCOMPATIBILITIES: none known
## PROVENANCE: created Steve Dunbar, based on documentation in
## "The markovchain Package: A Package for Easily Handling Discrete
##  Markov Chains in R" by G. Spedicato, et al.
## BUGS AND LIMITATIONS: none known
## FEATURES AND POTENTIAL IMPROVEMENTS:
## 1.  Make a random start state.
## AUTHOR:  Steve Dunbar
## VERSION: Version 1.0 as of Tue Jan 22 07:43:53 CST 2019
## KEYWORDS: Markov chain, stable distribution

