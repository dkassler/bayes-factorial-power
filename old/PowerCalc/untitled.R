# Created on 4-4-18 to replace the old power framework with something properly
# reproducible. These functions are meant to be stable between versions,
# allowing us to vary only a set of (easily documented) arguments to produce all
# our charts and graphs for the paper.

makeExperimentData <- function(
  N,                 # number of participants
  config,            # integer(4) giving #levels for arms
  cauchyScale = 5,
  diff = .25         # max - min treatment
) {
  #Assign level sizes
  stopifnot(length(config) == 4)
  nA <- config[1]
  nB <- config[2]
  nC <- config[3]
  nD <- config[4]
}