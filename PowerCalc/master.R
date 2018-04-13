# Master script for power calculations. Running this script should reproduce all
# figures and data used in the paper.

setwd(here::here('PowerCalc'))

library(tidyverse)
library(rstan)
library(doParallel)

#' Load functions in `lib` folder.

for (file in dir('lib')) {
  message(sprintf('Sourcing %s...', file))
  source(file.path('lib', file))
}

#' Run functions in `src` folder.

for (file in dir('src')) {
  message(sprintf('Sourcing %s...', file))
  source(file.path('src', file))
}
