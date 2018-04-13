setwd(here::here('PowerCalc'))

library(tidyverse)
library(rstan)
library(doParallel)

#' Load functions in `lib` folder.

for (file in dir('lib')) {
  message(sprintf('Sourcing %s...', file))
  source(file.path('lib', file))
}
