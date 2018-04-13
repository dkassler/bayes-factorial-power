# This script manages calling the simulations code in a more cleaned up way than
# the previous all in one approach. The sim function is sourced in from its own
# script. This script deals with setting up the designs and calling code

#source("powerSim4varNonAdaptive_simfunction.R")
source("PowerCalc/powerSim4varNonAdaptive_simfunction.R")

library(tidyverse)

# Parameters --------------------------------------------------------------

nSims <- 1000
design_params <- list(
  nArms = list(
    c(2,2,2,2), #16
    #c(2,2,2,3), #24
    c(2,2,3,3), #36
    ##c(2,2,3,4), #48
    #c(2,3,3,3), #54
    ##c(2,2,3,5), #60
    ##c(2,3,3,4), #72
    c(3,3,3,3), #81
    ##c(2,3,3,5), #90
    #c(3,3,3,4), #108
    ##c(3,3,3,5), #135
    ##c(2,4,5,5), #160
    c(3,3,4,4), #144
    ##c(3,3,4,5), #180
    #c(3,4,4,4), #192
    ##c(2,4,5,5), #200
    ##c(3,3,5,5), #225
    c(4,4,4,4), #256
    #c(4,4,4,5), #320
    c(4,4,5,5)#, #400
    #c(4,5,5,5), #500
    #c(5,5,5,5),  #625
  ),
  N = seq(500, 5000, by = 500)
)


# Build # arms ------------------------------------------------------------

fnum <- expand.grid(data_frame(nA = 2:10, nB = nA, nC = nA, nD = nA))

design <- fnum %>% 
  mutate(prod = nA * nB * nC * nD) %>% 
  filter(prod <= 400) %T>% 
  qplot(data = ., x = prod, geom = "histogram") %>%
  mutate(bin = cut(prod, 20)) %>% 
  group_by(bin) %>% 
  sample_n(2, replace = TRUE) %>% 
  ungroup() %>% 
  mutate(N = sample(500:10000, nrow(.), replace = TRUE))

library(doParallel)

cl <- makeCluster(42)
registerDoParallel(cl)
MDEs <- foreach(i = seq(nrow(design))) {
  sim(
    N = design[i, "N"],
    nA = design[i, "nA"],
    nB = design[i, "nB"],
    nC = design[i, "nC"],
    nD = design[i, "nD"]
  )
}

design %>% 
  mutate(MDE = MDEs) %>% 
  saveRDS(sprintf("output_C_%s.Rds", format(Sys.time(), "%Y%m%d%H%M%S")))

