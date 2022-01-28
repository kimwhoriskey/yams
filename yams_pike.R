#### Fitting YAMS to hald pike data ####
# author: Kim Whoriskey and Henrik Baktoft


setwd("~/Desktop/yams")
rm(list=ls())	
library(data.table)
library(ggplot2)
library(yaps)
library(TMB)
library(dplyr) 
library(gridExtra)
library(mgcv)
library(sp)
source('yapsfunctions.R')
source('issmfunctions.R')
source("simtrack.R")
fmf <- wesanderson::wes_palette("FantasticFox1", type = "discrete")


# TMB likelihoods
compile("yaps_ssm.cpp")
dyn.load(dynlib("yaps_ssm"))
compile("yaps_hmm.cpp")
dyn.load(dynlib("yaps_hmm"))


# bring in sync data to extract hydro positions
sync <- readRDS('data/sync.RDS')
hydros <- data.table::data.table(sync$pl$TRUE_H)
colnames(hydros) <- c('hx','hy','hz')

# bring in time of arrival data for model fitting (see yams_prepdata.R)
toalist <- readRDS('output/toalist.RDS')

# extract toa and burst interval vector
toa <- toalist$toa
bivec <- toalist$seq
# 20,000 pings, so we have to split the data up

# manually input rbi_min and _max, which are 10 and 30 for this dataset
# burst interval mins and maxes
rbi_min <- 10
rbi_max <- 30









########################################################
### run all of the mods for the pike 



# first split the data into five groups
modidx <- rep(1:5, each=5000) # 32539
toamats <- split(data.table(toa), factor(modidx)) %>% lapply(as.matrix)
bivecs <- split(bivec, factor(modidx))



# going to run models four times
# have starting values for hmm at log(1), thats svhmm1
# have starting values of hmm to match ssm, ie at log(60), svhmm60
# each of these for two states and three states, m2 vs m3
# all done with a timescale of 60


# 2 states, hmm at log(1)

mods <- list()
for(i in 1:length(toamats)){
  mods[[i]] <- tryCatch(fit_issm(inp_init = getmodinp(hydros,
                                                      toamats[[i]], sdInits=1,
                                                      rbi_min=rbi_min, rbi_max=rbi_max, biTable=bivecs[[i]],
                                                      m = 1, b = NULL, timescale=60, dstart=60),
                                 inp_ssm = getmodinp(hydros,
                                                     toamats[[i]], sdInits=1,
                                                     rbi_min=rbi_min, rbi_max=rbi_max, biTable=bivecs[[i]],
                                                     m=2, b = NULL, timescale=60, dstart=60),
                                 maxsteps=10,
                                 m=2,
                                 logdstart = rep(log(1), 2), # this sets the starting values of the hmm
                                 fixsteps=TRUE,
                                 allowfc=TRUE,
                                 jiggle_fc=0.01,
                                 initssmargs = list(inner_control = list(maxit=1500),
                                                    mapping = list(working_A = factor(matrix(NA)))),
                                 ssmargs = list(inner_control = list(maxit=1500),
                                                mapping = list(#logD_xy = factor(c(NA, NA)),
                                   working_A = factor(matrix(NA, nrow=2, ncol=2))),
                                   optimizer='nlminb'),
                                 timescalehmm=60,
                                 setinitstatdist=1), error=function(e)e)
  paste('finished mod', i)
}
saveRDS(mods, file="output/pikemodsm2ts60svhmm1.RDS")


# 2 states, hmm at log(60)


mods <- list()
for(i in 1:length(toamats)){
  mods[[i]] <- tryCatch(fit_issm(inp_init = getmodinp(hydros,
                                                      toamats[[i]], sdInits=1,
                                                      rbi_min=rbi_min, rbi_max=rbi_max, biTable=bivecs[[i]],
                                                      m = 1, b = NULL, timescale=60, dstart=60),
                                 inp_ssm = getmodinp(hydros,
                                                     toamats[[i]], sdInits=1,
                                                     rbi_min=rbi_min, rbi_max=rbi_max, biTable=bivecs[[i]],
                                                     m=2, b = NULL, timescale=60, dstart=60),
                                 maxsteps=10,
                                 m=2,
                                 logdstart = rep(log(60), 2), # this sets the starting values of the hmm
                                 fixsteps=TRUE,
                                 allowfc=TRUE,
                                 jiggle_fc=0.01,
                                 initssmargs = list(inner_control = list(maxit=1500),
                                                    mapping = list(working_A = factor(matrix(NA)))),
                                 ssmargs = list(inner_control = list(maxit=1500),
                                                mapping = list(#logD_xy = factor(c(NA, NA)),
                                   working_A = factor(matrix(NA, nrow=2, ncol=2))),
                                   optimizer='nlminb'),
                                 timescalehmm=60,
                                 setinitstatdist=1), error=function(e)e)
  paste('finished mod', i)
}
saveRDS(mods, file="output/pikemodsm2ts60svhmm60.RDS")



# 3 states, hmm at log(1)

mods <- list()
for(i in 1:length(toamats)){
  mods[[i]] <- tryCatch(fit_issm(inp_init = getmodinp(hydros,
                                                      toamats[[i]], sdInits=1,
                                                      rbi_min=rbi_min, rbi_max=rbi_max, biTable=bivecs[[i]],
                                                      m = 1, b = NULL, timescale=60, dstart=60),
                                 inp_ssm = getmodinp(hydros,
                                                     toamats[[i]], sdInits=1,
                                                     rbi_min=rbi_min, rbi_max=rbi_max, biTable=bivecs[[i]],
                                                     m=3, b = NULL, timescale=60, dstart=60),
                                 maxsteps=10,
                                 m=3,
                                 logdstart = rep(log(1), 3), # this sets the starting values of the hmm
                                 fixsteps=TRUE,
                                 allowfc=TRUE,
                                 jiggle_fc=0.01,
                                 initssmargs = list(inner_control = list(maxit=1500),
                                                    mapping = list(working_A = factor(matrix(NA)))), 
                                 ssmargs = list(inner_control = list(maxit=1500),
                                                mapping = list(#logD_xy = factor(c(NA, NA)),
                                   working_A = factor(matrix(NA, nrow=3, ncol=3))),
                                   optimizer='nlminb'),
                                 timescalehmm=60, 
                                 setinitstatdist=1), error=function(e)e)
  paste('finished mod', i)
}
saveRDS(mods, file="output/pikemodsm3ts60svhmm1.RDS")




# 3 states, hmm at log(60)


mods <- list()
for(i in 1:length(toamats)){
  mods[[i]] <- tryCatch(fit_issm(inp_init = getmodinp(hydros,
                                                      toamats[[i]], sdInits=1,
                                                      rbi_min=rbi_min, rbi_max=rbi_max, biTable=bivecs[[i]],
                                                      m = 1, b = NULL, timescale=60, dstart=60),
                                 inp_ssm = getmodinp(hydros,
                                                     toamats[[i]], sdInits=1,
                                                     rbi_min=rbi_min, rbi_max=rbi_max, biTable=bivecs[[i]],
                                                     m=3, b = NULL, timescale=60, dstart=60),
                                 maxsteps=10,
                                 m=3,
                                 logdstart = rep(log(60), 3), # this sets the starting values of the hmm
                                 fixsteps=TRUE,
                                 allowfc=TRUE,
                                 jiggle_fc=0.01,
                                 initssmargs = list(inner_control = list(maxit=1500),
                                                    mapping = list(working_A = factor(matrix(NA)))), 
                                 ssmargs = list(inner_control = list(maxit=1500),
                                                mapping = list(#logD_xy = factor(c(NA, NA)),
                                   working_A = factor(matrix(NA, nrow=3, ncol=3))),
                                   optimizer='nlminb'),
                                 timescalehmm=60, 
                                 setinitstatdist=1), error=function(e)e)
  paste('finished mod', i)
}
saveRDS(mods, file="output/pikemodsm3ts60svhmm60.RDS")










