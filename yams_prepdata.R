#### Prepping hald pike data for YAMS ####
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


# bring in data
detections <- readRDS('data/dets.RDS')
# ts is the time step
# tag is the tag id
# epo is the UNIX epoch time
# number of seconds since jan 1, 1970 
# frac is fraction of seconds
# serial is the receiver
# sensor val and sensor unit are paired
# sometimes it is depth, sometimes it is ping_num, which 
# is associated only with sync tags - ping_num is the known ping idx of the 
# sync tag, but it can be detected multiple times if it bleeds into other frequencies


hydros <- readRDS('data/hydros.RDS')
# serial is the id of the receiver
# x is the easting
# y is the northing
# z is depth
# sync_tag is the id of the sync tag if present or NA if not present 
# sync tag deployed at that receiver
# idx is just an index used instead of the serial number to designate each receiver

# sync model to estimate how the receivers have drifted (clocks in time, 
# and possibly locations in space)
sync <- readRDS('data/sync.RDS')

# applySync completes the "synchronization"
# you will see that it's only a couple of seconds in difference, 
# but that's important because SoS is about 1500 m/s, which means that a one 
# second difference means we can be off by more than a km!
synceddets <- applySync(toa=detections, hydros=hydros, sync_model=sync)
# adds eposync which is what goes into yaps (for the toa)

hydros <- data.table::data.table(sync$pl$TRUE_H)
# two of the hydros had uncertain positions 
# so we need to use the predicted positions from the sync model
colnames(hydros) <- c('hx','hy','hz')


# these tags were programmed with a random burst interval (rbi) that was given to us
# the sequence of random numbers determining these rbis are included
# in the hald data set: 
seq <- readRDS('data/seq.RDS')
# set the ping type, can either be random known, or random burst, or fixed
ping_type <- 'pbi' # pbi denotes random but known burst interval
# sets the min and max bounds of the burst interval if that's the case
# do this in seconds
rbi_min <- 10 
rbi_max <- 30


aligneddat <- alignBurstSeq(synced_dat=synceddets, burst_seq=seq, 
                             seq_lng_min=25, 
                             rbi_min=rbi_min, rbi_max=rbi_max, 
                             plot_diag=TRUE) %>%
  filter(eposync<1562439511)
saveRDS(aligneddat, 'output/aligneddat.RDS')
# filter out the first 25000 pings
# figures out where the data lie on the known random burst sequence
# this function is only to be used when you have a random known burst interval
# this is how we determine which observations (times of arrival) came from 
# which unobserved random effects (times of pings)
# if you don't have this, then you have to figure out which detections 
# correspond to the same ping based on the temporal difference between detections
# adds a few things
# seq_epo is the expected time of ping 
# generates seq_ping_idx, which is the index of the ping sequence (i.e., ping ID)


# if you have a random known burst interval, you use alignBurstSeq and then 
# buildToaKnownSeq
toalist <- buildToaKnownSeq(seq=seq, aligned_dat=aligneddat, hydros=hydros)
str(toalist)
# toa is the time of arrival matrix
# hydrophones are the columns
# time is the rows

saveRDS(toalist, 'output/toalist.RDS')





