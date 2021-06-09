# post hoc for yams pike analysis


#### set up working directory

setwd("~/Desktop/yams")
rm(list=ls())	
library(data.table)
library(ggplot2)
library(yaps)
library(TMB)
library(tidyverse) # uffff - ugly...
library(gridExtra)
library(mgcv)
library(sp)
source('yapsfunctions.R')
source('issmfunctions.R')
source("simtrack.R")
fmf <- wesanderson::wes_palette("FantasticFox1", type = "discrete")



#### load mods and data

# bring in sync data to extract hydro positions
sync <- readRDS('data/sync.RDS')
hydros <- data.table::data.table(sync$pl$TRUE_H)
colnames(hydros) <- c('hx','hy','hz')

# bring in the lake
lake <- readRDS("data/lake.rds")

# bring in time of arrival data (see yams_prepdata.R)
toalist <- readRDS('output/toalist.RDS')
# extract toa and burst interval vector
toa <- toalist$toa
bivec <- toalist$seq

# bring in the model results
mods21 <- readRDS("output/m2svhmm1.RDS")
mods22 <- readRDS("output/m2svhmm60.RDS")
mods31 <- readRDS("output/m3svhmm1.RDS")
mods32 <- readRDS("output/m3svhmm60.RDS")





####### pick the best mod for each group
# compare nlls
cbind(sapply(mods21, function(x)min(x$ssm_nll$nll)),
      sapply(mods22, function(x)min(x$ssm_nll$nll)))

mods2winners <- apply(cbind(sapply(mods21, function(x)min(x$ssm_nll$nll)),
            sapply(mods22, function(x)min(x$ssm_nll$nll))), 
      1, 
      which.min)
mods2winners

cbind(sapply(mods31, function(x)min(x$ssm_nll$nll)),
      sapply(mods32, function(x)min(x$ssm_nll$nll)))

mods3winners <- apply(cbind(sapply(mods31, function(x)min(x$ssm_nll$nll)),
            sapply(mods32, function(x)min(x$ssm_nll$nll))), 
      1,
      which.min)
mods3winners

allmods2 <- list(mods21, mods22)
allmods3 <- list(mods31, mods32)

# get the best mods for 2 and 3
mods2 <- list()
for(i in 1:5) mods2[[i]] <- allmods2[[mods2winners[i]]][[i]]

mods3 <- list()
for(i in 1:5) mods3[[i]] <- allmods3[[mods3winners[i]]][[i]]



# now extract both the location and behavioural states
stateslist <- list()
stateslist[[1]] <- data.frame(do.call(rbind, lapply(mods2, function(x)x$lstates))) %>% 
  mutate(b=do.call(c, lapply(mods2, function(x)x$bstates)),
         top=do.call(c, lapply(mods2, function(x)x$ssm_results[[x$winner]]$top)))
stateslist[[2]] <- data.frame(do.call(rbind, lapply(mods3, function(x)x$lstates))) %>% 
  mutate(b=do.call(c, lapply(mods3, function(x)x$bstates)),
         top=do.call(c, lapply(mods3, function(x)x$ssm_results[[x$winner]]$top)))
# check sizes
lapply(stateslist, dim)


states <- do.call(rbind, stateslist) %>% 
  mutate(modid=rep(c('mod2', 'mod3'), each=25000))
head(states)





################### RESULTS ################### 



# study system
lake %>% ggplot(aes(x=X, y=Y)) + 
  geom_polygon(colour="black", fill=NA) +
  geom_point(data=hydros, aes(x=hx, y=hy)) + 
  coord_fixed() + 
  theme_classic() + 
  xlab("Easting") + 
  ylab('Northing')




# location and behavioural states, all mods
states %>% 
  mutate(modid=factor(case_when(modid=='mod2' ~ 'Two-State Model',
                         T~'Three-State Model'),
                      levels=c('Two-State Model', 'Three-State Model'))) %>% 
  ggplot(aes(x=x, y=y)) + 
  geom_polygon(data=lake, aes(x=X, y=Y), colour="black", size=0.2, fill=NA) +
  geom_path() + 
  geom_point(aes(col=factor(b)), size=0.5) + 
  geom_point(data=hydros, aes(x=hx, y=hy), size=0.8) + 
  theme_classic() +
  scale_color_manual(values=c(fmf[2], fmf[3], fmf[5])) + 
  facet_wrap(~modid, ncol=2, scales='fixed') + 
  theme(legend.position='none',
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
  coord_fixed() +
  xlab('Eastings') + 
  ylab('Northings')







##############
# check out the winning iterations
lapply(mods2, function(x)x$winner)
lapply(mods3, function(x)x$winner)
# so with more states it becomes a bit more complex



##############
# get the values of the parameter D

sapply(mods2, function(x)x$hmm_results[[x$winner]]$params[7:8,1])
sapply(mods3, function(x)x$hmm_results[[x$winner]]$params[13:15,1])
# ouf, two is not super comparable



##############
# get the stationary distributions, activity budgets!

sapply(mods2, function(x)x$hmm_results[[x$winner]]$statdist)
sapply(mods3, function(x)x$hmm_results[[x$winner]]$statdist)


##############
# get the diagonals of the generator matrices 
# units mins

# these are the generator matrices
lapply(mods2, function(x)-(x$hmm_results[[x$winner]]$genmat))
lapply(mods3, function(x)-(x$hmm_results[[x$winner]]$genmat))

# check that the 3x3 generators are correct
lapply(mods3, function(x) rowSums(matrix(tail(x$hmm_results[[x$winner]]$params, 9)[,1], nrow=3)))
# okay, rows sum to 0 with a rounding error, that's good


# these are the diagonals
sapply(mods2, function(x)-(diag(x$hmm_results[[x$winner]]$genmat)))
sapply(mods3, function(x)-(diag(x$hmm_results[[x$winner]]$genmat)))

# this is the holding time
sapply(mods2, function(x)(1/-(diag(x$hmm_results[[x$winner]]$genmat))))
sapply(mods3, function(x)(1/-(diag(x$hmm_results[[x$winner]]$genmat))))






##############
# get and plot the qq plots

# single group examples
plot(mods2[[1]])
plotqq(mods2[[1]])


# pool all the results together
pseudoslist <- list()
pseudoslist[[1]] <- do.call(rbind, lapply(mods2, function(x)x$hmm_results[[x$winner]]$pseudo))
pseudoslist[[2]] <- do.call(rbind, lapply(mods3, function(x)x$hmm_results[[x$winner]]$pseudo))
lapply(pseudoslist, dim) # -1 for each 
pseudos <- do.call(rbind, pseudoslist) %>% 
  as.data.frame() %>% 
  rename(x=V1, y=V2) %>% 
  mutate(modidx=factor(rep(c('2 States', '3 States'), each=dim(pseudoslist[[1]])[1]))) %>% 
  pivot_longer(cols=c('x', 'y')) %>% 
  arrange(modidx, name) %>% 
  as.data.frame() 
pseudos <- pseudos %>% mutate(pseudoidx=paste(pseudos$modidx, pseudos$name))

# calculate the qqlines for the plots
qqlines <- do.call(rbind,
                   lapply(split(pseudos, f=paste(pseudos$modidx, pseudos$name)),
                          function(x)get_qqline(x$value))) %>%
  as.data.frame() %>% 
  mutate(pseudoidx = row.names(.)) %>% 
  rename(slope='75%', yint='25%') %>% 
  mutate(pseudoidx=case_when(pseudoidx=='2 States x' ~ '2 States; Eastings',
                             pseudoidx=='2 States y' ~ '2 States; Northings',
                             pseudoidx=='3 States x' ~ '3 States; Eastings',
                             pseudoidx=='3 States y' ~ '3 States; Northings')) %>% 
  mutate(name=c('Eastings', 'Northings', 'Eastings', 'Northings'),
         modidx=factor(c('Two-State Model', 'Two-State Model', 'Three-State Model', 'Three-State Model'),
                       levels=c('Two-State Model', 'Three-State Model')))



# now plot with facet grid
pseudos %>% 
  mutate(name=case_when(name=='x' ~ 'Eastings',
                        name=='y' ~ 'Northings'),
         modidx=factor(case_when(modidx=='2 States' ~ 'Two-State Model',
                          modidx=='3 States' ~ 'Three-State Model'), 
                       levels=c('Two-State Model', 'Three-State Model'))) %>%
  mutate(modidx=factor(modidx, 
                       levels=c('Two-State Model', 'Three-State Model'))) %>% 
  ggplot(aes(sample=value)) +
  geom_qq() +
  theme_bw() +
  geom_abline(data=qqlines, aes(slope=slope, intercept=yint), col=fmf[4], lwd=1.1) +
  facet_grid(modidx~name) +
  xlab('Theoretical Quantiles') + 
  ylab('Sample Quantiles') + 
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))








###############
# plot the behavioural states over time


# need to calculate the datetimes that each of the five groups starts
modidx <- rep(1:5, each=5000)
toamats <- split(data.table(toa), factor(modidx)) %>% lapply(as.matrix)
bivecs <- split(bivec, factor(modidx))
dates <- sapply(lapply(toamats, range, na.rm=TRUE),
                function(x)as.character(lubridate::as_datetime(x))) %>% 
  t() %>% 
  as.data.frame() %>% 
  rename(start=V1, stop=V2) %>% 
  mutate(start=lubridate::as_datetime(start),
         stop=lubridate::as_datetime(stop)) 
dates$avg <- lubridate::as_datetime(cbind(as.numeric(dates$start), as.numeric(dates$stop)) %>% rowMeans)
dates$grp <- c('Group 1', 'Group 2', 'Group 3', 'Group 4', 'Group 5')






#### plot the speeds vs the states
# checking the interpretability of the states that we've discovered
states <- do.call(rbind, stateslist) %>% 
  mutate(modid=rep(c('mod2', 'mod3'), each=25000))
head(states)
states$date <- c(dates$start[1] + states$top[1:5000], 
             dates$start[2] + states$top[5001:10000],
             dates$start[3] + states$top[10001:15000],
             dates$start[4] + states$top[15001:20000],
             dates$start[5] + states$top[20001:25000],
             dates$start[1] + states$top[25001:30000], 
             dates$start[2] + states$top[30001:35000],
             dates$start[3] + states$top[35001:40000],
             dates$start[4] + states$top[40001:45000],
             dates$start[5] + states$top[45001:50000])
head(states)

# claculate the distance
d <- do.call(c, 
             lapply(split(states, states$modid),
                    function(y)c(NA,apply(cbind(diff(y$x), diff(y$y)), 1, function(x)sqrt(sum(x^2)))))
)
# calculate the time difference by using the time of ping
deltat <- do.call(c, lapply(split(states, states$modid), function(x)c(NA, diff(x$top))))
# calculate the speed
v <- d/deltat
v[seq(1, 50000, by=5000)] <- NA # otherwise will be negative
# add to dataframe 
states <- states %>% mutate(deltat=deltat, 
                            d=d, 
                            v=v)

# plot the speed over time
states %>% 
  mutate(modid=factor(case_when(modid=='mod2' ~ 'Two-State Model',
                         modid=='mod3' ~ 'Three-State Model'), 
                      levels=c("Two-State Model",
                               'Three-State Model'))) %>% 
  ggplot(aes(x=date, y=v, col=factor(b))) + 
  ylim(-0.05, 1.2) + 
  xlim(dates$start[1], dates$stop[5]) +
  geom_vline(xintercept=c(dates$start[1:5], dates$stop[5]),
             linetype='dotted') + 
  geom_text(data=dates[1:5,], aes(x=avg, y=0.95, label=grp),
            inherit.aes=FALSE) + 
  geom_point() + 
  facet_wrap(~modid, nrow=2) +
  theme_bw() +
  scale_color_manual(values=c(fmf[2], fmf[3], fmf[5])) + 
  theme(legend.position='none',
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
  xlab('Time') + 
  ylab('Observed Speed (m/s)') 



# do a distribution plot
# basically same as above but density form instead of longitudinal
# had to scale the density for visual interpretability
states %>% 
  mutate(grp=case_when(grp==1 ~ 'Group 1',
                       grp==2 ~ 'Group 2',
                       grp==3 ~ 'Group 3',
                       grp==4 ~ 'Group 4',
                       grp==5 ~ 'Group 5'),
         modid=factor(case_when(modid=='mod2'~ 'Two-State Model',
                                modid=='mod3'~ 'Three-State Model'),
                      levels=c('Two-State Model', 'Three-State Model'))) %>% 
  ggplot(aes(x=v)) + 
  geom_histogram(aes(y=stat(count) / sum(count)), 
                 fill='slategray') +
  geom_density(aes(x=v, y=15* stat(count) / sum(count), col=factor(b)),
               position='stack',
               lwd=1) + 
  scale_color_manual(values=c(fmf[2], fmf[3], fmf[5])) + 
  facet_grid(modid~ grp) + 
  theme_bw() + 
  theme(legend.position='none',
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
  xlab('Observed Velocity (m/s)') + 
  ylab('Density')




#########
# check pike speeds within each state within each group
# min and max in table form 

states$grp <- rep(rep(1:5, each=5000), 2)

states %>% 
  group_by(modid, b, grp) %>% 
  summarize(min=round(min(v, na.rm=TRUE),3), max=round(max(v, na.rm=TRUE),3)) %>% 
  as.data.frame() %>% 
  arrange(modid, grp, b)


# interquartile range instead
states %>% 
  group_by(modid, b, grp) %>% 
  summarize(min=quantile(v, 0.25, na.rm=TRUE), max=quantile(v, 0.75, na.rm=TRUE)) %>% 
  as.data.frame() %>% 
  arrange(modid, grp, b)
  

##########
# check out the number of switches for each model
# subtract one for the number of switches
sapply(mods2, function(x)length(rle(x$bstates)$values)-1)
sapply(mods3, function(x)length(rle(x$bstates)$values)-1)
# more switches for three states














