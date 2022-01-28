# bootstrap the 2 behaviour model 
#### set up working directory

setwd("~/Desktop/yams")
rm(list=ls())	
library(data.table)
library(ggplot2)
library(yaps)
library(TMB)
library(tidyverse) 
library(gridExtra)
library(mgcv)
library(sp)
source('yapsfunctions.R')
source('issmfunctions.R')
source("simtrack.R")
fmf <- wesanderson::wes_palette("FantasticFox1", type = "discrete")

# TMB functions
compile("yaps_ssm.cpp")
dyn.load(dynlib("yaps_ssm"))
compile("yaps_hmm.cpp")
dyn.load(dynlib("yaps_hmm"))


#### load mods and data

#### load mods and data

# bring in sync data to extract hydro positions
sync <- readRDS('data/sync.RDS')
hydros <- data.table::data.table(sync$pl$TRUE_H)
colnames(hydros) <- c('x','y','z')

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

# need this too
rbi_min <- 10
rbi_max <- 30


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









# pick a single mod to calculate the detection efficiency
# we'll pick the first group of the 2 state model
mod <- mods2[[1]]
modidx <- rep(1:5, each=5000)
toamats <- split(data.table(toa), factor(modidx)) %>% lapply(as.matrix)
bivecs <- split(bivec, factor(modidx))

#########################################################
########### estimate the detection efficiency ###########
#########################################################

# we need to estimate the detection efficiency for the simulation study
# i explored several options, this one worked by far the best
# to get the data in the right format, translate toa into 0s and 1s


toa01 <- ifelse(is.na(toamats[[1]]), 0, 1) # this is only for one subset
dets <- reshape2::melt(toa01 %>% t()) %>% rename(hydro=Var1, pingidx=Var2, det=value)
# var1 is row of original data frame / matrix (in our case hydro)
# var 2 is column of original data frame / matrix (in our case time)
# value is the value of the element
obj <- mod$ssm_results[[mod$winner]]$obj
# get the distances from the tmb object
dets$dist <- reshape2::melt(obj$report()$dist)$value
dets <- data.table(dets)


#######
# basic approach is to bin the data and calculating detection probabilities empirically
# then fit a binomial glm with weights
# use 5 m intervals
# assuming detection probability is constant through space and time
distidx <- seq(0, max(dets$dist), by=5)
dets$distidx <- findInterval(dets$dist, distidx)*5 # basically rounding to the nearest five
deteff <- dets %>% group_by(distidx) %>% arrange(distidx) %>% summarize(eff = mean(det), totdet = n())
head(deteff)
# find the largest distance without a 0 value
max(deteff[deteff$eff!=0,]$distidx)
# everything above 930 is = 0 (deteff that is) 
# we'll assume everything over 1000m has a zero probability
# that is pretty standard for acoustic telemetry
which(deteff$eff==1)
# we also have two rows with efficiency = 1
deteff <- deteff[deteff$distidx < 1000,]
deteff %>% 
  ggplot(aes(x=distidx, y=eff, col=totdet)) + 
  geom_point() + 
  xlim(0, 1000) + 
  theme_bw()
# fit the binomial with weights
de <- gam(eff~ s(distidx), data = deteff, family=binomial, weights = totdet)
deteff %>% ggplot(aes(x=distidx, y=eff, col=totdet)) +
  geom_point() +
  geom_line(aes(x=distidx, y=fitted(de)), col='tomato') +
  theme_bw()

# quick cross validation, sample 70% in inidx, the rest go into outidx
inidx <- sort(sample(1:nrow(deteff), size=round(.7*nrow(deteff))))
outidx <- (1:nrow(deteff))[!(1:nrow(deteff)) %in% inidx]
detefftrain <- deteff[inidx,]
detrain <- gam(eff~ s(distidx), data = detefftrain, family=binomial, weights = totdet)
outpreds <- predict(detrain, newdata=deteff[outidx,], type="response")
predpts <- rbind(
  cbind(deteff[outidx,c('distidx', 'eff')], type=rep("Observed", length(outidx))),
  cbind(deteff[outidx,'distidx'], eff=outpreds, type=rep("Predicted", length(outidx)))
)

# plot the training/testing results
detefftrain %>% ggplot(aes(x=distidx, y=eff, col=totdet)) +
  geom_point() +
  scale_color_gradient(low=fmf[3], high="navyblue", name="Total No. of\nDetections") +
  geom_line(aes(x=distidx, y=fitted(detrain)), col=fmf[5]) +
  ggnewscale::new_scale("col") + #whoa ggnewscale FOR THE WIN, man that was a b
  geom_point(data=predpts, aes(x=distidx, y=eff, col=factor(type)), size=2, inherit.aes = FALSE) +
  scale_color_manual(values=c(fmf[4], fmf[5]), name='Testing Data') +
  xlab("Distance (m)") +
  ylab("Detection Efficiency") +
  theme_bw()


# efficiency at 250 and 500 m
predict(de, newdata=data.frame(distidx=c(100, 250, 500, 750)), type='response')
# really cool
# these are the same regardless of using mod2 or mod3

















########################################
# run  simulation studies 



# need the lake data because we want to place the sim dat randomly within the lake
# have to input the centered hydros and lake to match the predictions from yaps
centeredlake <- lake[, c("X", "Y")] - hydros[, c("x", "y")][rep(1, nrow(lake))]
centeredhydros <- hydros[, c("x", "y")] - hydros[, c("x", "y")][rep(1, nrow(hydros))]
lakepoly <- Orcs::coords2Polygons(as.matrix(centeredlake), ID="lake")




# pick mod to base sim on
mod <- mods2[[1]]
# starting values: use second choice because it was better for the first set
boot <- bootstrap(mod,
                  startseed=42, # started at 42, so now do 72
                  nsims=60,
                  savepath='output/twostatesboot',
                  useSSMds = FALSE,
                  dstart=60,
                  simargs =  list(
                    hydropos = centeredhydros,
                    poly=lakepoly,
                    nloc = 5000,
                    bi = c(10,30),
                    demod = de,
                    me="t",
                    sos=1465,
                    multipath=NULL,
                    move="Wiener",
                    case=2,
                    timescale=60),
                  issm.args = list(maxsteps=5,
                                   fixsteps=TRUE,
                                   allowfc=FALSE,
                                   fcmax=3,
                                   jiggle_fc=0.01,
                                   initssmargs = list(mapping = list(working_A = factor(matrix(NA)))),
                                   ssmargs = list(mapping = list(working_A = factor(matrix(NA, nrow=2, ncol=2))),
                                                  optimizer='nlminb'),
                                   timescalehmm=60,
                                   logdstart = c(log(60), log(60)),#c(1,3),
                                   setinitstatdist=1)
)


mod <- mods3[[1]]
# starting values: use second choice because it was better for the first set
boot <- bootstrap(mod,
                  startseed=42,
                  nsims=60,
                  savepath='output/threestatesboot',
                  dstart=60,
                  simargs =  list(
                    hydropos = centeredhydros,
                    poly=lakepoly,
                    nloc = 5000,
                    bi = c(10,30),
                    demod = de,
                    me="t",
                    sos=1465,
                    multipath=NULL,
                    move="Wiener",
                    case=2,
                    timescale=60),
                  issm.args = list(maxsteps=7,
                                   m=3,
                                   fixsteps=TRUE,
                                   allowfc=FALSE,
                                   fcmax=3,
                                   jiggle_fc=0.03,
                                   initssmargs = list(mapping = list(working_A = factor(matrix(NA)))),
                                   ssmargs = list(mapping = list(working_A = factor(matrix(NA, nrow=3, ncol=3))),
                                                  inner_control = list(maxit=1500),
                                                  optimizer='nlminb'),
                                   timescalehmm=60,
                                   logdstart = c(log(60), log(60), log(60)),#c(1,3),
                                   setinitstatdist=1)
)




##################################
# look at simulation results


boot2 <- readRDS('output/twostatesboot.RDS')
boot3 <- readRDS('output/threestatesboot.RDS')






# check if any didn't work
sapply(boot2$mods, length) %>% length
sapply(boot3$mods, length) %>% length
# only 57 worked

# check convergence
sapply(boot2$mods, function(x)x$ssm_results[[x$winner]]$mess)
sapply(boot2$mods, function(x)x$ssm_results[[x$winner]]$opt$convergence)
fcidx2 <- which(sapply(boot2$mods, function(x)x$ssm_results[[x$winner]]$opt$convergence)>0)
length(fcidx2)
length(boot2$mods)-length(fcidx2)

sapply(boot3$mods, function(x)x$ssm_results[[x$winner]]$mess)
sapply(boot3$mods, function(x)x$ssm_results[[x$winner]]$opt$convergence)
fcidx3 <- which(sapply(boot3$mods, function(x)x$ssm_results[[x$winner]]$opt$convergence)>0)
length(fcidx3)
length(boot3$mods)-length(fcidx3)





# summary stats
summary(1-boot2$err.rate[-fcidx2])
summary(boot2$lon.rmse[-fcidx2])
summary(boot2$lat.rmse[-fcidx2])
summary(boot2$rmse)

summary(1-boot3$err.rate[-fcidx3])
summary(boot3$lon.rmse[-fcidx3])
summary(boot3$lat.rmse[-fcidx3])
summary(boot3$rmse)



# get bootstrap results ready for plotting 
bootdat <- data.frame(val=c(1-boot2$err.rate[-fcidx2], 1-boot3$err.rate[-fcidx3],
                            boot2$lon.rmse[-fcidx2], boot3$lon.rmse[-fcidx3],
                            boot2$lat.rmse[-fcidx2],boot3$lat.rmse[-fcidx3],
                            boot2$rmse[-fcidx2], boot3$rmse[-fcidx3]),
                      type=rep(c('bstateacc', 'Eastings', 'Northings', 'locstateacc'), each=length(boot2$err.rate[-fcidx2])+length(boot3$err.rate[-fcidx3])),
                      m=factor(rep(c('Two-State Model', 'Three-State Model'), times=c(length(boot2$err.rate[-fcidx2]),length(boot3$err.rate[-fcidx3]))),
                               levels=c('Two-State Model', 'Three-State Model')))




##########################
# look at state accuracy

bootdat %>% 
  filter(type=='bstateacc') %>% 
  ggplot(aes(x=m, y=val, group=m, fill=m)) + 
  geom_boxplot() + 
  ylim(0.6, 1) + 
  theme_bw() + 
  ylab('Proportion') + 
  ggtitle('Behavioural State Accuracy') + 
  scale_fill_manual(values=c('slategray1', 'cadetblue')) + 
  theme(legend.position = 'none', 
        axis.title.x=element_blank(),
        plot.title = element_text(hjust = 0.5))

# location state accuracy in each axis
bootdat %>% 
  filter(type %in% c('Eastings', 'Northings')) %>% 
  ggplot(aes(x=m, y=val, fill=m)) + 
  geom_boxplot() + 
  theme_bw() + 
  ylab('Root Mean Squared Error (m)') + 
  xlab('Location Axis') + 
  ggtitle('Location State Accuracy') + 
  facet_wrap(~type) + 
  scale_fill_manual(values=c('slategray1', 'cadetblue')) + 
  theme(legend.position = 'none',
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.title.x=element_blank(),
        plot.title = element_text(hjust = 0.5))

# location state accuracy total
bootdat %>% 
  filter(type=='locstateacc') %>% 
  ggplot(aes(x=m, y=val, fill=m)) + 
  geom_boxplot() + 
  theme_bw() + 
  ylab('Root Mean Squared Error (m)') + 
  xlab('Location Axis') + 
  ggtitle('Location State Accuracy') + 
  scale_fill_manual(values=c('slategray1', 'cadetblue')) + 
  theme(legend.position = 'none', 
        axis.title.x=element_blank(),
        plot.title = element_text(hjust = 0.5))





# let's look at the winners
sapply(boot2$mods, function(x)x$winner)
sapply(boot3$mods, function(x)x$winner)
# don't need many iterations  




##################################
# let's look at the parameters now 

boot2$true.pars
boot2$pars



# lets do some tables
boot2$pars %>% 
  select(-par, -true) %>% 
  select(-fcidx2) %>% 
  mutate(mean=rowMeans(.), 
         sd=sqrt(apply(., 1, var)),
         lower2.5=apply(., 1, quantile, 0.025),
         upper97.5=apply(., 1, quantile, 0.975),
         median=apply(., 1, median)) %>% 
  mutate(true=boot2$true.pars,
         par=as.character(boot2$pars$par)) %>% 
  select(par, true, mean, median, sd, lower2.5, upper97.5) %>% 
  xtable::xtable(digits=3) %>% 
  print(include.rownames=FALSE)


boot3$pars %>% 
  select(-par, -true) %>% 
  select(-fcidx3) %>% 
  mutate(mean=rowMeans(.), 
         sd=sqrt(apply(., 1, var)),
         lower2.5=apply(., 1, quantile, 0.025),
         upper97.5=apply(., 1, quantile, 0.975),
         median=apply(., 1, median)) %>% 
  mutate(true=boot3$true.pars,
         par=as.character(boot3$pars$par)) %>% 
  select(par, true, mean, median, sd, lower2.5, upper97.5) %>% 
  xtable::xtable(digits=3) %>% 
  print(include.rownames=FALSE)


  
  







