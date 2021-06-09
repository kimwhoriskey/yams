

gentrack <- function(startseed, 
                     nrowhydros = 6, ncolhydros = 6, 
                     mean_hydro_dist = 1,
                     hydropos=NULL, poly=NULL, 
                     nloc,  
                     bi, 
                     demod, #demodlink, 
                     move='wiener',
                     par, 
                     me = "t",
                     multipath=NULL,
                     sos=1465, 
                     timescale=1, 
                     case=2){
 
  
  set.seed(startseed)
  
  
  ##########################
  ### time of pings part ###
  ##########################
  
  # first simulate the burst interval sequence bivec
  # as per manufacturer it should be U(min, max)
  # it's the tag, not the receiver
  bivec <- runif(nloc-1, bi[1], bi[2])
  # top_rkbi = c(0, cumsum(bivec))/timescale
  top_rkbi = c(0, cumsum(bivec))
  
  
  ##########################
  ######  behaviours  ######
  ##########################

  
  #continuous time
  # alpha is based off of a different time scale, so we're simulating these switches
  # based on that time scale
  # the time scale it's based off of is called "timescale"
  m = nrow(par$alpha)
  B <-  factor(1:m)
  bt <- numeric()
  trueb <- sample(B, size=1, prob=rep(1/m, m))
  truebt <- rexp(1, rate=-par$alpha[trueb[1], trueb[1]]) #the mod is scaled to mins, but bivec is in secs, going to have to generalize this
  i=2
  while(sum(truebt)<sum(bivec/timescale)){
    probability <- par$alpha[as.numeric(trueb[i-1]), -as.numeric(trueb[i-1])]/sum(par$alpha[as.numeric(trueb[i-1]), -as.numeric(trueb[i-1])])
    thing <- B[B!=trueb[i-1]]
    trueb <- append(trueb, sample(thing, size=1, 
                          prob=probability))
    truebt <- append(truebt, rexp(1, rate=-par$alpha[as.numeric(trueb[i]), as.numeric(trueb[i])]))
    i=i+1
  }
  # do one more past
  probability <- par$alpha[as.numeric(tail(trueb, 1)), -as.numeric(tail(trueb[i-1],2)[1])]/sum(par$alpha[as.numeric(tail(trueb[i-1],2)[1]), -as.numeric(tail(trueb[i-1],2)[1])])
  thing <- B[B!=tail(trueb[i-1],2)[1]]
  trueb <- append(trueb, sample(thing, size=1, 
                                prob=probability))
  truebt <- append(truebt, rexp(1, rate=-par$alpha[as.numeric(tail(trueb, 1)), as.numeric(tail(trueb, 1))]))
 

  
  
  
  #####################################################
  ######  simulate tag drift and observed times  ######
  #####################################################
  

  # next simulate the tag drift
  # this follows a random walk
  tagdrift <- numeric(nloc) # in seconds
  # tagdrift[1] <- rnorm(1, 0, 4)
  # tagdrift[2] <- rnorm(1, tagdrift[1], 4)
  tagdrift[1] <- rnorm(1, 0, par$sigmabi)
  tagdrift[2] <- rnorm(1, tagdrift[1], par$sigmabi)
  # this second tagdrift has a big effect, if you put this to 4 like in the 
  # likelihood then the shape becomes very linear
  for(i in 3:nloc) tagdrift[i] <- 2*tagdrift[i-1] - tagdrift[i-2] + rnorm(1, 0, par$sigmabi)
  # next calculate the time of the ping of the random burst interval
  top <- top_rkbi + tagdrift # top_rkbi already has been divided by the timescale
  # tag drift is estimated on the scale of the observed times
  
  
  
  
  
  ##############################
  ######  true locations  ######
  ##############################
  
  # get full states for movement sim
  # case 1 - not accounting for snapshot principle, so figure out the behaviours at the 
  # observation times
  # this should line up with top_rkbi, not top, because top includes measurement error (tag drift)
  # case 2 - accounts for snapshot principle, 
  # augment the other simulation by adding in times and locationssimulate locations under trueb and then linearly interpolate
  # using top right now because that is consistent with the likelihood
  
  behavswitches <- cumsum(truebt)*timescale # on the scale of timescale, so multiply by timescale to get to observed scale
  
  if(case==1){
    # use -5 because of a teeny amount of tag drift
    intidx <- findInterval(top, c(-5,behavswitches)) # this should maybe be top_rkbi? 
    b <- trueb[intidx]
  } else if(case==2){
    fulltimes <- sort(c(top, behavswitches))
    intidx <- findInterval(fulltimes, c(-5,behavswitches))
    # when it changes, the behavswitches time is also included in the switch, which is good
    b <- trueb[intidx]
    b[length(b)] <- b[length(b)-1] # last b will be NA because it is the last switch
    # we remove this b and X anyways so i'll just set it to the previous value
    lateridx <- which(fulltimes %in% behavswitches) # locations to remove later
  }
  
  
  # pick times based on cases
  # times here need to be scaled based on the time scale
  if(case == 1){
    times = top/timescale # could do top_rkbi
  } else if (case==2){
    times = fulltimes/timescale
  }
  
  ### true locations part
  X <- matrix(0, nrow=length(b), ncol=2)
  # X[1,] <- mvtnorm::rmvnorm(1, c(0,0), 
  #                           matrix(c(10000,0,0,1000), nrow=2))
  if(move=="Wiener"){
    for(i in 2:length(b)) {
      X[i,] <- mvtnorm::rmvnorm(1, 
                                X[i-1,], 
                                diag(2*par$Dxy[b[i]]*(times[i] - times[i-1]), 2))
    }
  } 

  # now remove the points from the behavioural switches if it's case 2
  if(case==2){
    X <- X[-lateridx, ]
    b <- b[-lateridx]
  } 
  
 
  # need to simulate a hydrophone area that is relevant to the 
  # detection efficiency
  # take the mean distance of all pairs of hydrophones from original data
  # two options here
  # first, we can assume a regular grid that is spaced similarly to original hydro files
  # second, if we input a shapefile and hydrophone positions then we can simulate data within 
  # the shapefile, and the hydrophone positions are still needed for calculating the 
  # detection efficiency
  # HYDROPHONE POSITIONS AND SHAPEFILE NEED TO BE TRANSLATED BACK TO 0
  if(is.null(hydropos)){
    # regular grid
    hydros <- expand.grid(x=(1:ncolhydros)*mean_hydro_dist - (ncolhydros+1)*mean_hydro_dist/2, 
                          y=(1:nrowhydros)*mean_hydro_dist  - (nrowhydros+1)*mean_hydro_dist/2) 
  } else {
    # hydrophones, needed below for lhdists
    hydros <- hydropos
    # shapefile, ideally 
    sX <- sp::SpatialPoints(X)
    plot(poly)
    points(sX, col='blue')
    sXnew <- sX
    while(anyNA(over(sXnew, poly))){
      sXnew <- sp::SpatialPoints(sX@coords + sp::spsample(poly, n=1, type='random')@coords[rep(1,nrow(sX@coords)),])
    } 
    points(sXnew, col='seagreen1')
    X <- sXnew@coords
  }

  # now detection efficiency
  # need the distances between the track and all simulated receivers at any point in time
  # lhdists => location-hydro dists
  # eassiest to just predict based on the detection efficiency mod
  lhdists <- outer(X=1:nloc, Y=1:nrow(hydros), 
                   FUN=function(loc,h) { 
                     sqrt(  (hydros[h,]$x - X[loc,1])^2  +  (hydros[h,]$y - X[loc,2])^2 )  
                     }
                   )
  
   newdat <- data.frame(as.numeric(lhdists))
   names(newdat) <- names(demod$model)[2]
   p <- predict(demod, newdata=newdat, type="response") %>% matrix(nrow=nrow(lhdists), ncol=ncol(lhdists))
  
  # anything over 1000, the probability should automatically be assumed to be 0
  p[which(lhdists>1000)] <- 0

  # sample based on a binomial dist
  dets <- matrix(0, nrow=nrow(lhdists), ncol=ncol(lhdists))
  for(i in 1:nrow(lhdists)){
    for(j in 1:ncol(lhdists)){
      dets[i,j] <- rbinom(1, 1, p[i,j])
    }
  }

  # HOLY SHIT WE DID IT

  # measurement error distribution comes into play with the temporal differences between the 
  # time of ping and time of arrival 
  # can use a gaussian or t
  if(me == "gau"){
    toa <- (top + lhdists/sos + rnorm(length(lhdists), 0, scale)) * ifelse(dets==1, 1, NA)
  } else if(me == "t"){
    toa <- (top + lhdists/sos + par$scale*rt(length(lhdists), df=3)) * ifelse(dets==1, 1, NA)
  }
  # to account for multipaths, i.e., when the pings are bouncing off rocks etc, 
  # we can add some extra noise
  if(!is.null(multipath)) toa <- toa + runif(length(lhdists), 50, 500) * rbinom(length(lhdists), 1, multipath) / sos
 
  rslt <- list(args = match.call(), 
               states = b, 
               locs = X, 
               top = list(bivec = bivec, 
                          top_rkbi = top_rkbi,
                          tagdrift = tagdrift, 
                          top = top), 
               hydros = hydros, 
               deteff = list(lhdists = lhdists, p = p),
               dets = dets,
               toa = toa) 
  class(rslt) <- "yapssim"
  
  return(rslt)
}  




plot.yapssim <- function(x){
  
  plot(x$locs, type='l', asp=1, 
       main = "simulated track with states",
       xlab = "easting", ylab = "northing")
  points(x$locs, pch = 20, col=ifelse(x$states==2, "tomato", "cadetblue"))
  
  plot(x$hydros, col="navy", pch=20, asp=1, 
       main = "simulated track with hydros",
       xlab = "easting", ylab = "northing")
  points(x$locs, type='l')
  
  # plot(x$deteff$p~ x$deteff$lhdists, ylim = c(0, 1), 
       # main = "detection efficiency", xlab = 'distance (m)', ylab = "detection probability")
  
}

# saveRDS(de_m, file='de_m.rda')


# plot(sim)

# sim <- gentrack(startseed = 4242, 
#          nrowhydros = 6, ncolhydros = 6, 
#          mean_hydro_dist = mean_hydro_dist,
#          nloc = 2000, 
#          bi = c(30,60),
#          decoeff = fixef(de_m), 
#          par = list(alpha = matrix(c(0.95, 0.05, 0.03, 0.97), nrow=2, byrow=TRUE),
#                     sigmabi = 2.260329e-06,
#                     Dxy = c(0.09, 0.001),
#                     scale = exp(pl$logScale)),
#          me="t",
#          sos=1465, multipath=NULL)
# 
# 

# 
# 
# plot(top, type='l', col='cadetblue', lwd=3)
# points(top_rkbi, type='l', col='tomato', lwd=3) # not gonna see the difference
# 
# 
# 
# simdets_mixed <- matrix(rbinom(length(p), 1, p), nrow=nrow(p))
# nobs <- apply(simdets_mixed, 1, function(k) sum(k!=0))
# plot(nobs)
# lines(caTools::runmean(nobs, k=100), col="blue")
# 
# plot(simhydros, type='p', col='blue')
# points(X, type='l')
# points(X, type="p", cex=0.7, col=ifelse(b=="2", "tomato", "cadetblue"))
# apply(simdets, 1, FUN=function(k)sum(k==1)) %>% plot()
# 
# 
# 
#   # dim(toa)  
#   # plotToa(toa)
#   # plotToa(toa_mp)
#   # plotNobs(toa)
#   # 
# 
#   
#   
#   
#   
#   
#   