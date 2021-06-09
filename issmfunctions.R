##### issm fitting functions



###############################################################
##################### First Step: Fit HMM #####################
###############################################################

fit_hmm <- function(dat,
                    p_start = list(logD_xy = rep(log(1), 2),
                                   working_A = matrix(rep(log(1),2), nrow=2)),
                    silence=TRUE){
  
  # create TMB object
  obj <- MakeADFun(dat, p_start, DLL="yaps_hmm", silent=silence)
  
  # optimize TMB object
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  
  # likelihood value
  nll <- obj$fn(opt$par)
  
  # extract behavioural states
  b_states <- obj$report()$states
  
  # extract parameter estimates
  params_hmm <- summary(sdreport(obj))
  
  # get the stationary distribution
  stationary <- obj$report()$delta
  
  # get the generator matrix
  m <- length(p_start$logD_xy)
  genmat <- matrix(tail(params_hmm, m^2)[,1],nrow=m)
  
  # get the pseudoresiduals
  pseudo <- qnorm(obj$report()$pseudo[-1,])

  # return results
  rslt <- list(mess = opt$mess, obj=obj, opt=opt, nll=nll, params = params_hmm, 
               b_hat=b_states, statdist = stationary, genmat=genmat, pseudo = pseudo)
  return(rslt)
  
}














##############################################################################################
##################### Second Step: Fit SSM with fixed behavioural states #####################
##############################################################################################




fit_ssm <- function(dat, 
                    p_start,
                    hydro_translation_factors,
                    mapping=NULL, 
                    inner_control = list(maxit=500),
                    allowErr=FALSE, 
                    jiggle_err=0.01, 
                    optimizer="nlminb", 
                    optimizer_arguments=NULL,
                    silence=TRUE){

  #extract the names of the fixed variables
  fixed_names <- names(mapping)
  
  
  # set up an error vector so that we can keep track of how many errors occur
  err_num = 0 # set initial error value at 0, base an iterative loop on this
  err = 0
  class(err) <- "error"
  
  while("error"%in%class(err) & allowErr==FALSE){
    
    if(err_num>10) stop("encountered errors > 10 times while attempting to fit the SSM") #cut optimization if errors > 10
    
    if(err_num>0){

      # print a message letting us know we're working on errors
      cat("\n -------------------- \n working on errors \n err: ", err$message, " \n -------------------- \n")
      
      # jiggle starting parameters - add a little bit of random error to the ones that aren't fixed
      # exclude the random effects as well
      idx <- which(!(names(p_start) %in% c("x", fixed_names)))
      for(i in idx) p_start[[i]]<- p_start[[i]] + rnorm(length(p_start[[i]]), 0, jiggle_err)
    }

    # re-create TMB object, remove previous just in case
    if(exists("obj")) rm(obj)
    obj <- MakeADFun(dat, p_start, map=mapping, random=c("X", "Y", "tag_drift"), DLL="yaps_ssm",
                     inner.control = inner_control, silent=silence)
    argus <- append(list(obj$par, obj$fn, obj$gr), optimizer_arguments)
    
    # use trycatch loop again
    err <- tryCatch({
      opt <- do.call(optimizer, argus)
      srep <- summary(sdreport(obj)) # extract parameter estimates
    },  error=function(e){e} #capture errors
    )
    if("error" %in% class(err)) err_num = err_num+1 #count the errors
  } # close while loop
  
  
  # extract the likelihood value
  nll <- obj$fn(opt$par)
  
  # parameter estimates, take the random effects out of it
  params_ssm <- srep[!(rownames(srep) %in% c("X", "Y", "tag_drift")),]
  
  # extract location state estimates
  # if statements just in case we fix them
  if(!("X" %in% fixed_names)) x <- srep[row.names(srep) == "X",]
  if(!("Y" %in% fixed_names)) y <- srep[row.names(srep) == "Y",]
  if(!("tag_drift" %in% fixed_names)) td <- srep[row.names(srep) == "tag_drift", 1]
  
  if(!("X" %in% fixed_names)){
    l_states = data.frame(x=x[,1]+hydro_translation_factors[1], 
                          y=y[,1]+hydro_translation_factors[2])
  } else {
    l_states = NULL
  }
  
  top <- obj$report()$top
  
  
  rslt <- list(mess = opt$mess, obj=obj, opt=opt, nll=nll, 
               params=params_ssm, l_hat=l_states, x=x, y=y, td_hat=td, top=top)
  return(rslt)
  
}







######################################################################################
##################### Iteratively fit multiple SSM and HMM steps #####################
######################################################################################







fit_issm <- function(inp_init, 
                     inp_ssm, 
                     maxsteps=10, 
                     m=2, 
                     fixsteps=FALSE, 
                     allowfc=FALSE, 
                     fcmax=10,
                     jiggle_fc=0.01, 
                     logdstart = rep(log(1), 2),
                     initssmargs = list(mapping = list(working_A = factor(matrix(NA)))), 
                     ssmargs = list(mapping = list(working_A = factor(matrix(NA, 
                                                              nrow=nrow(p_start_ssm$working_A),
                                                              ncol=ncol(p_start_ssm$working_A))))), 
                     # ssm_map=NULL,
                     allsilent=TRUE,
                     timescalehmm=1,
                     optimizer="nlminb",
                     setinitstatdist=1){
  
 
  

  # set up accumulators
  # all of the results from each step
  hmm_results <- list()
  ssm_results <- list()
  # times it takes to fit each step
  hmm_time <- list() 
  ssm_time <- list()
  # convergence of the smm, cue a while loop on this if we don't allow false convergence
  ssm_conv <- numeric()
  fc_count <- numeric() # store the false convergence counts (how many at each step)
  
  # store the nll of the ssm, cue another while loop on this
  nll_ssm <- numeric()
  

 
  ###############################################################################
  ############ initial step, fit one behaviour YAPS to observed data ############
  ###############################################################################
  
  initmodargs <- append(list(dat=inp_init$datTmb,
                             p_start = inp_init$params,
                             hydro_translation_factors = c(inp_init$inp_params$Hx0,
                                                           inp_init$inp_params$Hy0),
                             silence=allsilent),
                        initssmargs)
  init_mod <- do.call(fit_ssm, initmodargs)

  cat("\n ------------------------------------------------------------ \n finished initial SSM, convergence was:",
      init_mod$mess, 
      "\n ------------------------------------------------------------ \n ")
  
  ################################################################
  ############ first step, fit HMM to observed data ############
  ################################################################

  
  # fit the hmm
  p_start_hmm = list(logD_xy = logdstart,
                     working_A = matrix(rep(log(1),m*m), nrow=m))
  hmm_time[[1]] <- system.time(hmm_results[[1]] <- fit_hmm(dat = list(X = init_mod$x[,1], 
                                                                      Y = init_mod$y[,1], 
                                                                      top = init_mod$top,
                                                                      timescale=timescalehmm, 
                                                                      setinitstatdist=setinitstatdist),
                                                           p_start = p_start_hmm,
                                                           silence=allsilent))
    # print out an update
  cat("\n ------------------------------------------------------------ \n finished HMM 1, convergence was:",
      hmm_results[[1]]$mess, 
      "\n ------------------------------------------------------------ \n ")
  
  ##############################################################################
  ########### second step, fit SSM to observed data with fixed states ###########
  ##############################################################################


  
  # get starting values for the parameters
  p_start_ssm <- inp_ssm$params
  # overwrite the A
  p_start_ssm$working_A = matrix(hmm_results[[1]]$params[row.names(hmm_results[[1]]$params) %in% 'working_A',1],nrow=m)
  # overwrite the logD_xy if they're fixed
  if("logD_xy" %in% names(ssmargs$mapping)){
    p_start_ssm$logD_xy <- hmm_results[[1]]$params[row.names(hmm_results[[1]]$params) == "logD_xy", 1] 
  }
  # doesn't work if we overwrite the logD_xy, should build this back in
  

  
  # change over the ssm_dat behavioural states
  ssm_dat <- inp_ssm$datTmb
  ssm_dat$b <- hmm_results[[1]]$b_hat-1
  ssm_dat$matexp <- hmm_results[[1]]$obj$report()$matexp

  # ssm_map <- list(logD_xy = rep(factor(NA),2) )
  # ssm_map <- list(working_A = factor(matrix(NA, 
  #                                           nrow=nrow(p_start_ssm$working_A),
  #                                           ncol=ncol(p_start_ssm$working_A))))
  # put together the argument list for the model
  ssmargs <- append(list(dat = ssm_dat,
                         p_start = p_start_ssm, 
                         hydro_translation_factors = c(inp_ssm$inp_params$Hx0,
                                                       inp_ssm$inp_params$Hy0),
                         # mapping = ssm_map,
                         silence=allsilent),
                    ssmargs)
  
  # exclude the fixed parameters from the jiggle 
  ssm_fixed_names <- names(ssmargs$mapping) # record the names of the fixed parameters
  pidx <- which(!(names(p_start_ssm) %in%  c("x", ssm_fixed_names))) #create an index to loop over, p for parameter

  
  # use a while loop acting on whether we allow fc and what the count of the false convergence is
  fc_count[1]=0
  ssm_conv[[1]]=0 #set it to 1 to initialize the loop
  runme = TRUE
  while(runme){
    
    if((sum(ssm_conv[[1]])>0) & (allowfc==FALSE) & (fc_count[1]>0)){
      # print a message letting us know we're working on false convergence
      cat("\n --------------------------- \n working on false convergence \n --------------------------- \n")
    
    # jiggle the starting parameters, jiggle_fc is the st dev for random normal error
    for(i in pidx) ssmargs$p_start[[i]] <- ssmargs$p_start[[i]] + rnorm(length(ssmargs$p_start[[i]]), 0, jiggle_fc)
    }
    
    # fit the ssm 
    ssm_time[[1]] <- system.time(ssm_results[[1]] <- do.call(fit_ssm, ssmargs))
                                                             
    # save the new convergence, cue on this
    ssm_conv[1] <- ifelse(ssmargs$optimizer=="nloptr", 
                          ssm_results[[1]]$opt$status, 
                          ssm_results[[1]]$opt$convergence)
    # keep track of nll if we use fixsteps=FALSE to tell the optimizer to stop right after the nll starts to increase
    nll_ssm[1] = ifelse(ssmargs$optimizer=="nloptr", 
                        ssm_results[[1]]$opt$objective, 
                        ssm_results[[1]]$nll)

    fc_count[1] = fc_count[1] + 1 # add to the count every time we come back up to the beginning of this loop
    if(fc_count[1]>=fcmax) warning(paste("had false convergence ", fcmax, "times while trying to get the SSM to fit"))
    # stop if greater than 10 times
    
    if(sum(ssm_conv[[1]])==0 | allowfc==TRUE | fc_count[1]>=fcmax) runme=FALSE
    
  } # close while loop
  

  # update on where we are
  cat("\n ------------------------------------------------------------ \n finished SSM 1, convergence was:",
      ssm_results[[1]]$mess, 
      "\n ------------------------------------------------------------ \n ")
  
  
  ###########################################
  ############ rest of the steps ############
  ###########################################
  
  # initialize the indexer for the rest of the iterations
  i=2
  
  # use a while loop, and if fixsteps = false, then stop the optimization when nll_diff becomes geq 0
  nll_diff = -1
  while(nll_diff<0){
    
    # fit HMM again
    # data are now the last locations
    hmm_time[[i]] <- system.time(hmm_results[[i]] <- fit_hmm(dat = list(X = ssm_results[[i-1]]$x[,1], 
                                                                        Y = ssm_results[[i-1]]$y[,1], 
                                                                        top = ssm_results[[i-1]]$top, 
                                                                        timescale=timescalehmm,
                                                                        setinitstatdist=setinitstatdist),
                                                             p_start = p_start_hmm,
                                                             silence = allsilent))
    # the data here are the last locations
    
    # update on where we are (iteratively) in the optimization
    cat("\n ------------------------------------------------------------ \n finished HMM", i,  ", convergence was:",
        hmm_results[[i]]$mess, 
        "\n ------------------------------------------------------------ \n ")
    
    # switch the generator matrix
    # could also switch the D values (if we ended up fixing them)
    ssmargs$p_start$working_A = matrix(hmm_results[[i]]$params[row.names(hmm_results[[i]]$params) %in% 'working_A',1],nrow=m)
    # p_start_ssm$logD_xy <- hmm_results[[i]]$params[row.names(hmm_results[[i]]$params) == "logD_xy", 1] #%>% sort(decreasing=TRUE)
    # overwrite the logD_xy if they're fixed
    if("logD_xy" %in% names(ssmargs$mapping)){
      ssmargs$p_start$logD_xy <- hmm_results[[i]]$params[row.names(hmm_results[[i]]$params) == "logD_xy", 1] 
    }
    # update the data for the ssm, changing the behavioural states and the matrix exponentials
    ssmargs$dat$b <- hmm_results[[i]]$b_hat-1
    ssmargs$dat$matexp <- hmm_results[[i]]$obj$report()$matexp
 
 
    # block for false convergence, same deal as above
    fc_count[i]=0
    ssm_conv[[i]]=0 #set it to 1 to initialize the loop
    runme = TRUE
    while(runme){
      
      if((sum(ssm_conv[[i]])>0) & (allowfc==FALSE) & (fc_count[i]>0)){
        # print a message letting us know we're working on false convergence
        cat("\n --------------------------- \n working on false convergence \n --------------------------- \n")
        
        # jiggle the starting parameters, jiggle_fc is the st dev for random normal error
        for(j in pidx) ssmargs$p_start[[j]] <- ssmargs$p_start[[j]] + rnorm(length(ssmargs$p_start[[j]]), 0, jiggle_fc)
      }
      
      # fit SSM
      ssm_time[[i]] <- system.time(ssm_results[[i]] <- do.call(fit_ssm, ssmargs))

      # save the new convergence, cue on this
      ssm_conv[i] <- ifelse(ssmargs$optimizer=="nloptr", 
                            ssm_results[[i]]$opt$status, 
                            ssm_results[[i]]$opt$convergence)
      # keep track of nll if we use fixsteps=FALSE to tell the optimizer to stop right after the nll starts to increase
      nll_ssm[i] = ifelse(ssmargs$optimizer=="nloptr", 
                          ssm_results[[i]]$opt$objective, 
                          ssm_results[[i]]$nll)
      
      fc_count[i] = fc_count[i] + 1 
      if(fc_count[i]>fcmax) warning(paste("had false convergence", fcmax, "times while trying to get the SSM to fit"))

      if(sum(ssm_conv[[i]])==0 | allowfc==TRUE | fc_count[i]>=fcmax) runme=FALSE
      
    } # close while loop
    
    
    # save the next nll value
    nll_ssm[i] = ssm_results[[i]]$nll
    # calculate the difference between this and the last nll
    nll_diff <- nll_ssm[i]-nll_ssm[i-1]
    # if we've fixed the number of iterations then set this difference back to a negative number
    if(fixsteps) nll_diff <- -1
    # if we've reached the max number of steps allowed then also set this back to a negative value
    if(length(ssm_results)==maxsteps) nll_diff <- 1
    
    # update on where we are in the iteration
    cat("\n ------------------------------------------------------------ \n finished SSM", i,", convergence was:",
        ssm_results[[i]]$mess, 
        "\n ------------------------------------------------------------ \n ")    
    # add to the iteration counter
    i=i+1
  } # close while loop
  
  
  #####################################
  ############ get results ############
  #####################################
  
  # nll values, convergence and messages, stuff that wasn't tracked in the for loops
  nll_hmm <- sapply(hmm_results, function(x)x$nll)
  hmm_conv <- sapply(hmm_results, function(x)x$opt$convergence)
  hmm_mess <- sapply(hmm_results, function(x)x$opt$message)
  ssm_mess <- sapply(ssm_results, function(x)x$opt$message)
  # put the results into a more convenient data frame
  ssm_nll <- data.frame(iter=1:length(ssm_results), mess = ssm_mess, conv = ssm_conv, nll = nll_ssm)
  hmm_nll <- data.frame(iter=1:length(ssm_results), mess = hmm_mess, conv = hmm_conv, nll = nll_hmm)
  
  # the best of the models based on proper convergence and minimum ssm nll
  winner = ssm_nll[ssm_nll$conv==0,'iter'][which.min(ssm_nll[ssm_nll$conv==0,'nll'])]
  if(length(winner)==0){
    winner=1
    warning("no winning iteration reached; setting winner to 1")
  }
  
  # bring up the best states
  lstates = ssm_results[[winner]]$l_hat
  bstates = hmm_results[[winner]]$b_hat
 
  
  # final results list
  rslts <- list(args = match.call(),
                hmm_results = hmm_results, ssm_results = ssm_results, 
                hmm_time = hmm_time, ssm_time = ssm_time,
                hmm_nll = hmm_nll, ssm_nll = ssm_nll, 
                fc_count = fc_count, winner = winner,
                lstates = lstates, 
                bstates = bstates)
  class(rslts) <- "issm" #set a class, super simple, just for printing/plotting/summary
  return(rslts)
  
}




############################################
############ set a print method ############
############################################

# self explanatory
print.issm <- function(mod){
  
  cat("----------------------------------------------------------- \n the best iteration is ... \n\n")
  print(mod$winner)
  cat("\n\n\n")
  
  cat("----------------------------------------------------------- \n the associated HMM parameters are ... \n\n")
  print(round(mod$hmm_results[[mod$winner]]$params,4))
  cat("\n\n\n")
  
  cat("----------------------------------------------------------- \n the associated SSM parameters are ... \n\n")
  print(round(mod$ssm_results[[mod$winner]]$params,4))
  cat("\n\n\n")
  
  cat("----------------------------------------------------------- \n the convergence information for the SSM is ... \n\n")
  print(mod$ssm_nll)
  cat("\n\n\n")
  

  
}





############################################
############ set a plot method #############
############################################

# this function gets the slope and intercept of a qqline for the residuals below
get_qqline <- function(vec){
  y <- quantile(vec[!is.na(vec)], c(0.25, 0.75))
  x <- qnorm(c(0.25, 0.75))
  m <- diff(y)/diff(x)
  b <- y[1] - m * x[1]
  return(c(m,b))
}


# plot the predicted states along with the pseudoresiduals
plot.issm <- function(mod, sim=NULL){
  
  fmf <- wesanderson::wes_palette("FantasticFox1", type = "discrete")
  
  if(!is.null(sim)){
    lplot <- ggplot() + 
      ylab("Northing") +
      xlab("Easting") +
      geom_path(data=as.data.frame(sim$locs), aes(x=coords.x1, y=coords.x2), col='navy', lwd=1.1) +
      geom_path(data=as.data.frame(mod$ssm_results[[mod$winner]]$l_hat), aes(x=x, y=y), col='cyan3', lwd=1.1) +
      ggtitle("True and Predicted Path") +
      theme_bw()
  }
  
  # plot the predicted locations and their predicted behavioural states
  bplot <- ggplot() + 
    ylab("Northing") +
    xlab("Easting") +
    geom_path(data=as.data.frame(mod$ssm_results[[mod$winner]]$l_hat), aes(x=x, y=y)) +
    geom_point(data=as.data.frame(mod$ssm_results[[mod$winner]]$l_hat), aes(x=x, y=y, col=as.factor(mod$hmm_results[[mod$winner]]$b_hat))) +
    scale_color_manual(values=c(fmf[2], fmf[3], fmf[5])) +
    guides(colour=guide_legend(title="b state")) +
    ggtitle("Predicted Path and Predicted States") +
    theme_bw()

  # take the residuals, gotta get rid of the rows with Inf
  pseudo <- as.data.frame(mod$hmm_results[[mod$winner]]$pseudo)
  names(pseudo) <- c("lon", "lat")
  pseudo <- pseudo[is.finite(rowSums(pseudo)),]

  # plot the residuals
  linetmp <- get_qqline(pseudo[,1]) # get qqline parameters
  xresidplot <- ggplot(data=pseudo, aes(sample=lon)) +
    geom_qq() +
    geom_abline(slope=linetmp[1], intercept=linetmp[2], col=fmf[4]) + 
    ggtitle("X Coordinate Residuals") +
    theme_bw()

  linetmp <- get_qqline(pseudo[,2])
  yresidplot <- ggplot(data=pseudo, aes(sample=lat)) +
    geom_qq() +
    geom_abline(slope=linetmp[1], intercept=linetmp[2], col=fmf[4]) +
    ggtitle("Y Coordinate Residuals") +
    theme_bw()

  # put it all together
  if(is.null(sim)){
    grid.arrange(bplot,
                 xresidplot, yresidplot,
                 widths=c(2,1),
                 layout_matrix=rbind(c(1,2),
                                     c(1,3)),
                 top = paste("winning iteration =", mod$winner))
  } else {
    grid.arrange(lplot, bplot,
                 xresidplot, yresidplot,
                 widths=c(2,1),
                 layout_matrix=rbind(c(1,3),
                                     c(2,4)),
                 top = paste("winning iteration =", mod$winner))
  }

}

# retrieve just the qqplots
plotqq <- function(mod){
  
  fmf <- wesanderson::wes_palette("FantasticFox1", type = "discrete")
  
  pseudo <- as.data.frame(mod$hmm_results[[mod$winner]]$pseudo)
  names(pseudo) <- c("lon", "lat")
  pseudo <- pseudo[is.finite(rowSums(pseudo)),]
  
  linetmp <- get_qqline(pseudo[,1]) # get qqline parameters
  xresidplot <- ggplot(data=pseudo, aes(sample=lon)) +
    geom_qq() +
    geom_abline(slope=linetmp[1], intercept=linetmp[2], col=fmf[4]) + 
    # ggtitle("X Coordinate Residuals") +
    theme_bw() + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) 
  
  linetmp <- get_qqline(pseudo[,2])
  yresidplot <- ggplot(data=pseudo, aes(sample=lat)) +
    geom_qq() +
    geom_abline(slope=linetmp[1], intercept=linetmp[2], col=fmf[4]) +
    # ggtitle("Y Coordinate Residuals") +
    theme_bw() + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) 
  
  rslt <- list(xqqplot=xresidplot, yqqplot=yresidplot)
  return(rslt)
}




bootstrap <- function(mod, 
                      useSSMds = TRUE,
                     startseed, 
                     nsims, 
                     savepath=NULL,
                     dstart=60, 
                     simargs = list(),
                     issm.args = list()
){
  
  
  hmmparnames <- row.names((mod$hmm_results[[mod$winner]]$params))
  ssmparnames <- row.names((mod$ssm_results[[mod$winner]]$params))
  m = length(which(hmmparnames=="logD_xy"))
  if(useSSMds){
    simpars = list(alpha = matrix(mod$hmm_results[[mod$winner]]$params[hmmparnames %in% "A",1],
                                  nrow=m, byrow=FALSE),
                   sigmabi = exp(mod$ssm_results[[mod$winner]]$params['logSigma_bi',1]),
                   Dxy = cumsum(exp(mod$ssm_results[[mod$winner]]$params[ssmparnames %in% 'logD_xy',1])),
                   scale = exp(mod$ssm_results[[mod$winner]]$params['logScale',1]))
  } else {
    simpars = list(alpha = matrix(mod$hmm_results[[mod$winner]]$params[hmmparnames %in% "A",1],
                                  nrow=m, byrow=FALSE),
                   sigmabi = exp(mod$ssm_results[[mod$winner]]$params['logSigma_bi',1]),
                   Dxy = mod$hmm_results[[mod$winner]]$params[hmmparnames %in% 'D_xy',1],
                   scale = exp(mod$ssm_results[[mod$winner]]$params['logScale',1]))
  }
  simargs$par = simpars 
  
 
  sims <- list()
  mods <- list()
  
  for(i in 1:nsims){
    
    seed = startseed-1+i
    
    simargs$startseed = seed
    
    sims[[i]] <- do.call(gentrack, simargs)
    
    inp_init <- getmodinp(sims[[i]]$hydros %>% rename(hx=x, hy=y) %>% data.table(), 
                           sims[[i]]$toa, sdInits=1, 
                           rbi_min=10, rbi_max=30, biTable=sims[[i]]$top$bivec, 
                           m=1, b = NULL,
                          dstart=dstart)
    
    inp_ssm <- getmodinp(sims[[i]]$hydros %>% rename(hx=x, hy=y) %>% data.table(), 
                          sims[[i]]$toa, sdInits=1, 
                          rbi_min=10, rbi_max=30, biTable=sims[[i]]$top$bivec, 
                          m=m, b = NULL, setre0=TRUE, 
                         dstart=dstart)
    
    
    issmargs <- append(list(inp_init = inp_init, 
                            inp_ssm = inp_ssm),
                       issm.args)
  
    mods[[i]] <- tryCatch({
      do.call(fit_issm, issmargs)
    },
    error=function(e)e)
    if(!is.null(savepath)) saveRDS(mods, file=paste(savepath, '_bootstrapmods.RDS', sep=''))
    
    # update on where we are in the iteration
    cat("\n ------------------------------------------------------------ \n finished sim study", i,
        "\n ------------------------------------------------------------ \n ")
  }
  
    
    nullidx <- which(sapply(mods, length) < 10)  
    if(length(nullidx)>0) for(i in 1:length(nullidx)) mods[[nullidx[i]+1-i]] <- NULL
    if(length(nullidx)>0) for(i in 1:length(nullidx)) sims[[nullidx[i]+1-i]] <- NULL
    
    # get trues
    true.pars <- c(simpars$Dxy, as.numeric(simpars$alpha), simpars$Dxy, simpars$sigmabi, simpars$scale)
    
    # get estimated
    pars <- data.frame(par = c(rep("hmm_D_xy", m), rep("A", m*m), rep("ssm_D_xy", m), 'sigmabi', 'scale'),
                       true = true.pars) %>% 
      cbind(., rbind(sapply(mods, function(x)x$hmm_results[[x$winner]]$params[hmmparnames %in% 'D_xy', 'Estimate']),
                  sapply(mods, function(x)x$hmm_results[[x$winner]]$params[hmmparnames %in% 'A', 'Estimate']),
                  apply(exp(sapply(mods, function(x)x$ssm_results[[x$winner]]$params[ssmparnames %in% 'logD_xy', 'Estimate'])), 2, cumsum),
                  exp(sapply(mods, function(x)x$ssm_results[[x$winner]]$params['logSigma_bi', 'Estimate'])),
                  exp(sapply(mods, function(x)x$ssm_results[[x$winner]]$params['logScale', 'Estimate']))
    ) %>% as.data.frame() )
    
    # parameters
    par.stats <- pars %>% 
      dplyr::select(-true, -par) %>% 
      mutate(mean = rowMeans(.),
             median = apply(., 1, median), 
             sd = apply(., 1, sd),
             lower2.5 = apply(., 1, quantile, probs=0.025),
             upper97.5 = apply(., 1, quantile, probs=0.975),
             meanbias = (rowMeans(.)-true.pars),
             rmse = sqrt(rowMeans((. - matrix(rep(true.pars, length(mods)), ncol=length(mods)))^2)),
             par=c(rep("hmm_D_xy", m), rep("A", m*m), rep("ssm_D_xy", m), 'sigmabi', 'scale'),
             true= true.pars) %>% 
      select(par, true, mean, median, sd, lower2.5, upper97.5, meanbias, rmse)
    
    # behavioural states
    # adjust for 0-1 loss with three states
    # err.rate = colMeans((sapply(mods, function(x)x$bstates) - sapply(sims, function(x)x$states))^2) #only accurate for two states
    statediff <- ifelse((sapply(mods, function(x)x$bstates) - sapply(sims, function(x)x$states))==0, 0, 1)
    err.rate  <- colMeans(statediff)
    
    #location states
    lon.rmse = sqrt(colMeans((sapply(mods, function(x)x$lstates[,1]) - sapply(sims, function(x)x$locs[,1]))^2)) 
    lat.rmse = sqrt(colMeans((sapply(mods, function(x)x$lstates[,2]) - sapply(sims, function(x)x$locs[,2]))^2)) 
    rmsemat = mapply(FUN=function(x,y) sqrt(rowSums((x-y)^2)),
                            lapply(mods, function(x)x$lstates),
                            lapply(sims, function(x)x$locs),
                            SIMPLIFY=TRUE)
    rmse = colMeans(rmsemat)
    
    # convergence 
    winnermess <- sapply(mods, function(x)x$ssm_results[[x$winner]]$mess)
    
    # results
    boot <- list(args = match.call(),
                 sims=sims, 
                 mods=mods,
                 true.pars=true.pars,
                 pars=pars, 
                 par.stats=par.stats, 
                 err.rate=err.rate, 
                 lon.rmse=lon.rmse,
                 lat.rmse=lat.rmse, 
                 rmse=rmse,
                 winnermess=winnermess)
    class(boot) <- 'yams_bootstrap'
    if(!is.null(savepath)) saveRDS(boot, file=paste(savepath, '.RDS', sep=''))
    
  
    return(boot)
    
    
  }






