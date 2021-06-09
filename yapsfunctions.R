#' Get prepared inp-object for use in TMB-call
#'
#' Wrapper-function to compile a list of input needed to run TMB
#' @param hydros Dataframe from simHydros() or Dataframe with columns hx and hy containing positions of the receivers. Translate the coordinates to get the grid centre close to (0;0).
#' @param toa TOA-matrix: matrix with receivers in rows and detections in columns. Make sure that the receivers are in the same order as in hydros, and that the matrix is very regular: one ping per column (inlude empty columns if a ping is not detected).
#' @param E_dist Which distribution to use in the model - "Gaus" = Gaussian, "Mixture" = mixture of Gaussian and t or "t" = pure t-distribution
#' @param n_ss Number of soundspeed estimates: one estimate per hour is usually enough
#' @param pingType Type of transmitter to simulate - either stable burst interval ('sbi'), random burst interval ('rbi') or random burst interval but where the random sequence is known a priori
#' @param rbi_min,rbi_max Minimum and maximum BI for random burst interval transmitters
#' @param sdInits If >0 initial values will be randomized around the normally fixed value using rnorm(length(inits), mean=inits, sd=sdInits)
#' @param ss_data_what What speed of sound (ss) data to be used. Default ss_data_what='est': ss is estimated by the model. Alternatively, if ss_data_what='data': ss_data must be provided and length(ss_data) == ncol(toa)
#' @param ss_data Vector of ss-data to be used if ss_data_what = 'est'. Otherwise ss_data <- 0 (default)
#' @param biTable Table of known burst intervals. Only used when pingType == "pbi". Default=NULL


#' @return List of input data ready for use in TMB-call
#' @export
getmodinp <- function(hydros, toa, sdInits=1, rbi_min=0, rbi_max=0, biTable=NULL, 
                      m, b=NULL, setre0=FALSE, timescale=60, dstart=1, setinitstatdist=1, working_A = matrix(0)){
  # m is the number of states
  # b is a state vector (if applicable)
	inp_params <- gettranslatevars(hydros, toa)
	datTmb <- getdat(hydros, toa, rbi_min, rbi_max, biTable, inp_params, m=m, b=b)
	datTmb$timescale <- timescale
	datTmb$setinitstatdist <- setinitstatdist
	datTmb$matexp <- matrix(0)
	params <- getpars(datTmb, m, scale=dstart, set0=setre0)
	params$working_A = working_A
	inits <- getpinits(sdInits, m)
	return(list(
		datTmb = datTmb,
		params= params,
		inits = inits,
		inp_params = inp_params
		)
	)
}



#' Get inits for use in TMB
#'
#' Compile a vector of initial values to use in TMB. One value for each estimated parameter (not random effects).
#' Should all be in a credible range.
#' @inheritParams getInp
#' @return Vector of initial values to use in TMB
#' @export
getpinits <- function(sdInits=1, m) {
  # set.seed(42)
	init_logD_xy <- rep(1, m)
	init_logSigma_bi <- -5
	init_logScale <- 1
	inits <- c(init_logD_xy, init_logSigma_bi,  init_logScale)

	inits <- stats::rnorm(length(inits), mean=inits, sd=sdInits)
	return(inits)
}



#' Get data for input to TMB
#'
#' Compile data for input to TMB.
#' @param inp_params Selection of parameters used to setup and run YAPS.
#' @inheritParams getInp
#'
#' @return List for use in TMB.
#' @export
getdat <- function(hydros, toa, rbi_min, rbi_max, biTable, inp_params, m, b){
	T0 <- inp_params$T0
	Hx0 <- inp_params$Hx0
	Hy0 <- inp_params$Hy0
	
	toa <- toa - T0
	
	# allowing slight out-of-bounds BIs
	rbi_min <- rbi_min - rbi_min * 0.05
	rbi_max <- rbi_max + rbi_max * 0.05
	
	# attempting to make sure toa is oriented correct
	if(!nrow(toa) == nrow(hydros)){
		toa <- t(toa)
	}
	
	# when m =1 this will all be 1s
	# when m>1 and no b's are set then the code will randomly sample
	# if b != NULL then it just keeps b the way it is
	if(is.null(b)) b <- sample(1:m,  size = ncol(toa), replace=TRUE)

	dat <- list(
		H = matrix(c(hydros$hx-Hx0, hydros$hy-Hy0), ncol=2),
		toa = toa,
		nh = nrow(hydros),
		np = ncol(toa),
		rbi_min = rbi_min,
		rbi_max = rbi_max,
		biTable = biTable,
		b = b - 1 # subtract 1 for compatibility with C++ indexing
	)

	return(dat)
}

#' Get params-list for use in TMB
#'
#' Compile a list of parameters for use in TMB.
#' @param datTmb Object obtained using getDatTmb()
#' @return List of params for use in TMB
#' @export
getpars <- function(datTmb, m, scale, set0 = FALSE){
	params_XY <- getxyinitsfromcoa(datTmb)
	if(set0){
	  X <- rep(0, ncol(datTmb$toa))
	  Y <- rep(0, ncol(datTmb$toa))
	  tag_drift <- rep(0, datTmb$np)
	} else {
	  X = params_XY$X + stats::rnorm(ncol(datTmb$toa), sd=10)
	  Y = params_XY$Y + stats::rnorm(ncol(datTmb$toa), sd=10)
	  tag_drift = stats::rnorm(datTmb$np, 0, 1e-2)
	}
	
	list(X = X,
	     Y = Y,
	     logD_xy = rep(log(scale), m),				#diffusivity of transmitter movement (D_xy in ms)
	     logSigma_bi = 0,			#sigma  burst interval (sigma_bi in ms)
	     logScale = 0,				#scale parameter for t-distribution
	     tag_drift = tag_drift
	)
}

#' Get initial values for X and Y based on Center Of Activity - i.e. hydrophones positions
#'
#' Attempts to give meaningful initial values for X and Y based on which hydros detected each ping
#' @inheritParams getInp
#' @noRd
getxyinitsfromcoa <- function(datTmb){
	toa <- datTmb$toa
	hydros <- datTmb$H

	toa_detect <- toa
	toa_detect[!is.na(toa_detect)] <- 1

	X <- zoo::na.approx(colMeans((toa_detect) * hydros[,1], na.rm=TRUE))
	Y <- zoo::na.approx(colMeans((toa_detect) * hydros[,2], na.rm=TRUE))
	
	return(list(X=X, Y=Y))

}

#' Get parameters for this specific data set
#'
#' Compile a list of relevant parameters (e.g. T0) to use later on
#' @inheritParams getInp
gettranslatevars <- function(hydros, toa){
	T0 <- min(toa, na.rm=TRUE)
		
	Hx0 <- hydros[1,hx]
	Hy0 <- hydros[1,hy]

	return(list(T0=T0, Hx0=Hx0, Hy0=Hy0))

}















getNobs <- function(toa){
	nobs <- apply(toa, 1, function(k) sum(!is.na(k)))
	return(nobs)

}
plotNobs <- function(toa){
	nobs <- getNobs(toa)
	plot(nobs)
	lines(caTools::runmean(nobs, k=50), col="red")
	lines(caTools::runmean(nobs, k=100), col="blue")
}

plotToa <- function(toa){
	best_h <- which.max(apply(toa, 2, function(k) sum(!is.na(k))))
	# finds the column (hydrophone) with the most number of detections
	# uses that as a reference column
	suppressWarnings(matplot(toa[,best_h], (toa[,best_h] - toa) * 1450))
	# plots the best_h toa column of toa against the time differences 
	# with all other columns multiplied by the speed of sound
	# so basically plotting the best hydrophone timestamp of each detection
	# vs the distances from the other receivers 
	# i'm guessing this doesn't do well with NAs which there will be
	# because sometimes the best hydrophone won't detect the animal
	# whereas other receivers will, hence the suppressWarnings()
	# 
}


downSampleToa <- function(toa, max_hydro){
	toa_down <- toa
	nobs <- apply(toa, 1, function(k) sum(!is.na(k)))
	plot(nobs)
	lines(runmean(nobs, k=100), col="red")

	while(sum(nobs > max_hydro) > 0){
		toa_down[nobs > max_hydro, ] <- toa_down[nobs >  max_hydro] * rbinom(length(toa_down[nobs >  max_hydro]), 1, .95)
		toa_down[toa_down == 0] <- NA
		nobs <- apply(toa_down, 1, function(k) sum(!is.na(k)))
		plot(nobs)
		lines(runmean(nobs, k=100), col="red")
	}
	nobs <- apply(toa_down, 1, function(k) sum(!is.na(k)))
	plot(nobs)
	lines(runmean(nobs, k=100), col="red")

	return(toa_down)
}


getSS <- function(temp, sal=0, depth=5){
	ss <- 1410 + 4.21 * temp - 0.037 * temp^2 + 1.10 * sal + 0.018 * depth
	return(ss)
}

buildToaKnownSeq <- function(seq, aligned_dat, hydros){
	toa <- matrix(nrow=length(seq), ncol=nrow(hydros))
	inp <- as.matrix(aligned_dat[,c('seq_ping_idx', 'hydro_idx', 'eposync')])
	toa[inp[, 1:2]] <- inp[,3]
	
	nobs <- apply(toa, 1, function(k) sum(!is.na(k)))
	
	first_det <- head(which(nobs >= 1), n=1)
	last_det <- tail(which(nobs >= 1), n=1)
	toa <- toa[first_det:last_det, ]
	seq <- 	seq[first_det:last_det]

	return(list(toa=toa, seq=seq))
}
