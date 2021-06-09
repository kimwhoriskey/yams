#include <TMB.hpp>
using namespace density;

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
Type softplus(Type x,Type epsilon)
{
  return 0.5*(x+sqrt(x*x+epsilon*epsilon));
}

template<class Type>
bool isFinite(Type x){
  return R_finite(asDouble(x));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  // data 
	DATA_ARRAY(H);			// Position of hydros
	DATA_ARRAY(toa);   		// Time of arrival at hydro. One row per buoy, one column per ping
	DATA_INTEGER(nh);     // number of hydros
	DATA_INTEGER(np);     // number of pings
	DATA_SCALAR(rbi_min); // random burst interval min
	DATA_SCALAR(rbi_max); // random burst interval max
	DATA_VECTOR(biTable); // known burst interval sequence 
	DATA_IVECTOR(b);     // behaviour
  DATA_SCALAR(timescale); // scaling the time differences
  DATA_INTEGER(setinitstatdist); // set the initial distribution for the markov chain
  DATA_ARRAY(matexp); // matrix exponentials from the hidden markov model
  

  
  
  // random effects
	PARAMETER_VECTOR(X);	//Position at time of ping
	PARAMETER_VECTOR(Y);	//Position at time of ping
	int n_obs = X.size();

	
	
	// parameters 
	
	PARAMETER_VECTOR(logD_xy);    		// Diffusivity of fish
	int m = logD_xy.size();
	vector<Type> D_xy_int = exp(logD_xy);
	vector<Type> D_xy = D_xy_int;
	for(int i = 1; i < m; ++i){
	  D_xy(i) = D_xy(i-1) + D_xy_int(i); // ensures they are increasing in order
	}
	REPORT(D_xy);

	PARAMETER(logSigma_bi);		// standard deviation for the tag drift
	Type sigma_bi = exp(logSigma_bi);

	PARAMETER(logScale);		// scale-parameter for t-dist for the time of arrival
	Type scale = exp(logScale);

	
	// random effects
	PARAMETER_VECTOR(tag_drift);

	

	
	
	
  // parameters - the generator matrix
  // if m = 1, doesn't enter likelihood at all
	PARAMETER_MATRIX(working_A); //working generator matrix, mxm
	matrix<Type> A(m, m);
	for(int i = 0; i < m; ++i){
	  for(int j = 0; j < m; ++j){
	    if(i != j){
	      A(i,j) = exp(working_A(i,j));  //ensures the off-diagonal rates are >0
	    } else if (i==j){
	      A(i,j) = 0.0;
	    }
	  }
	  A(i,i) = -A.row(i).sum();
	}
	REPORT(A);
	// solving for the stationary distribution
	// system of equations
	matrix<Type> system(m, m);
	for(int i=0; i < (m-1); i++) {
	  for(int j=0; j<m; j++) {
	    system(i,j) = A(j,i);
	  }
	}
	for(int j = 0; j<m; j++){
	  system(m-1,j) = 1;
	}
	// right hand side of the equation
	vector<Type> rhs(m);
	for(int i=0; i<(m-1); i++) {
	  rhs(i) = 0.0;
	}
	rhs(m-1) = 1.0;
	// take the inverse of the system matrix and then multiply by the vector of ones
	matrix<Type> sys_inverse = system.inverse();
	vector<Type> delta = sys_inverse*rhs;
	REPORT(delta);

	// set delta_row based on option setinitstatdist: 1 for the stat dist, 0 for uniform
	matrix<Type> delta_row(1,m);
	if(setinitstatdist==1){
	  for(int i=0; i<m; i++) delta_row(0,i) = delta(i);
	} else {
	  for(int i=0; i<m; i++) delta_row(0,i) = 1.0/m;
	}

	
	
	
	
	
	
	///////////////////////////////////
	///// likelihood calculations /////
	///////////////////////////////////
	

	array<Type> mu_toa(nh,np);  // expected time of arrival
	array<Type> dist(nh,np);	// distances between the hydrophones and the pings
	vector<Type> top(np); // time of ping
	vector<Type> top_pbi(np); // time of ping without drift

	
	
	// initialize nll
	Type nll = 0.0;

	
	
	/////////////////////////
  // calculate the time of ping from the drift and the programmed time of ping 
	top_pbi(0) = Type(0.0);
	for(int i = 1; i < np; ++i)	{
		top_pbi(i) = top_pbi(i-1) + biTable(i-1);
	}
	for(int i = 0; i < np; ++i)	{
		top(i) = top_pbi(i) + tag_drift(i);
	}
	REPORT(top); // report these random effects

	



	
	/////////////////////////
	// contribution from the behaviour
	// import the matexp as data because otherwise these calculations slow way down because 
	// of the dependence on top, which is a random effect, thus all of those derivatives are calculated
	Type ll_behav = 0.0;

	if(m>1){ 
	  // only evaluate this if we have more than one state
	  // first contribution comes from the stationary distribution
	  ll_behav += log(delta_row(0, b(0)));
	  // contributions from the rest of the behaviours
	  for(int i=1; i < b.size(); ++i){
	    ll_behav += log(matexp(b(i-1), b(i), i));
	    // previous state in rows, the next state in columns
	    // likelihood of this observed markov chain is just the correct entries of the log matrix added together
	  }
	} else {
	  // set to zero if there aren't any behaviours
	  ll_behav = 0.0;
	}
	nll -= ll_behav;

	
	
	
	/////////////////////////
	// contribution from the tag drift
	Type ll_tagdrift = 0.0; 
	ll_tagdrift += dnorm(tag_drift(0), Type(0.0), Type(4), true); // prior
	ll_tagdrift += dnorm(tag_drift(1), Type(0.0), Type(4), true); // prior
	for(int i = 2; i < np; ++i)	{
		ll_tagdrift += dnorm(tag_drift(i)-2*tag_drift(i-1)+tag_drift(i-2), Type(0), sigma_bi, true);
	}
	nll -= ll_tagdrift;
	

	
	

	/////////////////////////
	// contribution from the time of arrivals
  Type ll_toa = 0.0; 
	for(int i=0; i<np; ++i){ //iterate over pings
		for(int h=0; h<nh; ++h){ //iterate over hydros
		  
		  // distance from the hydro to the location, would make it hard to scale the locations
		  // although we could scale the speed of sound too
		  dist(h,i) = sqrt((H(h,0)-X(i))*(H(h,0)-X(i)) + (H(h,1)-Y(i))*(H(h,1)-Y(i)));

			if(!isNA(toa(h,i))){ //ignore NA's

			  // expected time of arrival based on the distance, the real time of ping, and the speed of sound
			  // here speed of sound is constant (1465) but can be replaced by random effects as well
				mu_toa(h,i) = top(i) +  dist(h,i)/1465;
				Type eps = toa(h,i) - mu_toa(h,i);

				// difference in actual time of arrival and expected one follows a scaled t distribution 
				ll_toa += log(dt(eps/scale, Type(3.0), false)/scale);

			}
		}
	}
	REPORT(dist);
	nll -= ll_toa; 

	
	
	
	/////////////////////////
	// contribution from the movement of the animal
	Type ll_loc = 0.0; 
	// starting locations, variance would also have to change if we scaled the hydros positions
	ll_loc += dnorm(X(0),Type(0),Type(10000),true); 
	ll_loc += dnorm(Y(0),Type(0),Type(10000),true);
  vector<Type> diffop(np-1);
	for(int i=1; i<np; ++i)	{
    diffop(i-1) = (top(i) - top(i-1))/timescale;
		ll_loc += dnorm(X(i), X(i-1),sqrt(2*D_xy[b[i]]*diffop(i-1)),true);
		ll_loc += dnorm(Y(i), Y(i-1),sqrt(2*D_xy[b[i]]*diffop(i-1)),true);
	}	
  nll -= ll_loc; 
	

	
	
	
	// a few last things to report
	REPORT(ll_behav);
	REPORT(ll_tagdrift);
	REPORT(ll_toa); 
	REPORT(ll_loc);
	REPORT(nll);
	
	
	return nll;
	
	
}
