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
	DATA_VECTOR(X); // animal positions, eastings
  DATA_VECTOR(Y); // animal positions, northings
  DATA_VECTOR(top); // time of ping
  DATA_SCALAR(timescale); //how much to scale the time by
  DATA_SCALAR(setinitstatdist); 
  int n_obs = X.size();
  
	
	/////////////////////////////
	//////// Paramaeters ////////
	/////////////////////////////

  // movement of animal
	PARAMETER_VECTOR(logD_xy);    		// Diffusivity of fish
	int m = logD_xy.size();
	vector<Type> D_xy_int = exp(logD_xy);
	vector<Type> D_xy = D_xy_int;
	for(int i = 1; i < m; ++i){
	  D_xy(i) = D_xy(i-1) + D_xy_int(i);
	}
  REPORT(D_xy);
	ADREPORT(D_xy);
  
  // the working generator matrix
	PARAMETER_MATRIX(working_A); //matrix of switching probabilities

  // actual generator matrix
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
	ADREPORT(A);

  ////////////////////////////////
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
  // set delta_row based on option: 1 for the stat dist, 0 for uniform 
  matrix<Type> delta_row(1,m);
  if(setinitstatdist==1){
    for(int i=0; i<m; i++) delta_row(0,i) = delta(i); 
  } else {
    for(int i=0; i<m; i++) delta_row(0,i) = 1.0/m; 
  }


	REPORT(delta); // stationary dist vec
	REPORT(delta_row); // this one is used in the likelihood


	
	



	//////////////////////////////////////////
	///////// Movement Probabilities /////////
	//////////////////////////////////////////



	//// probability density matrices
	vector<Type> TMP(2); //temporary variable

	array<Type> P_array(m, m, n_obs);
	array<Type> C_x_array(m, m, n_obs);
	array<Type> C_y_array(m, m, n_obs);
	for(int k=0; k<n_obs; k++) {    // k indexes time
	  for(int i=0; i<m; i++) {
	    for(int j=0; j<m; j++) { // i and j both index behavioural state
	      if( i==j ) { // diagonal entries
	        //process equation
	        if(k == 0){
	          P_array(i,j,k) = 0.0;
	            // dnorm(X(0),Type(0),Type(10000),true);
	            // dnorm(Y(0),Type(0),Type(10000),true);
	            C_x_array(i,j,k) = 0.0;
	            C_y_array(i,j,k) = 0.0;
	        } else {
              Type diffop = (top(k) - top(k-1))/timescale;
	            P_array(i,j,k) = exp(dnorm(X(k), X(k-1),sqrt(2*D_xy(i)*diffop), true))*
	              exp(dnorm(Y(k), Y(k-1),sqrt(2*D_xy(i)*diffop), true));
	            C_x_array(i,j,k) = pnorm(X(k), X(k-1),sqrt(2*D_xy(i)*diffop));
	            C_y_array(i,j,k) = pnorm(Y(k), Y(k-1),sqrt(2*D_xy(i)*diffop));

	        }
	      } else {
	        P_array(i,j,k) = 0.0; //off-diagonals are zero
	        C_x_array(i,j,k) = 0.0;
	        C_y_array(i,j,k) = 0.0;
	      }
	    }
	  }
	} // matrix form for ease of forward algorithm

	REPORT(P_array);


	
	
	

	/////////////////////////////////////
	///////// Viterbi Algorithm /////////
	/////////////////////////////////////



	array<Type> viterbi(m, n_obs);    // Columns index time, rows index state
	array<Type> forward_max(m);    // Vector of most likely states
	vector<int> states(n_obs);
  array<Type> matexp(m, m, n_obs);

	int min = 0;
	int max = n_obs;

	//calculate the matrix exponentials once
	for(int k=1; k<n_obs; k++){
	  matrix<Type> matxprod = A*(top(k)-top(k-1))/timescale;
	  matrix<Type> matxexp = expm(matxprod);
	  for(int i = 0; i<m; i++){
	    for(int j = 0; j<m; j++){
	      matexp(i,j,k) = matxexp(i,j);
	    }
	  }
	}
	REPORT(matexp);

//starting state likelihoods
for(int i=0; i<m; i++) {
  viterbi(i,min) = log(delta_row(0, i));
}

// forward probabilities
for(int k=(min+1); k<max; k++) {
  for(int i=0; i<m; i++) {
    for(int j=0; j<m; j++) {
        forward_max(j) = viterbi(j,k-1) + log(matexp(i,j,k)); // transition probabilities
    }
    viterbi(i,k) = forward_max.maxCoeff() + log(P_array(i,i,k));    // Choose most likely transition
  }
}

REPORT(viterbi);

Eigen::Index max_row; 
Eigen::Index max_col;
Type foo;

// Get the last column of viterbi matrix for backwards probabilities
matrix<Type> column = viterbi.col(max-1);  
foo = column.maxCoeff(&max_row, &max_col);    

states(max-1) = max_row + 1;    

// Backwards recursion
for(int k=n_obs-2; k>=0; k--) {
  for(int i=0; i<m; i++) { // index state
      column(i,0) = viterbi(i,k) +  log(matexp(i, max_row,k));    // backwards Viterbi coefficient
  }
  foo = column.maxCoeff(&max_row, &max_col);    // sends index of largest coefficient to max_row
  // figures out which of the previous states led to the next one
  states(k) = max_row + 1;    // Eigen::Index starts at 0, want it to start at 1
}

REPORT(states);




	//////////////////////////////
	///////// Likelihood /////////
	//////////////////////////////



	Type ll = 0.0;

	// Forward Probabilities as row vector
	matrix<Type> alpha(1,m);
	// let's store them all 
	matrix<Type> pseudo(n_obs, 2);

	// calculate likelihood and pseudoresiduals
	for(int k=0; k<n_obs; k++) {
	  if( k == 0 ) {
	    pseudo(k) = delta_row.sum();
	    alpha = delta_row;
	  } else {
	        pseudo(k,0) = ( alpha * matexp.col(k).matrix() * ( C_x_array.col(k).matrix())).sum() / alpha.sum();   // pseudoresidual for the kth longitude
	        pseudo(k,1) = ( alpha * matexp.col(k).matrix() * ( C_y_array.col(k).matrix())).sum() / alpha.sum();   // pseudoresidual for the kth latitude
	      alpha = alpha * matexp.col(k).matrix() * ( P_array.col(k).matrix() ); // Add k'th observation to forward probability
	  }
	  ll += log(alpha.sum());  // add log of forward probability to recursive likelihood
	  alpha = alpha/alpha.sum();    // rescale forward probabilities to sum to 1
	}


	REPORT(pseudo);
	return -ll; 


}
