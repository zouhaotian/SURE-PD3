data{
  int n; // number of subjects 
  int nobs; // number of observations 
  int P; 
  int items; // number of total items
  int id[nobs]; // subject ID
  matrix[nobs, P] time_spline;
  real time[nobs]; // observed time
  
  int ncat[items]; // number of categories for each item
  int sum_ncat; // sum of categories for each item
  int delta_ncat; // sum_ncat - items
  
  int Y[nobs, items]; // items response
  int obs_cat[nobs, items]; // observed category for Y_{ik}
  vector[2] zero;
}
parameters{
  real a_temp[items]; // difficulty parameter of first category for each item: vectorized form 
  real<lower=0> delta[delta_ncat]; // difference of a[l+1]-a[l], because a must be monotone increasing for one item 
  //for example, a[2] = a[1] + delta[1], delta[1] ~ normal(0, 10) T[0, ];
  real<lower=0> b[items]; // discrimination parameter
  vector[P] beta;
  real<lower=0> sigma_u1;
  real<lower=-1, upper=1> rho_01; // correlation coefficients
  vector[2] U[n]; // random intercept and slope
}
transformed parameters{
  cov_matrix[2] Sigma_u;
  real a[sum_ncat]; // difficulty parameter
  real location[delta_ncat]; // item location parameter
  real theta[nobs]; // theta_i for non-tremor items
  real<lower=0, upper=1> psi_prob[sum_ncat]; // P(Y_{kl}<=l | theta_i(t))
  real<lower=0, upper=1> cat_prob[sum_ncat]; // P(Y_{kl}= l | theta_i(t))
  real<lower=0, upper=1> obs_cat_prob[nobs, items]; // observed categorical probability for Y_{ik}
  
  { // variables within this block cannot be seen outside
    
  int which_cat = 0; // determine which category
  int which_delta = 0;
  int which_location = 0;
  
  for (i in 1:nobs){
    which_cat = 0;
    which_delta = 0;
    which_location = 0;
    theta[i] = time_spline[i]*beta + U[id[i], 1] + U[id[i], 2]*time[i];
    
    for (k in 1:items){
      for (l in 1:(ncat[k]-1)){ // ncat[1]=3
        which_cat+=1;
        which_location+=1;
        if (l==1){
          a[which_cat] = a_temp[k];
          psi_prob[which_cat] = inv_logit(a[which_cat] - b[k]*theta[i]);
          cat_prob[which_cat] = psi_prob[which_cat];
        } else {
          which_delta+=1;
          a[which_cat] = a[which_cat-1] + delta[which_delta];
          psi_prob[which_cat] = inv_logit(a[which_cat] - b[k]*theta[i]);
          cat_prob[which_cat] = psi_prob[which_cat] - psi_prob[which_cat-1];
        }
        location[which_location] = a[which_cat];
      }
      which_cat+=1;
      psi_prob[which_cat] = 1; // last category cumulative probability is 1
      cat_prob[which_cat] = 1 - psi_prob[which_cat-1];
      obs_cat_prob[i, k] = cat_prob[obs_cat[i, k]];    
    }
  }
  
  } // end of this block
  
  Sigma_u[1, 1] = 1;
  Sigma_u[1, 2] = rho_01*sigma_u1;
  Sigma_u[2, 1] = Sigma_u[1, 2];
  Sigma_u[2, 2] = sigma_u1*sigma_u1;
  
}
model{
  for (i in 1:nobs){
    target+=sum(log(obs_cat_prob[i])); // define the log likelihood of categorical distribution
  }
  a_temp ~ normal(0, 10);
  for (i in 1:delta_ncat){
    delta[i] ~ normal(0, 10) T[0, ];
  }
  b ~ uniform(0, 10);
  beta ~ normal(0, 10); 
  sigma_u1 ~ inv_gamma(0.1, 0.1);
  U ~ multi_normal(zero, Sigma_u);
  rho_01 ~ uniform(-1, 1);
}
generated quantities{
  real ll[nobs]; // log likelihood for i-th visit
  for (i in 1:nobs){
    ll[i] = sum(log(obs_cat_prob[i]));
  }
}
