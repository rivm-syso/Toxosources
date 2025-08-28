functions {
  real si_model( real lambda, real age ){
    real p;
    p = 1-exp(-lambda * age );
    return p;
  }
  
  real sis_model( real lambda, real gamma, real age){
    real p;
    p = - lambda * expm1(-(gamma+lambda)*age) / (lambda+gamma);
    return p;
  }
  
  // Sum to zero constraints: 
  // https://discourse.mc-stan.org/t/test-soft-vs-hard-sum-to-zero-constrain-choosing-the-right-prior-for-soft-constrain/3884/32
  vector Q_sum_to_zero_QR(int N) {
    vector [2*N] Q_r;
  
    for(i in 1:N) {
      Q_r[i] = -sqrt((N-i)/(N-i+1.0));
      Q_r[i+N] = inv_sqrt((N-i) * (N-i+1));
    }
    return Q_r;
  }
  
  vector sum_to_zero_QR(vector x_raw, vector Q_r) {
    int N = num_elements(x_raw) + 1;
    vector [N] x;
    real x_aux = 0;
  
    for(i in 1:N-1){
      x[i] = x_aux + x_raw[i] * Q_r[i];
      x_aux = x_aux + x_raw[i] * Q_r[i+N];
    }
    x[N] = x_aux;
    return x;
  }
}

data {
  // https://github.com/stan-dev/stan/issues/2377
  int<lower=0, upper=1> use_sis_model;

  int<lower=1>       n;
  int<lower=1>       n_region;
  int<lower=1>       n_tot[n];
  int<lower=0>       n_pos[n];
  
  int<lower=1, upper=n_region>      region[n];

  vector<lower=0>[n] age_low;
  vector<lower=0>[n] age_high;
  vector<lower=0>[n] age;
  
}

transformed data{
  real sigma_region_raw = inv_sqrt(1 - inv(n_region));
  vector[2*n_region] Q_region = Q_sum_to_zero_QR(n_region);
}

parameters {
  real  lambda_baseline;
  vector[n_region-1]      lambda_region_raw;
  real                    gamma_baseline[use_sis_model];
  real<lower=0>           sigma_lambda_region;  
  vector<lower=0, upper=1>[n] age_raw;
  vector<lower=0, upper=1>[n] weight;
}


transformed parameters{
  vector[n_region]  lambda_region = sum_to_zero_QR(lambda_region_raw, Q_region );
}

model {
  vector[n] p;
  real      lambda[n];
  real      gamma[use_sis_model, n];
  vector[n] agedist;
 
  age_raw ~ beta_proportion( (age-age_low)./(age_high-age_low), 1 ); # Arno: hier stond 10

  agedist = age_low + (age_high - age_low) .* age_raw;
  
  // Centered parametrisation: sigma * N(0,1) instead of N(0,sigma)
  // Complicated expression for the SD below ensures sum-to-one.
  // https://mc-stan.org/docs/2_27/stan-users-guide/parameterizing-centered-vectors.html
  lambda_baseline         ~ normal( 0, 3 );
  
  lambda_region_raw       ~ normal(0, sigma_region_raw);
  sigma_lambda_region     ~ std_normal(  );

  if( use_sis_model ){
    gamma_baseline[use_sis_model] ~ normal( 0, 3);
  
    for( j in 1:n ){
      gamma[use_sis_model,j] = exp( gamma_baseline[use_sis_model]);
    }
  }
  
  for( j in 1:n )
    lambda[j] = exp(lambda_baseline + sigma_lambda_region *lambda_region[region[j]] );
                     
  if( use_sis_model ){    
    for( j in 1:n )
      p[j] = sis_model( lambda[j], 
                        gamma[1,j], //first index has size one
                        agedist[j] ); 
  }else{
    for( j in 1:n )
      p[j] = si_model( lambda[j], agedist[j] );
  }
  
  // Version with weighting of the observations, because some populations
  // are tested multiple times they should be downweighted.
  for( i in 1:n ){
    target += weight[i] * binomial_lpmf(n_pos[i]| n_tot[i], p[i] );
  }
}

