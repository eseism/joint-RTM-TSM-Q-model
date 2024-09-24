data {
  int<lower=0> n_R;
  int<lower=0> n_T;
  int<lower=0> n_S;
  int<lower=0> p;
  int<lower=0> r;
  vector[n_R] Y_R;
  vector[n_T] Y_T;
  matrix[n_R, p] X_R;
  matrix[n_T, p] X_T;
  vector[n_S] Y_S;
  matrix[n_T, r] Phi_T;
  matrix[n_S, r] Phi_S;
}
parameters {
  vector<lower=0.001, upper=0.1>[p] beta; 
  real<lower=0> sigma_R;
  real<lower=0> sigma_T;
  real<lower=0> sigma_S;
  real<lower=0> sigma_beta;
  vector<lower=-5, upper=5>[r] eta;
}
model {
  beta ~ normal(0.005, sqrt(sigma_beta));
  sum(eta) ~ normal(0,0.001*r);
  sigma_beta ~ inv_gamma(1, 1);
  sigma_R ~ inv_gamma(1, 1);
  sigma_T ~ inv_gamma(1, 1);
  sigma_S ~ inv_gamma(1, 1);
  Y_R ~ normal(X_R * beta, sqrt(sigma_R));
  Y_T ~ normal(X_T * beta + Phi_T*eta, sqrt(sigma_T));
  Y_S ~ normal(Phi_S*eta, sqrt(sigma_S));
}

