data {
  int<lower=1> N;             // Number of observations
  int<lower=1> K;             // Number of clusters
  int<lower=1> D;             // Number of predictors
  matrix[N, D] X;             // Design matrix
  int<lower=1, upper=K> y[N]; // Response labels
}

parameters {
  matrix[D, K-1] beta_raw;    // Coefficients for clusters 2 through 6
}

transformed parameters {
  matrix[D, K] beta;
  // Enforce sum-to-zero across clusters for each predictor
  for (d in 1:D) {
    beta[d] = append_row(beta_raw[d]', -sum(beta_raw[d]))';
  }
}

model {
  // Prior: Weakly informative Normal priors
  to_vector(beta_raw) ~ normal(0, 2);

  y ~ categorical_logit_glm(X, rep_vector(0, K), beta);
}
