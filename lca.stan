data {
  int<lower=1> N;                 // 2053
  int<lower=1> K;                 // number of classes
  int<lower=1, upper=3> Y1[N];    // item 1 (3 categories)
  int<lower=1, upper=3> Y2[N];    // item 2 (3 categories)
  int<lower=1, upper=4> Y3[N];    // item 3 (4 categories)
  int<lower=1, upper=3> Y4[N];    // item 4 (3 categories)
  int<lower=1, upper=3> Y5[N];    // item 5 (3 categories)
}

parameters {
  simplex[K] pi;                 // class proportions

  simplex[3] theta1[K];
  simplex[3] theta2[K];
  simplex[4] theta3[K];
  simplex[3] theta4[K];
  simplex[3] theta5[K];
}

model {

  // ---- Priors ----
  pi ~ dirichlet(rep_vector(1.0, K));

  for (k in 1:K) {
    theta1[k] ~ dirichlet(rep_vector(1.0, 3));
    theta2[k] ~ dirichlet(rep_vector(1.0, 3));
    theta3[k] ~ dirichlet(rep_vector(1.0, 4));
    theta4[k] ~ dirichlet(rep_vector(1.0, 3));
    theta5[k] ~ dirichlet(rep_vector(1.0, 3));
  }

  // ---- Likelihood (marginalized over classes) ----
  for (i in 1:N) {

    vector[K] log_prob;

    for (k in 1:K) {
      log_prob[k] =
          log(pi[k])
        + log(theta1[k][Y1[i]])
        + log(theta2[k][Y2[i]])
        + log(theta3[k][Y3[i]])
        + log(theta4[k][Y4[i]])
        + log(theta5[k][Y5[i]]);
    }

    target += log_sum_exp(log_prob);
  }
}

generated quantities {

  vector[N] log_lik;

  for (i in 1:N) {

    vector[K] log_prob;

    for (k in 1:K) {
      log_prob[k] =
          log(pi[k])
        + log(theta1[k][Y1[i]])
        + log(theta2[k][Y2[i]])
        + log(theta3[k][Y3[i]])
        + log(theta4[k][Y4[i]])
        + log(theta5[k][Y5[i]]);
    }

    log_lik[i] = log_sum_exp(log_prob);
  }
}
