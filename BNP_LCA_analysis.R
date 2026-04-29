# latent class analysis with vaccine survey data
# wave 1 (2024)
library(poLCA)
library(dplyr)
library(nnet)
library(ggplot2)
library(effects)
library(tidyr)
library(purrr)

set.seed(123)

load(here::here("../data_clean_wave1.RData"))

# model formula (all variables must be factors or integers)
# no covariates
form0 <- cbind(vacc_last_yr, vacc_next_yr, mmr_doses, vacc_eff, concerned) ~ 1

data0 <- data[, c(
  "vacc_last_yr", "vacc_next_yr", "mmr_doses", "vacc_eff", "concerned",
  "Age", "Gender", "social_grade", "ethnicity_grp", "marital_status", "marital_collapsed", "house_tenure",
  "working_status", "parent_guardian18", "parent_guardian", "children_in_household")]

# check NAs
all(sapply(data0, FUN = \(x) table(is.na(x))) == nrow(data0))

data_complete <- na.omit(data[, c(
  "vacc_last_yr", "vacc_next_yr", "mmr_doses", "vacc_eff", "concerned",
  "Age", "Gender", "social_grade", "ethnicity_grp", "marital_status", "marital_collapsed", "house_tenure",
  "working_status", "parent_guardian18", "parent_guardian", "children_in_household")])

data_null <- data_complete

library(forcats)
data_null$vacc_next_yr <- fct_collapse(data_null$vacc_next_yr,
                                       "Likely" = c("Very likely", "Fairly likely"),
                                       "Unlikely" = c("Not very likely", "Not at all likely"),
                                       "Unknown" = "Don't know/ unsure")

data_null$mmr_doses <- fct_collapse(data_null$mmr_doses,
                                    "Yes" = "Yes, I have",
                                    "No" = "No, I haven't",
                                    "Don't know/ unsure" = "Don't know/ unsure",
                                    "Not applicable" = c(
                                      "Not parent",
                                      "Not applicable - My oldest child is still not old enough to get the MMR vaccine"
                                    )
)


data_null$vacc_eff <- fct_collapse(data_null$vacc_eff,
                                   "Agree"              = c("Strongly agree", "Tend to agree"),
                                   "Disagree"           = c("Strongly disagree", "Tend to disagree"),
                                   "Neutral" = c(
                                     "Neither agree nor disagree",
                                     "Don't know"
                                   )
)

data_null$concerned <- fct_collapse(data_null$concerned,
                                    "Agree"              = c("Strongly agree", "Tend to agree"),
                                    "Disagree"           = c("Strongly disagree", "Tend to disagree"),
                                    "Neutral" = c(
                                      "Neither agree nor disagree",
                                      "Don't know"
                                    )
)

indicator_data <- data_null[, c("vacc_last_yr", "vacc_next_yr",
                                "mmr_doses", "vacc_eff", "concerned")]
Y <- indicator_data %>%
  mutate(across(everything(), ~ as.integer(as.factor(.))))
Y = as.matrix(Y)
# 1. Basic Dimensions
n <- nrow(Y)        # Number of samples
p <- ncol(Y)        # Number of items
d_j <- c(3, 3, 4, 3, 3) # Number of categories per item
alpha <- 0.05#0.3460165
# 2. Initialize z (Each sample in its own component)
#z <- 1:n
k = 2
z <- sample(1:k, n, replace = TRUE)
# 3. Initialize PSI
# Since d_j varies, we use a list of matrices
PSI <- vector("list", p)

for(j in 1:p) {
  # Create the matrix
  PSI[[j]] <- matrix(NA, nrow = k, ncol = d_j[j])
  # IMPORTANT: Initialize rownames to match z (1 to n)
  rownames(PSI[[j]]) <- as.character(1:k)
  # Sample Dirichlet for each cluster
  for(i in 1:k) {
    PSI[[j]][i, ] <- extraDistr::rdirichlet(1, rep(1, d_j[j]))
  }
}
# Assume beta_jc is 1 for all j, c (flat prior)
beta <- lapply(d_j, function(d) rep(1, d))

calculate_class_likelihood <- function(i, h, Y, x) {
  p <- ncol(Y)

  # Initialize likelihood
  likelihood <- 1

  for (j in 1:p) {
    # Get the observed category for individual i, item j
    observed_category <- Y[i, j]

    # Select the probability of that category for class h
    # This represents: beta_{jhv}^1 where v = Y[i,j]
    prob_val <- x[[j]][h, observed_category]

    # Multiply into the product
    likelihood <- likelihood * prob_val
  }

  return(likelihood)
}

# Function to calculate the marginal likelihood for a NEW cluster
# Y: Data matrix (n x p)
# i: Index of the individual
# beta_prior: List of length p, where beta_prior[[j]] is a vector of length d_j
calculate_new_cluster_marginal <- function(i, Y, beta_prior) {
  p <- ncol(Y)

  # Initialize likelihood
  marginal_likelihood <- 1

  for (j in 1:p) {
    # Calculate sum of betas for item j
    sum_beta_j <- gamma(sum(beta_prior[[j]]))
    prod_beta_j <- prod(gamma(beta_prior[[j]]))
    beta_sum_gamma_j = gamma(sum(beta_prior[[j]]) + 1)


    # Get the observed category for individual i, item j
    observed_category <- Y[i, j]

    # Extract beta for the observed category
    beta_jc <- beta_prior[[j]][observed_category]
    beta_prior_temp = beta_prior[[j]]
    beta_prior_temp[observed_category] = beta_prior_temp[observed_category] + 1
    prod_beta_temp_gamma = prod(gamma(beta_prior_temp))

    # Multiply into the product across all items p
    marginal_likelihood <- marginal_likelihood *
      (  (prod_beta_temp_gamma * sum_beta_j) / (beta_sum_gamma_j * prod_beta_j) )
  }

  return(marginal_likelihood)
}


####################
n_iter <- 100000  # Total iterations

# --- Storage for results ---
history_z <- matrix(NA, nrow = n_iter, ncol = n)
history_K <- numeric(n_iter)

history_log_post <- numeric(n_iter) # <--- NEW: Storage for log-posterior
history_PSI <- vector("list", n_iter) # <--- NEW: List to store PSI objects
for (iter in 1:n_iter) {

  # 1. Update Cluster Assignments
  for (i in 1:n) {
    z_others <- z[-i]
    active_labels <- unique(z_others)
    k_current <- length(active_labels)

    n_h_table <- table(factor(z_others, levels = active_labels))
    n_h_minus_i <- as.numeric(n_h_table)

    all_h_likelihoods <- sapply(active_labels, function(h) {
      calculate_class_likelihood(i, h, Y, PSI)
    })

    z_old_weights <- (n_h_minus_i / (n + alpha - 1)) * all_h_likelihoods
    z_new_weight  <- (alpha / (n + alpha - 1)) * calculate_new_cluster_marginal(i, Y, beta)

    total_w <- sum(z_old_weights) + z_new_weight
    probs <- c(z_old_weights, z_new_weight) / total_w

    sample_idx <- sample(1:(k_current + 1), 1, prob = probs)

    if (sample_idx <= k_current) {
      z[i] <- active_labels[sample_idx]
    } else {
      new_label <- max(z) + 1
      z[i] <- new_label
      for (quest in 1:p) {
        ones <- rep(0, length(beta[[quest]]))
        ones[Y[i, quest]] <- 1
        new_param <- extraDistr::rdirichlet(1, beta[[quest]] + ones)
        PSI[[quest]] <- rbind(PSI[[quest]], new_param)
        rownames(PSI[[quest]])[nrow(PSI[[quest]])] <- as.character(new_label)
      }
    }
  }

  # 2. Cleanup & Relabeling
  active_final <- unique(z)
  for (quest in 1:p) {
    PSI[[quest]] <- PSI[[quest]][rownames(PSI[[quest]]) %in% as.character(active_final), , drop = FALSE]
  }

  new_map <- match(z, active_final)
  for (quest in 1:p) {
    PSI[[quest]] <- PSI[[quest]][match(as.character(active_final), rownames(PSI[[quest]])), , drop = FALSE]
    rownames(PSI[[quest]]) <- as.character(1:length(active_final))
  }
  z <- new_map

  # 3. Update PSI Parameters (Post-Assignment Update)
  K_current <- length(unique(z))
  for (h in 1:K_current) {
    indices_in_h <- which(z == h)
    for (quest in 1:p) {
      counts_h_quest <- table(factor(Y[indices_in_h, quest], levels = 1:d_j[quest]))
      c_hj <- as.numeric(counts_h_quest)
      PSI[[quest]][h, ] <- extraDistr::rdirichlet(1, beta[[quest]] + c_hj)
    }
  }

  # 4. Save progress
  history_z[iter, ] <- z
  history_K[iter] <- K_current

  # A. Log-likelihood of the data: log P(Y | z, PSI)
  log_lik <- 0
  for (quest in 1:p) {
    # For each observation, find its cluster 'h' and its response value
    # Then look up the probability in PSI[[quest]][h, response]
    cluster_assignments <- z
    responses <- Y[, quest]

    # Extract probabilities for each person based on their cluster
    probs <- PSI[[quest]][cbind(cluster_assignments, responses)]
    log_lik <- log_lik + sum(log(probs))
  }

  history_log_post[iter] <- log_lik

  history_PSI[[iter]] <- PSI

  # Optional: Print progress every 10 iterations
  if (iter %% 10 == 0) cat("Iteration:", iter, "| Active Clusters:", K_current, "\n")
}


#after burnin
history_K = history_K[3000:10000]
history_PSI = history_PSI[3000:10000]
history_log_post = history_log_post[3000:10000]

table(history_K)
K_mode = 6
# 1. Find the iterations where K was 6
idc_mode <- which(history_K == K_mode)

# 2. Find which of *those* iterations had the highest log-posterior
# we use idc_mode[which.max(...)] to get the original iteration number
best_iter_idx <- idc_mode[which.max(history_log_post[idc_mode])]

# 3. Extract the PSI object from that specific iteration
map_psi <- history_PSI[[best_iter_idx]]

library(gtools)
# 1. Setup
idc_6 <- which(history_K == 6)
ref_psi <- map_psi  # Using  MAP estimate
all_perms <- permutations(n = 6, r = 6) # Matrix of 720 possible orderings

aligned_history_PSI_6 <- list()
for (i in seq_along(idc_6)) {
  iter_idx <- idc_6[i]
  current_psi_list <- history_PSI[[iter_idx]]

  best_dist <- Inf
  best_perm <- 1:6

  # 3. Find the best permutation for this iteration
  for (p_idx in 1:nrow(all_perms)) {
    perm <- all_perms[p_idx, ]
    current_dist <- 0

    # Calculate distance  across all questions
    for (quest in 1:length(current_psi_list)) {
      # Permute the rows of the current PSI matrix
      permuted_mat <- current_psi_list[[quest]][perm, ]
      current_dist <- current_dist + sqrt(sum((permuted_mat - ref_psi[[quest]])^2))
    }

    if (current_dist < best_dist) {
      best_dist <- current_dist
      best_perm <- perm
    }
  }

  # 4. Save the aligned version
  aligned_psi <- lapply(current_psi_list, function(mat) mat[best_perm, ])
  aligned_history_PSI_6[[i]] <- aligned_psi
}


# 1. Setup
n_people <- nrow(Y)
n_sims <- length(aligned_history_PSI_6)
simulated_z <- matrix(NA, nrow = n_people, ncol = n_sims)

# 2. Loop through each iteration's aligned PSI
for (s in 1:n_sims) {
  current_psi_list <- aligned_history_PSI_6[[s]]

  # Calculate likelihood for each person in each of the 6 clusters
  log_lik_matrix <- matrix(0, nrow = n_people, ncol = 6)

  for (h in 1:6) {
    for (quest in 1:ncol(Y)) {
      # Probabilities for this specific cluster's profile
      probs_h <- current_psi_list[[quest]][h, Y[, quest]]
      log_lik_matrix[, h] <- log_lik_matrix[, h] + log(probs_h)
    }
  }

  # Convert log-likelihoods to probabilities (Softmax)
  # We assume a uniform prior here as K is fixed at 6
  probs_matrix <- t(apply(log_lik_matrix, 1, function(x) {
    exp_x <- exp(x - max(x))
    exp_x / sum(exp_x)
  }))

  # Sample a cluster for each person based on these probabilities
  simulated_z[, s] <- apply(probs_matrix, 1, function(p) {
    sample(1:6, size = 1, prob = p)
  })
}



# 1. Calculate the log-likelihood for each person being in each of the 6 clusters
n <- nrow(Y)
K <- 6
p <- ncol(Y)
log_post_probs <- matrix(0, nrow = n, ncol = K)

for (h in 1:K) {
  for (quest in 1:p) {
    # Get probabilities for all people for this question if they were in cluster h
    probs <- map_psi[[quest]][h, Y[, quest]]
    log_post_probs[, h] <- log_post_probs[, h] + log(probs)
  }
}

# 2. Convert log-probabilities to normalized posterior probabilities
# We subtract the max log-prob for numerical stability
adj_log_probs <- t(apply(log_post_probs, 1, function(x) x - max(x)))
post_probs <- exp(adj_log_probs) / rowSums(exp(adj_log_probs))
# Replace 0s with a tiny number to avoid log(0) = -Inf
post_probs_safe <- ifelse(post_probs == 0, 1e-10, post_probs)
# Total Entropy
entropy_val <- -sum(post_probs_safe * log(post_probs_safe))
# Normalized Entropy (0 to 1 scale)
# 1 is perfect separation, 0 is total overlap
relative_entropy <- 1 - (entropy_val / (n * log(K)))

############################
class_probs = matrix(0,dim(simulated_z)[1],K)
for(i in 1:dim(simulated_z)[1]){
  class_probs[i,c(as.numeric(names(table(simulated_z[i,]))))] =
    table(simulated_z[i,])/sum(table(simulated_z[i,]))
}

row_modes = apply(class_probs,1,which.is.max)
miscl_probs = matrix(0,K,K)
for(i in 1:K){
  for(j in 1:K){
    miscl_probs[i,j] = sum(class_probs[which(row_modes == i),j])/sum(class_probs[,j])
  }
}

miscl_probs = array(NA,c(K,K,dim(simulated_z)[2]))
for(l in 1:K){
  for(j in 1:K){
    for(iter in 1:dim(simulated_z)[2]){
      miscl_probs[l,j,iter] = sum(simulated_z[which(row_modes == l),iter] == j)/
        sum(simulated_z[,iter] == j)
    }
  }
}


###########################################
df_bayesian_cov <- data_complete[,c(6,7,8,9,11,12,13,14,16)]
# 1. Select covariates
X_data <- df_bayesian_cov

# 2. Create the design matrix (including intercept)
# This automatically handles the dummy coding for categorical factors
X <- model.matrix(~ ., data = X_data)

library(DescTools)
row_modes = c()
for(i in 1:dim(simulated_z)[1]){
  row_modes[i] = Mode(simulated_z[i,])
}
# Apply to each row (1 indicates rows)
# row_modes <- apply(simulated_z, 1, get_mode)

# 3. Prepare the data list for Stan
stan_data <- list(
  N = nrow(X),              # Number of observations
  K = 6,                    # Number of clusters
  D = ncol(X),              # Number of predictors (including intercept)
  X = X,                    # The design matrix
  y = row_modes # Cluster labels (1 to 6)
)

library(rstan)

# Compile and sample
fit <- stan(
  file = "two_step_mult.stan", # Or file = "model.stan"
  data = stan_data,
  chains = 1,
  iter = 1500,
  warmup = 1000,
  cores = parallel::detectCores()
)
#save(fit, file = "fit_mult_zero_constr_MAP.RData")
# Extract posterior samples of beta
beta_post <- rstan::extract(fit, "beta")$beta

odds_ratio_beta = array(NA,c(dim(beta_post)[1],dim(beta_post)[2],dim(beta_post)[3]))
for(i in 1:dim(beta_post)[3]){
  for(j in 1:dim(beta_post)[2]){
    odds_ratio_beta[,j,i] = exp(beta_post[,j,i])
  }
}

mean_log_or <- apply(log(odds_ratio_beta), c(2,3), mean)
lower_log_or <- apply(log(odds_ratio_beta), c(2,3), quantile, probs = 0.025)
upper_log_or <- apply(log(odds_ratio_beta), c(2,3), quantile, probs = 0.975)
cov_names <- colnames(X)

# Convert matrix to long dataframe
df_summary <- as.data.frame(mean_or) %>%
  mutate(covariate = cov_names) %>%
  pivot_longer(
    cols = -covariate,
    names_to = "cluster",
    values_to = "mean_odds"
  )
df_summary$cluster <- paste("Cluster", as.numeric(gsub("V", "", df_summary$cluster)))
# Plot
ggplot(df_summary, aes(x = mean_odds, y = factor(covariate, levels = rev(cov_names)))) +
  geom_vline(size = 1, xintercept = 0, linetype = "dashed", color = "red") +
  geom_point(size = 2.5, color = "blue") +
  facet_wrap(~ cluster, ncol = 6) +
  labs(
    x = "Posterior mean log odds ratio (cluster vs rest)",
    y = "Covariate",
    title = "Posterior mean log odds ratio for each cluster and covariate"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.x = unit(0.5, "cm")
  )

library(rstan)
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Extract and Reconstruct the Beta Matrix
# ---------------------------------------------------------
post_beta <- extract(fit)$beta
n_iter <- dim(post_beta_raw)[1]
D      <- ncol(X)
K      <- 6

post_beta[,,1]
# Reconstruct full beta [D, K] using posterior means
# We include the first column as zeros (your reference Cluster 1)
beta_means <- matrix(0, nrow = D, ncol = K)
beta_means[, 2:K] <- apply(post_beta_raw, c(2, 3), mean)
colnames(beta_means) <- paste0("Cluster_", 1:K)
rownames(beta_means) <- colnames(X)

# 2. Function to calculate Average Marginal Effects (AME)
# ---------------------------------------------------------
calc_ame <- function(target_cov_idx, design_matrix, beta_mat) {
  # Baseline probabilities for all observations
  xb_base <- design_matrix %*% beta_mat
  # Softmax calculation: exp(xb) / rowSums(exp(xb))
  probs_base <- exp(xb_base) / rowSums(exp(xb_base))

  # Probabilities if we increase the covariate by 1 unit
  X_alt <- design_matrix
  X_alt[, target_cov_idx] <- X_alt[, target_cov_idx] + 1
  xb_alt <- X_alt %*% beta_mat
  probs_alt <- exp(xb_alt) / rowSums(exp(xb_alt))

  # AME is the average change in probability across the sample
  return(colMeans(probs_alt - probs_base))
}

# 3. Calculate AMEs for all predictors
# ---------------------------------------------------------
# Map the function across all columns of X
ame_list <- lapply(1:D, function(d) calc_ame(d, X, beta_means))
ame_matrix <- do.call(rbind, ame_list)
rownames(ame_matrix) <- colnames(X)
colnames(ame_matrix) <- paste0("Cluster_", 1:K)

library(dplyr)
library(tidyr)
library(ggplot2)

# --- 1. Prepare the AME Data ---
ame_long <- as.data.frame(ame_matrix) %>%
  mutate(Covariate = rownames(.)) %>%
  filter(Covariate != "(Intercept)") %>%
  pivot_longer(cols = starts_with("Cluster"),
               names_to = "Cluster",
               values_to = "AME")

# --- 2. Generate Individual Plots ---
# We loop through each unique cluster name
cluster_names <- unique(ame_long$Cluster)

for (cl in cluster_names) {

  # Filter data for just this cluster
  plot_data <- ame_long %>%
    filter(Cluster == cl) %>%
    # Sort by absolute impact for the plot
    mutate(abs_ame = abs(AME)) %>%
    arrange(desc(abs_ame))

  # Create the plot
  p <- ggplot(plot_data, aes(x = reorder(Covariate, AME), y = AME, fill = AME > 0)) +
    geom_bar(stat = "identity", width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = c("firebrick3", "dodgerblue4"),
                      labels = c("Lower Prob.", "Higher Prob.")) +
    theme_minimal(base_size = 12) +
    labs(
      title = paste("Covariate Contributions:", cl),
      subtitle = "Average Marginal Effect (Absolute change in probability)",
      x = NULL,
      y = "Change in Probability",
      fill = "Direction"
    ) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )

  # Display the plot in your R console/Plots pane
  print(p)

  # Optional: Save each plot as a PNG
  # ggsave(filename = paste0("AME_", cl, ".png"), plot = p, width = 8, height = 6)
}



post = extract(fit)

# 1. Calculate Mean and 95% Credible Intervals
beta_samples <- post$beta[,,3] # Dimensions: [iterations, predictors]

means <- apply(beta_samples, 2, mean)
lower <- apply(beta_samples, 2, quantile, probs = 0.025)
upper <- apply(beta_samples, 2, quantile, probs = 0.975)

# 2. Create a data frame for plotting
plot_data <- data.frame(
  predictor = 1:ncol(beta_samples),
  mean = means,
  lower = lower,
  upper = upper
)

# 3. Filter for "Significant" predictors (CI does not cross 0)
sig_data <- plot_data[plot_data$lower > 0 | plot_data$upper < 0, ]

# Sort by mean for better visualization
sig_data <- sig_data[order(sig_data$mean), ]
sig_data$predictor_factor <- factor(1:nrow(sig_data), labels = sig_data$predictor)

colnames(X)[sig_data$predictor]
library(ggplot2)

ggplot(sig_data, aes(x = mean, y = predictor_factor)) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  geom_pointrange(aes(xmin = lower, xmax = upper),
                  color = "royalblue", size = 0.3, alpha = 0.7) +
  labs(title = "Posterior Coefficients (Cluster 1)",
       subtitle = "Only predictors with 95% CI not overlapping zero are shown",
       x = "Effect Size (Mean & 95% Credible Interval)",
       y = "Predictor Index") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 6)) # Shrink text if many variables remain





##################################

# Parameters
S <- dim(simulated_z)[2]  # 5000 draws
N <- dim(simulated_z)[1]  # 2053 individuals
K <- 6                  # 6 Clusters

# Initialize the weight matrix (Individuals as rows, Clusters as columns)
w_nk <- matrix(0, nrow = N, ncol = K)

# Loop through each individual to count their cluster frequency
for (i in 1:N) {
  # Create a factor with fixed levels 1:K to ensure 0s are counted for empty clusters
  draws_for_i <- factor(simulated_z[i, ], levels = 1:K)

  # Tabulate and divide by total draws (S)
  w_nk[i, ] <- as.numeric(table(draws_for_i)) / S
}

# Verification: Every row must sum to 1
if(all(abs(rowSums(w_nk) - 1) < 1e-9)) {
  message("Success: All individual weights sum to 1.")
}



df_bayesian_cov <- data_complete[,6:16]
# 1. Select your covariates
X_data <- df_bayesian_cov

# 2. Create the design matrix (including intercept)
# This automatically handles the dummy coding for your categorical factors
X <- model.matrix(~ ., data = X_data)
# 3. Prepare the data list for Stan
stan_data <- list(
  N = nrow(X),              # Number of observations
  K = 6,                    # Number of clusters
  D = ncol(X),              # Number of predictors (including intercept)
  X = X,                    # The design matrix
  w = w_nk # Cluster labels (1 to 6)
)

library(rstan)

# Compile and sample
fit <- stan(
  file = "soft_two_step.stan", # Or file = "model.stan"
  data = stan_data,
  chains = 2,
  iter = 1000,
  cores = parallel::detectCores()
)



Y <- indicator_data %>%
  mutate(across(everything(), ~ as.integer(as.factor(.))))
Y = as.matrix(Y)
N <- nrow(Y)
J <- ncol(Y)
R <- as.vector(apply(Y,2,max))
library(rstan)

stan_data <- list(
  N = N,
  J = J,
  K = 3,
  R = R,
  Y1 = Y[,1],
  Y2 = Y[,2],
  Y3 = Y[,3],
  Y4 = Y[,4],
  Y5 = Y[,5]
)

fits <- list()

for(k in 1:10){

  stan_data$K <- k

  fits[[k]] <- stan(
    file = "lca.stan",
    data = stan_data,
    chains = 2,
    iter = 1000,
    warmup = 500,
    seed = 123,
    control = list(adapt_delta = 0.9),
    cores = parallel::detectCores()
  )

}


library(rstan)
library(gtools)
library(matrixStats)

# --- 1. SETUP ---
K_val <- 3
fit_3 <- fits[[K_val]]
params_3 <- rstan::extract(fit_3)

n_samples <- dim(params_3$pi)[1]
n_people <- nrow(Y)

# --- 2. FIND MAP REFERENCE ---
# Using the log-posterior to find the reference iteration
map_idx <- which.max(params_3$lp__)
ref_theta <- list(
  params_3$theta1[map_idx, , ],
  params_3$theta2[map_idx, , ],
  params_3$theta3[map_idx, , ],
  params_3$theta4[map_idx, , ],
  params_3$theta5[map_idx, , ]
)

# --- 3. INITIALIZE STORAGE ---
all_perms <- gtools::permutations(n = K_val, r = K_val)
simulated_z_3 <- matrix(NA, nrow = n_people, ncol = n_samples)

# List to store the aligned theta matrices for each iteration
aligned_thetas_history_3 <- vector("list", n_samples)

# --- 4. ALIGNMENT & ALLOCATION LOOP ---
pb <- txtProgressBar(min = 0, max = n_samples, style = 3)

for (s in 1:n_samples) {

  # Current sample theta values (unaligned)
  current_thetas <- list(
    params_3$theta1[s, , ],
    params_3$theta2[s, , ],
    params_3$theta3[s, , ],
    params_3$theta4[s, , ],
    params_3$theta5[s, , ]
  )

  best_dist <- Inf
  best_perm <- 1:K_val

  # Find the permutation that minimizes distance to the MAP reference
  for (p in 1:nrow(all_perms)) {
    perm <- all_perms[p, ]
    d <- 0
    for (j in 1:5) {
      d <- d + sum((current_thetas[[j]][perm, ] - ref_theta[[j]])^2)
    }

    if (d < best_dist) {
      best_dist <- d
      best_perm <- perm
    }
  }

  # --- 5. SAVE ALIGNED PARAMETERS ---
  # Reorder rows of each theta matrix based on the best permutation
  aligned_thetas_history_3[[s]] <- lapply(current_thetas, function(mat) mat[best_perm, ])

  # --- 6. ALLOCATION SAMPLING ---
  # Align mixing proportions
  pi_aligned <- params_3$pi[s, best_perm]
  log_probs <- matrix(0, nrow = n_people, ncol = K_val)

  for (k in 1:K_val) {
    log_probs[, k] <- log(pi_aligned[k]) +
      log(params_3$theta1[s, best_perm[k], Y[,1]]) +
      log(params_3$theta2[s, best_perm[k], Y[,2]]) +
      log(params_3$theta3[s, best_perm[k], Y[,3]]) +
      log(params_3$theta4[s, best_perm[k], Y[,4]]) +
      log(params_3$theta5[s, best_perm[k], Y[,5]])
  }

  # Draw a cluster assignment for each individual
  simulated_z_3[, s] <- apply(log_probs, 1, function(lp) {
    p <- exp(lp - matrixStats::logSumExp(lp))
    sample(1:K_val, size = 1, prob = p)
  })

  setTxtProgressBar(pb, s)
}
close(pb)


# --- 1. SETUP ---
K_val <- 6
fit_6 <- fits[[K_val]]
params_6 <- rstan::extract(fit_6)

n_samples <- dim(params_6$pi)[1]
n_people <- nrow(Y)

# --- 2. FIND MAP REFERENCE ---
map_idx <- which.max(params_6$lp__)
ref_theta <- list(
  params_6$theta1[map_idx, , ],
  params_6$theta2[map_idx, , ],
  params_6$theta3[map_idx, , ],
  params_6$theta4[map_idx, , ],
  params_6$theta5[map_idx, , ]
)

# --- 3. INITIALIZE STORAGE ---
all_perms <- gtools::permutations(n = K_val, r = K_val)
simulated_z_6 <- matrix(NA, nrow = n_people, ncol = n_samples)

# List to store the aligned theta matrices for each iteration
# This will be a list of lists: [sample][item_number]
aligned_thetas_history <- vector("list", n_samples)

# --- 4. ALIGNMENT & ALLOCATION LOOP ---
pb <- txtProgressBar(min = 0, max = n_samples, style = 3)

for (s in 1:n_samples) {

  # Current sample theta values (unaligned)
  current_thetas <- list(
    params_6$theta1[s, , ],
    params_6$theta2[s, , ],
    params_6$theta3[s, , ],
    params_6$theta4[s, , ],
    params_6$theta5[s, , ]
  )

  best_dist <- Inf
  best_perm <- 1:K_val

  # Find the permutation that minimizes distance to the MAP reference
  for (p in 1:nrow(all_perms)) {
    perm <- all_perms[p, ]
    d <- 0
    for (j in 1:5) {
      d <- d + sum((current_thetas[[j]][perm, ] - ref_theta[[j]])^2)
    }

    if (d < best_dist) {
      best_dist <- d
      best_perm <- perm
    }
  }

  # --- 5. SAVE ALIGNED PARAMETERS ---
  # Reorder rows of each theta matrix based on the best permutation
  aligned_thetas_history[[s]] <- lapply(current_thetas, function(mat) mat[best_perm, ])

  # --- 6. ALLOCATION SAMPLING ---
  pi_aligned <- params_6$pi[s, best_perm]
  log_probs <- matrix(0, nrow = n_people, ncol = K_val)

  for (k in 1:K_val) {
    log_probs[, k] <- log(pi_aligned[k]) +
      log(params_6$theta1[s, best_perm[k], Y[,1]]) +
      log(params_6$theta2[s, best_perm[k], Y[,2]]) +
      log(params_6$theta3[s, best_perm[k], Y[,3]]) +
      log(params_6$theta4[s, best_perm[k], Y[,4]]) +
      log(params_6$theta5[s, best_perm[k], Y[,5]])
  }

  simulated_z_6[, s] <- apply(log_probs, 1, function(lp) {
    p <- exp(lp - matrixStats::logSumExp(lp))
    sample(1:K_val, size = 1, prob = p)
  })

  setTxtProgressBar(pb, s)
}
close(pb)




# 1. Define a helper function to find the mode (MAP) for each individual
get_map_allocation <- function(allocations_matrix) {
  apply(allocations_matrix, 1, function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  })
}

# 2. Calculate the MAP allocation for each individual
map_k3 <- get_map_allocation(simulated_z_3)
map_k6 <- get_map_allocation(simulated_z_6)

# 3. Generate the Cross-Table (Confusion Matrix)
map_cross_table <- table(MAP_K6 = map_k6, MAP_K3 = map_k3)

# 4. View the table
print(map_cross_table)

# 5. Convert to proportions
map_prop_table <- prop.table(map_cross_table, margin = 1)
print(round(map_prop_table, 3))


