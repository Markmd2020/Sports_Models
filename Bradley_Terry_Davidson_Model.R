# ============================================================
# 1. Load packages
# ============================================================

#Load libararies
library(dplyr)
library(goalmodel)
library(ggplot2)
library(MASS)
library(readr)
library(worldfootballR) 
library(implied)
library(rstan)
library(tidyverse)
library(MCMCpack)



# ============================================================
# 2. Stan model embedded directly inside R
# ============================================================
stan_code <- "
data {
    int<lower=0> N; 
    int<lower=0> P; 
    int<lower=1> team1[N]; 
    int<lower=1> team2[N]; 
    int<lower=1, upper=3> results[N]; 
    real<lower=0> nu_prior_rate;
    vector[P] alpha;
}

parameters {
    simplex[P] ratings;
    real<lower=0> nu;
}

model {
    ratings ~ dirichlet(alpha);
    nu ~ exponential(nu_prior_rate);

    for (i in 1:N){
        real r1 = ratings[team1[i]];
        real r2 = ratings[team2[i]];
        real nu_term = nu * sqrt(r1 * r2);
        real denom = r1 + r2 + nu_term;

        vector[3] p;
        p[1] = r1 / denom;       // home win
        p[3] = nu_term / denom;  // draw
        p[2] = r2 / denom;       // away win

        target += categorical_lpmf(results[i] | p);
    }
}
"

# ============================================================
# 3. Generate synthetic data
# ============================================================
set.seed(135)

P <- 6
N <- 150

alpha <- rep(2, P)
ratings_true <- as.numeric(rdirichlet(1, alpha))
nu_true <- 0.7

team1 <- sample(1:P, N, replace = TRUE)
team2 <- sample(1:P, N, replace = TRUE)

same <- team1 == team2
while(any(same)) {
  team2[same] <- sample(1:P, sum(same), replace = TRUE)
  same <- team1 == team2
}

compute_probs <- function(t1, t2, ratings, nu) {
  r1 <- ratings[t1]
  r2 <- ratings[t2]
  nu_term <- nu * sqrt(r1 * r2)
  denom <- r1 + r2 + nu_term
  c(r1/denom, r2/denom, nu_term/denom)
}

results <- integer(N)
for (i in 1:N) {
  p <- compute_probs(team1[i], team2[i], ratings_true, nu_true)
  results[i] <- sample(1:3, 1, prob = p)
}

stan_data <- list(
  N = N,
  P = P,
  team1 = team1,
  team2 = team2,
  results = results,
  alpha = alpha,
  nu_prior_rate = 1
)

# ============================================================
# 4. Fit the Stan model
# ============================================================
fit <- stan(
  model_code = stan_code,
  data = stan_data,
  iter = 2000,
  chains = 4
)

# ============================================================
# 5. Print posterior summaries
# ============================================================
print(fit, pars = c("ratings", "nu"))

# Posterior means
post <- rstan::extract(fit)

ratings_hat <- colMeans(post$ratings)
nu_hat <- mean(post$nu)

cat("\nTrue ratings:\n")
print(round(ratings_true, 3))

cat("\nEstimated ratings:\n")
print(round(ratings_hat, 3))

cat("\nTrue nu:", nu_true, "\n")
cat("Estimated nu:", round(nu_hat, 3), "\n")

# ============================================================
# 6. Posterior predictive example
# ============================================================
example_probs <- compute_probs(1, 2, ratings_hat, nu_hat)

cat("\nPosterior predictive probabilities for Team 1 vs Team 2:\n")
names(example_probs) <- c("Home win", "Away win", "Draw")
print(round(example_probs, 3))

