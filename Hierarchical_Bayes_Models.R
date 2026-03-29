#######Hierarchical Bayesian Models#######

#Load libraries
library(brms)
library(tidyverse)
library(posterior)
library(bayesplot)
library(tidybayes)
library(loo)

#Generate synthetic data
set.seed(135)

n_teams <- 20
teams   <- paste0("Team_", seq_len(n_teams))

# True latent parameters (unknown in real life)
true_attack  <- rnorm(n_teams, 0, 0.35)
true_defense <- rnorm(n_teams, 0, 0.35)
true_home_adv <- 0.20

# Schedule: random pairings
n_matches <- 600
home_id <- sample(seq_len(n_teams), n_matches, replace = TRUE)
away_id <- sample(seq_len(n_teams), n_matches, replace = TRUE)
# Avoid self-matches
same <- home_id == away_id
while (any(same)) {
  away_id[same] <- sample(seq_len(n_teams), sum(same), replace = TRUE)
  same <- home_id == away_id
}

home_team <- teams[home_id]
away_team <- teams[away_id]

# Poisson rates for goals
log_lambda_home <- true_home_adv + true_attack[home_id] - true_defense[away_id]
log_lambda_away <- true_attack[away_id] - true_defense[home_id]

lambda_home <- exp(log_lambda_home)
lambda_away <- exp(log_lambda_away)

home_goals <- rpois(n_matches, lambda_home)
away_goals <- rpois(n_matches, lambda_away)

matches <- tibble(
  match_id   = seq_len(n_matches),
  home_team  = factor(home_team),
  away_team  = factor(away_team),
  home_goals = home_goals,
  away_goals = away_goals
)

matches %>% glimpse()

#Complete Pooling: No Teams Effects
#This would be high bias
m_pool <- brm(
  home_goals ~ 1,
  data   = matches,
  family = poisson(),
  prior  = c(
    prior(normal(0, 1.0), class = "Intercept")
  ),
  chains = 4, cores = 4, iter = 2000, seed = 135
)
summary(m_pool)

#No Pooling:Separate Teams Parameters which leads to high variance
m_nopool <- brm(
  home_goals ~ 1 + home_team + away_team,
  data   = matches,
  family = poisson(),
  prior  = c(
    prior(normal(0, 1.0), class = "Intercept"),
    prior(normal(0, 0.5), class = "b")  # regularization, but still no pooling
  ),
  chains = 4, cores = 4, iter = 2000, seed = 135
)

summary(m_nopool)

#Partial Pooling:Hierarchical (Multilevel) Team Effects
priors_hier <- c(
  prior(normal(0, 1.0), class = "Intercept"),
  # SD priors control how much team-to-team variation is plausible
  prior(exponential(1.0), class = "sd")
)

m_hier <- brm(
  home_goals ~ 1 + (1 | home_team) + (1 | away_team),
  data   = matches,
  family = poisson(),
  prior  = priors_hier,
  chains = 4, cores = 4, iter = 2500, seed = 135
)

summary(m_hier)

#Model both scores jointly
bf_home <- bf(home_goals ~ 1 + (1 | home_team) + (1 | away_team))
bf_away <- bf(away_goals ~ 1 + (1 | away_team) + (1 | home_team))


m_mv <- brm(
  bf_home + bf_away + set_rescor(FALSE),
  data   = matches,
  family = poisson(),
  prior  = c(
    prior(normal(0, 1), class = "Intercept", resp = "homegoals"),
    prior(normal(0, 1), class = "Intercept", resp = "awaygoals"),
    
    prior(exponential(1), class = "sd", group = "home_team", resp = "homegoals"),
    prior(exponential(1), class = "sd", group = "away_team", resp = "homegoals"),
    prior(exponential(1), class = "sd", group = "home_team", resp = "awaygoals"),
    prior(exponential(1), class = "sd", group = "away_team", resp = "awaygoals")
  ),
  chains = 4, cores = 4, iter = 3000, seed = 135
)

summary(m_mv)

#Diagnostics and Posterior Predictive Checks
# Convergence diagnostics
m_hier %>% summary()

# Posterior predictive checks: do simulated goals resemble observed?
pp_check(m_hier, type = "hist", ndraws = 100)


#Build a negative binomial model to handle overdispersion
m_hier_nb <- brm(
  home_goals ~ 1 + (1 | home_team) + (1 | away_team),
  data   = matches,
  family = negbinomial(),
  prior  = c(
    prior(normal(0, 1.0), class = "Intercept"),
    prior(exponential(1.0), class = "sd"),
    prior(exponential(1.0), class = "shape")  # NB dispersion
  ),
  chains = 4, cores = 4, iter = 2500, seed = 5
)
summary(m_hier_nb)

pp_check(m_hier_nb, type = "hist", ndraws = 100)


#Visualising Partial Pooling
#Shrinkage Effect
# Extract team effects from both models
re_hier <- ranef(m_hier)$home_team[, , "Intercept"] %>%
  as_tibble(.name_repair = "minimal") %>%
  setNames(c("estimate", "est_error", "q2.5", "q97.5")) %>%
  mutate(team = rownames(ranef(m_hier)$home_team[, , "Intercept"]))

# No pooling fixed effects: home_team coefficients (approx comparison)
fix_nopool <- fixef(m_nopool) %>% as.data.frame() %>% rownames_to_column("term") %>%
  filter(str_starts(term, "home_team")) %>%
  mutate(team = str_remove(term, "home_team")) %>%
  transmute(team, estimate = Estimate, q2.5 = Q2.5, q97.5 = Q97.5)

# Join and compare
comp <- re_hier %>%
  left_join(fix_nopool, by = "team", suffix = c("_hier", "_nopool"))

comp %>%
  ggplot(aes(x = estimate_nopool, y = estimate_hier)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(
    x = "No pooling estimate (fixed effects)",
    y = "Partial pooling estimate (hierarchical)",
    title = "Shrinkage: hierarchical estimates pull extreme values toward the mean"
  )

#Derive predictions from posterior distributions
# Create a small set of future fixtures (example)
new_matches <- tibble(
  home_team = factor(c("Team_1", "Team_2", "Team_3"), levels = levels(matches$home_team)),
  away_team = factor(c("Team_4", "Team_5", "Team_6"), levels = levels(matches$away_team))
)
# Posterior expected goals (lambda) for home_goals model
epred <- posterior_epred(m_hier, newdata = new_matches, ndraws = 2000)
# epred is draws x rows. Summarize mean and interval per match:
pred_summary <- apply(epred, 2, function(x) {
  c(mean = mean(x), q10 = quantile(x, 0.10), q90 = quantile(x, 0.90))
}) %>% t() %>% as_tibble()
bind_cols(new_matches, pred_summary)

#Simulated Goal Counts
yrep <- posterior_predict(m_hier, newdata = new_matches, ndraws = 2000)

sim_summary <- apply(yrep, 2, function(x) {
  c(mean = mean(x), q10 = quantile(x, 0.10), q90 = quantile(x, 0.90))
}) %>% t() %>% as_tibble()

bind_cols(new_matches, sim_summary)

#Model Comparison
loo_pool   <- loo(m_pool)
loo_nopool <- loo(m_nopool)
loo_hier   <- loo(m_hier)
loo_compare(loo_pool, loo_nopool, loo_hier)

# Example: informative intercept prior based on typical goals per match
# If average home goals ~ 1.4, then log(1.4) ~ 0.336
priors_sporty <- c(
  prior(normal(log(1.4), 0.5), class = "Intercept"),
  prior(exponential(1.0), class = "sd")
)

m_hier2 <- brm(
  home_goals ~ 1 + (1 | home_team) + (1 | away_team),
  data   = matches,
  family = poisson(),
  prior  = priors_sporty,
  chains = 4, cores = 4, iter = 2500, seed = 135
)

summary(m_hier2)

#Add days of rest as covariate for the model
matches2 <- matches %>%
  mutate(rest_diff = rnorm(n(), 0, 1))  # placeholder for real engineered feature

m_cov <- brm(
  home_goals ~ 1 + rest_diff + (1 | home_team) + (1 | away_team),
  data   = matches2,
  family = poisson(),
  prior  = c(
    prior(normal(0, 1.0), class = "Intercept"),
    prior(normal(0, 0.3), class = "b"),   # effect size prior for covariate
    prior(exponential(1.0), class = "sd")
  ),
  chains = 4, cores = 4, iter = 2500, seed = 135
)
summary(m_cov)