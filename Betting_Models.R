#Sports Betting Models

#Load libraries
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(lubridate)
library(ggplot2)
library(scales)
library(rstanarm)
library(posterior)

#Generate sythetic Data
simulate_league <- function(
    n_teams = 12,
    n_games = 2000,
    home_adv = 0.25,
    sigma_strength = 0.9,
    bookmaker_margin = 0.05
) {
  teams <- paste0("T", sprintf("%02d", 1:n_teams))
  strength <- rnorm(n_teams, 0, sigma_strength)
  names(strength) <- teams
  
  team <- sample(teams, n_games, replace = TRUE)
  opponent <- sample(teams, n_games, replace = TRUE)
  while(any(team == opponent)) {
    idx <- which(team == opponent)
    opponent[idx] <- sample(teams, length(idx), replace = TRUE)
  }
  
  home <- rbinom(n_games, 1, 0.5)
  date <- as.Date("2022-01-01") + sample(0:900, n_games, replace = TRUE)
  
  lin <- (strength[team] - strength[opponent]) + home_adv * home
  p_true <- plogis(lin)
  
  p_mkt <- pmin(pmax(p_true + rnorm(n_games, 0, 0.03), 0.02), 0.98)
  p_raw <- pmin(pmax(p_mkt * (1 + bookmaker_margin), 0.02), 0.995)
  odds_decimal <- 1 / p_raw
  
  result <- rbinom(n_games, 1, p_true)
  
  tibble(
    date = date,
    team = team,
    opponent = opponent,
    home = home,
    result = result,
    p_true = p_true,
    odds_decimal = odds_decimal
  ) %>%
    arrange(date)
}

matches <- simulate_league()
glimpse(matches)
summary(matches$odds_decimal)

#Bayesian Foundation
# Prior: Beta(a, b)
a0 <- 10
b0 <- 10

# Observed: k wins out of n
k <- 12
n <- 20

# Posterior: Beta(a0 + k, b0 + n - k)
a1 <- a0 + k
b1 <- b0 + (n - k)

# Posterior mean and interval
post_mean <- a1 / (a1 + b1)
post_ci <- qbeta(c(0.05, 0.95), a1, b1)

c(mean = post_mean, lo90 = post_ci[1], hi90 = post_ci[2])

# Visualize prior vs posterior
grid <- seq(0, 1, length.out = 400)

df_beta <- bind_rows(
  tibble(p = grid, density = dbeta(grid, a0, b0), dist = "prior"),
  tibble(p = grid, density = dbeta(grid, a1, b1), dist = "posterior")
)

ggplot(df_beta, aes(p, density, color = dist)) +
  geom_line() +
  labs(
    title = "Beta Prior vs Posterior",
    x = "Win probability p",
    y = "Density"
  )

#Bayesian Win Rate Conjugate Model
# Build per-team counts
team_stats <- matches %>%
  group_by(team) %>%
  summarise(
    n = n(),
    k = sum(result),
    .groups = "drop"
  )

# Choose a prior: mildly skeptical
a0 <- 8
b0 <- 8

team_post <- team_stats %>%
  mutate(
    a = a0 + k,
    b = b0 + (n - k),
    post_mean = a/(a+b),
    lo90 = qbeta(0.05, a, b),
    hi90 = qbeta(0.95, a, b)
  ) %>%
  arrange(desc(post_mean))

team_post %>% slice(1:10)
ggplot(team_post, aes(reorder(team, post_mean), post_mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = lo90, ymax = hi90), width = 0.2) +
  coord_flip() +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(
    title = "Team Win Probability Posteriors (Beta-Binomial)",
    x = "Team",
    y = "Posterior mean (90% interval)"
  )

#Posterior Sampling
# Suppose you want to bet Team X at odds d
team_x <- team_post$team[1]
d <- 1.95

a <- team_post$a[team_post$team == team_x]
b <- team_post$b[team_post$team == team_x]

#Expected Value Function
ev_decimal <- function(p, d) {
  # expected profit per unit stake for decimal odds d
  # win: profit = d - 1; lose: -1
  p*(d - 1) - (1 - p)
}

p_draws <- rbeta(20000, a, b)
ev_draws <- ev_decimal(p_draws, d)

quantile(ev_draws, c(0.05, 0.5, 0.95))

# Probability that EV is positive (posterior probability of an edge)
mean(ev_draws > 0)

#Sequential Updating
update_beta <- function(a, b, y) {
  # y is 1 for win, 0 for loss
  c(a = a + y, b = b + (1 - y))
}

sequential_posterior <- function(y, a0 = 8, b0 = 8) {
  a <- a0; b <- b0
  out <- vector("list", length(y))
  for (i in seq_along(y)) {
    ab <- update_beta(a, b, y[i])
    a <- ab["a"]; b <- ab["b"]
    out[[i]] <- c(
      i = i,
      a = a, b = b,
      mean = a/(a+b),
      lo90 = qbeta(0.05, a, b),
      hi90 = qbeta(0.95, a, b)
    )
  }
  out
}

# Pick one team and track posterior over time
team_pick <- sample(unique(matches$team), 1)
y_seq <- matches %>% 
  filter(team == team_pick) %>% 
  arrange(date) %>% 
  pull(result)

seq_df <- sequential_posterior(y_seq, a0 = 8, b0 = 8) %>%
  mutate(team = team_pick)

ggplot(seq_df, aes(i, mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = lo90, ymax = hi90), alpha = 0.2) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(
    title = paste("Sequential Posterior for", team_pick),
    x = "Game index",
    y = "Posterior mean (90% band)"
  )

#Empirical Bayes Function
fit_beta_mom <- function(p_hat, w = NULL) {
  # method of moments for Beta parameters
  # p_hat: observed proportions, w optional weights (e.g., n games)
  if (is.null(w)) w <- rep(1, length(p_hat))
  
  m <- weighted.mean(p_hat, w)
  v <- weighted.mean((p_hat - m)^2, w)
  
  # clamp variance to avoid weirdness
  v <- max(v, 1e-6)
  
  # Beta variance: m(1-m)/(a+b+1)
  t <- m*(1-m)/v - 1
  a <- max(0.5, m*t)
  b <- max(0.5, (1-m)*t)
  c(a = a, b = b)
}

p_hat <- team_stats$k / team_stats$n
prior_ab <- fit_beta_mom(p_hat, w = team_stats$n)
prior_ab

a0 <- prior_ab["a"]
b0 <- prior_ab["b"]

team_post_eb <- team_stats %>%
  mutate(
    a = a0 + k,
    b = b0 + (n - k),
    post_mean = a/(a+b),
    lo90 = qbeta(0.05, a, b),
    hi90 = qbeta(0.95, a, b)
  ) %>%
  arrange(desc(post_mean))

# Compare naive sample rates vs EB posterior means
compare_df <- team_stats %>%
  mutate(sample_rate = k/n) %>%
  left_join(team_post_eb %>% select(team, post_mean), by = "team") %>%
  pivot_longer(cols = c(sample_rate, post_mean), names_to = "type", values_to = "p")

ggplot(compare_df, aes(p, fill = type)) +
  geom_histogram(bins = 20, alpha = 0.6, position = "identity") +
  labs(
    title = "Shrinkage Effect: Sample Rates vs EB Posterior Means",
    x = "Probability",
    y = "Count"
  )

#Bayesian Logistic Regression
# Use rolling win rate as a proxy strength feature (in real data: Elo, xG, etc.)
matches_feat <- matches %>%
  arrange(date) %>%
  group_by(team) %>%
  mutate(
    games_played = row_number() - 1,
    roll_winrate = ifelse(games_played == 0, NA_real_,
                          cummean(lag(result)))
  ) %>%
  ungroup() %>%
  mutate(
    roll_winrate = ifelse(is.na(roll_winrate), mean(result), roll_winrate),
    # center features
    roll_winrate_c = roll_winrate - mean(roll_winrate),
    home_c = home - mean(home)
  )
str(matches_feat$date)
unique(matches_feat$date)

# Train-test split by time (avoid leakage)
cut_date <- quantile(matches_feat$date, 0.8)
train <- matches_feat %>% filter(date <= "2023-11-25")
test  <- matches_feat %>% filter(date >  "2023-11-25")

nrow(train); nrow(test)
summary(train)

fit_bayes_logit <- stan_glm(
  result ~ home_c ,
  data = train,
  family = binomial(link = "logit"),
  prior = normal(location = 0, scale = 1, autoscale = TRUE),
  prior_intercept = normal(0, 2.5),
  chains = 2,
  iter = 800,
  refresh = 0
)

print(fit_bayes_logit)

# Posterior predictive probabilities on test
p_test <- posterior_linpred(fit_bayes_logit, newdata = test, transform = TRUE)
# p_test is draws x observations; take posterior mean probability
p_hat <- colMeans(p_test)

test_pred <- test %>%
  mutate(p_model = p_hat)

test_pred %>%
  summarise(
    mean_p = mean(p_model),
    mean_y = mean(result)
  )

#Elo rating
elo_update <- function(Ra, Rb, y, k = 20) {
  # y: 1 if A wins, 0 if B wins
  Ea <- 1 / (1 + 10^((Rb - Ra)/400))
  Ra_new <- Ra + k*(y - Ea)
  Rb_new <- Rb + k*((1 - y) - (1 - Ea))
  c(Ra = Ra_new, Rb = Rb_new, Ea = Ea)
}

run_elo <- function(df, init = 1500, k = 20) {
  teams <- sort(unique(c(df$team, df$opponent)))
  R <- setNames(rep(init, length(teams)), teams)
  
  out <- vector("list", nrow(df))
  for (i in seq_len(nrow(df))) {
    a <- df$team[i]
    b <- df$opponent[i]
    y <- df$result[i]
    
    upd <- elo_update(R[a], R[b], y, k = k)
    R[a] <- upd["Ra"]
    R[b] <- upd["Rb"]
    
    out[[i]] <- tibble(
      date = df$date[i],
      team = a,
      opponent = b,
      result = y,
      elo_team = R[a],
      elo_opp = R[b],
      p_elo = upd["Ea"]
    )
  }
  bind_rows(out)
}

elo_df <- run_elo(matches %>% arrange(date), k = 18)
elo_df %>% slice(1:5)
summary(elo_df)

# Join Elo features back to matches
matches_elo <- matches %>%
  arrange(date) %>%
  mutate(row_id = row_number()) %>%
  left_join(
    elo_df %>% mutate(row_id = row_number()) %>% select(row_id, p_elo),
    by = "row_id"
  )

summary(matches_elo$p_elo)

#Skellam Distribution For Score Models
# Add synthetic goals to the toy dataset (for demonstration)
add_goals <- function(df) {
  # tie lambda to latent p_true to create plausible relationship
  # this is just for tutorial purposes
  base <- 1.2
  spread <- (df$p_true - 0.5) * 2.0
  
  lambda_home <- pmax(0.2, base + 0.4 + spread)
  lambda_away <- pmax(0.2, base - 0.1 - spread)
  
  df %>%
    mutate(
      goals_team = rpois(n(), lambda_home),
      goals_opp  = rpois(n(), lambda_away),
      total_goals = goals_team + goals_opp,
      goal_diff = goals_team - goals_opp
    )
}

matches_goals <- add_goals(matches)
summary(matches_goals$total_goals)

# Poisson model for "team goals" using simple predictors
# In real data: team/opponent effects, home advantage, pace, etc.
cut_date <- "2023-11-25"
train_g <- matches_goals %>% filter(date <= cut_date)
test_g  <- matches_goals %>% filter(date >  cut_date)

fit_pois <- glm(
  goals_team ~ home + p_true,
  data = train_g,
  family = poisson()
)

summary(fit_pois)

# Predict lambdas on test
lambda_hat <- predict(fit_pois, newdata = test_g, type = "response")
test_g2 <- test_g %>% mutate(lambda_team = lambda_hat)
summary(test_g2$lambda_team)

