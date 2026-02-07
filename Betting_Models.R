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

