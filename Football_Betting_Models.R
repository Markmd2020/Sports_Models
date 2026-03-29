############Football Betting Model######

#Load libraries
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(lubridate)
library(readr)
library(ggplot2)
library(janitor)
library(glue)
library(worldfootballR)



epl_url <- fb_league_urls(country = "ENG", gender = "M", 
                            season_end_year = 2024,tier = "1st")

matches <- fb_match_results(epl_url,season_end_year = 2024,gender = "M") %>%
  janitor::clean_names()

head(matches)

#Generate sythetic data
set.seed(123)

# Number of synthetic matches
n_matches <- 50

# Some fake Premier League teams
teams <- c(
  "Arsenal", "Aston Villa", "Brentford", "Brighton",
  "Chelsea", "Crystal Palace", "Everton", "Fulham",
  "Liverpool", "Manchester City", "Manchester United",
  "Newcastle", "Nottingham Forest", "Tottenham",
  "West Ham", "Wolves"
)

# Generate synthetic matches
n_matches <- 50

teams <- c(
  "Arsenal", "Aston Villa", "Brentford", "Brighton",
  "Chelsea", "Crystal Palace", "Everton", "Fulham",
  "Liverpool", "Manchester City", "Manchester United",
  "Newcastle", "Nottingham Forest", "Tottenham",
  "West Ham", "Wolves"
)

matches <- data.frame(
  match_id   = seq_len(n_matches),
  season     = "2025-2026",
  date       = as.Date("2025-08-01") + sample(0:250, n_matches, replace = TRUE),
  home_team  = sample(teams, n_matches, replace = TRUE),
  away_team  = sample(teams, n_matches, replace = TRUE),
  home_goals = rpois(n_matches, lambda = 1.6),
  away_goals = rpois(n_matches, lambda = 1.2),
  venue      = "Stadium X",
  stringsAsFactors = FALSE
)

# Avoid home_team == away_team by resampling those rows
same_team_idx <- which(matches$home_team == matches$away_team)
if (length(same_team_idx) > 0) {
  matches$away_team[same_team_idx] <- sample(
    teams[teams != matches$home_team[same_team_idx]],
    length(same_team_idx),
    replace = TRUE
  )
}

# Add result column from home perspective
matches$result <- with(matches, ifelse(
  home_goals > away_goals, "H",
  ifelse(home_goals < away_goals, "A", "D")
))

# Inspect
head(matches)


# Make sure we have standardized column names
matches <- matches %>%
  transmute(
    date = as.Date(date),
    season = if_else(month(date) >= 7, year(date) + 1L, year(date)), # football season heuristic
    home = as.character(home_team),
    away = as.character(away_team),
    hg = as.integer(home_goals),
    ag = as.integer(away_goals)
  ) %>%
  filter(!is.na(date), !is.na(home), !is.na(away), !is.na(hg), !is.na(ag)) %>%
  arrange(date)

# Basic sanity checks
stopifnot(all(matches$hg >= 0), all(matches$ag >= 0))

long <- matches %>%
  mutate(match_id = row_number()) %>%
  tidyr::pivot_longer(
    cols = c(home, away),
    names_to = "side",
    values_to = "team"
  ) %>%
  mutate(
    opp = if_else(side == "home", "away", "home"),
    goals = if_else(side == "home", hg, ag),
    conceded = if_else(side == "home", ag, hg),
    is_home = as.integer(side == "home")
  ) %>%
  select(match_id, date, season, team, opp, is_home, goals, conceded)
head(long)

#Add form features
# Rolling averages for goals scored/conceded over last N matches (per team)
add_form_features <- function(df, n = 5) {
  df %>%
    arrange(team, date, match_id) %>%
    group_by(team) %>%
    mutate(
      gf_roll = zoo::rollapplyr(goals, width = n, FUN = mean, fill = NA, partial = TRUE),
      ga_roll = zoo::rollapplyr(conceded, width = n, FUN = mean, fill = NA, partial = TRUE)
    ) %>%
    ungroup()
}

#Train Test split
# Time-based split (e.g., last 20% of matches as test)
n_total <- nrow(matches)
cut_idx <- floor(n_total * 0.80)

train <- matches %>% slice(1:cut_idx)
test  <- matches %>% slice((cut_idx + 1):n_total)

# Ensure consistent factor levels
teams <- sort(unique(c(matches$home, matches$away)))
train <- train %>% mutate(home = factor(home, levels = teams), away = factor(away, levels = teams))
test  <- test  %>% mutate(home = factor(home, levels = teams), away = factor(away, levels = teams))

# Fit models
home_mod <- glm(hg ~ 1 + home + away, data = train, family = poisson())
away_mod <- glm(ag ~ 1 + away + home, data = train, family = poisson())

summary(home_mod)
summary(away_mod)

# Build a modeling frame in the classic attack/defense form
# We model home goals:
# log(lambda_home) = home_adv + attack_home - defense_away
# And away goals:
# log(lambda_away) = attack_away - defense_home

# Create team factors
train2 <- train %>%
  mutate(
    home = factor(home, levels = teams),
    away = factor(away, levels = teams)
  )

# We'll encode attack and defense as separate factors by prefixing labels
mk_attack <- function(team) factor(paste0("att_", team), levels = paste0("att_", teams))
mk_def    <- function(team) factor(paste0("def_", team), levels = paste0("def_", teams))

train_home <- train2 %>%
  transmute(
    goals = hg,
    is_home = 1L,
    att = mk_attack(home),
    def = mk_def(away)
  )

train_away <- train2 %>%
  transmute(
    goals = ag,
    is_home = 0L,
    att = mk_attack(away),
    def = mk_def(home)
  )

train_long <- bind_rows(train_home, train_away)
head(train_long)

# Fit a single Poisson model with:
# goals ~ is_home + att + def
# Note: to reflect "- defense" we can include def and allow coefficients to learn direction;
# For stricter structure you can re-code defense sign, but this works well in practice.
ad_mod <- glm(goals ~ is_home + att + def, data = train_long, family = poisson())
summary(ad_mod)


#Predict expect goals for each match
predict_lambdas <- function(df, model, teams) {
  df2 <- df %>%
    mutate(
      home = factor(home, levels = teams),
      away = factor(away, levels = teams),
      att_home = factor(paste0("att_", home), levels = paste0("att_", teams)),
      def_away = factor(paste0("def_", away), levels = paste0("def_", teams)),
      att_away = factor(paste0("att_", away), levels = paste0("att_", teams)),
      def_home = factor(paste0("def_", home), levels = paste0("def_", teams))
    )
  
  # home lambda
  new_home <- df2 %>%
    transmute(is_home = 1L, att = att_home, def = def_away)
  
  # away lambda
  new_away <- df2 %>%
    transmute(is_home = 0L, att = att_away, def = def_home)
  
  lam_home <- predict(model, newdata = new_home, type = "response")
  lam_away <- predict(model, newdata = new_away, type = "response")
  
  df2 %>%
    mutate(lambda_home = lam_home, lambda_away = lam_away)
  
}
test_pred <- predict_lambdas(test, ad_mod, teams)
head(test_pred)
  
  #Dixon Cole Adjustments
  # Dixon-Coles tau adjustment for low-score dependence
  tau_dc <- function(x, y, lam_x, lam_y, rho) {
    # x = home goals, y = away goals
    # rho is the dependence parameter
    if (x == 0 && y == 0) return(1 - (lam_x * lam_y * rho))
    if (x == 0 && y == 1) return(1 + (lam_x * rho))
    if (x == 1 && y == 0) return(1 + (lam_y * rho))
    if (x == 1 && y == 1) return(1 - rho)
    return(1)
  }
  
  # Scoreline probability matrix up to max_goals
  score_matrix <- function(lam_h, lam_a, rho = 0, max_goals = 10) {
    xs <- 0:max_goals
    ys <- 0:max_goals
    
    ph <- dpois(xs, lam_h)
    pa <- dpois(ys, lam_a)
    
    # outer product for independent probabilities
    P <- outer(ph, pa)
    
    # apply DC tau correction
    for (i in seq_along(xs)) {
      for (j in seq_along(ys)) {
        P[i, j] <- P[i, j] * tau_dc(xs[i], ys[j], lam_h, lam_a, rho)
      }
    }
    
    # renormalize
    P / sum(P)
  }
  
# Example
P_ex <- score_matrix(lam_h = 1.4, lam_a = 1.1, rho = 0.05, max_goals = 8)
round(P_ex[1:5,1:5], 4)

# Estimate rho by maximizing log-likelihood on train set given lambdas
train_pred <- predict_lambdas(train, ad_mod, teams)

dc_loglik <- function(rho, df, max_goals = 10) {
  # clamp rho to a reasonable range to avoid numerical issues
  rho <- max(min(rho, 0.3), -0.3)
  
  ll <- 0
  for (k in seq_len(nrow(df))) {
    lam_h <- df$lambda_home[k]
    lam_a <- df$lambda_away[k]
    hg <- df$hg[k]
    ag <- df$ag[k]
    
    P <- score_matrix(lam_h, lam_a, rho = rho, max_goals = max_goals)
    
    # if score exceeds max_goals, treat as tiny prob (or increase max_goals)
    if (hg > max_goals || ag > max_goals) {
      ll <- ll + log(1e-12)
    } else {
      ll <- ll + log(P[hg + 1, ag + 1] + 1e-15)
    }
  }
  ll
}

opt <- optimize(
  f = function(r) -dc_loglik(r, train_pred, max_goals = 10),
  interval = c(-0.2, 0.2)
)

rho_hat <- opt$minimum
rho_hat

#From scorelines to 1X2 probabilities
p1x2_from_matrix <- function(P) {
  max_g <- nrow(P) - 1
  xs <- 0:max_g
  ys <- 0:max_g
  
  p_home <- 0
  p_draw <- 0
  p_away <- 0
  
  for (i in seq_along(xs)) {
    for (j in seq_along(ys)) {
      if (xs[i] > ys[j]) p_home <- p_home + P[i, j]
      if (xs[i] == ys[j]) p_draw <- p_draw + P[i, j]
      if (xs[i] < ys[j]) p_away <- p_away + P[i, j]
    }
  }
  
  tibble(p_home = p_home, p_draw = p_draw, p_away = p_away)
}

predict_1x2 <- function(df, rho = 0, max_goals = 10) {
  out <- vector("list", nrow(df))
  for (k in seq_len(nrow(df))) {
    P <- score_matrix(df$lambda_home[k], df$lambda_away[k], rho = rho, max_goals = max_goals)
    out[[k]] <- p1x2_from_matrix(P)
  }
  bind_rows(out)
}

test_1x2 <- bind_cols(
  test_pred,
  predict_1x2(test_pred, rho = rho_hat, max_goals = 10)
)

test_1x2 %>%
  select(date, home, away, hg, ag, lambda_home, lambda_away, p_home, p_draw, p_away) %>%
  head(10)

#Calculating Optimal Odds
set.seed(135)

# Number of synthetic matches
n_matches <- 50

# Some fake Premier League teams
teams <- c(
  "Arsenal", "Aston Villa", "Brentford", "Brighton",
  "Chelsea", "Crystal Palace", "Everton", "Fulham",
  "Liverpool", "Manchester City", "Manchester United",
  "Newcastle", "Nottingham Forest", "Tottenham",
  "West Ham", "Wolves"
)

# Generate synthetic odds dataset
odds <- data.frame(
  date       = as.Date("2025-08-01") + sample(0:250, n_matches, replace = TRUE),
  home       = sample(teams, n_matches, replace = TRUE),
  away       = sample(teams, n_matches, replace = TRUE),
  odds_home  = round(runif(n_matches, 1.3, 3.5), 2),
  odds_draw  = round(runif(n_matches, 2.8, 4.2), 2),
  odds_away  = round(runif(n_matches, 1.8, 5.0), 2),
  stringsAsFactors = FALSE
)

# Ensure home != away
same_team_idx <- which(odds$home == odds$away)
if (length(same_team_idx) > 0) {
  odds$away[same_team_idx] <- sample(
    teams[teams != odds$home[same_team_idx]],
    length(same_team_idx),
    replace = TRUE
  )
}

head(odds)

odds <- odds%>%
  janitor::clean_names() %>%
  mutate(date = as.Date(date)) %>%
  transmute(
    date,
    home = as.character(home),
    away = as.character(away),
    o_home = as.numeric(odds_home),
    o_draw = as.numeric(odds_draw),
    o_away = as.numeric(odds_away)
  )

head(odds)

df <- test_1x2 %>%
  mutate(home = as.character(home), away = as.character(away)) %>%
  left_join(odds, by = c("date","home","away"))

# Implied probs (no vig removal yet)
df <- df %>%
  mutate(
    imp_home = 1 / o_home,
    imp_draw = 1 / o_draw,
    imp_away = 1 / o_away,
    overround = imp_home + imp_draw + imp_away
  )

# Simple vig removal by normalization
df <- df %>%
  mutate(
    mkt_home = imp_home / overround,
    mkt_draw = imp_draw / overround,
    mkt_away = imp_away / overround
  )

# EV per 1 unit stake
ev <- function(p, o) p*(o - 1) - (1 - p)

df <- df %>%
  mutate(
    ev_home = ev(p_home, o_home),
    ev_draw = ev(p_draw, o_draw),
    ev_away = ev(p_away, o_away)
  )

df %>%
  select(date, home, away, p_home, p_draw, p_away, o_home, o_draw, o_away, ev_home, ev_draw, ev_away) %>%
  head(10)

#Pick bets with thresholds
# Practical filters: require at least some edge and avoid tiny probabilities
EDGE_MIN <- 0.02   # 2% EV edge
P_MIN    <- 0.05   # avoid extreme longshots unless you model them well

df_bets <- df %>%
  mutate(
    pick = case_when(
      ev_home == pmax(ev_home, ev_draw, ev_away, na.rm = TRUE) ~ "H",
      ev_draw == pmax(ev_home, ev_draw, ev_away, na.rm = TRUE) ~ "D",
      TRUE ~ "A"
    ),
    p_pick = case_when(pick == "H" ~ p_home, pick == "D" ~ p_draw, TRUE ~ p_away),
    o_pick = case_when(pick == "H" ~ o_home, pick == "D" ~ o_draw, TRUE ~ o_away),
    ev_pick = case_when(pick == "H" ~ ev_home, pick == "D" ~ ev_draw, TRUE ~ ev_away)
  ) %>%
  filter(!is.na(o_pick)) %>%
  filter(p_pick >= P_MIN, ev_pick >= EDGE_MIN)

df_bets %>% count(pick)

#Backtesting
# Practical filters: require at least some edge and avoid tiny probabilities
EDGE_MIN <- 0.02   # 2% EV edge
P_MIN    <- 0.05   # avoid extreme longshots unless you model them well

df_bets <- df %>%
  mutate(
    pick = case_when(
      ev_home == pmax(ev_home, ev_draw, ev_away, na.rm = TRUE) ~ "H",
      ev_draw == pmax(ev_home, ev_draw, ev_away, na.rm = TRUE) ~ "D",
      TRUE ~ "A"
    ),
    p_pick = case_when(pick == "H" ~ p_home, pick == "D" ~ p_draw, TRUE ~ p_away),
    o_pick = case_when(pick == "H" ~ o_home, pick == "D" ~ o_draw, TRUE ~ o_away),
    ev_pick = case_when(pick == "H" ~ ev_home, pick == "D" ~ ev_draw, TRUE ~ ev_away)
  ) %>%
  filter(!is.na(o_pick)) %>%
  filter(p_pick >= P_MIN, ev_pick >= EDGE_MIN)

df_bets %>% count(pick)

#Kelly Staking
# Practical filters: require at least some edge and avoid tiny probabilities
EDGE_MIN <- 0.02   # 2% EV edge
P_MIN    <- 0.05   # avoid extreme longshots unless you model them well

df_bets <- df %>%
  mutate(
    pick = case_when(
      ev_home == pmax(ev_home, ev_draw, ev_away, na.rm = TRUE) ~ "H",
      ev_draw == pmax(ev_home, ev_draw, ev_away, na.rm = TRUE) ~ "D",
      TRUE ~ "A"
    ),
    p_pick = case_when(pick == "H" ~ p_home, pick == "D" ~ p_draw, TRUE ~ p_away),
    o_pick = case_when(pick == "H" ~ o_home, pick == "D" ~ o_draw, TRUE ~ o_away),
    ev_pick = case_when(pick == "H" ~ ev_home, pick == "D" ~ ev_draw, TRUE ~ ev_away)
  ) %>%
  filter(!is.na(o_pick)) %>%
  filter(p_pick >= P_MIN, ev_pick >= EDGE_MIN)

df_bets %>% count(pick)

#Plot bankroll curves
plot_df <- df_bets %>%
  select(date, br_flat, br_kelly) %>%
  pivot_longer(cols = c(br_flat, br_kelly), names_to = "strategy", values_to = "bankroll")

ggplot(plot_df, aes(x = date, y = bankroll, group = strategy)) +
  geom_line() +
  labs(x = NULL, y = "Bankroll", title = "Backtest Bankroll: Flat vs Fractional Kelly")

#Calibration diagnostics
# Create a probability matrix and truth labels
df_eval <- df %>%
  filter(!is.na(p_home), !is.na(p_draw), !is.na(p_away)) %>%
  mutate(
    truth = case_when(hg > ag ~ "H", hg == ag ~ "D", TRUE ~ "A"),
    truth = factor(truth, levels = c("H","D","A"))
  )

# Log loss (manual)
log_loss_1x2 <- function(pH, pD, pA, y) {
  eps <- 1e-15
  p <- ifelse(y=="H", pH, ifelse(y=="D", pD, pA))
  -mean(log(pmax(p, eps)))
}

ll <- log_loss_1x2(df_eval$p_home, df_eval$p_draw, df_eval$p_away, df_eval$truth)
ll

#Reliability plot
# Example: calibration for HOME-win probability
calib_home <- df_eval %>%
  mutate(bin = ntile(p_home, 10)) %>%
  group_by(bin) %>%
  summarise(
    p_mean = mean(p_home),
    freq = mean(truth == "H"),
    n = n(),
    .groups = "drop"
  )

ggplot(calib_home, aes(x = p_mean, y = freq)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "Predicted P(Home win)", y = "Observed frequency", title = "Calibration (Home win)")

#Create production pipeline
# 1) Load historical matches (from worldfootballR or your CSV)
matches <- readr::read_csv("matches.csv") %>% janitor::clean_names()

# 2) Train up to cutoff date (e.g., yesterday)
cutoff_date <- Sys.Date() - 1

hist <- matches %>%
  mutate(date = as.Date(date)) %>%
  filter(date <= cutoff_date) %>%
  transmute(date, home = home_team, away = away_team, hg = home_goals, ag = away_goals) %>%
  arrange(date)

teams <- sort(unique(c(hist$home, hist$away)))

# 3) Fit attack/defense model
train_long <- bind_rows(
  hist %>%
    transmute(goals = hg, is_home = 1L,
              att = factor(paste0("att_", home), levels = paste0("att_", teams)),
              def = factor(paste0("def_", away), levels = paste0("def_", teams))),
  hist %>%
    transmute(goals = ag, is_home = 0L,
              att = factor(paste0("att_", away), levels = paste0("att_", teams)),
              def = factor(paste0("def_", home), levels = paste0("def_", teams)))
)

ad_mod <- glm(goals ~ is_home + att + def, data = train_long, family = poisson())

# 4) Predict for upcoming fixtures (you need a fixtures table)
fixtures <- readr::read_csv("fixtures.csv") %>%
  janitor::clean_names() %>%
  mutate(date = as.Date(date)) %>%
  transmute(date, home = home_team, away = away_team)

# Compute lambdas
fixtures2 <- fixtures %>%
  mutate(
    home = factor(home, levels = teams),
    away = factor(away, levels = teams)
  )

fixtures_pred <- predict_lambdas(
  df = fixtures2 %>% mutate(hg = 0L, ag = 0L), # placeholders
  model = ad_mod,
  teams = teams
) 

# 5) Estimate rho (optional) on recent history only (faster)
hist_pred <- hist %>% mutate(home = factor(home, levels=teams), away=factor(away, levels=teams))
hist_pred <- predict_lambdas(hist_pred, ad_mod, teams)

opt <- optimize(
  f = function(r) -dc_loglik(r, hist_pred %>% mutate(lambda_home=lambda_home, lambda_away=lambda_away),
                             max_goals = 10),
  interval = c(-0.2, 0.2)
)
rho_hat <- opt$minimum

# 6) Convert to 1X2
fixtures_1x2 <- bind_cols(
  fixtures_pred,
  predict_1x2(fixtures_pred, rho = rho_hat, max_goals = 10)
) %>%
  select(date, home, away, lambda_home, lambda_away, p_home, p_draw, p_away)

write_csv(fixtures_1x2, "model_probs.csv")