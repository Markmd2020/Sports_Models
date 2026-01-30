#####VolleyBall Analytics Examples#####

#Load libraries
library(tidyverse)
library(janitor)
library(broom)
library(stringr)
library(tidymodels)
library(lme4)
library(brms) 

#Set seed to ensure reproducibility
set.seed(135)

# Number of synthetic rows
n <- 5000

# Helper vectors
teams <- c("A", "B")
skills <- c("Serve", "Reception", "Set", "Attack", "Block", "Dig", "Freeball")
evaluations <- c("Error", "Poor", "OK", "Good", "Perfect")

# Generate synthetic dataset
events <- data.frame(
  match_id   = sample(1:40, n, replace = TRUE),
  set_no     = sample(1:5, n, replace = TRUE),
  rally_id   = sample(1:200, n, replace = TRUE),
  team       = sample(teams, n, replace = TRUE),
  opponent   = sample(teams, n, replace = TRUE),
  player     = paste0("Player_", sample(1:30, n, replace = TRUE)),
  skill      = sample(skills, n, replace = TRUE),
  evaluation = sample(evaluations, n, replace = TRUE),
  start_zone = sample(c(1:9, NA), n, replace = TRUE, prob = c(rep(0.1, 9), 0.1))
)

# Introduce some missing players to mimic real data issues
missing_idx <- sample(1:n, size = 0.05*n)
events$player[missing_idx] <- NA

# Introduce some empty opponent fields
empty_idx <- sample(1:n, size = 0.03*n)
events$opponent[empty_idx] <- ""

head(events)

#Data Quality Checks
# Basic validation
stopifnot(all(c("match_id","set_no","rally_id","team","skill","evaluation") %in% names(events)))
# Remove obvious duplicates (same match/set/rally/team/player/skill)
events <- events %>%
  distinct(match_id, set_no, rally_id, team, player, skill, evaluation, .keep_all = TRUE)
# Ensure opponent field exists
events <- events %>%
  mutate(opponent = if_else(is.na(opponent) | opponent == "",
                            NA_character_, opponent))
# Quick data quality report
quality_report <- list(
  n_rows = nrow(events),
  n_matches = n_distinct(events$match_id),
  missing_player = mean(is.na(events$player) | events$player == ""),
  missing_zone = mean(is.na(events$start_zone)),
  skill_counts = events %>% count(skill, sort = TRUE)
)
quality_report

#Run Rally Event function
derive_rally_context <- function(df) {
  df %>%
    group_by(match_id, set_no, rally_id) %>%
    mutate(
      serving_team = team[which(skill == "serve")[1]],
      receiving_team = setdiff(unique(team), serving_team)[1],
      phase = case_when(
        team == receiving_team ~ "sideout",
        team == serving_team   ~ "transition",
        TRUE ~ NA_character_
      ) %>% factor(levels = c("sideout","transition"))
    ) %>%
    ungroup()
}

derive_rally_context <- function(df) {
  df %>%
    group_by(match_id, set_no, rally_id) %>%
    mutate(
      serving_team = team[which(skill == "serve")[1]],
      receiving_team = setdiff(unique(team), serving_team)[1],
      phase = case_when(
        team == receiving_team ~ "sideout",
        team == serving_team   ~ "transition",
        TRUE ~ NA_character_
      ) %>% factor(levels = c("sideout","transition"))
    ) %>%
    ungroup()
}

#Apply function
derive_rally_context(events)

#Serve target heatmap
serve_zones <- events %>%
  filter(skill == "serve") %>%
  count(team, end_zone, name = "serves") %>%
  group_by(team) %>%
  mutate(pct = serves / sum(serves)) %>%
  ungroup()
ggplot(serve_zones, aes(x = factor(end_zone), y = pct)) +
  geom_col() +
  facet_wrap(~ team) +
  labs(
    title = "Serve Target Distribution by Zone",
    x = "End Zone (Serve Target)",
    y = "Share of Serves"
  )

#Plotting shot chart
shot_chart <- events %>%
  filter(skill == "attack") %>%
  mutate(
    outcome = case_when(
      tolower(evaluation) %in% eval_map$attack$kill ~ "kill",
      tolower(evaluation) %in% eval_map$attack$error ~ "error",
      tolower(evaluation) %in% eval_map$attack$blocked ~ "blocked",
      TRUE ~ "in_play"
    )
  )

ggplot(shot_chart, aes(x = factor(end_zone), fill = outcome)) +
  geom_bar(position = "fill") +
  facet_wrap(~ player) +
  labs(
    title = "Attack Outcome Mix by Target Zone (End Zone)",
    x = "Target Zone",
    y = "Share"
  )

#Create a logistic regression model to predict rallies
# Create a rally-level modeling table
rally_model_df <- events %>%
  group_by(match_id, set_no, rally_id) %>%
  summarise(
    serving_team = team[which(skill == "serve")[1]],
    receiving_team = setdiff(unique(team), serving_team)[1],
    pass_eval = evaluation[which(skill == "pass" & team == receiving_team)[1]],
    pass_score = case_when(
      tolower(pass_eval) %in% eval_map$pass$perfect ~ 3,
      tolower(pass_eval) %in% eval_map$pass$positive ~ 2,
      tolower(pass_eval) %in% eval_map$pass$negative ~ 1,
      tolower(pass_eval) %in% eval_map$pass$error ~ 0,
      TRUE ~ NA_real_
    ),
    serve_zone = end_zone[which(skill == "serve")[1]],
    point_won_by = first(na.omit(point_won_by)),
    .groups = "drop"
  ) %>%
  filter(!is.na(pass_score), !is.na(serve_zone)) %>%
  mutate(
    sideout_success = point_won_by == receiving_team
  )

# Baseline xSO model
xso_fit <- glm(
  sideout_success ~ pass_score + factor(serve_zone),
  data = rally_model_df,
  family = binomial()
)

tidy(xso_fit)
summary(xso_fit)

rally_model_df <- rally_model_df %>%
  mutate(xSO = predict(xso_fit, type = "response"))

rally_model_df %>%
  group_by(receiving_team) %>%
  summarise(
    actual_SO = mean(sideout_success),
    expected_SO = mean(xSO),
    delta = actual_SO - expected_SO,
    .groups = "drop"
  ) %>%
  arrange(desc(delta))

#Simple set model from probability score differential
# If you have event-level score columns, you can build a win probability model.
# Here we illustrate a simple logistic model from score differential and set number.

wp_df <- events %>%
  filter(!is.na(score_team), !is.na(score_opp)) %>%
  mutate(score_diff = score_team - score_opp) %>%
  group_by(match_id, set_no, rally_id) %>%
  summarise(
    team = first(team),
    score_diff = first(score_diff),
    point_won_by = first(na.omit(point_won_by)),
    .groups = "drop"
  ) %>%
  mutate(won_point = point_won_by == team)

wp_fit <- glm(won_point ~ score_diff + factor(set_no), data = wp_df, family = binomial())
wp_df <- wp_df %>%
  mutate(win_prob_point = predict(wp_fit, type = "response"))

#ELO based models
# Minimal Elo example (team-level). You can replace with your season match table.
matches <- tibble(
  match_id = c("m1","m2","m3"),
  date = as.Date(c("2025-09-01","2025-09-05","2025-09-10")),
  home = c("Team A","Team B","Team A"),
  away = c("Team B","Team C","Team C"),
  winner = c("Team A","Team C","Team A")
)

elo_update <- function(r_home, r_away, home_won, k = 20) {
  p_home <- 1 / (1 + 10^((r_away - r_home)/400))
  s_home <- ifelse(home_won, 1, 0)
  r_home_new <- r_home + k * (s_home - p_home)
  r_away_new <- r_away + k * ((1 - s_home) - (1 - p_home))
  list(home = r_home_new, away = r_away_new, p_home = p_home)
}

teams <- sort(unique(c(matches$home, matches$away)))
ratings <- setNames(rep(1500, length(teams)), teams)

elo_log <- vector("list", nrow(matches))

for (i in seq_len(nrow(matches))) {
  m <- matches[i,]
  rH <- ratings[[m$home]]
  rA <- ratings[[m$away]]
  upd <- elo_update(rH, rA, home_won = (m$winner == m$home))
  ratings[[m$home]] <- upd$home
  ratings[[m$away]] <- upd$away
  elo_log[[i]] <- tibble(match_id = m$match_id, p_home = upd$p_home,
                         home = m$home, away = m$away,
                         winner = m$winner,
                         r_home_pre = rH, r_away_pre = rA,
                         r_home_post = upd$home, r_away_post = upd$away)
}

bind_rows(elo_log) %>% arrange(match_id)
tibble(team = names(ratings), elo = as.numeric(ratings)) %>% arrange(desc(elo))

#Markov Chain Model For Rally Sequences
# Build simple sequences per rally: skill chain for receiving team until point ends
rally_sequences <- events %>%
  arrange(match_id, set_no, rally_id) %>%
  group_by(match_id, set_no, rally_id) %>%
  summarise(
    serving_team = team[which(skill == "serve")[1]],
    receiving_team = setdiff(unique(team), serving_team)[1],
    seq = paste(skill, collapse = "-"),
    point_won_by = first(na.omit(point_won_by)),
    .groups = "drop"
  )

# Count bigrams (transitions) from sequences
extract_bigrams <- function(seq_str) {
  tokens <- str_split(seq_str, "-", simplify = TRUE)
  tokens <- tokens[tokens != ""]
  if (length(tokens) < 2) return(tibble(from = character(), to = character()))
  tibble(from = tokens[-length(tokens)], to = tokens[-1])
}

transitions <- rally_sequences %>%
  mutate(bigrams = map(seq, extract_bigrams)) %>%
  select(match_id, bigrams) %>%
  unnest(bigrams) %>%
  count(from, to, name = "n") %>%
  group_by(from) %>%
  mutate(p = n / sum(n)) %>%
  ungroup() %>%
  arrange(from, desc(p))

transitions

wp_fit %>% broom::tidy()

rally_model_df

#Tidymodels workflow

df <- rally_model_df %>%
  mutate(
    serve_zone = factor(serve_zone),
    receiving_team = factor(receiving_team)
  )

split <- initial_split(df, prop = 0.8, strata = sideout_success)
train <- training(split)
test  <- testing(split)

rec <- recipe(sideout_success ~ pass_score + serve_zone, data = train) %>%
  step_impute_median(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors())

model <- logistic_reg() %>%
  set_engine("glm")

wf <- workflow() %>%
  add_recipe(rec) %>%
  add_model(model)

fit <- wf %>% fit(data = train)

pred <- predict(fit, test, type = "prob") %>%
  bind_cols(test %>% select(sideout_success))

roc_auc(pred, truth = sideout_success, .pred_TRUE)
accuracy(predict(fit, test) %>% bind_cols(test), truth = sideout_success, estimate = .pred_class)

# Example: include receiving_team as a random intercept
xso_glmm <- glmer(
  sideout_success ~ pass_score + factor(serve_zone) + (1 | receiving_team),
  data = rally_model_df,
  family = binomial()
)
summary(xso_glmm)

#Bayesian Model
#Generate synthetic dataset
rally_model_df <- data.frame(
  sideout_success = rbinom(n, 1, 0.55),                 # typical sideout rate
  pass_score       = sample(1:5, n, replace = TRUE),    # 1–5 pass quality
  serve_zone       = sample(1:9, n, replace = TRUE),    # zones 1–9
  receiving_team   = factor(sample(paste0("Team", 1:12), 
                                   n, replace = TRUE))  # 12 synthetic teams
)


bayes_fit <- brm(
  sideout_success ~ pass_score + factor(serve_zone) + (1 | receiving_team),
  data = rally_model_df,
  family = bernoulli(),
  chains = 2, cores = 2, iter = 1500,
  seed = 135
)
summary(bayes_fit)
posterior_summary(bayes_fit) 

#Bootstrap Confidence Intervals
bootstrap_ci <- function(x, B = 2000, conf = 0.95) {
  n <- length(x)
  boots <- replicate(B, mean(sample(x, n, replace = TRUE)))
  alpha <- (1 - conf) / 2
  quantile(boots, probs = c(alpha, 1 - alpha), na.rm = TRUE)
}
so_ci <- rallies %>%
  mutate(sideout_success = point_won_by == receiving_team) %>%
  group_by(receiving_team) %>%
  summarise(
    so = mean(sideout_success),
    ci_low = bootstrap_ci(sideout_success)[1],
    ci_high = bootstrap_ci(sideout_success)[2],
    n = n(),
    .groups = "drop"
  )
so_ci
