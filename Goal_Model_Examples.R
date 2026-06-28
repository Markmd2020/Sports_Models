#Load libraries

#Load libararies
library(dplyr)
library(goalmodel)
library(ggplot2)
library(MASS)
library(readr)

#Set seed to ensure reproducibility
set.seed(135)

teams <- c(
  "Arsenal","Aston Villa","Brentford","Brighton","Chelsea","Crystal Palace",
  "Everton","Fulham","Liverpool","Manchester City","Manchester United",
  "Newcastle","Nottingham Forest","Tottenham","West Ham","Wolves"
)

fixtures <- expand.grid(
  home = teams,
  away = teams,
  stringsAsFactors = FALSE
) %>% 
  filter(home != away)

n_matches <- nrow(fixtures)

fixtures <- fixtures %>%
  mutate(
    date   = as.Date("2025-08-01") + sample(0:250, n_matches, replace = TRUE),
    goals1 = rpois(n_matches, 1.4),
    goals2 = rpois(n_matches, 1.2)
  )

m_dc <- goalmodel(
  goals1 = fixtures$goals1,
  goals2 = fixtures$goals2,
  team1  = fixtures$home,
  team2  = fixtures$away,
  dc = TRUE
)

summary(m_dc)
any(fixtures$home == "Arsenal" | fixtures$away == "Arsenal")
any(fixtures$home == "Chelsea" | fixtures$away == "Chelsea")

predict_goals(m_dc, "Arsenal", "Chelsea")


# 1. Get expected goals for the matchup
pg <- predict_goals(m_dc, team1 = "Arsenal", team2 = "Chelsea")

str(pg)
# Typically a 1x2 matrix or vector: c(expg1, expg2)

# 2. Extract expg1 and expg2
pgmat <- pg[[1]]

# home goals: multiply each row index by row sums
expg1 <- sum(rowSums(pgmat) * (0:(nrow(pgmat)-1)))

# away goals: multiply each column index by column sums
expg2 <- sum(colSums(pgmat) * (0:(ncol(pgmat)-1)))


# 3. Compute 1X2 probabilities
p1x2(expg1, expg2)

predict_btts(m_dc, "Arsenal", "Chelsea")
predict_ou(m_dc, "Arsenal", "Chelsea", 2.5)

score_predictions(m_dc, fixtures)

#League Table Predictions
league_table(
  goals1 = fixtures$goals1,
  goals2 = fixtures$goals2,
  team1  = fixtures$home,
  team2  = fixtures$away
)


#Plot predictions
ppg <- predict_goals(m_dc, "Arsenal", "Chelsea")

pg <- predict_goals(m_dc, "Arsenal", "Chelsea")

# Extract the matrix from the list
probmat <- pg[[1]]

# Truncate to 5x5 (goals 0 to 4)
probmat5 <- probmat[1:5, 1:5]

df_plot <- expand.grid(
  home_goals = 0:4,
  away_goals = 0:4
)

df_plot$prob  <- c(probmat5)
df_plot$label <- sprintf("%.3f", df_plot$prob)


ggplot(df_plot, aes(home_goals, away_goals, fill = prob, label = label)) +
  geom_tile() +
  geom_text() +
  scale_fill_continuous(low = "white", high = "orange") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.ticks = element_blank())

#Negative Bimonail Model Workflow
fixtures$goals1 <- rnegbin(nrow(fixtures), mu = 1.4, theta = 1.0)
fixtures$goals2 <- rnegbin(nrow(fixtures), mu = 1.2, theta = 1.0)

# Fit NB model (DC must be FALSE)
m_nb <- goalmodel(
  goals1 = fixtures$goals1,
  goals2 = fixtures$goals2,
  team1  = fixtures$home,
  team2  = fixtures$away,
  model  = "negbin",
  dc     = FALSE
)

# Fit NB model
m_nb_dc <- goalmodel(
  goals1 = fixtures$goals1,
  goals2 = fixtures$goals2,
  team1  = fixtures$home,
  team2  = fixtures$away,
  model  = "negbin",
  dc     = FALSE
)

# Predict scoreline probabilities
pg <- predict_goals(m_nb_dc, "Arsenal", "Chelsea")
probmat <- pg[[1]]

# Truncate to 5x5
probmat5 <- probmat[1:5, 1:5]

# Expected goals
expg1_nb <- sum(rowSums(probmat5) * 0:4)
expg2_nb <- sum(colSums(probmat5) * 0:4)

# 1X2 probabilities
p1x2(expg1_nb, expg2_nb)

# Plot
df_plot <- expand.grid(
  home_goals = 0:4,
  away_goals = 0:4
)

df_plot$prob  <- c(probmat5)
df_plot$label <- sprintf("%.3f", df_plot$prob)


ggplot(df_plot, aes(home_goals, away_goals, fill = prob, label = label)) +
  geom_tile() +
  geom_text() +
  scale_fill_continuous(low = "white", high = "orange") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.ticks = element_blank())

#Workflow where data is retrieved online
library(readr)
library(dplyr)
library(goalmodel)
library(ggplot2)

# 1. Download Premier League data
url <- "https://www.football-data.co.uk/mmz4281/2324/E0.csv"
raw <- read_csv(url)

# 2. Prepare dataset
df <- raw %>%
  transmute(
    date   = as.Date(Date, format = "%d/%m/%Y"),
    home   = HomeTeam,
    away   = AwayTeam,
    goals1 = FTHG,
    goals2 = FTAG
  ) %>%
  filter(!is.na(goals1), !is.na(goals2))

# 3. Fit Dixon–Coles model
m_dc <- goalmodel(
  goals1 = df$goals1,
  goals2 = df$goals2,
  team1  = df$home,
  team2  = df$away,
  dc     = TRUE
)

# 4. Predict scoreline probabilities
pg <- predict_goals(m_dc, "Arsenal", "Chelsea")
probmat <- pg[[1]]

# 5. Truncate to 5x5
probmat5 <- probmat[1:5, 1:5]

# 6. Expected goals
expg1 <- sum(rowSums(probmat5) * 0:4)
expg2 <- sum(colSums(probmat5) * 0:4)

# 7. 1X2 probabilities
p1x2(expg1, expg2)

# 8. Plot
df_plot <- expand.grid(
  home_goals = 0:4,
  away_goals = 0:4
)

df_plot$prob  <- c(probmat5)
df_plot$label <- sprintf("%.3f", df_plot$prob)

ggplot(df_plot, aes(home_goals, away_goals, fill = prob, label = label)) +
  geom_tile() +
  geom_text() +
  scale_fill_continuous(low = "white", high = "orange") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.ticks = element_blank())
