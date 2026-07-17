#Load libraries

#Load libararies
library(dplyr)
library(goalmodel)
library(footBayes)

#World Cup Final Model
expg1 <- 2.02*0.95

# away goals: multiply each column index by column sums
expg2 <- 1.86*0.72
p1x2(expg1, expg2)
spain_prob <- 0.5102634/(0.5102634  +  0.2675816)
arg_prob <- 1-spain_prob
#Probability of both teams to score
pbtts(expg1,expg2)

#Dixon Cole From First Principles

#Part 1
# Build match-specific lambdas from xG scored and conceded (neutral venue)
build_lambdas_neutral <- function(xG_for_A, xG_conc_A,
                                  xG_for_B, xG_conc_B,
                                  mu) {
  # Attack strengths
  att_A <- xG_for_A / mu
  att_B <- xG_for_B / mu
  
  # Defence strengths (higher = more goals conceded, i.e. weaker defence)
  def_A <- xG_conc_A / mu
  def_B <- xG_conc_B / mu
  
  # Neutral venue: no home-advantage multiplier
  lambda1 <- mu * att_A * def_B
  lambda2 <- mu * att_B * def_A
  
  list(lambda1 = lambda1,
       lambda2 = lambda2,
       att_A = att_A, def_A = def_A,
       att_B = att_B, def_B = def_B)
}

#Part 2
# Dixon–Coles tau adjustment
dc_tau <- function(x, y, lambda1, lambda2, rho) {
  if (x == 0 && y == 0) {
    1 - (lambda1 * lambda2 * rho)
  } else if (x == 0 && y == 1) {
    1 + lambda1 * rho
  } else if (x == 1 && y == 0) {
    1 + lambda2 * rho
  } else if (x == 1 && y == 1) {
    1 - rho
  } else {
    1
  }
}

# Bivariate DC probability for a single scoreline
dc_score_prob <- function(x, y, lambda1, lambda2, rho) {
  p_indep <- dpois(x, lambda1) * dpois(y, lambda2)
  p_indep * dc_tau(x, y, lambda1, lambda2, rho)
}

# Full probability matrix
dc_prob_matrix <- function(lambda1, lambda2, rho, max_goals = 10) {
  scores <- 0:max_goals
  n <- length(scores)
  P <- matrix(0, nrow = n, ncol = n,
              dimnames = list(paste0("T1_", scores),
                              paste0("T2_", scores)))
  
  for (i in seq_along(scores)) {
    for (j in seq_along(scores)) {
      x <- scores[i]; y <- scores[j]
      P[i, j] <- dc_score_prob(x, y, lambda1, lambda2, rho)
    }
  }
  
  P / sum(P)
}

# Aggregate to match outcomes
dc_match_outcomes <- function(lambda1, lambda2, rho, max_goals = 10) {
  P <- dc_prob_matrix(lambda1, lambda2, rho, max_goals)
  scores <- 0:max_goals
  
  p_T1_win <- 0; p_draw <- 0; p_T2_win <- 0
  
  for (i in seq_along(scores)) {
    for (j in seq_along(scores)) {
      if (scores[i] > scores[j]) {
        p_T1_win <- p_T1_win + P[i, j]
      } else if (scores[i] == scores[j]) {
        p_draw   <- p_draw + P[i, j]
      } else {
        p_T2_win <- p_T2_win + P[i, j]
      }
    }
  }
  
  list(p_T1_win = p_T1_win,
       p_draw   = p_draw,
       p_T2_win = p_T2_win,
       prob_matrix = P)
}

#WC 2026 Predictions
# Example team stats (per match averages)
xG_for_Spain  <- 2.02
xG_conc_Spain <- 0.72
xG_for_Arg  <- 1.86
xG_conc_Arg <- 0.95

mu <- 1.35   # league average goals per team per match
rho <- -0.1  # DC correlation parameter

# Build lambdas from scored/conceded, neutral venue
lam <- build_lambdas_neutral(xG_for_Spain, xG_conc_Spain,
                             xG_for_Arg, xG_conc_Arg,
                             mu)

lam$lambda1
lam$lambda2

# Dixon–Coles outcome probabilities
res <- dc_match_outcomes(lam$lambda1, lam$lambda2, rho, max_goals = 8)
res$p_T1_win
res$p_draw
res$p_T2_win 
