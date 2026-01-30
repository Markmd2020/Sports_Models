######Empirical  Bayes Examples For Basebal####

#Load libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(Lahman)
library(stats4)
library(VGAM)

#Simulation of hits examples
num_trials <- 10e6
simulations <- data_frame(
  true_average = rbeta(num_trials, 81, 219),
  hits = rbinom(num_trials, 300, true_average)
)

hit_100 <- simulations %>%
  filter(hits == 100)

#Plotting simulation
simulations %>%
  filter(hits %in% c(60, 80, 100)) %>%
  ggplot(aes(true_average, color = factor(hits))) +
  geom_density() +
  labs(x = "True average of players with H hits / 300 at-bats",
       color = "H")

# Filter out pitchers
career <- Batting %>%
  filter(AB > 0) %>%
  anti_join(Pitching, by = "playerID") %>%
  group_by(playerID) %>%
  summarize(H = sum(H), AB = sum(AB)) %>%
  mutate(average = H / AB)

# Include names along with the player IDs
career <- Master %>%
  tbl_df() %>%
  dplyr::select(playerID, nameFirst, nameLast) %>%
  unite(name, nameFirst, nameLast, sep = " ") %>%
  inner_join(career, by = "playerID") %>%
  dplyr::select(-playerID)

career_filtered <- career %>%
  filter(AB > 500)

# log-likelihood function
ll <- function(alpha, beta) {
  x <- career_filtered$H
  total <- career_filtered$AB
  -sum(VGAM::dbetabinom.ab(x, total, alpha, beta, log = TRUE))
}

# maximum likelihood estimation
m <- mle(ll, start = list(alpha = 1, beta = 10), method = "L-BFGS-B",
         lower = c(0.0001, .1))
ab <- coef(m)
alpha0 <- ab[1]
beta0 <- ab[2]

#Empirical Bayes Estimate
career_eb <- career %>%
  mutate(eb_estimate = (H + alpha0) / (AB + alpha0 + beta0))

#Credible Intervals
career <- Batting %>%
  filter(AB > 0) %>%
  anti_join(Pitching, by = "playerID") %>%
  group_by(playerID) %>%
  summarize(H = sum(H), AB = sum(AB)) %>%
  mutate(average = H / AB)

career <- Master %>%
  tbl_df() %>%
  dplyr::select(playerID, nameFirst, nameLast) %>%
  unite(name, nameFirst, nameLast, sep = " ") %>%
  inner_join(career, by = "playerID")

#Prior Distribution
# values estimated by maximum likelihood in Chapter 3
alpha0 <- 101.4
beta0 <- 287.3
career_eb <- career %>%
  mutate(eb_estimate = (H + alpha0) / (AB + alpha0 + beta0))

#Posterior Distribution
career_eb <- career_eb %>%
  mutate(alpha1 = alpha0 + H,
         beta1 = beta0 + AB- H)

#Plotting Credible Intervals
yankee_1998 <- c("brosisc01", "jeterde01", "knoblch01", "martiti02",
                 "posadjo01", "strawda01", "willibe02")

yankee_1998_career <- career_eb %>%
  filter(playerID %in% yankee_1998)

yankee_beta <- yankee_1998_career %>%
  crossing(x = seq(.18, .33, .0002)) %>%
  ungroup() %>%
  mutate(density = dbeta(x, alpha1, beta1))

ggplot(yankee_beta, aes(x, density, color = playerID)) +
  geom_line() +
  stat_function(fun = function(x) dbeta(x, alpha0, beta0),
                lty = 2, color = "black") +
  labs(x = "Batting average",
       color = "Player")  

yankee_1998_career <- yankee_1998_career %>%
  mutate(low = qbeta(.025, alpha1, beta1),
         high = qbeta(.975, alpha1, beta1))

yankee_1998_career %>%
  mutate(name = reorder(playerID, eb_estimate)) %>%
  ggplot(aes(eb_estimate, playerID)) +
  geom_point() +
  geom_errorbarh(aes(xmin = low, xmax = high)) +
  geom_vline(xintercept = alpha0 / (alpha0 + beta0), color = "red", lty = 2) +
  xlab("Estimated batting average (w/ 95% interval)") +
  ylab("Player") 