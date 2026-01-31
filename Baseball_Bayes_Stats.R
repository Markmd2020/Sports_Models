######Empirical  Bayes Examples For Baseball####

#Load libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(Lahman)
library(stats4)
library(VGAM)
library(gamlss)
library(broom)
library(splines)

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

career <- Batting %>%
  filter(AB > 0) %>%
  anti_join(Pitching, by = "playerID") %>%
  group_by(playerID) %>%
  summarize(H = sum(H), AB = sum(AB)) %>%
  mutate(average = H / AB)

# values estimated by maximum likelihood in Chapter 3
alpha0 <- 101.4
beta0 <- 287.3
career_eb <- career %>%
  mutate(eb_estimate = (H + alpha0) / (AB + alpha0 + beta0),
         alpha1 = H + alpha0,
         beta1 = AB- H + beta0)

career_eb %>%
  filter(playerID == "Hank Aaron")

#Calculating Posterior Error Probability
pbeta(.3, 3850, 8818)

career_eb <- career_eb %>%
  mutate(PEP = pbeta(.3, alpha1, beta1))

#Identify top 100 players
top_players <- career_eb %>%
  arrange(PEP) %>%
  head(100)

sum(top_players$PEP)

#Calculate FDR
mean(top_players$PEP)

#Repeat the analysis for top 50 players
sorted_PEP <- career_eb %>%
  arrange(PEP)
mean(head(sorted_PEP$PEP, 50))

#Calculate Q-Values
#The term q-value was first defined by John Storey (Storey, 2002)
#as an analogue to the p-value for controlling FDRs in multiple testing

career_eb <- career_eb %>%
  arrange(PEP) %>%
  mutate(qvalue = cummean(PEP))

hall_of_fame <- career_eb %>%
  filter(qvalue < .05)

strict_hall_of_fame <- career_eb %>%
  filter(qvalue < .01)

#Bayesian A/B Testing
# For each player, update the beta prior based on the evidence
# to get posterior parameters alpha1 and beta1
career_eb <- career %>%
  mutate(eb_estimate = (H + alpha0) / (AB + alpha0 + beta0)) %>%
  mutate(alpha1 = H + alpha0,
         beta1 = AB- H + beta0) %>%
  arrange(desc(eb_estimate))

# while we're at it, save them as separate objects too for later:
hamil <- career_eb %>% filter(playerID == "hamilbi01")
gehri <- career_eb %>% filter(playerID == "gehrilo01")

two_players <- bind_rows(hamil, gehri)

#Simulation of posterior draws
hamil_simulation <- rbeta(1e6, hamil$alpha1, hamil$beta1)
gehri_simulation <- rbeta(1e6, gehri$alpha1, gehri$beta1)
sim <- mean(hamil_simulation > gehri_simulation)

#Using integration  to run A/B tests
d <- .00002
limits <- seq(.29, .33, d)

sum(outer(limits, limits, function(x, y) {
  (x > y) *
    dbeta(x, hamil$alpha1, hamil$beta1) *
    dbeta(y, gehri$alpha1, gehri$beta1) *
    d ^ 2
}))

#Applying closed form solution
h <- function(alpha_a, beta_a,
              alpha_b, beta_b) {
  j <- seq.int(0, round(alpha_b)- 1)
  log_vals <- (lbeta(alpha_a + j, beta_a + beta_b)- log(beta_b + j)-
               lbeta(1 + j, beta_b)- lbeta(alpha_a, beta_a))
  1- sum(exp(log_vals))
}

h <- function(alpha_a, beta_a,
alpha_b, beta_b) {
j <- seq.int(0, round(alpha_b)- 1)
log_vals <- (lbeta(alpha_a + j, beta_a + beta_b)- log(beta_b + j) -
lbeta(1 + j, beta_b)- lbeta(alpha_a, beta_a))
1- sum(exp(log_vals))
}

h(hamil$alpha1, hamil$beta1,gehri$alpha1, gehri$beta1)

#Apply close form approximation
h_approx <- function(alpha_a, beta_a, alpha_b, beta_b) {
  u1 <- alpha_a / (alpha_a + beta_a)
  u2 <- alpha_b / (alpha_b + beta_b)
  var1 <- (alpha_a * beta_a) /
    ((alpha_a + beta_a) ^ 2 * (alpha_a + beta_a + 1))
  var2 <- (alpha_b * beta_b) /
    ((alpha_b + beta_b) ^ 2 * (alpha_b + beta_b + 1))
  pnorm(0, u2- u1, sqrt(var1 + var2))
}

h_approx(hamil$alpha1, hamil$beta1,gehri$alpha1, gehri$beta1)

#Run Chi-Square test
prop.test(two_players$H, two_players$AB)

#Calculate credible intervals
credible_interval_approx <- function(a, b, c, d) {
  u1 <- a / (a + b)
  u2 <- c / (c + d)
  var1 <- a * b / ((a + b) ^ 2 * (a + b + 1))
  var2 <- c * d / ((c + d) ^ 2 * (c + d + 1))
  mu_diff <- u2- u1
  sd_diff <- sqrt(var1 + var2)
  data_frame(posterior = pnorm(0, mu_diff, sd_diff),
             estimate = mu_diff,
             conf.low = qnorm(.025, mu_diff, sd_diff),
             conf.high = qnorm(.975, mu_diff, sd_diff))
}

credible_interval_approx(hamil$alpha1, hamil$beta1,gehri$alpha1, gehri$beta1)

#Beta Binomial Regression

# values estimated by maximum likelihood in Chapter 3
alpha0 <- 101.4
beta0 <- 287.3
prior_mu <- alpha0 / (alpha0 + beta0)
# for each player, update the beta prior based on the evidence
# to get posterior parameters alpha1 and beta1
career_eb <- career %>%
  mutate(eb_estimate = (H + alpha0) / (AB + alpha0 + beta0)) %>%
  mutate(alpha1 = H + alpha0,
         beta1 = AB- H + beta0) %>%
  arrange(desc(eb_estimate))

#Fit model across all players
fit <- gamlss(cbind(H, AB- H) ~ log(AB),
              data = career_eb,
              family = BB(mu.link = "identity"))

td <- tidy(fit)

mu <- fitted(fit, parameter = "mu")
sigma <- fitted(fit, parameter = "sigma")
head(mu)
head(sigma)

career_eb_wAB <- career_eb %>%
  dplyr::select(playerID, H, AB, original_eb = eb_estimate) %>%
  mutate(mu = mu,
         alpha0 = mu / sigma,
         beta0 = (1- mu) / sigma,
         alpha1 = alpha0 + H,
         beta1 = beta0 + AB- H,
         new_eb = alpha1 / (alpha1 + beta1))

#Empirical Bayesian Hierarchical Modelling

#Right and left hand batters
career %>%
  count(bats)

# relevel to set right-handed batters as the baseline
career2 <- career %>%
  filter(!is.na(bats)) %>%
  mutate(bats = relevel(bats, "R"))

fit2 <- gamlss(cbind(H, AB- H) ~ log(AB) + bats,
               data = career2,
               family = BB(mu.link = "identity"))

tidy(fit2)

#Add splines to allow the model vary over time
fit3 <- gamlss(cbind(H, AB- H) ~ 0 + ns(year, df = 5)
               + bats + log(AB),
               data = career2,
               family = BB(mu.link = "identity"))

#Adding an interaction term
fit4 <- gamlss(cbind(H, AB- H) ~ 0 + ns(year, 5) * bats
               + log(AB),
               data = career2,
               family = BB(mu.link = "identity"))

#Calculate posterior distributions
players <- crossing(year = c(1915, 1965, 2015),
                    bats = c("L", "R"),
                    H = 30,
                    AB = 100)

players_posterior <- players %>%
  mutate(mu = predict(fit4, what = "mu", newdata = players),
         sigma = predict(fit4, what = "sigma",
                         newdata = players, type = "response"),
         alpha0 = mu / sigma,
         beta0 = (1- mu) / sigma,
         alpha1 = alpha0 + H,
         beta1 = beta0 + AB- H)

#Mixture Models and Expectation-Maximisation
