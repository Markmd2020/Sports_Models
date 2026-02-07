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
library(purrr)
library(ebbr)

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

#Generate sythentic dataset

# number of observations
n <- 300

career <- data.frame(
  year = sample(1990:2020, n, replace = TRUE),
  AB   = sample(200:700, n, replace = TRUE),          # at-bats
  bats = factor(sample(c("L", "R", "S"), n, replace = TRUE))
)

# Generate synthetic hit probabilities using a smooth function of year + bats
# This ensures the BB model has realistic structure
p_base <- 0.20 + 0.05 * sin((career$year - 1990) / 30 * pi)
p_bats <- ifelse(career$bats == "L", 0.02,
                 ifelse(career$bats == "R", -0.01, 0))

p <- pmin(pmax(p_base + p_bats + rnorm(n, 0, 0.02), 0.05), 0.40)

# Hits (H) must be <= AB
career$H <- rbinom(n, size = career$AB, prob = p)

# Inspect
head(career)

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

summary(players_posterior)

#Mixture Models and Expectation-Maximisation
# We'll fit the clusters only with players that have had at least 20 at-bats
starting_data <- career %>%
  filter(AB >= 20) %>%
  mutate(cluster = factor(sample(c("A", "B"), n(), replace = TRUE)))

fit_bb_mle <- function(x, n) {
  # dbetabinom.ab is the likelihood function for a beta-binomial
  # using n, alpha and beta as parameters
  ll <- function(alpha, beta) {-sum(dbetabinom.ab(x, n, alpha, beta, log = TRUE))
  }
  m <- stats4::mle(ll, start = list(alpha = 3, beta = 10),
                   method = "L-BFGS-B", lower = c(0.001, .001))
  ab <- stats4::coef(m)
  data_frame(alpha = ab[1], beta = ab[2])
}

fit_bb_mle(starting_data$H, starting_data$AB)

#Build Mixture Model
fits <- starting_data %>%
  group_by(cluster) %>%
  do(fit_bb_mle(.$H, .$AB)) %>%
  ungroup()

crosses <- starting_data %>%
  select(-cluster) %>%
  crossing(fits) %>%
  mutate(likelihood = VGAM::dbetabinom.ab(H,AB,alpha,beta))

#Assign players to clusters
assignments<-starting_data%>%
  select(-cluster)%>%
  crossing(fits) %>%
  mutate(likelihood = VGAM::dbetabinom.ab(H,AB,alpha,beta))%>%
  group_by(playerID) %>%
  top_n(1,likelihood) %>%
  ungroup()

#Iterate the EM algorithm
iterate_em <- function(state, ...) {
  # maximization
  fits <- state$assignments %>%
    group_by(cluster) %>%
    do(fit_bb_mle(.$H, .$AB)) %>%
    ungroup()
  # expectation
  assignments <- state$assignments %>%
    select(playerID:average) %>%
    crossing(fits) %>%
  mutate(likelihood = VGAM::dbetabinom.ab(H, AB, alpha, beta)) %>%
    group_by(playerID) %>%
    top_n(1, likelihood) %>%
    ungroup()
  list(assignments = assignments, fits = fits)
}
 
init <- list(assignments = starting_data)
iterations <- accumulate(1:5, iterate_em, .init = init)

fit_iterations %>%
  crossing(x = seq(.001, .4, .001)) %>%
  mutate(density = dbeta(x, alpha, beta)) %>%
  ggplot(aes(x, density, color = iteration, group = iteration)) +
  geom_line() +
  facet_wrap(~ cluster)

final_parameters <- last(iterations)$fits
career_likelihoods <- career %>%
  filter(AB > 20) %>%
  crossing(final_parameters) %>%
  mutate(likelihood = VGAM::dbetabinom.ab(H, AB, alpha, beta)) %>%
  group_by(playerID) %>%
  mutate(posterior = likelihood / sum(likelihood))

career_assignments <- career_likelihoods %>%
  top_n(1, posterior) %>%
  ungroup()

eb_shrinkage <- career_likelihoods %>%
mutate(shrunken_average = (H + alpha) / (AB + alpha + beta)) %>%
  group_by(playerID) %>%
  summarize(shrunken_average = sum(posterior * shrunken_average))

#Multinomial and Direchlet

#include the "bats"(handedness)and"year"columnforlater

pitchers <- Pitching %>%
  group_by(playerID) %>%
  summarize(gamesPitched = sum(G)) %>%
  filter(gamesPitched> 3)

hit_types <-Batting%>%
  filter(AB > 0)%>%
  anti_join(pitchers, by= "playerID")%>%
  rename(Double = X2B,Triple= X3B)%>%
  group_by(playerID) %>%
  summarize_each(funs(sum(.,na.rm=TRUE)),AB,H,Double,Triple,HR)%>%
  inner_join(player_names, by="playerID")%>%
  transmute(playerID,name,AB,H,
            Single = H-Double-Triple-HR,
            Double,Triple,HR,
            NonHit = AB-H)

#Multinomial Distribution Example
rmultinom(3, 100, c(.2, .2, .2, .2, .2))
rmultinom(3, 100, c(.05, .05, .2, .2, .5))

#Using a Direchlet prior
VGAM::rdiric(3, c(1, 1, 1, 1, 1))

#Fitting a Direchlet Multinomial Model
hit_500 <- hit_types %>%
  filter(AB >= 500)
hit_matrix <- hit_500 %>%
  select(Single, Double, Triple, HR, NonHit) %>%
  as.matrix()
dm_fit <- DirichletMultinomial::dmn(hit_matrix, 1)

#Extract model parameters
tidy.DMN <- function(x, ...) {
  ret <- as.data.frame(x@fit)
  tbl_df(fix_data_frame(ret, c("conf.low", "estimate", "conf.high")))
}
dm_params <- tidy(dm_fit)

#Empirical Bayes shrinkage of slugging percentage
par_total <- sum(dm_params$estimate)
par <- dm_params %>%
  select(term, estimate) %>%
  spread(term, estimate)

w<-c(1:4,0)
slugging_mean<-sum(w* dm_params$estimate)/ sum(dm_params$estimate)

hit_types <-hit_types%>%
  mutate(slugging= (Single+ 2*Double+ 3*Triple+ 4*HR)/ AB)

hit_types_eb<-hit_types%>%
  mutate(slugging_eb= ((Single+ par$Single) +
                         (Double+par$Double)* 2+
                         (Triple+par$Triple)* 3+
                         (HR+par$HR)* 4)/
           (AB + par_total))

#Hierarchical Modelling with Built In Function
career %>%
  filter(AB >= 10) %>%
  ggplot(aes(AB, H / AB)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_log10()

#Perform empirical Bayes shrinkage 
eb_career_ab <- career %>%
  add_ebb_estimate(H, AB, method = "gamlss",
                   mu_predictors = ~ log10(AB))

#Simulation Examples
prior <- career %>%
  ebb_fit_prior(H, AB)

alpha0 <- tidy(prior)$alpha
beta0 <- tidy(prior)$beta

rbeta(10, alpha0, beta0)

career_sim <- career %>%
  mutate(p = rbeta(n(), alpha0, beta0),
         H = rbinom(n(), AB, p))

career_sim_eb <- career_sim %>%
  add_ebb_estimate(H, AB)

career_sim_gathered <- career_sim_eb %>%
  rename(Shrunken = .fitted, Raw = .raw) %>%
  gather(type, estimate, Shrunken, Raw)

career_sim_gathered %>%
  filter(AB >= 10) %>%
  ggplot(aes(p, estimate, color = AB)) +
  geom_point() +
  geom_abline(color = "red") +
  geom_smooth(method = "lm", color = "white", lty = 2, se = FALSE) +
  scale_color_continuous(trans = "log",
                         breaks = c(10, 100, 1000, 10000)) +
  facet_wrap(~ type) +
  labs(x = "True batting average (p)",
       y = "Raw or shrunken batting average")

#Model Evaluation 
career_sim_gathered %>%
  group_by(type) %>%
  summarize(mse = mean((estimate- p) ^ 2))

metric_by_bin <- career_sim_gathered %>%
  group_by(type, AB = 10 ^ (round(log10(AB)))) %>%
  summarize(mse = mean((estimate- p) ^ 2))

ggplot(metric_by_bin, aes(AB, mse, color = type)) +
  geom_line() +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Number of at-bats (AB)",
       y = "Mean-squared-error within this bin")

#Check how often credible intervals capture the actual batting average
career_sim_eb %>%
  sample_n(20) %>%
  mutate(playerID = reorder(playerID, .fitted)) %>%
  ggplot(aes(.fitted, playerID)) +
  geom_point() +
  geom_point(aes(x = p), color = "red") +
  geom_errorbarh(aes(xmin = .low, xmax = .high)) +
  theme(axis.text.y = element_blank()) +
  labs(x = "Estimated batting average (w/ 95% credible interval)",
       y = "Player")

#Evaluate coverage 
career_sim_eb %>%
  summarize(coverage = mean(.low <= p & p <= .high))

# fit the prior once
sim_prior <- ebb_fit_prior(career_sim, H, AB)

# find the coverage probability for each level
estimate_by_cred_level <- data_frame(level = seq(.5, .98, .02)) %>%
  unnest(map(level, ~ augment(sim_prior, career_sim, cred_level = .)))

estimate_by_cred_level %>%
  group_by(level) %>%
  mutate(cover = .low <= p & p <= .high) %>%
  summarize(coverage = mean(cover)) %>%
  ggplot(aes(level, coverage)) +
  geom_line() +
  geom_abline(color = "red", lty = 2) +
  labs(x = "Level of credible interval",
       y = "Probability credible interval contains the true value")

#FDR Control
pt <- career_sim_eb %>%
  add_ebb_prop_test(.3, sort = TRUE) 

# Control for FDR of 10%
hall_of_fame <- pt %>%
  filter(.qvalue <= .1)
nrow(hall_of_fame)
mean(hall_of_fame$p < .3)

pt %>%
  mutate(true_fdr = cummean(p < .3)) %>%
  ggplot(aes(.qvalue, true_fdr)) +
  geom_line() +
  geom_abline(color = "red") +
  labs(x = "q-value threshold",
       y = "True FDR at this q-value threshold")

#Beta Binomial Regression
bb_reg <- career %>%
  ebb_fit_prior(H, AB, method = "gamlss", mu_predictors = ~ log10(AB))
tidy(bb_reg)

career_sim_ab <- augment(bb_reg, career) %>%
  select(playerID, AB,
         true_alpha0 = .alpha0,
         true_beta0 = .beta0) %>%
  mutate(p = rbeta(n(), true_alpha0, true_beta0),
         H = rbinom(n(), AB, p))

#Assess model performance
career_ab_prior <- career_sim_ab %>%
  ebb_fit_prior(H, AB, method = "gamlss", mu_predictors = ~ log10(AB))
tidy(bb_reg)
tidy(career_ab_prior)

career_flat_prior <- career_sim_ab %>%
  ebb_fit_prior(H, AB)

data_frame(method = c("Flat prior", "Prior depending on AB"),
           model = list(career_flat_prior, career_ab_prior)) %>%
  unnest(map(model, augment, data = career_sim_ab)) %>%
  ggplot(aes(p, .fitted, color = AB)) +
  geom_point() +
  scale_color_continuous(trans = "log", breaks = c(1, 10, 100, 1000)) +
  geom_abline(color = "red") +
  facet_wrap(~ method) +
  labs(x = "True batting average (p)",
       y = "Shrunken batting average estimate")

#Running multiple simulations
prior <- ebb_fit_prior(career, H, AB)
alpha0 <- tidy(prior)$alpha
beta0 <- tidy(prior)$beta

sim_replications <- career %>%
  crossing(replication = 1:50) %>%
  mutate(p = rbeta(n(), alpha0, beta0),
         H = rbinom(n(), AB, p))

#Fit multiple priors
sim_replication_models <- sim_replications %>%
  nest(-replication) %>%
  mutate(prior = map(data, ~ ebb_fit_prior(., H, AB)))

sim_replication_priors <- sim_replication_models %>%
  unnest(map(prior, tidy), .drop = TRUE)


sim_replication_au <- sim_replication_models %>%
  unnest(map2(prior, data, augment))

sim_replication_mse <- sim_replication_au %>%
  rename(Raw = .raw, Shrunken = .fitted) %>%
  gather(type, estimate, Raw, Shrunken) %>%
  group_by(type, replication) %>%
  summarize(mse = mean((estimate- p) ^ 2))

sim_replication_intervals <- sim_replication_models %>%
  crossing(cred_level = c(seq(.5, .9, .05), .95)) %>%
  unnest(pmap(list(prior, data, cred_level = cred_level), augment)) %>%
  select(replication, cred_level, p, .low, .high)

#Simulating varying observations
varying_size_sim <- career %>%
  select(-H) %>%
  crossing(size = c(30, 100, 300, 1000, 3000, 10000),
           replication = 1:50) %>%
  group_by(size, replication) %>%
  sample_frac(1, replace = TRUE) %>%
  filter(row_number() <= size) %>%
  ungroup()

varying_size_priors <- varying_size_sim %>%
  mutate(p = rbeta(n(), alpha0, beta0),
         H = rbinom(n(), AB, p)) %>%
  nest(-size,-replication) %>%
  mutate(prior = map(data, ~ ebb_fit_prior(., H, AB)))