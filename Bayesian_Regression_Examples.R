#Bayesian Regression Example

#Load libraries
library(brms)
library(tidyverse)
library(tidybayes)
library(bayesplot)
library(posterior)
library(loo)
library(rstanarm)

#Create a simple dataset
set.seed(135)
n <- 120
advertising_spend <- rnorm(n, mean = 15, sd = 4)
sales <- 20 + 3.5 * advertising_spend + rnorm(n, mean = 0, sd = 8)

df <- data.frame(
  advertising_spend = advertising_spend,
  sales = sales
)

head(df)

#Data Exploration
summary(df)

ggplot(df, aes(x = advertising_spend, y = sales)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal()

lm_model <- lm(sales ~ advertising_spend, data = df)
summary(lm_model)

#Choosing priors
priors <- c(
  prior(normal(0, 20), class = "Intercept"),
  prior(normal(0, 10), class = "b"),
  prior(student_t(3, 0, 10), class = "sigma")
)


#Fit Bayesian Linear Regression Model
bayes_model <- brm(
  formula = sales ~ advertising_spend,
  data = df,
  prior = priors,
  family = gaussian(),
  chains = 4,
  iter = 4000,
  warmup = 2000,
  seed = 135
)

#Model Summary
summary(bayes_model)

#Extracting posterior draws
draws <- as_draws_df(bayes_model)
head(draws)

mean(draws$b_advertising_spend > 0)

#Plot posterior distributions
plot(bayes_model)

#Visualise intervals
mcmc_areas(
  as.array(bayes_model),
  pars = c("b_Intercept", "b_advertising_spend", "sigma")
)

#Checking convergence
mcmc_trace(
  as.array(bayes_model),
  pars = c("b_Intercept", "b_advertising_spend", "sigma")
)

#Posterior Predictive Checks
pp_check(bayes_model)

pp_check(bayes_model, type = "dens_overlay")
pp_check(bayes_model, type = "hist")
pp_check(bayes_model, type = "scatter_avg")

tidy_draws <- bayes_model %>%
  spread_draws(b_Intercept, b_advertising_spend, sigma)

head(tidy_draws)

tidy_draws %>%
  ggplot(aes(x = b_advertising_spend)) +
  geom_density(fill = "steelblue", alpha = 0.4) +
  theme_minimal()

#Generate predictions for new data
new_customers <- data.frame(
  advertising_spend = c(10, 15, 20, 25)
)

predict(bayes_model, newdata = new_customers)

fitted(bayes_model, newdata = new_customers)

#Display conditional effects
conditional_effects(bayes_model)

#Alternative approach with STAN model
rstanarm_model <- stan_glm(
  sales ~ advertising_spend,
  data = df,
  family = gaussian(),
  chains = 4,
  iter = 4000,
  seed = 135
)

print(rstanarm_model)

#Leave One Out CrossValidation
loo_result <- loo(bayes_model)
print(loo_result)