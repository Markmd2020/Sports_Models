######Empirical  Bayes Examples For Basebal####

#Load libraries
library(dplyr)

#Simulation of hits examples
num_trials <- 10e6
simulations <- data_frame(
  true_average = rbeta(num_trials, 81, 219),
  hits = rbinom(num_trials, 300, true_average)
)