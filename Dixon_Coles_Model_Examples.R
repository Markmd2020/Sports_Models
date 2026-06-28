#######Dixon Cole Models#######

##Load libraries
library(goalmodel)

# Expected goals by the the two oppsing teams.
expg1 <- 1.1
expg2 <- 1.9

# The upper limit of how many goals to compute probabilities for.
maxgoal <- 10

# The "S" matrix, which can also be computed by predict_goals().
# Assuming the independent Poisson model.
probmat <- dpois(0:maxgoal, expg1) %*% t(dpois(0:maxgoal, expg2))

# Compute the BTTS probability using the formula from the paper.
prob_btts <- 1 - (sum(probmat[2:nrow(probmat),1]) +
                    sum(probmat[1,2:ncol(probmat)]) + probmat[1,1]) 

#2nd Approach
prob_btts_2 <- (1 - dpois(0, lambda = expg1)) * (1 - dpois(0, lambda = expg2))

#Dixon Cole Adjustment

# Dixon-Coles adjustment parameter.
rho <- 0.13

# 1-1 probability for the independent Poisson model.
p11 <- dpois(1, lambda = expg1) * dpois(1, lambda = expg2) 

# Add DC adjusted 1-1 probability, subtract unadjusted 1-1 probability.
dc_correction <- (p11 * (1-rho)) - p11

# Apply the corrections
prob_btts_dc  <- prob_btts_2 + dc_correction
 
