# A comparison of Alonzo, Beta, and Lehmann ROC regression methods 
# Covariate adjusted data
# Normal, Extreme Value, and Weibull

########################################
##    Load the necessary packages     ##
########################################
library(tidyverse)
library(survival)
library(survminer)
library(PRROC)
library(quantreg)
library(betareg)

########################################
##      THE SIMULATION IS BELOW       ##
########################################
# What is the number of simulations?
nsim <- 1000

# Results will be stored in these matrices
matrix_list_normal <- list(auc_alonzo = matrix(NA, nsim, 4),
                    auc_beta = matrix(NA, nsim, 4),
                    auc_lehmann = matrix(NA, nsim, 4))

matrix_list_exval <- list(auc_alonzo = matrix(NA, nsim, 4),
                           auc_beta = matrix(NA, nsim, 4),
                           auc_lehmann = matrix(NA, nsim, 4))

matrix_list_weib <- list(auc_alonzo = matrix(NA, nsim, 4),
                          auc_beta = matrix(NA, nsim, 4),
                          auc_lehmann = matrix(NA, nsim, 4))

# Maximum Youden Index
youden_list_normal <- list(youden_alonzo = matrix(NA, nsim, 4), 
                          youden_beta = matrix(NA, nsim, 4),
                          youden_lehmann = matrix(NA, nsim, 4))

youden_list_exval <- list(youden_alonzo = matrix(NA, nsim, 4), 
                          youden_beta = matrix(NA, nsim, 4),
                          youden_lehmann = matrix(NA, nsim, 4))

youden_list_weib <- list(youden_alonzo = matrix(NA, nsim, 4), 
                          youden_beta = matrix(NA, nsim, 4),
                          youden_lehmann = matrix(NA, nsim, 4))

# What I think are the true AUC values
# These are the x values we are considering
x2 <- c(0.2, 0.4, 0.6, 0.8)
# For normal case
true_normal_auc <- pnorm((0.5 + 4*x2)/1.5)
# For extreme value case
aucfunction <- function(x){
  roc_true_exval <- function(fpr){
    1 - exp(-exp(log(-log(1-fpr)) + ((2 - 1.5) + (4 - 0)*x)/1.5))
  }
  
  # Evaluate that function across the FPR's
  integrate(f = roc_true_exval, lower = 0.001, upper = 0.999)$value
}

# Integrate all four simultaneously!
true_exval_auc <- sapply(x2, aucfunction)

# For the Weibull case
truetheta <- 3/7 - 0.1*x2
true_weib_auc <- 1/(1+truetheta)

for(counter in 865:nsim){
  print(counter)
#### Normal Case ####
# Generate data (no repeated measures -- yet)
  data_frame(x = runif(100)) %>% 
    mutate(dis = sapply(x, function(x) rnorm(1, 2 + 4*x, (1.5)))) %>% 
    mutate(ref = sapply(x, function(x) rnorm(1, 1.5 + 0*x, (1.5))))  %>% 
    gather(group, obs, -x) %>% 
    mutate(ind = as.factor(case_when(
      (group == "ref") ~ 0,
      (group == "dis") ~ 1
    ))) -> df1

## ALONZO METHOD ##
# Find the quantiles
quant_matrix <- quantreg::rq(obs ~ x, data = df1[df1$group == "ref",], 
                             tau = seq(0.01, 0.99, by = 0.01))$fitted.values
quants <- colMeans(quant_matrix)

# Rearrange the data into the necessary format
df1 %>%
  slice(rep(1:n(), each = 99)) %>% 
  mutate(fpr = rep(seq(0.01, 0.99, by = 0.01), 200)) %>% 
  mutate(q = rep(quants, 200)) %>% 
  mutate(phi_inv = qnorm(fpr)) %>% 
  mutate(u = obs >= q) %>% 
  filter(group == "dis") -> df2

# Find the coefficients for the Alonzo model 
probit_coeff <- glm(u ~ phi_inv + x, data = df2, 
                    family = binomial(link = "probit"))$coefficients

# This is just so R shows us the results
probit_coeff

# Calculate some ROC values
data_frame(x = rep(c(0.2, 0.4, 0.6, 0.8),200), 
           fpr = rep(seq(0.001, 0.999, by = 0.005), 4)) %>% 
  mutate(phi_inv = qnorm(fpr)) %>% 
  mutate(ROC = pnorm(1*probit_coeff[1] + abs(probit_coeff[2])*phi_inv +
                       probit_coeff[3]*x)) %>% 
  mutate(Youden = abs(ROC - fpr)) -> df3  

# Youden Index Section
df3 %>% 
  filter(x == 0.20) %>% 
  select(Youden) %>% 
  max() -> youden_list_normal$youden_alonzo[counter,1]
df3 %>% 
  filter(x == 0.40) %>% 
  select(Youden) %>% 
  max() -> youden_list_normal$youden_alonzo[counter,2]
df3 %>% 
  filter(x == 0.60) %>% 
  select(Youden) %>% 
  max() -> youden_list_normal$youden_alonzo[counter,3]
df3 %>% 
  filter(x == 0.80) %>% 
  select(Youden) %>% 
  max() -> youden_list_normal$youden_alonzo[counter,4]
  
# ROC finding function
big_function <- function(x){
  roc_f <- function(fpr){
    pnorm(1*probit_coeff[1] + probit_coeff[2]*qnorm(fpr) +
            probit_coeff[3]*x)
  }

  
# Evaluate that function across the FPR's
integrate(f = roc_f, lower = 0.001, upper = 0.999)$value
}

# Integrate all four simultaneously!
Alonzo_AUC <- sapply(c(0.2, 0.4, 0.6, 0.8), big_function)

# Save the results
matrix_list_normal$auc_alonzo[counter,] <- Alonzo_AUC

## BETA REGRESSION METHOD ##
# Calculate the placement values
pv <- sapply(df1[df1$group == "dis",]$obs,function(y) (1 - mean(y >= quants)))

df1 %>% 
  filter(group == "dis") %>% 
  dplyr::select(c(x, obs)) %>% 
  mutate(pv = pv) -> df4

# The dependent variable must be in (0,1)
# Change 0 to 0.0001 and 1 to 0.9999, if necessary
if(any(df4$pv == 0)){
df4[df4$pv == 0,]$pv <- 0.0001
}
if(any(df4$pv == 1)){
df4[df4$pv == 1,]$pv <- 0.9999
}

# Beta regression with the placement values
beta_fit <- betareg(pv ~ x, data = df4)

# The coefficients returned by beta regression
beta0 <- beta_fit$coefficients$mean[1] # Intercept
beta1 <- beta_fit$coefficients$mean[2] # Slope
phi <- beta_fit$coefficients$precision # Precision 

# Calculate some ROC values
data_frame(x = rep(c(0.2, 0.4, 0.6, 0.8),200), 
           fpr = rep(seq(0.001, 0.999, by = 0.005), 4)) %>% 
  mutate(m = 1/(1 + exp(-beta0 - beta1*x))) %>% 
  mutate(ROC = pbeta(fpr, shape1 = m*phi, shape2 = phi*(1-m))) %>% 
  mutate(Youden = abs(ROC - fpr)) -> df5 

# Youden Index Section
df5 %>% 
  filter(x == 0.20) %>% 
  select(Youden) %>% 
  max() -> youden_list_normal$youden_beta[counter,1]
df5 %>% 
  filter(x == 0.40) %>% 
  select(Youden) %>% 
  max() -> youden_list_normal$youden_beta[counter,2]
df5 %>% 
  filter(x == 0.60) %>% 
  select(Youden) %>% 
  max() -> youden_list_normal$youden_beta[counter,3]
df5 %>% 
  filter(x == 0.80) %>% 
  select(Youden) %>% 
  max() -> youden_list_normal$youden_beta[counter,4]


# Finding the ROC
big_function2 <- function(x){
  roc_beta <- function(fpr){
    mu <- 1/(1 + exp(-beta0 - beta1*x))
    a <- mu * phi 
    b <- phi * (1 - mu)
    
    pbeta(fpr, shape1 = a, shape2 = b)
  }
  
  # Evaluate that function across the FPR's
  integrate(f = roc_beta, lower = 0.001, upper = 0.999)$value
}

# Integrate all four simultaneously!
beta_AUC <- sapply(c(0.2, 0.4, 0.6, 0.8), big_function2)

# Save the results
matrix_list_normal$auc_beta[counter,] <- beta_AUC

## LEHMANN METHOD ##
# Find the survival curves
S <- Surv(df1$obs)
fit <- survfit(S ~ df1$group)

df_surv <- data_frame(time = fit$time, 
                      survival = fit$surv,
                      group = rev(df1$group))

# Check the proportional hazards assumption
coxx <- coxph(formula = Surv(obs) ~ as.factor(ind) + x 
              + as.factor(ind)*x,
              data = df1)
coxx_px <- cox.zph(coxx)

# Cox proportional hazards model to estimate beta 
# Which can be used to estimate theta
beta1 <- unname(coxx$coefficients[1])
beta3 <- unname(coxx$coefficients[3])

# Calculate some ROC values
data_frame(x = rep(c(0.2, 0.4, 0.6, 0.8),200), 
           fpr = rep(seq(0.001, 0.999, by = 0.005), 4))  %>% 
  mutate(theta = exp(beta1 + beta3*x)) %>% 
  mutate(ROC = fpr^theta) %>% 
  mutate(Youden = abs(ROC - fpr)) -> df6

# Youden Index Section
df6 %>% 
  filter(x == 0.20) %>% 
  select(Youden) %>% 
  max() -> youden_list_normal$youden_lehmann[counter,1]
df6 %>% 
  filter(x == 0.40) %>% 
  select(Youden) %>% 
  max() -> youden_list_normal$youden_lehmann[counter,2]
df6 %>% 
  filter(x == 0.60) %>% 
  select(Youden) %>% 
  max() -> youden_list_normal$youden_lehmann[counter,3]
df6 %>% 
  filter(x == 0.80) %>% 
  select(Youden) %>% 
  max() -> youden_list_normal$youden_lehmann[counter,4]


# Finding the ROC
big_function3 <- function(x){
  roc_lehmann <- function(fpr){
    theta <- exp(beta1 + beta3*x)
    
    fpr^theta
  }
  
  # Evaluate that function across the FPR's
  integrate(f = roc_lehmann, lower = 0.001, upper = 0.999)$value
}

# Integrate all four simultaneously!
Lehmann_AUC <- sapply(c(0.2, 0.4, 0.6, 0.8), big_function3)

# Save the results
matrix_list_normal$auc_lehmann[counter,] <- Lehmann_AUC
###########################################################################
#### EXTREME VALUE CASE ####
data_frame(x = runif(100)) %>% 
  mutate(dis = sapply(x, function(x) 2 + 4*x - 1.5*log(rexp(1)))) %>% 
  mutate(ref = sapply(x, function(x) 1.5 + 0*x - 1.5*log(rexp(1))))  %>% 
  gather(group, obs, -x) %>% 
  mutate(ind = as.factor(case_when(
    (group == "ref") ~ 0,
    (group == "dis") ~ 1
  ))) -> df1

## ALONZO METHOD ##
# Find the quantiles
quant_matrix <- quantreg::rq(obs ~ x, data = df1[df1$group == "ref",], 
                             tau = seq(0.01, 0.99, by = 0.01))$fitted.values
quants <- colMeans(quant_matrix)

# Rearrange the data into the necessary format
df1 %>%
  slice(rep(1:n(), each = 99)) %>% 
  mutate(fpr = rep(seq(0.01, 0.99, by = 0.01), 200)) %>% 
  mutate(q = rep(quants, 200)) %>% 
  mutate(phi_inv = qnorm(fpr)) %>% 
  mutate(u = obs >= q) %>% 
  filter(group == "dis") -> df2

# Find the coefficients for the Alonzo model 
probit_coeff <- glm(u ~ phi_inv + x, data = df2, 
                    family = binomial(link = "probit"))$coefficients

# This is just so R shows us the results
probit_coeff

# Calculate some ROC values
data_frame(x = rep(c(0.2, 0.4, 0.6, 0.8),200), 
           fpr = rep(seq(0.001, 0.999, by = 0.005), 4)) %>% 
  mutate(phi_inv = qnorm(fpr)) %>% 
  mutate(ROC = pnorm(1*probit_coeff[1] + abs(probit_coeff[2])*phi_inv +
                       probit_coeff[3]*x)) %>% 
  mutate(Youden = abs(ROC - fpr)) -> df3

# Youden Index Section
df3 %>% 
  filter(x == 0.20) %>% 
  select(Youden) %>% 
  max() -> youden_list_exval$youden_alonzo[counter,1]
df3 %>% 
  filter(x == 0.40) %>% 
  select(Youden) %>% 
  max() -> youden_list_exval$youden_alonzo[counter,2]
df3 %>% 
  filter(x == 0.60) %>% 
  select(Youden) %>% 
  max() -> youden_list_exval$youden_alonzo[counter,3]
df3 %>% 
  filter(x == 0.80) %>% 
  select(Youden) %>% 
  max() -> youden_list_exval$youden_alonzo[counter,4]


# ROC finding function
big_function <- function(x){
  roc_f <- function(fpr){
    pnorm(1*probit_coeff[1] + probit_coeff[2]*qnorm(fpr) +
            probit_coeff[3]*x)
  }
  
  # Evaluate that function across the FPR's
  integrate(f = roc_f, lower = 0.001, upper = 0.999)$value
}

# Integrate all four simultaneously!
Alonzo_AUC <- sapply(c(0.2, 0.4, 0.6, 0.8), big_function)

# Save the results
matrix_list_exval$auc_alonzo[counter,] <- Alonzo_AUC

## BETA METHOD ##
# Calculate the placement values
pv <- sapply(df1[df1$group == "dis",]$obs,function(y) (1 - mean(y >= quants)))


df1 %>% 
  filter(group == "dis") %>% 
  dplyr::select(c(x, obs)) %>% 
  mutate(pv = pv) -> df4

# The dependent variable must be in (0,1)
# Change 0 to 0.0001 and 1 to 0.9999, if necessary
if(any(df4$pv == 0)){
  df4[df4$pv == 0,]$pv <- 0.0001
}
if(any(df4$pv == 1)){
  df4[df4$pv == 1,]$pv <- 0.9999
}

# Beta regression with the placement values
beta_fit <- betareg(pv ~ x, data = df4)

# The coefficients returned by beta regression
beta0 <- beta_fit$coefficients$mean[1] # Intercept
beta1 <- beta_fit$coefficients$mean[2] # Slope
phi <- beta_fit$coefficients$precision # Precision 

# Calculate some ROC values
data_frame(x = rep(c(0.2, 0.4, 0.6, 0.8),200), 
           fpr = rep(seq(0.001, 0.999, by = 0.005), 4)) %>% 
  mutate(m = 1/(1 + exp(-beta0 - beta1*x))) %>% 
  mutate(ROC = pbeta(fpr, shape1 = m*phi, shape2 = phi*(1-m))) %>% 
  mutate(Youden = abs(ROC - fpr)) -> df5 

# Youden Index Section
df5 %>% 
  filter(x == 0.20) %>% 
  select(Youden) %>% 
  max() -> youden_list_exval$youden_beta[counter,1]
df5 %>% 
  filter(x == 0.40) %>% 
  select(Youden) %>% 
  max() -> youden_list_exval$youden_beta[counter,2]
df5 %>% 
  filter(x == 0.60) %>% 
  select(Youden) %>% 
  max() -> youden_list_exval$youden_beta[counter,3]
df5 %>% 
  filter(x == 0.80) %>% 
  select(Youden) %>% 
  max() -> youden_list_exval$youden_beta[counter,4]


# Finding the ROC
big_function2 <- function(x){
  roc_beta <- function(fpr){
    mu <- 1/(1 + exp(-beta0 - beta1*x))
    a <- mu * phi 
    b <- phi * (1 - mu)
    
    pbeta(fpr, shape1 = a, shape2 = b)
  }
  
  # Evaluate that function across the FPR's
  integrate(f = roc_beta, lower = 0.001, upper = 0.999)$value
}

# Integrate all four simultaneously!
beta_AUC <- sapply(c(0.2, 0.4, 0.6, 0.8), big_function2)

# Save the results
matrix_list_exval$auc_beta[counter,] <- beta_AUC

## LEHMANN METHOD ##
# Find the survival curves
S <- Surv(df1$obs)
fit <- survfit(S ~ df1$group)

df_surv <- data_frame(time = fit$time, 
                      survival = fit$surv,
                      group = rev(df1$group))

# Check the proportional hazards assumption
coxx <- coxph(formula = Surv(obs) ~ as.factor(ind) + x 
              + as.factor(ind)*x,
              data = df1)
coxx_px <- cox.zph(coxx)

# Cox proportional hazards model to estimate beta 
# Which can be used to estimate theta
beta1 <- unname(coxx$coefficients[1])
beta3 <- unname(coxx$coefficients[3])

# Calculate some ROC values
data_frame(x = rep(c(0.2, 0.4, 0.6, 0.8),200), 
           fpr = rep(seq(0.001, 0.999, by = 0.005), 4))  %>% 
  mutate(theta = exp(beta1 + beta3*x)) %>% 
  mutate(ROC = fpr^theta) %>% 
  mutate(Youden = abs(ROC - fpr)) -> df6

# Youden Index Section
df6 %>% 
  filter(x == 0.20) %>% 
  select(Youden) %>% 
  max() -> youden_list_exval$youden_lehmann[counter,1]
df6 %>% 
  filter(x == 0.40) %>% 
  select(Youden) %>% 
  max() -> youden_list_exval$youden_lehmann[counter,2]
df6 %>% 
  filter(x == 0.60) %>% 
  select(Youden) %>% 
  max() -> youden_list_exval$youden_lehmann[counter,3]
df6 %>% 
  filter(x == 0.80) %>% 
  select(Youden) %>% 
  max() -> youden_list_exval$youden_lehmann[counter,4]

# Finding the ROC
big_function3 <- function(x){
  roc_lehmann <- function(fpr){
    theta <- exp(beta1 + beta3*x)
    
    fpr^theta
  }
  
  # Evaluate that function across the FPR's
  integrate(f = roc_lehmann, lower = 0.001, upper = 0.999)$value
}

# Integrate all four simultaneously!
Lehmann_AUC <- sapply(c(0.2, 0.4, 0.6, 0.8), big_function3)

# Save the results
matrix_list_exval$auc_lehmann[counter,] <- Lehmann_AUC
###########################################################################
#### WEIBULL CASE ####
# Generate data (No repeated measures -- yet)
# Misc
theta <- 3/7
s <- 1
# Scale parameters
lambda_r <- 1
lambda_d <- theta^1/s
# Shape parameters
nu_r <- 1
nu_d <- 1
# Put it together and make data
df1 <- data_frame(x = runif(100),
                  uni = runif(100),
                  dis = (-log(uni)/(lambda_d*exp(0 + 0.6*-x)))^(1/nu_d),
                  ref = (-log(uni)/(lambda_r*exp(0 + 0.5*-x)))^(1/nu_r))

df1 %>% 
  gather(group, obs, -c(x,uni)) %>% 
  select(-uni) %>% 
  mutate(ind = as.factor(case_when(
    (group == "ref") ~ 0,
    (group == "dis") ~ 1
  ))) -> df1

## ALONZO ##
# Find the quantiles
quant_matrix <- quantreg::rq(obs ~ x, data = df1[df1$group == "ref",], 
                             tau = seq(0.01, 0.99, by = 0.01))$fitted.values
quants <- colMeans(quant_matrix)

# Rearrange the data into the necessary format
df1 %>%
  slice(rep(1:n(), each = 99)) %>% 
  mutate(fpr = rep(seq(0.01, 0.99, by = 0.01), 200)) %>% 
  mutate(q = rep(quants, 200)) %>% 
  mutate(phi_inv = qnorm(fpr)) %>% 
  mutate(u = obs >= q) %>% 
  filter(group == "dis") -> df2

# Find the coefficients for the Alonzo model 
probit_coeff <- glm(u ~ phi_inv + x, data = df2, 
                    family = binomial(link = "probit"))$coefficients

# This is just so R shows us the results
probit_coeff

# Calculate some ROC values
data_frame(x = rep(c(0.2, 0.4, 0.6, 0.8),200), 
           fpr = rep(seq(0.001, 0.999, by = 0.005), 4)) %>% 
  mutate(phi_inv = qnorm(fpr)) %>% 
  mutate(ROC = pnorm(1*probit_coeff[1] + abs(probit_coeff[2])*phi_inv +
                       probit_coeff[3]*x)) %>% 
  mutate(Youden = abs(ROC - fpr)) -> df3

# Youden Index Section
df3 %>% 
  filter(x == 0.20) %>% 
  select(Youden) %>% 
  max() -> youden_list_weib$youden_alonzo[counter,1]
df3 %>% 
  filter(x == 0.40) %>% 
  select(Youden) %>% 
  max() -> youden_list_weib$youden_alonzo[counter,2]
df3 %>% 
  filter(x == 0.60) %>% 
  select(Youden) %>% 
  max() -> youden_list_weib$youden_alonzo[counter,3]
df3 %>% 
  filter(x == 0.80) %>% 
  select(Youden) %>% 
  max() -> youden_list_weib$youden_alonzo[counter,4]

# ROC finding function
big_function <- function(x){
  roc_f <- function(fpr){
    pnorm(1*probit_coeff[1] + probit_coeff[2]*qnorm(fpr) +
            probit_coeff[3]*x)
  }
  
  # Evaluate that function across the FPR's
  integrate(f = roc_f, lower = 0.001, upper = 0.999)$value
}

# Integrate all four simultaneously!
Alonzo_AUC <- sapply(c(0.2, 0.4, 0.6, 0.8), big_function)

# Save the results
matrix_list_weib$auc_alonzo[counter,] <- Alonzo_AUC

## BETA METHOD ## 
# Calculate the placement values
pv <- sapply(df1[df1$group == "dis",]$obs,function(y) (1 - mean(y >= quants)))

df1 %>% 
  filter(group == "dis") %>% 
  dplyr::select(c(x, obs)) %>% 
  mutate(pv = pv) -> df4

# The dependent variable must be in (0,1)
# Change 0 to 0.0001 and 1 to 0.9999, if necessary
if(any(df4$pv == 0)){
  df4[df4$pv == 0,]$pv <- 0.0001
}
if(any(df4$pv == 1)){
  df4[df4$pv == 1,]$pv <- 0.9999
}

# Beta regression with the placement values
beta_fit <- betareg(pv ~ x, data = df4)

# The coefficients returned by beta regression
beta0 <- beta_fit$coefficients$mean[1] # Intercept
beta1 <- beta_fit$coefficients$mean[2] # Slope
phi <- beta_fit$coefficients$precision # Precision 

# Calculate some ROC values
data_frame(x = rep(c(0.2, 0.4, 0.6, 0.8),200), 
           fpr = rep(seq(0.001, 0.999, by = 0.005), 4)) %>% 
  mutate(m = 1/(1 + exp(-beta0 - beta1*x))) %>% 
  mutate(ROC = pbeta(fpr, shape1 = m*phi, shape2 = phi*(1-m))) %>% 
  mutate(Youden = abs(ROC - fpr)) -> df5 

# Youden Index Section
df5 %>% 
  filter(x == 0.20) %>% 
  select(Youden) %>% 
  max() -> youden_list_weib$youden_beta[counter,1]
df5 %>% 
  filter(x == 0.40) %>% 
  select(Youden) %>% 
  max() -> youden_list_weib$youden_beta[counter,2]
df5 %>% 
  filter(x == 0.60) %>% 
  select(Youden) %>% 
  max() -> youden_list_weib$youden_beta[counter,3]
df5 %>% 
  filter(x == 0.80) %>% 
  select(Youden) %>% 
  max() -> youden_list_weib$youden_beta[counter,4]


# Finding the ROC
big_function2 <- function(x){
  roc_beta <- function(fpr){
    mu <- 1/(1 + exp(-beta0 - beta1*x))
    a <- mu * phi 
    b <- phi * (1 - mu)
    
    pbeta(fpr, shape1 = a, shape2 = b)
  }
  
  # Evaluate that function across the FPR's
  integrate(f = roc_beta, lower = 0.001, upper = 0.999)$value
}

# Integrate all four simultaneously!
beta_AUC <- sapply(c(0.2, 0.4, 0.6, 0.8), big_function2)

# Save the results
matrix_list_weib$auc_beta[counter,] <- beta_AUC

## LEHMANN METHOD ## 
# Find the survival curves
S <- Surv(df1$obs)
fit <- survfit(S ~ df1$group)

df_surv <- data_frame(time = fit$time, 
                      survival = fit$surv,
                      group = rev(df1$group))


# Check the proportional hazards assumption
coxx <- coxph(formula = Surv(obs) ~ as.factor(ind) + x 
              + as.factor(ind)*x,
              data = df1)
coxx_px <- cox.zph(coxx)

# Cox proportional hazards model to estimate beta 
# Which can be used to estimate theta
beta1 <- unname(coxx$coefficients[1])
beta3 <- unname(coxx$coefficients[3])

# Calculate some ROC values
data_frame(x = rep(c(0.2, 0.4, 0.6, 0.8),200), 
           fpr = rep(seq(0.001, 0.999, by = 0.005), 4))  %>% 
  mutate(theta = exp(beta1 + beta3*x)) %>% 
  mutate(ROC = fpr^theta) %>% 
  mutate(Youden = abs(ROC - fpr)) -> df6

# Youden Index Section
df6 %>% 
  filter(x == 0.20) %>% 
  select(Youden) %>% 
  max() -> youden_list_weib$youden_lehmann[counter,1]
df6 %>% 
  filter(x == 0.40) %>% 
  select(Youden) %>% 
  max() -> youden_list_weib$youden_lehmann[counter,2]
df6 %>% 
  filter(x == 0.60) %>% 
  select(Youden) %>% 
  max() -> youden_list_weib$youden_lehmann[counter,3]
df6 %>% 
  filter(x == 0.80) %>% 
  select(Youden) %>% 
  max() -> youden_list_weib$youden_lehmann[counter,4]


# Finding the ROC
big_function3 <- function(x){
  roc_lehmann <- function(fpr){
    theta <- exp(beta1 + beta3*x)
    
    fpr^theta
  }
  
  # Evaluate that function across the FPR's
  integrate(f = roc_lehmann, lower = 0.001, upper = 0.999)$value
}

# Integrate all four simultaneously!
Lehmann_AUC <- sapply(c(0.2, 0.4, 0.6, 0.8), big_function3)

# Save the results
matrix_list_weib$auc_lehmann[counter,] <- Lehmann_AUC
}

########################################
##          Calculate MSE             ##
########################################
## NORMAL CASE ##
sq_err_normal <- list(alonzo = matrix(NA, nsim, 4),
                     beta = matrix(NA, nsim, 4),
                     lehmann = matrix(NA, nsim, 4))
for(i in 1:3){
  for(j in 1:nsim){
    sq_err_normal[[i]][j,] <- (matrix_list_normal[[i]][j,] - true_normal_auc)^2
  }
}

# MSE
normal_mse <- matrix(NA, nrow = 3, ncol = 4)
rownames(normal_mse) <- c("Alonzo","Beta","Lehmann")
colnames(normal_mse) <- c("x = 0.2", "x = 0.4", "x = 0.6", "x = 0.8")

for(i in 1:3){
  normal_mse[i,] <- colMeans(sq_err_normal[[i]])
}

## EXTREME VALUE CASE ##
sq_err_exval <- list(alonzo = matrix(NA, nsim, 4),
                              beta = matrix(NA, nsim, 4),
                              lehmann = matrix(NA, nsim, 4))
for(i in 1:3){
  for(j in 1:nsim){
  sq_err_exval[[i]][j,] <- (matrix_list_exval[[i]][j,] - true_exval_auc)^2
  }
}

# MSE
exval_mse <- matrix(NA, nrow = 3, ncol = 4)
rownames(exval_mse) <- c("Alonzo","Beta","Lehmann")
colnames(exval_mse) <- c("x = 0.2", "x = 0.4", "x = 0.6", "x = 0.8")

for(i in 1:3){
  exval_mse[i,] <- colMeans(sq_err_exval[[i]])
}

## WEIBULL CASE ##
sq_err_weib <- list(alonzo = matrix(NA, nsim, 4),
                     beta = matrix(NA, nsim, 4),
                     lehmann = matrix(NA, nsim, 4))
for(i in 1:3){
  for(j in 1:nsim){
    sq_err_weib[[i]][j,] <- (matrix_list_weib[[i]][j,] - true_weib_auc)^2
  }
}

# MSE
weib_mse <- matrix(NA, nrow = 3, ncol = 4)
rownames(weib_mse) <- c("Alonzo","Beta","Lehmann")
colnames(weib_mse) <- c("x = 0.2", "x = 0.4", "x = 0.6", "x = 0.8")

for(i in 1:3){
  weib_mse[i,] <- colMeans(sq_err_weib[[i]])
}

# Youden Index Section
# Weibull 
weib_true_maxyouden <- vector()
for(m in 1:4){
  data_frame(fpr = seq(0.01, 0.99, by = 0.01), x = rep(x2[m], 99)) %>% 
    mutate(roc = fpr^(3/7 - 0.1*x)) %>% 
    mutate(youden = abs(roc - fpr)) %>% 
    select(youden) %>% 
    max() -> weib_true_maxyouden[m]
}

# Extreme Value
exval_true_maxyouden <- vector()
for(m in 1:4){
data_frame(fpr = seq(0.01, 0.99, by = 0.01), x = rep(x2[m], 99)) %>% 
  mutate(roc = 1 - exp(-exp(log(-log(1 - fpr))) - (0.5 + x)/1.5)) %>% 
  mutate(youden = abs(roc - fpr)) %>% 
  select(youden) %>% 
  max() -> exval_true_maxyouden[m]
}

# Normal
normal_true_maxyouden <- vector()
for(m in 1:4){
  a <- ((2 + 4*x2[m]) - (1.5 + 0*x2[m]))/1.5
  b <- 1.5/1.5
data_frame(fpr = seq(0.01, 0.99, by = 0.01)) %>% 
  mutate(roc = pnorm(a + b*qnorm(fpr))) %>% 
  mutate(youden = abs(roc - fpr)) %>% 
  select(youden) %>% 
  max() -> normal_true_maxyouden[m]
}

## NORMAL CASE ##
sq_err_normal_youden <- list(alonzo = matrix(NA, nsim, 4),
                      beta = matrix(NA, nsim, 4),
                      lehmann = matrix(NA, nsim, 4))
for(i in 1:3){
  for(j in 1:nsim){
    sq_err_normal_youden[[i]][j,] <- (youden_list_normal[[i]][j,] - normal_true_maxyouden)^2
  }
}

# MSE
normal_mse_youden <- matrix(NA, nrow = 3, ncol = 4)
rownames(normal_mse_youden) <- c("Alonzo","Beta","Lehmann")
colnames(normal_mse_youden) <- c("x = 0.2", "x = 0.4", "x = 0.6", "x = 0.8")

for(i in 1:3){
  normal_mse_youden[i,] <- colMeans(sq_err_normal_youden[[i]])
}

## EXTREME VALUE CASE ##
sq_err_exval_youden <- list(alonzo = matrix(NA, nsim, 4),
                             beta = matrix(NA, nsim, 4),
                             lehmann = matrix(NA, nsim, 4))
for(i in 1:3){
  for(j in 1:nsim){
    sq_err_exval_youden[[i]][j,] <- (youden_list_exval[[i]][j,] - exval_true_maxyouden)^2
  }
}

# MSE
exval_mse_youden <- matrix(NA, nrow = 3, ncol = 4)
rownames(exval_mse_youden) <- c("Alonzo","Beta","Lehmann")
colnames(exval_mse_youden) <- c("x = 0.2", "x = 0.4", "x = 0.6", "x = 0.8")

for(i in 1:3){
  exval_mse_youden[i,] <- colMeans(sq_err_exval_youden[[i]])
}

## WEIBULL CASE ##
sq_err_weib_youden <- list(alonzo = matrix(NA, nsim, 4),
                             beta = matrix(NA, nsim, 4),
                             lehmann = matrix(NA, nsim, 4))
for(i in 1:3){
  for(j in 1:nsim){
    sq_err_weib_youden[[i]][j,] <- (youden_list_weib[[i]][j,] - weib_true_maxyouden)^2
  }
}

# MSE
weib_mse_youden <- matrix(NA, nrow = 3, ncol = 4)
rownames(weib_mse_youden) <- c("Alonzo","Beta","Lehmann")
colnames(weib_mse_youden) <- c("x = 0.2", "x = 0.4", "x = 0.6", "x = 0.8")

for(i in 1:3){
  weib_mse_youden[i,] <- colMeans(sq_err_weib_youden[[i]])
}


## SIMULATION RESULTS ##
results <- list(Weibull = list(weib_mse, weib_mse_youden),
     ExtremeValue = list(exval_mse, exval_mse_youden),
     Normal = list(normal_mse, normal_mse_youden))
