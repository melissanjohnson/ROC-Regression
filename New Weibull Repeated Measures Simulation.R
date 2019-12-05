library(tidyverse)
library(glmmTMB)
library(betareg)
library(survival)
library(survminer)
library(PRROC)
############################################
## Function to Sample from Clayton Copula ##
############################################
copula <- function(nsubject,ntime,clayton_theta){
  # Sampling function
  cop <- function(ntime, clayton_theta){
    u <- vector()
    u[1] <- runif(1)
    s <- 0
    for(j in 2:ntime){
      s <- s + u[(j-1)]^(-clayton_theta) - 1
      v <- runif(1)
      u[j] <- ((1 + s)*v^(-clayton_theta/(1 + (j-1)*clayton_theta)) - s)^-(1/clayton_theta)
    }
    u
  }
  # Take more than one sample
  t(replicate(nsubject, cop(ntime,clayton_theta)))
}

########################################
##      SIMULATION STARTS HERE        ##
########################################
# Which case?
ss <- 1
# How many observations per subject?
ntime <- 5
# How many subjects (per group, not total)?
nsub <- 40
# How many iterations of the simulation?
nsim <- 1000
# What were the true AUC values I was aiming for?
# I calculated these numbers elsewhere
trueAUC <- list(c(0.50, 0.50, 0.50, 0.50, 0.50),
                c(0.50, 0.50, 0.50, 0.50, 0.70),
                c(0.50, 0.50, 0.50, 0.70, 0.70),
                c(0.50, 0.50, 0.70, 0.70, 0.70),
                c(0.50, 0.70, 0.70, 0.70, 0.70))
# What Youden was I aiming for?
truetheta <- list(c(1, 1, 1, 1, 1),
                  c(1, 1, 1, 1,3/7),
                  c(1, 1, 1, 3/7,3/7),
                  c(1, 1, 3/7, 3/7,3/7),
                  c(1, 3/7, 3/7, 3/7,3/7))


trueYouden <- list(vector(length = ntime),vector(length = ntime),
                   vector(length = ntime),vector(length = ntime),
                   vector(length = ntime))
for(i in 1:ntime){
  tibble(x = seq(0.01, 0.99, by = 0.01)) %>% 
    mutate(y = x^truetheta[[ss]][i]) %>% 
    mutate(youden_index = y - x) %>% 
    select(youden_index) %>% 
    max() -> trueYouden[[ss]][i]
}

# A list to hold the generated AUCs
AUC <- list(Lehmann = matrix(NA, nrow = nsim, ncol = ntime),
            Beta = matrix(NA, nrow = nsim, ncol = ntime))
# A list ot hold the generated Youden indeces
Youden <- list(Lehmann = matrix(NA, nrow = nsim, ncol = ntime),
               Beta = matrix(NA, nrow = nsim, ncol = ntime))

for(iter in 1:nsim){
########################################
##           Generate Data            ##
########################################
# This is the theta for the Clayton copula. 
# Completely unrelated to the Weibull theta!!
ctheta <- 1
# Random intercept
# Aka the starting point for each subject
b0 <- runif(nsub, 1, 4)
# This is the part that determines the AUC at each time
# This is 1 + k*gamma_k all combined into one step!
gammak <- list(c(0,0, 0, 0, 0),
               c(0,0, 0,0, 7/3),
               c(0,0,0, 7/3, 7/3),
               c(0,0, 7/3, 7/3, 7/3),
               c(0,7/3, 7/3, 7/3, 7/3))
# Generate dependent uniforms from the copula
u1 <- copula(nsub, ntime, ctheta)
# Turn those into dependent Weibull errors for the diseased group
epsilon1 <- qweibull(u1, 1, 2.5)
# Generate more dependent uniforms from the copula
u0 <- copula(nsub, ntime, ctheta)
# Turn those into dependent Weibull errors for the reference group
epsilon0 <- qweibull(u0, 1, 2.5)
# These are empty matrices that I am going to put the observations in
y0 <- matrix(NA, nrow = nsub, ncol = ntime)
y1 <- matrix(NA, nrow = nsub, ncol = ntime)
# 1 + gamma_k*k + epsilon_ijk
for(i in 1:nsub){
  y1[i,] <- gammak[[ss]] + epsilon1[i,]
}
# Now add the intercept
# This was done in two parts because of the dimensions of my vectors
# I was just lazy and couldn't think of a better way to do it...
# There is really no significance to separating it!
for(i in 1:ntime){
  y1[,i] <- y1[,i] + b0
}
# b0i + epsilon_ijk for the reference group data 
for(i in 1:ntime){
  y0[,i] <- b0 + epsilon0[,i]
}
# Put everything in a nice, labeled data frame
tibble(ref = c(y0), 
       dis = c(y1)) %>% 
  mutate(time = rep(1:ntime, each = nsub)) %>% 
  mutate(sub = rep(1:nsub, ntime)) %>% 
  gather(group, obs, -c(time, sub)) -> df


########################################
##       Lehmann ROC Regression       ##
########################################
# A vector to hold the AUC at each time
AUC1 <- vector()
# A for loop to find the ROC/AUC at each time
for(i in 1:ntime){
  df %>% 
    mutate(ind = as.factor(case_when(
      group == "ref" ~ 0,
      group == "dis" ~ 1
    )))-> tempdf
  
  # What patient?
  ID <- tempdf$sub
  
  # Diseased (1) or reference (0)?
  D <- tempdf$ind
  
  # Turn matrices into vectors
  y <- tempdf$obs
  
  # Time
  k <- as.numeric(tempdf$time)
  
  #### Lehmann Regression #####
  # Lehmann ROC using Cox PH Regression
  fit <- coxph(Surv(y) ~ D + D*k + cluster(ID))
  
  # Estimated (Lehmann) theta
  theta_est <- exp(unname(coef(fit)[1] + coef(fit)[3]*i))
  
  # Store the AUC
  AUC1[i] <- 1/(1 + theta_est)
  # Youden Index Lehmann
  tibble(x = seq(0.01, 0.99, by = 0.025)) %>% 
    mutate(y = x^theta_est) %>% 
    mutate(youden_index = y - x) %>% 
    select(youden_index) %>% 
    max() -> Youden$Lehmann[iter,i]
   }
AUC$Lehmann[iter,] <- AUC1


########################################
##       Beta ROC Regression          ##
########################################
# A vector to store the AUC at each time 
AUC2 <- vector()
# A for loop to find the AUC/ROC at each time
for(i in 1:ntime){
  df %>% 
    filter(time == i) -> tempdf2
  
  quant_matrix <- quantreg::rq(obs ~ 1, data = tempdf2[tempdf2$group == "ref",], 
                               tau = seq(0.01, 0.99, by = 0.025))$fitted.values
  quants <- colMeans(quant_matrix)
  
  # Calculate the placement values
  pv <- sapply(tempdf2[tempdf2$group == "dis",]$obs,function(y) (1 - mean(y >= quants)))
  
  # Select the diseased ones and make a new data set with them
  tempdf2 %>% 
    filter(group == "dis") %>% 
    select(c(obs, sub)) %>% 
    mutate(pv = pv) -> df2
  df2
  
  # The dependent variable must be in (0,1)
  # Change 0 to 0.0001 and 1 to 0.9999, if necessary
  if(any(df2$pv == 0)){
    df2[df2$pv == 0,]$pv <- 0.0001
  }
  if(any(df2$pv == 1)){
    df2[df2$pv == 1,]$pv <- 0.9999
  }
  
  # Beta regression with the placement values
  beta_fit <- betareg(pv ~ 1, data = df2)
  
  # The coefficients returned by beta regression
  beta0 <- beta_fit$coefficients$mean[1] # Intercept
  phi <- beta_fit$coefficients$precision # Precision 
  
  # ROC function
  roc_beta <- function(fpr){
    mu <- 1/(1 + exp(-beta0))
    a <- mu * phi 
    b <- phi * (1 - mu)
    
    pbeta(fpr, shape1 = a, shape2 = b)
  } 
# Store the AUC
AUC2[i] <- as.numeric(integrate(roc_beta, lower = 0.001, upper = 0.999)[1])
# Youden Index Beta
tibble(x = seq(0.01, 0.99, by = 0.025)) %>% 
  mutate(y = roc_beta(x)) %>% 
  mutate(youden_index = y - x) %>% 
  select(youden_index) %>% 
  max() -> Youden$Beta[iter,i]
}
AUC$Beta[iter,] <- AUC2
}

########################################
##          Calculate MSE             ##
########################################
# List to store the squared errors
SE <- list(Lehmann = matrix(NA, nrow = nsim, ncol = ntime),
           Beta = matrix(NA, nrow = nsim, ncol = ntime))
for(i in 1:ntime){
  SE$Lehmann[,i] <- (AUC$Lehmann[,i] - trueAUC[[ss]][i])^2
  SE$Beta[,i] <- (AUC$Beta[,i] - trueAUC[[ss]][i])^2
}

# Matrix to store the MSE
MSE <- matrix(NA, nrow = 2, ncol = ntime)
for(i in 1:ntime){
  MSE[1,i] <- mean(SE$Lehmann[,i])
  MSE[2,i] <- mean(SE$Beta[,i])
}

rownames(MSE) <- c("Lehmann", "Beta")
colnames(MSE) <- paste("Time",1:ntime)

########################################
##       Calculate MSE (Youden)       ##
########################################
# List to store the squared errors
SE2 <- list(Lehmann = matrix(NA, nrow = nsim, ncol = ntime),
           Beta = matrix(NA, nrow = nsim, ncol = ntime))
for(i in 1:ntime){
  SE2$Lehmann[,i] <- (Youden$Lehmann[,i] - trueYouden[[ss]])^2
  SE2$Beta[,i] <- (Youden$Beta[,i] - trueYouden[[ss]])^2
}

# Matrix to store the MSE
MSE2 <- matrix(NA, nrow = 2, ncol = ntime)
for(i in 1:ntime){
  MSE2[1,i] <- mean(SE2$Lehmann[,i])
  MSE2[2,i] <- mean(SE2$Beta[,i])
}

rownames(MSE2) <- c("Lehmann", "Beta")
colnames(MSE2) <- paste("Time",1:ntime)

########################################
##        Summary of Results          ##
########################################
result_table <- list(MSE_for_AUC = MSE, 
     MSE_for_Youden = MSE2)

result_table
    
########################################
##   These are plots for the paper    ##
########################################
# The observations across time for each subject in the treatment group
df %>% 
  filter(group == "ref") %>% 
  ggplot() + geom_line(aes(time, obs, color = as.factor(sub))) +
  guides(color = FALSE) +
  theme_bw() +
  labs(x = "k", y = "y")
# One individual in the treatment group, and one in the control group over time
df %>% 
  filter(sub == 1) %>% 
  ggplot() + geom_line(aes(time, obs, color = as.factor(group))) +
  theme_bw() +
  labs(x = "k", y = "y") +
  scale_color_discrete("Group", labels = c("Treatment", "Control"))
# Densities at time 1
df %>% 
  filter(time == 1) %>% 
  ggplot() + geom_density(aes(obs, color = group)) +
  theme_bw() +
  labs(x = "y", y = "") +
  scale_color_discrete("Group", labels = c("Treatment", "Control"))
# Densities at time 2
df %>% 
  filter(time == 2) %>% 
  ggplot() + geom_density(aes(obs, color = group)) +
  theme_bw() +
  labs(x = "y", y = "") +
  scale_color_discrete("Group", labels = c("Treatment", "Control"))
# Densities at time 3
df %>% 
  filter(time == 3) %>% 
  ggplot() + geom_density(aes(obs, color = group)) +
  theme_bw() +
  labs(x = "y", y = "") +
  scale_color_discrete("Group", labels = c("Treatment", "Control"))
# Densities at time 4
df %>% 
  filter(time == 4) %>% 
  ggplot() + geom_density(aes(obs, color = group)) +
  theme_bw() +
  labs(x = "y", y = "") +
  scale_color_discrete("Group", labels = c("Treatment", "Control"))
# Densities at time 5
df %>% 
  filter(time == 5) %>% 
  ggplot() + geom_density(aes(obs, color = group)) +
  theme_bw() +
  labs(x = "y", y = "") +
  scale_color_discrete("Group", labels = c("Treatment", "Control"))


tl <- vector()
for(i in 1:ntime){
  df %>% 
    filter(time == i)  %>% 
    mutate(ind = as.factor(case_when(
      group == "ref" ~ 0,
      group == "dis" ~ 1
    )))-> tempdf
  
  # What patient?
  ID <- tempdf$sub
  
  # Diseased (1) or reference (0)?
  D <- tempdf$ind
  
  # Turn matrices into vectors
  y <- tempdf$obs
  
  #### Lehmann Regression #####
  # Lehmann ROC using Cox PH Regression
  fit <- coxph(Surv(y) ~ D + cluster(ID))
  
  # Estimated (Lehmann) theta
  tl[i] <- exp(unname(coef(fit)))
}

ll <- list()
for(i in 1:ntime){
  df %>% 
    filter(time == i) -> tempdf2
  
  quant_matrix <- quantreg::rq(obs ~ 1, data = tempdf2[tempdf2$group == "ref",], 
                               tau = seq(0.01, 0.99, by = 0.01))$fitted.values
  quants <- colMeans(quant_matrix)
  
  # Calculate the placement values
  pv <- sapply(tempdf2[tempdf2$group == "dis",]$obs,function(y) (1 - mean(y >= quants)))
  
  # Select the diseased ones and make a new data set with them
  tempdf2 %>% 
    filter(group == "dis") %>% 
    select(c(obs, sub)) %>% 
    mutate(pv = pv) -> df2
  df2
  
  # The dependent variable must be in (0,1)
  # Change 0 to 0.0001 and 1 to 0.9999, if necessary
  if(any(df2$pv == 0)){
    df2[df2$pv == 0,]$pv <- 0.0001
  }
  if(any(df2$pv == 1)){
    df2[df2$pv == 1,]$pv <- 0.9999
  }
  
  # Beta regression with the placement values
  beta_fit <- betareg(pv ~ 1, data = df2)
  
  # The coefficients returned by beta regression
  beta0 <- beta_fit$coefficients$mean[1] # Intercept
  phi <- beta_fit$coefficients$precision # Precision 
  
  
  # ROC function
  roc_beta <- function(fpr){
    mu <- 1/(1 + exp(-beta0))
    a <- mu * phi 
    b <- phi * (1 - mu)
    
    pbeta(fpr, shape1 = a, shape2 = b)
  } 
  
  ll[[i]] <- beta_fit$coefficients
}

roc_beta1 <- function(fpr){
  beta0 <- ll[[1]]$mean
  phi <- ll[[1]]$precision
  mu <- 1/(1 + exp(-beta0))
  a <- mu * phi 
  b <- phi * (1 - mu)
  
  pbeta(fpr, shape1 = a, shape2 = b)
  
}

# Time 1
tibble(fpr = seq(0.001, 0.999, by = 0.001)) %>% 
  mutate(rocL = fpr^tl[1]) %>% 
  mutate(rocB = roc_beta1(fpr)) %>% 
  mutate(rocT = fpr) %>% 
  gather(method, roc, - fpr) %>% 
  ggplot() + geom_line(aes(fpr, roc, color = method)) +
  theme_bw() +
  labs(x = "x", y = "y") +
  scale_color_discrete("Method", labels = c("Beta", "Lehmann", "True"))

# Time 2
roc_beta2 <- function(fpr){
  beta0 <- ll[[2]]$mean
  phi <- ll[[2]]$precision
  mu <- 1/(1 + exp(-beta0))
  a <- mu * phi 
  b <- phi * (1 - mu)
  
  pbeta(fpr, shape1 = a, shape2 = b)
  
}

tibble(fpr = seq(0.001, 0.999, by = 0.001)) %>% 
  mutate(rocL = fpr^tl[2]) %>% 
  mutate(rocB = roc_beta2(fpr)) %>% 
  mutate(rocT = fpr^truetheta[2]) %>% 
  gather(method, roc, - fpr) %>% 
  ggplot() + geom_line(aes(fpr, roc, color = method)) +
  theme_bw() +
  labs(x = "x", y = "y") +
  scale_color_discrete("Method", labels = c("Beta", "Lehmann", "True"))

# Time 3
roc_beta3 <- function(fpr){
  beta0 <- ll[[3]]$mean
  phi <- ll[[3]]$precision
  mu <- 1/(1 + exp(-beta0))
  a <- mu * phi 
  b <- phi * (1 - mu)
  
  pbeta(fpr, shape1 = a, shape2 = b)
  
}

tibble(fpr = seq(0.001, 0.999, by = 0.001)) %>% 
  mutate(rocL = fpr^tl[3]) %>% 
  mutate(rocB = roc_beta3(fpr)) %>% 
  mutate(rocT = fpr^truetheta[3]) %>% 
  gather(method, roc, - fpr) %>% 
  ggplot() + geom_line(aes(fpr, roc, color = method)) +
  theme_bw() +
  labs(x = "x", y = "y") +
  scale_color_discrete("Method", labels = c("Beta", "Lehmann", "True"))

# Time 4
roc_beta4 <-function(fpr){
  beta0 <- ll[[4]]$mean
  phi <- ll[[4]]$precision
  mu <- 1/(1 + exp(-beta0))
  a <- mu * phi 
  b <- phi * (1 - mu)
  
  pbeta(fpr, shape1 = a, shape2 = b)
  
}

tibble(fpr = seq(0.001, 0.999, by = 0.001)) %>% 
  mutate(rocL = fpr^tl[4]) %>% 
  mutate(rocB = roc_beta4(fpr)) %>% 
  mutate(rocT = fpr^truetheta[4]) %>% 
  gather(method, roc, - fpr) %>% 
  ggplot() + geom_line(aes(fpr, roc, color = method)) +
  theme_bw() +
  labs(x = "x", y = "y") +
  scale_color_discrete("Method", labels = c("Beta", "Lehmann", "True"))

# Time 5
roc_beta5 <- function(fpr){
  beta0 <- ll[[5]]$mean
  phi <- ll[[5]]$precision
  mu <- 1/(1 + exp(-beta0))
  a <- mu * phi 
  b <- phi * (1 - mu)
  
  pbeta(fpr, shape1 = a, shape2 = b)
  
}

tibble(fpr = seq(0.001, 0.999, by = 0.001)) %>% 
  mutate(rocL = fpr^tl[5]) %>% 
  mutate(rocB = roc_beta5(fpr)) %>% 
  mutate(rocT = fpr^truetheta[5]) %>% 
  gather(method, roc, - fpr) %>% 
  ggplot() + geom_line(aes(fpr, roc, color = method)) +
  theme_bw() +
  labs(x = "x", y = "y") +
  scale_color_discrete("Method", labels = c("Beta", "Lehmann", "True"))




















