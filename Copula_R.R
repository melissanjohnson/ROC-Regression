# Function to sample from a Weibull with copulas
copula <- function(d,n,theta){
  # Sampling function
  cop <- function(n, theta){
    u <- vector()
    u[1] <- runif(1)
    s <- 0
    for(j in 2:n){
      s <- s + u[(j-1)]^(-theta) - 1
      v <- runif(1)
      u[j] <- ((1 + s)*v^(-theta/(1 + (j-1)*theta)) - s)^-(1/theta)
    }
    u
  }
  # Take more than one sample
  replicate(d, cop(n,theta))
}

# Sample using the copula
sample_matrix <- copula(5, 100, 2)

# If I had multiple samples, I could turn them all 
# into a sample from the normal this way
# I could do this with the other two distributions too
normal_samples <- apply(sample_matrix, 2, FUN = qnorm)


# I am going to use one sample
# And turn it into a sample from all 3 distributions
# That way I can look at graphs and stuff
df <- as_data_frame(copula(1,100,2))
colnames(df) <- "x"

df %>% 
  mutate(w = qweibull(x,1)) %>% # Sample from Weibull
  mutate(ev = -1.5*log(x)) %>% # Sample from Extreme Value
  mutate(n = qnorm(x))  -> df # Sample from Normal

# Plot the three distributions on the same graph
df %>% 
  gather(dist, val, -x) %>% 
  ggplot() +
    geom_density(aes(val, color = dist)) +
    labs(x = "x", y = "Density", 
         title = "Normal, Weibull, and Extreme Value Distributions
         Sampled Using the Clayton Copula") +
    scale_color_discrete(name ="Distribution", 
                         labels = c("Extreme Value", "Normal","Weibull")) +
    theme_bw()  
  

