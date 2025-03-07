library(tidyverse)
library(gridExtra)

tests = tibble()


beta.statistics = function(a, b) {
  mean_value = a / (a + b)
  variance = (a * b) / (((a + b)^2) * (a + b + 1))
  skewness = 2 * (b - a) * sqrt(a + b + 1) / ((a + b + 2) * sqrt(a * b))
  excess_kurtosis = 6 * (((a - b)^2) * (a + b + 1) - a * b * (a + b + 2)) / ((a * b) * (a + b + 2) * (a + b + 3))
  
  
  return(tibble(a = a, b = b, mean = mean_value, variance = variance, skewness = skewness, excess_kurtosis = excess_kurtosis))
}

beta.total = tibble()

beta25 = beta.statistics(2,5)
beta55 = beta.statistics(5,5)
beta52 = beta.statistics(5,2)
beta.5.5 = beta.statistics(.5,.5)



beta.total = beta25 |>
  rbind(beta55) |>
  rbind(beta52) |>
  rbind(beta.5.5) |>
  mutate(Method = "Formula")


beta.graph = function(a,b){
  q1.fig.dat <- tibble(x = seq(-0.25, 1.25, length.out=1000))|>   # generate a grid of points
    mutate(beta.pdf = dbeta(x, a, b),                      # compute the beta PDF
           norm.pdf = dnorm(x,                                    # Gaussian distribution with
                            mean = a/(a+b),            # same mean and variance
                            sd = sqrt((a*b)/((a+b)^2*(a+b+1)))))
  ggplot(data= q1.fig.dat)+                                              # specify data
    geom_line(aes(x=x, y=beta.pdf)) +                 # plot beta dist
    geom_hline(yintercept=0)+                                            # plot x axis
    theme_bw()+                                                          # change theme
    xlab("x")+                                                           # label x axis
    ylab("Density")+                                                     # label y axis
    scale_color_manual("", values = c("black", "grey"))+                 # change colors
    theme(legend.position = "bottom")                                    # move legend to bottom
}

l = beta.graph(2,5)
m = beta.graph(0.5,0.5)
n = beta.graph(5,5)
o = beta.graph(5,2)

grid.arrange(l,m,n,o, ncol = 2)



beta.moment = function(a, b, k, centered) {
  # Define the integrand function for the raw moment
  integrand = function(x) {(x^k) * dbeta(x, a, b)}
  # Compute the expected value (raw moment)
  E = integrate(integrand, 0, 1)$value
  
  # Check if the moment is centered
  if (identical(centered, FALSE)) {
    return(E)  # Return the raw moment
  } else {
    # Define the integrand for calculating the mean
    integrand_mean = function(x) {x * dbeta(x, a, b)}
    mean = integrate(integrand_mean, 0, 1)$value
    
    # Define the integrand for the centered moment
    integrand_centered = function(x) {((x - mean)^k) * dbeta(x, a, b)}
    Ec = integrate(integrand_centered, 0, 1)$value
    
    return(Ec)  # Return the centered moment
  }
}

beta.summarize.noncenter = function(a,b){
  a = a 
  b = b
  mean = beta.moment(a,b,1,F)
  variance = beta.moment(a,b,2,T)
  skewness = (beta.moment(a,b,3,T))/(beta.moment(a,b,2,T)^(3/2))
  excess_kurtosis = beta.moment(a,b,4,T)/((beta.moment(a,b,2,T))^2) - 3
  return(tibble(a = a, b = b, mean = mean, variance = variance, skewness = skewness, excess_kurtosis = excess_kurtosis))
}

derived25 = beta.summarize.noncenter(2,5)
derived55 = beta.summarize.noncenter(5,5)
derived52 = beta.summarize.noncenter(5,2)
derived.5.5 = beta.summarize.noncenter(.5,.5)






#TASK 3



#BETA A = 2 B = 5
library(tidyverse)
library(e1071)
set.seed(7272) # Set seed so we all get the same results.
sample.size <- 500 # Specify sample details
alpha <- 2
beta <- 5
x = seq(-0.25, 1.25, length.out=1000)
beta.sample <- rbeta(n = sample.size,  # sample size
                     shape1 = alpha,   # alpha parameter
                     shape2 = beta)    # beta parameter
beta.df = tibble(sample = beta.sample)

beta.pdf = tibble(x = x, beta = dbeta(x, alpha, beta))

summary.statistics.of.beta.25 = tibble(beta.sample) |>
  summarize(sample.mean = mean(beta.sample),
            sample.variance = var(beta.sample),
            sample.sd = sd(beta.sample),
            sample.skewness = skewness(beta.sample),
            sample.excess.kurtosis = kurtosis(beta.sample))


plot1 = ggplot(beta.df)+
  geom_histogram(aes(x = sample, y = after_stat(density)), bins = 34, fill = "purple")+
  geom_density(aes(x=sample, color = "Density of Sample"))+
  geom_hline(yintercept = 0)+
  xlab("Sample")+
  ylab("Density")+
  theme_minimal()+
  geom_line(data = beta.pdf, aes(x = x, y = beta, color = "Beta (2,5)"))+
  guides(color = guide_legend(title = NULL))
  

#BETA A = 5 B = 5
library(tidyverse)
library(e1071)
set.seed(7272) # Set seed so we all get the same results.
sample.size <- 500 # Specify sample details
alpha <- 5
beta <- 5
x = seq(-0.25, 1.25, length.out=1000)
beta.sample <- rbeta(n = sample.size,  # sample size
                     shape1 = alpha,   # alpha parameter
                     shape2 = beta)    # beta parameter
beta.df = tibble(sample = beta.sample)

beta.pdf = tibble(x = x, beta = dbeta(x, alpha, beta))

summary.statistics.of.beta.55 = tibble(beta.sample) |>
  summarize(sample.mean = mean(beta.sample),
            sample.variance = var(beta.sample),
            sample.sd = sd(beta.sample),
            sample.skewness = skewness(beta.sample),
            sample.excess.kurtosis = kurtosis(beta.sample))


plot2 = ggplot(beta.df)+
  geom_histogram(aes(x = sample, y = after_stat(density)), bins = 34, fill = "purple")+
  geom_density(aes(x=sample, color = "Density of Sample"))+
  geom_hline(yintercept = 0)+
  xlab("Sample")+
  ylab("Density")+
  theme_minimal()+
  geom_line(data = beta.pdf, aes(x = x, y = beta, color = "Beta (5,5)"))+
  guides(color = guide_legend(title = NULL))
            
#BETA A = 5 B = 2
library(tidyverse)
library(e1071)
set.seed(7272) # Set seed so we all get the same results.
sample.size <- 500 # Specify sample details
alpha <- 5
beta <- 2
x = seq(-0.25, 1.25, length.out=1000)
beta.sample <- rbeta(n = sample.size,  # sample size
                     shape1 = alpha,   # alpha parameter
                     shape2 = beta)    # beta parameter
beta.df = tibble(sample = beta.sample)

beta.pdf = tibble(x = x, beta = dbeta(x, alpha, beta))

summary.statistics.of.beta.52 = tibble(beta.sample) |>
  summarize(sample.mean = mean(beta.sample),
            sample.variance = var(beta.sample),
            sample.sd = sd(beta.sample),
            sample.skewness = skewness(beta.sample),
            sample.excess.kurtosis = kurtosis(beta.sample))


plot3 = ggplot(beta.df)+
  geom_histogram(aes(x = sample, y = after_stat(density)), bins = 34, fill = "purple")+
  geom_density(aes(x=sample, color = "Density of Sample"))+
  geom_hline(yintercept = 0)+
  xlab("Sample")+
  ylab("Density")+
  theme_minimal()+
  geom_line(data = beta.pdf, aes(x = x, y = beta, color = "Beta (5,2)"))+
  guides(color = guide_legend(title = NULL))
 
#BETA A = 0.5 B = 0.5
library(tidyverse)
library(e1071)
set.seed(7272) # Set seed so we all get the same results.
sample.size <- 500 # Specify sample details
alpha <- 0.5
beta <- 0.5
x = seq(-0.25, 1.25, length.out=1000)
beta.sample <- rbeta(n = sample.size,  # sample size
                     shape1 = alpha,   # alpha parameter
                     shape2 = beta)    # beta parameter
beta.df = tibble(sample = beta.sample)

beta.pdf = tibble(x = x, beta = dbeta(x, alpha, beta))

summary.statistics.of.beta.0.5.0.5 = tibble(beta.sample) |>
  summarize(sample.mean = mean(beta.sample),
            sample.variance = var(beta.sample),
            sample.sd = sd(beta.sample),
            sample.skewness = skewness(beta.sample),
            sample.excess.kurtosis = kurtosis(beta.sample))


plot4 = ggplot(beta.df)+
  geom_histogram(aes(x = sample, y = after_stat(density)), bins = 34, fill = "purple")+
  geom_density(aes(x=sample, color = "Density of Sample"))+
  geom_hline(yintercept = 0)+
  xlab("Sample")+
  ylab("Density")+
  theme_minimal()+
  geom_line(data = beta.pdf, aes(x = x, y = beta, color = "Beta (0.5,0.5)"))+
  guides(color = guide_legend(title = NULL))

grid.arrange(plot1,plot2,plot3,plot4,ncol = 2)


#STEP 4############
library(cumstats)

set.seed(7272)
beta.sample = rbeta(n=500, shape1 = 2, shape2 = 5)

ggplot(tibble(beta.sample))+
  geom_line(x = 1:500, y = cummean(beta.sample))

