library(tidyverse)
library(gridExtra)
library(e1071)
library(patchwork)
library(cumstats)
library(nleqslv)

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
  return(tibble(a = a, b = b, mean = mean, variance = variance, skewness = skewness, excess_kurtosis = excess_kurtosis, Method = "Derived"))
}

derived25 = beta.summarize.noncenter(2,5)
derived55 = beta.summarize.noncenter(5,5)
derived52 = beta.summarize.noncenter(5,2)
derived.5.5 = beta.summarize.noncenter(.5,.5)

beta.total = beta.total |>
  rbind(derived25) |>
  rbind(derived55) |>
  rbind(derived52) |>
  rbind(derived.5.5) |>
  unique()






#TASK 3



#BETA A = 2 B = 5
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

summary.statistics.of.beta25 = tibble(beta.sample) |>
  summarize(sample.mean = mean(beta.sample),
            sample.variance = var(beta.sample),
            sample.sd = sd(beta.sample),
            sample.skewness = skewness(beta.sample),
            sample.excess.kurtosis = kurtosis(beta.sample))


plot1 = ggplot(beta.df)+
  geom_histogram(aes(x = sample, y = after_stat(density)), bins = 30, fill = "purple")+
  geom_density(aes(x=sample, color = "Density of Sample"))+
  geom_hline(yintercept = 0)+
  xlab("Sample")+
  ylab("Density")+
  theme_minimal()+
  geom_line(data = beta.pdf, aes(x = x, y = beta, color = "Beta (2,5)"))+
  guides(color = guide_legend(title = NULL))
  

#BETA A = 5 B = 5
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

summary.statistics.of.beta55 = tibble(beta.sample) |>
  summarize(sample.mean = mean(beta.sample),
            sample.variance = var(beta.sample),
            sample.sd = sd(beta.sample),
            sample.skewness = skewness(beta.sample),
            sample.excess.kurtosis = kurtosis(beta.sample))


plot2 = ggplot(beta.df)+
  geom_histogram(aes(x = sample, y = after_stat(density)), bins = 30, fill = "purple")+
  geom_density(aes(x=sample, color = "Density of Sample"))+
  geom_hline(yintercept = 0)+
  xlab("Sample")+
  ylab("Density")+
  theme_minimal()+
  geom_line(data = beta.pdf, aes(x = x, y = beta, color = "Beta (5,5)"))+
  guides(color = guide_legend(title = NULL))
            
#BETA A = 5 B = 2
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

summary.statistics.of.beta52 = tibble(beta.sample) |>
  summarize(sample.mean = mean(beta.sample),
            sample.variance = var(beta.sample),
            sample.sd = sd(beta.sample),
            sample.skewness = skewness(beta.sample),
            sample.excess.kurtosis = kurtosis(beta.sample))


plot3 = ggplot(beta.df)+
  geom_histogram(aes(x = sample, y = after_stat(density)), bins = 30, fill = "purple")+
  geom_density(aes(x=sample, color = "Density of Sample"))+
  geom_hline(yintercept = 0)+
  xlab("Sample")+
  ylab("Density")+
  theme_minimal()+
  geom_line(data = beta.pdf, aes(x = x, y = beta, color = "Beta (5,2)"))+
  guides(color = guide_legend(title = NULL))
 
#BETA A = 0.5 B = 0.5
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
  geom_histogram(aes(x = sample, y = after_stat(density)), bins = 30, fill = "purple")+
  geom_density(aes(x=sample, color = "Density of Sample"))+
  geom_hline(yintercept = 0)+
  xlab("Sample")+
  ylab("Density")+
  theme_minimal()+
  geom_line(data = beta.pdf, aes(x = x, y = beta, color = "Beta (0.5,0.5)"))+
  guides(color = guide_legend(title = NULL))

(plot1 + plot2) / (plot3 + plot4)

total.summary = summary.statistics.of.beta.0.5.0.5 |>
  rbind(summary.statistics.of.beta25) |>
  rbind(summary.statistics.of.beta55) |>
  rbind(summary.statistics.of.beta52)

#STEP 4############

set.seed(7272)
beta.sample = rbeta(n=500, shape1 = 2, shape2 = 5)

# Calculate cumulative statistics
cum.statistics = tibble(
  cum_mean = cummean(beta.sample),
  cum_skew = cumskew(beta.sample),
  cum_kurt = cumkurt(beta.sample) - 3,  # Subtract 3 for excess kurtosis
  cum_var = cumvar(beta.sample)
) |>
  mutate(n = row_number())


mean.plot = ggplot(cum.statistics, aes(x=n,y=cum_mean))+
  geom_line(color = "blue")+
  geom_hline(yintercept = beta25$mean, linetype = "dashed", color = "red")+
  labs(x = "Sample Size", y = "Mean")

skew.plot = ggplot(cum.statistics, aes(x=n,y=cum_skew))+
  geom_line(color = "blue")+
  geom_hline(yintercept = beta25$skewness, linetype = "dashed", color = "red")+
  labs(x = "Sample Size", y = "Skewness")

kurt.plot = ggplot(cum.statistics, aes(x=n,y=cum_kurt))+
  geom_line(color = "blue")+
  geom_hline(yintercept = beta25$excess_kurtosis, linetype = "dashed", color = "red")+
  labs(x = "Sample Size", y = "Kurtosis")

var.plot = ggplot(cum.statistics, aes(x=n,y=cum_var))+
  geom_line(color = "blue")+
  geom_hline(yintercept = beta25$variance, linetype = "dashed", color = "red")+
  labs(x = "Sample Size", y = "Variance")


(mean.plot + skew.plot) / (kurt.plot + var.plot)



final_df <- tibble()

for (i in 2:50) {
  set.seed(7272 + i)
  

  beta_sample <- rbeta(n = 500, shape1 = 2, shape2 = 5)
  
  # Compute cumulative statistics
  simulation_cum_statistics <- tibble(
    n = seq_along(beta_sample),
    cum_mean = cummean(beta_sample),  
    cum_var = cumvar(beta_sample),  
    cum_skew = cumskew(beta_sample),  
    cum_kurt = cumkurt(beta_sample)-3,  
    simulation = i  
  )
  
  final_df = final_df |>
    bind_rows(simulation_cum_statistics) |>
    unique()
  
}

#PLOTTING

total.mean.plot = ggplot(final_df, aes(x=n,y=cum_mean,color = factor(simulation)))+
  geom_line()+
  geom_hline(yintercept = beta25$mean, linetype = "dashed", color = "red")+
  labs(x = "Sample Size", y = "Mean")+
  guides(color = "none")+
  theme_minimal()

total.skew.plot = ggplot(final_df, aes(x=n,y=cum_skew,color = factor(simulation)))+
  geom_line()+
  geom_hline(yintercept = beta25$skewness, linetype = "dashed", color = "red")+
  labs(x = "Sample Size", y = "Skewness")+
  guides(color = "none")+
  theme_minimal()

total.kurt.plot = ggplot(final_df, aes(x=n,y=cum_kurt,color = factor(simulation)))+
  geom_line()+
  geom_hline(yintercept = beta25$excess_kurtosis, linetype = "dashed", color = "red")+
  labs(x = "Sample Size", y = "Kurtosis")+
  guides(color = "none")+
  theme_minimal()

total.var.plot = ggplot(final_df, aes(x=n,y=cum_var, color = factor(simulation)))+
  geom_line()+
  geom_hline(yintercept = beta25$variance, linetype = "dashed", color = "red")+
  labs(x = "Sample Size", y = "Variance")+
  guides(color = "none")+
  theme_minimal()

(total.mean.plot + total.skew.plot) / (total.kurt.plot + total.var.plot)


#STEP 5

simulated.data = tibble()
statistics.dist = tibble()

for(i in 1:1000){
  
  set.seed(7272 + i)
  
  
  beta_sample <- rbeta(n = 500, shape1 = 2, shape2 = 5)
  
  # Compute cumulative statistics
  simulated.data <- tibble(
    mean = mean(beta_sample),  
    var = var(beta_sample),  
    skew = skewness(beta_sample),  
    kurt = kurtosis(beta_sample)-3,  
    simulation = i  
  )
  
  statistics.dist = statistics.dist |>
    bind_rows(simulated.data) |>
    unique()
}



#PLOTTING THIS DATA

mean.dist = ggplot(statistics.dist, aes(x=mean,y=after_stat(density)))+
  geom_histogram(fill = "blue")+
  theme_minimal()+
  labs(x= "Mean", y = "Occurances")+
  geom_density(aes(x=mean))

var.dist = ggplot(statistics.dist, aes(x=var,y=after_stat(density)))+
  geom_histogram(fill = "blue")+
  theme_minimal()+
  labs(x= "Variance", y = "Occurances")+
  geom_density(aes(x=var))

kurt.dist = ggplot(statistics.dist)+
  geom_histogram(aes(x=kurt, y = after_stat(density)),fill = "blue",bins=28)+
  theme_minimal()+
  labs(x= "Excess Kurtosis", y = "Occurances")+
  geom_density(aes(x=kurt))

skew.dist = ggplot(statistics.dist, aes(x=skew,y=after_stat(density)))+
  geom_histogram(fill = "blue",bins=26)+
  theme_minimal()+
  labs(x= "Skewness", y = "Occurances")+
  geom_density(aes(x=skew))

(mean.dist + skew.dist) / (kurt.dist + var.dist)




#TASK 6 LAB 8!!!!!

death.data = read_csv("realreal.csv") |>
  mutate("Deaths" = get("2022")) |>
  select("Country Name", "Deaths") |>
  mutate(prop.death = get("Deaths")/1000) |>
  filter(!is.na(prop.death))

mombeta = function(par,data){
  alpha = par[1]
  beta  = par[2]
  
  EX =  mean(data, na.rm = T)
  EX2 = mean(data^2, na.rm = T)
  
  m1 = (alpha)/(alpha + beta)
  m2 = ((alpha+1)*(alpha)) / ((alpha + beta + 1)*(alpha + beta))
  
  return(c(EX - m1, EX2 - m2))
}
 

moms = nleqslv(x = c(0.5,0.5),
        fn = mombeta,
        data = death.data$prop.death) 


alpha.hat.mom <- moms$x[1]
beta.hat.mom <- moms$x[2]

llbeta <- function(data, par, neg=F){
  alpha <- par[1]
  beta <- par[2]
  
  loglik <- sum(log(dbeta(x=data, shape1=alpha, shape2=beta)),na.rm = T)
  
  return(ifelse(neg, -loglik, loglik))
}

(mles <- optim(par = c(0.5,0.5),
               fn = llbeta,
               data = death.data$prop.death,
               neg=T))
alpha.hat.mle <- mles$par[1]
beta.hat.mle <- mles$par[2]

mleplot <- tibble(x = seq(-0, 0.02, length.out = 1000)) |> 
  mutate(y = dbeta(x, alpha.hat.mle, beta.hat.mle))

momplot <- tibble(x = seq(-0, 0.02, length.out = 1000)) |> 
  mutate(y = dbeta(x, alpha.hat.mom, beta.hat.mom))

                   
                 

ggplot()+
  geom_histogram(data = death.data, aes(x=prop.death, y = after_stat(density),color = "lightgrey"), bins = 36)+
  geom_line(data = mleplot, aes(x=x,y=y,color = "MLE Estimate"))+
  geom_line(data = momplot, aes(x=x,y=y,color = "MOM Estimate"))+
  scale_color_manual(values = c("MLE Estimate" = "blue",
                                "MOM Estimate" = "red"))+
  theme_minimal()+
  xlab("Proportion of Deaths")+
  ylab("Density")


#TASK 8


simulation = tibble()

for(i in 1:1000){
  set.seed(7272 + i)
  
  current.sample = rbeta(266, 8, 950)
  
  mles.looped <- optim(par = c(0.5,0.5),
                 fn = llbeta,
                 data = current.sample,
                 neg=T)
  alpha.hat.mle.store <- mles.looped$par[1]
  beta.hat.mle.store <- mles.looped$par[2]
  
  
  moms.looped = nleqslv(x = c(0.5,0.5),
                 fn = mombeta,
                 data = current.sample) 
  
  alpha.hat.mom.store <- moms.looped$x[1]
  beta.hat.mom.store <- moms.looped$x[2]
  
  
  in.simulation = tibble(mle.alpha = alpha.hat.mle.store,
                         mle.beta = beta.hat.mle.store,
                         mom.alpha = alpha.hat.mom.store,
                         mom.beta = beta.hat.mom.store,
                         simulation = i)
  
  simulation = simulation |>
    bind_rows(in.simulation) |>
    unique()
}


mle.alpha.plot = ggplot(simulation)+
  geom_density(aes(x=mle.alpha))+
  theme_minimal()+
  xlab("MLE Alphas")+
  ylab("Density")+
  geom_hline(yintercept = 0)

mle.beta.plot = ggplot(simulation)+
  geom_density(aes(x=mle.beta))+
  theme_minimal()+
  xlab("MLE Betas")+
  ylab("Density")+
  geom_hline(yintercept = 0)

mom.alpha.plot = ggplot(simulation)+
  geom_density(aes(x=mom.alpha))+
  theme_minimal()+
  xlab("MOM Alphas")+
  ylab("Density")+
  geom_hline(yintercept = 0)

mom.beta.plot = ggplot(simulation)+
  geom_density(aes(x=mom.beta))+
  theme_minimal()+
  xlab("MOM Betas")+
  ylab("Density")+
  geom_hline(yintercept = 0)

(mle.alpha.plot + mle.beta.plot) / (mom.alpha.plot + mom.beta.plot)


a.theta = 8
b.theta = 950

names = c("alpha.mle", "alpha.mom", "beta.mle", "beta.mom")

bias = c(mean(simulation$mle.alpha) - a.theta,
mean(simulation$mom.alpha) - a.theta,
mean(simulation$mle.beta) - b.theta,
mean(simulation$mom.beta) - b.theta)


presicion = c(1/var(simulation$mle.alpha),
1/var(simulation$mom.alpha),
1/var(simulation$mle.beta),
1/var(simulation$mom.beta))

MSE = c(var(simulation$mle.alpha) + (mean(simulation$mle.alpha) - a.theta)^2,
var(simulation$mom.alpha) + (mean(simulation$mom.alpha) - a.theta)^2,
var(simulation$mle.beta) + (mean(simulation$mle.beta) - b.theta)^2,
var(simulation$mom.beta) + (mean(simulation$mom.beta) - b.theta)^2)


predictors = tibble(names,bias,presicion,MSE)

