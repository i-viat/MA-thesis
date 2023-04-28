library(tidyverse)
library(rjags)
library(ggpubr)
library(glue) 



data = read.csv("C://Users/User/Desktop/MA_Thesis/OSF_upload/data/dataset_for_analysis.csv") # specify path to the dataset



### Data transformation

# Replace NA with 0 in protests
data["pr_quan_avg"][is.na(data["pr_quan_avg"])] <- 0

# Log the value of social expenditure per capita
data$soc_sp_pc_log = log(data$soc_sp_pc_ria_avg)
data$soc_sp_pc_19_log = log(data$soc_sp_pc_ria_19)
data$soc_sp_pc_20_log = log(data$soc_sp_pc_ria_20)
data$soc_sp_pc_21_log = log(data$soc_sp_pc_ria_21)

# Skip Leningard oblast as an outlier
data_84 = data %>% filter(region != 'Leningrad')



## Dispersion test

dispersion_test <- function(x) 
{
  res <- 1-2 * abs((1 - pchisq((sum((x - mean(x))^2)/mean(x)), length(x) - 1))-0.5)
  
  cat("Dispersion test of count data:\n",
      length(x), " data points.\n",
      "Mean: ",mean(x),"\n",
      "Variance: ",var(x),"\n",
      "Probability of being drawn from Poisson distribution: ", 
      round(res, 3),"\n", sep = "")
  
  invisible(res)
}

dispersion_test(na.omit(data_84$pb_count_sum))



### Defining model templates

run_jags = function(data_raw, d, string_mod_final, inter_cov_1, inter_cov_2){
  params = c("int")
  for (i in seq_len(nrow(d))){
    params <- c(params, d[i, 1])  
  }
  params <- c(params, "r")
  if(inter_cov_1 != FALSE & inter_cov_2 != FALSE){
    params <- c(params, "b_int")
  }
  set.seed(1)
  print(params)
  mod = jags.model(textConnection(string_mod_final), data=as.list(na.omit(data_raw)), n.chains=3)
  update(mod, 1e3)
  mod_sim = coda.samples(model=mod, variable.names=params, n.iter=3e4)
  print(summary(mod_sim))
  print(dic.samples(mod, n.iter=1e3))
  return(mod_sim)
}

neg_bin_likelihood = function(d, dv, inter_cov_1, inter_cov_2){
  string_mod = glue(" model {{
    for (i in 1:length({dv})) {{
        {dv}[i] ~ dnegbin(p[i],r)
        p[i] <- r/(r+lambda[i]) 
        log(lambda[i]) = int")
  
  for (i in seq_len(nrow(d))){
    string_mod = paste(string_mod,' + ', d[i, 1], '*', d[i, 2], '[i]', sep='')
  }
  
  if(inter_cov_1 != FALSE){
    string_mod = paste(string_mod,' + ', 'b_int*', inter_cov_1, '[i]*', inter_cov_2, '[i]', sep='')
  }
  return(string_mod)
}

zero_infl_likelihood = function(d, dv, inter_cov_1, inter_cov_2){
  string_mod = glue(" model {{
     for (i in 1:length({dv})) {{
        {dv}[i] ~ dnegbin(p[i], r)
        p[i] <- r/(r+(1-zero[i])*lambda.count[i]) - 1e-10*zero[i]
        lambda.count[i] = exp(int")
  for (i in seq_len(nrow(d))){
    string_mod = paste(string_mod,' + ', d[i, 1], '*', d[i, 2], '[i]', sep='')
  }
  
  if(inter_cov_1 != FALSE){
    string_mod = paste(string_mod,' + ', 'b_int*', inter_cov_1, '[i]*', inter_cov_2, '[i]', sep='')
  }
  string_mod = paste(string_mod, ')')
  
  # Zero-Inflation
  zero_inf = "zero[i] ~ dbern(pi[i]) pi[i] <- ilogit(int"
  
  for (i in seq_len(nrow(d))){
    zero_inf = paste(zero_inf,' + ', d[i, 1], '*', d[i, 2], '[i]', sep='')
  }
  if(inter_cov_1 != FALSE){
    zero_inf = paste(zero_inf,' + ', 'b_int*', inter_cov_1, '[i]*', inter_cov_2, '[i]', sep='')
  }
  zero_inf = paste(zero_inf, ')', sep='')
  string_mod = paste(string_mod, zero_inf)
  return(string_mod)
}

priors = function(d, string_mod, inter_cov_1, inter_cov_2){
  string_mod_full = paste(string_mod, '} int ~ dnorm(0.0, 1.0/10e1)')
  
  for (i in seq_len(nrow(d))){
    string_mod_full = paste(string_mod_full, d[i, 1], '~ dnorm(0.0, 1.0/10e1)')
  }
  
  if(inter_cov_1 != FALSE){
    string_mod_full = paste(string_mod_full, 'b_int ~ dnorm(0.0, 1.0/10e1)')
  }
  string_mod_final = paste(string_mod_full, ' r ~ dunif(0, 50) }')
  return(string_mod_final)
}

run_model = function(dv, vector_coef, vector_cov, dataset, inter_cov_1=FALSE, inter_cov_2=FALSE, zero_infl=FALSE){
  
  jags_cov = c(vector_cov, dv)
  
  data_raw = dataset %>% dplyr::select(all_of(jags_cov))
  
  d <- data.frame(coef = vector_coef,
                  cov = vector_cov)
  
  if(zero_infl==FALSE){
    string_mod = neg_bin_likelihood(d, dv, inter_cov_1, inter_cov_2)
  } else {
    string_mod = zero_infl_likelihood(d, dv, inter_cov_1, inter_cov_2)
  }
  
  string_mod_final = priors(d, string_mod, inter_cov_1, inter_cov_2)
  
  return(run_jags(data_raw, d, string_mod_final, inter_cov_1, inter_cov_2))
}



### Main analysis

## Model 1

mod_1 = run_model('pb_count_sum', c("b_ins", "b_elt", "b_pro"), 
                  c('inst_qual_17', 'elites', 'pr_quan_avg'), data_84)

## Model 2

mod_1 = run_model('pb_count_sum', c("b_ins", "b_elt", "b_pro", "b_eff", "b_hdi", "b_tci"),  
                  c('inst_qual_17', 'elites', 'pr_quan_avg', 'total_eff', 'hdi_16', 'tci_avg'), data_84)

## Model 3

mod_1 = run_model('pb_count_sum', c("b_bal", "b_trf", "b_soc"),
                  c('budg_bal_avg',  'cbf_budg_tr_share_avg', 'soc_sp_pc_log'), data_84)

## Model 4

mod_1 = run_model('pb_count_sum', c("b_bal", "b_trf", "b_soc", "b_eff", "b_hdi", "b_tci"),  
                  c('budg_bal_avg', 'cbf_budg_tr_share_avg', 'soc_sp_pc_log', 'total_eff', 'hdi_16', 'tci_avg'))

## Model 5

mod_1 = run_model('pb_count_sum', c("b_ins", "b_elt", "b_pro", "b_bal", "b_trf", "b_soc"), 
                  c('inst_qual_17', 'elites', 'pr_quan_avg', 'budg_bal_avg', 'cbf_budg_tr_share_avg', 'soc_sp_pc_log'), data_84)

## Model 6

full_model_coef = c("b_ins", "b_elt", "b_pro", "b_bal", "b_trf", "b_soc", "b_eff", "b_hdi", "b_tci")
full_model_cov = c('inst_qual_17', 'elites', 'pr_quan_avg', 'budg_bal_avg', 'cbf_budg_tr_share_avg', 
                   'soc_sp_pc_log', 'total_eff', 'hdi_16', 'tci_avg')

mod_6 = run_model('pb_count_sum', full_model_coef,  full_model_cov, data_84)

## Model 7

mod_7 = run_model('pb_count_sum', full_model_coef,  full_model_cov, data_84, 
                  inter_cov_1='inst_qual_17', inter_cov_2='elites')

## Model 8

mod_8 = run_model('pb_count_sum', full_model_coef,  full_model_cov, data_84, 
                  inter_cov_1='inst_qual_17', inter_cov_2='pr_quan_avg')

## Model 9

mod_9 = run_model('pb_count_sum', full_model_coef,  full_model_cov, data_84, 
                  inter_cov_1='elites', inter_cov_2='pr_quan_avg')

## Model 10

mod_10 = run_model('pb_count_sum', full_model_coef,  full_model_cov, data_84, 
                   inter_cov_1='budg_bal_avg', inter_cov_2='cbf_budg_tr_share_avg')

## Model 11

mod_11 = run_model('pb_count_sum', full_model_coef,  full_model_cov, data_84,
                   inter_cov_1='budg_bal_avg', inter_cov_2='soc_sp_pc_log')

## Model 12

mod_12 = run_model('pb_count_sum', full_model_coef,  full_model_cov, data_84,
                   inter_cov_1='cbf_budg_tr_share_avg', inter_cov_2='soc_sp_pc_log')

## Model 13

mod_13 = run_model('pb_count_sum', full_model_coef,  full_model_cov, data_84,
                   inter_cov_1='inst_qual_17', inter_cov_2='budg_bal_avg')

## Model 14

mod_14 = run_model('pb_count_sum', full_model_coef,  full_model_cov, data_84,
                   inter_cov_1='inst_qual_17', inter_cov_2='cbf_budg_tr_share_avg')

## Model 15

mod_15 = run_model('pb_count_sum', full_model_coef,  full_model_cov, data_84,
                   inter_cov_1='inst_qual_17', inter_cov_2='soc_sp_pc_log')

## Model 16

mod_16 = run_model('pb_count_sum', full_model_coef,  full_model_cov, data_84, 
                   inter_cov_1='elites', inter_cov_2='budg_bal_avg')

## Model 17

mod_17 = run_model('pb_count_sum', full_model_coef,  full_model_cov, data_84,
                   inter_cov_1='elites', inter_cov_2='cbf_budg_tr_share_avg')

## Model 18

mod_18 = run_model('pb_count_sum', full_model_coef,  full_model_cov, data_84, 
                   inter_cov_1='elites', inter_cov_2='soc_sp_pc_log')

## Model 19

mod_19 = run_model('pb_count_sum', full_model_coef,  full_model_cov, data_84,
                   inter_cov_1='pr_quan_avg', inter_cov_2='budg_bal_avg')

## Model 20

mod_20 = run_model('pb_count_sum', full_model_coef,  full_model_cov, data_84,
                   inter_cov_1='pr_quan_avg', inter_cov_2='cbf_budg_tr_share_avg')

## Model 21

mod_21 = run_model('pb_count_sum', full_model_coef,  full_model_cov, data_84, 
                   inter_cov_1='pr_quan_avg', inter_cov_2='soc_sp_pc_log')



### Plotting marginal effects - model 17

int.mcmc = as.mcmc(do.call(rbind, mod_17))
int.mcmc.mat <- as.matrix(int.mcmc)
int.mcmc.dat <- as.data.frame(int.mcmc.mat)

# Simulate the range of the moderating variable:
x2.sim <- seq(min(data_84$cbf_budg_tr_share_avg), max(data_84$cbf_budg_tr_share_avg), by = 0.1) 

# Calculate conditional effect of X1 across the range of X2, first using the Bayesian estimates:
int.sim <- matrix(rep(NA, 
                      nrow(int.mcmc.dat)*length(x2.sim)), 
                  nrow = nrow(int.mcmc.dat))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- int.mcmc.dat$b_elt + int.mcmc.dat$b_int * x2.sim[i]
}


bayes.c.eff.mean <- apply(int.sim, 2, mean)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.025)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.975)))


# Combine both estimates into one dataframe for plotting:
plot.dat <- data.frame(x2.sim, 
                       bayes.c.eff.mean, 
                       bayes.c.eff.lower, 
                       bayes.c.eff.upper)


ggplot(data=plot.dat, aes(x=x2.sim, y=bayes.c.eff.mean)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), alpha = 0.1) +
  xlab("Average Share of Federal Fiscal Transfers Received in the Regional Budget") +
  ylab("Posterior Estimate for Regional Elite Openness")



### Distributions visualization (Appendix 1)

plot_hist = function(variable, name){
  data = data %>% filter(data[,variable] < 500) # only for projects cost
  ggplot(data = data, aes(x=data[,variable])) + 
    geom_histogram(color="black", fill="white") + 
    xlab(paste(glue('Project Count by Region\n({name})'))) + 
    ylab('') +
    geom_vline(aes(xintercept=median(na.omit(data[,variable]))),
               color="red", linetype="dashed", size=1)
}  

p1 = plot_hist('pb_count_19', '2019')
p2 = plot_hist('pb_count_20', '2020')
p3 = plot_hist('pb_count_21', '2021')
p4 = plot_hist('pb_count_mean', 'average 2019–2021')
p5 = plot_hist('pb_count_sum', 'total 2019–2021')

ggarrange(p1, p2, p3, p4, p5, ncol = 3, nrow = 2)




### Diagnostics of models 6 and 17 (Appendix 2)

plot(mod_6)
gelman.diag(mod_6)
autocorr.diag(mod_6)
autocorr.plot(mod_6)

plot(mod_17)
gelman.diag(mod_17)
autocorr.diag(mod_17)
autocorr.plot(mod_17)



### Annual cross-sections - zero-inflated negative binomial regressions (Appendix 3)

## 2019

data_19 = data %>% filter(region != 'Kirov') # remove an outlier

ggplot(data = data_19, aes(x = pb_count_19)) +
  geom_histogram()

dispersion_test(na.omit(data_19$pb_count_19))


mod_infl = run_model('pb_count_19', c("b_ins", "b_elt", "b_pro", "b_bal", "b_trf", "b_soc"), 
                     c('inst_qual_17', 'elites', 'pr_quan_avg', "budg_bal_19", "cbf_budg_tr_share_19", "soc_sp_pc_19_log"),
                     data_19, zero_infl=TRUE)

full_model_cov_19 = c('inst_qual_17', 'elites', 'pr_quan_avg', 'budg_bal_19', 'cbf_budg_tr_share_19',
                      'soc_sp_pc_19_log', 'total_eff', 'hdi_16', 'tci_19')

mod_infl = run_model('pb_count_19', full_model_coef,  full_model_cov_19,
                     data_19, zero_infl=TRUE)

mod_infl = run_model('pb_count_19', full_model_coef,  full_model_cov_19, inter_cov_1='elites',
                     inter_cov_2='cbf_budg_tr_share_19', data_19, zero_infl=TRUE)



## 2020

data_20 = data %>% filter(region != 'Leningrad') # remove an outlier

ggplot(data = data_20, aes(x = pb_count_20)) +
  geom_histogram()

dispersion_test(na.omit(data_20$pb_count_20))

mod_infl = run_model('pb_count_20', c("b_ins", "b_elt", "b_pro", "b_bal", "b_trf", "b_soc"), 
                     c('inst_qual_17', 'elites', 'pr_quan_avg', "budg_bal_20", "cbf_budg_tr_share_20", "soc_sp_pc_20_log"),
                     data_20, zero_infl=TRUE)

full_model_cov_20 = c('inst_qual_17', 'elites', 'pr_quan_avg', 'budg_bal_20', 'cbf_budg_tr_share_20', 'soc_sp_pc_20_log', 'total_eff', 'hdi_16', 'tci_20')

mod_infl = run_model('pb_count_20', full_model_coef,  full_model_cov_20,
                     data_20, zero_infl=TRUE)

mod_infl = run_model('pb_count_20', full_model_coef,  full_model_cov_20, inter_cov_1='elites',
                     inter_cov_2='cbf_budg_tr_share_20', data_20, zero_infl=TRUE)



## 2021

data_21 = data %>% filter(region != 'Leningrad' & region != 'Krasnodar' & region != 'Udmurtia') # remove outliers

ggplot(data = data_21, aes(x = pb_count_21)) +
  geom_histogram()

dispersion_test(na.omit(data_21$pb_count_21))


mod_infl = run_model('pb_count_21', c("b_ins", "b_elt", "b_pro", "b_bal", "b_trf", "b_soc"), 
                     c('inst_qual_17', 'elites', 'pr_quan_avg', "budg_bal_21", "cbf_budg_tr_share_21", "soc_sp_pc_21_log"),
                     data_21, zero_infl=TRUE)

full_model_cov_21 = c('inst_qual_17', 'elites', 'pr_quan_avg', 'budg_bal_21', 'cbf_budg_tr_share_21', 'soc_sp_pc_21_log', 'total_eff', 'hdi_16', 'tci_21')

mod_infl = run_model('pb_count_21', full_model_coef,  full_model_cov_21,
                     data_21, zero_infl=TRUE)

mod_infl = run_model('pb_count_21', full_model_coef,  full_model_cov_21, inter_cov_1='elites',
                     inter_cov_2='cbf_budg_tr_share_21', data_21, zero_infl=TRUE)
