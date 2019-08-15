library(ggthemes); library(stringr);library(reshape2);library(readabs);library(dplyr); library(ggplot2); library(zoo); library(raustats); library(Quandl); library(rstan)
options(mc.cores = parallel::detectCores())

library(httr)
httr::set_config(config(ssl_verifypeer = FALSE))

#consumer price index - qtly original
CPI <- read_abs(series_id = "A2325846C") %>%
  mutate(date = zoo::as.yearmon(date)) %>% 
  mutate(CPI = 100*log(value/lag(value, 1))) %>% 
  dplyr::select(date, CPI)

#real gdp millions chained dollars - qtly seasonally adj
GDPA <- read_abs(series_id = "A2304402X") %>%
  mutate(date = zoo::as.yearmon(date)) %>%   
  mutate(GDPA = 100*log(value/lag(value, 1))) %>% 
  dplyr::select(date, GDPA)

#Trimmed Mean Inflation - Qtly Growth
CPIT <- Quandl("RBA/G01_GCPIOCPMTMQP") %>% 
  mutate(Date = zoo::as.yearmon(Date)) %>% 
  arrange(Date) %>%
  rename(date = Date) %>%
  rename(value = "Quarterly trimmed mean inflation. Units: Per cent; Series ID: GCPIOCPMTMQP") %>%
  mutate(CPIT = value) %>%
  dplyr::select(date,CPIT)

#Inflation Expectations - Break-even 10-year inflation rate (converted to qtly)
CPIE <- Quandl("RBA/G03_GBONYLD") %>% 
  mutate(Date = zoo::as.yearmon(Date)) %>% 
  arrange(Date) %>%
  rename(date = Date) %>%
  rename(value = "Break-even 10-year inflation rate. Units: Per cent ; Series ID: GBONYLD") %>%
  mutate(CPIE = ((1+value/100)^(1/4)-1)*100) %>%
  dplyr::select(date, CPIE)

#Nominal Unit Labor Costs - Qtly Seasonally Adjusted
ULC <- read_abs(series_id = "A2433068T") %>%
  mutate(date = zoo::as.yearmon(date)) %>%   
  mutate(ULC = 100*log(value/lag(value, 1))) %>% 
  dplyr::select(date, ULC)

#Unemployment Rate Persons - Monthly Seasonally Adjusted
UNR <- read_abs(series_id = "A84423050A") %>%
  mutate(date = as.Date(date)) %>% 
  mutate(date = zoo::as.yearmon(date)) %>%   
  mutate(UNR = (value+lag(value,1)+lag(value,2))/3) %>%
  dplyr::select(date, UNR)

#Import Deflator - Qtly Seasonally Adjusted
IMPD <- read_abs(series_id = "A2303729J") %>%
  mutate(date = zoo::as.yearmon(date)) %>% 
  mutate(IMPD = 100*log(value/lag(value, 4))) %>% 
  dplyr::select(date, IMPD)

# Join series together, and set missing values to -999 (Stan doesn't like NAs)
full_data <- list(CPI, GDPA, ULC, IMPD, CPIT, CPIE, UNR) %>%
  Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by="date"), .)

#Pick Sample
est_data <- subset(full_data, date>"June 1986")

# Subset Data for Stan
stan_data <- est_data[-nrow(est_data),-1]

# run model ---------------------------------------------------------------
data_list <- list(T = nrow(stan_data),
                  J = ncol(stan_data),
                  Y = stan_data)

# Compile The Model
compiled_model <- stan_model(file = "NAIRU_calibrate.stan")

sampled_model <- sampling(compiled_model, data = data_list, iter = 1000, cores = 4)

summarised_state <- as.data.frame(sampled_model) %>% 
  select(contains("NAIRU")) %>%
  melt() %>% 
  group_by(variable) %>% 
  summarise(median = median(value),
            lowera = quantile(value, 0.05),
            uppera = quantile(value, 0.95),
            lowerb = quantile(value, 0.15),
            upperb = quantile(value, 0.85)) 

summarised_state$date<-seq(as.Date("1986/12/1"), as.Date("2019/06/1"), by = "quarter")
summarised_state <- summarised_state %>%
  mutate(date = zoo::as.yearmon(date))

UR <- full_data[,-(2:7)]
summarised_state <- left_join(summarised_state,UR,by="date")

summarised_state %>% 
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin = lowera, ymax = uppera), fill = "orange", alpha = 0.3) +
  geom_line(aes(y = median), colour="red") +
  geom_line(aes(y = UNR)) +
  ggthemes::theme_economist() +
  ggtitle("NAIRU Estimate")

# Print estimated parameters from the model
#print(sampled_model, pars = c("tau", "delta_pt", "beta_pt", "phi_pt", "gamma_pt", "lambda_pt", "alpha_pt", "eps_pt","delta_pu", "beta_pu", "gamma_pu", "lambda_pu", "eps_pu"))



