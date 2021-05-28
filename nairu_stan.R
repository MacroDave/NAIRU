library(ggthemes)
library(reshape2)
library(readabs)
library(dplyr)
library(ggplot2)
library(zoo)
library(rstan)
library(readrba)
library(lubridate)

options(mc.cores = parallel::detectCores())

#---------------------------------------------------------------------------------------------------------
		#Download Most Recent ABS and RBA Data
#---------------------------------------------------------------------------------------------------------
#Import Data from ABS Website
abs_5206 <- read_abs(series_id = c("A2304402X", "A2302915V"))
abs_6202 <- read_abs(series_id = c("A84423043C", "A84423047L"))
abs_6457 <- read_abs(series_id = c("A2298279F"))
abs_1364 <- read_abs(series_id = c("A2454521V", "A2454517C"))
rba_g3 <- read_rba(series_id = c("GBONYLD")) 
rba_g1 <- read_rba(series_id = c("GCPIOCPMTMQP")) 

#---------------------------------------------------------------------------------------------------------		
#CLEANUP ABS SPREADSHEETS
#---------------------------------------------------------------------------------------------------------
#5206.0 Australian National Accounts: National Income, Expenditure and Product
R_5206 <- abs_5206 %>%
  filter(series_id %in% c("A2304402X", "A2302915V")) %>% 
  mutate(date = zoo::as.yearqtr(date)) %>% 
  dplyr::select(date, series_id, value) 
R_5206 <- distinct(R_5206,date,series_id, .keep_all= TRUE)
R_5206 <- dcast(R_5206, date ~ series_id)
R_5206 <- R_5206 %>%
  mutate(NULC = A2302915V/A2304402X) %>%
  mutate(DLNULC = 100*(log(NULC)-log(lag(NULC,1)))) %>%
  select(date,DLNULC)

#6457.0 International Trade Price Indexes, Australia
R_6457 <- abs_6457 %>%
  filter(series_id %in% c("A2298279F")) %>%
  mutate(date = zoo::as.yearqtr(date)) %>%
  mutate(dl4pmcg = 100*(log(value)-log(lag(value,4)))) %>%
  dplyr::select(date, dl4pmcg)

#6202.0 Labour Force, Australia - Monthly
R_6202 <- abs_6202 %>%
  filter(series_id %in% c("A84423043C", "A84423047L")) %>%
  select(date, series_id, value)
R_6202 <- distinct(R_6202,date,series_id, .keep_all= TRUE)
R_6202 <- dcast(R_6202, date ~ series_id)
R_6202 <- R_6202 %>% group_by(date=floor_date(date, "quarter")) %>%
  summarize(A84423043C=mean(A84423043C), A84423047L=mean(A84423047L)) %>%
  mutate(date = zoo::as.yearqtr(date))
R_6202 <- R_6202 %>%
  mutate(LUR = 100*(1-A84423043C/A84423047L)) %>%
  select(date, LUR)

#Import RBA Data#
#Trimmed-Mean Inflation
R_g1 <- rba_g1 %>%
  filter(series_id %in% c("GCPIOCPMTMQP")) %>%
  mutate(date = zoo::as.yearqtr(date)) %>%
  rename(DLPTM = value) %>%
  select(date, DLPTM)

#Bond-market inflation expectations
R_g3 <- rba_g3 %>%
  filter(series_id %in% c("GBONYLD")) %>%
  mutate(date = zoo::as.yearqtr(date)) %>%
  mutate(pie_bondq = ((1+value/100)^(1/4)-1)*100) %>%
  select(date, pie_bondq)

NAIRU_data <- list(R_5206, R_6457, R_6202, R_g1, R_g3) %>%
  Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by="date"), .)

#Pick Sample
est_data <- NAIRU_data %>%
  filter(date>"1986q2" & date<"2020q1")

# Subset Data for Stan
stan_data <- est_data[,-1]

# run model ---------------------------------------------------------------
data_list <- list(T = nrow(stan_data),
                  J = ncol(stan_data),
                  Y = stan_data)

# Compile The Model
compiled_model <- stan_model(file = "NAIRU_est.stan")

sampled_model <- sampling(compiled_model, data = data_list, iter = 2000, control = list(max_treedepth = 15, adapt_delta = 0.99))

summarised_state <- as.data.frame(sampled_model) %>% 
  select(contains("NAIRU")) %>%
  melt() %>% 
  group_by(variable) %>% 
  summarise(median = median(value),
            lowera = quantile(value, 0.05),
            uppera = quantile(value, 0.95),
            lowerb = quantile(value, 0.15),
            upperb = quantile(value, 0.85)) %>%
  mutate(date = as.Date(est_data$date)) %>%
  mutate(date = zoo::as.yearqtr(date)) %>%
  mutate(LUR = est_data$LUR)

summarised_state %>% 
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin = lowera, ymax = uppera), fill = "orange", alpha = 0.3) +
  geom_line(aes(y = median), colour="red") +
  geom_line(aes(y = LUR)) +
  ggthemes::theme_economist() +
  ggtitle("NAIRU Estimate")

# Print estimated parameters from the model
print(sampled_model, pars = c("tau", "delta_pt", "beta_pt", "phi_pt", "gamma_pt", "lambda_pt", "alpha_pt", "eps_pt","delta_pu", "beta_pu", "gamma_pu", "lambda_pu", "eps_pu"))



