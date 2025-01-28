# Source Github code
{
  source("https://raw.githubusercontent.com/stsegan/VSLiteRVol/master/R/VSLite.R")
  source("https://raw.githubusercontent.com/stsegan/VSLiteRVol/master/R/compute.gE.R")
  source("https://raw.githubusercontent.com/stsegan/VSLiteRVol/master/R/dataset_doc.R")
  source("https://raw.githubusercontent.com/stsegan/VSLiteRVol/master/R/daylength.factor.from.lat.R")
  source("https://raw.githubusercontent.com/stsegan/VSLiteRVol/master/R/leakybucket.monthly.R")
  source("https://raw.githubusercontent.com/stsegan/VSLiteRVol/master/R/leakybucket.submonthly.R")
  source("https://raw.githubusercontent.com/stsegan/VSLiteRVol/master/R/std.ramp.R")
  source("https://raw.githubusercontent.com/stsegan/VSLiteRVol/master/R/test.R")
  source("https://raw.githubusercontent.com/stsegan/VSLiteRVol/master/R/param_est.r")
  source("https://raw.githubusercontent.com/stsegan/VSLiteRVol/master/R/sample_thresh_pars.r")
}

# Necessary Packages
{
library(cowplot)
library(wesanderson)
library(fagus)
library(gridExtra)
library(tidyverse)
library(strucchange)
library(dplR)
library(fGarch)
library(rugarch)
library(mgcv)
library(forecast)
library(aTSA)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
}

#### 0. Pick site ####
rw <- read.rwl("C:/Users/User/Local Documents/R_Data/.rwl/UK/NET.rwl") 



#### 1. Filter for site climate data ####
clim <- read_csv("C:/Users/User/Local Documents/R_Data/.csv & .xlsx/Climate/fagus_climate.csv")

clim <- clim %>% filter(site_id == "NET") %>% 
  select(month, year, tmp, pre)
lat <- 51.57 # found in fagus_meta
syear <- 1901
eyear <- 2016



#### 2. Format climate data ####

# temp
{
  tmp <- clim %>% 
    select(month, tmp) %>% 
    group_by(month) %>% 
    arrange(.by_group = TRUE) %>% 
    mutate(id = row_number()) %>% 
    pivot_wider(names_from = month, values_from = tmp)
  
  
  tmp <- tmp %>% 
    mutate(year = c(syear:eyear)) %>% 
    relocate(year) %>% 
    select(-c(id, year))
  tmp <- as.matrix(tmp)
  tmp <- t(tmp)
}

# prec
{
  pre <- clim %>% 
    select(month, pre) %>% 
    group_by(month) %>% 
    arrange(.by_group = TRUE) %>% 
    mutate(id = row_number()) %>% 
    pivot_wider(names_from = month, values_from = pre)
  
  pre <- pre %>% 
    mutate(year = c(syear:eyear)) %>% 
    relocate(year) %>% 
    select(-c(id, year))
  pre <- as.matrix(pre)
  pre <- t(pre)
}


#### 2.5 Alterations to climate data #### 
clim_reg <- read_csv2("C:/Users/User/Local Documents//R_Data/.csv & .xlsx/Climate/REG_Climate.csv")
clim_reg[,1]
clim_reg <- clim_reg %>% 
  mutate(num = c(1:1440))

# Option 1
clim_reg <- clim_reg %>% 
  mutate(trend = 0.5 + (0.001 * (num * sin(num) ^ 2)))

# Option 2
#out <- lm(tmax ~ num * sin(num), data = clim_reg)
#out_fv <- fitted.values(out)
#out_fv <- 0.1 * out_fv # 1 + (out_fv - mean(out_fv))
#clim_reg <- clim_reg %>% 
#  mutate(trend = out_fv)

plot(x = clim_reg$year, y = clim_reg$trend, 
     type = "p")

clim_reg <- clim_reg %>%
  mutate(trend_tmax = tmax * trend)
clim_reg <- clim_reg %>% 
  mutate(trend_prec = prec * trend)


clim_reg %>% 
  ggplot(aes(x = num, y = trend_tmax)) +
  geom_point(col = "darkblue", alpha = 0.4) + 
  geom_point(aes(x = num, y = tmax), colour = "lightblue3", alpha = 0.4) +
  geom_smooth(method = "loess", aes(y = trend_tmax), se = F, colour = "darkblue", lwd = 1.5, span = 2) +
  geom_smooth(method = "loess", aes(y = tmax), se = F, colour = "lightblue", lwd = 1.5, span = 2) +
  labs(
    title = "Observed vs. Artificial Monthly Max. Temperature",
    x = "Month",
    y = "Maximum Temperature"
  ) +
  theme_cowplot()

clim_reg %>% 
  ggplot(aes(x = num, y = trend_prec)) +
  geom_point(col = "darkblue", alpha = 0.4) + 
  geom_point(aes(x = num, y = prec), colour = "lightblue3", alpha = 0.4) +
  geom_smooth(method = "loess", aes(y = trend_prec), se = F, colour = "darkblue", lwd = 1.5, span = 2) +
  geom_smooth(method = "loess", aes(y = prec), se = F, colour = "lightblue", lwd = 1.5, span = 2) +
  labs(
    title = "Observed vs. Artificial Monthly Precipitation",
    x = "Month",
    y = "Monthly Precipitation (mm)"
  ) +
  theme_cowplot()




#### 3. VSLite runs ####

# Formatting 
vs_list <- list()
trw_list <- list()
k <- seq(0.1, 10, by = 0.1)

# Run
for(i in 1:length(k)){
  
  vs_list[[i]] <- VSLite(syear = 1901, eyear = 2016, phi = lat, Te = tmp, 
                         Pr = pre, k = k[i], m = 1)
  
  trw_list[[i]] <- t(as.data.frame(vs_list[[i]]$trw))
}

# Formatting Output
trw_df <- as.data.frame(trw_list)
colnames(trw_df) <- k
rownames(trw_df) <- c(syear:eyear)

#### 4. Response curve graph #### 

RespCur <- as.data.frame(vs_list[[i]]$gT) 

RespCur <- RespCur%>% 
  mutate(Year = rownames(RespCur))

RC_longer <-  pivot_longer(RespCur, cols = -Year, names_to = "Series", values_to = "RC")

ggplot(RC_longer, aes(x = Year, y = RC, color = Series, group = Series)) +
  geom_line(stat = "smooth",method = "loess", alpha = 0.4, aes(group = Series), se = F, show.legend = F) +
  theme_cowplot(font_size = 10, rel_small = 0.75, rel_large = 1.25) +
  scale_x_discrete(
    breaks = seq(1900, 2040, by = 20)
  ) +
  labs(title = "Sigmoidal Ensemble Response Curve", x = "Time", y = "Response Curve")

#### 5. Create .rwl of artificial ring-widths ####

write.tucson(rwl.df = trw_df, fname = "Nettlebed Sigmoid.rwl")
trw_rwl <- read.rwl("C:/Users/User/Local Documents/R_Data/.rwl/UK/Nettlebed Sigmoid.rwl")

#### 6. Artificial Ring-Width plot ####

rw_df <- as.data.frame(trw_rwl) 

rw_df <- rw_df %>% 
  mutate(Year = rownames(rw_df))

rw_longer <- pivot_longer(rw_df, cols = -Year, names_to = "Series", values_to = "RW")

ggplot(rw_longer, aes(x = Year, y = RW, color = Series, group = Series)) +
  geom_line(alpha = 0.20, show.legend = F) +
  theme_cowplot(font_size = 10, rel_small = 0.75, rel_large = 1.25) +
  scale_x_discrete(
    breaks = seq(1900, 2040, by = 20)
  ) +
  labs(title = "Sigmoidal Ensemble Ring Width Plot, Nettlebed UK", x = "Time", y = "Ring Width (mm)")

#### 7. GARCH Run #### 

# GARCH function
do_garch <- function(x, detrending = "gam", ...) {
  
  get_aic <- function(x) {
    if (any(class(x) == "error")) {
      return(NA_real_)
    }
    infocriteria(x)[1]
  }
  
  .names <- colnames(x)
  
  # step 1: detrending
  if (detrending == "gam") {
    xd <- detrend_gam(x)
  } else {
    xd <- detrend(x, method = "Spline", ...)
  }
  
  out <- list()
  out$rwi <- xd
  
  # step 2: ARIMA model
  n_series <- ncol(xd)
  years <- as.numeric(rownames(xd))
  arima_models <- list()
  garch_models <- list()
  arima_output <- list()
  volatility_output <- list()
  for (i in 1:n_series) {
    series <- xd[ ,i]
    .name <- .names[i]
    series_df <- data.frame(years, series) |> na.omit()
    start_year <- min(series_df$years)
    series_ts <- ts(series_df$series, start = start_year)  
    aa <- forecast::auto.arima(series_ts, stationary = FALSE)
    arma_ar <- aa$arma[1]
    arma_i <- aa$arma[2]
    arma_ma <- aa$arma[3]
    arma_cal <- arima(series_ts, c(arma_ar, arma_i, arma_ma))
    mlt <- aTSA::arch.test(arma_cal, FALSE)
    arima_output[[i]] <- data.frame(year = series_df$years,
                                    resid = residuals(arma_cal))
    colnames(arima_output[[i]])[2] <- .name
    arima_models[[i]] <- aa
    if (any(mlt[, 3] < 0.05)) { # PQ only
      garch_needed <- TRUE
    } else {
      garch_needed <- FALSE
    }
    
    # step 3: if necessary, fit GARCH to individual series
    if (garch_needed) {
      ps <- qs <- 1:3
      garch_coefs <- expand.grid(ps, qs)
      n <- nrow(garch_coefs)
      fit_ugarch <- function(i){
        tryCatch(ugarchfit(
          ugarchspec(variance.model = list(
            model = "sGARCH", garchOrder = c(garch_coefs[i, 1],
                                             garch_coefs[i, 2])),
            mean.model = list(armaOrder = c(arma_ar, arma_ma))),
          data = series_ts), error = function(e) e)
      }
      garch_fits <- lapply(1:n, fit_ugarch)
      aics <- sapply(garch_fits, get_aic)
      fit_classes <- sapply(garch_fits, class)
      fit_success <- !sapply(fit_classes, function(x) any(x == "error"))
      # if (!all(fit_success)) browser()
      if (any(fit_success)) {
        best_index <- which.min(aics)
        best_order <- garch_coefs[best_index, ]
        best_model <- garch_fits[[best_index]]
        vol <- best_model@fit$sigma
        volatility_output[[i]] <- data.frame(year = series_df$years,
                                             vol = vol)
        colnames(volatility_output[[i]])[2] <- .name
        garch_models[[i]] <- best_model  
      } else {
        volatility_output[[i]] <- data.frame(year = series_df$years,
                                             vol = NA_real_)
        colnames(volatility_output[[i]])[2] <- .name  
      }
      
    } else {
      volatility_output[[i]] <- data.frame(year = series_df$years,
                                           vol = NA_real_)
      colnames(volatility_output[[i]])[2] <- .name
    }
  }
  out$arima_models <- arima_models
  arima_output <- Reduce(merge, arima_output)
  rownames(arima_output) <- arima_output$year
  arima_output$year <- NULL
  out$arima <- arima_output
  volatility_output <- Reduce(merge, volatility_output)
  rownames(volatility_output) <- volatility_output$year
  volatility_output$year <- NULL
  out$volatility <- volatility_output
  out$garch_models <- garch_models
  return(out)
}

# Run
vs_garch <- do_garch(x = trw_rwl, detrending = "Spline", nyrs = 32)
vs_vol <- vs_garch$volatility
vs_vol <- vs_vol[, !sapply(vs_vol, anyNA)]
vs_vol$Year <- rownames(vs_vol)

#### 8. Volatility Plot ####

vs_garch_long <- pivot_longer(vs_vol, cols = -Year, names_to = "Series", values_to = "Volatility")

ggplot(vs_garch_long, aes(x = Year, y = Volatility, color = Series, group = Series)) +
  geom_line(alpha = 0.20, show.legend = F) +
  theme_cowplot(font_size = 10, rel_small = 0.75, rel_large = 1.25) +
  scale_x_discrete(
    breaks = seq(1900, 2040, by = 20)
  ) +
  labs(title = "Sigmoidal VSLite Volatility Plot, Nettlebed UK", x = "Time", y = "Volatility")

