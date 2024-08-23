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

# 0. Pick site
view(fagus_meta)
rw <- read.rwl("C:/Users/User/Local Documents/R_Data/.rwl/UK/Lady Park.rwl") 

# 1. Filter for site climate data
clim <- read_csv("C:/Users/User/Local Documents/R_Data/.csv & .xlsx/Climate/fagus_climate.csv")
clim <- clim %>% filter(site_id == "UK10") %>% 
  select(month, year, tmp, pre)
lat <- 51.83
syear <- 1901
eyear <- 2016

# Format climate data
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

# VSLite runs
vs_list <- list()
trw_list <- list()
k <- seq(0.1, 10, by = 0.1)

for(i in 1:length(k)){
  
vs_list[[i]] <- VSLite(syear = 1901, eyear = 2016, phi = lat, Te = tmp, 
                      Pr = pre, k = k[i])

trw_list[[i]] <- t(as.data.frame(vs_list[[i]]$trw))
}

trw_df <- as.data.frame(trw_list)
colnames(trw_df) <- k
rownames(trw_df) <- c(syear:eyear)



# Create rwl. file for VSLite runs
write.tucson(rwl.df = trw_df, fname = "Lady Park VSLite linear.rwl")
trw_rwl <- read.rwl("C:/Users/User/Local Documents/R_Data/.rwl/Lady Park VSLite linear.rwl")

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
    if (any(mlt[, 3] < 0.05)) { # PQ only ATM
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

# GARCH run
vs_garch <- do_garch(x = trw_rwl, detrending = "Spline", nyrs = 32)
vs_vol <- vs_garch$volatility
vs_vol <- vs_vol[, !sapply(vs_vol, anyNA)]
vs_vol$Year <- rownames(vs_vol)

# vs_vol <- vs_vol %>% 
#  select(-'66')

vs_garch_long <- pivot_longer(vs_vol, cols = -Year, names_to = "Series", values_to = "Volatility")

ggplot(vs_garch_long, aes(x = Year, y = Volatility, color = Series, group = Series)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Linear VSLite Volatility Plot", x = "Time", y = "Volatility")



# Need to: 
# NEW 
# 1. Figure out why lines aren't showing on ggplot. 
# 2. Once sorted, run for quadratic and linear ramps. 
# 3. Need to try and find a way of how to present this data. 
# Maybe something we talk about with Christian.



# 1. Decide if I want to detrend outputted VSLite RW or not. 
# 2. Figure out how to loop GARCH model application.
# 3. Look at RW values. Are they okay or do they need to be scaled?
# 4. Figure out what to do with linear ramp. Can't scale in a similar way to 
# quadratic and sigmoid as endpoint will change. Is this a problem? Should this
# be incorporated to sigmoid/quadratic? 
  





