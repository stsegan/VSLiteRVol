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
library(grid)
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
library(progress)
}

#### 1. Filter for site climate data ####
clim <- read_csv("C:/Users/User/Local Documents/R_Data/.csv & .xlsx/Climate/fagus_climate.csv")
# view(fagus_meta)
clim <- clim %>% filter(site_id == "NET") %>% 
  select(month, year, tmp, pre)
lat <- 51.57 # found in fagus_meta
syear <- head(clim$year, n = 1)
eyear <- tail(clim$year, n = 1)



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

clim <- read_csv("C:/Users/User/Local Documents/R_Data/.csv & .xlsx/Climate/fagus_climate.csv")

clim <- clim %>% filter(site_id == "NET") %>% 
  select(month, year, tmp, pre)
clim <- clim %>% 
  mutate(num = c(1:length(year)))


# Trend
clim <- clim %>% 
  mutate(trend = 0.5 + (0.001 * (num * sin(num) ^ 2)))


# Plot
plot_trend <- ggplot(clim, aes(x = num, y = trend)) +
  geom_line(color = "dodgerblue4") +
  labs(
    title = "",
    x = "Monthly observation",
    y = "Trend Value"
  ) +
  theme_cowplot(font_size = 10, rel_small = 0.75, rel_large = 1.25)



# Induce trend on data
clim <- clim %>%
  mutate(trend_tmp = tmp * trend)
clim <- clim %>% 
  mutate(trend_pre = pre * trend)

# Format Climate Data

# trended temp
{
  trend_tmp <- clim %>% 
    select(month, trend_tmp) %>% 
    group_by(month) %>% 
    arrange(.by_group = TRUE) %>% 
    mutate(id = row_number()) %>% 
    pivot_wider(names_from = month, values_from = trend_tmp)
  
  
  trend_tmp <- trend_tmp %>% 
    mutate(year = c(syear:eyear)) %>% 
    relocate(year) %>% 
    select(-c(id, year))
  trend_tmp <- as.matrix(trend_tmp)
  trend_tmp <- t(trend_tmp)
  }

# trended prec
{
  trend_pre <- clim %>% 
    select(month, trend_pre) %>% 
    group_by(month) %>% 
    arrange(.by_group = TRUE) %>% 
    mutate(id = row_number()) %>% 
    pivot_wider(names_from = month, values_from = trend_pre)
  
  trend_pre <- trend_pre %>% 
    mutate(year = c(syear:eyear)) %>% 
    relocate(year) %>% 
    select(-c(id, year))
  trend_pre <- as.matrix(trend_pre)
  trend_pre <- t(trend_pre)
}



plot_tmp <- clim %>% 
  ggplot(aes(x = num, y = trend_tmp)) +
  geom_point(col = "darkorchid3", alpha = 0.4) + 
  geom_point(aes(x = num, y = tmp), colour = "goldenrod2", alpha = 0.4) +
  labs(
    title = "",
    x = "Month",
    y = "Mean Temperature (\u00B0C)"
  ) +
  theme_cowplot(font_size = 10, rel_small = 0.75, rel_large = 1.25)

plot_pre <- clim %>% 
  ggplot(aes(x = num)) +
  geom_point(aes(y = trend_pre, color = "Artificial"), alpha = 0.4) + 
  geom_point(aes(y = pre, color = "Observed"), alpha = 0.4) +
  labs(
    title = "",
    x = "Month",
    y = "Total Monthly Precipitation (mm)"
  ) +
  scale_color_manual(values = c("Artificial" = "darkorchid3", "Observed" = "goldenrod2"),
                     name = "Climate Variable") +
  theme_cowplot(font_size = 10, rel_small = 0.75, rel_large = 1.25)


lay <- rbind(c(1, 2),c(1, 2),
             c(3, 3))

obs_art_clim <- grid.arrange(plot_tmp, plot_pre, plot_trend, nrow = 2, layout_matrix = lay, 
             top = textGrob("Observed vs. Artificial Monthly Average Climate"
                            ,gp=gpar(fontsize = 16,font = 1)))

# ggsave(filename = "Observed vs. Artificial Monthly Average Climate.png", obs_art_clim, 
#       dpi = 500, device = "png", width = 10, height = 6,
#       path = "C:/Users/User/Local Documents/R_Data/Graphs/GARCH Paper")





#### 3. VSLite runs ####

# Formatting 
vs_list <- list()
trw_list <- list()
m <- seq(0.02, 10, by = 0.02)


# Run - for non-climate run, need to change ramp function in VSLite.R. 
# For climate run, change Te and Pr to trend_tmp and trend_pre, with Linear ramp.
for(i in 1:length(m)){
  
  vs_list[[i]] <- VSLite(syear = 1901, eyear = 2016, phi = lat, Te = tmp, 
                         Pr = pre, m = m[i], k = 2)
  
  trw_list[[i]] <- t(as.data.frame(vs_list[[i]]$trw))
}

# Formatting Output
trw_df <- as.data.frame(trw_list)
colnames(trw_df) <- m
trw_transformed <- apply(trw_df, 2, function(x) x + abs(min(x)) + 0.001)
trw_df <- as.data.frame(trw_transformed)


trw_df <- trw_df %>% 
  mutate(Year = c(syear:eyear))

#### 4. Response curve graph #### 

RespCur <- as.data.frame(vs_list[[i]]$gT) 

RespCur <- RespCur%>% 
  mutate(Month = rownames(RespCur))

RC_longer <-  pivot_longer(RespCur, cols = -Month, names_to = "Series", values_to = "RC")

vs_rc <- ggplot(RC_longer, aes(x = Month, y = RC, color = Series, group = Series)) +
  geom_line(stat = "smooth",method = "loess", alpha = 0.4, aes(group = Series), se = F, show.legend = F) +
  theme_cowplot(font_size = 10, rel_small = 0.75, rel_large = 1.25) +
  scale_x_discrete(breaks = seq(1, 12, 1)) +
  labs(title = "", x = "Month", y = "gT")



#### 5. Artificial Ring-Width plot ####

rw_longer <- pivot_longer(trw_df, cols = -Year, names_to = "Series", values_to = "RW")

vs_rwi <- ggplot(rw_longer, aes(x = Year, y = RW, color = Series, group = Series)) +
  geom_line(alpha = 0.1, show.legend = F) +
  theme_cowplot(font_size = 10, rel_small = 0.75, rel_large = 1.25) +
  scale_x_continuous(breaks = seq(1900, 2020, by = 20)) +
  labs(title = "", x = "Year", y = "Ring Width Index (RWI)")


#### 6. GARCH Run #### 

# Data formatting 
rownames(trw_df) <- trw_df$Year
trw_df <- trw_df %>% 
  select(-Year)

# GARCH function
detrend_gam <- function(rwl) {
  year <- as.numeric(rownames(rwl))
  tree_names <- colnames(rwl)
  out <- list()
  for (i in tree_names) {
    x <- rwl[, i]
    .data <- na.omit(data.frame(year = year, x = x))
    gam_model <- tryCatch(
      gam(x ~ s(year), data = .data),
      error = function(e) e)
    if (any(class(gam_model) == "error")) stop()
    .data$res <- .data$x
    .data$tree <- i
    out[[i]] <- .data
  }
  d <- bind_rows(out) %>% select(-x) %>% 
    arrange(year) %>% 
    pivot_wider(names_from = "tree",
                values_from = "res", values_fill = NA) %>% 
    data.frame()
  rownames(d) <- d$year
  d$year <- NULL
  d  
}
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
  
  pb <- txtProgressBar(min = 0, max = n_series, style = 3)
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
    # Update the progress bar
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
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
vs_garch <- do_garch(x = trw_df)
vs_vol <- vs_garch$volatility
vs_vol <- vs_vol[, !sapply(vs_vol, anyNA)]
vs_vol$Year <- rownames(vs_vol)
vs_vol <- vs_vol %>% 
  relocate(Year)

#### 7. Volatility Plot ####

vs_garch_long <- pivot_longer(vs_vol, cols = -Year, names_to = "Series", values_to = "Volatility")

vs_volplot <- ggplot(vs_garch_long, aes(x = Year, y = Volatility, color = Series, group = Series)) +
  geom_line(alpha = 0.2, show.legend = F) +
  theme_cowplot(font_size = 10, rel_small = 0.75, rel_large = 1.25) +
  scale_x_discrete(breaks = seq(1900, 2040, by = 20)) +
  labs(title = "", x = "Time", y = "Volatility")


lay <- rbind(c(1, 2),c(1, 2),
             c(3, 3))

vs_plots <- grid.arrange(vs_rwi, vs_volplot, vs_rc, nrow = 2, layout_matrix = lay, 
                            top = textGrob("Quadratic VSLite Plots (Lady Park, UK)"
                                           ,gp=gpar(fontsize = 20,font = 1)))

ggsave(filename = "Quadratic VSLite Volatility Plots, Lady Park UK.png", vs_plots, 
       dpi = 500, device = "png", width = 10, height = 6,
       path = "C:/Users/User/Local Documents/R_Data/Graphs/GARCH Paper")

#### 8. Analyse Time-Varying Volatility

# Formatting 
vs_vol <- vs_vol %>% 
  select(-Year)

trw_df_subset <- trw_df %>% 
  select(all_of(names(vs_vol)))

# Assuming your data frames are named 'ring_widths' and 'volatility_trends'
# Set the volatility threshold
volatility_threshold <- 1.5  # Adjust this value as needed

# Create a data frame with high volatility series
high_volatility_df <- trw_df_subset %>%
  select(which(apply(vs_vol, 2, function(col) any(col > volatility_threshold))))

# Create a data frame with non-high-volatility series
non_high_volatility_df <- trw_df_subset %>%
  select(-names(high_volatility_df))

# Function to calculate coefficient of variation
calc_cv <- function(x) {
  (sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)) * 100
}

# Calculate CV for each series in high volatility dataframe
high_vol_cv <- apply(high_volatility_df, 2, calc_cv)

# Calculate CV for each series in non-high-volatility dataframe
non_high_vol_cv <- apply(non_high_volatility_df, 2, calc_cv)

# Calculate average CV for each dataframe
avg_cv_high <- mean(high_vol_cv, na.rm = TRUE)
avg_cv_non_high <- mean(non_high_vol_cv, na.rm = TRUE)

# Print results
cat("Average CV for high volatility series:", round(avg_cv_high, 2), "%\n")
cat("Average CV for non-high-volatility series:", round(avg_cv_non_high, 2), "%\n")

# Calculate the difference in average CV
cv_difference <- avg_cv_high - avg_cv_non_high
cat("Difference in average CV:", round(cv_difference, 2), "percentage points\n")


# Find the lowest CV in the high volatility dataframe
highest_cv_high <- max(high_vol_cv, na.rm = TRUE)
highest_cv_high_series <- names(high_vol_cv)[which.max(high_vol_cv)]

# Find the highest CV in the non-high-volatility dataframe
lowest_cv_non_high <- min(non_high_vol_cv, na.rm = TRUE)
lowest_cv_non_high_series <- names(non_high_vol_cv)[which.min(non_high_vol_cv)]

# Print results
cat("Highest CV in high volatility series:", round(highest_cv_high, 2), "% (Series:", highest_cv_high_series, ")\n")
cat("Lowest CV in non-high-volatility series:", round(lowest_cv_non_high, 2), "% (Series:", lowest_cv_non_high_series, ")\n")

# Calculate the difference between these CVs
cv_difference <- highest_cv_high - lowest_cv_non_high
cat("Difference between highest high-volatility CV and lowest non-high-volatility CV:", round(cv_difference, 2), "percentage points\n")
