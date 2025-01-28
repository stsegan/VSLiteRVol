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
rw <- read.rwl("C:/Users/User/Local Documents/R_Data/.rwl/UK/NET.rwl") 

# 1. Filter for site climate data
clim <- read_csv("C:/Users/User/Local Documents/R_Data/.csv & .xlsx/Climate/fagus_climate.csv")

clim <- clim %>% filter(site_id == "NET") %>% 
  select(month, year, tmp, pre)
lat <- 51.57
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

for(i in 1:length(m)){
  
vs_list[[i]] <- VSLite(syear = 1901, eyear = 2016, phi = lat, Te = tmp, 
                      Pr = pre, k = 2, m = m[i])

trw_list[[i]] <- t(as.data.frame(vs_list[[i]]$trw))
}

RespCur <- as.data.frame(vs_list[[i]]$gT)

RespCur <- RespCur %>% 
  mutate(Year = rownames(RespCur))

RC_longer <-  pivot_longer(RespCur, cols = -Year, names_to = "Series", values_to = "RC")

ggplot(RC_longer, aes(x = Year, y = RC, color = Series, group = Series)) +
  geom_line(stat = "smooth",method = "loess", alpha = 0.4, aes(group = Series), se = F, show.legend = F) +
  theme_cowplot(font_size = 10, rel_small = 0.75, rel_large = 1.25) +
  scale_x_discrete(
    breaks = seq(1900, 2040, by = 20)
  ) +
  labs(title = "Sigmoidal Ensemble Response Curve", x = "Time", y = "Response Curve")
