
library(tidyverse)
library(MARSS)


#### function definitions ####

BB <- matrix(c(0.8, 0.2, -0.1, 0.6), 2, 2)

CC <- matrix(c(0.5, 0.4, -0.2, -0.3), 2, 2)

LL <- function(B, C) {
  ll <- matrix(NA, nrow(B), ncol(C))
  for(i in 1:ncol(C)) {
    for(j in 1:nrow(B)) {
      numer <- B
      numer[,j] <- C[,i]
      ll[j,i] <- det(numer) / det(B)
    }
  }
  return(ll)
}

LL(BB, CC)

data(package = "MARSS")

head(ivesDataByWeek)

lwa_1 <- lakeWAplanktonRaw %>%
  as_tibble() %>%
  filter(Year <= 1975) %>%
  apply(MARGIN = 2, quantile, na.rm = TRUE)



covars_raw <- lakeWAplanktonRaw %>%
  as_tibble() %>%
  filter(Year >= 1976) %>%
  select(c("Temp", "TP")) %>%
  ts(start = c(1976, 1), frequency = 12)

plot.ts(covars_raw, las = 1,
        ylab = "")

plot(covars_raw[,"Temp"], covars_raw[,"TP"])

R2 <- lm(Temp ~ TP, data = covars_raw) %>%
  summary
  

VIF <- 1 / (1 - R2$r.squared)


lwa <- lakeWAplanktonTrans %>%
  as_tibble() %>%
  filter(Year >= 1976)

covars <- lwa %>%
  select(c("Temp", "TP"))

plank <- lwa %>%
  select(-c("Year", "Month", "Temp", "TP", "pH"))

