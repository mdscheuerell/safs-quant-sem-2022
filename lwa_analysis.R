
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
  select(-c("Year", "Month", "Temp", "TP", "pH", "Unicells", "Other.algae", "Neomysis")) %>%
  t() %>%
  zscore(mean.only = TRUE)

nn <- nrow(plank)

CC <- matrix(list(0), nrow = ncol(plank), ncol = 2)
CC[,1] <- colnames(plank)
CC[c(1, 2, 3, 4),2] <- colnames(plank)[c(1, 2, 3, 4)]
CC

## rownames(plank)
# [1,] "1"  "Cryptomonas"            
# [2,] "2"  "Diatoms"                
# [3,] "3"  "Greens"                 
# [4,] "4"  "Bluegreens"             
# [5,] "5"  "Conochilus"             
# [6,] "6"  "Cyclops"                
# [7,] "7"  "Daphnia"                
# [8,] "8"  "Diaptomus"              
# [9,] "9"  "Epischura"              
# [10,] "10" "Leptodora"              
# [11,] "11" "Non.daphnid.cladocerans"
# [12,] "12" "Non.colonial.rotifers" 

BB <- matrix(list(0), nn, nn)
for(rr in 1:nn) {
  for(cc in 1:nn) {
    BB[rr, cc] <- paste(rr, cc, sep = ",")
  }
}
BB[seq(4), seq(9,10)] <- BB[seq(9,10), seq(4)] <- 0
BB

model_list <- list(
  B = BB,
  U = "zero",
  Q = "diagonal and unequal",
  C = "unconstrained",
  c = t(covars),
  Z = "identity",
  A = "zero",
  R = "diagonal and equal"
)

control_list <- list(
  maxit = 5000
)


lwa_fit <- MARSS(plank, model = model_list, control = control_list, method = "BFGS")





