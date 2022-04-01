
library(tidyverse)
library(MARSS)


#### function definitions ####

BB <- matrix(c(0.8, 0.2, -0.1, 0.6), 2, 2)

CC <- matrix(c(0.5, 0.4, -0.2, -0.3), 2, 2)

dpiB <- function(B) {
  2 * det(B) * t(solve(B))
}


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




plot.ts(covars_raw, las = 1,
        ylab = "")

# plot(covars_raw[,"Temp"], covars_raw[,"TP"])
# 
# R2 <- lm(Temp ~ TP, data = covars_raw) %>%
#   summary
#   
# 
# VIF <- 1 / (1 - R2$r.squared)


lwa <- lakeWAplanktonRaw %>%
  as_tibble() %>%
  filter(Year >= 1977) %>%
  mutate(lg_phyto = log(Diatoms + Greens),
         sm_phyto = log(Unicells + Cryptomonas),
         non_Daphnia = log(Conochilus + Cyclops + Diaptomus +
           Non.daphnid.cladocerans + Non.colonial.rotifers),
         Daphnia = log(Daphnia + 0.001)) %>%
  select(c("Year", "Month", "Temp", "TP",
           "lg_phyto", "sm_phyto", "Daphnia", "non_Daphnia"))

p1 <- seq(1977, 1982)
p2 <- seq(1989, 1994)


colnames(lwa)

covars_raw <- lakeWAplanktonRaw %>%
  as_tibble() %>%
  filter(Year >= 1977) %>%
  select(c("Temp", "TP")) %>%
  ts(start = c(1977, 1), frequency = 12)

colnames(covars_raw) <- c("Temp (C)", "P (mg/L)")

png("covars.png", width = 8, height = 5, units = "in",
    res = 300)
par(mai = c(0.9, 1, 0.1, 1),
    omi = c(0, 0, 0, 0))
plot.ts(covars_raw, lwd = 2,
        cex.axis = 1.2, cex.lab = 1.5,
        main = "", yax.flip = TRUE, col = c("dodgerblue"))
dev.off()


plank_raw <- lwa %>%
  select(-c("Year", "Month", "Temp", "TP"))  %>%
  ts(start = c(1977, 1), frequency = 12)

colnames(plank_raw) <- c("Lg phyto", "Sm phyto", "Daphnia", "non-Daphnia")

png("plank.png", width = 8, height = 6, units = "in",
    res = 300)
par(mai = c(0.9, 1, 0.1, 1),
    omi = c(0, 0, 0, 0))
plot.ts(plank_raw, lwd = 2,
        cex.axis = 1.5, cex.lab = 1.4,
        main = "", yax.flip = TRUE, col = c("dodgerblue"))
dev.off()


#### for models ####

plank <- lwa %>%
  select(-c("Month", "Temp", "TP"))
  
plank_1 <- plank %>%
  filter(Year >= 1977 & Year <= 1982) %>%
  select(-c("Year")) %>%
  as.matrix() %>%
  zscore() %>%
  t()

plank_2 <- plank %>%
  filter(Year >= 1989 & Year <= 1994) %>%
  select(-c("Year")) %>%
  as.matrix() %>%
  zscore() %>%
  t()

nn <- nrow(plank_1)


covars <- lakeWAplanktonTrans %>%
  as_tibble() %>%
  filter(Year >= 1977) %>%
  select(c("Year", "Temp", "TP"))

cov_1 <- covars %>%
  filter(Year >= 1977 & Year <= 1982) %>%
  select(c("Temp", "TP")) %>%
  as.matrix() %>%
  t()

cov_2 <- covars %>%
  filter(Year >= 1989 & Year <= 1994) %>%
  select(c("Temp", "TP")) %>%
  as.matrix() %>%
  t()


CC <- matrix(list(0), nrow = nn, ncol = 2)
CC[,1] <- paste0("T_", rownames(plank_1))
CC[c(1, 2), 2] <- paste0("P_", rownames(plank_1)[c(1, 2)])
CC

# BB <- matrix(list(0), nn, nn)
# for(rr in 1:nn) {
#   for(cc in 1:nn) {
#     BB[rr, cc] <- paste(rr, cc, sep = ",")
#   }
# }
# BB[seq(3), seq(8,9)] <- BB[seq(8,9), seq(3)] <- 0
# BB

model_list <- list(
  B = "unconstrained",
  U = "zero",
  Q = "diagonal and unequal",
  C = CC,
  c = cov_1,
  Z = "identity",
  A = "zero",
  R = "diagonal and equal"
)

control_list <- list(
  maxit = 5000
)


lwa_fit_1 <- MARSS(plank_1, model = model_list, control = control_list, method = "BFGS")

saveRDS(lwa_fit_1, file = "lwa_fit_1.rds")


model_list$c <- cov_2

lwa_fit_2 <- MARSS(plank_2, model = model_list, control = control_list, method = "BFGS")

saveRDS(lwa_fit_2, file = "lwa_fit_2.rds")




lwa_fit <- readRDS(file = "lwa_fit.rds")

B_fit <- coef(lwa_fit, type = "matrix")$B

C_fit <- coef(lwa_fit, type = "matrix")$C


pi_B <- det(B_fit)^2
d_pi_B <- dpiB(B_fit)

react <- log(max(svd(B_fit)$d))

ret_mu <- max(abs(eigen(B_fit)$values))

ret_sig <- max(abs(eigen(B_fit %x% B_fit)$values))

LL(B_fit, C_fit)

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

