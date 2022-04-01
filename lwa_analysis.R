
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
  t() %>%
  zscore()

plank_2 <- plank %>%
  filter(Year >= 1989 & Year <= 1994) %>%
  select(-c("Year")) %>%
  as.matrix() %>%
  t() %>%
  zscore()

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
  R = "diagonal and equal",
  tinitx = 1
)

control_list <- list(
  maxit = 5000,
  allow.degen = TRUE
)


lwa_fit_1 <- MARSS(plank_1, model = model_list, control = control_list, method = "BFGS")

saveRDS(lwa_fit_1, file = "lwa_fit_1.rds")


model_list$c <- cov_2

lwa_fit_2 <- MARSS(plank_2, model = model_list, control = control_list, method = "BFGS")

saveRDS(lwa_fit_2, file = "lwa_fit_2.rds")



# lwa_fit_1 <- MARSS(plank_1, model = model_list, control = control_list)
# 
# model_list$c <- cov_2
# lwa_fit_2 <- MARSS(plank_2, model = model_list, control = control_list)


#### period 1 ####

# lwa_fit_1 <- readRDS(file = "lwa_fit.rds")

B_fit_1 <- coef(lwa_fit_1, type = "matrix")$B

C_fit_1 <- coef(lwa_fit_1, type = "matrix")$C


(pi_B_1 <- det(B_fit_1)^2)
d_pi_B_1 <- dpiB(B_fit_1)

(react_1 <- log(max(svd(B_fit_1)$d)))

(ret_mu_1 <- max(abs(eigen(B_fit_1)$values)))

(ret_sig_1 <- max(abs(eigen(B_fit_1 %x% B_fit_1)$values)))

LL(B_fit_1, C_fit_1)


#### period 2 ####

# lwa_fit_2 <- readRDS(file = "lwa_fit.rds")

B_fit_2 <- coef(lwa_fit_2, type = "matrix")$B

C_fit_2 <- coef(lwa_fit_2, type = "matrix")$C


(pi_B_2 <- det(B_fit_2)^2)
d_pi_B_2 <- dpiB(B_fit_2)

(react_2 <- log(max(svd(B_fit_2)$d)))

(ret_mu_2 <- max(abs(eigen(B_fit_2)$values)))

(ret_sig_2 <- max(abs(eigen(B_fit_2 %x% B_fit_2)$values)))

LL(B_fit_2, C_fit_2)


#### Ives et al (2003) ####

Bi <- matrix(
  c(0.5, 0, 0, 0,
    0, 0.6, -0.02, 0,
    0, 0, 0.77, 0,
    0, 0.1, 0, 0.55),
  4, 4, byrow = TRUE
)

(pi_B_i <- det(Bi)^2)
(d_pi_B_i <- dpiB(Bi))

(react_i <- log(max(svd(Bi)$d)))

(ret_mu_i <- max(abs(eigen(Bi)$values)))

(ret_sig_i <- max(abs(eigen(Bi %x% Bi)$values)))

