#========================================================#
# Quantitative ALM, Financial Econometrics & Derivatives 
# ML/DL using R, Python, Tensorflow by Sang-Heon Lee 
#
# https://kiandlee.blogspot.com
#——————————————————–#
# Vector Error Correction Model and Cointegration
#========================================================#
graphics.off()  # clear all graphs
rm(list = ls()) # remove all files from your workspace
library(urca)  # ca.jo, ur.df, finland
library(vars)  # vec2var
library(tsDyn) # VECM
#========================================================
# Data
#========================================================
# level data : 1958q1 – 1984q2
data(finland)
lev <- finland
lev_lag <- dplyr::lag(lev,1)
lev_lag <- lev_lag[-1,]
nr_lev <- nrow(lev)

# the sample period
yq <- expand.grid(1:4, 1958:1984)[1:nr_lev,]
colnames(yq) <- c("q", "yyyy")
rownames(yq) <- NULL

# quarterly centered dummy variables
yq$Q1 <- (yq$q==1)-1/4
yq$Q2 <- (yq$q==2)-1/4
yq$Q3 <- (yq$q==3)-1/4
dum_season <- yq[,-c(1,2)]

# 1st differenced data
dif <- as.data.frame(diff(as.matrix(lev), lag = 1))

#========================================================
# Cointegration Test
#========================================================
#———————————————-
# Johansen Cointegration Procedure
#———————————————-
# ecdet  = ‘none’  for no intercept 
#          ‘const’ for constant term
#          ‘trend’ for trend variable 
#          in cointegration
# type   =  eigen or trace test
# K      =  lag order of VAR
# spec   = “transitory” or “longrun”
# season = centered seasonal dummy (4:quarterly)
# dumvar = another dummy variables
#———————————————-
coint_ca.jo <- ca.jo(
  lev, ecdet = "none", type  = "eigen", K = 2, 
  spec = "transitory", season = 4, dumvar = NULL)
summary(coint_ca.jo) #we take the first two columns cause they have the largest eigenvalues
beta <- coint_ca.jo@V[,c(1:2)]
alpha <- coint_ca.jo@W[,c(1:2)]
beta 
alpha
#========================================================
# VECM model estimation
#========================================================
#————————————————
# VECM estimation
#————————————————
# VECM(data, lag, r = 1, 
#      include = c(“const”, “trend”, “none”, “both”),
#      beta = NULL, estim = c(“2OLS”, “ML”), 
#      LRinclude = c(“none”, “const”,”trend”, “both”), 
#      exogen = NULL)
#————————————————
VECM_tsDyn <- VECM(lev, lag=1, r=2,
                   estim = "ML",
                   LRinclude = "none",
                   beta = beta, #if betas (cointegrating relationship) are the same, then ETC coefficients are the same as alpha matrix!
                   exogen = dum_season)
VECM_tsDyn
#————————————————
# restricted VECM -> input for r
#————————————————
cajorls_ca.jo <- cajorls(coint_ca.jo, r=2)
cajorls_ca.jo #same shit as VECM_tsDyn!!!
#————————————————
# the VAR representation of a VECM from ca.jo
#————————————————
# vec2var: Transform a VECM to VAR in levels - this function zero restricts some beta coefficients - this is why it gives different results
# ca.jo is transformed to a VAR in level
# r : The cointegration rank 
#————————————————
vec2var_ca.jo <- vec2var(coint_ca.jo, r=2)
vec2var_ca.jo
#we could get the VECM parameters from VAR representation with a little manipulation
#——————————
# AR(1)
#——————————
VECM_tsDyn$coefficients[,4:7]
t(cajorls_ca.jo$rlm$coefficients[7:10,])
-vec2var_ca.jo$A$A2 
#——————————
# ECT : long-run impact matrix
#——————————
VECM_tsDyn$coefficients[,1:2] #this is different here cause I manually added a different beta matrix in VECM() function!
t(cajorls_ca.jo$rlm$coefficients[1:2,])
(vec2var_ca.jo$A$A1 + vec2var_ca.jo$A$A2 - diag(4))[,1:2]
#some other helpful shit https://www.r-bloggers.com/2021/12/some-interesting-issues-in-vecm-using-r/
