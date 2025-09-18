setwd("/Users/nejcperme/Desktop/Faks/Podiplomski - 2.letnik/Monetary Economics 2/SVAR")
library(readxl)
library(tseries)
#library(quantmod)
library(ggplot2)
library(vars)
library(data.table)
library(jtools)
library(dplyr)
library(zoo)
library(forecast)
library(reshape2)
#install.packages(vars)


#load data
# Original time series are from the Journal of Applied Econometrics (JAE) data archive 
#(see data from Weber, C.E. (1995). Cyclical Output, cyclical unemployment, and Okun's coefficient: A new approach, Journal of Applied Econometrics, Vol. 10, pp. 433-335.).
Blanchard_Quah <- read_excel("Blanchard-Quah.xlsx")
dat <- Blanchard_Quah

#Transform data
DQ <- ts(dat$DQ, start = c(1948, 1), freq = 4)  # 1st differences of 100*log(real GNP) is taken, where "real GNP" is from the Weber data
U <- ts(dat$U , start = c(1948, 1), freq = 4)  # Unemployment rate


# Plot data
data <- cbind(DQ, U)
colnames(data) <- c("DQ", "U")
plot.ts(data)

#Unit root tests

adf.DQ <- ur.df(DQ, type = "drift", selectlags = c("AIC"))
summary(adf.DQ)
adf.U <- ur.df(U, type = "drift", selectlags = c("AIC"))
summary(adf.U)

#VAR lag length selection
info.var <- VARselect(data, lag.max = 8, type = "const")
info.var$selection

#VAR estimation
var.est1 <- VAR(data, p = 2,type = "const")

summary(var.est1)
svar.BQ=BQ(var.est1 )
svar.BQ$LRIM
Amats <- Acoef(var.est1)
P <- var.est1$p
Ident <- diag(var.est1$K)
mat1 <- matrix(0, var.est1$K, var.est1$K)
mat2 <- mat1
for(i in 1:P){
  mat1 <- mat1 - Amats[[i]] # O_k - A1 - A2
  mat2 <- mat2 - t(Amats[[i]]) #same but transposed
}
mat1 <- Ident + mat1 #I-A1-A2 - which inversed is the infinite sum moving average!!
mat2 <- Ident + mat2 #same but transposed
df <- summary(var.est1$varresult[[1]])$df[2]
SigmaU <- crossprod(resid(var.est1)) / df
eval <- solve(mat1) %*% SigmaU %*% solve(mat2) #(I-B)^-1PP'((I-B)^-1)'
lrim <- t(chol(eval)) #Cholesky! - chol(CC')=C long run impulse responses
colnames(lrim) <- colnames(data)
rownames(lrim) <- colnames(data)
cim <- mat1 %*% lrim # (I-A1-A2 )*C contemporaneous impulse responses
colnames(cim) <- colnames(lrim)
rownames(cim) <- colnames(lrim)
# Compute structural IRFs
irf <- Phi(var.est1, nstep=39)
rownames(irf) <- c("DQ", "U")
for(i in 1:dim(irf)[3]){
  irf[,,i] <- irf[,,i] %*% cim
}

# Extract responses to both shocks
irf_BQ <- array(0, dim=c(2, 40, 2))  # variables × horizons × shocks
for(s in 1:2){  # loop over shocks
  for(i in 1:dim(irf)[3]){
    irf_BQ[, i, s] <- irf[1:2, s, i]
  }
}

# Convert to tidy data frame
df_list <- list()
for(s in 1:2){
  tmp <- as.data.frame(t(irf_BQ[,,s]))
  colnames(tmp) <- c("Output Growth Rate", "Unemployment")
  tmp$Period <- 1:40
  tmp$Shock <- paste0("Shock ", s)
  df_list[[s]] <- tmp
}
irf_df <- do.call(rbind, df_list)

# Long format
irf_long <- reshape2::melt(irf_df, id.vars = c("Period", "Shock"),
                           variable.name = "Variable", value.name = "Response")

# Plot: facets by Variable, separate plots for each shock
p <- ggplot(irf_long, aes(x = Period, y = Response)) +
  geom_line(color = "blue") +
  facet_grid(Variable ~ Shock, scales = "free_y") +
  theme_bw() +
  labs(title = "Impulse Responses to Shocks",
       x = "Period", y = "") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, color="black"),
    panel.grid = element_blank(),
    strip.text = element_text(size = 12, face = "bold")
  )

# Show plot
print(p)

# Save figure
ggsave("IRFs_two_shocks.png", p, width = 8, height = 6, dpi = 300)
one.BQ <- irf(svar.BQ, response = NULL, impulse = NULL, 
              n.ahead = 40, ortho = TRUE, boot = TRUE, runs=1000)
plot(one.BQ)
#SVAR impulse response from U, negative spending measure will increase the unemployment rate (only in the first few quarters), growth goes down, then resumes but a long run effect on output is equal to 0!!! so the cumulative impulse response
#would go to zero. FIrstly it is negative (the level response of Q). Output goes from negative to a through and goes back up to the 0 level. 
#SVAR impulse responses from DQ - permnanet productivity shock. Positive effect on growth rate - output would be positive go up then dip slightly and then go up again. For unemployment, it looks like skill biased technological change. 
