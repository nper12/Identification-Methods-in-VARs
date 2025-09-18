
library(readxl)
library(tseries)
library(quantmod)
library(ggplot2)
library(vars)
library(data.table)
library(jtools)
library(dplyr)
library(zoo)
library(forecast)
library(tidyr)
library(dplyr)
library(tibble)
library(patchwork)
#install.packages("stargazer")
library(stargazer)
library(tidyverse)
library(svars)

data <- read_excel("Data/Cholesky.xlsx")

#Transform data
y <- ts(log(data$y)*100, start = c(2000, 1), freq = 12)  # output
p <- ts(log(data$hcpi)*100, start = c(2000, 1), freq = 12)  # consumer inflation
i <- ts(data$i, start = c(2000, 1), freq = 12)
carbon <- ts(data$carbon, start=c(2000,1), freq=12)
u <- ts(data$u, start=c(2000,1), freq=12)

# Plot data
data <- cbind(carbon, p, u,y, i)
data <- as.data.frame(data)
colnames(data) <- c("carbon","p","u", "y", "i")
plot.ts(data)

#VAR estimation
var.est1 <- vars::VAR(data, p = 3,type = "cons")

# Setting up the SVAR with Cholesky decomposition
a.mat <- diag(5)
diag(a.mat) <- NA
a.mat[2, 1] <- NA
a.mat[3, 1] <- NA
a.mat[3, 2] <- NA
a.mat[4, 1] <- NA
a.mat[4, 2] <- NA
a.mat[4, 3] <- NA
a.mat[5, 1] <- NA
a.mat[5, 2] <- NA
a.mat[5, 3] <- NA
a.mat[5, 4] <- NA
print(a.mat)

b.mat <- diag(5)
diag(b.mat) <- NA
print(b.mat)

svar.one <- SVAR(var.est1, Amat = a.mat, Bmat=b.mat, max.iter = 10000, 
                 hessian = TRUE)
#v svar ne mors dt identity matrixa za B recimo, pol sam matriko vn vrzs!!!
svar.one$A
svar.one$B

# Compute and store IRFs
par(mfrow = c(4, 1), mar = c(2.2, 2.2, 1, 1), cex = 0.6)
irf_y <- irf(svar.one, response = "y", impulse = "carbon", 
             n.ahead = 36, ortho = TRUE, boot = TRUE)
irf_p <- irf(svar.one, response = "p", impulse = "carbon", 
             n.ahead = 36, ortho = TRUE, boot = TRUE)
irf_i <- irf(svar.one, response = "i", impulse = "carbon", 
             n.ahead = 36, ortho = TRUE, boot = TRUE)
irf_u <- irf(svar.one, response = "u", impulse = "carbon", 
             n.ahead = 36, ortho = TRUE, boot = TRUE)
irf_carbon <- irf(svar.one, response = "carbon", impulse = "carbon", 
                  n.ahead = 36, ortho = TRUE, boot = TRUE)

plot(irf_y)
plot(irf_p)
plot(irf_i)
plot(irf_u)
plot(irf_carbon)

# Compute and store IRFs
two.int.pi <- irf(svar.one, response = "p", impulse = "i", 
                  n.ahead = 36, ortho = TRUE, boot = TRUE)
two.int.y <- irf(svar.one, response = "y", impulse = "i", 
                 n.ahead = 36, ortho = TRUE, boot = TRUE)
two.int.i <- irf(svar.one, response = "i", impulse = "i", 
                 n.ahead = 36, ortho = TRUE, boot = TRUE)
two.int.u <- irf(svar.one, response = "u", impulse = "i", 
                 n.ahead = 36, ortho = TRUE, boot = TRUE)
# Plot each IRF
plot(two.int.pi, main = "Response of p to monetary policy shock")
plot(two.int.y, main = "Response of y to monetary policy shock")
plot(two.int.i, main = "Response of i to monetary policy shock")
plot(two.int.u, main = "Response of u to monetary policy shock")

#---------------------- Getting the right dataframes ----------------------------
lagged_data <- data %>%
  mutate(across(everything(),
                .fns = list(
                  lag1 = ~lag(.x, 1),
                  lag2 = ~lag(.x, 2),
                  lag3 = ~lag(.x, 3)
                ),
                .names = "{.col}_{fn}"))
lagged_data$cons <- 1 
clean_frame <- na.omit(lagged_data)
Y <- as.matrix(clean_frame[, c("carbon","p", "u", "y", "i")])
regressors <- grep("_lag", names(clean_frame), value=TRUE)
X <- as.matrix(clean_frame[, c(regressors, "cons")])
Y <- t(Y)
X <- t(X)
#-------------------------------- OLS - matrix ------------------------------------------
B_hat <- Y%*%t(X)%*%solve(X%*%t(X)) #OLS matrix (X'X)^-1 X'Y
U_hat <- Y-B_hat%*%X #residual matrix
T <- ncol(X) #nobs
Sigma_u <- (U_hat %*% t(U_hat))/ (T-nrow(X)) #variance covariance matrix of residuals
invXtX <- solve(X%*% t(X))
vcov_B <- kronecker(Sigma_u, invXtX) #variance covariance matrix of regression coefficient
std_err <- sqrt(diag(vcov_B))
se_matrix <- matrix(std_err, nrow = nrow(B_hat), ncol = ncol(B_hat), byrow=TRUE)
colnames(se_matrix) <- rownames(X)
rownames(se_matrix) <- rownames(B_hat)

#VAR estimation
var.est1 <- vars::VAR(data, p = 3,type = "const")
results <- var.est1

# Extract standard errors for each equation
se_list <- lapply(results$varresult, function(eq_summary) {
  coef_table <- summary(eq_summary)$coefficients
  coef_table[, "Std. Error"]
})
coef <- var.est1$varresult$carbon$coefficients
coef
#-------------------------------- OLS - eq-by-eq ------------------------------------------
carbon_reg <- lm(carbon ~ carbon_lag1 + carbon_lag2+carbon_lag3+p_lag1+p_lag2+p_lag3+u_lag1+u_lag2+u_lag3+y_lag1+y_lag2+y_lag3+i_lag1+i_lag2+i_lag3, data=lagged_data)
p_reg <- lm(p ~ carbon_lag1 + carbon_lag2+carbon_lag3+p_lag1+p_lag2+p_lag3+u_lag1+u_lag2+u_lag3+y_lag1+y_lag2+y_lag3+i_lag1+i_lag2+i_lag3, data=lagged_data)
y_reg <- lm(y ~ carbon_lag1 + carbon_lag2+carbon_lag3+p_lag1+p_lag2+p_lag3+u_lag1+u_lag2+u_lag3+y_lag1+y_lag2+y_lag3+i_lag1+i_lag2+i_lag3, data=lagged_data)
u_reg <- lm(u ~ carbon_lag1 + carbon_lag2+carbon_lag3+p_lag1+p_lag2+p_lag3+u_lag1+u_lag2+u_lag3+y_lag1+y_lag2+y_lag3+i_lag1+i_lag2+i_lag3, data=lagged_data)
i_reg <- lm(i ~ carbon_lag1 + carbon_lag2+carbon_lag3+p_lag1+p_lag2+p_lag3+u_lag1+u_lag2+u_lag3+y_lag1+y_lag2+y_lag3+i_lag1+i_lag2+i_lag3, data=lagged_data)
vcov_carbon_reg <- vcov(carbon_reg)
vcov_p_reg <- vcov(p_reg)
vcov_y_reg <- vcov(y_reg)
vcov_u_reg <- vcov(u_reg)
vcov_i_reg <- vcov(i_reg)
se_carbon_reg <- sqrt(diag(vcov_carbon_reg))
se_p_reg <- sqrt(diag(vcov_p_reg))
se_y_reg <- sqrt(diag(vcov_y_reg))
se_u_reg <- sqrt(diag(vcov_u_reg))
se_i_reg <- sqrt(diag(vcov_i_reg))

#--------------- Manual structural shock from a recursive VAR  -----------------
residuals_var <- residuals(var.est1)
residuals_var
T <- nrow(residuals_var)
vcov <- t(residuals_var)%*%residuals_var/(T-15) #need also small sample adjustment
vcov(var.est1)
P <- chol(vcov) #t(P)=solve(svar.one$A)%*%svar.one$B
P <- t(P) 
P
diagonal_elements <- diag(P)
D <- array(0,dim=c(5,5))
diag(D) <- diagonal_elements
standardized_matrix <- P%*%solve(D)
solve(standardized_matrix) #v tem primeru sta enaka A in B!
svar.one$A
structural_matrix <- solve(svar.one$A)%*%svar.one$B
structural_matrix
#----------------------- STRUCTURAL SHOCKS ---------------------------
structural_shocks <- solve(structural_matrix)%*%t(residuals_var)
#oziroma
structural_shocks <- solve(P)%*%t(residuals_var)
structural_shocks <- as.data.frame(t(structural_shocks))
plot(structural_shocks$carbon)
plot(structural_shocks$i)
carbon_shock <- ts(structural_shocks$carbon, start=c(2000,3), frequency=12)
monetary_shock <- ts(structural_shocks$i, start=c(2000,3), frequency=12)
plot(carbon_shock)
carbon_plot <- autoplot(carbon_shock, color="blue") + theme_bw()+ theme(panel.grid = element_blank(),  
                                                                        axis.text = element_text(color = "black", size =13),   plot.title = element_text(hjust = 0.5, size = 15),
                                                                        axis.title.x = element_text(size = 15),   panel.border = element_rect(color = "black", linewidth = 1)) + ggtitle("Carbon policy shock") + ylab("") + xlab("Year")
monetary_plot <- autoplot(monetary_shock, color="blue") + theme_bw()+ theme(panel.grid = element_blank(),  
                                                                            axis.text = element_text(color = "black", size =13),   plot.title = element_text(hjust = 0.5, size = 15),
                                                                            axis.title.x = element_text(size = 15),   panel.border = element_rect(color = "black", linewidth = 1)) + ggtitle("Monetary policy shock")+ ylab("")+ xlab("Year")
shock_plot <- carbon_plot + monetary_plot
shock_plot
plot(shock_plot)
ggsave("Output/shock_plot.pdf", plot = shock_plot, width = 10, height = 3, dpi = 600)
#-------------------------- IRFs ročno -------------------------------
irf <- Phi(var.est1,nstep=39)
nstep <- abs(as.integer(39)) 
K <- 5 #določi koliko maš y
p <- 3 #določi koliko je endogenih lagov
A <- as.array(Acoef(var.est1)) #vsi ocenjeni koeficienti v A matrikah
if(nstep >= p){
  As <- array(0, dim = c(K, K, nstep + 1))
  for(i in (p + 1):(nstep + 1)){ #create nstep matrices filled with zeros
    As[, , i] <- matrix(0, nrow = K, ncol = K)
  }
} else {
  As <- array(0, dim = c(K, K, p))
}
for(i in 1:p){
  As[, , i] <- A[[i]] #fill the first p matrices with actual coefficients
}  
Phi <- array(0, dim=c(K, K, nstep + 1))
Phi[, ,1] <- diag(K)
Phi[, , 2] <- Phi[, , 1] %*% As[, , 1]
if (nstep > 1) {
  for (i in 3:(nstep + 1)) {
    tmp1 <- Phi[, , 1] %*% As[, , i-1]
    tmp2 <- matrix(0, nrow = K, ncol = K)
    idx <- (i - 2):1
    for (j in 1:(i - 2)) {
      tmp2 <- tmp2 + Phi[, , j+1] %*% As[, , idx[j]]
    }
    Phi[, , i] <- tmp1 + tmp2
  }
}
return(Phi)
}
Phi #equivalence
irf #equivalence
#----------------------------- Structural impulse responses z
rownames(irf) <- c("carbon", "p", "u", "y", "i")
for(i in 1: dim(irf)[3]){
  irf[, , i] <- irf[, , i] %*% P #this is how you get one unit shock impulse responses! - one unit shock to y!!! in our case carbon shock
}
irf_carbon = array(0, dim=c(5,40))
for(i in 1:dim(irf)[3]){
  irf_carbon[,i] <- irf[1:5,1,i]
}
# Loop through each row and plot
for (i in 1:nrow(irf_carbon)) {  # Corrected syntax
  plot(irf_carbon[i, ], type = "o", col = "blue", main = paste("Plot for Row", i), xlab = "Index", ylab = "Value")
}
irf
irf_carbon <- as.data.frame(t(irf_carbon))  # already done in your code
colnames(irf_carbon) <- c("Carbon", "Inflation", "Unemployment", "Output", "InterestRate")
irf_carbon$Period <- 1:nrow(irf_carbon)

# Convert to long format for faceting
irf_long <- irf_carbon %>%
  pivot_longer(
    cols = c("Carbon", "Inflation", "Unemployment", "Output", "InterestRate"),
    names_to = "Variable",
    values_to = "Response"
  )

# Faceted plot
ggplot(irf_long, aes(x = Period, y = Response)) +
  geom_line(color = "blue") +
  facet_wrap(~ Variable, scales = "free_y", ncol = 1) +
  theme_bw() +
  labs(title = "Impulse Responses to Carbon Shock",
       x = "Period", y = "") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, color="black"),
    panel.grid = element_blank(),
    strip.text = element_text(size = 12, face = "bold")
  )