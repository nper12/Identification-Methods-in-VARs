library(readxl)
library(boot)
library(tsDyn)
library(vars)
library(repr)
# Get original large dataset of BBE (2005): 120 series
data = read_excel("./Data/bbe_data.xlsx")
head(data)
# Standardizing data = all variables with mean 0 and standard deviation 1. 
# This step is crucial in PC analysis
data_s = scale(data[,2:121], center = TRUE, scale = TRUE)
pc_all = prcomp(data_s, center=FALSE, scale.=FALSE, rank. = 3) 
# summary(pc_all)
C = pc_all$x # saving the principal components
slow_vars = c("IP", "LHUR", "PUNEW", "IPP", "IPF", "IPC", "IPCD", "IPCN", "IPE", "IPI", "IPM",
              "IPMD", "IPMND", "IPMFG", "IPD", "IPN", "IPMIN", "IPUT", "IPXMCA", "PMI", "PMP",
              "GMPYQ", "GMYXPQ", "LHEL", "LHELX", "LHEM", "LHNAG", "LHU680", "LHU5", "LHU14",
              "LHU15", "LHU26", "LPNAG", "LP", "LPGD", "LPMI", "LPCC", "LPEM", "LPED", "LPEN",
              "LPSP", "LPTU", "LPT", "LPFR", "LPS", "LPGOV", "LPHRM", "LPMOSA", "PMEMP", "GMCQ",
              "GMCDQ", "GMCNQ", "GMCSQ", "GMCANQ", "PWFSA", "PWFCSA", "PWIMSA", "PWCMSA", "PSM99Q",
              "PU83", "PU84", "PU85", "PUC", "PUCD", "PUS", "PUXF", "PUXHS", "PUXM", "LEHCC", "LEHM")
data_slow = data_s[, slow_vars]
pc_slow = prcomp(data_slow, center=FALSE, scale.=FALSE, rank. = 3)
F_slow = pc_slow$x
# Next clean the PC of all space from the observed Y
reg = lm(C ~ F_slow + data_s[,"FYFF"])
#summary(reg)
F_hat = C - data.matrix(data_s[,"FYFF"])%*%reg$coefficients[5,] # cleaning and saving F_hat

data_var = data.frame(F_hat, "FYFF" = data_s[,"FYFF"])
var = VAR(data_var, p = 13)
#summary(var)

irf_point = irf(var, n.ahead = 48, impulse = "FYFF", response = "FYFF", boot = FALSE)

# Shock size of 25 basis points
impulse_sd = 0.25/sd(data$FYFF)
scale = impulse_sd/(irf_point$irf$FYFF[1]) # position of FYFF response at step 0


# Computing Loading Factors
reg_loadings = lm(data_s ~ 0 + F_hat + data_s[,"FYFF"])
loadings = reg_loadings$coefficients
# head(reg_loadings$coefficients)
#summary(reg_loadings)


#### BOOTSTRAPING ########

R = 500 # Number of simulations
nvars = 120 # Number of variables
nsteps = 49 # numbers of steps

IRFs = array(c(0,0,0), dim = c(nsteps,nvars,R))

var = lineVar(data_var, lag = 13, include = "const")
for(j in 1:R){    
  data_boot = VAR.boot(var, boot.scheme ="resample")
  var_boot = VAR(data_boot, lag = 13)
  irf1 = irf(var_boot, n.ahead = 48, impulse = "FYFF", boot = FALSE)
  for(i in 1:nvars){
    IRFs[,i,j] = (irf1$irf$FYFF %*% matrix(loadings[1:4, i]))*scale
  }
} ## Boot simulations done

# Extract the quantiles of IRFs we are interested: 90% confidence intervals in BBE
Upper = array(c(0,0), dim = c(nsteps, nvars))
for(k in 1:nsteps){
  for(i in 1:nvars){
    Upper[k,i] = quantile(IRFs[k,i,], probs = c(0.95))[1]
  }
}
Lower = array(c(0,0), dim = c(nsteps, nvars))
for(k in 1:nsteps){
  for(i in 1:nvars){
    Lower[k,i] = quantile(IRFs[k,i,], probs = c(0.05))[1]
  }
}
IRF = array(c(0,0), dim = c(nsteps, nvars))
for(k in 1:nsteps){
  for(i in 1:nvars){
    IRF[k,i] = quantile(IRFs[k,i,], probs = c(0.5))[1]
  }
}
rm(var_boot)
rm(IRFs)

# Select the Variables you are Interested in
# List of variables we are interested: FYFF, IP, CPI
variables = c(grep("^FYFF$", colnames(data_s)), grep("^IP$", colnames(data_s)), grep("^PUNEW$", colnames(data_s)),
              grep("^FYGM3$", colnames(data_s)), grep("^FYGT5$", colnames(data_s)), grep("^FMFBA$", colnames(data_s)),
              grep("^FM2$", colnames(data_s)), grep("^EXRJAN$", colnames(data_s)), grep("^PMCP$", colnames(data_s)),
              grep("^IPXMCA$", colnames(data_s)), grep("^GMCQ$", colnames(data_s)), grep("^GMCDQ$", colnames(data_s)), 
              grep("^GMCNQ$", colnames(data_s)), grep("^LHUR$", colnames(data_s)), grep("^PMEMP$", colnames(data_s)), 
              grep("^LEHM$", colnames(data_s)), grep("^HSFR$", colnames(data_s)), grep("^PMNO$", colnames(data_s)), 
              grep("^FSDXP$", colnames(data_s)), grep("^HHSNTN$", colnames(data_s))
)

transf_code = c(1, 5, 5,
                1, 1, 5,
                5, 5, 1,
                1, 5, 5,
                5, 1, 1,
                5, 1, 1,
                1, 1 
)

variable_names = c("Fed Funds Rate", "Industrial Production", "CPI",
                   "3m Treasury Bills", "5y Treasury Bonds", "Monetary Base",
                   "M2", "Exchange Rate Yen", "Commodity Price Index",
                   "Capacity Util Rate", "Personal Consumption", "Durable Cons",
                   "Nondurable Cons", "Unemployment", "Employment",
                   "Avg Hourly Earnings", "Housing Starts", "New Orders",
                   "Dividends", "Consumer Expectations" 
)

# Replicating Figure II in BBE (2005) - 3 Factors and Y = FYFF
# Change plot size to 15 x 10
options(repr.plot.width=12, repr.plot.height=8)

par(mfrow=c(5,4), 
    mar = c(2, 2, 2, 2))

for(i in variables){
  index = which(variables == i)
  if(transf_code[index] == 5){
    plot(cumsum(IRF[,i]), type ='l',lwd=2, main = variable_names[index],
         ylab= "", xlab="Steps", ylim=range(cumsum(Lower[,i]),cumsum(Upper[,i])),
         cex.main=1.8, cex.axis=1.3)
    lines(cumsum(Upper[,i]), lty=2, col="red")
    lines(cumsum(Lower[,i]), lty=2, col="red")
    abline(h=0)
  }
  else{
    plot(IRF[,i], type ='l',lwd=2, main = variable_names[index],
         ylab= "", xlab="Steps", ylim=range((Lower[,i]),(Upper[,i])),
         cex.main=1.8, cex.axis=1.3)
    lines((Upper[,i]), lty=2, col="red")
    lines((Lower[,i]), lty=2, col="red")
    abline(h=0)  
  }
}

options(repr.plot.width=12, repr.plot.height=6)

par(mfrow=c(1,2), 
    mar = c(2, 2, 2, 2))

plot(data_s[, variables[2]], type = "l", ylab = "", main = "Industrial Production")
lines(fitted(reg_loadings)[,variables[2]], lty=2, col="red")
legend(300,6, legend=c("Data", "Predicted by PC"), lty = 1:2, col = c("black", "red"), 
       box.lty=0)
plot(data_s[, variables[3]], type = "l", ylab = "", main = "CPI")
lines(fitted(reg_loadings)[,variables[3]], lty=2, col="red")

# Get the VAR point estimates
hor = 60
var = VAR(data_var, p = 13)
irf_point = irf(var, n.ahead = hor, boot = FALSE)

# Get IRFs for all of X we are interested in, Dimensions: (hor, key_nvars)
# Find loadings
results = summary(reg_loadings) # the warning comes because of FYFF

key_nvars = length(variables)
irf_X_pc1 = array(c(0,0), dim=c(hor+1, key_nvars))
irf_X_pc2 = array(c(0,0), dim=c(hor+1, key_nvars))
irf_X_pc3 = array(c(0,0), dim=c(hor+1, key_nvars))
irf_X_fyff = array(c(0,0), dim=c(hor+1, key_nvars))

for(i in 1:key_nvars){
  irf_X_pc1[,i] = irf_point$irf$PC1 %*% matrix(loadings[1:4, variables[i]])
  irf_X_pc2[,i] = irf_point$irf$PC2 %*% matrix(loadings[1:4, variables[i]])
  irf_X_pc3[,i] = irf_point$irf$PC3 %*% matrix(loadings[1:4, variables[i]])
  irf_X_fyff[,i] = (irf_point$irf$FYFF) %*% matrix(loadings[1:4, variables[i]])
}

# Get the IRFs squared and accumulate them
psi2_pc1 = array(0, dim = key_nvars)
psi2_pc2 = array(0, dim = key_nvars)
psi2_pc3 = array(0, dim = key_nvars)
psi2_fyff = array(0, dim = key_nvars)

for(i in 1:key_nvars){
  for(j in 1:hor){
    psi2_pc1[i] = psi2_pc1[i] + irf_X_pc1[j,i]^2
    psi2_pc2[i] = psi2_pc2[i] + irf_X_pc2[j,i]^2
    psi2_pc3[i] = psi2_pc3[i] + irf_X_pc3[j,i]^2
    psi2_fyff[i] = psi2_fyff[i] + irf_X_fyff[j,i]^2    
  }
}

var_total= array(0, dim = key_nvars)
var_fac= array(0, dim = key_nvars)
var_e= array(0, dim = key_nvars)

for(i in 1:key_nvars){
  var_fac[i] = psi2_pc1[i] + psi2_pc2[i] + psi2_pc3[i] + psi2_fyff[i]
  var_total[i] = psi2_pc1[i] + psi2_pc2[i] + psi2_pc3[i] + psi2_fyff[i] + results[[variables[i]]]$sigma^2
  var_e[i] = results[[variables[i]]]$sigma^2
}

table = data.frame("PC1" = round((psi2_pc1),3), "PC2" = round((psi2_pc2),3), "PC3" = round((psi2_pc3),3), "FYFF" = round((psi2_fyff),3),
                   "Factor_Y_total" = round(var_fac,3) ,"e" = round((var_e),3), "Total" = round(var_total,3))
row.names(table) = variable_names
table

# Replicating Table I in BBE (2005) - 3 Factors and Y = FYFF
r2 = array(0, dim = key_nvars)
for(i in 1:key_nvars){
  r2[i] = results[[variables[i]]]$r.squared
}

table2 = data.frame("Variables" = variable_names, "Contribution" = round((psi2_fyff/var_total),3), "R-squared" = round(r2,3))
table2