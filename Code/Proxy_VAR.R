library(vars)
library(tidyr)
library(dplyr)
library(ggplot2)
library(writexl)
library(readxl)
library(xlsx)
library(extraDistr)
library(lpirfs)
GKdata <- read_excel("Data/Proxy.xlsx")
#------------------------------ Data Wrangling ---------------------------------
var <- VAR(GKdata[, c("logip", "logcpi", "gs1", "ebp")], p = 12, type = "const") #creates a VAR with 12 lags involving a constant
res <- data.frame(stats::residuals(var)) #retains the residual matrix
p <- var$p #gets out the number of lags in VAR
seriesnames <- colnames(res) #get the names of the variables form the residual matrix
origorder <- seriesnames
if ("gs1" %in% seriesnames) { #reorder the list so that the dependent variable is first
  # order dependent first
  seriesnames <- seriesnames[seriesnames != "gs1"]
  seriesnames <- c("gs1", seriesnames) # Order the dependent variable first
} else {
  stop(paste("The series you are trying to instrument (", dependent, ") is not a series in the residual dataframe.", sep =""))
} #this gives you the original order in the var-cov - i dont know why seriesnames is not ordered okays
res[, "ff4_tc"] <- GKdata$ff4_tc[(p+1):length(GKdata$ff4_tc)] #you have to trim the column of first 12 observations, because of the VAR structure
# put together matrix of residuals
u <- as.matrix(res[, seriesnames])
# Now restrict to just the sample for the instrument (if necessary)
validrows <- !is.na(res[, "ff4_tc"])
u <- u[validrows,]
# Useful constants
T <- nrow(u)
k <- ncol(u)
# Some necessary parts of the covariance matrix
gamma <- (1 / (T - k*p - 1)) * t(u) %*% u
gamma_11 <- gamma[1,1]
gamma_21 <- matrix(gamma[2:nrow(gamma), 1], c(k-1,1))
gamma_22 <- matrix(gamma[2:nrow(gamma), 2:nrow(gamma)], c(k-1,k-1))
#---------------------------------First stage IV--------------------------------------------
firststage <- lm(res$gs1 ~ ff4_tc, res)
res[names(predict(firststage)), "fs"] <- stats::predict(firststage)
coefs <- rep(0, k) #k=4 in case of GK
names(coefs) <- seriesnames #the list is named
for (i in 1:k) {
  s <- seriesnames[i]
  if (s != "gs1") {
    secondstage <- stats::lm(stats::as.formula(paste(s, " ~ fs", sep = "")), res)
    coefs[i] <- secondstage$coefficients["fs"]
  } else {
    coefs[i] <- 1
  }
} #coef are S21/s11
s21_on_s11 <- matrix(coefs[2:k], c(k-1,1)) #3x1 matrix of coefficient (ratios)
Q <- (s21_on_s11 * gamma_11) %*% t(s21_on_s11) - (gamma_21 %*% t(s21_on_s11) + s21_on_s11 %*% t(gamma_21)) + gamma_22
s12s12 <- t(gamma_21 - s21_on_s11 * gamma_11) %*% solve(Q) %*% (gamma_21 - s21_on_s11 * gamma_11)
s11_squared <- gamma_11 - s12s12 #če ne bi mel s11 bi pomenil, da so to normalizirani koeficienti (če je s11 = 1)!
sp <- as.numeric(sqrt(s11_squared)) #s11
shockcolumn <- sp * coefs[origorder] #S21
shockcolumn
#---------------------------------IRFs-------------------------------------------
ma_representation <- Phi(var, 50)
ma_representation
Acoef(var)
irfs <- apply(ma_representation, 3, function(x) x %*% shockcolumn) #takes the MA coefficients and multiplies shockcol with each of the matrices. 
#The 3 is there because ma_representation is a three-dimensional array, this means that we have stacked matrices in an object.
irfs <- as.data.frame(t(irfs))
irfs#makes a dataframe from a matrix and transposes the matrix for 51 horizons
colnames(irfs) <- names(shockcolumn) #adds names to irfs
irfs <- mutate(irfs, horizon = 0:50) #adds horizon variable to irfs dataframe
irfs <- gather(irfs, key = variable, value = response, -horizon) #reshapse the dataframe, we get a long format, where
#horizon stays the same, a new column variable hold the names of variable in original dataframe and values are the ones that pertained to the each variable.
ggplot(irfs, aes(x = horizon, y = response, group = variable, color = variable)) + geom_line()


#--------------------------------- proving equivalence VARX ---------------------------
data_varx <- na.omit(GKdata)
data_varx <- data_varx[,-1]

#--------------------------creating the dataframes for OLS estimation
lagged_data <- data_varx %>%
  mutate(across(!ff4_tc,
                .fns = setNames(
                  lapply(1:12, function(k) ~lag(.x, k)),
                  paste0("lag", 1:12)
                ),
                .names = "{.col}_{fn}"))

# Add constant
lagged_data$cons <- 1 

# Drop NA rows (due to lags)
clean_frame <- na.omit(lagged_data)

# Dependent variables (endogenous system: 4 vars)
Y <- as.matrix(clean_frame[, c("logcpi", "logip", "gs1", "ebp")])

# Regressors: all lagged terms + exogenous ff4_tc contemporaneous + constant
regressors <- grep("_lag", names(clean_frame), value = TRUE)
X <- as.matrix(clean_frame[, c(regressors, "ff4_tc", "cons")])

# Transpose for later use
Y <- t(Y)
X <- t(X)

#---------------------------------------- OLS estimation
B_hat <- Y%*%t(X)%*%solve(X%*%t(X)) #OLS matrix YX'(XX')^-1

# number of endogenous variables
endog_vars <- c("logcpi", "logip", "gs1", "ebp")
p <- 12   # number of lags

# List to store A matrices
A_mats <- vector("list", p)

for (lag in 1:p) {
  # Build exact column names for this lag
  cols <- paste0(endog_vars, "_lag", lag)
  
  # Find which column positions exist in B_hat
  col_pos <- match(cols, colnames(B_hat))
  
  # Extract the columns from B_hat
  A_mats[[lag]] <- B_hat[, col_pos, drop = FALSE]
  
  # Assign row and column names
  rownames(A_mats[[lag]]) <- endog_vars
  colnames(A_mats[[lag]]) <- cols
}

# Exogenous contemporaneous matrix
B_mat <- B_hat[, "ff4_tc", drop = FALSE]
rownames(B_mat) <- endog_vars
colnames(B_mat) <- "ff4_tc"

# Constant vector
c_vec <- B_hat[, "cons", drop = FALSE]
rownames(c_vec) <- endog_vars
colnames(c_vec) <- "cons"

A_mats
B_mat

#----------------------- reduced form IRF with recursive formula -------------------
nstep <- abs(as.integer(50)) 
K <- 4 #določi koliko maš y
p <- 12 #določi koliko je endogenih lagov
if(nstep >= p){
  As <- array(0, dim = c(K, K, nstep + 1))
  for(i in (p + 1):(nstep + 1)){ #create nstep matrices filled with zeros
    As[, , i] <- matrix(0, nrow = K, ncol = K)
  }
} else {
  As <- array(0, dim = c(K, K, p))
}
for(i in 1:p){
  As[, , i] <- A_mats[[i]] #fill the first p matrices with actual coefficients
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
Phi #equivalence

#---------------------------- structural IRF --------------------
rownames(Phi) <- c("logcpi", "logip", "gs1", "ebp")
for(i in 1: dim(Phi)[3]){
  Phi[, , i] <- Phi[, , i] %*% B_mat  #this is how you get the partially identified structural coefficients from proxy VAR in a VARX framework! - gs1 response is normalized to 1!
}
irf_monetary = array(0, dim=c(4,51))
for(i in 1:dim(Phi)[3]){
  irf_monetary[,i] <- Phi[1:4,1,i]
}
# Loop through each row and plot
for (i in 1:nrow(irf_monetary)) {  # Corrected syntax
  plot(irf_monetary[i, ], type = "o", col = "blue", main = paste("Plot for Row", i), xlab = "Index", ylab = "Value")
}
irf
irf_monetary <- as.data.frame(t(irf_monetary))  # already done in your code
colnames(irf_monetary) <- c("Inflation", "Output", "Yield", "Premium")
irf_monetary$Period <- 1:nrow(irf_monetary)

# Convert to long format for faceting
irf_long <- irf_monetary %>%
  pivot_longer(
    cols = c("Inflation", "Output", "Yield", "Premium"),
    names_to = "Variable",
    values_to = "Response"
  )

# Faceted plot
ggplot(irf_long, aes(x = Period, y = Response)) +
  geom_line(color = "blue") +
  facet_wrap(~ Variable, scales = "free_y", ncol = 1) +
  theme_bw() +
  labs(title = "Impulse Responses to Monetary Shock",
       x = "Period", y = "") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, color="black"),
    panel.grid = element_blank(),
    strip.text = element_text(size = 12, face = "bold")
  )

#------------------------------------ structural 
#---------------- reduced and structural IRF with companion matrix - yes you can do it all at once!!!!----------
matpow <- function(A, n) {
  if (n == 0) return(diag(nrow(A)))
  if (n == 1) return(A)
  res <- A
  for (i in 2:n) res <- res %*% A
  return(res)
} #if n=0 it returns a diagonal matrix, if n=1 it returns A, if n is more, it defines A as res and mutliplies it by A and and then again multiplies res with A etc (loop)
nstep <- 50   # number of steps for IRFs
K <- 4        # number of endogenous variables
p <- 12        # number of lags of y
M <- 1        # number of exogenous variables
s <- 1        # number of lags of x
endog_vars <- c("logcpi", "logip", "gs1", "ebp")

# List to store A matrices
A_y <- vector("list", p)

A_x <- NULL # if you have lags of exogenous variable!

for (lag in 1:p) {
  # Build exact column names for this lag
  cols <- paste0(endog_vars, "_lag", lag)
  
  # Find which column positions exist in B_hat
  col_pos <- match(cols, colnames(B_hat))
  
  # Extract the columns from B_hat
  A_y[[lag]] <- B_hat[, col_pos, drop = FALSE]
  
  # Assign row and column names
  rownames(A_y[[lag]]) <- endog_vars
  colnames(A_y[[lag]]) <- cols
}
A_x <- NULL

# Exogenous contemporaneous matrix
B <- B_hat[, "ff4_tc", drop = FALSE]
rownames(B) <- endog_vars
colnames(B) <- "ff4_tc"

# Build companion matrix A_comp
Kp <- K*p
Ms <- M*s
dimA <- Kp + Ms
A_comp <- matrix(0, nrow = dimA, ncol = dimA)

# Endogenous lags
for(i in 1:p){
  A_comp[1:K, ((i-1)*K+1):(i*K)] <- A_y[[i]]
}

# Lagged exogenous coefficients (only if s > 0)
if (s > 0 && !is.null(A_x)){
  for(i in 1:s){
    A_comp[1:K, Kp + ((i-1)*M+1):(i*M)] <- A_x[[i]]
  }
}

# Shift matrices for endogenous lags
if(p > 1){
  A_comp[(K+1):Kp, 1:(K*(p-1))] <- diag(K*(p-1))
}

# Shift matrices for exogenous lags
if(s > 1){
  A_comp[(Kp+M+1):(Kp+Ms), Kp + 1:(M*(s-1))] <- diag(M*(s-1))
}


# ------------------------------- Build B0 (dimA x M) 
B0 <- matrix(0, nrow = dimA, ncol = M)
B0[1:K, 1:M] <- B  # contemporaneous effect
if(s > 0){
  B0[(Kp+1):(Kp+M), 1:M] <- diag(M) # identity for xt in state vector
}
# Selection matrix J (K x dimA)
J <- matrix(0, nrow = K, ncol = dimA)
J[1:K, 1:K] <- diag(K)

# Compute IRFs to exogenous shocks
Phi <- array(0, dim = c(K, M, nstep+1))
Phi[,,1] <- J %*% B0

if(nstep > 0){
  for(h in 1:nstep){
    Phi[,,h+1] <- J %*% matpow(A_comp, h) %*% B0
  }
}
Phi

#----------------------------- structural shocks ------------------------------
res <- res[,1:4]
head(res)
head(gamma)
head(shockcolumn)
shockcolumn <- as.data.frame(shockcolumn)
shockcolumn <- t(shockcolumn)
shockcolumn <- shockcolumn[,seriesnames]
shockcolumn <- as.data.frame(shockcolumn)
shockcolumn <- as.matrix(shockcolumn)
shockcolumn <- t(shockcolumn)
res <- res[,seriesnames]
res <- as.matrix(res)
res

#one way to get structural shocks
structural_shock <- res%*%solve(gamma)%*%t(shockcolumn) 
shockcolumn <- as.matrix(shockcolumn)
#second way to get structural shocks
factor <- t((shockcolumn%*%solve(gamma)))%*%solve(shockcolumn%*%solve(gamma)%*%t(shockcolumn))
structural_shock_2 <- t(factor)%*%t(res)
structural_shock <- structural_shock[116:384, ]
structural_shock <- as.data.frame(structural_shock)
GKdata_subset <- GKdata[128:396, ]
GKdata_subset <- GKdata_subset[,2:5]

#check with lps!!!
lp <- lpirfs::lp_lin_iv(endog_data=GKdata_subset, shock=structural_shock, lags_endog_lin = 12, lags_exog = 2, hor=50, trend=0, confint=1.96, use_nw = TRUE)
plot(lp)

#različno od ramey ker ona sam uzame ff4_tc, kar je pr ns bascially nečisti proxy za šok. Strukturne šoke potem dobiš posebi vn!
#-------------------------- Wild Bootstrap ----------------------
# Initialize a list to store bootstrap results
bootstrap_results <- vector("list", 3000)

# Extract the estimated residuals and coefficients from the original VAR model
residuals_original <- residuals(var)
A_hat_list <- Bcoef(var)  # Coefficient matrices A1 to A12, stored as a list
constant_term <- A_hat_list[,ncol(A_hat_list)]
A_hat_list <- A_hat_list[,-ncol(A_hat_list)]

# Get dimensions of the original data
n_periods <- nrow(residuals_original)  # Number of time periods (384)
n_vars <- ncol(residuals_original)  # Number of variables (4 in this case)
p <- 12  # Number of lags in the VAR model
new_proxy <- rep(NA, n_periods)

# Original dataset has 396 rows, but we are focused on the first 384 rows + 12 original
n_total <- 396

# Loop over 5000 bootstrap iterations
for (iter in 1:3000) {
  
  # Step 1: Initialize the simulated data
  simulated_data <- matrix(NA, nrow = n_total, ncol = n_vars)
  colnames(simulated_data) <- colnames(GKdata[, c("logip", "logcpi", "gs1", "ebp")])
  
  # Step 2: Plug the first 12 rows with the original data
  simulated_data[1:12, ] <- as.matrix(GKdata[1:12, c("logip", "logcpi", "gs1", "ebp")])
  
  # Step 3: For each time period (13 to 396)
  for (t in (p + 1):n_total) {
    # Generate a single Rademacher multiplier for the current time period
    rademacher_multiplier <- rsign(1)  # Generate one Rademacher value
    
    # Apply the Rademacher multiplier to the corresponding residuals
    new_residuals <- rademacher_multiplier * residuals_original[t - p, ]
    
    # Adjust the proxy value with the same Rademacher multiplier
      original_proxy_value <- GKdata$ff4_tc[t]  # Original proxy value
      new_proxy[t-p] <- rademacher_multiplier * original_proxy_value  # Adjusted proxy value
        
        # Initialize a vector to store the sum of all lag effects for this period
        new_values <- rep(0, n_vars)
        
        # Apply all 12 lags, summing the contributions from each lag
        for (lag in 1:12) {
          # Get the corresponding lagged values from the simulated data (rows t-1 to t-12)
          lagged_values <- as.matrix(simulated_data[t - lag, ])
          
          # Get the coefficients for this lag (for all variables: 4 columns)
          start_idx <- (lag - 1) * n_vars + 1      # Start index for the coefficient matrix
          end_idx <- lag * n_vars                  # End index for the coefficient matrix
          A_lag <- as.matrix(A_hat_list[, start_idx:end_idx]) # Coefficients for this lag (4 columns)
          
          # Multiply the lagged row values by the corresponding columns in A_lag
          new_values <- new_values + A_lag %*% lagged_values
        }
        
        # Add the constant term for each variable
        new_values <- new_values + constant_term
        
        # Add the newly generated residuals for the current period
        new_values <- new_values + new_residuals
        
        # Store the new values in the simulated data matrix
        simulated_data[t, ] <- new_values
      }
  
  # Step 6: Re-estimate the VAR model on the simulated data
  simulated_data_df <- as.data.frame(simulated_data)
  var_boot <- VAR(simulated_data_df, p = 12, type = "const")
  
  # Step 7: Continue with proxy VAR steps as before
  res_boot <- data.frame(stats::residuals(var_boot))
  p_boot <- var_boot$p
  seriesnames_boot <- colnames(res_boot)
  origorder_boot <- seriesnames_boot
  
  if ("gs1" %in% seriesnames_boot) {
    seriesnames_boot <- seriesnames_boot[seriesnames_boot != "gs1"]
    seriesnames_boot <- c("gs1", seriesnames_boot)  # Dependent variable first
  } else {
    stop(paste("The series you are trying to instrument (", dependent, ") is not in the residual dataframe.", sep = ""))
  }
  
  #add new_proxy to residual matrix
  res_boot[, "new_proxy"] <- new_proxy
  
  # Perform the proxy VAR steps
  u_boot <- as.matrix(res_boot[, seriesnames_boot])
  validrows <- !is.na(res_boot[, "new_proxy"])
  u_boot <- u_boot[validrows, ]
  
  T_boot <- nrow(u_boot)
  k_boot <- ncol(u_boot)
  gamma_boot <- (1 / (T_boot - k_boot * p_boot - 1)) * t(u_boot) %*% u_boot
  gamma_11_boot <- gamma_boot[1, 1]
  gamma_21_boot <- matrix(gamma_boot[2:nrow(gamma_boot), 1], c(k_boot - 1, 1))
  gamma_22_boot <- matrix(gamma_boot[2:nrow(gamma_boot), 2:nrow(gamma_boot)], c(k_boot - 1, k_boot - 1))
  
  firststage_boot <- lm(res_boot$gs1 ~ new_proxy, res_boot)
  res_boot[names(predict(firststage_boot)), "fs"] <- stats::predict(firststage_boot)
  
  coefs_boot <- rep(0, k_boot)
  names(coefs_boot) <- seriesnames_boot
  for (i in 1:k_boot) {
    s <- seriesnames_boot[i]
    if (s != "gs1") {
      secondstage_boot <- stats::lm(stats::as.formula(paste(s, " ~ fs", sep = "")), res_boot)
      coefs_boot[i] <- secondstage_boot$coefficients["fs"]
    } else {
      coefs_boot[i] <- 1
    }
  }
  
  s21_on_s11_boot <- matrix(coefs_boot[2:k_boot], c(k_boot - 1, 1))
  Q_boot <- (s21_on_s11_boot * gamma_11_boot) %*% t(s21_on_s11_boot) - (gamma_21_boot %*% t(s21_on_s11_boot) + s21_on_s11_boot %*% t(gamma_21_boot)) + gamma_22_boot
  s12s12_boot <- t(gamma_21_boot - s21_on_s11_boot * gamma_11_boot) %*% solve(Q_boot) %*% (gamma_21_boot - s21_on_s11_boot * gamma_11_boot)
  s11_squared_boot <- gamma_11_boot - s12s12_boot
  sp_boot <- as.numeric(sqrt(s11_squared_boot))
  
  shockcolumn_boot <- sp_boot * coefs_boot[origorder_boot]
  ma_representation_boot <- Phi(var_boot, 50)
  
  irfs_boot <- apply(ma_representation_boot, 3, function(x) x %*% shockcolumn_boot)
  irfs_boot <- as.data.frame(t(irfs_boot))
  
  # Store the IRFs in the bootstrap results list
  bootstrap_results[[iter]] <- irfs_boot
}

# Convert the list of results to a three-dimensional array
irfs_array <- array(unlist(bootstrap_results), dim = c(nrow(irfs_boot), ncol(irfs_boot), 3000))

# Calculate the mean IRFs across all bootstrap iterations
mean_irfs <- apply(irfs_array, c(1, 2), mean)

# Calculate upper and lower confidence bands
alpha <- 0.05  # for 95% confidence interval
lower_band <- apply(irfs_array, c(1, 2), quantile, probs = alpha / 2)
upper_band <- apply(irfs_array, c(1, 2), quantile, probs = 1 - alpha / 2)

# Combine results into a data frame for easier visualization
results_df <- data.frame(
  time = 1:nrow(mean_irfs),
  mean_irfs = list(mean_irfs),
  lower_band = list(lower_band),
  upper_band = list(upper_band)
)
colnames(lower_band) <- names(shockcolumn_boot)
colnames(upper_band) <- names(shockcolumn_boot)
lower_band <- as.data.frame(lower_band)
upper_band <- as.data.frame(upper_band)
irfs <- as.data.frame(irfs)
horizon <- c(0:50)
horizon <- as.data.frame(horizon)
lower_band <- cbind(horizon, lower_band)
upper_band <- cbind(horizon, upper_band)
colnames(lower_band)[-1] <- paste0(colnames(lower_band)[-1], "_lower")
colnames(upper_band)[-1] <- paste0(colnames(upper_band)[-1], "_upper")
combined_data <- merge(irfs, lower_band, by = "horizon")
combined_data <- merge(combined_data, upper_band, by = "horizon")
combined_data

# Only keep unique variables for reshaping
vars <- c("logip", "logcpi", "gs1", "ebp")

long_data <- combined_data %>%
  pivot_longer(
    cols = all_of(c("response", paste0(vars, "_lower"), paste0(vars, "_upper"))),
    names_to = "temp",
    values_to = "value"
  ) %>%
  mutate(
    type = case_when(
      temp == "response" ~ "irf",
      grepl("_lower$", temp) ~ "lower",
      grepl("_upper$", temp) ~ "upper"
    ),
    variable = case_when(
      temp == "response" ~ variable,  # keep existing variable column
      TRUE ~ sub("_(lower|upper)$", "", temp)
    )
  ) %>%
  select(horizon, variable, type, value)

ggplot(long_data, aes(x = horizon, y = value, group = type)) +
  # IRF line
  geom_line(data = subset(long_data, type == "irf"), color = "black", linewidth = 1) +
  # Confidence interval edges
  geom_line(data = subset(long_data, type %in% c("lower", "upper")), 
            aes(group = type), linetype = "dashed", color = "black", linewidth = 0.8) +
  facet_wrap(~ variable, scales = "free_y") +
  theme_minimal() +
  labs(title = "Impulse Response Functions with Confidence Intervals",
       x = "Horizon",
       y = "Response")
#-------------------------- Moving Block Bootstrap ----------------------------
# --- Inputs (placeholders, replace with your actual objects)
nboot <- 3000
residuals_original <- residuals(var)
T_est <- nrow(residuals_original)     # sample size
nvar  <- ncol(residuals_original)     # number of variables in VAR
p     <- 12                 # VAR lags (if needed later)
proxy <- GKdata$ff4_tc[(p+1):length(GKdata$ff4_tc)] %>%
  replace(is.na(.), 0) %>% as.matrix()
proxy
BlockSize <- round(5.03 * T_est^(0.25))
nBlock    <- ceiling(T_est / BlockSize)
horizon <- 0:50
irf_boot_list <- vector("list", nboot)
Y_init <- as.matrix(GKdata[1:p, c("logip","logcpi","gs1","ebp")])
A_hat_list <- Bcoef(var)  # Coefficient matrices A1 to A12, stored as a list
constant_term <- A_hat_list[,ncol(A_hat_list)]
A_hat_list <- A_hat_list[,-ncol(A_hat_list)]
Y_init
#creating matrices where to put blocks!
VARBlocks   <- array(0, dim = c(BlockSize, nvar, T_est - BlockSize + 1))
ProxyBlocks <- array(0, dim = c(BlockSize, ncol(proxy), T_est - BlockSize + 1))

#putting in the blocks
for (j in 1:(T_est - BlockSize + 1)) {
  VARBlocks[ , , j]   <- residuals_original[j:(BlockSize + j - 1), ]
  ProxyBlocks[ , , j] <- proxy[j:(BlockSize + j - 1), ]
}

#centering the residuals and proxies
VARcentering <- matrix(0, nrow = BlockSize, ncol = nvar)
for (j in 1:BlockSize) {
  VARcentering[j, ] <- colMeans(residuals_original[j:(T_est - BlockSize + j), , drop = FALSE])
}
VARcentering <- VARcentering[rep(1:BlockSize, length.out = nBlock * BlockSize), ]
VARcentering <- VARcentering[1:T_est, , drop = FALSE]
Proxycentering <- matrix(0, nrow = BlockSize, ncol = ncol(proxy))
for (j in 1:BlockSize) {
  subProxy <- proxy[j:(T_est - BlockSize + j), , drop = FALSE]
  nonzero_rows_sub <- subProxy[, 1] != 0
  nonzero_rows_all <- proxy[, 1] != 0
  Proxycentering[j, ] <- colMeans(subProxy[nonzero_rows_sub, , drop = FALSE]) -
    colMeans(proxy[nonzero_rows_all, , drop = FALSE])
}
Proxycentering <- Proxycentering[rep(1:BlockSize, length.out = nBlock * BlockSize), ] %>% as.matrix()
Proxycentering <- Proxycentering[1:T_est, , drop = FALSE]

#Part 2 - random draws out of the matrices

for (b in 1:nboot) {

# Randomly choose block indices
index <- ceiling((T_est - BlockSize + 1) * runif(nBlock))

# Bootstrap VAR residuals
bootU <- matrix(0, nrow = nBlock * BlockSize, ncol = ncol(residuals_original))
for (kk in 1:nBlock) {
  bootU[(1 + BlockSize * (kk - 1)):(BlockSize * kk), ] <- VARBlocks[ , , index[kk]]
}
bootU <- bootU[1:T_est, , drop = FALSE]

# Bootstrap proxies
bootProxy <- matrix(0, nrow = nBlock * BlockSize, ncol = ncol(proxy))
for (kk in 1:nBlock) {
  bootProxy[(1 + BlockSize * (kk - 1)):(BlockSize * kk), ] <- ProxyBlocks[ , , index[kk]]
}
bootProxy <- bootProxy[1:T_est, , drop = FALSE]

# Centering both
bootU <- bootU - VARcentering
for (kk in 1:ncol(proxy)) {
  nonzero_idx <- bootProxy[, kk] != 0
  bootProxy[nonzero_idx, kk] <- bootProxy[nonzero_idx, kk] - Proxycentering[nonzero_idx, kk]
}

#new sample construction and estimation
Y_boot <- matrix(NA, nrow = T_est + p, ncol = nvar)
colnames(Y_boot) <- c("logip", "logcpi", "gs1", "ebp")
Y_boot[1:p, ] <- Y_init
Y_boot
for (t in (p + 1):(T_est + p)) {
  lagged_vals <- as.vector(t(Y_boot[(t - 1):(t - p), ]))
  Y_boot[t, ] <- constant_term + 
    as.numeric(A_hat_list %*% lagged_vals) + 
    bootU[t - p, ]
}
Y_boot

# Step 2: re-estimate VAR on bootstrapped data
var_boot <- VAR(as.data.frame(Y_boot), p = p, type = "const")

# Step 3: prepare residuals and proxy for proxy VAR steps
res_boot <- data.frame(residuals(var_boot))
res_boot$ff4_tc <- bootProxy[,1]  # assuming one proxy variable
valid_idx <- !is.na(res_boot$ff4_tc)
res_boot_valid <- res_boot[valid_idx, ]
res_boot_valid

# Step 4: proxy VAR identification
# First stage: instrument
first_stage <- lm(gs1 ~ ff4_tc, data = res_boot_valid)
res_boot_valid$fs <- predict(first_stage)

# Second stage: get coefficients
var_names <- colnames(residuals_original)  # "logip", "logcpi", "gs1", "ebp"

# Second stage: regress each variable on fs, except gs1
coefs_boot <- sapply(var_names, function(v) {
  if (v == "gs1") {
    return(1)  # gs1 coefficient is 1
  } else {
    lm(as.formula(paste(v, "~ fs")), data = res_boot_valid)$coefficients["fs"]
  }
})

coefs_boot

# 4. Compute s11, s21, scaling as in your code
u_boot <- as.matrix(res_boot_valid[, !(colnames(res_boot_valid) %in% c("ff4_tc","fs"))])
gamma_boot <- t(u_boot) %*% u_boot / (nrow(u_boot) - ncol(u_boot)*p - 1)
gamma_11 <- gamma_boot[1,1]
gamma_21 <- gamma_boot[2:ncol(u_boot), 1, drop = FALSE]
gamma_22 <- gamma_boot[2:ncol(u_boot), 2:ncol(u_boot)]

s21_on_s11 <- matrix(coefs_boot[-1], ncol = 1)
Q_boot <- (s21_on_s11 * gamma_11) %*% t(s21_on_s11) - (gamma_21 %*% t(s21_on_s11) + s21_on_s11 %*% t(gamma_21)) + gamma_22
s12s12_boot <- t(gamma_21 - s21_on_s11 * gamma_11) %*% solve(Q_boot) %*% (gamma_21 - s21_on_s11 * gamma_11)
s11_squared <- gamma_11 - s12s12_boot
sp_boot <- sqrt(s11_squared)

shock_vector <- c(sp_boot) * coefs_boot

# 5. MA representation
ma_rep <- Phi(var_boot, max(horizon))  # or use your VAR estimate
irfs_boot <- sapply(1:length(horizon), function(h) ma_rep[,,h] %*% shock_vector)
irf_boot_list[[b]] <- t(irfs_boot)
}

# ------------------- Aggregate IRFs 
# Convert list of IRFs to 3D array: time x variables x bootstrap iterations
irf_array <- array(unlist(irf_boot_list),
                   dim = c(nrow(irf_boot_list[[1]]),
                           ncol(irf_boot_list[[1]]),
                           length(irf_boot_list)))

# Compute mean and confidence bands
mean_irfs  <- apply(irf_array, c(1,2), mean)
lower_band <- apply(irf_array, c(1,2), quantile, probs = 0.025)
upper_band <- apply(irf_array, c(1,2), quantile, probs = 0.975)

# ------------------- Combine into a single dataframe 
irfs_df <- data.frame(horizon = horizon, mean_irfs)
lower_df <- data.frame(horizon = horizon, lower_band)
upper_df <- data.frame(horizon = horizon, upper_band)

# Rename columns to follow consistent pattern var_type
colnames(irfs_df)[-1]   <- paste0(colnames(irfs_df)[-1], "_mean")
colnames(lower_df)[-1]  <- paste0(colnames(lower_df)[-1], "_lower")
colnames(upper_df)[-1]  <- paste0(colnames(upper_df)[-1], "_upper")

# Merge all into one dataframe
combined <- Reduce(function(x, y) merge(x, y, by = "horizon"), list(irfs_df, lower_df, upper_df))

# ------------------- Convert to long format for plotting
library(tidyr)
long_data <- combined %>%
  pivot_longer(
    cols = -horizon,
    names_to = c("var", "type"),
    names_sep = "_",
    values_to = "value"
  )

# Convert combined data to long format for plotting
long_data_wide <- combined %>%
  pivot_longer(
    cols = -horizon,
    names_to = c("var", "type"),
    names_sep = "_"
  ) %>%
  pivot_wider(
    names_from = type,
    values_from = value
  )

# Map internal variable names to desired plot labels
var_labels <- c("X1" = "logip", "X2" = "logcpi", "X3" = "gs1", "X4" = "ebp")
long_data_wide <- long_data_wide %>%
  mutate(var = recode(var, !!!var_labels))

# Plot: IRF line + upper/lower CI as black dashed lines
ggplot(long_data_wide, aes(x = horizon)) +
  # IRF line
  geom_line(aes(y = mean), color = "black", linewidth = 1) +
  # Upper CI
  geom_line(aes(y = upper), color = "black", linetype = "dashed", linewidth = 0.8) +
  # Lower CI
  geom_line(aes(y = lower), color = "black", linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~ var, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Impulse Response Functions (MBB)",
    x = "Horizon",
    y = "Response"
  ) +
  theme(legend.position = "none")
 #these are wider!! Courtesy of Kanzig github! - oil and news shocks
#-------------------------- robustness check -------------------
data_varx <- na.omit(GKdata)
data_varx <- data_varx[,-c(1,6)]
var
var <- VAR(data_varx, p = 12, type = "const") #creates a VAR with 12 lags involving a constant
res <- data.frame(stats::residuals(var)) #retains the residual matrix
p <- var$p #gets out the number of lags in VAR
seriesnames <- colnames(res) #get the names of the variables form the residual matrix
origorder <- seriesnames
if ("gs1" %in% seriesnames) { #reorder the list so that the dependent variable is first
  # order dependent first
  seriesnames <- seriesnames[seriesnames != "gs1"]
  seriesnames <- c("gs1", seriesnames) # Order the dependent variable first
} else {
  stop(paste("The series you are trying to instrument (", dependent, ") is not a series in the residual dataframe.", sep =""))
} #this gives you the original order in the var-cov - i dont know why seriesnames is not ordered okays
# Remove NA values from ff4_tc
ff4_tc_clean <- na.omit(GKdata$ff4_tc)
# Trim the first p observations to align with residuals
ff4_tc_clean <- ff4_tc_clean[(p+1):length(ff4_tc_clean)]

# Trim the first p observations for ff4_tc
res$ff4_tc <- ff4_tc_clean
# put together matrix of residuals
u <- as.matrix(res[, seriesnames])
# Now restrict to just the sample for the instrument (if necessary)
validrows <- !is.na(res[, "ff4_tc"])
u <- u[validrows,]
# Useful constants
T <- nrow(u)
k <- ncol(u)
# Some necessary parts of the covariance matrix
gamma <- (1 / (T - k*p - 1)) * t(u) %*% u
gamma_11 <- gamma[1,1]
gamma_21 <- matrix(gamma[2:nrow(gamma), 1], c(k-1,1))
gamma_22 <- matrix(gamma[2:nrow(gamma), 2:nrow(gamma)], c(k-1,k-1))
#---------------------------------First stage IV--------------------------------------------
firststage <- lm(res$gs1 ~ ff4_tc, res)
res[names(predict(firststage)), "fs"] <- stats::predict(firststage)
coefs <- rep(0, k) #k=4 in case of GK
names(coefs) <- seriesnames #the list is named
for (i in 1:k) {
  s <- seriesnames[i]
  if (s != "gs1") {
    secondstage <- stats::lm(stats::as.formula(paste(s, " ~ fs", sep = "")), res)
    coefs[i] <- secondstage$coefficients["fs"]
  } else {
    coefs[i] <- 1
  }
} #coef are S21/s11
coefs
s21_on_s11 <- matrix(coefs[2:k], c(k-1,1)) #3x1 matrix of coefficient (ratios)
Q <- (s21_on_s11 * gamma_11) %*% t(s21_on_s11) - (gamma_21 %*% t(s21_on_s11) + s21_on_s11 %*% t(gamma_21)) + gamma_22
s12s12 <- t(gamma_21 - s21_on_s11 * gamma_11) %*% solve(Q) %*% (gamma_21 - s21_on_s11 * gamma_11)
s11_squared <- gamma_11 - s12s12 #če ne bi mel s11 bi pomenil, da so to normalizirani koeficienti (če je s11 = 1)!
sp <- as.numeric(sqrt(s11_squared)) #s11 - note that this can be interpreted as standard deviation of the strucutral shock!!!!!
coefs <- coefs[origorder]
shockcolumn <- sp * coefs[origorder] #S21
shockcolumn
#---------------------------------IRFs-------------------------------------------
ma_representation <- Phi(var, 50)
ma_representation
irfs <- apply(ma_representation, 3, function(x) x %*% shockcolumn) #takes the MA coefficients and multiplies shockcol with each of the matrices. 
#The 3 is there because ma_representation is a three-dimensional array, this means that we have stacked matrices in an object.
irfs <- as.data.frame(t(irfs))
irfs#makes a dataframe from a matrix and transposes the matrix for 51 horizons
colnames(irfs) <- names(shockcolumn) #adds names to irfs
irfs <- mutate(irfs, horizon = 0:50) #adds horizon variable to irfs dataframe
irfs <- gather(irfs, key = variable, value = response, -horizon) #reshapse the dataframe, we get a long format, where
#horizon stays the same, a new column variable hold the names of variable in original dataframe and values are the ones that pertained to the each variable.
ggplot(irfs, aes(x = horizon, y = response, group = variable, color = variable)) + geom_line()

#Note that VARX does not allow you to identify absolute responses to underlying true structural shocks, only relative. These are the same as in the case of Proxy-VAR. 
#What you get when estimating the B matrix you get the raw structural coefficients that relate to how proxy impacts the dependent variables. You say unit shock in the proxy or one standard deviation of the proxy.
#If you normalize the policy rate to rise by 100bps you actually get the true monetary policy shock, but conditional on the set response in the rate!!!! 
#Absolute responses of variable to a structural underlying shock (not proxy) are only possible to get in a proxy VAR setting cause you have the explicit mapping!!!
