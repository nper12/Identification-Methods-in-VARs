library(vars)
library(tidyr)
library(dplyr)
library(ggplot2)
library(writexl)
library(readxl)
library(xlsx)
library(extraDistr)
library(lpirfs)
#------------------------------ Data Wrangling ---------------------------------
GKdata <- read_excel("/Users/nejcperme/Desktop/GKdata.xlsx")
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
}
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
structural_shock <- res%*%solve(gamma)%*%t(shockcolumn) #tko se dobi strukturni šok vn!!! - thanks to the guy the myth the legend Kaezing
shockcolumn <- as.matrix(shockcolumn)
structural_shock <- as.data.frame(structural_shock)
shock <- firststage$fitted.values
structural_shock <- structural_shock[116:384, ]
GKdata_subset <- GKdata[128:396, ]
GKdata_subset <- GKdata_subset[,1:4]
GKdata_subset_cpi <- GKdata_subset[,1]
GKdata_subset_cpi <- as.data.frame(GKdata_subset_cpi)
#GKdata_subset$structural_shock <- structural_shock$V1
check <- lm(logip ~ structural_shock,data= GKdata_subset)
summary(check)
lp <- lpirfs::lp_lin_iv(endog_data=GKdata_subset, shock=structural_shock, lags_endog_lin = 2, lags_exog = 2, hor=50, trend=0, confint=1.96, use_nw = TRUE)
plot(lp)

#različno od ramey ker ona sam uzame ff4_tc, kar je pr ns bascially nečisti proxy za šok. Strukturne šoke potem dobiš posebi vn!
#-------------------------- WILD BOOTSTRAP ----------------------
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
colnames(irfs) <- names(shockcolumn_boot)
colnames(lower_band) <- names(shockcolumn_boot)
colnames(upper_band) <- names(shockcolumn_boot)
lower_band <- as.data.frame(lower_band)
upper_band <- as.data.frame(upper_band)
irfs <- as.data.frame(irfs)
horizon <- c(0:50)
horizon <- as.data.frame(horizon)
lower_band <- cbind(horizon, lower_band)
upper_band <- cbind(horizon, upper_band)
irfs <- cbind(horizon, irfs)
combined_data <- merge(irfs, lower_band, by = "horizon", suffixes = c("_irf", "_lower"))
combined_data <- merge(combined_data, upper_band, by = "horizon", suffixes = c("", "_upper"))
combined_data <- combined_data[,-c(10,11,12,13)]

# Convert combined_data from wide to long format
long_data <- combined_data %>%
  pivot_longer(
    cols = -horizon, 
    names_to = c("name", ".value"),  # Use "name" for variable names, ".value" for response types
    names_sep = "_"  # Separate based on underscore
  )

# Check the structure of long_data
str(long_data)

# Plot using ggplot2 with facets
library(ggplot2)

ggplot(long_data, aes(x = horizon)) +
  geom_line(aes(y = irf, color = name), size = 1) +  # Use 'name' for line colors
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = name), alpha = 0.5) +  # Use 'name' for fill colors
  facet_wrap(~ name, scales = "free_y") +  # Separate plot for each variable
  theme_minimal() +
  labs(title = "Impulse Response Functions",
       x = "Horizon",
       y = "Response") +
  scale_color_manual(values = c("black", "black", "black", "black")) +  # Customize colors
  scale_fill_manual(values = c("black", "black", "black", "black")) +    # Fill colors for ribbons
  theme(legend.position = "bottom")


