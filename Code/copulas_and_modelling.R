library(quantmod)
library(tidyverse)
library(MCMCpack)
library(gridExtra)
library(VineCopula)
library(nleqslv)
library(forecast)
library(tseries)
library(FinTS)

### ====================== Data Prep ====================== ###
RNGversion("4.3.0")
set.seed(2026)

start_date <- as.Date("2016-06-01")
end_date <- as.Date("2026-06-06")

tickers <- c("NVDA", "AMD")
getSymbols(tickers, src = "yahoo", from = start_date, to = end_date, auto.assign
           = TRUE)
Adj_NVDA <- Ad(NVDA)
Adj_AMD <- Ad(AMD)
X1 <- (log(Adj_NVDA) - log(lag(Adj_NVDA, 1)))[-1]
X1 <- as.numeric(coredata(X1))
X2 <- (log(Adj_AMD) - log(lag(Adj_AMD, 1)))[-1]
X2 <- as.numeric(coredata(X2))

L1_NVDA = 1 - exp(X1)
L2_AMD = 1 - exp(X2)

quantiles = c(0.95, 0.99, 0.995)

n = 10000

### ====================== VAR and ES ====================== ###


VAR = function(loss, alpha) {
  return(as.numeric(quantile(loss, probs = alpha)))
}

ES_alpha = function(loss, alpha) {
  out = numeric(3)
  VaR = as.numeric(quantile(loss, probs = alpha))
  for (i in 1:3) {
    out[i] = mean(loss[loss > VaR[i]])
  }
  return(out)
}

### ====================== Copula Tasks ====================== ###


### Task 1 ###
RNGversion("4.3.0")
set.seed(2026)

VAR_ES_normal = data.frame(Alpha = quantiles, L1_VAR = VAR(L1_NVDA, quantiles), L2_VAR = VAR(L2_AMD, quantiles), 
           L1_ES = ES_alpha(L1_NVDA, quantiles), L2_ES = ES_alpha(L2_AMD, quantiles))

write.csv(round(VAR_ES_normal, 4), file = 'Results/normal_risk_metrics.csv')

### Task 2 ###
RNGversion("4.3.0")
set.seed(2026)

kendals_tau = function(l1, l2) {
  concordance_count = 0
  discordance_count = 0
  n = length(L1_NVDA)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      prod = (l1[i] - l1[j]) * (l2[i] - l2[j])
      if (prod > 0) {
        concordance_count  = concordance_count + 1
      } else if (prod < 0) {
        discordance_count = discordance_count + 1
      }
    }
  }
  return((concordance_count / choose(n, 2)) - (discordance_count / choose(n, 2)))
}


spearmans_rho <- function(l1, l2) {
  cdf_l1 <- ecdf(l1)(l1)
  cdf_l2 <- ecdf(l2)(l2)
  return(cor(cdf_l1, cdf_l2, method = "pearson"))
}

kendals_tau(L1_NVDA, L2_AMD)
spearmans_rho(L1_NVDA, L2_AMD)



### Task 3: Independence Copula ###

RNGversion("4.3.0")
set.seed(2026)

uniforms = matrix(c(runif(2*n)), nrow = n, ncol = 2)

X1_ind = quantile(X1, probs = uniforms[,1], names = FALSE)
X2_ind = quantile(X2, probs = uniforms[,2], names = FALSE)

L_ind = 0.5 * (1 - exp(X1_ind)) + 0.5 * (1 - exp(X2_ind))

independence_copula_summary = data.frame(Alpha = quantiles, Value_At_Risk = VAR(L_ind, quantiles), Expected_Shortfall = ES_alpha(L_ind, quantiles))

write.csv(round(independence_copula_summary,4), file = 'Results/independence_copula_eq_port.csv')

### Task 4: Gaussian Copula ###

RNGversion("4.3.0")
set.seed(2026)

#X~N_(0, rho) -> X = sqrt(rho) * Z
rho_hat = cor(X1, X2, method = 'pearson')
correlation_matrix = matrix(c(1,rho_hat, rho_hat, 1), ncol = 2, byrow = TRUE)

Z0 <- matrix(rnorm(2* n), nrow = n, ncol = 2)
A = t(chol(correlation_matrix))

normals = data.frame(Z1 = numeric(n), Z2 = numeric(n))
for (i in 1:n) {
  x = A %*% matrix(Z0[i, ], ncol = 1)
  normals[i, 1] = x[1,1]
  normals[i, 2] = x[2,1]
}

U1 <- pnorm(normals$Z1)
U2 <- pnorm(normals$Z2)

ggplot(data.frame(U1 = U1, U2 = U2), aes(x = U1, y = U2)) +
  geom_point(size = 0.5, colour = 'blue') + 
  theme_minimal() + 
  labs(title = 'Figure 1: U1 vs U2 (Gaussian Copula)')


X1_gaussian = as.numeric(quantile(X1, probs = U1))
X2_gaussian = as.numeric(quantile(X2, probs = U2))

L_gaussian = 0.5 * (1 - exp(X1_gaussian)) + 0.5 * (1 - exp(X2_gaussian))

gaussian_summary = data.frame(Alpha = quantiles, Value_At_Risk = VAR(L_gaussian, quantiles), Expected_Shortfall = ES_alpha(L_gaussian, quantiles))
write.csv(round(gaussian_summary,4), file = 'Results/gaussian_summary.csv')

### Task 5: T Copula ###
RNGversion("4.3.0")
set.seed(2026)

#Z = sqrt(W) * A * Z where A = sigma^1/2
rho_hat <- cor(X1, X2, method = "pearson")
sigma <- matrix(c(1, rho_hat, rho_hat, 1), nrow = 2, byrow = TRUE)
A <- t(chol(sigma))

generate_W <- function(v) {
  v / rchisq(n, df = v)
}

generate_T_copula <- function(v) {
  W <- generate_W(v)
  Z <- matrix(rnorm(2 * n), ncol = 2)
  Y <- Z %*% t(A)
  T_var <- Y * sqrt(W)
  U <- data.frame(
    U1 = pt(T_var[, 1], df = v),
    U2 = pt(T_var[, 2], df = v)
  )
  
  return(U)
}

RNGversion("4.3.0")
set.seed(2026)

u_3df = generate_T_copula(3)

RNGversion("4.3.0")
set.seed(2026)

u_10df = generate_T_copula(10)

RNGversion("4.3.0")
set.seed(2026)

u_10000df = generate_T_copula(10000)

plt1 = ggplot(u_3df, aes(x = U1, y = U2)) + 
  geom_point(color = 'blue', size = 0.5) + 
  theme_minimal() + 
  labs(title = 'T Copula 3 df: U1 vs U2')

plt2 = ggplot(u_10df, aes(x = U1, y = U2)) + 
  geom_point(color = 'brown', size = 0.5) + 
  theme_minimal() + 
  labs(title = 'T Copula 10 df: U1 vs U2')

plt3 = ggplot(u_10000df, aes(x = U1, y = U2)) + 
  geom_point(color = 'orange', size = 0.5) + 
  theme_minimal() + 
  labs(title = 'T Copula 10000 df: U1 vs U2')

grid.arrange(plt1, plt2, plt3, layout_matrix = rbind(c(1,2), c(3,3)), top = "Figure 2: U1 vs U2 (T-Copula)")

X1_3df = quantile(X1, probs = u_3df[, 1], names = FALSE) 
X1_10df = quantile(X1, probs = u_10df[, 1] , names = FALSE) 
X1_10000df = quantile(X1, probs = u_10000df[, 1], names = FALSE) 

X2_3df = quantile(X2, probs = u_3df[, 2], names = FALSE) 
X2_10df = quantile(X2, probs = u_10df[, 2] , names = FALSE) 
X2_10000df = quantile(X2, probs = u_10000df[, 2], names = FALSE) 

L_3df = 0.5 * (1 - exp(X1_3df)) + 0.5 * (1 - exp(X2_3df))
L_10df = 0.5 * (1 - exp(X1_10df)) + 0.5 * (1 - exp(X2_10df))
L_10000df = 0.5 * (1 - exp(X1_10000df)) + 0.5 * (1 - exp(X2_10000df))

T_copula_summary = data.frame(Alpha = quantiles, VAR_3df = VAR(L_3df, quantiles), VAR_10df = VAR(L_10df, quantiles), VAR_10000df = VAR(L_10000df, quantiles),
           ES_3df = ES_alpha(L_3df, quantiles), ES_10df = ES_alpha(L_10df, quantiles), ES_10000df = ES_alpha(L_10000df, quantiles))

write.csv(round(T_copula_summary,4), file = 'Results/tcopula_summary.csv')

### Task 6: Comonotonicity and Countermonotonicity ###

RNGversion("4.3.0")
set.seed(2026)

U_s = runif(n)

Comonotonic_Copula = data.frame(U1 = U_s, U2 = U_s)

plt4 = ggplot(Comonotonic_Copula, aes(x = U1, y = U2)) + 
  geom_point(color = 'orange', size = 0.5) + 
  theme_minimal() + 
  labs(title = 'Comonotonic Copula: U1 vs U2')


Countermonotonic_Copula = data.frame(U1 = U_s, U2 = 1 - U_s)

plt5 = ggplot(Countermonotonic_Copula, aes(x = U1, y = U2)) + 
  geom_point(color = 'steelblue', size = 0.5) + 
  theme_minimal() + 
  labs(title = 'Countermonotonic Copula: U1 vs U2')


X1_comono = quantile(X1, probs = Comonotonic_Copula[,1], names = FALSE)
X2_comono = quantile(X2, probs = Comonotonic_Copula[,2], names = FALSE)

X1_countermono = quantile(X1, probs = Countermonotonic_Copula[,1], names = FALSE)
X2_countermono = quantile(X2, probs = Countermonotonic_Copula[,2], names = FALSE)

L_comono = 0.5 * (1 - exp(X1_comono)) + 0.5 * (1 - exp(X2_comono))
L_countermono = 0.5 * (1 - exp(X1_countermono)) + 0.5 * (1 - exp(X2_countermono))

mono_summary = data.frame(Alpha = quantiles, VAR_comono = VAR(L_comono, quantiles) , VAR_countermono = VAR(L_countermono, quantiles), 
                          ES_comono = ES_alpha(L_comono, quantiles) , ES_countermono = ES_alpha(L_countermono, quantiles))


write.csv(round(mono_summary,4), file = 'Results/mono_summary.csv')

### Task 7 ###

RNGversion("4.3.0")
set.seed(2026)

probability = function(a, b) {
  out = numeric(3)
  u = ecdf(a)(a)
  v = ecdf(b)(b)
  for (i in 1:3) {
    alpha = quantiles[i]
    out[i] = sum((u <= 1 - alpha) & (v <= 1 - alpha)) / sum(v <= 1 - alpha)
  }
  return(out)
}


prob_table = data.frame(original = probability(X1, X2), 
                        independence = probability(X1_ind, X2_ind), 
                        gaussian = probability(X1_gaussian, X2_gaussian),
                        t3df = probability(X1_3df, X2_3df), 
                        t10df = probability(X1_10df, X2_10df), 
                        t10000df = probability(X1_10000df, X2_10000df),
                        comono = probability(X1_comono, X2_comono), 
                        countermono = probability(X1_countermono, X2_countermono))  

write.csv(round(prob_table, 4), file = 'Results/prob_metrics.csv')


L_original = 0.5 * L1_NVDA + 0.5*L2_AMD
original_df = data.frame(Percentiles = quantiles, VaR = VAR(L_original, quantiles), ES = ES_alpha(L_original, quantiles))
write.csv(round(original_df, 4), file = 'Results/truevals.csv')

### Task 8: BiCopSelect ###

RNGversion("4.3.0")
set.seed(2026)

test = BiCopSelect(ecdf(X1)(X1), ecdf(X2)(X2), selectioncrit = 'logLik', method = 'mle')

summary(test)
test$familyname
test$par
test$par2
test$logLik
test$AIC
test$BIC
test$tau


### ====================== Time Series Tasks ====================== ###


### Task 1: Time Series and ACF Plots ###

RNGversion("4.3.0")
set.seed(2026)

S_t_data = data.frame(NVDA_St = nvda_adj, AMD_St = amd_adj)
S_t_long = S_t_data %>%
  pivot_longer(cols = everything(), 
               names_to = 'Stock',
               values_to = 'Prices') %>%
  group_by(Stock) %>%
  mutate(t = row_number())

S_t = ggplot(S_t_long, aes(x = t, y = Prices, color = Stock)) +
  geom_line(linewidth = 1, alpha = 0.6) +
  theme_minimal() + 
  labs(y = 'S_t', x = 't', title = 'NVDA and AMD Stock Price Processes')


X_t_data = data.frame(NVDA_Xt = X1, AMD_Xt = X2)
X_t_long = X_t_data %>%
  pivot_longer(cols = everything(), 
               names_to = 'Stock',
               values_to = 'Prices') %>%
  group_by(Stock) %>%
  mutate(t = row_number())

X_t = ggplot(X_t_long, aes(x = t, y = Prices, color = Stock)) +
  geom_line(linewidth = 1, alpha = 0.6) +
  theme_minimal() + 
  labs(y = 'X_t', x = 't', title = 'NVDA and AMD Log Return Processes')

grid.arrange(S_t, X_t, layout_matrix = rbind(c(1,1), c(2,2)), top = 'Figure 3: NVDA and AMD Processes')


acf_extended = function(x, up_to) {
  return(as.vector(acf(x, lag.max = up_to, plot = FALSE)$acf))
}


nvda_acf_data = data.frame(X_i = acf_extended(X1, up_to = 15), S_i = acf_extended(nvda_adj, up_to = 15))
nvda_acf_data_long = nvda_acf_data %>%
  pivot_longer(cols = everything(),
               names_to = 'Process',
               values_to = 'Values') %>%
  group_by(Process) %>%
  mutate(t = row_number() - 1)

plt6 = ggplot(nvda_acf_data_long,
       aes(x = t, y = Values, color = Process)) +
  geom_line() +
  geom_point() + 
  theme_minimal() + 
  labs(title = 'ACF of NVDA Stock Price and Log Return Processes',
       y = 'ACF')


amd_acf_data = data.frame(X_i = acf_extended(X2, up_to = 15), S_i = acf_extended(amd_adj, up_to = 15))
amd_acf_data_long = amd_acf_data %>%
  pivot_longer(cols = everything(),
               names_to = 'Process',
               values_to = 'ACF') %>%
  group_by(Process) %>%
  mutate(t = row_number() - 1)

plt7 = ggplot(amd_acf_data_long,
       aes(x = t, y = ACF, color = Process)) +
  geom_line() +
  geom_point() + 
  theme_minimal() + 
  labs(title = 'ACF of AMD Stock Price and Log Return Processes')

grid.arrange(plt6, plt7, layout_matrix = rbind(c(1,1), c(2,2)), top = 'Figure 4: ACF of NVIDIA and AMD Processes')

### Task 2: Moving Average Process of Xt order of 3 ###

RNGversion("4.3.0")
set.seed(2026)

#Cov(X_t) = sigma(1 + theta_1^2 + theta_2^2 + theta_3^2)
#Cov(X_t, X_t - 1) = sigma(theta_1^2 + theta_1*theta_2 + theta_2 * theta_3)
#Cov(X_t,X_t-2) = σε^2(θ2+θ1θ3)
#Cov(X_t,X_t-3) = σε^2 x θ3,

#Sample Mean: Mu



ACF_X1_up_to_3 = nvda_acf_data[2:4 ,1]
ACF_X2_up_to_3 = amd_acf_data[2:4, 1]

#z = (theta_1, theta_2, theta_3)
RNGversion("4.3.0")
set.seed(2026)

f_x1 = function(z) {
  theta_1 = z[1]
  theta_2 = z[2]
  theta_3 = z[3]
  
  denominator = 1 + theta_1^2 + theta_2^2 + theta_3^2
  return(c(
    (theta_1 + theta_1 * theta_2 + theta_2 * theta_3) /denominator - ACF_X1_up_to_3[1],
    (theta_2 + theta_1 * theta_3) / denominator -  ACF_X1_up_to_3[2],
    theta_3 / denominator - ACF_X1_up_to_3[3]
  ))
}
  
f_x2 = function(z) {
  theta_1 = z[1]
  theta_2 = z[2]
  theta_3 = z[3]
  
  denominator = 1 + theta_1^2 + theta_2^2 + theta_3^2
  return(c(
    (theta_1 + theta_1 * theta_2 + theta_2 * theta_3) /denominator - ACF_X2_up_to_3[1],
    (theta_2 + theta_1 * theta_3) / denominator -  ACF_X2_up_to_3[2],
    theta_3 / denominator - ACF_X2_up_to_3[3]
  ))
}

test_invertibility <- function(theta) {
  poly_coeff <- c(theta[3], theta[2], theta[1], 1)
  roots <- polyroot(poly_coeff)
  root_modulus <- abs(roots)
  if (all(root_modulus > 1)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

starts <- data.frame(
  theta1 = runif(1000, -2, 2),
  theta2 = runif(1000, -2, 2),
  theta3 = runif(1000, -2, 2)
)  

solutions_x1 = data.frame(NA, ncol = 3, nrow = 1000)
colnames(solutions_x1) = c('theta_1', 'theta_2', 'theta_3')

for (i in 1:1000) {
  solutions_x1[i, ] = nleqslv(as.numeric(starts[i, ]), f_x1)$x
}

invertable_x1 = logical(1000)
for (i in 1:1000) {
  theta = as.numeric(solutions_x1[i, ])
  invertable_x1[i] = test_invertibility(theta = theta)
}

solutions_x1 = solutions_x1[invertable_x1, ] %>%
  mutate(euclidean = sqrt(theta_1^2 + theta_2^2 + theta_3^2))

thetas_x1 = solutions_x1[which.min(solutions_x1$euclidean), ]

solutions_x2 = data.frame(NA, ncol = 3, nrow = 1000)
colnames(solutions_x2) = c('theta_1', 'theta_2', 'theta_3')

for (i in 1:1000) {
  solutions_x2[i, ] = nleqslv(as.numeric(starts[i, ]), f_x2)$x
}

invertable_x2 = logical(1000)
for (i in 1:1000) {
  theta = as.numeric(solutions_x2[i, ])
  invertable_x2[i] = test_invertibility(theta = theta)
}

solutions_x2 <- solutions_x2[invertable_x2, ] %>%
  mutate(euclidean = sqrt(theta_1^2 + theta_2^2 + theta_3^2))

thetas_x2 <- solutions_x2[which.min(solutions_x2$euclidean), ]


theta_1_x1 = as.numeric(thetas_x1[1])
theta_2_x1 = as.numeric(thetas_x1[2])
theta_3_x1 = as.numeric(thetas_x1[3])


theta_1_x2 = as.numeric(thetas_x2[1])
theta_2_x2 = as.numeric(thetas_x2[2])
theta_3_x2 = as.numeric(thetas_x2[3])


sigma_x1 = sqrt(var(X1) / (1 + theta_1_x1 ^ 2 + theta_2_x1 ^ 2 + theta_3_x1 ^ 2))
sigma_x2 = sqrt(var(X2) / (1 + theta_1_x2 ^ 2 + theta_2_x2 ^ 2 + theta_3_x2 ^ 2))

mu_x1 = mean(X1)
mu_x2 = mean(X2)


parameter_estimation_results <- data.frame(
  Series = c("X1", "X2"),
  mu = c(mu_x1, mu_x2),
  sigma = c(sigma_x1, sigma_x2),
  theta_1 = c(theta_1_x1, theta_1_x2),
  theta_2 = c(theta_2_x1, theta_2_x2),
  theta_3 = c(theta_3_x1, theta_3_x2)
)

write.csv(parameter_estimation_results, file = 'Results/parameter_estimation.csv')

test_invertibility(theta = as.numeric(parameter_estimation_results[1, 4:6])) #Double checking invertability
test_invertibility(theta = as.numeric(parameter_estimation_results[2, 4:6]))


### Task 3 ###

RNGversion("4.3.0")
set.seed(2026)

fit_1 = auto.arima(X1, ic = 'aic')
fit_2 = auto.arima(X2, ic = 'aic')

coef(fit_1)
coef(fit_2)

### Task 4: Box Jenkins Predict ###

RNGversion("4.3.0")
set.seed(2026)

predict_x1 = function() {
  #MA(2) = MU + Z_t + THETA_1*Z_t-1 + THETA_2 * Z_t-2
  #X_t+1 = mu + theta_1*z_t + theta_2 * z_t-1
  #X_t+2 = mu + theta_2 * z_t
  #X_t+3 = mu
  #X_t+4 = mu
  #X_t+5 = mu
  mu = as.numeric(coef(fit_1)['intercept'])
  theta_1 = as.numeric(coef(fit_1)['ma1'])
  theta_2 = as.numeric(coef(fit_1)['ma2'])
  z = resid(fit_1)
  n = length(z)
  z_t = z[n]
  z_tsub1 = z[n-1]

  return(c(mu + theta_1 * z_t + theta_2 * z_tsub1,
                   mu + theta_2 * z_t, mu, mu, mu))
}


predict_x2 = function() {
  #ARMA(2,2) = mu + psi_1*x_t-1 + psi_2 * x_t-2 + z_t + theta_1 * z_t-1 + theta_2 * z_t-2
  #X_t+1 = mu + phi_1*x_t + phi_2 * x_t-1 + theta_1 * z_t + theta_2 * z_t-1
  
  z = resid(fit_2)
  n = length(z)
  z_t = z[n]
  z_tsub1 = z[n-1]
  mu = as.numeric(coef(fit_2)['intercept'])
  theta_1 = as.numeric(coef(fit_2)['ma1'])
  theta_2 = as.numeric(coef(fit_2)['ma2'])
  phi_1 = as.numeric(coef(fit_2)['ar1'])
  phi_2 = as.numeric(coef(fit_2)['ar2'])
  x_t = X1[length(X1)]
  x_tsub1 = X1[length(X1)-1]

  x_tilda1 = mu + phi_1 * x_t + phi_2 * x_tsub1 + theta_1 * z_t + theta_2*z_tsub1
  x_tilda2 = mu + phi_1 * x_tilda1 + phi_2 * x_t + theta_2 * z_t
  x_tilda3 = mu + phi_1 * x_tilda2 + phi_2 *x_tilda1
  x_tilda4 = mu + phi_1 * x_tilda3 + phi_2 * x_tilda2
  x_tilda5 = mu + phi_1 * x_tilda4 + phi_2 * x_tilda3
  
  return(c(x_tilda1, x_tilda2, x_tilda3, x_tilda4, x_tilda5))
}

 
x1_results = data.frame(time = c(1,2,3,4,5), expected_nvda_logreturn = predict_x1())
x2_results = data.frame(time = c(1,2,3,4,5), expected_amd_logreturn = predict_x2())


# actual stock prices (note some days are not available, so i used the nearest 5 available days of stock data)

start_new = as.Date('2026-06-05')
end_new = as.Date('2026-06-15')


getSymbols(tickers, src = "yahoo", from = start_new, to = end_new, auto.assign
           = TRUE)

future_nvda = Ad(NVDA)
future_X1 = as.vector(na.omit(diff(log(future_nvda$NVDA.Adjusted))))

future_amd = Ad(AMD)
future_X2 = as.vector(na.omit(diff(log(future_amd$AMD.Adjusted))))


mse = function(predicted, actual) {
  return(mean((actual - predicted)^2))
}

absolute_error = function(predicted, actual) {
  return(mean(abs(actual - predicted)))
}


# comparison and test metrics

x1_results['actual'] = future_X1 
x1_results_long = x1_results %>%
  pivot_longer(cols = c('expected_nvda_logreturn', 'actual'), values_to = 'log_return', names_to = 'Actual_and_Expected')

plt8 = ggplot(x1_results_long, aes(x = time, y = log_return, color = Actual_and_Expected)) + 
  geom_line() + 
  geom_point(color = 'black') + 
  theme_minimal() + 
  labs(y = 'Log Return', title = 'Actual NVDA vs Box-Jenkins Predicted Log Returns') +
  scale_color_manual(
    labels = c(
      "expected_nvda_logreturn" = "Predicted",
      "actual" = "Actual"
    ),
    values = c(
      "expected_nvda_logreturn" = "blue",
      "actual" = "red"
    )
  ) 



x2_results['actual'] = future_X2 

x2_results_long = x2_results %>%
  pivot_longer(cols = c('expected_amd_logreturn', 'actual'), values_to = 'log_return', names_to = 'Actual_and_Expected')
plt9 = ggplot(x2_results_long, aes(x = time, y = log_return, color = Actual_and_Expected)) + 
  geom_line() + 
  geom_point(color = 'black') + 
  theme_minimal() + 
  labs(y = 'Log Return', title = 'Actual AMD vs Box-Jenkins Predicted Log Returns') +
  scale_color_manual(
    labels = c(
      "expected_amd_logreturn" = "Predicted",
      "actual" = "Actual"
    ),
    values = c(
      "expected_amd_logreturn" = "blue",
      "actual" = "red"
    )
  ) 

grid.arrange(plt8, plt9, layout_matrix = rbind(c(1,1), c(2,2)), top = 'Figure 5: Prediction Estimates vs Actual Experience')


prediction_results <- data.frame(
  Stock = c("NVDA", "AMD"),
  MeanSquaredError = c(
    mse(
      predicted = x1_results$expected_nvda_logreturn,
      actual = x1_results$actual
    ),
    mse(
      predicted = x2_results$expected_amd_logreturn,
      actual = x2_results$actual
    )
  ),
  AbsoluteError = c(
    absolute_error(
      predicted = x1_results$expected_nvda_logreturn,
      actual = x1_results$actual
    ),
    absolute_error(
      predicted = x2_results$expected_amd_logreturn,
      actual = x2_results$actual
    )
  )
)

prediction_results
write.csv(prediction_results, file = 'Results/prediction_results.csv')


### Task 5: Testing Goodness of Fit ###

RNGversion("4.3.0")
set.seed(2026)

checkresiduals(fit_1)
checkresiduals(fit_2)

#Null: residuals are not autocorrelated
#Alternative Residuals: at least one autocorrelation is non zero
Box.test(
  resid(fit_1),
  lag = 15,
  type = "Ljung-Box",
  fitdf = length(coef(fit_1))
)

#The residuals from your NVDA ARIMA(0,0,2) model do not behave like white noise. 
#There is evidence that some autocorrelation remains in the residuals, 
#suggesting that the model has not fully captured the time-series dependence structure.

Box.test(
  resid(fit_2),
  lag = 15,
  type = "Ljung-Box",
  fitdf = length(coef(fit_2))
)

#The residuals from the AMD ARIMA(2,0,2) model are consistent with a white noise process.
#There is no significant evidence of remaining autocorrelation, 
#suggesting that the model has adequately captured the dependence structure in the AMD log-return series.


ArchTest(resid(fit_1), lags = 12)
ArchTest(resid(fit_2), lags = 12)

#both models exhibit time-varying conditional variance,

jarque.bera.test(resid(fit_1))
jarque.bera.test(resid(fit_2))

#Evidence that are not normal

t.test(resid(fit_1), mu = 0)
t.test(resid(fit_2), mu = 0)

#Not fully consistent with white noise, not normal, residuals are correlated, heteroskedasticity but do have zero mean



