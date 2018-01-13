library(ggplot2)
library(lubridate)
library(fpp2)
library(reshape2)
library(rstan)
library(BRRR)
library(bridgesampling)
set.seed(123)

# Load data
setwd("/Users/skyehersh/Documents/fall_2017/stats/")

# Original data
df.raw <- read.csv('mauna_loa/co2_mlo.csv', header = TRUE)

# Data w/ interpolations + future timesteps 
df.new <- read.csv('mauna_loa/final_co2.csv', header = TRUE)
df.orig <- df.new[1:3093, ]
df.futr <- df.new[3094:nrow(df.new), ]

# Time series for visual analysis. I start the time 
# series at 1959 (idx 26) b/c many vals missing 1958.

ppm.time.series <- ts(df.orig[26:nrow(df.orig), ]$ppm, 
                      start = decimal_date(ymd(df.orig$date[26])), 
                      frequency = 365.25/7)

autoplot(ppm.time.series) +
  labs(title = expression(paste("Atmospheric CO"[2]," concentration")), 
       subtitle = "1958 to present",
       x = "Year", 
       y = expression(paste("CO"[2], " concentration (ppm)")))  + 
  theme(text = element_text(family = 'serif'))

ggsave('original_data.png')

print("Number of datapoints per year:")
table(year(df.raw$date))

### 1961 plot
qplot(df.raw$day[127:178], df.raw$ppm[127:178]) +
  geom_point(color = "turquoise3") +
  labs(title = expression(paste("Atmospheric CO"[2]," concentration, 1961")), 
       x = "Date in 1961 (m/d)", 
       y = expression(paste("CO"[2], " concentration (ppm)"))) + 
  scale_x_continuous(labels = c('1/7', '4/1', '7/1', '9/30', '12/30'), 
                     breaks = c(1015, 1099, 1190, 1281, 1365)) +
  theme(text = element_text(family = 'serif'),
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))

ggsave('1961.png')

### 1964

# Binary vector indicating which
# dates in 1964 were interpolated
interp <- c(rep("Original data", 3), 
            rep("Interpolated data", 18), 
            rep("Original data", 2), 
            rep("Interpolated data", 2), 
            rep("Original data", 26))

# Dataframe for 1964
temp <- data.frame(interp, 
                   cbind(df.orig$day[286:336], 
                         df.orig$ppm[286:336]))

ggplot(data = temp, aes(X1, X2)) +
  geom_point(aes(color = interp)) + 
  labs(title = expression(paste("Atmospheric CO"[2]," concentration, 1964")), 
       x = "Date in 1964 (m/d)", 
       y = expression(paste("CO"[2], " concentration (ppm)"))) + 
  scale_x_continuous(labels = c('1/4', '3/28', '6/27', '9/26', '12/26'), 
                     breaks = c(2107, 2191, 2282, 2373, 2464)) +
  theme(text = element_text(family = 'serif'), 
        legend.title = element_blank(),
        legend.position = c(0.8, 0.8),
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))

rm(temp)

ggsave('1964.png')

print("Data points per year, after interpolating and generating future dates")
table(year(df.new$date))

# Pick idx at which to split into train / test
split  <- floor((nrow(df.orig)*0.8))

# Split into training and testing sets
train  <- df.orig[1:split, ]
test   <- df.orig[(split + 1):nrow(df.orig), ]

autoplot(stl(ppm.time.series, s.window=365.25/7, robust=TRUE))
ggsave('components.png')

### Model decomposition

trendpattern = filter(df.orig$scaled_ppm, 
                      filter = rep(1/52.18, 52.18)
)

seasonals = df.orig$scaled_ppm - trendpattern

components <- data.frame(1:nrow(df.orig), 
                         trendpattern, 
                         df.orig$scaled_ppm, 
                         seasonals)

names(components) <- c('x', 
                       'trendpattern', 
                       'scaled.ppm', 
                       'seasonals')

# Indices for formatting x-ticks
breaks = floor(seq(1, nrow(df.orig), length.out = 10))

# Long term trend plot
ggplot(components, aes(x)) +
  geom_path(aes(y = scaled.ppm)) +
  geom_path(aes(y = trendpattern, col = 'Trend'), lwd = 1.3) +
  theme(text = element_text(family = 'serif'),
        legend.title = element_blank(),
        legend.position = c(0.9, 0.1)) +
  labs(title = expression(paste("Atmospheric CO"[2],", moving average annual trend")),
       x = "Year",
       y = expression(paste("CO"[2], " concentration (ppm), normalized"))) +
  scale_x_continuous(labels = year((df.orig$date))[breaks],
                     breaks = breaks)

ggsave('moving_averages.png')

### Polynomial fitting

linear  <- lm(train$scaled_ppm ~ train$scaled_day)
quad    <- lm(train$scaled_ppm ~ poly(train$scaled_day, 2))
cubic   <- lm(train$scaled_ppm ~ poly(train$scaled_day, 3))
quartic <- lm(train$scaled_ppm ~ poly(train$scaled_day, 4))

linear.pred   <- predict(linear, data.frame(x=train$scaled_day))
quad.pred     <- predict(quad, data.frame(x=train$scaled_day))
cubic.pred    <- predict(cubic, data.frame(x=train$scaled_day))
quartic.pred  <- predict(quartic, data.frame(x=train$scaled_day))

fits <- data.frame(train$scaled_day, 
                   train$scaled_ppm,
                   linear.pred,
                   quad.pred,
                   cubic.pred,
                   quartic.pred)

library(reshape2)
long.fits <- melt(fits, id = 'train.scaled_day')

# For x ticks
train.idxs <- seq(1, nrow(train), length.out = 10)
train.lbls <- year(train$date[train.idxs])
train.brks <- train$scaled_day[train.idxs]

ggplot(data = long.fits, aes(x = train.scaled_day, y = value, color = variable)) + 
  geom_line() +
  scale_colour_manual(values = c('black', 'green', 'blue', 'brown', 'hotpink'),
                      labels = c("Original data", 
                                 "Linear fit", 
                                 "Quadratic fit", 
                                 "Cubic fit", 
                                 "Quartic fit")) + 
  labs(title = "Polynomial fits on training data", 
       x = "Year", 
       y = expression(paste("CO"[2], " concentration (ppm), normalized"))) +
  scale_x_continuous(labels = train.lbls, 
                     breaks = train.brks) +
  theme(text = element_text(family = 'serif'),
        legend.title = element_blank(), 
        legend.position = c(0.85, 0.3))

ggsave('polynomial_fits.png')

### Peaks and troughs
mays <- which(month(df.orig$date) == 5)
septs <- which(month(df.orig$date) == 9)

# Indices to select x-tick labels
idxs <- floor(seq(1, 700, length.out = 10))

ggplot(data = df.orig, aes(day, ppm)) +
  geom_line() +
  coord_cartesian(xlim = c(1, 5000), ylim = c(300, 340)) +
  labs(title = expression(paste("CO"[2]," seasonality pattern")), 
       subtitle = "May colored pink, September colored blue",
       caption = "Closeup on 1958 to 1971",
       x = "Year", 
       y = expression(paste("CO"[2], " concentration (ppm)")))  + 
  scale_x_continuous(labels = year(df.orig$date)[idxs], 
                     breaks = df.orig$day[idxs]) +
  theme(text = element_text(family = 'serif')) +
  geom_vline(xintercept = df.orig$day[mays], lwd = 0.08, color = 'violetred3', alpha = 0.6) +
  geom_vline(xintercept = df.orig$day[septs], lwd = 0.08, color = 'blue', alpha = 0.6)

ggsave('peaks_troughs.png')

### Seasonal plotting
ggseasonplot(ppm.time.series, 
             ylab = expression(paste("CO"[2], " concentration (ppm), normalized"))) +
  labs(title = expression(paste("CO"[2], " concentration, seasonal plot"))) +
  scale_x_continuous(labels = c(1, 13, 26, 39, 52), breaks = c(0, 0.25, 0.5, 0.75, .98)) +
  theme(text = element_text(family = 'serif'))

ggsave('seasonal_plot.png')

# Seasonal trend
ggplot(components, aes(x)) +
  geom_path(aes(y = seasonals)) + 
  labs(title = expression(paste("Atmospheric CO"[2],", seasonal trend")), 
       x = "Year") +
  geom_hline(aes(yintercept = 0, col = "Long-term trend")) +
  theme(text = element_text(family = 'serif'),
        axis.title.y = element_blank(),
        legend.title = element_blank()) +
  scale_x_continuous(labels = year((df.orig$date))[breaks], 
                     breaks = breaks)

ggsave('seasonal_trend.png')

# Trend model prior sampling
n.samples = 8

c0 <- rnorm(n.samples, 0.02, 0.01)
c1 <- rnorm(n.samples, 0.45, 0.1)
c2 <- rnorm(n.samples, 0.5, 0.1)

# Make matrix in which each row will be a sample prior function
quad.samples <- matrix(rep(NA, nrow(train) * n.samples), 
                       ncol=nrow(train))

for(i in 1:n.samples) {
  quad.samples[i,] <- (c0[i] + 
                         c1[i] * train$scaled_day + 
                         c2[i] * train$scaled_day^2
  )
}


# Visualize the samples
quad.trends <- data.frame(t(rbind(train$scaled_day, train$scaled_ppm, quad.samples)))
names(quad.trends) <- c('day', 'ppm', seq(1, n.samples))

plt.quad <- melt(quad.trends, id = 'day')

ggplot(data = plt.quad, aes(x = day, y = value, color = variable)) +
  geom_line() +
  labs(title = "Long-term trend prior samples: quadratic model",
       subtitle = expression(paste("Samples of c"[0], " + c"[1], "t"[i], " + c"[2],
                                   "t"[i]^2)),
       x = "Time, normalized (t)", 
       y = expression(paste("CO"[2], " concentration (ppm), normalized"))) +
  theme(legend.position = "none",
        text = element_text(family = 'serif'))

ggsave('quad_prior_samples.png')

c0 <- tanh(rnorm(n.samples, -0.7, 0.04))
c1 <- 1 / (1 + exp(-rnorm(n.samples, 0.4, 0.05)))
c2 <- exp(rnorm(n.samples, 0.95, 0.05))

# Make matrix in which each row will be a sample prior function
exp.samples <- matrix(rep(NA, nrow(train) * n.samples), 
                      ncol=nrow(train))

for(i in 1:n.samples) {
  exp.samples[i,] <- (c0[i] + c1[i] * c2[i] ^ train$scaled_day)
}

# Visualize the samples
exp.trends <- data.frame(t(rbind(train$scaled_day, train$scaled_ppm, exp.samples)))
names(exp.trends) <- c('day', 'ppm', seq(1, n.samples))

plt.exp <- melt(exp.trends, id = 'day')

ggplot(data = plt.exp, aes(x = day, y = value, color = variable)) +
  geom_line() +
  labs(title = "Long-term trend prior samples: exponential model",
       subtitle = expression(paste("Samples of c"[0], " + c"[1],"c"[2]^"t")),
       x = "Time, normalized (t)", 
       y = expression(paste("CO"[2], " concentration (ppm), normalized"))) +
  theme(legend.position = "none",
        text = element_text(family = 'serif'))

ggsave('exp_prior_samples.png')

# Stan data (for both models)
stan.data <- list(
  N = nrow(train),                    # number training data points
  t = train$scaled_day,               # rescaled training dates/timesteps
  N_test = nrow(test),                # number test timesteps
  t_test = test$scaled_day,           # rescaled test timesteps
  N_future = nrow(df.futr),           # number of future timesteps
  t_future = df.futr$scaled_day,        # rescaled future dates/timesteps
  ppm = train$scaled_ppm,               # rescaled ppm
  k = (2 * pi * (tail(df.orig$day, 1) / # periodic frequency
                   365.25))
)

# THE QUADRATIC MODEL
quad.model = "
data {
int N;                   // number of training timesteps
real t[N];               // train timesteps
int N_test;              // number of test timesteps
real t_test[N_test];     // test timesteps
int N_future;            // number of future timesteps 
real t_future[N_future]; // future timesteps 
real ppm[N];             // observed ppm values
real k;                   // denominator of period
}

parameters {
real c0_prime;      // intercept
real c1_prime;      // linear component
real c2_prime;      // quadratic component
real A_prime;       // amplitude
real phi_prime;     // phase shift 
real sigma_prime;   // variance
}

transformed parameters {
real c0;
real c1;
real c2;
real A;
real phi;
real sigma;

c0 = c0_prime;      
c1 = c1_prime;   
c2 = c2_prime;                              
A = exp(A_prime);         
phi = phi_prime;                                        
sigma = exp(sigma_prime); 
}

model {
// Priors
c0_prime ~ normal(0.02, 0.01);   
c1_prime ~ normal(0.45, 0.1);
c2_prime ~ normal(0.5, 0.1);
A_prime ~ normal(-3.5, 0.5);     
phi_prime ~ normal(0, pi() / 2);           
sigma_prime ~ normal(-4.5, 1); 

// Likelihood 
for (n in 1:N){
ppm[n] ~ normal(c0 + 
c1 * t[n] + 
c2 * square(t[n]) + 
A  * cos((t[n] * k) + phi),
sigma); 
}
}

generated quantities {
// Posterior samples
real ppm_pred[N];            // posterior predicted, train set
real ppm_test_pred[N_test];  // posterior predicted, test set
real ppm_future[N_future];  // posterior predicted, proj. future

for (n in 1:N) {
ppm_pred[n] = normal_rng(c0 +                        // intercept
c1 * t[n] +                 // linear
c2 * square(t[n]) +         // quadratic
A  * cos((t[n] * k) + phi), // seasonal
sigma);                     // error
}

for (n in 1:N_test) {
ppm_test_pred[n] = normal_rng(c0 + 
c1 * t_test[n] + 
c2 * square(t_test[n]) +
A  * cos((t_test[n] * k) + phi), 
sigma);
}

for (n in 1:N_future) {
ppm_future[n] = normal_rng(c0 + 
c1 * t_future[n] + 
c2 * square(t_future[n]) + 
A  * cos((t_future[n] * k) + phi),
sigma);
}
}
"

# THE EXP MODEL
exp.model = "
data {
int N;                   // number of training timesteps
real t[N];               // train timesteps
int N_test;              // number of test timesteps
real t_test[N_test];     // test timesteps
int N_future;            // number of future timesteps 
real t_future[N_future]; // future timesteps 
real ppm[N];             // observed ppm values
real k;                   // denominator of period
}

parameters {
real c0_prime;      // intercept
real c1_prime;      // scaling constant
real c2_prime;      // exponentiated constant
real A_prime;       // amplitude
real phi_prime;     // phase shift 
real sigma_prime;   // variance
}

transformed parameters {
real c0;
real c1;
real c2;
real A;
real phi;
real sigma;

c0 = tanh(c0_prime);      
c1 = 1 / (1 + exp(-c1_prime));   
c2 = exp(c2_prime); 
A = exp(A_prime);         
phi = phi_prime;                                        
sigma = exp(sigma_prime); 
}

model {
// Priors
c0_prime ~ normal(-0.7, 0.01);   
c1_prime ~ normal(0.4, 0.05);
c2_prime ~ normal(0.95, 0.05);
A_prime ~ normal(-3.2, 5);     
phi_prime ~ normal(0, pi() / 2);           
sigma_prime ~ normal(-4.5, 1); 

// Likelihood 
for (n in 1:N){
ppm[n] ~ normal(c0 + 
c1 * pow(c2, t[n]) + 
A  * cos((t[n] * k) + phi),
sigma); 
}
}

generated quantities {
// Posterior samples
real ppm_pred[N];            // posterior predicted, train set
real ppm_test_pred[N_test];  // posterior predicted, test set
real ppm_future[N_future];   // posterior predicted, proj. future

for (n in 1:N) {
ppm_pred[n] = normal_rng(c0 +                        // intercept
c1 * pow(c2, t[n]) + 
A  * cos((t[n] * k) + phi), // seasonal
sigma);                     // error
}

for (n in 1:N_test) {
ppm_test_pred[n] = normal_rng(c0 + 
c1 * pow(c2, t_test[n]) +
A  * cos((t_test[n] * k) + phi), 
sigma);
}

for (n in 1:N_future) {
ppm_future[n] = normal_rng(c0 + 
c1 * pow(c2, t_future[n]) + 
A  * cos((t_future[n] * k) + phi),
sigma);
}
}
"

quad.fit <- stan(
  model_code = quad.model, 
  data = stan.data,  
  chains = 3,           
  warmup = 1000,     
  iter = 2000,     
  cores = 2,            
  refresh = 1000,         
  control = list(adapt_delta = 0.999)
); skrrrahh(4)

exp.fit <- stan(
  model_code = exp.model, 
  data = stan.data,  
  chains = 3,           
  warmup = 1000,     
  iter = 2000,     
  cores = 2,            
  refresh = 1000,         
  control = list(adapt_delta = 0.999)
); skrrrahh(4)

### Stan results
print("QUADRATIC MODEL:")
print(quad.fit, par=c('c0', 'c1', 'c2', 'A', 'phi', 'sigma'), probs=c(.05, 0.95))
quad.samples <- extract(quad.fit)

print("EXPONENTIAL MODEL:")
print(exp.fit, par=c('c0', 'c1', 'c2', 'A', 'phi', 'sigma'), 
      probs=c(.05, 0.95))
exp.samples <- extract(exp.fit)

# Autocorrelation plots
c0.acf <- acf(quad.samples$c0, plot = FALSE)
c0.acf.df <- with(c0.acf, data.frame(lag, acf))

ggplot(data = c0.acf.df, aes(lag, acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  theme(text = element_text(family = 'serif'))
ggsave('auto1.png')

c1.acf <- acf(quad.samples$c1, plot = FALSE)
c1.acf.df <- with(c1.acf, data.frame(lag, acf))

ggplot(data = c1.acf.df, aes(lag, acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  theme(text = element_text(family = 'serif'))
ggsave('auto2.png')

c2.acf <- acf(quad.samples$c2, plot = FALSE)
c2.acf.df <- with(c2.acf, data.frame(lag, acf))

ggplot(data = c2.acf.df, aes(lag, acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  theme(text = element_text(family = 'serif'))
ggsave('auto3.png')

A.acf <- acf(quad.samples$A, plot = FALSE)
A.acf.df <- with(A.acf, data.frame(lag, acf))

ggplot(data = A.acf.df, aes(lag, acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  theme(text = element_text(family = 'serif'))
ggsave('auto4.png')

phi.acf <- acf(quad.samples$phi, plot = FALSE)
phi.acf.df <- with(phi.acf, data.frame(lag, acf))

ggplot(data = phi.acf.df, aes(lag, acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  theme(text = element_text(family = 'serif'))
ggsave('auto5.png')

sigma.acf <- acf(quad.samples$sigma, plot = FALSE)
sigma.acf.df <- with(sigma.acf, data.frame(lag, acf))

ggplot(data = sigma.acf.df, aes(lag, acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  theme(text = element_text(family = 'serif'))
ggsave('auto6.png')

inverse_rescale <- function(x, original.max, original.min) {
  # Will transform posterior predictions back to original scale
  x * (original.max - original.min) + original.min
}

### Quadratic model intervals 

# Posterior predicted sample interval, training set
ppm.pred      <- apply(quad.samples$ppm_pred, 2, 
                       quantile, probs=c(0.025, 0.975))
# Posterior predicted sample interval, test set
ppm.test.pred <- apply(quad.samples$ppm_test_pred, 2,
                       quantile, probs=c(0.025, 0.975))
# Posterior predicted sample interval, future projection
ppm.future    <- apply(quad.samples$ppm_future, 2, 
                       quantile, probs=c(0.025, 0.975))

# Mean posterior predictions, test set
mean.test <- inverse_rescale(apply(quad.samples$ppm_test_pred, 2, mean),
                             max(df.orig$ppm),
                             min(df.orig$ppm))

# Extract upper and lower interval bounds; rescale to original scale
ppm.pred.lower <- inverse_rescale(ppm.pred[1,], 
                                  max(df.orig$ppm), 
                                  min(df.orig$ppm))
ppm.pred.upper <- inverse_rescale(ppm.pred[2,], 
                                  max(df.orig$ppm), 
                                  min(df.orig$ppm))

ppm.test.pred.lower <- inverse_rescale(ppm.test.pred[1,], 
                                       max(df.orig$ppm), 
                                       min(df.orig$ppm))
ppm.test.pred.upper <- inverse_rescale(ppm.test.pred[2,], 
                                       max(df.orig$ppm), 
                                       min(df.orig$ppm))

ppm.future.lower <- inverse_rescale(ppm.future[1,], 
                                    max(df.orig$ppm), 
                                    min(df.orig$ppm))
ppm.future.upper <- inverse_rescale(ppm.future[2,], 
                                    max(df.orig$ppm), 
                                    min(df.orig$ppm))

### Exponential model intervals

# Posterior predicted sample interval, training set
ppm.pred2      <- apply(exp.samples$ppm_pred, 2, 
                        quantile, probs=c(0.025, 0.975))
# Posterior predicted sample interval, test set
ppm.test.pred2 <- apply(exp.samples$ppm_test_pred, 2, 
                        quantile, probs=c(0.025, 0.975))
# Posterior predicted sample interval, future projection
ppm.future2    <- apply(exp.samples$ppm_future, 2, 
                        quantile, probs=c(0.025, 0.975))

# Mean posterior predictions, test set
mean.test2 <- inverse_rescale(apply(exp.samples$ppm_test_pred, 2, mean),
                              max(df.orig$ppm), 
                              min(df.orig$ppm))  

# Extract upper and lower interval bounds; rescale to original scale
ppm.pred.lower2 <- inverse_rescale(ppm.pred2[1,], 
                                   max(df.orig$ppm), 
                                   min(df.orig$ppm))
ppm.pred.upper2 <- inverse_rescale(ppm.pred2[2,], 
                                   max(df.orig$ppm), 
                                   min(df.orig$ppm))

ppm.test.pred.lower2 <- inverse_rescale(ppm.test.pred2[1,], 
                                        max(df.orig$ppm), 
                                        min(df.orig$ppm))
ppm.test.pred.upper2 <- inverse_rescale(ppm.test.pred2[2,], 
                                        max(df.orig$ppm), 
                                        min(df.orig$ppm))

ppm.future.lower2 <- inverse_rescale(ppm.future2[1,], 
                                     max(df.orig$ppm), 
                                     min(df.orig$ppm))
ppm.future.upper2 <- inverse_rescale(ppm.future2[2,], 
                                     max(df.orig$ppm), 
                                     min(df.orig$ppm))

### Train and test confidence interval plots

# Indices to select x-tick labels
idxs <- floor(seq(1, nrow(df.new), length.out = 15))

# Quadratic plot
ggplot(data = df.new, aes(day, ppm)) +
  coord_cartesian(xlim = c(1, 21791), ylim = c(300, 420)) +
  labs(title = "Posterior predictions with quadratic model",
       subtitle = "Confidence intervals for training and test date predictions",
       x = "Year", 
       y = expression(paste("CO"[2], " concentration (ppm)"))) +
  scale_x_continuous(labels = year(df.new$date)[idxs], 
                     breaks = df.new$day[idxs]) +
  geom_ribbon(data = train,
              aes(ymin = ppm.pred.lower,
                  ymax = ppm.pred.upper,
                  fill = '95% CI: Training set'),
              alpha = 0.5) +
  geom_ribbon(data = test,
              aes(ymin = ppm.test.pred.lower,
                  ymax = ppm.test.pred.upper,
                  fill = '95% CI: Test set'),
              alpha = 0.5) +
  geom_line(lwd = 0.3) +
  theme(text = element_text(family = 'serif'),
        legend.title = element_blank(),
        legend.position = c(0.8, 0.2))

ggsave('quad_intervals.png')

# Exponential plot
ggplot(data = df.new, aes(day, ppm)) +
  coord_cartesian(xlim = c(1, 21791), ylim = c(300, 420)) +
  labs(title = "Posterior predictions with exponential model",
       subtitle = "Confidence intervals for training and test date predictions",
       x = "Year", 
       y = expression(paste("CO"[2], " concentration (ppm)"))) +
  scale_x_continuous(labels = year(df.new$date)[idxs], 
                     breaks = df.new$day[idxs]) +
  geom_ribbon(data = train,
              aes(ymin = ppm.pred.lower2,
                  ymax = ppm.pred.upper2,
                  fill = '95% CI: Training set'),
              alpha = 0.5) +
  geom_ribbon(data = test,
              aes(ymin = ppm.test.pred.lower2,
                  ymax = ppm.test.pred.upper2,
                  fill = '95% CI: Test set'),
              alpha = 0.5) +
  geom_line(lwd = 0.3) +
  theme(text = element_text(family = 'serif'),
        legend.title = element_blank(),
        legend.position = c(0.8, 0.2))

ggsave('exp_intervals.png')

# Root mean square error
rmse <- function(estimate, obs){
  sqrt(mean((estimate - obs)^2))
}

rmse.quad <- rmse(mean.test, test$ppm)
rmse.exp <- rmse(mean.test2, test$ppm)

# Stan data (for both models)
stan.data <- list(
  N = nrow(df.orig),                    # number training data points
  t = df.orig$scaled_day,               # rescaled training dates/timesteps
  N_future = nrow(df.futr),             # number of future timesteps
  t_future = df.futr$scaled_day,        # rescaled future dates/timesteps
  ppm = df.orig$scaled_ppm,               # rescaled ppm
  k = (2 * pi * (tail(df.orig$day, 1) / # periodic frequency
                   365.25))
)

# THE QUADRATIC MODEL
quad.model = "
data {
int N;                    // number of timesteps
real t[N];                // timesteps
int N_future;             // number of future timesteps 
real t_future[N_future];  // future timesteps 
real ppm[N];              // observed ppm values
real k;                   // denominator of period
}

parameters {
real c0_prime;      // intercept
real c1_prime;      // linear component
real c2_prime;      // quadratic component
real A_prime;       // amplitude
real phi_prime;     // phase shift 
real sigma_prime;   // variance
}

transformed parameters {
real c0;
real c1;
real c2;
real A;
real phi;
real sigma;

c0 = c0_prime;      
c1 = c1_prime;   
c2 = c2_prime;                              
A = exp(A_prime);         
phi = phi_prime;                                        
sigma = exp(sigma_prime); 
}

model {
// Priors
c0_prime ~ normal(0.02, 0.01);   
c1_prime ~ normal(0.45, 0.1);
c2_prime ~ normal(0.5, 0.1);
A_prime ~ normal(-3.5, 0.5);     
phi_prime ~ normal(0, pi() / 2);           
sigma_prime ~ normal(-4.5, 1); 

// Likelihood 
for (n in 1:N){
ppm[n] ~ normal(c0 + 
c1 * t[n] + 
c2 * square(t[n]) + 
A  * cos((t[n] * k) + phi),
sigma); 
}
}

generated quantities {
// Posterior samples
real ppm_pred[N];            // posterior predicted, train set
real ppm_future[N_future];  // posterior predicted, proj. future

for (n in 1:N) {
ppm_pred[n] = normal_rng(c0 +                        // intercept
c1 * t[n] +                 // linear
c2 * square(t[n]) +         // quadratic
A  * cos((t[n] * k) + phi), // seasonal
sigma);                     // error
}

for (n in 1:N_future) {
ppm_future[n] = normal_rng(c0 + 
c1 * t_future[n] + 
c2 * square(t_future[n]) + 
A  * cos((t_future[n] * k) + phi),
sigma);
}
}
"

# THE EXP MODEL
exp.model = "
data {
int N;                   // number of timesteps
real t[N];               // timesteps
int N_future;            // number of future timesteps 
real t_future[N_future]; // future timesteps 
real ppm[N];             // observed ppm values
real k;                   // denominator of period
}

parameters {
real c0_prime;      // intercept
real c1_prime;      // scaling constant
real c2_prime;      // exponentiated constant
real A_prime;       // amplitude
real phi_prime;     // phase shift 
real sigma_prime;   // variance
}

transformed parameters {
real c0;
real c1;
real c2;
real A;
real phi;
real sigma;

c0 = tanh(c0_prime);      
c1 = 1 / (1 + exp(-c1_prime));   
c2 = exp(c2_prime); 
A = exp(A_prime);         
phi = phi_prime;                                        
sigma = exp(sigma_prime); 
}

model {
// Priors
c0_prime ~ normal(-0.7, 0.01);   
c1_prime ~ normal(0.4, 0.05);
c2_prime ~ normal(0.95, 0.05);
A_prime ~ normal(-3.45, 0.58);     
phi_prime ~ normal(0, pi() / 2);           
sigma_prime ~ normal(-4.5, 1); 

// Likelihood 
for (n in 1:N){
ppm[n] ~ normal(c0 + 
c1 * pow(c2, t[n]) + 
A  * cos((t[n] * k) + phi),
sigma); 
}
}

generated quantities {
// Posterior samples
real ppm_pred[N];            // posterior predicted, train set
real ppm_future[N_future];   // posterior predicted, proj. future

for (n in 1:N) {
ppm_pred[n] = normal_rng(c0 +                        // intercept
c1 * pow(c2, t[n]) + 
A  * cos((t[n] * k) + phi), // seasonal
sigma);                     // error
}


for (n in 1:N_future) {
ppm_future[n] = normal_rng(c0 + 
c1 * pow(c2, t_future[n]) + 
A  * cos((t_future[n] * k) + phi),
sigma);
}
}
"

quad.fit.full <- stan(
  model_code = quad.model, 
  data = stan.data,  
  chains = 3,           
  warmup = 1000,     
  iter = 2000,     
  cores = 2,            
  refresh = 1000,         
  control = list(adapt_delta = 0.999)
); skrrrahh(4)

exp.fit.full <- stan(
  model_code = exp.model, 
  data = stan.data,  
  chains = 3,           
  warmup = 1000,     
  iter = 2000,     
  cores = 2,            
  refresh = 1000,         
  control = list(adapt_delta = 0.999)
); skrrrahh(4)

quad.ml <- bridge_sampler(quad.fit.full, silent = TRUE)
exp.ml <-  bridge_sampler(exp.fit.full, silent = TRUE)

BF <- bayes_factor(quad.ml, exp.ml)

print("Quadratic model:")
print(quad.ml)

print("Exponential model:")
print(exp.ml)
print(BF)

quad.fit.full.samples <- extract(quad.fit.full)
exp.fit.full.samples <- extract(exp.fit.full)

# Posterior predicted sample interval, future projections
quad.future <- apply(quad.fit.full.samples$ppm_future, 2, 
                     quantile, probs=c(0.025, 0.975))

exp.future <- apply(exp.fit.full.samples$ppm_future, 2, 
                    quantile, probs=c(0.025, 0.975))

# Generate lower and upper bounds
quad.future.lower <- inverse_rescale(quad.future[1,], 
                                     max(df.orig$ppm), 
                                     min(df.orig$ppm))
quad.future.upper <- inverse_rescale(quad.future[2,], 
                                     max(df.orig$ppm), 
                                     min(df.orig$ppm))

exp.future.lower <- inverse_rescale(exp.future[1,], 
                                    max(df.orig$ppm), 
                                    min(df.orig$ppm))
exp.future.upper <- inverse_rescale(exp.future[2,], 
                                    max(df.orig$ppm), 
                                    min(df.orig$ppm))

# Indices to select x-tick labels
idxs <- floor(seq(1, nrow(df.new), length.out = 15))

# Plot quadratic forecast
ggplot(data = df.new, aes(day, ppm)) +
  geom_line(lwd = 0.3) + 
  coord_cartesian(xlim = c(1, 36442), ylim = c(300, 530)) +
  labs(title = "Posterior predictions with quadratic model: 40 year forecast",
       x = "Year", 
       y = expression(paste("CO"[2], " concentration (ppm)"))) +
  scale_x_continuous(labels = year(df.new$date)[idxs], 
                     breaks = df.new$day[idxs]) +
  geom_ribbon(data = df.futr,
              aes(ymin = quad.future.lower,
                  ymax = quad.future.upper,
                  fill = '95% CI: Forecast')) +
  geom_vline(xintercept = tail(df.orig$day, 1)) +
  annotate("text", x = 18500, y = 310, 
           label = "Original data", family = 'serif', col = 'gray50') +
  theme(text = element_text(family = 'serif'),
        legend.title = element_blank(),
        legend.position = c(0.7, 0.09))

ggsave('quad_forecast.png')

# Plot exponential forecast
ggplot(data = df.new, aes(day, ppm)) +
  geom_line(lwd = 0.3) + 
  coord_cartesian(xlim = c(1, 36442), ylim = c(300, 550)) +
  labs(title = "Posterior predictions with exponential model: 40 year forecast",
       x = "Year", 
       y = expression(paste("CO"[2], " concentration (ppm)"))) +
  scale_x_continuous(labels = year(df.new$date)[idxs], 
                     breaks = df.new$day[idxs]) +
  geom_ribbon(data = df.futr,
              aes(ymin = exp.future.lower,
                  ymax = exp.future.upper,
                  fill = '95% CI: Forecast')) +
  geom_vline(xintercept = tail(df.orig$day, 1)) +
  annotate("text", x = 18500, y = 310, 
           label = "Original data", family = 'serif', col = 'gray50') +
  theme(text = element_text(family = 'serif'),
        legend.title = element_blank(),
        legend.position = c(0.7, 0.09))

ggsave('exp_forecast.png')

# Confidence interval for CO2 on Jan. 5, 2058
print("Lower bound:")
print(tail(quad.future.lower, 1))
print("Upper bound:")
print(tail(quad.future.upper, 1))

# 450 ppm plot

# Plot
ggplot(data = df.new, aes(day, ppm)) +
  geom_line(lwd = 0.3) + 
  coord_cartesian(xlim = c(27560, 28200), ylim = c(440, 460)) +
  labs(title = "Zoom in on 95% confidence interval, future forecast",
       subtitle = "Lower and upper bound intersections with 450 ppm",
       x = "Date", 
       y = expression(paste("CO"[2], " concentration (ppm)"))) +
  scale_x_continuous(labels = ymd(df.futr$date)[seq(810, 910, 11)], 
                     breaks = df.futr$day[seq(810, 910, 11)]) +
  geom_ribbon(data = df.futr,
              aes(ymin = quad.future.lower,
                  ymax = quad.future.upper),
              alpha = 0.5,
              fill = 'indianred2') +
  geom_vline(aes(xintercept = df.futr$day[847], col = "Feb. 18, 2034")) +
  geom_vline(aes(xintercept = df.futr$day[903], col = "Mar. 17, 2035")) +
  geom_hline(yintercept = 450) +
  theme(text = element_text(family = 'serif'),
        legend.title = element_blank(),
        legend.position = c(0.1, 0.8),
        axis.text.x = element_text(angle = 40, vjust = 0.5))

ggsave('zoom-in.png')
