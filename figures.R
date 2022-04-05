## Script for evaluation of aggregation methods of NN methods

#### Housekeeping ####
rm(list=ls())
gc()

#### Settings ####
# Load package
library(lubridate)
library(scoringRules)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)

# Path for figures
pdf_path <- "plots/"

# Path of R-data
data_r_path <- "data/"

# Path of network ensemble data
data_ss_ens_path <- "-"
data_cs_ens_path <- "-"

# Path of aggregated network data
data_ss_agg_path <- "-"
data_cs_agg_path <- "-"

# Load functions
source(file = paste0(get_wd(), "/fn_basic.R"))
source(file = paste0(get_wd(), "/fn_eval.R"))

#### Initialize ####
# Vector for plotting on [0,1]
x_plot <- seq(0, 1, 0.001)

# Vector for plotting on [0,50]
x_plot50 <- seq(0, 50, 0.01)

# Evaluation measures
sr_eval <- c("crps", "me", "lgt", "cov")

# Skill scores
sr_skill <- c("crps")

# Network types
nn_vec <- c("drn", "bqn", "hen")

# Names of aggregation methods
agg_names <- c(
  "lp" = "Linear Pool",
  "vi" = "Vincentization",
  "vi-a" = "Vincentization (a)",
  "vi-w" = "Vincentization (w)",
  "vi-aw" = "Vincentization (a, w)"
)

# Names of aggregation methods
agg_abr <- c(
  "lp" = expression("LP"),
  "vi" = expression("V[0]^\"=\""),
  "vi-a" = expression("V[a]^\"=\""),
  "vi-w" = expression("V[0]^w"),
  "vi-aw" = expression("V[a]^w")
)
agg_abr <- c(
  "lp" = "LP",
  "vi" = "V[0]^\"=\"",
  "vi-a" = "V[a]^\"=\"",
  "vi-w" = "V[0]^w",
  "vi-aw" = "V[a]^w"
)

# Aggregation methods
agg_meths <- names(agg_names)

# Methods with coefficient estimation
coeff_meths <- c("vi-a", "vi-w", "vi-aw")

# Get colors
cols <- brewer.pal(n = 8,
                   name = "Dark2")

# Colors of aggregation methods
agg_col <- c(
  "lp" = cols[5],
  "vi" = cols[6],
  "vi-a" = cols[3],
  "vi-w" = cols[1],
  "vi-aw" = cols[4],
  "ens" = cols[8],
  "opt" = cols[2]
)

# Line types of aggregation methods
agg_lty <- c(
  "lp" = 2,
  "vi" = 1,
  "vi-a" = 1,
  "vi-w" = 1,
  "vi-aw" = 1,
  "ens" = 4,
  "opt" = 4
)

# Function to get legend
get_legend <- function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#### Section 2: Example aggregation ####
# Line width for plotting
lwd_plot <- 3

# Sizes
cex_main <- 1.3
cex_lab <- 1.2
cex_axis <- 1.2

# Evaluations for plotting
n_plot <- 1e+3

# Aggregation methods to plot
agg_meths_plot <- c("lp", "vi", "vi-a", "vi-w", "vi-aw")

# Number of forecasts to aggregate
n <- 2

# Lower and upper boundaries
l <- 0.5
u <- 13.5

# Parameters of individual distributions
mu_1 <- 7
mu_2 <- 10
sd_1 <- 1
sd_2 <- 1

# Data frame of distributions
df_plot <- data.frame(y = numeric(length = n_plot))

# Add probabilities and densities
for(temp in c("1", "2", agg_meths_plot)){
  df_plot[[paste0("p_", temp)]] <- df_plot[[paste0("d_", temp)]] <- 
    numeric(length = n_plot) }

# Vector to plot on
df_plot[["y"]] <- seq(from = l,
                      to = u,
                      length.out = n_plot)

# Calculate probabilities and densities of individual forecasts
for(temp in 1:2){
  df_plot[[paste0("p_", temp)]] <- pnorm(q = df_plot[["y"]],
                                         mean = get(paste0("mu_", temp)),
                                         sd = get(paste0("sd_", temp)))
  df_plot[[paste0("d_", temp)]] <- dnorm(x = df_plot[["y"]],
                                         mean = get(paste0("mu_", temp)),
                                         sd = get(paste0("sd_", temp)))
}

# LP
df_plot[["p_lp"]] <- rowMeans(df_plot[,paste0("p_", 1:2)])
df_plot[["d_lp"]] <- rowMeans(df_plot[,paste0("d_", 1:2)])

# For-Loop over Vincentization approaches
for(temp in agg_meths_plot[grepl("vi", agg_meths_plot, fixed = TRUE)]){
  # Intercept
  if(is.element(temp, c("vi", "vi-w"))){ a <- 0 }
  else if(is.element(temp, c("vi-a", "vi-aw"))){ a <- -6 }
  
  # Weights
  if(is.element(temp, c("vi", "vi-a"))){ w <- 1/n }
  else if(is.element(temp, c("vi-w", "vi-aw"))){ w <- 1/n + 0.15 }
  
  # Calculate mean and standard deviation
  mu_vi <- a + w*sum(c(mu_1, mu_2))
  sd_vi <- w*sum(c(sd_1, sd_2))
  
  # Calculate probabilities and densities
  df_plot[[paste0("p_", temp)]] <- pnorm(q = df_plot[["y"]],
                                         mean = mu_vi,
                                         sd = sd_vi)
  df_plot[[paste0("d_", temp)]] <- dnorm(x = df_plot[["y"]],
                                         mean = mu_vi,
                                         sd = sd_vi)
}

# Name of PDF
file_pdf <- paste0(pdf_path, "/aggregation_methods.pdf")

# Start PDF
pdf(file = file_pdf,
    height = 8,
    width = 25,
    pointsize = 28)

# Set margins
par(mfrow = c(1, 3),
    oma = c(4, 1, 0, 1),
    mar = c(2, 2, 3, 1))

## PDF
# Empty plot
plot(x = 0,
     y = 0,
     type = "n",
     xlab = "y",
     ylab = "f(y)",
     main = "Probability density function (PDF)",
     # main = "PDF",
     cex.axis = cex_axis,
     cex.lab = cex_lab,
     cex.main = cex_main,
     ylim = c(0, max(df_plot[,grepl("d_", colnames(df_plot), fixed = TRUE)])),
     xlim = c(l, u))

# Draw individual PDFs
for(i in 1:n){
  lines(x = df_plot[["y"]],
        y = df_plot[[paste0("d_", i)]],
        col = agg_col["ens"],
        lty = agg_lty["ens"],
        cex = lwd_plot,
        pch = lwd_plot,
        lwd = lwd_plot) 
}

# For-Loop over aggregation methods
for(temp in agg_meths_plot){
  lines(x = df_plot[["y"]],
        y = df_plot[[paste0("d_", temp)]],
        col = agg_col[temp],
        lty = agg_lty[temp],
        lwd = lwd_plot)
}

## CDF
# Empty plot
plot(x = 0,
     y = 0,
     type = "n",
     xlab = "y",
     ylab = "F(y)",
     main = "Cumulative distribution function (CDF)",
     # main = "CDF",
     cex.axis = cex_axis,
     cex.lab = cex_lab,
     cex.main = cex_main,
     ylim = c(0, 1),
     xlim = c(l, u))

# Draw individual CDFs
for(i in 1:n){
  lines(x = df_plot[["y"]],
        y = df_plot[[paste0("p_", i)]],
        col = agg_col["ens"],
        lty = agg_lty["ens"],
        lwd = lwd_plot) 
}

# For-Loop over aggregation methods
for(temp in agg_meths_plot){
  lines(x = df_plot[["y"]],
        y = df_plot[[paste0("p_", temp)]],
        col = agg_col[temp],
        lty = agg_lty[temp],
        lwd = lwd_plot)
}

## Quantile functions
# Empty plot
plot(x = 0,
     y = 0,
     type = "n",
     xlab = "p",
     ylab = "Q(p)",
     main = "Quantile function",
     cex.axis = cex_axis,
     cex.lab = cex_lab,
     cex.main = cex_main,
     xlim = c(0, 1),
     ylim = c(l, u))

# Draw individual quantile functions
for(i in 1:n){
  lines(y = df_plot[["y"]],
        x = df_plot[[paste0("p_", i)]],
        col = agg_col["ens"],
        lty = agg_lty["ens"],
        lwd = lwd_plot) 
}

# For-Loop over aggregation methods
for(temp in agg_meths_plot){
  lines(y = df_plot[["y"]],
        x = df_plot[[paste0("p_", temp)]],
        col = agg_col[temp],
        lty = agg_lty[temp],
        lwd = lwd_plot)
}

# Set margins
par(fig = c(0, 1, 0, 1),
    oma = c(0, 0, 0, 0),
    mar = c(0, 0, 0, 0),
    new = TRUE)

# Empty Plot
plot(x = 0,
     y = 0,
     type = 'l',
     bty = 'n',
     axes = FALSE)

# Add Legend
legend(x = "bottom",
       bty = "n",
       horiz = TRUE,
       inset = 0.01,
       xpd = TRUE,
       cex = 1.5,
       legend = as.expression(parse(text = c("F[1]", "F[2]", agg_abr))),
       col = c(rep(agg_col["ens"], 2), agg_col[agg_meths]),
       lty = c(rep(agg_lty["ens"], 2), agg_lty[agg_meths]),
       lwd = lwd_plot + 2)

# End PDF
dev.off()

#### Section 2: Aggregation of HEN ####
# Line width for plotting
lwd_plot <- 3

# Sizes
cex_main <- 1.7
cex_lab <- 1.7
cex_axis <- 1.7

# Number of forecasts to aggregate
n <- 2

# Number of bins
n_bins <- 8

# Parameters of underlying normal distribution
mu_1 <- n_bins/2 - 0.75
mu_2 <- n_bins/2 + 0.75
sd_1 <- 1
sd_2 <- 1

# Define binning
bin_edges <- seq(from = 0,
                 to = n_bins,
                 by = 1)

# Via normal distribution
p_bins <- rbind(dnorm(x = bin_edges[-1],
                      mean = mu_1,
                      sd = sd_1),
                dnorm(x = bin_edges[-1],
                      mean = mu_2,
                      sd = sd_2))

# Standardize probabilities
p_bins <- p_bins/rowSums(p_bins)

# Get accumulated sums
p_cum <- t(apply(p_bins, 1, cumsum))

# Get accumulated sums of aggregated forecast    
p_cum_sort <- unique(c(0, round(sort(as.vector(p_cum)), 4)))

# Generate corresponding bin probabilities
p_bins_f <- diff(p_cum_sort)

# Generate bin edges for each forecast
bin_edges_f <- rowMeans(sapply(1:n, function(i_sim){
  quant_hd(tau = p_cum_sort,
           probs = p_bins[i_sim,],
           bin_edges = bin_edges) }))

# Number of bins for aggregated forecast
n_bins_f <- length(bin_edges_f) - 1

# Name of PDF
file_pdf <- paste0(pdf_path, "/aggregation_hen.pdf")

# Start PDF
pdf(file = file_pdf,
    height = 9,
    width = 28,
    pointsize = 25)

# Plot together in two plots
par(mfrow = c(1, 3))

### PDF
# Empty plot
plot(x = 0,
     y = 0,
     xlim = range(bin_edges),
     ylim = c(0, max(c(p_bins, p_bins_f))),
     type = "n",
     main = "Probability density function (PDF)",
     # main = "PDF",
     cex.axis = cex_axis,
     cex.lab = cex_lab,
     cex.main = cex_main,
     xlab = "y",
     ylab = "f(y)")

# Axis

# For-Loop over individual forecasts
for(j in 1:n){ 
  # Standardize bin probabilities with length
  p_tilde <- p_bins[j,]/diff(bin_edges)
  
  # Color
  col_alpha <- col2rgb(agg_col["ens"])
  
  # For-Loop over bins
  for(i in 2:length(bin_edges)){
    # Rectangle for each bin
    rect(xleft = bin_edges[i-1],
         xright = bin_edges[i],
         ybottom = 0,
         ytop = p_tilde[i-1],
         col = rgb(col_alpha[1], col_alpha[2], col_alpha[3], 
                   max = 255, alpha = 50 + 50*(j-1)),
         lwd = 0,
         border = NA)
    
    # Draw borders
    lines(x = bin_edges[c(i-1, i-1, i)],
          y = c(0, p_tilde)[c(i-1, i, i)],
          lwd = lwd_plot,
          col = agg_col["ens"])
  }
  
  # close last box
  lines(x = bin_edges[rep((n_bins + 1), 2)],
        y = c(p_tilde[n_bins], 0),
        lwd = lwd_plot,
        col = agg_col["ens"])
}

# V_0^=
# Standardize bin probabilities with length
p_tilde <- p_bins_f/diff(bin_edges_f)

# Color
col_alpha <- col2rgb(agg_col["vi"])

# For-Loop over bins
for(i in 2:length(bin_edges_f)){
  # Rectangle for each bin
  rect(xleft = bin_edges_f[i-1],
       xright = bin_edges_f[i],
       ybottom = 0,
       ytop = p_tilde[i-1],
       col = rgb(col_alpha[1], col_alpha[2], col_alpha[3], 
                 max = 255, alpha = 100),
       lwd = 0,
       border = NA)
  
  # Draw borders
  lines(x = bin_edges_f[c(i-1, i-1, i)],
        y = c(0, p_tilde)[c(i-1, i, i)],
        lwd = lwd_plot,
        # lty = agg_lty["vi"],
        col = agg_col["vi"])
}

# close last box
lines(x = bin_edges_f[rep((n_bins_f + 1), 2)],
      y = c(p_tilde[n_bins_f], 0),
      lwd = lwd_plot,
      col = agg_col["vi"])

# LP
# Standardize bin probabilities with length
p_tilde <- colMeans(p_bins)/diff(bin_edges)

# Color
col_alpha <- col2rgb(agg_col["lp"])

# For-Loop over bins
for(i in 2:length(bin_edges)){
  # Rectangle for each bin
  rect(xleft = bin_edges[i-1],
       xright = bin_edges[i],
       ybottom = 0,
       ytop = p_tilde[i-1],
       col = rgb(col_alpha[1], col_alpha[2], col_alpha[3], 
                 max = 255, alpha = 100),
       lwd = 0,
       border = NA)
  
  # Draw borders
  lines(x = bin_edges[c(i-1, i-1, i)],
        y = c(0, p_tilde)[c(i-1, i, i)],
        lwd = lwd_plot,
        # lty = agg_lty["lp"],
        col = agg_col["lp"])
}

# close last box
lines(x = bin_edges[rep((n_bins + 1), 2)],
      y = c(p_tilde[n_bins], 0),
      lwd = lwd_plot,
      col = agg_col["lp"])

### CDF
# Empty plot
plot(x = 0,
     y = 0,
     type = "n",
     xlab = "y",
     ylab = "F(y)",
     main = "Cumulative distribution function (CDF)",
     # main = "CDF",
     cex.axis = cex_axis,
     cex.lab = cex_lab,
     cex.main = cex_main,
     ylim = c(0, 1),
     xlim = range(bin_edges))

# Draw vertical lines at bin edges
abline(v = bin_edges,
       lty = 2,
       col = "lightgrey")

# Draw individual CDFs
for(i in 1:n){
  lines(x = bin_edges,
        y = c(0, p_cum[i,]),
        col = agg_col["ens"],
        # lty = agg_lty["ens"],
        lwd = lwd_plot)
  points(x = bin_edges,
         y = c(0, p_cum[i,]),
         col = agg_col["ens"],
         pch = 4,
         lwd = lwd_plot)
}

# Draw vertical lines at bin edges for LP
for(i in 2:(n_bins + 1)){
  lines(x = rep(bin_edges[i], 2),
        y = p_cum[,(i-1)],
        lwd = lwd_plot - 1,
        lty = 2,
        col = agg_col["lp"]) }

# V_0^=
lines(x = bin_edges_f,
      y = p_cum_sort,
      col = agg_col["vi"],
      # lty = agg_lty["vi"],
      lwd = lwd_plot)
points(x = bin_edges_f,
       y = p_cum_sort,
       col = agg_col["vi"],
       pch = 4,
       lwd = lwd_plot)

# LP
lines(x = bin_edges,
      y = c(0, colMeans(p_cum)),
      col = agg_col["lp"],
      # lty = agg_lty["lp"],
      lwd = lwd_plot)
points(x = bin_edges,
       y = c(0, colMeans(p_cum)),
       col = agg_col["lp"],
       pch = 4,
       lwd = lwd_plot)

# Legend
legend(x = "bottomright",
       cex = 1.8,
       legend = as.expression(parse(text = c("F[1]", "F[2]", agg_abr[c("lp", "vi")]))),
       col = c(rep(agg_col["ens"], 2), agg_col[c("lp", "vi")]),
       # lty = c(rep(agg_lty["ens"], 2), agg_lty[c("lp", "vi")]),
       lwd = lwd_plot)

### Quantile functions
# Empty plot
plot(x = 0,
     y = 0,
     type = "n",
     xlab = "p",
     ylab = "Q(p)",
     main = "Quantile function",
     xlim = c(0, 1),
     cex.axis = cex_axis,
     cex.lab = cex_lab,
     cex.main = cex_main,
     ylim = range(bin_edges))

# Draw vertical lines at bin edges
abline(v = p_cum_sort,
       lty = 2,
       col = "lightgrey")

# Draw individual CDFs
for(i in 1:n){
  lines(y = bin_edges,
        x = c(0, p_cum[i,]),
        col = agg_col["ens"],
        # lty = agg_lty["ens"],
        lwd = lwd_plot)
  points(y = bin_edges,
         x = c(0, p_cum[i,]),
         col = agg_col["ens"],
         pch = 4,
         lwd = lwd_plot)
}

# Draw vertical lines at bin edges for V_0^=
for(i in 2:length(p_cum_sort)){
  lines(x = rep(p_cum_sort[i], 2),
        y = sapply(1:n, function(l){ 
          quant_hd(tau = p_cum_sort[i],
                   probs = p_bins[l,],
                   bin_edges = bin_edges) }),
        lty = 2,
        lwd = lwd_plot - 1,
        col = agg_col["vi"]) }

# V_0^=
lines(y = bin_edges_f,
      x = p_cum_sort,
      col = agg_col["vi"],
      # lty = agg_lty["lp"],
      lwd = lwd_plot)
points(y = bin_edges_f,
       x = p_cum_sort,
       col = agg_col["vi"],
       pch = 4,
       lwd = lwd_plot)

# LP
lines(y = bin_edges,
      x = c(0, colMeans(p_cum)),
      col = agg_col["lp"],
      # lty = agg_lty["lp"],
      lwd = lwd_plot)
points(y = bin_edges,
       x = c(0, colMeans(p_cum)),
       col = agg_col["lp"],
       pch = 4,
       lwd = lwd_plot)

# End PDF
dev.off()

#### Simu: Load data ####
# Load scores
load(file = paste0(data_r_path, "nn_agg_ss_scores.RData"))

#### Simu: Initialize ####
# Models considered
scenario_vec <- 1:4

# Number of simulations
n_sim <- 50

# Vector of ensemble members
n_ens_vec <- 2*(1:20)

#### Simu: Score panel ####
# Vector of scores/quantities to plot
score_vec <- c("crps", "crpss", "me", "lgt", "cov", "a", "w")

# For-Loop over scenarios
for(i_scenario in scenario_vec){
  #### Initialization ####
  # Only scenario
  df_sc <- subset(df_scores, (model == i_scenario))
  
  #### Calculate quantities ####
  # List for matrices
  df_plot <- data.frame(nn = character(),
                        metric = character(),
                        agg = character(),
                        n_ens = integer(),
                        score = numeric(),
                        stringsAsFactors = FALSE)
  
  # For-Loop over quantities of interest
  for(temp_sr in score_vec){
    # Consider special case CRPSS
    if(temp_sr == "crpss"){ temp_out <- "crps"
    }else{ temp_out <- temp_sr }
    
    # Get optimal score of scenario for CRPSS
    if(!is.element(temp_sr, c("a", "w"))){
      s_opt <- mean(subset(df_sc, (type == "ref"))[[temp_out]]) }
    
    # For-Loop over network variants
    for(temp_nn in nn_vec){
      # Only network type
      df_nn <- subset(df_sc, (nn == temp_nn))
      
      # For-Loop over ensemble sizes and aggregation methods
      for(i_ens in n_ens_vec){ for(temp_agg in c("opt", "ens", agg_meths)){
        # Skip Ensemble for skill
        if((temp_sr == "crpss") & (temp_agg == "ens")){ next 
        }else if((temp_sr == "crpss") & (temp_agg == "opt")){ next
        }else if((temp_sr == "a") & (temp_agg == "opt")){ next
        }else if((temp_sr == "w") & (temp_agg == "opt")){ next }
        
        # Get row of data frame
        i <- nrow(df_plot) + 1
        
        # Fill in data frame
        df_plot[i, "nn"] <- temp_nn
        df_plot[i, "metric"] <- temp_sr
        df_plot[i, "n_ens"] <- i_ens
        df_plot[i, "agg"] <- temp_agg
        
        # Reference: Average score of ensemble members
        s_ref <- mean(subset(df_nn, (n_rep <= i_ens) & (type == "ind"))[[temp_out]])
        
        # Special case: Average ensemble score
        if(temp_agg == "ens"){ df_plot[i, "score"] <- s_ref
        }else if(temp_agg == "opt"){ df_plot[i, "score"] <- s_opt
        }else{
          # Read out score
          df_plot[i, "score"] <- mean(subset(df_nn, (n_ens == i_ens) &
                                               (type == temp_agg))[[temp_out]])
          
          # Special case: CRPSS
          if(temp_sr == "crpss"){
            # Calculate skill
            df_plot[i, "score"] <- 100*(s_ref - df_plot[i, "score"])/(s_ref - s_opt)
          }
          
          # Relative weight difference to equal weights in %
          if(temp_sr == "w"){ df_plot[i, "score"] <- 100*(i_ens*df_plot[i, "score"] - 1) }
        }
      }}
    }
  }
  
  #### Prepare data ####
  # Remove negative CRPSS smaller than 5%
  df_plot <- subset(df_plot, (metric != "crpss") | (score >= -5))
  
  # Rename networks
  df_plot[["nn"]] <- c("drn" = "DRN",
                       "bqn" = "BQN",
                       "hen" = "HEN")[df_plot[["nn"]]]
  
  # Rename metrics
  score_labels <- c("a" = "Intercept",
                    "cov" = "PI coverage in %",
                    "crps" = "CRPS",
                    "crpss" = "CRPSS in %",
                    "lgt" = "PI length",
                    "me" = "Bias",
                    "w" = "Rel. weight diff. in %")
  
  # Order metrics
  df_plot[["metric"]] <- factor(df_plot[["metric"]],
                                levels = score_vec,
                                labels = score_labels[score_vec])
  
  # Intercept values
  hline_vec0 <- c(
    # "crpss" = 0,
    "me" = 0,
    "cov" = 100*19/21,
    "a" = 0,
    "w" = 0
  )
  
  # Minimal lower/upper limit: Per Hand for each scenario
  if(i_scenario == 1){ 
    # Upper
    u_vec <- c(
      "crps" = 1.15,
      "me" = 0.01,
      "a" = 0.01,
      "w" = 1
    )
    
    # Lower
    l_vec <- c(
      "crps" = 0.97,
      "a" = -0.001
    )
  }else if(i_scenario == 2){ 
    # Upper
    u_vec <- c(
      "crps" = 2.32,
      "a" = 0.07,
      "w" = 1.1
    )
    
    # Lower
    l_vec <- c(
      "crps" = 2.19,
      "a" = -0.07,
      w = -0.1
    )
  }else if(i_scenario == 3){ 
    # Upper
    u_vec <- c(
      "crps" = 0.7,
      "a" = 0.01,
      "w" = 2
    )
    
    # Lower
    l_vec <- c(
      "crps" = 0.65,
      "a" = -0.01
    )
  }else if(i_scenario == 4){ 
    # Upper
    u_vec <- c(
      "crps" = 0.64,
      "a" = 0.01,
      "w" = 0.5
    )
    
    # Lower
    l_vec <- c(
      "crps" = 0.42,
      "a" = -0.07
    )
  }else{ u_vec <- l_vec <- c() }
  
  # For-Loop over upper and lower
  for(temp_limit in c("l", "u")){
    # Upper or lower vector
    temp_vec <- get(paste0(temp_limit, "_vec"))
    
    # Check if limits are given
    if(length(temp_vec) == 0){ next }
    
    # Get upper values
    temp_df <- data.frame(metric = character(length = length(temp_vec)),
                         value = numeric(length = length(temp_vec)),
                         stringsAsFactors = FALSE)
    
    # Set intercepts
    temp_df[,1] <- names(temp_vec)
    temp_df[,2] <- temp_vec
    
    # Order metrics
    temp_df[["metric"]] <- factor(temp_df[["metric"]],
                                 levels = score_vec,
                                 labels = score_labels[score_vec])
    
    # Assign
    assign(x = paste0(temp_limit, "_plot"),
           value = temp_df)
    rm(temp_df, temp_vec)
  }
  
  # Legend labels
  leg_labels <- c("ens" = "DE", 
                  "opt" = "F^\"*\"", 
                  agg_abr)
  
  #### PDF ####
  # For-Loop over networks
  for(temp_nn in c("DRN", "BQN", "HEN")){
    ## CRPSS hline depending on method
    # Draw no skill only if negative skill is included
    if(any(subset(df_plot, (nn == temp_nn) &
                  (metric == score_labels[["crpss"]]))[["score"]] <= 1)){ 
      hline_vec <- c(hline_vec0, "crpss" = 0) }
    else{ hline_vec <- hline_vec0 }
    
    # Get intercept values
    hline_plot <- data.frame(metric = character(length = length(hline_vec)),
                             value = numeric(length = length(hline_vec)),
                             stringsAsFactors = FALSE)
    
    # Set intercepts
    hline_plot[,1] <- names(hline_vec)
    hline_plot[,2] <- hline_vec
    
    # Order metrics
    hline_plot[["metric"]] <- factor(hline_plot[["metric"]],
                                     levels = score_vec,
                                     labels = score_labels[score_vec])
    
    # Make plot
    temp_plot <- ggplot(subset(df_plot, nn == temp_nn), 
                          aes(x = n_ens, y = score, color = agg)) + 
      geom_line(size = 1,
                aes(linetype = agg)) +
      scale_color_manual(values = agg_col[names(leg_labels)],
                         labels = as.expression(parse(text = leg_labels))) +
      scale_linetype_manual(values = agg_lty[names(leg_labels)],
                         labels = as.expression(parse(text = leg_labels))) +
      facet_grid(rows = vars(metric), 
                 cols = vars(nn), 
                 scales = "free") +
      geom_hline(data = hline_plot, aes(yintercept = value), 
                 linetype = "dashed") +
      theme_bw() +
      theme(axis.text.x = element_text(size = 10))
    
    # Add upper and lower
    if(length(l_vec) > 0){ temp_plot <- temp_plot + 
      geom_hline(data = l_plot, aes(yintercept = value), 
                 linetype = "blank") }
    if(length(u_vec) > 0){ temp_plot <- temp_plot + 
      geom_hline(data = u_plot, aes(yintercept = value), 
                 linetype = "blank") }
    
    # Legend only once
    if(temp_nn == "DRN"){
      # make legend
      temp_plot <- temp_plot + 
        theme(legend.position = "top", 
              legend.justification = "center",
              legend.text = element_text(size = rel(1.1)), 
              legend.title = element_blank()) +
        guides(colour = guide_legend(nrow = 1))
      
      # Get legend
      my_legend <- get_legend(temp_plot) 
    }
    
    # No legend in plot
    temp_plot <- temp_plot + theme(legend.position = "none")
    
    # Text on y-axis
    if(temp_nn != "HEN"){ temp_plot <- temp_plot + theme(strip.text.y = element_blank()) }
    
    # Different y-labels
    if(temp_nn == "DRN"){ temp_plot <- temp_plot + ylab("Value") 
    }else{ temp_plot <- temp_plot + ylab(element_blank()) }
    
    # Different x-labels
    if(temp_nn == "BQN"){ temp_plot <- temp_plot + xlab("Network ensemble size") 
    }else{ temp_plot <- temp_plot + xlab("") }
    
    # Assign name
    assign(x = paste0("plot_", temp_nn),
           value = temp_plot)
    rm(temp_plot)
  }
  
  # Plot together
  pdf_plot <- grid.arrange(my_legend, 
                           plot_DRN, plot_BQN, plot_HEN, 
                           ncol = 3, nrow = 2,
                           layout_matrix = rbind(c(1, 1, 1),
                                                 c(2, 3, 4)),
                           heights = c(0.2, 3))
  
  # Save as PDF
  ggsave(filename = paste0("ss_gg_panel_model", i_scenario, ".pdf"),
         path = pdf_path,
         plot = pdf_plot,
         width = 10,
         height = 1.9*length(score_vec),
         scale = 0.8)
}

#### Simu: CRPSS boxplot panel ####
# For-Loop over scenarios
for(i_scenario in scenario_vec){
  #### Initialization ####
  # Only scenario
  df_sc <- subset(df_scores, model == i_scenario)
  
  #### Calculate quantities ####
  # List for matrices
  df_plot <- data.frame(nn = character(),
                        agg = character(),
                        n_ens = integer(),
                        crpss = numeric(),
                        stringsAsFactors = FALSE)
  
  # For-Loop over network variants, ensemble sizes and aggregation methods
  for(temp_nn in nn_vec){ for(i_ens in n_ens_vec){ for(temp_agg in agg_meths){
    # Get subset of scores data
    df_sub <- subset(df_sc, (n_ens == i_ens) & 
                       (type == temp_agg) & (nn == temp_nn))
    
    # Get row of data frame
    i <- nrow(df_plot) + 1:nrow(df_sub)
    
    # Fill in data frame
    df_plot[i, "n_ens"] <- i_ens
    df_plot[i, "nn"] <- temp_nn
    df_plot[i, "agg"] <- temp_agg
    
    # Reference: Average score of ensemble members
    df_plot[i, "crpss"] <- 100*df_sub[["crpss"]]
  }}}
  
  #### Prepare data ####
  # Minimal lower/upper limit: Per Hand for each scenario
  if(i_scenario == 1){ lu_vec <- c(-20, 40)
  }else{ lu_vec <- range(df_plot[["crpss"]]) }
  
  # Rename networks
  df_plot[["nn"]] <- c("drn" = "DRN",
                       "bqn" = "BQN",
                       "hen" = "HEN")[df_plot[["nn"]]]
  
  # Change order of networks
  df_plot[["nn"]] <- factor(df_plot[["nn"]], 
                            levels = c("DRN", "BQN", "HEN"))
  
  # Rename aggregation methods
  df_plot[["agg"]] <- agg_abr[df_plot[["agg"]]]
  
  # Colors for aggregation methods
  agg_col_gg <- agg_col
  names(agg_col_gg) <- agg_abr[names(agg_col)]
  
  #### PDF ####
  # Make plot
  pdf_plot <- 
    ggplot(df_plot, 
           aes(x = n_ens, y = crpss, group = n_ens, 
               col = agg, fill = agg)) + 
    facet_grid(rows = vars(nn),
               cols = vars(agg),
               labeller = label_parsed,
               scales = "free") +
    geom_boxplot(alpha = 0.5, 
                 outlier.size = 0.6,
                 na.rm = TRUE) +
    scale_color_manual(values = agg_col_gg[unique(df_plot[["agg"]])]) +
    scale_fill_manual(values = agg_col_gg[unique(df_plot[["agg"]])]) +
    coord_cartesian(ylim = lu_vec) +
    xlab("Network ensemble size") + 
    ylab("CRPSS in %")  +
    theme_bw() +
    theme(axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15)) + 
    theme(legend.position = "none") +
    guides(colour = guide_legend(nrow = 1))
  
  # Add intercept
  if(min(lu_vec) < 0.1){ pdf_plot <- pdf_plot + geom_hline(yintercept = 0, 
                                                         linetype = "dashed",
                                                         color = agg_col["ens"]) }
  
  # Save as PDF
  ggsave(filename = paste0("ss_gg_panel_crpss_boxplots_model", i_scenario, ".pdf"),
         path = pdf_path,
         plot = pdf_plot,
         width = 20,
         height = 15,
         scale = 0.7)
}

#### Simu: PIT histograms ####
# Network ensemble size
n_ens <- 2

# Number of bins in histogram
n_bins <- 21

# For-Loop over scenarios
for(i_scenario in scenario_vec){
  #### Initialization ####
  # Only scenario
  df_sc <- subset(df_scores, (model == i_scenario))
  
  # Index vector for validation and testing
  i_valid <- 1:unique(subset(df_sc, type != "ref")[["n_valid"]])
  i_test <- max(i_valid) + 1:unique(subset(df_sc, type != "ref")[["n_test"]])
  
  #### Get PIT values ####
  # List for PIT values
  pit_ls <- list()
  
  # For-Loop over network variants
  for(temp_nn in nn_vec){
    #### PIT values of ensemble member ####
    # Vector for PIT values
    temp_pit <- c()

    # For-Loop over ensemble member and simulations
    for(i_rep in 1:n_ens){ for(i_sim in 1:n_sim){
      # Load ensemble member
      load(file = paste0(data_ss_ens_path, "model", i_scenario,
                         "/model", i_scenario, "_", temp_nn,
                         "_sim", i_sim, "_ens", i_rep, ".RData"))

      # Read out
      temp_pit <- c(temp_pit, pred_nn[["scores"]][["pit"]][i_test])
    }}

    # Save PIT
    pit_ls[[paste0(temp_nn, "_ens")]] <- temp_pit

    #### Aggregation methods ####
    # For-Loop over aggregation methods
    for(temp_agg in agg_meths){
      # Vector for PIT values
      temp_pit <- c()

      # For-Loop over simulations
      for(i_sim in 1:n_sim){
        # Load aggregated forecasts
        load(file = paste0(data_ss_agg_path, "model", i_scenario,
                           "/model", i_scenario, "_", temp_nn, "_sim", i_sim,
                           "_", temp_agg, "_ens", n_ens, ".RData"))

        # Read out PIT-values
        temp_pit <- c(temp_pit, pred_agg[["scores"]][["pit"]])
      }

      # Save PIT
      pit_ls[[paste0(temp_nn, "_", temp_agg)]] <- temp_pit
    }
  }
  
  #### Calculate histograms ####
  # List for matrices
  df_plot <- data.frame(nn = character(),
                        agg = character(),
                        breaks = numeric(),
                        pit = numeric(),
                        stringsAsFactors = FALSE)
  
  # For-Loop over network variants and aggregation methods
  for(temp_nn in nn_vec){ for(temp_agg in c("ens", agg_meths)){
    # Calculate histogram and read out values (see pit function)
    temp <- hist(pit_ls[[paste0(temp_nn, "_", temp_agg)]],
                 breaks = (0:n_bins)/n_bins,
                 plot = FALSE)
    
    # Get row of data frame
    i <- nrow(df_plot) + 1:(n_bins + 1)
      
    # Fill in data frame
    df_plot[i, "nn"] <- temp_nn
    df_plot[i, "agg"] <- temp_agg
    
    # Breaks
    df_plot[i, "breaks"] <- temp[["breaks"]]
    
    # Density
    df_plot[i, "pit"] <- c(0, temp[["density"]])
  }}
  
  #### Prepare data ####
  # Rename networks
  df_plot[["nn"]] <- c("drn" = "DRN",
                       "bqn" = "BQN",
                       "hen" = "HEN")[df_plot[["nn"]]]
  
  # Change order of networks
  df_plot[["nn"]] <- factor(df_plot[["nn"]], 
                            levels = c("DRN", "BQN", "HEN"))
  
  # Order metrics
  df_plot[["agg"]] <- factor(df_plot[["agg"]],
                             levels = c("ens", agg_meths),
                             labels = c("DE",
                                        agg_abr[agg_meths]))
  
  #### PDF ####
  # Limit of y-axis
  if(i_scenario == 1){ y_lim <- 1.5
  }else if(i_scenario == 2){ y_lim <- 2.5
  }else if(i_scenario == 3){ y_lim <- 1.5
  }else if(i_scenario == 4){ y_lim <- 2 }

  # Make plot
  pdf_plot <- ggplot(df_plot, 
                     aes(x = breaks, y = pit)) + 
    geom_rect(aes(xmin = breaks, xmax = lead(breaks),
                  ymin = 0, ymax = lead(pit)),
                  color = "black", fill = "black", alpha = 0.3) +
    facet_grid(cols = vars(agg), 
               rows = vars(nn), 
               labeller = label_parsed,
               scales = "free")  +
    ylab("Density") +
    xlab("PIT") +
    theme_bw() +
    theme(axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          strip.text.x = element_text(size = 12),
          strip.text.y = element_text(size = 12)) + 
    geom_hline(aes(yintercept = 1), 
               linetype = "dashed") +
    geom_hline(aes(yintercept = y_lim), 
               linetype = "blank")

  # Save as PDF
  ggsave(filename = paste0("ss_gg_pit_ens", n_ens, "_model", i_scenario, ".pdf"),
         path = pdf_path,
         plot = pdf_plot,
         width = 15,
         height = 5,
         scale = 0.8)
}



#### Case: Load data ####
# From general parameters
load(file = paste0(data_r_path, "general_pars.RData"))

# Load scores
load(file = paste0(data_r_path, "nn_agg_cs_scores_rep.RData"))

# Treat Quasi-Vincentization as V^w_0
df_scores <- subset(df_scores, !((type == "vi") & (nn == "drn")))

# Rename "vi-q"
df_scores[["type"]][df_scores[["type"]] == "vi-q"] <- "vi"

#### Case: Initialize ####
# Initialization hour
init_hr <- 0

# Vector of forecast steps
temp_step_vec <- c(0, 6, 12, 18)

# Number of simulations carried out (network runs)
n_sim <- 100

# Number of aggregation repetitions
n_rep <- 20

# Vector of ensemble members
n_ens_vec <- 2*(1:20)

#### Case: Score panel ####
# Vector of scores/quantities to plot
score_vec <- c("crps", "crpss", "me", "lgt", "cov", "a", "w")

# List for matrices
df_plot <- data.frame(nn = character(),
                      metric = character(),
                      agg = character(),
                      n_ens = integer(),
                      score = numeric(),
                      stringsAsFactors = FALSE)

# For-Loop over quantities of interest
for(temp_sr in score_vec){
  # Consider special case CRPSS
  if(temp_sr == "crpss"){ temp_out <- "crps"
  }else{ temp_out <- temp_sr }
  
  # For-Loop over network variants
  for(temp_nn in nn_vec){
    # Only network type
    df_nn <- subset(df_scores, (nn == temp_nn))
    
    # For-Loop over ensemble sizes and aggregation methods
    for(i_ens in n_ens_vec){ for(temp_agg in c("ens", agg_meths)){
      # Skip Ensemble for skill and no estimation
      if((temp_sr == "crpss") & (temp_agg == "ens")){ next 
      }else if((temp_sr == "a") & !is.element(temp_agg, c("vi-a", "vi-aw"))){ next
      }else if((temp_sr == "w") & !is.element(temp_agg, c("vi-w", "vi-aw"))){ next }
      
      # Get row of data frame
      i <- nrow(df_plot) + 1
      
      # Fill in data frame
      df_plot[i, "nn"] <- temp_nn
      df_plot[i, "metric"] <- temp_sr
      df_plot[i, "n_ens"] <- i_ens
      df_plot[i, "agg"] <- temp_agg
      
      # Subset of repetitions
      if(temp_agg == "ens"){ df_ens <- subset(df_nn, (n_ens == i_ens))
      }else{ df_ens <- subset(df_nn, (n_ens == i_ens) & (type == temp_agg)) }
      
      # Reference: Average score of ensemble members
      if(!is.element(temp_sr, c("a", "w"))){ 
        s_ref <- mean(df_ens[[paste0(temp_out, "_ref")]]) }
      
      # Special case: Average ensemble score
      if(temp_agg == "ens"){ df_plot[i, "score"] <- s_ref
      }else{
        # Read out score
        df_plot[i, "score"] <- mean(df_ens[[temp_out]])
        
        # Special case: CRPSS
        if(temp_sr == "crpss"){ df_plot[i, "score"] <- 100*(1 - df_plot[i, "score"]/s_ref) }
        
        # Relative weight difference to equal weights in %
        if(temp_sr == "w"){ df_plot[i, "score"] <- 100*(i_ens*df_plot[i, "score"] - 1) }
      }
    }}
  }
}

# Remove negative CRPSS smaller than 5%
df_plot <- subset(df_plot, (metric != "crpss") | (score >= -5))

# Rename networks
df_plot[["nn"]] <- c("drn" = "DRN",
                     "bqn" = "BQN",
                     "hen" = "HEN")[df_plot[["nn"]]]

# Rename metrics
score_labels <- c("a" = "Intercept in m/s",
                  "cov" = "PI coverage in %",
                  "crps" = "CRPS in m/s",
                  "crpss" = "CRPSS in %",
                  "lgt" = "PI length in m/s",
                  "me" = "Bias in m/s",
                  "w" = "Rel. weight diff. in %")

# Order metrics
df_plot[["metric"]] <- factor(df_plot[["metric"]],
                              levels = score_vec,
                              labels = score_labels[score_vec])

# Intercept values
hline_vec <- c(
  "me" = 0,
  "cov" = 100*19/21,
  "a" = 0,
  "w" = 0
)

# Draw no skill only if negative skill is included
if(any(subset(df_plot, metric == score_labels[["crpss"]])[["score"]] <= 0)){ 
  hline_vec <- c(hline_vec, "crpss" = 0) }

# Get intercept values
hline_plot <- data.frame(metric = character(length = length(hline_vec)),
                         value = numeric(length = length(hline_vec)),
                         stringsAsFactors = FALSE)

# Set intercepts
hline_plot[,1] <- names(hline_vec)
hline_plot[,2] <- hline_vec

# Order metrics
hline_plot[["metric"]] <- factor(hline_plot[["metric"]],
                                 levels = score_vec,
                                 labels = score_labels[score_vec])

# Minimal lower/upper limit: Per Hand
u_vec <- c(
  "crps" = 0.865,
  "me" = 0.05,
  "a" = 0,
  "w" = 4
)

# Lower
l_vec <- c(
  "crps" = 0.83,
  "me" = -0.1,
  "a" = -0.3,
  "w" = -1
)

# For-Loop over upper and lower
for(temp_limit in c("l", "u")){
  # Upper or lower vector
  temp_vec <- get(paste0(temp_limit, "_vec"))
  
  # Check if limits are given
  if(length(temp_vec) == 0){ next }
  
  # Get upper values
  temp_df <- data.frame(metric = character(length = length(temp_vec)),
                        value = numeric(length = length(temp_vec)),
                        stringsAsFactors = FALSE)
  
  # Set intercepts
  temp_df[,1] <- names(temp_vec)
  temp_df[,2] <- temp_vec
  
  # Order metrics
  temp_df[["metric"]] <- factor(temp_df[["metric"]],
                                levels = score_vec,
                                labels = score_labels[score_vec])
  
  # Assign
  assign(x = paste0(temp_limit, "_plot"),
         value = temp_df)
  rm(temp_df, temp_vec)
}

# For-Loop over networks
for(temp_nn in c("DRN", "BQN", "HEN")){
  # Make plot
  temp_plot <- ggplot(subset(df_plot, nn == temp_nn), 
                      aes(x = n_ens, y = score, color = agg)) + 
    geom_line(size = 1,
              aes(linetype = agg)) +
    scale_color_manual(values = agg_col[c("ens", agg_meths)],
                       labels = as.expression(parse(text = c("ens" = "DE", agg_abr)))) +
    scale_linetype_manual(values = agg_lty[c("ens", agg_meths)],
                          labels = as.expression(parse(text = c("ens" = "DE", agg_abr)))) +
    facet_grid(rows = vars(metric), 
               cols = vars(nn), 
               scales = "free") +
    geom_hline(data = hline_plot, aes(yintercept = value), 
               linetype = "dashed") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 10))
  
  # Add upper and lower
  if(length(l_vec) > 0){ temp_plot <- temp_plot + 
    geom_hline(data = l_plot, aes(yintercept = value), 
               linetype = "blank") }
  if(length(u_vec) > 0){ temp_plot <- temp_plot + 
    geom_hline(data = u_plot, aes(yintercept = value), 
               linetype = "blank") }
  
  # Legend only once
  if(temp_nn == "DRN"){
    # make legend
    temp_plot <- temp_plot + 
      theme(legend.position = "top", 
            legend.justification = "center", 
            legend.text = element_text(size = rel(1.1)), 
            legend.title = element_blank()) +
      guides(colour = guide_legend(nrow = 1))
    
    # Get legend
    my_legend <- get_legend(temp_plot) 
  }
  
  # No legend in plot
  temp_plot <- temp_plot + theme(legend.position = "none")
  
  # Text on y-axis
  if(temp_nn != "HEN"){ temp_plot <- temp_plot + theme(strip.text.y = element_blank()) }
  
  # Different y-labels
  if(temp_nn == "DRN"){ temp_plot <- temp_plot + ylab("Value") 
  }else{ temp_plot <- temp_plot + ylab(element_blank()) }
  
  # Different x-labels
  if(temp_nn == "BQN"){ temp_plot <- temp_plot + xlab("Network ensemble size") 
  }else{ temp_plot <- temp_plot + xlab("") }
  
  # Assign name
  assign(x = paste0("plot_", temp_nn),
         value = temp_plot)
  rm(temp_plot)
}

# Plot together
pdf_plot <- grid.arrange(my_legend, 
                         plot_DRN, plot_BQN, plot_HEN, 
                         ncol = 3, nrow = 2,
                         layout_matrix = rbind(c(1, 1, 1),
                                               c(2, 3, 4)),
                         heights = c(0.2, 3))

# Save as PDF
ggsave(filename = paste0("cs_gg_panel.pdf"),
       path = pdf_path,
       plot = pdf_plot,
       width = 10,
       height = 1.9*length(score_vec),
       scale = 0.8)

#### Case: Stepwise score panels ####
# Vector of scores/quantities to plot
score_vec <- c("crps", "crpss", "me", "lgt", "cov", "a", "w")

# For-Loop over forecast steps
for(temp_step in temp_step_vec){
  #### Initialization ####
  # Only forecast step
  df_step <- subset(df_scores, (fc_step == temp_step))
  
  #### Calculate quantities ####
  # List for matrices
  df_plot <- data.frame(nn = character(),
                        metric = character(),
                        agg = character(),
                        n_ens = integer(),
                        score = numeric(),
                        stringsAsFactors = FALSE)
  
  # For-Loop over quantities of interest
  for(temp_sr in score_vec){
    # Consider special case CRPSS
    if(temp_sr == "crpss"){ temp_out <- "crps"
    }else{ temp_out <- temp_sr }
    
    # For-Loop over network variants
    for(temp_nn in nn_vec){
      # Only network type
      df_nn <- subset(df_step, (nn == temp_nn))
      
      # For-Loop over ensemble sizes and aggregation methods
      for(i_ens in n_ens_vec){ for(temp_agg in c("ens", agg_meths)){
        # Skip Ensemble for skill and no estimation
        if((temp_sr == "crpss") & (temp_agg == "ens")){ next 
        }else if((temp_sr == "a") & !is.element(temp_agg, c("vi-a", "vi-aw"))){ next
        }else if((temp_sr == "w") & !is.element(temp_agg, c("vi-w", "vi-aw"))){ next }
        
        # Get row of data frame
        i <- nrow(df_plot) + 1
        
        # Fill in data frame
        df_plot[i, "nn"] <- temp_nn
        df_plot[i, "metric"] <- temp_sr
        df_plot[i, "n_ens"] <- i_ens
        df_plot[i, "agg"] <- temp_agg
        
        # Subset of repetitions
        if(temp_agg == "ens"){ df_ens <- subset(df_nn, (n_ens == i_ens))
        }else{ df_ens <- subset(df_nn, (n_ens == i_ens) & (type == temp_agg)) }
        
        # Reference: Average score of ensemble members
        if(!is.element(temp_sr, c("a", "w"))){ 
          s_ref <- mean(df_ens[[paste0(temp_out, "_ref")]]) }
        
        # Special case: Average ensemble score
        if(temp_agg == "ens"){ df_plot[i, "score"] <- s_ref
        }else{
          # Read out score
          df_plot[i, "score"] <- mean(df_ens[[temp_out]])
          
          # Special case: CRPSS
          if(temp_sr == "crpss"){ df_plot[i, "score"] <- 100*(1 - df_plot[i, "score"]/s_ref) }
          
          # Relative weight difference to equal weights in %
          if(temp_sr == "w"){ df_plot[i, "score"] <- 100*(i_ens*df_plot[i, "score"] - 1) }
        }
      }}
    }
  }
  
  #### Prepare data ####
  # Remove negative CRPSS smaller than 5%
  df_plot <- subset(df_plot, (metric != "crpss") | (score >= -5))
  
  # Rename networks
  df_plot[["nn"]] <- c("drn" = "DRN",
                       "bqn" = "BQN",
                       "hen" = "HEN")[df_plot[["nn"]]]
  
  # Rename metrics
  score_labels <- c("a" = "Intercept in m/s",
                    "cov" = "PI coverage in %",
                    "crps" = "CRPS in m/s",
                    "crpss" = "CRPSS in %",
                    "lgt" = "PI length in m/s",
                    "me" = "Bias in m/s",
                    "w" = "Rel. weight diff. in %")
  
  # Order metrics
  df_plot[["metric"]] <- factor(df_plot[["metric"]],
                                levels = score_vec,
                                labels = score_labels[score_vec])
  
  # Intercept values
  hline_vec <- c(
    # "crpss" = 0,
    "me" = 0,
    "cov" = 100*19/21,
    "a" = 0,
    "w" = 0
  )
  
  # Draw no skill only if negative skill is included
  if(any(subset(df_plot, metric == score_labels[["crpss"]])[["score"]] <= 0)){ 
    hline_vec <- c(hline_vec, "crpss" = 0) }
  
  # Get intercept values
  hline_plot <- data.frame(metric = character(length = length(hline_vec)),
                           value = numeric(length = length(hline_vec)),
                           stringsAsFactors = FALSE)
  
  # Set intercepts
  hline_plot[,1] <- names(hline_vec)
  hline_plot[,2] <- hline_vec
  
  # Order metrics
  hline_plot[["metric"]] <- factor(hline_plot[["metric"]],
                                   levels = score_vec,
                                   labels = score_labels[score_vec])
  
  # Minimal lower/upper limit: Per Hand for each step
  # if(temp_step >= 0){ 
  # Upper
  u_vec <- c(
    "me" = 0.08,
    "a" = 0,
    "w" = 6
  )
  
  # Lower
  l_vec <- c(
    "me" = -0.2,
    "a" = -0.5,
    "w" = -1.5
  )
  
  # CRPS
  if(temp_step == 0){ 
    u_vec <- c(u_vec, "crps" = 0.785)
    l_vec <- c(l_vec, "crps" = 0.745)
  }else if(temp_step == 6){ 
    u_vec <- c(u_vec, "crps" = 0.82)
    l_vec <- c(l_vec, "crps" = 0.79)
  }else if(temp_step == 12){ 
    u_vec <- c(u_vec, "crps" = 0.885)
    l_vec <- c(l_vec, "crps" = 0.85)
  }else if(temp_step == 18){ 
    u_vec <- c(u_vec, "crps" = 0.98)
    l_vec <- c(l_vec, "crps" = 0.945)
  }
  
  # For-Loop over upper and lower
  for(temp_limit in c("l", "u")){
    # Upper or lower vector
    temp_vec <- get(paste0(temp_limit, "_vec"))
    
    # Check if limits are given
    if(length(temp_vec) == 0){ next }
    
    # Get upper values
    temp_df <- data.frame(metric = character(length = length(temp_vec)),
                          value = numeric(length = length(temp_vec)),
                          stringsAsFactors = FALSE)
    
    # Set intercepts
    temp_df[,1] <- names(temp_vec)
    temp_df[,2] <- temp_vec
    
    # Order metrics
    temp_df[["metric"]] <- factor(temp_df[["metric"]],
                                  levels = score_vec,
                                  labels = score_labels[score_vec])
    
    # Assign
    assign(x = paste0(temp_limit, "_plot"),
           value = temp_df)
    rm(temp_df, temp_vec)
  }
  
  #### PDF ####
  # For-Loop over networks
  for(temp_nn in c("DRN", "BQN", "HEN")){
    # Make plot
    temp_plot <- ggplot(subset(df_plot, nn == temp_nn), 
                        aes(x = n_ens, y = score, color = agg)) + 
      geom_line(size = 1,
                aes(linetype = agg)) +
      scale_color_manual(values = agg_col[c("ens", agg_meths)],
                         labels = as.expression(parse(text = c("ens" = "DE", agg_abr)))) +
      scale_linetype_manual(values = agg_lty[c("ens", agg_meths)],
                            labels = as.expression(parse(text = c("ens" = "DE", agg_abr)))) +
      facet_grid(rows = vars(metric), 
                 cols = vars(nn), 
                 scales = "free") +
      geom_hline(data = hline_plot, aes(yintercept = value), 
                 linetype = "dashed") +
      theme_bw() +
      theme(axis.text.x = element_text(size = 10))
    
    # Add upper and lower
    if(length(l_vec) > 0){ temp_plot <- temp_plot + 
      geom_hline(data = l_plot, aes(yintercept = value), 
                 linetype = "blank") }
    if(length(u_vec) > 0){ temp_plot <- temp_plot + 
      geom_hline(data = u_plot, aes(yintercept = value), 
                 linetype = "blank") }
    
    # Legend only once
    if(temp_nn == "DRN"){
      # make legend
      temp_plot <- temp_plot + 
        theme(legend.position = "top", 
              legend.justification = "center", 
              legend.text = element_text(size = rel(1.1)), 
              legend.title = element_blank()) +
        guides(colour = guide_legend(nrow = 1))
      
      # Get legend
      my_legend <- get_legend(temp_plot) 
    }
    
    # No legend in plot
    temp_plot <- temp_plot + theme(legend.position = "none")
    
    # Text on y-axis
    if(temp_nn != "HEN"){ temp_plot <- temp_plot + theme(strip.text.y = element_blank()) }
    
    # Different y-labels
    if(temp_nn == "DRN"){ temp_plot <- temp_plot + ylab("Value") 
    }else{ temp_plot <- temp_plot + ylab(element_blank()) }
    
    # Different x-labels
    if(temp_nn == "BQN"){ temp_plot <- temp_plot + xlab("Network ensemble size") 
    }else{ temp_plot <- temp_plot + xlab("") }
    
    # Assign name
    assign(x = paste0("plot_", temp_nn),
           value = temp_plot)
    rm(temp_plot)
  }
  
  # Plot together
  pdf_plot <- grid.arrange(my_legend, 
                           plot_DRN, plot_BQN, plot_HEN, 
                           ncol = 3, nrow = 2,
                           layout_matrix = rbind(c(1, 1, 1),
                                                 c(2, 3, 4)),
                           heights = c(0.2, 3))
  
  # Save as PDF
  ggsave(filename = paste0("cs_gg_panel_step", temp_step, ".pdf"),
         path = pdf_path,
         plot = pdf_plot,
         width = 10,
         height = 1.9*length(score_vec),
         scale = 0.8)
}

#### Case: Stepwise CRPSS boxplot panel ####
# For-Loop over forecast steps
for(temp_step in temp_step_vec){
  #### Initialization ####
  # Only forecast step
  df_step <- subset(df_scores, (fc_step == temp_step))
  
  #### Calculate quantities ####
  # List for matrices
  df_plot <- data.frame(nn = character(),
                        agg = character(),
                        n_ens = integer(),
                        crpss = numeric(),
                        stringsAsFactors = FALSE)
  
  # For-Loop over network variants, ensemble sizes and aggregation methods
  for(temp_nn in nn_vec){ for(i_ens in n_ens_vec){ for(temp_agg in agg_meths){
    # Get subset of scores data
    df_sub <- subset(df_step, (n_ens == i_ens) & 
                       (type == temp_agg) & (nn == temp_nn))
    
    # Get row of data frame
    i <- nrow(df_plot) + 1:nrow(df_sub)
    
    # Fill in data frame
    df_plot[i, "n_ens"] <- i_ens
    df_plot[i, "nn"] <- temp_nn
    df_plot[i, "agg"] <- temp_agg
    
    # Reference: Average score of ensemble members
    df_plot[i, "crpss"] <- 100*df_sub[["crpss"]]
  }}}
  
  #### Prepare data ####
  # Minimal lower/upper limit
  lu_vec <- range(df_plot[["crpss"]])
  
  # Rename networks
  df_plot[["nn"]] <- c("drn" = "DRN",
                       "bqn" = "BQN",
                       "hen" = "HEN")[df_plot[["nn"]]]
  
  # Change order of networks
  df_plot[["nn"]] <- factor(df_plot[["nn"]], 
                            levels = c("DRN", "BQN", "HEN"))
  
  # Rename aggregation methods
  df_plot[["agg"]] <- agg_abr[df_plot[["agg"]]]
  
  # Colors for aggregation methods
  agg_col_gg <- agg_col
  names(agg_col_gg) <- agg_abr[names(agg_col)]
  
  #### PDF ####
  # Make plot
  pdf_plot <- ggplot(df_plot, 
                     aes(x = n_ens, y = crpss, group = n_ens, 
                         col = agg, fill = agg)) + 
    facet_grid(rows = vars(nn),
               cols = vars(agg),
               labeller = label_parsed,
               scales = "free") +
    geom_boxplot(alpha = 0.5, 
                 outlier.size = 0.6,
                 na.rm = TRUE) +
    scale_color_manual(values = agg_col_gg[unique(df_plot[["agg"]])]) +
    scale_fill_manual(values = agg_col_gg[unique(df_plot[["agg"]])]) +
    coord_cartesian(ylim = lu_vec) +
    xlab("Network ensemble size") + 
    ylab("CRPSS in %") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15)) + 
    theme(legend.position = "none") +
    guides(colour = guide_legend(nrow = 1))
  
  # Add intercept
  if(min(lu_vec) < 0){ pdf_plot <- pdf_plot + geom_hline(yintercept = 0, 
                                                         linetype = "dashed",
                                                         color = agg_col["ens"]) }
  
  # Save as PDF
  ggsave(filename = paste0("cs_gg_panel_crpss_boxplots_step", temp_step, ".pdf"),
         path = pdf_path,
         plot = pdf_plot,
         width = 20,
         height = 15,
         scale = 0.7)
}
