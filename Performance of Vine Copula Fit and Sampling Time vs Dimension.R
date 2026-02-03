# ============================================================
# Performance of Vine Copula: Fit and Sampling Time vs Dimension
# ============================================================

# Install/load required packages
# install.packages("VineCopula")
library(VineCopula)
library(ggplot2)
library(gridExtra)

set.seed(123)

# -----------------------------
# Parameters
# -----------------------------
dims <- 3:20         # Dimensions to test
n_sample <- 1000      # Number of samples
copula_family <- 1    # Gaussian copula (family=1)

# Prepare a results data frame
results <- data.frame(
  Dimension = dims,
  FitTime_sec = NA,
  SampleTime_sec = NA
)

# -----------------------------
# Loop over dimensions
# -----------------------------
for (d in dims) {
  
  cat("Running dimension:", d, "\n")
  
  # Simulate uniform data for Vine copula fitting
  data <- matrix(runif(n_sample * d), nrow = n_sample, ncol = d)
  
  # Record vine copula fitting time
  fit_time <- system.time({
    vine_fit <- RVineStructureSelect(
      data,
      familyset = copula_family,   # Single copula family
      selectioncrit = "AIC"
    )
  })
  results$FitTime_sec[d - 2] <- fit_time["elapsed"]
  
  # Record time to generate one sample
  sample_time <- system.time({
    sample <- RVineSim(1, vine_fit)
  })
  results$SampleTime_sec[d - 2] <- sample_time["elapsed"]
}

# -----------------------------
# Convert sample time to milliseconds for plotting
# -----------------------------
results$SampleTime_ms <- results$SampleTime_sec * 1000

# -----------------------------
# Plotting setup
# -----------------------------
theme_paper <- theme(
  panel.background = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line(color="black"),
  axis.ticks = element_line(color="black"),
  axis.text = element_text(family="Arial", size=12, color="black"),
  axis.title = element_text(family="Arial", size=12, color="black"),
  panel.border = element_rect(color="black", fill=NA, size=1)
)

# -----------------------------
# Plot 1: Fit Time
# -----------------------------
p1 <- ggplot(results, aes(x=Dimension, y=FitTime_sec)) +
  geom_line(color="steelblue", size=1.2) +
  geom_point(color="steelblue", size=3) +
  labs(x="Dimension", y="Fit Time (s)") +
  theme_paper

# -----------------------------
# Plot 2: Sample Time
# -----------------------------
p2 <- ggplot(results, aes(x=Dimension, y=SampleTime_ms)) +
  geom_line(color="darkorange", size=1.2, linetype="dashed") +
  geom_point(color="darkorange", size=3, shape=17) +
  labs(x="Dimension", y="Sample Time (ms)") +
  theme_paper

# -----------------------------
# Arrange plots vertically
# -----------------------------
grid.arrange(p1, p2, ncol = 1)
