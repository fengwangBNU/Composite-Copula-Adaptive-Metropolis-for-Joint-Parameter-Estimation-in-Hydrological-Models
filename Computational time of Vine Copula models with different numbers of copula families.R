# ============================================================
#Computational time of Vine Copula models with different numbers of copula families
# ============================================================



# ============================================================
# Preparation: data and marginal transformation
# ============================================================

set.seed(123)

n_sample <- 1000   # number of samples
d <- 5             # dimensionality (HYMOD parameters)

DFN <- x[1:n_sample, ]

# Transform marginals to U(0,1) using empirical normal scores
norm_x <- matrix(NA, n_sample, d)
for (t in 1:d) {
  norm_x[, t] <- pnorm(DFN[, t],
                       mean = mean(DFN[, t]),
                       sd   = sd(DFN[, t]))
}
# ============================================================
# Function: fit vine copula and record computational time
# ============================================================

fit_vine_with_timing <- function(norm_x, DFN, familyset) {
  
  # ---- Copula fitting ----
  t_start <- Sys.time()
  
  RVM <- CDVineCondFit(
    norm_x,
    Nx            = 5,
    familyset     = familyset,
    treecrit      = "BIC",
    type          = "CVine",
    selectioncrit = "AIC"
  )
  
  t_fit_end <- Sys.time()
  
  # ---- One-step sampling from fitted copula ----
  u_new <- CDVineCondSim(RVM, runif(1))
  
  y_new <- numeric(d)
  for (t in 1:d) {
    y_new[t] <- qnorm(u_new[t],
                      mean = mean(DFN[, t]),
                      sd   = sd(DFN[, t]))
  }
  
  t_sample_end <- Sys.time()
  
  return(list(
    fit_time    = as.numeric(difftime(t_fit_end, t_start, units = "secs")),
    sample_time = as.numeric(difftime(t_sample_end, t_fit_end, units = "secs"))
  ))
}


# ============================================================
# Candidate copula family sets
# ============================================================

family_sets <- list(
  "1 family"  = c(1),
  "2 families" = c(1, 3),
  "3 families" = c(1, 3, 4),
  "4 families" = c(1, 3, 4, 5),
  "5 families" = c(1, 3, 4, 5, 6),
  "10 families" = c(1, 2, 3, 4, 5, 6, 7, 13, 23, 24),
  "All (43)"    = NA
)


# ============================================================
# Run experiment
# ============================================================

results <- data.frame(
  FamilySet   = character(),
  FitTime_s  = numeric(),
  SampleTime_s = numeric(),
  stringsAsFactors = FALSE
)

for (name in names(family_sets)) {
  
  timing <- fit_vine_with_timing(
    norm_x  = norm_x,
    DFN     = DFN,
    familyset = family_sets[[name]]
  )
  
  results <- rbind(results, data.frame(
    FamilySet    = name,
    FitTime_s   = timing$fit_time,
    SampleTime_s = timing$sample_time
  ))
}

print(results)



# install.packages("ggplot2")
# install.packages("gridExtra")

library(ggplot2)
library(gridExtra)

# # Convert 'FamilySet' to a factor to maintain the correct order in plots
results$FamilySet <- factor(results$FamilySet, levels = results$FamilySet)

# Define a theme for clean, paper-style plots
theme_paper <- theme(
  panel.background = element_blank(),        # Remove internal background
  panel.grid = element_blank(),              # Remove internal grid lines
  axis.line = element_blank(),               # Remove axis lines inside the panel
  axis.ticks = element_line(color="black"),  # Keep ticks
  axis.text = element_text(family="Arial", size=12, color="black"),
  axis.title = element_text(family="Arial", size=12, color="black"),
  plot.background = element_blank(),         # No extra background
  panel.border = element_rect(color="black", fill=NA, size=1) # Four border lines
)

# --------------------------
# Subplot 1: FitTime
# --------------------------
p1 <- ggplot(results, aes(x=FamilySet, y=FitTime_s, group=1)) +
  geom_line(color="steelblue", size=1.2) +
  geom_point(color="steelblue", size=3) +
  labs(x="", y="Fit Time (s)") +
  theme_paper +
  theme(axis.text.x = element_blank())  # Remove x-axis text for top subplot

# --------------------------
# Subplot 2: SampleTime
# --------------------------
# Convert SampleTime to milliseconds for better visualization
results$SampleTime_ms <- results$SampleTime_s * 1000

p2 <- ggplot(results, aes(x=FamilySet, y=SampleTime_ms, group=1)) +
  geom_line(color="darkorange", size=1.2, linetype="dashed") +
  geom_point(color="darkorange", size=3, shape=17) +
  labs(x="Number of Copula Families", y="Sample Time (ms)") +
  theme_paper +
  theme(axis.text.x = element_text(angle=45, hjust=1))

# --------------------------
# Arrange subplots vertically
# --------------------------
grid.arrange(p1, p2, ncol=1)
# --------------------------
