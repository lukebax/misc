# =========================
# MEDIAN • HS 95% CI + BOOTSTRAP SE (BCa) — TEMPLATE SCRIPT
# =========================
# Purpose
#   • Compute a robust central tendency (the sample median) and:
#       (1) a 95% confidence interval (CI) via the Hettmansperger–Sheather (HS) interpolation method; and
#       (2) a standard error (SE) for the median via BCa bootstrap (square it to get a variance).
#
# Why these methods?
#   • HS CI: nonparametric, distribution-free, and small-n friendly. It interpolates between order statistics,
#     smoothing away the discreteness of “exact” rank intervals, typically giving near-nominal 95% coverage
#     with shorter, more stable intervals at small n. This CI is for visualisation / display.
#   • Bootstrap SE: the BCa bootstrap is a principled, nonparametric way to quantify the sampling variability
#     of the median (handles skew/heavy tails). This SE (and its square) is for downstream inverse-variance weighting.
#
# Small-n guidance
#   • Two-sided 95% distribution-free median CIs are sensible from n ≥ 6 (discreteness dominates below this).
#   • Bootstrap SE for the median is computable from n ≥ 3 and stabilizes by about n ≈ 6–10.
#   • Ties/discrete data can make any nonparametric CI a bit “chunky”; HS helps by smoothing endpoints.
#
# Outputs
#   • Console: n, sample median, HS 95% CI, bootstrap SE, and variance (SE^2).
#   • Plot: one figure with two side-by-side panels:
#       Left  = raw data (faded) + median (black dot) + HS 95% CI (green).
#       Right = single green bar showing the bootstrap SE, with its numeric value annotated above the bar.
#
# Adapting this template to real data
#   • Replace the “DATA BLOCK — simulated” with the CSV stub (uncomment and point to your file/column).
#   • Keep the “CORE ANALYSIS BLOCK” unchanged; it works for any numeric vector x.
#   • The “CONSOLE SUMMARY” and “PLOT BLOCK” are optional; keep or remove as needed.
#
# Dependencies
#   • quantileCI  (HS CI): install from R-universe or GitHub.
#       install.packages('quantileCI',
#                        repos = c('https://hoehleatsu.r-universe.dev','https://cloud.r-project.org'))
#     # or: remotes::install_github('hoehleatsu/quantileCI')
#   • boot        (bootstrap / BCa)
#   • ggplot2     (graphics)
#   • patchwork   (side-by-side plot layout)
# =========================


# 0) PACKAGE SETUP -----------------------------------------------------------
suppressPackageStartupMessages({
  library(quantileCI)   # HS interpolation CI for the median
  library(boot)         # bootstrap machinery (boot, boot.ci)
  library(ggplot2)      # plotting
  library(patchwork)    # layout for side-by-side panels
})


# 1) CONFIG BLOCK — USER KNOBS (EDIT HERE) ----------------------------------
reproducible <- FALSE    # TRUE → identical runs; FALSE → new random draw each run
seed         <- 12345    # only used when reproducible = TRUE
n            <- 6        # default small-n; change as needed
alpha        <- 0.05     # 95% CI
B            <- 10000    # bootstrap resamples (stable SE for small n)

# Plot colors / cosmetics
col_green   <- "#00733E" # readable green for “recommended” elements
col_grey_pt <- "#AFAFAF" # grey for raw points
pt_alpha    <- 0.6       # transparency for raw data points
jitter_w    <- 0.06      # horizontal jitter width for raw points


# 2) DATA BLOCK — CHOOSE ONE (SIMULATED by default; CSV STUB provided) ------
## A) Simulated demo data (default)
if (reproducible) set.seed(seed)
x <- rt(n, df = 3)       # heavy-tailed sample (Student-t with 3 df)

## B) Real data from CSV (uncomment and edit path/column)
# library(readr)
# dat <- readr::read_csv("path/to/your_file.csv")
# x   <- dat$your_numeric_column
# n   <- length(x)  # optional: overwrite n to reflect your data length


# 3) CORE ANALYSIS BLOCK — NO HELPER FUNCTIONS (PACKAGE CALLS ONLY) ---------
## 3.1) Point estimate: sample median
med <- median(x)

## 3.2) HS 95% CI for the median (returns c(lower, upper))
hs_ci <- quantileCI::median_confint_hs(x, conf.level = 1 - alpha, interpolate = TRUE)
hs_lower <- hs_ci[1]
hs_upper <- hs_ci[2]

## 3.3) Bootstrap SE of the sample median (BCa bootstrap → replicates; SE = sd of medians)
boot_obj <- boot::boot(
  data = x,
  statistic = function(d, i) median(d[i]),
  R = B
)
se_boot <- sd(boot_obj$t[, 1])   # bootstrap SE (principled, from replicates)
var_boot <- se_boot^2            # variance for inverse-variance weighting


# 4) CONSOLE SUMMARY BLOCK — OPTIONAL ---------------------------------------
cat("\n=== Median (HS CI) + Bootstrap SE (BCa) ===\n")
cat(sprintf("n = %d   alpha = %.2f   B = %d\n", length(x), alpha, B))
cat(sprintf("Median (sample): %g\n", med))
cat(sprintf("HS 95%% CI (median): [%g, %g]\n", hs_lower, hs_upper))
cat(sprintf("Bootstrap SE (median): %g   Variance: %g\n", se_boot, var_boot))


# 5) PLOT BLOCK — OPTIONAL (SIDE-BY-SIDE PANELS) ----------------------------
## Left: raw data (faded) + median dot + HS 95% CI (green vertical line with caps)
df_raw <- data.frame(
  xj = 1 + runif(length(x), -jitter_w, jitter_w),
  y  = x
)

p_left <-
  ggplot(df_raw, aes(x = xj, y = y)) +
  geom_point(color = col_grey_pt, alpha = pt_alpha) +
  
  # Median dot (single item; use annotate to avoid inheriting df_raw rows)
  annotate("point", x = 1, y = med, color = "black", size = 3) +
  
  # HS CI vertical line (single item)
  annotate("segment",
           x = 1, xend = 1, y = hs_lower, yend = hs_upper,
           color = col_green, linewidth = 1.1) +
  
  # HS CI caps (single items)
  annotate("segment",
           x = 1 - 0.02, xend = 1 + 0.02, y = hs_lower, yend = hs_lower,
           color = col_green, linewidth = 1.1) +
  annotate("segment",
           x = 1 - 0.02, xend = 1 + 0.02, y = hs_upper, yend = hs_upper,
           color = col_green, linewidth = 1.1) +
  
  scale_x_continuous(limits = c(0.8, 1.2)) +
  labs(
    title = "Raw data, median • HS 95% CI",
    subtitle = "HS CI (green) is for display",
    x = NULL, y = "Value"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_blank(),
        panel.grid.minor = element_blank())

## Right: single green bar for bootstrap SE with numeric label above the bar
df_se <- data.frame(method = factor("Bootstrap SE", levels = "Bootstrap SE"),
                    se = se_boot)

p_right <-
  ggplot(df_se, aes(x = method, y = se)) +
  geom_col(fill = col_green, width = 0.6) +
  geom_text(aes(label = sprintf("SE = %.3g", se), y = se * 1.04),
            vjust = 0, color = "black") +
  labs(
    title = "Standard error of the median (BCa bootstrap)",
    subtitle = "BCa SE (variance = SE^2) is for inverse-variance weighting",
    x = NULL, y = "SE"
  ) +
  expand_limits(y = se_boot * 1.12) +  # headroom so label doesn’t clip
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

## Combine panels horizontally (patchwork)

## Combine panels horizontally (patchwork)
p_right <- p_right + theme(plot.margin = margin(5.5, 30, 5.5, 5.5))
p_combined <- p_left | p_right
p_combined

## Save the figure:
ggsave("median_ci_se.png", plot = p_combined, width = 10, height = 4.5, dpi = 300)

