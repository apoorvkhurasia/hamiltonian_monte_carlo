source("hmc.R")
library(latex2exp)
library(ggplot2)

# This only needs to run once and is needed to match the fonts used in slides
library(extrafont)

loadfonts(device = "postscript")
euler_theme <- theme_bw() +
    theme(
        text = element_text(size = 16, family = "Palatino"),
        legend.text = element_text(
            size = 10,
            face = "bold", family = "Palatino"
        ),
        legend.position = "bottom"
    )
theme_set(euler_theme)

# The log of 2-D Rosenbrock density function
log_rosenbrock <- function(x) {
    -((1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2) / 20
}

# The 2-D Rosenbrock density function itself
rosenbrock <- function(x) {
    exp(log_rosenbrock(x))
}

# The gradient of `log_rosenbrock`
grad_log_rosenbrock <- function(x) {
    return(-c(
        2 * (1 - x[1]) + 400 * x[1] * (x[2] - x[1]^2),
        -200 * (x[2] - x[1]^2)
    ) / 20)
}

# Use a simple mass (identity matrix)
mass <- diag(1, 2, 2)

# Run the simulation
res <- draw_samples_hmc(
    log_rosenbrock,
    grad_log_rosenbrock,
    2,
    10000,
    mass,
    integration_step_size = 0.05,
    num_steps_integration = 30
)

# Create a data frame from drawn samples
sample_df <- as.data.frame(res$samples)
names(sample_df) <- c("x1", "x2")

mesh <- expand.grid(
    x1 = seq(-4, 4, length.out = 200),
    x2 = seq(-2, 8, length.out = 200)
)
density <- cbind(mesh, prob = rosenbrock(mesh))
names(density) <- c("x1", "x2", "prob")
total_weight <- sum(density$prob)
density$prob <- density$prob / total_weight

# Scatter plot
ggplot() +
    geom_contour(
        data = density,
        aes(x = x1, y = x2, z = prob), alpha = 0.3
    ) +
    geom_point(
        data = sample_df, aes(x1, x2),
        alpha = 0.2, color = "black", size = 0.1
    ) +
    xlab(TeX("$x_1$")) +
    ylab(TeX("$x_2$")) +
    xlim(-4, 4) +
    ylim(-1, 8)

# Chain analysis
par(mfrow = c(2, 6))
ts.plot(res$samples[, 1],
    xlab = "Iterations",
    ylab = TeX("$x_1$")
)
acf(res$samples[, 1], main = "")
hist(res$samples[, 1], prob = TRUE, main = "", xlab = TeX("$x_1$"))

ts.plot(res$samples[, 2],
    xlab = "Iterations",
    ylab = TeX("$x_2$")
)
acf(res$samples[, 2], main = "")
hist(res$samples[, 2], prob = TRUE, main = "", xlab = TeX("$x_2$"))
cat(sprintf("Acceptance ratio = %.2f", res$accept_ratio))
