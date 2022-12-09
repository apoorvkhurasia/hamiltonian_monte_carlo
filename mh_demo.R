library(latex2exp)
library(ggplot2)

# This only needs to run once and is needed to match the fonts used in slides
library(extrafont)

loadfonts(device = "postscript")
euler_theme <- theme_bw() +
    theme(
        text = element_text(size = 16, family = "Palatino"),
        legend.text = element_text(size = 10,
                                   face = "bold", family = "Palatino"),
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

burn_in <- 10000
num_samples <- 10000
total_iterations <- burn_in + num_samples
samples <- matrix(0, total_iterations, 2)
samples[1, ] <- c(0, 0)
rejections <- 0
for (i in 2:total_iterations) {
    proposal <- rnorm(2, samples[i - 1][1], sd=0.1)
    alpha <- min(1, rosenbrock(proposal)/rosenbrock(samples[i-1,]))
    if (runif(1) < alpha) {
        samples[i, ] <- proposal
    } else {
        samples[i, ] <- samples[i - 1, ]
        rejections <- rejections + 1
    }
}

mh_df <- as.data.frame(samples[(burn_in + 1):total_iterations,])
names(mh_df) <- c("x1", "x2")

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
        data = mh_df, aes(x1, x2),
        alpha = 0.2, color = "black", size = 0.1
    ) +
    xlab(TeX("$x_1$")) +
    ylab(TeX("$x_2$")) +
    xlim(-4, 4) +
    ylim(-1, 8)

# Chain analysis
par(mfrow = c(2, 6))
ts.plot(samples[, 1],
        xlab = "Iterations",
        ylab = TeX("$x_1$")
)
acf(samples[, 1], main = "")
hist(samples[, 1], prob = TRUE, main = "", xlab = TeX("$x_1$"))

ts.plot(samples[, 2],
        xlab = "Iterations",
        ylab = TeX("$x_2$")
)
acf(samples[, 2], main = "")
hist(samples[, 2], prob = TRUE, main = "", xlab = TeX("$x_2$"))

cat(sprintf("Acceptance ratio = %.2f",
            1 - rejections / num_samples))