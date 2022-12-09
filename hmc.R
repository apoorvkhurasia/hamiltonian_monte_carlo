#' Performs leapfrog integration
#' @param pos: The starting position vector
#' @param mom: The starting momentum vector
#' @param step_size: The step size to use
#' @param num_steps: Total number of steps needed for integration
#' @param gradient: The gradient function to use
#' @return The position and momentum vectors at the end of
#'         integration
leapfrog_integration <- function(pos,
                                 mom,
                                 step_size,
                                 num_steps,
                                 gradient) {
    for (i in 1:num_steps) {
        mom <- mom - step_size * gradient(pos) / 2
        pos <- pos + step_size * mom
        mom <- mom - step_size * gradient(pos) / 2
    }
    return(list(pos = pos, mom = mom))
}

#' Returns the kinetic energy given the mass and momentum.
#' @param inverse_mass_matrix: The inverse of the mass matrix
#' @param mom: The momentum vector
kinetic_energy <- function(inverse_mass_matrix, mom) {
    return(t(mom) %*% inverse_mass_matrix %*% mom) / 2
}

#' Draws samples from a distribution using HMC. Also
#' returns the acceptance ratio for diagnostics.
#' @param log_density The log of the density of the distribution
#'                       This function must return a vector of the
#'                       same length as the number of dimensions in
#'                       the distribution.
#' @param grad_log_density The gradient of `potential_func`
#'                       This function must return a vector of the
#'                       same length as the number of dimensions in
#'                       the distribution.
#' @param num_dimensions The number of dimensions in the distribution
#' @param num_samples The number of samples needed
#' @param mass_matrix The mass matrix to use
#'                    (must be symmetric and positive definite)
#' @param integrator The step size for integration (default 0.01)
#' @param num_steps_integration The number of steps to use in integration
#'                              (default 20)
draw_samples_hmc <- function(log_density,
                             grad_log_density,
                             num_dimensions,
                             num_samples,
                             mass_matrix,
                             integration_step_size,
                             num_steps_integration) {
    inv_mass <- solve(mass_matrix)
    samples <- matrix(0, nrow = num_samples, ncol = num_dimensions)
    rejections <- 0
    for (i in 2:num_samples) {
        pos <- samples[i - 1, ]
        mom <- t(chol(mass_matrix)) %*% rnorm(num_dimensions, 0, 1)
        initial_energy <- log_density(pos) + kinetic_energy(inv_mass, mom)
        z <- leapfrog_integration(
            pos, mom,
            integration_step_size,
            num_steps_integration,
            grad_log_density
        )
        pos <- z$pos
        mom <- z$mom
        ending_energy <- log_density(pos) + kinetic_energy(inv_mass, mom)
        boltzmann <- exp(-initial_energy + ending_energy)
        alpha <- min(1, boltzmann)
        # Accept if we fit the Boltzmann distribution
        if (runif(1) < alpha) {
            samples[i, ] <- pos
        } else {
            rejections <- rejections + 1
        }
    }
    return(list(samples = samples, accept_ratio = 1 - rejections / num_samples))
}
