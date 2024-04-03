#' @importFrom graphics legend abline text contour
NULL

#' Plot contour for a bivariate normal distribution
#' 
#' @param mean1 Mean of the first normal distribution (on x-axis).
#' @param sd1 Standard deviation of thefirst normal distribution.
#' @param mean2 Mean of the second normal distribution (on y-axis).
#' @param sd2 Standard deviation of the second normal distribution.
#' @param cor12 Correlation in the bivariate normal distribution.
#' @param cov12 Covariance in the bivariate normal distribution. If not input,
#'     compute the covariance using the correlation and the standard deviations.
#' @param density Density level, i.e., probability enclosed by the ellipse.
#' @param length_out Number of values on the x-axis and on the y-axis to be
#'               evaluated; default to 101.
#' @param bty Argument passed to the `contour` function.
#' @param ... Additional arguments passed to \code{\link[graphics]{contour}}
#'
#' @return A plot showing the contour of the bivariate normal distribution on
#'     a two-dimensional space.
#' @examples
#' \dontrun{
#' contour_bvnorm(
#'   0.5, 1, 0.57, 1.03, cov12 = 0.8,
#'   xlab = bquote("Latent Composite" ~ (eta)),
#'   ylab = bquote("Observed Composite" ~ (italic(Z))),
#'   lwd = 2, col = "red", xlim = c(-3.0, 3.5),
#'   ylim = c(-2.97, 3.67)
#' )
#' }
contour_bvnorm <- function(mean1 = 0, sd1 = 1, mean2 = 0, sd2 = 1,
                           cor12 = 0, cov12 = NULL,
                           density = .95, length_out = 101,
                           bty = "L", ...) {
  # Error handling
  stopifnot(cor12 >= -1, cor12 <= 1)
  
  
  if (is.null(cov12)) cov12 <- cor12 * sd1 * sd2
  x_seq <- mean1 + seq(-3, 3, length.out = length_out) * sd1
  y_seq <- mean2 + seq(-3, 3, length.out = length_out) * sd2
  z <- outer(x_seq, y_seq, .bvnorm_kernel, mu_x = mean1, mu_y = mean2,
             sd_x = sd1, sd_y = sd2, cov_xy = cov12)
  contour(x_seq, y_seq, z, levels = qchisq(density, 2), drawlabels = FALSE,
          bty = bty, ...)
}

#' Contour plots for multiple groups.
#' 
#' \code{plot.PartInv} plots the contours for any number of groups.
#' 
#' @param x PartInv output.
#' @param labels By default, c("Reference", "Focal_1", ..., "Focal_g") where `g`
#'     is the number of groups.
#' @param which_result Whether to plot the partial or the strict invariance plot.
#' @param custom_colors Optional argument for specifying the colors of the 
#'     ellipses. 
#' @param ... Additional arguments.
#'@export
plot.PartInv <- function(x, labels = x[["labels"]],
                         which_result = c("pi", "mi"),
                         custom_colors = NULL, ...) {
    which_result <- match.arg(which_result)
    if (which_result == "pi") {
        plot_dat <- x$bivar_data
        cut_xi <- x$cutpt_xi
        cut_z <- x$cutpt_z
        summ <- x$summary
        title <- c("Partial Measurement Invariance")
    } else if (which_result == "mi") {
        summ <- x$summary_mi
        if (is.null(summ)) {
            stop("Strict invariance results not found. ",
                 "Please include `show_mi_result = TRUE` when running `PartInv()`.")
        }
        plot_dat <- x$bivar_data_mi
        cut_xi <- x$cutpt_xi_mi
        cut_z <- x$cutpt_z_mi
        title <- c("Strict Measurement Invariance")
    }

    # Determine the ranges of the y and x axis 
    x_lim <- y_lim <- c()
    n_g <- length(plot_dat$mn_xi) # number of groups
    
    for (i in seq_len(n_g)) {
      x_lim <- c(x_lim, c(plot_dat$mn_xi[i] + c(-3, 3) * plot_dat$sd_xi[i]))
      y_lim <- c(y_lim, plot_dat$mn_z[i] + c(-3, 3) * plot_dat$sd_z[i])
    }
    x_lim <- range(x_lim); y_lim <- range(y_lim)
    
    colorlist <-  c('#e6194b', '#4363d8', '#3cb44b', '#ffe119', '#f58231', 
                    '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
                    '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', 
                    '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', 
                    '#ffffff', '#000000') 
                    #https://sashamaps.net/docs/resources/20-colors/
    ltylist <- rep(c('twodash', 'longdash', 'dotdash', 'dashed', 'dotted'), 
                   length.out = n_g)
    if (!is.null(custom_colors)) { colorlist <- custom_colors }

    # Plot the ellipse for the reference group
    contour_bvnorm(plot_dat$mn_xi[1], plot_dat$sd_xi[1],
                   plot_dat$mn_z[1], plot_dat$sd_z[1],
                   cov12 = plot_dat$cov_z_xi[1],
                   xlab = bquote("Latent Composite" ~ (eta)),
                   ylab = bquote("Observed Composite" ~ (italic(Z))),
                   lwd = 2, col = colorlist[1], xlim = x_lim, ylim = y_lim,
                   main = title)
    # Add on the ellipses for the focal groups
    for (i in 2:n_g) {
      contour_bvnorm(plot_dat$mn_xi[i], plot_dat$sd_xi[i],
                     plot_dat$mn_z[i], plot_dat$sd_z[i],
                     cov12 = plot_dat$cov_z_xi[i],
                     add = TRUE, lwd = 2, col = colorlist[i], 
                     lty = ltylist[i]
                     )
    }
     legend("topleft", labels, lty = c("solid", ltylist[2:n_g]), 
            col = colorlist[1:n_g])
     abline(h = cut_z, v = cut_xi)
     x_cord <- rep(cut_xi + c(.8, -.8) * plot_dat$sd_xi[1], 2)
     y_cord <- rep(cut_z + c(.8, -.8) * plot_dat$sd_z[1], each = 2)
     text(x_cord, y_cord, c("A", "B", "D", "C"))
     if (n_g > 20) {
       warning("If you would like to plot the contours of more than 20 groups, 
               please provide a list of 20 color names.")
     }
}