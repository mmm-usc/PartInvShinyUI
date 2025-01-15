#' @title 
#' Unnest list elements. 
#' 
#' @name 
#' unnest_list
#'
#' @description
#' \code{unnest_list} takes in a list object and returns an unnested list
#'  of lists
#'  
#' @param ins list object (e.g., lavaan CFA fit)
#' 
#' @return The output will be a list of lists.
#'
#' @export
unnest_list <- function(ins) {
        nms <- names(ins[[1]])
        num_gr <- length(ins)
        out <- rep(list(vector(mode = "list", length = num_gr)), length(nms))
        names(out) <- nms
        for (i in seq_along(ins)) {
            for (nm in nms) {
                out[[nm]][[i]] <- ins[[i]][[nm]]
            }
        }
        out
    }

#' @title 
#' Convert acronym of composite index to full form. 
#' 
#' @name 
#' lab_cai
#'
#' @description
#' \code{lab_cai} takes in a string (PS, SR, SE, or SP) and returns the string 
#' for the full form.
#'  
#' @param cai A two-letter string indicating the classification accuracy index.
#' 
#' @return The output will be a string.
lab_cai <- function(cai) {
  out <- ""
  if(cai == "PS") {
    out <- "Proportion selected"
  }
  if(cai == "SR") {
    out <- "Success ratio"
  }
  if(cai == "SE") {
    out <- "Sensitivity"
  }
  if(cai == "SP") {
    out <- "Specificity"
  }
  out
}

#' @title 
#' Extract and format parameter values for `PartInv`.
#' 
#' @name 
#' format_cfa_partinv
#'
#' @description
#' \code{format_cfa_partinv} takes in a lavaan CFA fit object and a component
#'  and returns the necessary inputs for PartInv in a list
#'  
#' @param obj lavaan CFA output
#' @param comp a string indicating the lavaan object component of interest
#'             e.g., "se", "est"
#'  
#' @return The output will be a list of 5 elements:
#'    \item{nu}{A list of length `g` containing `1 x n` measurement intercept
#'     vectors where `g` is the number of groups and `n` is the number of items
#'     in the scale.}
#'    \item{alpha}{A list of length `g` containing `1 x d` latent factor mean
#'     vectors where `g` is the number of groups and `d` is the number of 
#'     latent dimensions.}
#'    \item{lambda}{A list of length `g` containing `n x d` factor loading
#'     matrices where `g` is the number of groups, `d` is the number of 
#'     latent dimensions, and `n` is the number of items in the scale.}
#'    \item{psi}{A list of length `g` containing `d x d` latent factor
#'     variance-covariance matrices where `g` is the number of groups and `d` 
#'     is the number of latent dimensions.}
#'    \item{theta}{A list of length `g` containing `1 x n` vectors or `n x n`
#'     matrices of unique factor variances and covariances, where `g` is the
#'     number of groups and `n` is the number of items in the scale.}
#'     
#' @export
format_cfa_partinv <- function(obj, comp) {
  ins <- lavaan::lavInspect(obj, what = comp)
  num_gr <- length(ins)
  
  psi_matrices <- lambda_matrices <- alpha_list <- nu_list <- 
    theta_list <- vector(mode = "list", length = num_gr)
  
  # Extract and format the parameters for each group
  for (i in seq_len(num_gr)) {
    psi_matrices[[i]] <- ins[[i]]$psi
    lambda_matrices[[i]] <- ins[[i]]$lambda
    alpha_list[[i]] <- ins[[i]]$alpha
    nu_list[[i]] <- ins[[i]]$nu
    theta_list[[i]] <- ins[[i]]$theta
  }
  return(list("lambda" = lambda_matrices, 
              "theta" = theta_list, 
              "psi" = psi_matrices, 
              "nu" = nu_list, 
              "alpha" = alpha_list))
}


#' Compute the mean, standard deviation, and covariance of latent and observed
#' variables.
#'
#' \code{mn_sd_cov} is a helper function that computes the mean, 
#' standard deviation, and covariance of latent and observed variables for each 
#' group.
#' @param weights_item A vector of item weights.
#' @param weights_latent A vector of latent factor weights.
#' @param alpha A list of length `g` containing `1 x d` latent factor mean 
#'    vectors where `g` is the number of groups and `d` is the number of latent 
#'    dimensions. The first element is assumed to belong to the reference group.
#' @param psi A list of length `g` containing `d x d` latent factor 
#'    variance-covariance matrices where `g` is the number of groups and `d` is 
#'    the number of latent dimensions. The first element is assumed to belong 
#'    to the reference group.
#' @param lambda A list of length `g` containing `n x d` factor loading matrices 
#'    where `g` is the number of groups, `d` is the number of latent dimensions, 
#'    and `n` is the number of items in the scale. The first element is assumed 
#'    to belong to the reference group.
#' @param nu A list of length `g` containing `1 x n` measurement intercept
#'    vectors where `g` is the number of groups and `n` is the number of items 
#'    in the scale. The first element is assumed to belong to the reference 
#'    group.
#' @param theta A list of length `g` containing `1 x n` vectors or `n x n` 
#'    matrices of unique factor variances and covariances, where `g` is the 
#'    number of groups and `n` is the number of items in the scale. The first 
#'    element is assumed to belong to the reference group.
#' @return The output will be a list of 5 elements:
#'    \item{mn_z}{Mean of the observed variable.}
#'    \item{sd_z}{Standard deviation of the observed variable.}
#'    \item{mn_xi}{Mean of the latent variable.}
#'    \item{sd_xi}{Standard deviation of the latent variable.}
#'    \item{cov_z_xi}{Covariance of the latent and observed variables.}
mn_sd_cov <- function(weights_item, weights_latent, alpha, psi, lambda, nu, 
                      theta) {
  mn_z <- sd_z <- mn_xi <- sd_xi <- cov_z_xi <- NULL
  for (i in seq_along(alpha)) {
    mn_z[i] <- c(crossprod(weights_item, nu[[i]] + lambda[[i]] %*% alpha[[i]]))
    sd_z[i] <- c(sqrt(crossprod(weights_item, lambda[[i]] %*% psi[[i]] %*% 
                                  t(lambda[[i]]) + theta[[i]]) %*%
                        weights_item))
    names(alpha) <- names(psi) <- NULL
    mn_xi[i] <- c(crossprod(weights_latent, alpha[[i]]))
    sd_xi[i] <- c(sqrt(crossprod(weights_latent, psi[[i]]) %*% weights_latent))
    cov_z_xi[i] <- c(crossprod(weights_item, lambda[[i]] %*% psi[[i]]) %*% 
                    weights_latent)
  } 
  return(list(mn_z = mn_z, sd_z = sd_z, mn_xi = mn_xi, sd_xi = sd_xi,
              cov_z_xi = cov_z_xi))
}
#' Compute summary statistics.
#'
#' \code{compute_cai} computes classification accuracy indices.
#' @param weights_item A vector of item weights.
#' @param weights_latent A vector of latent factor weights.
#' @param alpha A list of length `g` containing `1 x d` latent factor mean 
#'    vectors where `g` is the number of groups and `d` is the number of latent 
#'    dimensions. The first element is assumed to belong to the reference group.
#' @param psi A list of length `g` containing `d x d` latent factor 
#'    variance-covariance matrices where `g` is the number of groups and `d` is 
#'    the number of latent dimensions. The first element is assumed to belong 
#'    to the reference group.
#' @param lambda A list of length `g` containing `n x d` factor loading matrices 
#'    where `g` is the number of groups, `d` is the number of latent dimensions, 
#'    and `n` is the number of items in the scale. The first element is assumed 
#'    to belong to the reference group.
#' @param nu A list of length `g` containing `1 x n` measurement intercept
#'    vectors where `g` is the number of groups and `n` is the number of items 
#'    in the scale. The first element is assumed to belong to the reference 
#'    group.
#' @param theta A list of length `g` containing `1 x n` vectors or `n x n` 
#'    matrices of unique factor variances and covariances, where `g` is the 
#'    number of groups and `n` is the number of items in the scale. The first 
#'    element is assumed to belong to the reference group.
#' @param pmix List of length `g` containing the mixing proportions of each 
#'    group.
#' @param propsel Proportion of selection. If missing, computed using `cut_z`.
#' @param labels A character vector with `g` elements to label the reference
#'    and focal groups on the plot, where `g` is the number of groups.
#' @param cut_z Pre-specified cutoff score on the observed composite. This
#'    argument is ignored when `propsel` has input.
#' @param is_mi Whether summary statistics should be computed for strict vs. 
#'    partial measurement invariance. `FALSE` by default (partial).
#' @return The output will be a list of 5 elements:
#'    \item{propsel}{Proportion selected.}
#'    \item{cutpt_xi}{Cut point on the latent variable.}
#'    \item{cutpt_z}{Cut point on the observed variable.}
#'    \item{summary}{Summary statistics.}
#'    \item{bivar_data}{The mean, standard deviation, and covariance of latent 
#'     and observed variables for each group.}
compute_cai <- function(weights_item, weights_latent, alpha, psi, lambda, nu, 
                        theta, pmix, propsel, labels, cut_z = NULL, num_g,
                        is_mi = FALSE) {
  
  lst <- mn_sd_cov(weights_item, weights_latent, alpha, psi, lambda, nu, theta)

  if (!is.null(propsel)) {  # if there is an input for selection proportion
    # compute the cut score using qnormmix based on input selection proportion
    cut_z <- qnormmix_mult(propsel, means = lst$mn_z, sds = lst$sd_z,
                           pmix = pmix, lower.tail = FALSE)
  } else if (!is.null(cut_z) && is.null(propsel)) {
    # compute the selection proportion using pnormmix based on the cutoff value
    propsel <- pnormmix_mult(cut_z, lst$mn_z, lst$sd_z, pmix = pmix,
                             lower.tail = FALSE)
  }
  
  # compute the threshold for the latent variable based on the selection 
  # proportion provided by the user/computed using cut_z
  cut_xi <- qnormmix_mult(propsel, lst$mn_xi, lst$sd_xi, pmix = pmix,
                          lower.tail = FALSE)
  
  # computing summary statistics
  #CAIs <- matrix(ncol = ifelse(is_mi, num_g, num_g + num_g - 1) , nrow = 8) 
  CAIs <- matrix(ncol = (num_g + num_g - 1), nrow = 8) 
  
  for (i in seq_along(1:num_g)) {
    CAIs[, i] <- .partit_bvnorm(cut_xi, cut_z, lst$mn_xi[[i]], lst$sd_xi[[i]],
                                lst$mn_z[[i]], lst$sd_z[[i]],
                                cov12 = lst$cov_z_xi[[i]])
  }
  
  # Store mean, sd, cov values for the obs/latent variables
  zf_par <- list(mn_xi = lst$mn_xi, sd_xi = lst$sd_xi, mn_z = lst$mn_z,
                 sd_z = lst$sd_z, cov_z_xi = lst$cov_z_xi)
   dat <- data.frame(CAIs,
                    row.names = c("A (true positive)", "B (false positive)",
                                  "C (true negative)", "D (false negative)",
                                  "Proportion selected", "Success ratio",
                                  "Sensitivity", "Specificity"))
   nms <- labels #c("Reference", paste0("Focal_", 1:(num_g - 1)))
   
 ### if (!is_mi) {
    # selection indices for the focal group if its distribution matches the
    # distribution of the reference group (Efocal)
    mn_z_Ef <- sd_z_Ef <- cov_z_xi_Ef <- vector(mode = "list")
    for (i in 2:num_g) {
      mn_z_Ef[i - 1] <- c(crossprod(weights_item, nu[[i]] + lambda[[i]]
                                    %*% alpha[[1]]))
      sd_z_Ef[i - 1] <- c(sqrt(crossprod(weights_item, lambda[[i]] %*% psi[[1]]
                                         %*% t(lambda[[i]]) + theta[[i]])
                               %*% weights_item))
      cov_z_xi_Ef[i - 1] <- c(crossprod(weights_item, lambda[[i]] %*% psi[[1]]) 
                              %*% weights_latent)
      
      CAIs[, i + num_g - 1] <- .partit_bvnorm(cut_xi, cut_z, lst$mn_xi[[1]],
                                              lst$sd_xi[[1]], mn_z_Ef[[i - 1]],
                                              sd_z_Ef[[i - 1]],
                                              cov12 = cov_z_xi_Ef[[i - 1]])
      nms <- c(labels, paste0("E_R(", labels[2:length(labels)], ")"))
    }
###  }
   dat <- 
     data.frame(CAIs, row.names = c("A (true positive)", "B (false positive)",
                                    "C (true negative)", "D (false negative)",
                                    "Proportion selected", "Success ratio",
                                    "Sensitivity", "Specificity"))
   names(dat) <- nms
   
   list(propsel = propsel, cutpt_xi = cut_xi, cutpt_z = cut_z, summary = dat,
        bivar_data = zf_par)
}

#' Distribution function (pdf) of a mixture of two normal distributions. 
#' 
#' \code{pnormmix} returns the cumulative probability of q or \eqn{1 - q} on 
#'  the mixture normal distribution.
#' 
#' @param q A vector of quantiles.
#' @param mean1 Mean of the first normal distribution.
#' @param sd1 Standard deviation of the first normal distribution.
#' @param mean2 Mean of the second normal distribution.
#' @param sd2 Standard deviation of the second normal distribution.
#' @param pmix1 Mixing proportion for the first distribution. Should be a 
#'   number in the range (0, 1).
#' @param lower.tail A logical scalar; if TRUE (default), probabilities are 
#' \eqn{P[X <= x]}; otherwise, \eqn{P[X > x]}.
#' @return The output will be the cumulative probability of q or \eqn{1 - q} on 
#'   the mixture normal distribution.
#' @examples
#' \dontrun{
#' pnormmix(1, 0, 3.1, 1.7, 3.1, lower.tail = FALSE)
#' }

pnormmix <- function(q, mean1 = 0, sd1 = 1, mean2 = 0, sd2 = 1, pmix1 = 0.5,
                     lower.tail = TRUE) {
  stopifnot(pmix1 > 0, pmix1 < 1)
  as.vector(c(pmix1, 1 - pmix1) %*%
              sapply(q, pnorm, mean = c(mean1, mean2), sd = c(sd1, sd2),
                     lower.tail = lower.tail))
}


pnormmix_mult <- function(q, means = 0, sds = 1, pmix = NULL, 
                          lower.tail = TRUE) {

  stopifnot("Provide mixing proportions between 0 and 1." = all(pmix > 0,
                                                                pmix < 1))
  as.vector(pmix %*%
              sapply(q, pnorm, mean = means, sd = sds, lower.tail = lower.tail))
}

#' Quantile function of a mixture of two normal distributions. 
#' 
#' \code{qnormmix} returns the quantile corresponding to \eqn{p} or \eqn{1 - q}  
#' on the mixture normal distribution.
#' 
#' @param p A vector of probabilities.
#' @param mean1 Mean of the first normal distribution.
#' @param sd1 Standard deviation of the first normal distribution.
#' @param mean2 Mean of the second normal distribution.
#' @param sd2 Standard deviation of the second normal distribution.
#' @param pmix1 Mixing proportion for the first distribution. Should be a 
#'   number in the range (0, 1).
#' @param lower.tail A logical scalar; if TRUE (default), probabilities are 
#' \eqn{P[X <= x]}; otherwise, \eqn{P[X > x]}. 
#' @return The output will be the quantile corresponding to p or \eqn{1 - q} on 
#' the mixture normal distribution.
#' @examples
#' \dontrun{
#' qnormmix(0.8, 0, 3.1, 1.7, 0.5, lower.tail = FALSE)
#' }

qnormmix <- function(p, mean1 = 0, sd1 = 1, mean2 = 0, sd2 = 1, pmix1 = 0.5,
                     lower.tail = TRUE) {
  stopifnot(pmix1 > 0, pmix1 < 1, p >= 0, p <= 1)
  f <- function(x) {
    (pnormmix(x, mean1, sd1, mean2, sd2, pmix1, lower.tail) - p)^2
    }
  start <- as.vector(c(pmix1, 1 - pmix1) %*%
                       sapply(p, qnorm, c(mean1, mean2), c(sd1, sd2),
                              lower.tail = lower.tail))
  nlminb(start, f)$par
}

qnormmix_mult <- function(p, means = c(0), sds = 1, pmix = NULL, 
                          lower.tail = TRUE) {
  stopifnot("Provide mixing proportions between 0 and 1." = 
              all(pmix > 0, pmix < 1,p >= 0, p <= 1))
  f <- function(x) (pnormmix_mult(x, means, sds, pmix, lower.tail) - p)^2
  start <- as.vector(pmix %*% sapply(p, qnorm, means, sds, 
                              lower.tail = lower.tail))
  nlminb(start, f)$par
}

#' Helper function for computing the kernel for bivariate normal density. 
#' 
#' \code{.bvnorm_kernel} returns the kernel for bivariate normal density
#' 
#' @param x A normal distribution
#' @param mu_x Mean of the normal distribution x.
#' @param y A normal distribution
#' @param sd_x Standard deviation of the normal distribution x.
#' @param mu_y Mean of the normal distribution y.
#' @param sd_y Standard deviation of the normal distribution y.
#' @param cov_xy covariance between x and y
#' 
#' @return The output will be the quantile corresponding to p or \eqn{1 - q} on 
#' the mixture normal distribution.
#' @examples
#' \dontrun{
#' .bvnorm_kernel(x = -2.50, y = -2.52, mu_x = 1, mu_y = 0.57, sd_x = 1,
#' sd_y = 1.03, cov_xy = 0.8)
#' }

.bvnorm_kernel <- function(x, y, mu_x = 0, mu_y = 0, sd_x = 1, sd_y = 1, 
                           cov_xy = 0) {
  
  cor <- cov_xy / sd_x / sd_y
  numer <- (x - mu_x)^2 / sd_x^2 + (y - mu_y)^2 / sd_y^2 - 
    2 * cor * (x - mu_x) * (y - mu_y) / sd_x / sd_y
  numer / (1 - cor^2)
}

#' Computing summary statistics from a selection approach
#' 
#' \code{.partit_bvnorm} returns a table of selection accuracy indices
#' 
#' @param cut1 Cut score based on the latent score
#' @param cut2 Cut score based on the observed score
#' @param mean1 Mean of first normal distribution (on x-axis).
#' @param sd1 Standard deviation of first normal distribution.
#' @param mean2 Mean of second normal distribution (on y-axis).
#' @param sd2 Standard deviation of second normal distribution.
#' @param cor12 Correlation in the bivariate normal.
#' @param cov12 Covariance in the bivariate normal. If not input, compute the
#'          covariance using the correlation and the standard deviations.
#' 
#' @return A table of selection accuracy indices
#' @examples 
#' \dontrun{
#' .partit_bvnorm(cut1 = 2, cut2 = 2, mean1 = 0, sd1 = 1, 
#'                mean2 = 1.53, sd2 = 0.89, cov12 = 0.80)
#' }

.partit_bvnorm <- function(cut1, cut2, mean1 = 0, sd1 = 1, mean2 = 0, sd2 = 1, 
                           cor12 = 0, cov12 = cor12 * sd1 * sd2) {
  Sigma <- matrix(c(sd1^2, cov12, cov12, sd2^2), nrow = 2)
  C <- mnormt::pmnorm(c(cut1, cut2), c(mean1, mean2), Sigma)
  B <- pnorm(cut1, mean1, sd1) - C
  D <- pnorm(cut2, mean2, sd2) - C
  A <- 1 - B - C - D
  propsel <- A + B
  success_ratio <- A / propsel
  sensitivity <- A / (A + D)
  specificity <- C / (C + B)
  c(A, B, C, D, propsel, success_ratio, sensitivity, specificity)
}


is_symmetric_posdef <- function(x, tol = 1e-08) {
  # Borrow from matrixcalc::is.positive.definite()
  if (!isSymmetric(x)) return(FALSE)
  eigenvalues <- eigen(x, only.values = TRUE)$values
  all(eigenvalues >= tol)
}
