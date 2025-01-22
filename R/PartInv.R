#' @importFrom stats qchisq pnorm qnorm nlminb
#' @importFrom mnormt pmnorm
NULL

#' Classification Accuracy Analysis with MCAA
#'
#' \code{PartInv, PartInvMulti_we} evaluates partial measurement invariance
#' using the multidimensional classification accuracy analysis framework
#' (Lai & Zhang, 2022), which is an extension of Millsap & Kwok's (2004)
#' approach.
#'
#' @param cfa_fit CFA model output from lavaan.
#' @param propsel Proportion of selection. If missing, computed using `cut_z`.
#' @param cut_z Pre-specified cutoff score on the observed composite. Ignored 
#'     when `propsel` has input.
#' @param weights_item A vector of item weights.
#' @param weights_latent A vector of latent factor weights.
#' @param alpha A list of length `g` containing `1 x d` latent factor mean
#'     vectors where `g` is the number of groups and `d` is the number of latent
#'     dimensions. The first element is assumed to belong to the reference group.
#' @param psi A list of length `g` containing `d x d` latent factor
#'     variance-covariance matrices where `g` is the number of groups and `d` is
#'     the number of latent dimensions. The first element is assumed to belong
#'     to the reference group.
#' @param lambda A list of length `g` containing `n x d` factor loading matrices
#'     where `g` is the number of groups, `d` is the number of latent dimensions,
#'     and `n` is the number of items. The first element is assumed
#'     to belong to the reference group.
#' @param nu A list of length `g` containing `1 x n` measurement intercept
#'     vectors where `g` is the number of groups and `n` is the number of items
#'    . The first element is assumed to belong to the reference
#'     group.
#' @param theta A list of length `g` containing `1 x n` vectors or `n x n`
#'     matrices of unique factor variances and covariances, where `g` is the
#'     number of groups and `n` is the number of items. The first
#'     element is assumed to belong to the reference group.
#' @param pmix List of length `g` containing the mixing proportions of each
#'     group. If `NULL`, defaults to `1/g` for each group (i.e., equal sizes).
#' @param plot_contour Logical; whether the contour of the populations should be
#'     plotted; `TRUE` by default.
#' @param show_mi_result If \code{TRUE}, perform selection accuracy analysis
#'     for both the input parameters and the implied parameters based on a
#'     strict invariance model, with common parameter values as weighted
#'     averages of the input values using `pmix`.
#' @param labels A character vector with `g` elements to label the reference
#'     and focal groups on the plot, where `g` is the number of groups. If `NULL`
#'     groups are labeled as 'Reference' and 'Focal_1' through 'Focal_(g-1)'.
#' @param custom_colors Optional argument for specifying group colors.
#' @param reference Optional argument for specifying the reference group.
#' Currently only functional when cfa_fit is provided. If providing parameter 
#' estimates instead, reorder estimates such that the first estimates belong to 
#' the reference group.
#' @param ... Other arguments for \code{\link[graphics]{contour}}.
#' @param alpha_r,alpha_f,nu_r,nu_f,Theta_r,Theta_f,psi_r,psi_f,lambda_r,lambda_f,phi_r,phi_f,tau_r,tau_f,kappa_r,kappa_f,pmix_ref
#'     Deprecated; included for backward compatibility. With two groups, '_r' 
#'     and '_f' suffixes refer to the reference group and the focal group.
#' @return The output will be a list of six elements and a plot if
#'     \code{plot_contour == TRUE}:
#'         \item{propsel}{Proportion selected.}
#'         \item{cutpt_xi}{Cut point on the latent scale (xi).}
#'         \item{cutpt_z}{Cut point on the observed scale (Z).}
#'         \item{summary}{A `8 x (g + g - 1)` table, with columns representing
#'             the reference and `g - 1` focal groups, and  the expected results
#'             if the latent distribution of `g - 1` focal group match the
#'             reference group. The rows represent probabilities of true
#'             positive (A), false positive (B), true negative (C), false
#'             negative (D); proportion selected, success ratio, sensitivity,
#'             and specificity.}
#'         \item{bivardata}{List of length `5` containing `1 x g` vectors of
#'             latent and observed means, standard deviations, and covariances
#'             computed for each groups.}
#'         \item{ai_ratio}{A list of length `g - 1` containing the Adverse
#'             Impact (AI) ratio computed for each focal group. A result less
#'             than 80% may be considered evidence of adverse impact.}
#'      If \code{show_mi_result = TRUE}, the returned list will have the
#'      additional elements below:
#'          \item{propsel_mi}{Proportion selected under strict invariance.}
#'          \item{cutpt_xi_mi}{Cut point on the latent scale (xi) under strict
#'             invariance.}
#'          \item{cutpt_z_mi}{Cut point on the observed scale (Z) under strict
#'             invariance.}
#'         \item{summary_mi}{A `8 x (g + g - 1)` table, with columns
#'             representing the reference and `g - 1` focal groups and the
#'             expected results if the latent distributions of `g - 1` focal
#'             groups match the reference group, under strict invariance. The
#'             rows represent probabilities of true positive (A), false positive
#'             (B), true negative (C), false negative (D); proportion selected,
#'             success ratio, sensitivity, and specificity.}
#'          \item{bivardata_mi}{List of length `5` containing `1 x g` vectors of
#'             latent and observed means, standard deviations, and covariances
#'             computed for each group under strict invariance.}
#'          \item{labels}{List of labels.}
#'          \item{functioncall}{Function call to PartInv.}   
#' @examples
#' set.seed(7)  
#' cols <- c("salmon1", "lightgreen", "skyblue1", "pink")
#' sim_m <-
#'   "f =~ c(1, .7, 1) * x1 + c(.2, 1.1, 1) * x2 + 1 * x3 + 1 * x4 + 1 * x5
#'    f ~~ c(1, 1.3, 1.5) * f
#'    f ~  c(0, .5, 1) * 1
#'    x1 ~ c(0, .3, 0) * 1
#'    x3 ~ c(.3, 0, -.3) * 1
#'    x1 ~~ c(1, .5, 1) * x1"
#' dat_sim <- lavaan::simulateData(sim_m, sample.nobs = c(80, 100, 110))
#' fit_sim <- lavaan::cfa(model = sim_m, data = dat_sim, group = "group")
#' PartInv(cfa_fit = fit_sim, propsel = .05, plot_contour = TRUE, 
#'         custom_colors = cols[1:3], show_mi_result = TRUE)
#'
#' library(lavaan)
#' HS <- HolzingerSwineford1939
#' HS$sex <- as.factor(HS$sex)
#' HS.model <- ' visual  =~ x1 + x2 + x3
#'               textual =~ x4 + x5 + x6
#'               speed   =~ x7 + x8 + x9 '
#' fit <- cfa(HS.model, data = HS, group = "sex")
#' PartInv(fit, propsel = .7, plot_contour = TRUE, show_mi_result = TRUE,
#'         labels = c("M", "F"))
#' # Two groups, single dimension
#' PartInv(cfa_fit = NULL, propsel = .30, weights_item = rep(1, 4),
#'         alpha = list(0, 0), psi = list(1, 1), 
#'         lambda = list(rep(1, 4), c(1, 1, .7, .4)),
#'         nu = list(rep(1, 4), rep(1, 4)), theta = list(diag(1, 4), diag(1, 4)),
#'         labels = c("Female", "Male"), show_mi_result = TRUE, 
#'         plot_contour = TRUE)
#' # Two groups, two dimensions
#' l_mat <- matrix(0, nrow = 5, ncol = 2)
#' l_mat[1:2, 1] <- c(.3, .7); l_mat[3:5, 2] <- c(.4, .8, .5)
#' PartInv(propsel = .05, weights_latent = c(.5, .5),
#'         alpha = list(c(0, 0), c(-.3, .1)),
#'         psi = list(matrix(c(1, .5, .5, 1), nrow = 2),
#'                    matrix(c(1, .5, .5, 1), nrow = 2)),
#'         lambda = list(l_mat, l_mat),
#'         nu = list(c(.3, .3, .01, .3, .1), c(.3, -.05, .3, -.03, .1)),
#'         theta = list(diag(1, 5), c(1, .9, .8, .8, 1)),
#'         plot_contour = TRUE, show_mi_result = TRUE)
#' # Multiple groups, multiple dimensions
#' l_mat <- matrix(c(.7, .8, -.4, .7, .6, .5, .8, -.3, .6, .6, .6, .7, .6, 
#'                   .6, .6), nrow = 15, ncol = 1)
#' nu_mat <- matrix(c(3.6, 3.1, 2.7, 2.9, 2.5, 2.1, 3.5, 2.6, 3.2, 2.8,
#'                    3.5, 3.3, 2.5, 3.4, 2.5), nrow = 15, ncol = 1)
#' nu_mat1 <- nu_mat2 <- nu_mat3 <- nu_mat
#' nu_mat1[c(1, 14:15), 1] <- c(3.9, 3.8, 2.8)
#' nu_mat2[c(9:10, 12, 15), 1] <- c(3.6, 3.2, 3.5, 2.8); nu_mat3[15] <- 2.9
#' th_mat <- (c(.6, .6, .8, .6, .8, .5, .7, 1.1, .8, .9, .6, .7, .8, .5, .9))
#' th_mat1 <- th_mat2 <- th_mat3 <- th_mat
#' th_mat1[1] <- .4; th_mat2[c(6, 14)] <- c(.9, .7); th_mat3[14] <- .7
#' PartInv(propsel = .25, pmix = rep(1/4, 4),
#'         alpha = list(0, -.7, -1.1, -1.1), psi = list(1, 1.2, 1.3, 1.3),
#'         nu = list(nu_mat, nu_mat1, nu_mat2, nu_mat3),
#'         lambda = list(l_mat, l_mat, l_mat, l_mat),
#'         theta = list(th_mat, th_mat1, th_mat2, th_mat3),
#'         plot_contour = TRUE, show_mi_result = TRUE,
#'         labels = c("G1", "G2", "G3", "G4"),
#'         custom_colors = cols, reference = "G2")
#' @export
PartInv <- function(cfa_fit = NULL,
                    propsel = NULL, cut_z = NULL,
                    weights_item = NULL, weights_latent = NULL,
                    alpha = NULL, psi = NULL, lambda = NULL, theta = NULL, nu = NULL,
                    pmix = NULL,
                    pmix_ref = NULL, plot_contour = FALSE,
                    show_mi_result = FALSE,
                    labels = NULL,
                    custom_colors = NULL,
                    reference = NULL,
                    kappa_r = NULL, kappa_f = kappa_r, alpha_r = NULL, alpha_f = alpha_r, phi_r = NULL, phi_f = phi_r, psi_r = NULL, psi_f = psi_r, lambda_r = NULL, lambda_f = lambda_r, tau_r = NULL, tau_f = tau_r, nu_r = NULL, nu_f = nu_r, Theta_r = NULL, Theta_f = Theta_r,
                    ...) {

  functioncall <- match.call()
  # make adjustments for formatting and backward compatibility
  pl <- prep_params(
    cfa_fit, propsel, cut_z, weights_item, weights_latent, alpha, psi, lambda, 
    theta, nu, pmix, pmix_ref, plot_contour, labels, n_dim = NULL, 
    n_i_per_dim = NULL, delete_items = NULL, delete_one_cutoff = NULL, alpha_r, 
    alpha_f, phi_r, phi_f, psi_r, psi_f, lambda_r, lambda_f, tau_r, tau_f, 
    kappa_r, kappa_f, nu_r, nu_f, Theta_r, Theta_f, reference, custom_colors,
    PartInv_fit = NULL)
  
  alpha <- pl$alpha
  psi <- pl$psi
  lambda <- pl$lambda
  nu <- pl$nu
  theta <- pl$theta
  pmix <- pl$pmix
  n_i <- pl$n_i
  num_g <- pl$num_g
  d <- pl$d
  weights_latent <- pl$weights_latent
  weights_item <- pl$weights_item
  labels <- pl$labels

  params <- list(weights_item = weights_item, weights_latent = weights_latent, 
                 alpha = alpha, psi = psi, lambda = lambda, nu = nu, 
                 theta = theta, pmix = pmix, propsel = propsel, labels = labels, 
                 cut_z = cut_z, num_g = num_g)
  out <- do.call(compute_cai, params)

  if (out$propsel <= 0.01) warning("Proportion selected is 1% or less.")

  ai_ratio <- as.data.frame(out$summary[5, (num_g + 1):(num_g + num_g - 1)] /
    out$summary[5, 1])
  names(ai_ratio) <- labels[-1]
  row.names(ai_ratio) <- c("")
  out[["ai_ratio"]] <- ai_ratio

  if (show_mi_result) {
    lambda_avg <- .weighted_average_list(lambda, weights = pmix)
    nu_avg <- .weighted_average_list(nu, weights = pmix)
    theta_avg <- .weighted_average_list(theta, weights = pmix)

    params[["lambda"]] <- replicate(num_g, lambda_avg, simplify = FALSE)
    params[["nu"]] <- replicate(num_g, nu_avg, simplify = FALSE)
    params[["theta"]] <- replicate(num_g, theta_avg, simplify = FALSE)
    
    out_mi <- do.call(compute_cai, params)   
    colnames(out_mi$summary) <- labels
    names(out_mi) <- paste0(names(out_mi), "_mi")
    
    # calculate AI ratio for strict invariance   
    ai_ratio_mi <- as.data.frame(
      out_mi$summary[5, (num_g + 1):(num_g + num_g - 1)] / out_mi$summary[5, 1])
    names(ai_ratio_mi) <- labels[-1]
    row.names(ai_ratio_mi) <- c("")
    
    # remove the E_R(Focal) columns from summary_mi
    out_mi$summary_mi <- out_mi$summary_mi[, 1:num_g]
    out <- c(out, out_mi)
    out[["ai_ratio_mi"]] <- ai_ratio_mi
  }

  if (plot_contour) {
    plot.PartInv(out, labels = labels, which_result = "pi",
                 custom_colors = custom_colors, ...)
    if (show_mi_result == TRUE) {
      plot.PartInv(out, labels = labels, which_result = "mi", 
                   custom_colors = custom_colors, ...)
    }
  }
  out[["labels"]] <- labels
  out[["functioncall"]] <- functioncall
  class(out) <- "PartInv"
  out
}

#' @rdname PartInv
#' @export
PartInvMulti_we <- function(...)
{
  .Deprecated("PartInv")
  # PartInv(...)
}

