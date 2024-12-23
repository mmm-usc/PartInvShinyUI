#' @importFrom stats qchisq pnorm qnorm nlminb
#' @importFrom mnormt pmnorm
NULL

#' Evaluating selection accuracy for two or more groups based on the MCAA
#' Framework
#'
#' \code{PartInv, PartInvMulti_we} evaluates partial measurement invariance
#' using the multidimensional classification accuracy analysis
#' (Lai & Zhang, 2022), which is an extension of Millsap & Kwok's (2004)
#' approach.
#'
#' @param cfa_fit CFA model output from lavaan.
#' @param propsel Proportion of selection. If missing, computed using `cut_z`.
#' @param cut_z Pre-specified cutoff score on the observed composite. This
#'     argument is ignored when `propsel` has input.
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
#'     and `n` is the number of items in the scale. The first element is assumed
#'     to belong to the reference group.
#' @param nu A list of length `g` containing `1 x n` measurement intercept
#'     vectors where `g` is the number of groups and `n` is the number of items
#'     in the scale. The first element is assumed to belong to the reference
#'     group.
#' @param theta A list of length `g` containing `1 x n` vectors or `n x n`
#'     matrices of unique factor variances and covariances, where `g` is the
#'     number of groups and `n` is the number of items in the scale. The first
#'     element is assumed to belong to the reference group.
#' @param pmix List of length `g` containing the mixing proportions of each
#'     group. If `NULL`, defaults to `1/g` for each group (i.e., the populations
#'     have equal size).
#' @param plot_contour Logical; whether the contour of the populations should be
#'     plotted; `TRUE` by default.
#' @param show_mi_result If \code{TRUE}, perform selection accuracy analysis
#'     for both the input parameters and the implied parameters based on a
#'     strict invariance model, with common parameter values as weighted
#'     averages of the input values using `pmix`.
#' @param labels A character vector with `g` elements to label the reference
#'     and focal groups on the plot, where `g` is the number of groups. If not
#'     provided, groups are labeled automatically as 'Reference' (for the first
#'     group) and 'Focal_1' through 'Focal_(g-1)', where `g` is the number of
#'     groups.
#' @param custom_colors Optional argument for specifying colors
#' @param ... Other arguments passed to the \code{\link[graphics]{contour}}
#'     function.
#' @param alpha_r,alpha_f,nu_r,nu_f,Theta_r,Theta_f,psi_r,psi_f,lambda_r,lambda_f,phi_r,phi_f,tau_r,tau_f,kappa_r,kappa_f,pmix_ref
#'     Deprecated; included only for backward compatibility. When comparing two
#'     groups, parameters with the '_r' suffix refer to the reference group while
#'     parameters with the '_f' suffix refer to the focal group.
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
#' # Two groups, using cfa fit object
#' library(lavaan)
#' HS <- HolzingerSwineford1939
#' HS$sex <- as.factor(HS$sex)
#' HS.model <- ' visual  =~ x1 + x2 + x3
#'               textual =~ x4 + x5 + x6
#'               speed   =~ x7 + x8 + x9 '
#' 
#' fit <- cfa(HS.model, data = HS, group = "sex")
#' PartInv(fit, propsel = .7, plot_contour = TRUE, show_mi_result = TRUE)
#' # Two groups, single dimension
#' PartInv(propsel = .30,
#'         weights_item = c(1, 1, 1, 1),
#'         weights_latent = 1,
#'         alpha = list(0, 0),
#'         psi = list(1, 1),
#'         lambda = list(c(1, 1, 1, 1), c(1, 1, 1, 1)),
#'         nu = list(c(1, 1, 1, 2), c(1, 1, 1, 1)),
#'         theta = list(diag(1, 4), diag(1, 4)),
#'         labels = c("Female", "Male"),
#'         show_mi_result = TRUE)
#' # Two groups, two dimensions
#' lambda_matrix <- matrix(0, nrow = 5, ncol = 2)
#' lambda_matrix[1:2, 1] <- c(.322, .655)
#' lambda_matrix[3:5, 2] <- c(.398, .745, .543)
#' PartInv(propsel = .05,
#'         weights_latent = c(0.5, 0.5),
#'         alpha = list(c(0, 0), c(-0.3, 0.1)),
#'         psi = list(matrix(c(1, 0.5, 0.5, 1), nrow = 2),
#'                    matrix(c(1, 0.5, 0.5, 1), nrow = 2)),
#'         lambda = list(lambda_matrix, lambda_matrix),
#'         nu = list(c(.225, .25, .010, .30, .125),
#'                   c(.225, -.05, .240, -.025, .125)),
#'         theta = list(diag(1, 5), c(1, .95, .80, .75, 1)),
#'         plot_contour = TRUE, show_mi_result = TRUE)
#' # Multiple groups, multiple dimensions
#' lambda_matrix <- matrix(0, nrow = 15, ncol = 1)
#' lambda_matrix[1:15, 1] <- c(0.68, 0.79, -0.39, 0.74, 0.59, 0.46, 0.78, -0.30,
#'                             0.59, 0.59, 0.64, 0.66, 0.59, 0.63, 0.64);
#' nu_matrix <- nu_matrix1 <- nu_matrix2 <- nu_matrix3 <-
#'   matrix(0, nrow = 15, ncol = 1)
#' nu_matrix[1:15, 1] <- c(3.6, 3.1, 2.7, 2.9, 2.5, 2.1, 3.45, 2.62, 3.2, 2.84,
#'                         3.51, 3.26, 2.45, 3.39, 2.47);
#' nu_matrix1[1:15, 1] <- c(3.9, 3.1, 2.7, 2.9, 2.5, 2.1, 3.45, 2.62, 3.2, 2.84,
#'                          3.51, 3.26, 2.45, 3.76, 2.81);
#' nu_matrix2[1:15, 1] <- c(3.6, 3.1, 2.7, 2.9, 2.5, 2.1, 3.45, 2.62, 3.6, 3.18,
#'                          3.51, 3.54, 2.45, 3.39, 2.81);
#' nu_matrix3[1:15, 1] <- c(3.6, 3.1, 2.7, 2.6, 2.5, 2.1, 3.45, 2.62, 3.2, 2.84,
#'                          3.51, 3.26, 2.45, 3.39, 2.81);
#' theta_matrix <- c(0.35, 0.62, 0.83, 0.61, 0.81, 0.87, 0.39, 1.05, 0.84, 0.92,
#'                   0.36, 0.66, 0.8, 0.66, 0.9);
#' theta_matrix1 <- c(0.61, 0.62, 0.83, 0.61, 0.81, 0.5, 0.7, 1.05, 0.84, 0.92,
#'                    0.61, 0.66, 0.8, 0.54, 0.9);
#' theta_matrix2 <- c(0.61, 0.62, 0.826, 0.61, 0.81, 0.87, 0.5, 1.05, 0.84,
#'                    0.92, 0.61, 0.66, 0.8, 0.66, 0.9);
#' theta_matrix3 <- c(0.61, 0.62, 0.826, 0.61, 0.81, 0.5, 0.7, 1.05, 0.84, 0.92,
#'                    0.61, 0.66, 0.8, 0.66, 0.9);
#' PartInv(propsel = 0.25, pmix = c(1/4, 1/4, 1/4, 1/4),
#'         alpha = list(0, -0.70, -1.05, -1.10), psi = list(1, 1.2, 1.29, 1.3),
#'         nu = list(nu_matrix, nu_matrix1, nu_matrix2, nu_matrix3),
#'         lambda = list(lambda_matrix, lambda_matrix, lambda_matrix,
#'                       lambda_matrix),
#'         theta = list(theta_matrix, theta_matrix1, theta_matrix2,
#'                      theta_matrix3),
#'         plot_contour = TRUE, show_mi_result = TRUE,
#'         labels = c("Group 1", "Group 2", "Group 3", "Group 4"),
#'         custom_colors = c("salmon1", "lightgreen", "skyblue1", "pink")
#'         )
#' @export
PartInv <- function(cfa_fit = NULL, propsel = NULL, cut_z = NULL,
                    weights_item = NULL,
                    weights_latent = NULL,
                    alpha = NULL, psi = NULL, lambda = NULL, theta = NULL, nu = NULL,
                    pmix = NULL,
                    pmix_ref = NULL, plot_contour = FALSE,
                    show_mi_result = FALSE,
                    labels = NULL,
                    custom_colors = NULL,
                    kappa_r = NULL, kappa_f = kappa_r,
                    alpha_r = NULL, alpha_f = alpha_r,
                    phi_r = NULL, phi_f = phi_r,
                    psi_r = NULL, psi_f = psi_r,
                    lambda_r = NULL, lambda_f = lambda_r,
                    tau_r = NULL, tau_f = tau_r,
                    nu_r = NULL, nu_f = nu_r,
                    Theta_r = NULL, Theta_f = Theta_r,
                    ...) {

  functioncall <- match.call()
  # for backward compatibility with different input names ####
  if (is.null(cfa_fit) && (is.null(nu) && !is.null(nu_r))) {
    nu <- vector(2, mode = "list")
    nu[[1]] <- nu_r
    nu[[2]] <- nu_f
  }
  if (is.null(cfa_fit) && (is.null(nu) && !is.null(tau_r))) {
    nu <- vector(2, mode = "list")
    nu[[1]] <- tau_r
    nu[[2]] <- tau_f
  }
  if (is.null(cfa_fit) && 
      ((is.null(alpha) || is.logical(alpha)) && !is.null(kappa_r))) {
    alpha <- vector(2, mode = "list")
    alpha[[1]] <- kappa_r
    alpha[[2]] <- kappa_f
  }
  if (is.null(cfa_fit) && 
      ((is.null(alpha) || is.logical(alpha)) && !is.null(alpha_r))) {
    alpha <- vector(2, mode = "list")
    alpha[[1]] <- as.numeric(alpha_r)
    alpha[[2]] <- as.numeric(alpha_f)
  }
  
  if (is.null(cfa_fit) && 
      ((is.null(psi) || is.logical(psi)) && !is.null(phi_r))) {
    psi <- vector(2, mode = "list")
    psi[[1]] <- phi_r
    psi[[2]] <- phi_f
  }
  if (is.null(cfa_fit) && 
      ((is.null(psi) || is.logical(psi)) && !is.null(psi_r))) {
    psi <- vector(2, mode = "list")
    psi[[1]] <- as.numeric(psi_r)
    psi[[2]] <- as.numeric(psi_f)
  }

  if (is.null(lambda) && !is.null(lambda_r)) {
    lambda <- vector(2, mode = "list")
    lambda[[1]] <- lambda_r
    lambda[[2]] <- lambda_f
  }
  if (is.null(theta) && !is.null(Theta_r)) {
    theta <- vector(2, mode = "list")
    theta[[1]] <- Theta_r
    theta[[2]] <- Theta_f
  }
  if (is.null(pmix) && !is.null(pmix_ref)) {
    pmix <- c(pmix_ref, 1 - pmix_ref) # assuming two groups
  }
  # end ####
  
  if (!is.null(cfa_fit)) {
    
    if (!any(unlist(lapply(list(nu, alpha, psi, lambda, theta), is.null)))) {
      message("Both cfa_fit and parameter estimates were provided. Defaulting
               to using cfa_fit.")
    }
    
    lav_cfa <- unnest_list(lavInspect(cfa_fit, "est"))
    # extract the parameter estimates from the cfa fit object
    alpha <- lapply(lav_cfa$alpha, FUN = c)
    nu <- lapply(lav_cfa$nu, FUN = c)
    theta <- lav_cfa$theta
    lambda <- lav_cfa$lambda
    psi <- lav_cfa$psi
  }
  
  
  stopifnot("theta, nu, and lambda must be lists. Consider using `format_cfa_partinv()`." =
              (all(is.list(theta) & is.list(nu) & is.list(lambda))))
  stopifnot("Number of groups as indicated in the lengths of parameters must 
              match." = length(alpha) == lengths(list(psi, lambda, nu, theta)))
  stopifnot(
    "Number of dimensions must match." =
      (((lengths(alpha) == dim(psi)[1]) & (dim(psi)[1] == dim(psi)[2]) &
          lengths(alpha) == unlist(lapply(lambda, ncol))))
    )

  num_g <- length(alpha)
  n <- length(nu[[1]])
  d <- length(alpha[[1]])

  if (is.null(pmix)) pmix <- as.matrix(c(rep(1 / num_g, num_g)), ncol = num_g)
  pmix <- as.vector(pmix)

  stopifnot(
    "Provide the correct number of mixing proportions." =
      length(pmix) == length(alpha)
  )
  
  
  if (length(weights_latent) == 1) weights_latent <- rep(1, d)
  
  if(length(alpha) == 1 & length(psi) == 1) {
    stop("Check whether alpha and psi have the correct dimensions.")
  }

  # If labels were not provided by the user or the number of labels provided or
  # the number of labels provided does not match num_g, define new labels
  if (is.null(labels) || (length(labels) != num_g)) {
    labels <- c("Reference", paste0("Focal_", 1:(num_g - 1)))
  }

  g <- c("r", paste0("f", 1:(num_g - 1)))
  names(alpha) <- paste("alpha", g, sep = "_")
  names(nu) <- paste("nu", g, sep = "_")
  names(lambda) <- paste("lambda", g, sep = "_")
  names(psi) <- paste("psi", g, sep = "_")
  names(theta) <- paste("theta", g, sep = "_")

  # Change any vector elements within the list theta into diagonal matrices
  theta <- lapply(seq_along(theta), function(x) {
    if (is.vector(theta[[x]])) {
      theta[[x]] <- diag(theta[[x]])
    } else {
      theta[[x]] <- theta[[x]] # necessary to ensure conformable arguments later
    }
  })

  alpha <- lapply(alpha, as.matrix)
  psi <- lapply(psi, matrix, nrow = d, ncol = d)

  if (is.null(weights_item)) weights_item <- rep(1, n)
  if (!is.null(weights_item) & length(weights_item) != n) {
    stop("Please provide a weights_item vector of the correct length.")
  }
  
  if (is.null(weights_latent)) weights_latent <- rep(1, d)
  if (!is.null(weights_latent) & length(weights_latent) != d) {
    stop("Please provide a weights_latent vector of the correct length.")
  }
  
  out <- compute_cai(weights_item, weights_latent, alpha, psi, lambda, nu,
    theta, pmix, propsel, labels, cut_z,
    is_mi = FALSE
  )

  if (out$propsel <= 0.01) warning("Proportion selected is 1% or less.")

  ai_ratio <- as.data.frame(out$summary[5, (num_g + 1):(num_g + num_g - 1)] /
    out$summary[5, 1])

#    message("Note: The first group is being used as the reference group. 
#Rearrange the inputs to designate a different group as the reference.\n\n")
  
  names(ai_ratio) <- labels[-1]
  row.names(ai_ratio) <- c("")
  out[["ai_ratio"]] <- ai_ratio

  if (show_mi_result) {
    pop_weights <- pmix
    lambda_average <- .weighted_average_list(lambda, weights = pop_weights)
    nu_average <- .weighted_average_list(nu, weights = pop_weights)
    theta_average <- .weighted_average_list(theta, weights = pop_weights)

    lambda_average_g <- nu_average_g <- theta_average_g <-
      vector(mode = "list", length = num_g)

    for (i in 1:num_g) {
      lambda_average_g[[i]] <- lambda_average
      nu_average_g[[i]] <- nu_average
      theta_average_g[[i]] <- theta_average
    }

    out_mi <- compute_cai(weights_item, weights_latent, alpha, psi,
      lambda_average_g, nu_average_g, theta_average_g,
      pmix, propsel, labels, cut_z,
      is_mi = TRUE
    )
    colnames(out_mi$summary) <- labels
    names(out_mi) <- paste0(names(out_mi), "_mi")
    out <- c(out, out_mi)
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
