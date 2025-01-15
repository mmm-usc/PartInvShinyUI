

#' Called in item_deletion_h.R and PartInv.R to ensure backward
#' compatibility with different input names, check dimension and object 
#' types, and prepare/format the variables. 
#'
#' \code{prep_params} 
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
#'     #' @param n_dim Number of dimensions, 1 by default. If the user does not supply
#'        a different value, proceeds with the assumption that the scale is
#'        unidimensional.
#' @param n_i_per_dim A vector containing the number of items in each
#'        dimension; `NULL` by default. If the user provides a value for `n_dim`
#'        that is \eqn{> 1} but leaves \code{n_i_per_dim = NULL}, assumes that
#'        the subscales have an equal number of items.
#' @param delete_items A vector; default to `NULL`. If the user does not
#'        input a vector of items, only the items determined to contain bias will
#'        be considered for deletion.
#' @param delete_one_cutoff (optional) User-specified cutoff to use in
#'        delete-one scenarios. `NULL` by default; if `NULL`, proportion
#'        selected under SFI and PFI when the full item set is used is passed
#'        onto calls to PartInv.
#' @param alpha_r,alpha_f,nu_r,nu_f,Theta_r,Theta_f,psi_r,psi_f,lambda_r,lambda_f,phi_r,phi_f,tau_r,tau_f,kappa_r,kappa_f,pmix_ref
#'        Deprecated; included only for backward compatibility. When comparing two
#'        groups, parameters with the '_r' suffix refer to the reference group while
#'        parameters with the '_f' suffix refer to the focal group.
prep_params <- function(cfa_fit, propsel, cut_z, weights_item, weights_latent,
     alpha, psi, lambda, theta, nu, pmix, pmix_ref, plot_contour,
     labels, n_dim, n_i_per_dim, delete_items, delete_one_cutoff, 
     alpha_r, alpha_f, phi_r, phi_f, psi_r, psi_f, lambda_r, lambda_f, tau_r, tau_f, 
     nu_r, nu_f, Theta_r, Theta_f, reference, custom_colors) {
 
  #### 'alpha', 'nu', 'lambda', 'theta', and 'psi' #### 
  if (is.null(cfa_fit)) {
    if (is.null(nu)) {
      nu <- vector(2, mode = "list")
      if(!is.null(nu_r)) {
        nu[[1]] <- nu_r; nu[[2]] <- nu_f
        }
      if (!is.null(tau_r)) {
        nu[[1]] <- tau_r; nu[[2]] <- tau_f
      }
    } # else, nu <- nu
    if (is.null(alpha)) { # '|| is.logical(alpha)) && !is.null(alpha_r)'
      alpha <- vector(2, mode = "list")
      if (!is.null(alpha_r)) {
      alpha[[1]] <- as.numeric(alpha_r); alpha[[2]] <- as.numeric(alpha_f)
      }
    } # else, alpha <- alpha
    if (is.null(psi)) { # '|| is.logical(psi)'
      psi <- vector(2, mode = "list")
      if(!is.null(phi_r)) {
        psi[[1]] <- phi_r; psi[[2]] <- phi_f
      }
      if (!is.null(psi_r)) {
        psi[[1]] <- as.numeric(psi_r); psi[[2]] <- as.numeric(psi_f)
      }
    } # else, psi <- psi
    if (is.null(lambda)) {
      lambda <- vector(2, mode = "list")
      if (!is.null(lambda_r)) {
        lambda[[1]] <- lambda_r; lambda[[2]] <- lambda_f
      }
    } # else, lambda <- lambda
    if (is.null(theta)) {
      theta <- vector(2, mode = "list")
      if (!is.null(Theta_r)) {
        theta[[1]] <- Theta_r; theta[[2]] <- Theta_f
      }
    } # else, theta <- theta
  } else { # if cfa_fit is not null
    if (!all(unlist(lapply(list(nu, alpha, psi, lambda, theta), is.null)))) {
      message("Both cfa_fit and estimates were provided. Defaulting to cfa_fit.")
    }
    # extract the parameter estimates from the cfa fit object
    lav_cfa <- unnest_list(lavInspect(cfa_fit, "est"))
    alpha <- lapply(lav_cfa$alpha, FUN = c)
    nu <- lapply(lav_cfa$nu, FUN = c)
    theta <- lav_cfa$theta
    lambda <- lav_cfa$lambda
    psi <- lav_cfa$psi
  }
    
  num_g <- length(alpha)
  n_i <- length(nu[[1]])
  d <- length(alpha[[1]])
  
  alpha <- lapply(alpha, FUN = as.matrix)
  psi <- lapply(psi, FUN = matrix, nrow = d, ncol = d)
  
  stopifnot("theta, nu, and lambda must be lists." =
              (all(is.list(theta) & is.list(nu) & is.list(lambda))))
  # change any vector elements within the list theta into diagonal matrices
  theta <- lapply(seq_along(theta), function(x) {
    if (is.vector(theta[[x]])) {
      theta[[x]] <- diag(theta[[x]])
    } else {
      theta[[x]] <- theta[[x]] # necessary to ensure conformable arguments later
    }
  })
  stopifnot("Number of groups as indicated in the estimates must match." = 
              length(alpha) == lengths(list(psi, lambda, nu, theta)))
  stopifnot("Number of dimensions must match." =
      (((lengths(alpha) == dim(psi)[1]) & (dim(psi)[1] == dim(psi)[2]) &
          lengths(alpha) == unlist(lapply(lambda, ncol)))))
  
  # if (length(alpha) == 1 & length(psi) == 1) {
  #   stop("Check whether alpha and psi have the correct dimensions.")
  # }
  
  #### 'weights_item' and 'weights_latent' ####
  if (is.null(weights_item)) {
    weights_item <- rep(1, n_i)
  } else {
    stopifnot("Please provide a weights_item vector of the correct length." = 
                length(weights_item) == n_i)
  }
  if (is.null(weights_latent)) {
    weights_latent <- rep(1, d)
  } else {
    stopifnot("Please provide a weights_latent vector of the correct length." = 
                length(weights_latent) == d)
  }
  
  #### 'pmix' ####
  if (is.null(pmix)) {
    if (!is.null(pmix_ref)) {
      if (num_g == 2) {
        pmix <- c(pmix_ref, 1 - pmix_ref) # two groups
        }
      } 
    if (!is.null(cfa_fit)) {
      cfa_sum <- summary(cfa_fit)
      pmix <- cfa_sum$data$nobs / sum(cfa_sum$data$nobs)
      }
    if (is.null(pmix_ref) && is.null(cfa_fit)) {
      message("Mixing proportions not provided (pmix). Assuming equal weights.")
      pmix <- rep(1, num_g) / num_g
      }
    } else { # if pmix not null
    stopifnot("Provide the correct number of mixing proportions." = 
                length(pmix) == num_g)
    }
  #### 'labels' #### 
  if (!is.null(labels)) { # 'labels' was provided
    if (length(labels) != num_g) {
      stop("The number of labels does not match the number of groups. Using defaults.")
      labels <- c("Reference", paste0("Focal_", 1:(num_g - 1)))
      lab_text <- "default"
    }
    lab_text <- "provided"
  } else {  # 'labels' is null
    if (!is.null(cfa_fit)) { # user supplied cfa_fit
      labels <- summary(cfa_fit)$data$group.label
      lab_text <- "cfa fit object"
    } else { # user did not supply cfa_fit
      labels <- c("Reference", paste0("Focal_", 1:(num_g - 1)))
      lab_text <- "default (Reference, Focal_1, Focal_2, ...)"
    }
  }
  
  #### set the reference group and reorder appropriately ####
  if (!is.null(reference)) {
    msg <- paste0("The reference group label string does not match any of the ", lab_text, " group labels.")
    if (!reference %in% labels) stop(msg)
    ind <- which(labels == reference)
    labels <- c(labels[[ind]], labels[-ind])
    alpha <- c(list(alpha[[ind]]), alpha[-ind])
    nu <- c(list(nu[[ind]]), nu[-ind])
    theta <- c(list(theta[[ind]]), theta[-ind])
    lambda <- c(list(lambda[[ind]]), lambda[-ind])
    psi <- c(list(psi[[ind]]), psi[-ind])
    pmix <- c(pmix[[ind]], pmix[-ind])
    custom_colors <- c(custom_colors[[ind]], custom_colors[-ind])
  }
  
  g_labs <- c("r", paste0("f", 1:(num_g - 1)))
  names(alpha) <- paste("alpha", g_labs, sep = "_")
  names(nu) <- paste("nu", g_labs, sep = "_")
  names(lambda) <- paste("lambda", g_labs, sep = "_")
  names(psi) <- paste("psi", g_labs, sep = "_")
  names(theta) <- paste("theta", g_labs, sep = "_")

  return(list("alpha" = alpha, "lambda" = lambda, "nu" = nu, "psi" = psi, 
              "theta" = theta, "pmix" = pmix, "num_g" = num_g, "n_i" = n_i, 
              "d" = d, "weights_item" = weights_item,
              "weights_latent" = weights_latent, "labels" = labels))
}
