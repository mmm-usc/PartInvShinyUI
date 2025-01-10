# This function is called in item_deletion_h.R and PartInv.R to ensure backward
# compatibility with different input names, check dimension and object 
# types, and prepare/format the variables. 
prep_params <- function(cfa_fit, propsel, cut_z, weights_item, weights_latent,
     alpha, psi, lambda, theta, nu, pmix, pmix_ref, plot_contour, show_mi_result,
     labels, n_dim, n_i_per_dim, user_specified_items, delete_one_cutoff, 
     alpha_r, alpha_f, phi_r, psi_r, psi_f, lambda_r, lambda_f, tau_r, tau_f, 
     nu_r, nu_f, Theta_r, Theta_f, func_called) {
  
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
    if ((is.null(alpha) || is.logical(alpha)) && !is.null(alpha_r)) {
      alpha <- vector(2, mode = "list")
      alpha[[1]] <- as.numeric(alpha_r); alpha[[2]] <- as.numeric(alpha_f)
    }
    if (is.null(psi) || is.logical(psi)) {
      psi <- vector(2, mode = "list")
      if(!is.null(phi_r)) {
        psi[[1]] <- phi_r; psi[[2]] <- phi_f
      }
      if (!is.null(psi_r)) {
        psi[[1]] <- as.numeric(psi_r); psi[[2]] <- as.numeric(psi_f)
      }
    } # else, psi <- psi

    if (is.null(lambda) && !is.null(lambda_r)) {
      lambda <- vector(2, mode = "list")
      lambda[[1]] <- lambda_r; lambda[[2]] <- lambda_f
    } # else, lambda <- lambda
    if (is.null(theta) && !is.null(Theta_r)) {
      theta <- vector(2, mode = "list")
      theta[[1]] <- Theta_r; theta[[2]] <- Theta_f
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
  
  if (is.null(pmix)) {
    if (!is.null(pmix_ref)) {
      pmix <- c(pmix_ref, 1 - pmix_ref) # assuming two groups
    }
    if (!is.null(cfa_fit)) {
      cfa_sum <- summary(cfa_fit)
      pmix <- cfa_sum$data$nobs / sum(cfa_sum$data$nobs)
    }
    if (is.null(cfa_fit) & is.null(pmix_ref)) {
      pmix <- rep(1, length(alpha)) / length(alpha)
      warning("Mixing proportions not provided (pmix). Assuming equal weights.")
    }
  }
  stopifnot("Provide the correct number of mixing proportions." = 
              length(pmix) == length(alpha))
  
  # Change any vector elements within the list theta into diagonal matrices
  theta <- lapply(seq_along(theta), function(x) {
    if (is.vector(theta[[x]])) {
      theta[[x]] <- diag(theta[[x]])
    } else {
      theta[[x]] <- theta[[x]] # necessary to ensure conformable arguments later
    }
  })
  
  num_g <- length(alpha)
  n_i <- length(nu[[1]])
  d <- length(alpha[[1]])
  
  alpha <- lapply(alpha, as.matrix)
  psi <- lapply(psi, matrix, nrow = d, ncol = d)
  
  stopifnot("theta, nu, and lambda must be lists. Consider using `format_cfa_partinv()`." =
              (all(is.list(theta) & is.list(nu) & is.list(lambda))))
  stopifnot("Number of groups as indicated in the lengths of parameters must 
              match." = length(alpha) == lengths(list(psi, lambda, nu, theta)))
  stopifnot("Number of dimensions must match." =
      (((lengths(alpha) == dim(psi)[1]) & (dim(psi)[1] == dim(psi)[2]) &
          lengths(alpha) == unlist(lapply(lambda, ncol))))
  )
  if(length(alpha) == 1 & length(psi) == 1) {
    stop("Check whether alpha and psi have the correct dimensions.")
  }
  
  if (length(weights_latent) == 1) weights_latent <- rep(1, d)
  
  if (is.null(weights_item)) weights_item <- rep(1, n_i)
  if (!is.null(weights_item) & length(weights_item) != n_i) {
    stop("Please provide a weights_item vector of the correct length.")
  }
  
  if (is.null(weights_latent)) weights_latent <- rep(1, d)
  if (!is.null(weights_latent) & length(weights_latent) != d) {
    stop("Please provide a weights_latent vector of the correct length.")
  }
  
  g <- c("r", paste0("f", 1:(num_g - 1)))
  names(alpha) <- paste("alpha", g, sep = "_")
  names(nu) <- paste("nu", g, sep = "_")
  names(lambda) <- paste("lambda", g, sep = "_")
  names(psi) <- paste("psi", g, sep = "_")
  names(theta) <- paste("theta", g, sep = "_")

  return(list("alpha" = alpha, "lambda" = lambda, "nu" = nu, "psi" = psi, 
              "theta" = theta, "pmix" = pmix, "num_g" = num_g, "n_i" = n_i, 
              "d" = d, "weights_item" = weights_item,
              "weights_latent" = weights_latent))
}