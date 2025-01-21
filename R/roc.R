#' @import lavaan

# helper function to compute AUC using SE and SP
auc <- function(SE, SP) {
  TPR <- SE
  FPR <- 1 - SP
  dFPR <- c(diff(FPR), 0)
  dTPR <- c(diff(TPR), 0)
  sum(TPR * dFPR) + sum(dTPR * dFPR) / 2
}


#' @title
#' Return a list containing sensitivity (SE) and specificity (SP) values
#'
#' @name
#' return_SE_SP
#'
#' @description
#' \code{return_SE_SP} takes in a CFA fit object, output from PartInv(), or 
#' vectors of model parameter estimates and returns a list object containing 
#' SE and SP values under partial and/or strict invariance.
#'
#' @param cfa_fit Output from cfa().
#' @param PartInv_fit Output from PartInv()
#' @param pmix Vector of mixing proportions.
#' @param from 0.01 by default.
#' @param to 0.9999 by default.
#' @param by 0.01 by default.
#' @param cutoffs_from  NULL by default.
#' @param cutoffs_to NULL by default.
#' @param mod_names c("partial", "strict") by default.
#' @param alpha A list of length `g` containing `1 x d` latent factor mean
#'     vectors where `g` is the number of groups and `d` is the number of latent
#'     dimensions. The first element is assumed to belong to the reference group. 
#'     NULL by default.
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
#' @param labels Group labels. NULL as default. 
#' @param ... Other arguments passed to the \code{\link[graphics]{contour}}
#'     function.
#' @return The output will be a list containing the SE and SP values under 
#'     the specified invariance conditions.
#' @examples
#' # when vectors of model parameter estimates are provided as input
#' lam_vec <- c(.68, .8, -.39, .74, .6, .5, .8, -.30, .6, .6, .64,.66, .6, .63)
#' nu_vec <- nu_vec1 <- nu_vec2 <- nu_vec3 <- 
#'  c(3.6, 3, 2.7, 3, 2.5, 2.1, 3.5, 2.6, 3.2, 2.8, 3.5, 3.3, 2.45, 3.4)
#' nu_vec1[1] <- 3.9; nu_vec1[12:14] <- c(3.6, 2.45, 3.7)
#' nu_vec2[8:9] <- c(3.6, 3.2); nu_vec2[12] <- 3.5
#' nu_vec3[4] <- 2.6
#' theta_vec <- theta_vec1 <- theta_vec2 <- theta_vec3 <- 
#'  c(.6, .6, .83, .6, .8, .87, .4, 1, .84, .9, .6, .66, .8, .66)
#' theta_vec1[1] <- .35; theta_vec1[6] <- .5; theta_vec1[11] <- .36
#' theta_vec2[3] <- .8; theta_vec2[3] <- .3
#' theta_vec3[3] <- .8; theta_vec3[6:7] <- c(.5, .7)
#' out <- return_SE_SP(pmix = rep(1/4, 4),
#'                    alpha = list(0, -.70, -1.05, -1.10),
#'                    psi = list(1, 1.2, 1.29, 1.3),
#'                    nu = list(nu_vec, nu_vec1, nu_vec2, nu_vec3),
#'                    lambda = list(lam_vec, lam_vec, lam_vec, lam_vec),
#'                    theta = list(theta_vec, theta_vec1, theta_vec2, theta_vec3),
#'                    plot_contour = FALSE, show_mi_result = TRUE,
#'                    mod_names = c("strict", "partial"))
#' str(out)
#' # when PartInv output is provided as input
#' PI_out <- PartInv(propsel = .30, weights_item = rep(1, 4),
#'           weights_latent = 1, alpha = list(0, 0), psi = list(1, 1),
#'           lambda = list(rep(1, 4), rep(1, 4)),
#'           nu = list(c(1, 1, 1, 2), rep(1, 4)),
#'           theta = list(diag(1, 4), diag(1, 4)),
#'           labels = c("Female", "Male"), show_mi_result = TRUE)
#' out2 <- return_SE_SP(PartInv_fit = PI_out, 
#'                     mod_names = c("strict", "partial"))
#' str(out2)
#' # when cfa output is provided as input
#' library(lavaan)
#' data("HolzingerSwineford1939", package = "lavaan")
#' HS <- HolzingerSwineford1939
#' HS$sex <- as.factor(HS$sex)
#' HS.model <- ' visual  =~ x1 + x2 + x3
#'               textual =~ x4 + x5 + x6
#'               speed   =~ x7 + x8 + x9 '
#' HS.fit <- cfa(HS.model, data = HS, group = "sex")
#' out3 <- return_SE_SP(cfa_fit = HS.fit, 
#'                     mod_names = c("strict", "partial"))
#' str(out3)
#' @export
return_SE_SP <- function(cfa_fit = NULL,
                         PartInv_fit = NULL,
                         pmix = NULL, 
                         from = 0.01,
                         to = 0.9999,
                         by = 0.01,
                         cutoffs_from = NULL,
                         cutoffs_to = NULL,
                         mod_names = c("partial", "strict"), 
                         alpha = NULL,
                         lambda = NULL,
                         theta = NULL,
                         psi = NULL,
                         nu = NULL, 
                         labels = NULL,
                         weights_item = NULL,
                         weights_latent = NULL,
                         ...) {
  # if PartInv_fit was provided, use it to extract relevant parameters.
  if (!is.null(PartInv_fit)) {
    if (is.null(labels)) labels <- PartInv_fit$labels
    psi <- eval(PartInv_fit$functioncall$psi)
    lambda <- eval(PartInv_fit$functioncall$lambda)
    theta <- eval(PartInv_fit$functioncall$theta)
    alpha <- eval(PartInv_fit$functioncall$alpha)
    nu <- eval(PartInv_fit$functioncall$nu)
    weights_latent <- eval(PartInv_fit$functioncall$weights_latent)
    weights_item <- eval(PartInv_fit$functioncall$weights_item)
    num_g <- length(alpha)
    n_i <- length(nu[[1]])
    d <- length(alpha[[1]])
    alpha <- lapply(alpha, FUN = as.matrix)
    psi <- lapply(psi, FUN = matrix, nrow = d, ncol = d)
    
    if (is.null(pmix)) {
      nbs <- eval(PartInv_fit$cfa_fit)
      if (!is.null(nbs)) {
        pmix <- nbs@data@nobs / sum(nbs@data@nobs)
      } else {
        message("Mixing proportions not provided (pmix). Assuming equal weights.")
        pmix <- rep(1, num_g) / num_g
      }
    } 
  } else {
    pl <- prep_params(
      cfa_fit, propsel, cut_z, weights_item, weights_latent, alpha, psi, lambda,
      theta, nu, pmix, pmix_ref = NULL, plot_contour = NULL, labels, n_dim = NULL,
      n_i_per_dim = NULL, delete_items = NULL, delete_one_cutoff = NULL, 
      alpha_r = NULL, alpha_f = NULL, phi_r = NULL, phi_f = NULL, psi_r = NULL, 
      psi_f = NULL, lambda_r = NULL, lambda_f = NULL, tau_r = NULL, tau_f = NULL,
      kappa_r = NULL, kappa_f = NULL, nu_r = NULL, nu_f = NULL, Theta_r = NULL, 
      Theta_f = NULL, reference = NULL, custom_colors = NULL, PartInv_fit)

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
  }
  
  propsels <- seq(from = from, to = to, by = by)
  use <- "propsels" # set to use proportion of selection by default
  rangeVals <- propsels

  # if the user provided the max and min cutoff values, update rangeVals with
  # a range of cutoffs
  if (!is.null(cutoffs_from) && !is.null(cutoffs_to)) {
    cutoffs <- seq(from = cutoffs_from, to = cutoffs_to, by = by)
    rangeVals <- cutoffs
    use <- "cutoffs" # update to use cutoffs instead
  }

  cai <- c("Sensitivity", "Specificity")
  ls_mat <- matrix(NA, ncol = length(rangeVals), nrow = num_g,
                   dimnames = list(labels, rangeVals))
  ls_names <- c(t(outer(cai, Y = mod_names, FUN = paste, sep = "_")))
  vals <- rep(list(ls_mat), length(ls_names))
  names(vals) <- ls_names

  # call PartInv with each propsel/cutoff and store CAI in the list of dfs
  for (p in seq_along(rangeVals)) {
    if (use == "cutoffs") {
      cut_z <- cutoffs[p]; propsel <- NULL
    } else {
      propsel <- propsels[p]; cut_z <- NULL
    }
    params <- list(weights_item = weights_item, weights_latent = weights_latent, 
                   alpha = alpha, psi = psi, lambda = lambda, nu = nu, 
                   theta = theta, pmix = pmix, propsel = propsel, 
                   labels = labels, cut_z = cut_z, num_g = num_g)
    pinv <- do.call(compute_cai, params)
    
    if ("strict" %in% mod_names) {
      lambda_avg <- .weighted_average_list(lambda, weights = pmix)
      nu_avg <- .weighted_average_list(nu, weights = pmix)
      theta_avg <- .weighted_average_list(theta, weights = pmix)
      
      params[["lambda"]] <- replicate(num_g, lambda_avg, simplify = FALSE)
      params[["nu"]] <- replicate(num_g, nu_avg, simplify = FALSE)
      params[["theta"]] <- replicate(num_g, theta_avg, simplify = FALSE)
      
      out_mi <- do.call(compute_cai, params)
      colnames(out_mi$summary) <- labels
      names(out_mi) <- paste0(names(out_mi), "_mi")
      out_mi$summary_mi <- out_mi$summary_mi[, 1:num_g]
      pinv <- c(pinv, out_mi)
    }
    
    num_comb <- length(cai) * length(mod_names) + 1 # for specifying the index
    ind <- 1

    while (ind < num_comb) {
      for (i in cai) { 
        for (j in seq_along(mod_names)) {
          vals[[ind]][, p] <-
            ifelse(rep(mod_names[j] == "partial", num_g),
                   as.numeric(pinv$summary[i, 1:num_g]),
                   as.numeric(pinv$summary_mi[i, 1:num_g]))
          ind <- ind + 1
        }
      }
    }
  }
  return(vals)
}


#' @title
#' Plot the receiver-operating characteristic (ROC) curve and compute AUC
#'
#' @name
#' roc_auc_PartInv
#'
#' @description
#' \code{roc_auc_PartInv} takes in the output from `return_SE_SP` and plots ROC
#' curves and computes AUCs for the specified groups under the specified
#' invariance conditions.
#'
#' @param cfa_fit Output from cfa(). NULL by default.
#' @param PartInv_fit Output from PartInv(). NULL by default.
#' @param return_SE_SP_output output of `return_SE_SP()`.
#' @param plot_mods vector of strings indicating invariance conditions.
#' @param plot_group vector of integers indicating the groups of interest.
#'        `NULL` by default
#' @param ... Other arguments passed to the \code{\link[graphics]{contour}}
#'     function.
#' @return A list of of length plot_mods
#'
#' @examples
#'# using the output from return_SE_SP() as input
#'lam_vec <- c(.68, .8, -.39, .74, .6, .5, .8, -.30, .6, .6, .64,.66, .6, .63)
#'nu_vec <- nu_vec1 <- nu_vec2 <- nu_vec3 <- 
#'  c(3.6, 3, 2.7, 3, 2.5, 2.1, 3.5, 2.6, 3.2, 2.8, 3.5, 3.3, 2.45, 3.4)
#'nu_vec1[1] <- 3.9; nu_vec1[12:14] <- c(3.6, 2.45, 3.7)
#'nu_vec2[8:9] <- c(3.6, 3.2); nu_vec2[12] <- 3.5
#'nu_vec3[4] <- 2.6
#'theta_vec <- theta_vec1 <- theta_vec2 <- theta_vec3 <- 
#'  c(.6, .6, .83, .6, .8, .87, .4, 1, .84, .9, .6, .66, .8, .66)
#'theta_vec1[1] <- .35; theta_vec1[6] <- .5; theta_vec1[11] <- .36
#'theta_vec2[3] <- .8; theta_vec2[3] <- .3
#'theta_vec3[3] <- .8; theta_vec3[6:7] <- c(.5, .7)
#'out <- return_SE_SP(pmix = rep(1/4, 4),
#'                    alpha = list(0, -.70, -1.05, -1.10),
#'                    psi = list(1, 1.2, 1.29, 1.3),
#'                    nu = list(nu_vec, nu_vec1, nu_vec2, nu_vec3),
#'                    lambda = list(lam_vec, lam_vec, lam_vec, lam_vec),
#'                    theta = list(theta_vec, theta_vec1, theta_vec2, theta_vec3),
#'                    plot_contour = FALSE, show_mi_result = TRUE,
#'                    mod_names = c("strict", "partial"))
#'roc_auc_PartInv(return_SE_SP_output = out, plot_mods = c("strict", "partial"))
#'roc_auc_PartInv(return_SE_SP_output = out, plot_mods = c("partial"))
#'
#'# when PartInv output is provided as input
#' PI_out <- PartInv(propsel = .30, weights_item = rep(1, 4),
#'                  weights_latent = 1, alpha = list(0, 0), psi = list(1, 1),
#'                  lambda = list(rep(1, 4), rep(1, 4)),
#'                  nu = list(c(1, 1, 1, 2), rep(1, 4)),
#'                  theta = list(diag(1, 4), diag(1, 4)),
#'                  labels = c("Female", "Male"), show_mi_result = TRUE)
#'roc_auc_PartInv(PartInv_fit = PI_out, plot_mods = c("partial"))
#'roc_auc_PartInv(PartInv_fit = PI_out, plot_mods = c("partial", "strict"))
#'
#'# when cfa output is provided as input
#'library(lavaan)
#'data("HolzingerSwineford1939", package = "lavaan")
#'HS <- HolzingerSwineford1939
#'HS$sex <- as.factor(HS$sex)
#'HS.model <- ' visual  =~ x1 + x2 + x3
#'              textual =~ x4 + x5 + x6
#'              speed   =~ x7 + x8 + x9 '
#'HS.fit <- cfa(HS.model, data = HS, group = "sex")
#'roc_auc_PartInv(cfa_fit = HS.fit, plot_mods = c("partial"))
#'roc_auc_PartInv(cfa_fit = HS.fit, plot_mods = c("partial", "strict"))

#' @export
roc_auc_PartInv <- function(cfa_fit = NULL, PartInv_fit = NULL,
                            return_SE_SP_output = NULL, 
                            plot_mods = c("partial", "strict"),
                            plot_group = NULL, ...) {
  out <- ""
  if(!is.null(cfa_fit) & is.null(return_SE_SP_output) & is.null(PartInv_fit)) {
    out <- return_SE_SP(cfa_fit, mod_names = plot_mods)
  }
  if(is.null(cfa_fit) & is.null(return_SE_SP_output) & !is.null(PartInv_fit)) {
    # since PartInv needs to be called repeatedly, we extract the relevant 
    # parameters from the cfa fit
    call_cfa <- eval(PartInv_fit$functioncall)
    out <- return_SE_SP(PartInv_fit = call_cfa, mod_names = plot_mods)
  }
  if(is.null(cfa_fit) & is.null(PartInv_fit) & !is.null(return_SE_SP_output)) {
    out <- return_SE_SP_output
  }
  if(sum(c(!is.null(cfa_fit),!is.null(return_SE_SP_output), 
           !is.null(PartInv_fit))) > 1) {
    out <- return_SE_SP_output
    warning('Please only provide one of PartInv_fit, cfa_fit and return_SE_SP_output.')
  }
 
  if(is.null(plot_group)) {
    plot_group <- c(1:length(rownames(out[1][[1]])))
  }

  #out <- return_SE_SP_output
  
  labx <- "FPR (1-SP)"
  laby <- "TPR (SE)"
  SEs <- SPs <- data.frame()
  AUCs <- c()

  par_ind <- grep("par", names(out))
  str_ind <- grep("str", names(out))
  # if partial invariance is indicated and `out` does contain partial inv. results
  if ("partial" %in% plot_mods && (!(is.integer(par_ind) && length(par_ind) == 0L))) {
    partial <- c()
    # for each group, extract SE and SP values, pad with 0 and 1, plot ROC
    for (g in plot_group) {
      SEs <- out[par_ind[1]][[1]][g,]
      SPs <- out[par_ind[2]][[1]][g,]
      labg <- ifelse(g == 1, "Reference", paste0("Group ", g))
      plot(y = c(0, SEs, 1), x = c(0, 1 - SPs,1), type = "l",
           ylim = c(0, 1), xlim = c(0, 1), xlab = labx, ylab = laby,
           main = paste0("ROC: ", labg," (Partial Invariance)"))
      # compute AUC for the group
      partial <- cbind(partial, auc(SEs, SPs))
    }
    colnames(partial) <- paste0("Group ", plot_group)
    AUCs[[grep("partial", plot_mods)]] <- partial
  }
  # if strict invariance is indicated and `out` does contain strict inv. results
  if ("strict" %in% plot_mods && (!(is.integer(str_ind) && length(str_ind) == 0L))) {
    strict <- c()
    # for each group, extract SE and SP values, pad with 0 and 1, plot ROC
    for (g in plot_group) {
      SEs <- out[str_ind[1]][[1]][g,]
      SPs <- out[str_ind[2]][[1]][g,]
      labg <- ifelse(g == 1, "Reference", paste0("Group ", g))
      plot(y = c(0, SEs, 1), x = c(0, 1 - SPs,1), type = "l",
           ylim = c(0,1), xlim = c(0,1), xlab = labx, ylab = laby,
           main = paste0("ROC: ", labg," (Strict Invariance)"))
      # compute AUC for the group
      strict <- cbind(strict, auc(SEs, SPs))
    }
    colnames(strict) <- paste0("Group ", plot_group)
    AUCs[[grep("strict", plot_mods)]] <- strict
  }

  names(AUCs) <- c(paste0("AUC under ", plot_mods, " invariance"))
  return(AUCs)
}

