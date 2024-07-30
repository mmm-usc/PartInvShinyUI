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
#'                   plot_contour = FALSE, show_mi_result = TRUE,
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
#'
#' # when cfa output is provided as input
#' library(lavaan)
#' #data("HolzingerSwineford1939", package = "lavaan")
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
return_SE_SP <- function(cfa_fit,
                         PartInv_fit,
                         pmix, 
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
                         nu = NULL, ...) {
  # if the user provided the output from cfa()
  if(!missing(cfa_fit) && missing(PartInv_fit)) {
    est <- format_cfa_partinv(cfa_fit, comp = "est")
    n_g <- cfa_fit@Data@ngroups # number of groups
    psi <- est$psi
    lambda <- est$lambda
    theta <- est$theta
    alpha <- est$alpha
    nu <- est$nu
  }
  # if the user provided the output from PartInv
  if(missing(cfa_fit) && !missing(PartInv_fit)) {
    n_g <- length(PartInv_fit$labels) # number of groups
    psi <- eval(PartInv_fit$functioncall$psi)
    lambda <- eval(PartInv_fit$functioncall$lambda)
    theta <- eval(PartInv_fit$functioncall$theta)
    alpha <- eval(PartInv_fit$functioncall$alpha)
    nu <- eval(PartInv_fit$functioncall$nu)
  }
  # if the user provided the relevant matrices of parameter estimates instead
  # of providing the cfa() or PartInv() outputs
  if(missing(cfa_fit) & missing(PartInv_fit) & !is.null(alpha) & !is.null(psi) & 
     !is.null(lambda) & !is.null(theta)  & !is.null(nu)) {
    n_g <- length(alpha) # number of groups
  }
  
  if(!missing(cfa_fit) && !missing(PartInv_fit) && 
     any(!is.null(c(alpha, lambda, theta, psi, nu)))) {
    warning('Please only provide the output of cfa(), the output of PartInv(), or vectors of parameter estimates.')
  }
  
  propsels <- seq(from = from, to = to, by = by)
  use <- "propsels" # set to use proportion of selection by default
  xl <- "Proportion of selection"
  rangeVals <- propsels

  # if the user provided the max and min cutoff values, update rangeVals with
  # a range of cutoffs
  if (!is.null(cutoffs_from) && !is.null(cutoffs_to)) {
    cutoffs <- seq(from = cutoffs_from, to = cutoffs_to, by = by)
    rangeVals <- cutoffs
    xl <- "Thresholds" # for the plots later
    use <- "cutoffs" # update to use cutoffs instead
  }

  # if the user did not provide labels, or provided the wrong number of labels,
  # and the cfa output was provided, extract relevant info from the output
  # to create labels
  if ((is.null(labels) || (length(labels) != n_g)) && !missing(cfa_fit)) {
    labels <- cfa_fit@Data@group.label
    labels <- paste(labels, c("(reference)", rep("(focal)", n_g - 1)))
  }
  # if the user did not provide labels, or provided the wrong number of labels,
  # but the cfa output was not provided, make up group labels with 'Group'
  if ((is.null(labels) || (length(labels) != n_g)) && missing(cfa_fit)) {
    labels <- paste(rep("Group"), seq(1:n_g))
    labels <- paste(labels, c("(reference)", rep("(focal)", n_g - 1)))
  }

  cai <- c("Sensitivity", "Specificity")
  ls_mat <- matrix(NA, ncol = length(rangeVals), nrow = n_g,
                   dimnames = list(labels, rangeVals))
  ls_names <- c(t(outer(cai, Y = mod_names, FUN = paste, sep = "_")))
  vals <- rep(list(ls_mat), length(ls_names))
  names(vals) <- ls_names

  # if pmix is missing, assume equal mixing proportions
  if (missing(pmix)) pmix <- as.matrix(c(rep(1 / n_g, n_g)), ncol = n_g)
  pmix <- as.vector(pmix)

  ylabs <- ""
  mains <- ""

  # call PartInv with each proportion of selection and store CAI in the list of
  # data frames
  for (p in seq_along(rangeVals)) {
    # if the user provided cutoff values
    if (use == "cutoffs") {
      suppressWarnings({
        pinv <- PartInv(cut_z = cutoffs[p],
                        psi = psi,
                        lambda = lambda,
                        theta = theta,
                        alpha = alpha,
                        nu = nu,
                        pmix = pmix,
                        plot_contour = FALSE,
                        labels = labels,
                        show_mi_result = TRUE)
      })
    }
    # if the user did not provide cutoff values, follow the default
    if (use == "propsels") {
      suppressWarnings({
        pinv <- PartInv(propsel = propsels[p],
                        psi = psi,
                        lambda = lambda,
                        theta = theta,
                        alpha = alpha,
                        nu = nu,
                        pmix = pmix,
                        plot_contour = FALSE,
                        labels = labels,
                        show_mi_result = TRUE)
      })
    }

    num_comb <- length(cai) * length(mod_names) + 1 # for specifying the
    # index within vals
    ind <- 1

    while (ind < num_comb) {
      for (i in cai) {
        # for each specified invariance condition
        for (j in seq_along(mod_names)) {
          # if the specified invariance condition is partial invariance,
          vals[[ind]][, p] <-
            ifelse(rep(mod_names[j] == "partial", n_g),
                   as.numeric(pinv$summary[i, 1:n_g]),
                   as.numeric(pinv$summary_mi[i, 1:n_g]))

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


