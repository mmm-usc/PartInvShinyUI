auc <- function(SE, SP){
  TPR <- SE
  FPR <- 1 - SP
  dFPR <- c(diff(FPR), 0)
  dTPR <- c(diff(TPR), 0)
  sum(TPR * dFPR) + sum(dTPR * dFPR)/2
}

#' @title 
#' Return a list containing sensitivity (SE) and specificity (SP) values
#' 
#' @name 
#' return_SE_SP
#'
#' @description
#' \code{return_SE_SP} takes in a CFA fit object or vectors of model parameter
#' estimates and returns a list object containing SE and SP values under partial
#' and/or strict invariance.
#'  
#' @param out list object
#' 
#' @return The output will be a list.
#' @export
return_SE_SP <- function(cfa_fit,
                         PartInv_fit,
                         pmix, 
                         from = 0.01,
                         to = 0.25,
                         by = 0.01,
                         cutoffs_from = NULL,
                         cutoffs_to = NULL,
                         mod_names = c("partial", "strict"), 
                         alpha = NULL,
                         lambda = NULL, 
                         theta = NULL,
                         psi = NULL,
                         nu = NULL, ...) {
  
  if(!missing(cfa_fit) && missing(PartInv_fit)) {
    est <- format_cfa_partinv(cfa_fit, comp = "est")
    n_g <- cfa_fit@Data@ngroups # number of groups
    psi <- est$psi
    lambda <- est$lambda
    theta <- est$theta
    alpha <- est$alpha
    nu <- est$nu
  }
  if(missing(cfa_fit) && !missing(PartInv_fit)) {
    n_g <- length(PartInv_fit$labels) # number of groups
    psi <- eval(PartInv_fit$functioncall$psi)
    lambda <- eval(PartInv_fit$functioncall$lambda)
    theta <- eval(PartInv_fit$functioncall$theta)
    alpha <- eval(PartInv_fit$functioncall$alpha)
    nu <- eval(PartInv_fit$functioncall$nu)
  }
  
  if(missing(cfa_fit) & missing(PartInv_fit) & !is.null(alpha) & !is.null(psi) & 
     !is.null(lambda) & !is.null(theta)  & !is.null(nu)) {
    n_g <- length(alpha) # number of groups
  }
  
  propsels <- seq(from = from, to = to, by = by)
  use <- "propsels"
  xl <- "Proportion of selection"
  rangeVals <- propsels
  
  if ((is.null(cutoffs_from) && !is.null(cutoffs_to)) ||
      (!is.null(cutoffs_from) && is.null(cutoffs_to))) {
  }
  
  # if the user provided the max and min cutoff values, update rangeVals with 
  # a range of cutoffs
  if (!is.null(cutoffs_from) && !is.null(cutoffs_to)) {
    cutoffs <- seq(from = cutoffs_from, to = cutoffs_to, by = by)
    rangeVals <- cutoffs
    xl <- "Thresholds" # for the plots later
    use <- "cutoffs"
  }
  
  
  # if the user did not provide labels, or provided the wrong number of labels,
  if ((is.null(labels) || (length(labels) != n_g)) && !missing(cfa_fit)) {
    labels <- cfa_fit@Data@group.label
    labels <- paste(labels, c("(reference)", rep("(focal)", n_g - 1)))
  }
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
    # if the user did not provide cutoff values
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
      for(i in cai) {
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
#' @param return_SE_SP_output output of `return_SE_SP()`.
#' @param plot_mods vector of strings indicating invariance conditions.
#' @param plot_group vector of integers indicating the groups of interest. 
#'        `NULL` by default.
#' 
#' @return A list of of length plot_mods
#'
#' @examples 
#' # Multiple groups, multiple dimensions
lambda_matrix <- matrix(0, nrow = 15, ncol = 1)
lambda_matrix[1:15, 1] <- c(0.68, 0.79, -0.39, 0.74, 0.59, 0.46, 0.78, -0.30,
                            0.59, 0.59, 0.64, 0.66, 0.59, 0.63, 0.64);
nu_matrix <- nu_matrix1 <- nu_matrix2 <- nu_matrix3 <-
  matrix(0, nrow = 15, ncol = 1)
nu_matrix[1:15, 1] <- c(3.6, 3.1, 2.7, 2.9, 2.5, 2.1, 3.45, 2.62, 3.2, 2.84,
                        3.51, 3.26, 2.45, 3.39, 2.47);
nu_matrix1[1:15, 1] <- c(3.9, 3.1, 2.7, 2.9, 2.5, 2.1, 3.45, 2.62, 3.2, 2.84,
                         3.51, 3.26, 2.45, 3.76, 2.81);
nu_matrix2[1:15, 1] <- c(3.6, 3.1, 2.7, 2.9, 2.5, 2.1, 3.45, 2.62, 3.6, 3.18,
                         3.51, 3.54, 2.45, 3.39, 2.81);
nu_matrix3[1:15, 1] <- c(3.6, 3.1, 2.7, 2.6, 2.5, 2.1, 3.45, 2.62, 3.2, 2.84,
                         3.51, 3.26, 2.45, 3.39, 2.81);
theta_matrix <- c(0.35, 0.62, 0.83, 0.61, 0.81, 0.87, 0.39, 1.05, 0.84, 0.92,
                  0.36, 0.66, 0.8, 0.66, 0.9);
theta_matrix1 <- c(0.61, 0.62, 0.83, 0.61, 0.81, 0.5, 0.7, 1.05, 0.84, 0.92,
                   0.61, 0.66, 0.8, 0.54, 0.9);
theta_matrix2 <- c(0.61, 0.62, 0.826, 0.61, 0.81, 0.87, 0.5, 1.05, 0.84,
                   0.92, 0.61, 0.66, 0.8, 0.66, 0.9);
theta_matrix3 <- c(0.61, 0.62, 0.826, 0.61, 0.81, 0.5, 0.7, 1.05, 0.84, 0.92,
                   0.61, 0.66, 0.8, 0.66, 0.9);
out <- return_SE_SP(from = 0.00001, to = 0.999999999999999, by = 0.001,
                    pmix = c(1/4, 1/4, 1/4, 1/4),
                    alpha = list(0, -0.70, -1.05, -1.10),
                    psi = list(1, 1.2, 1.29, 1.3),
                    nu = list(nu_matrix, nu_matrix1, nu_matrix2, nu_matrix3),
                    lambda = list(lambda_matrix, lambda_matrix, lambda_matrix,
                                  lambda_matrix),
                    theta = list(theta_matrix, theta_matrix1, theta_matrix2,
                                 theta_matrix3),
                    plot_contour = FALSE, show_mi_result = TRUE,
                    mod_names = c("strict", "partial"))
roc_auc_PartInv(out, plot_mods = c("strict", "partial"))
roc_auc_PartInv(out, plot_mods = c("partial"))
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
    out <- return_SE_SP(call_cfa, mod_names = plot_mods)
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
  
  labx <- "FPR (1-SP)"
  laby <- "TPR (SE)"
  SEs <- SPs <- data.frame()
  AUCs <- c()
  
  par_ind <- grep("par", names(out))
  str_ind <- grep("str", names(out))
  # if partial invariance is indicated and `out` does contain partial inv. results
  if("partial" %in% plot_mods && (!(is.integer(par_ind) && length(par_ind) == 0L))) {
    partial <- c()
    # for each group, extract SE and SP values, pad with 0 and 1, plot ROC
    for(g in plot_group){
      SEs <- out[par_ind[1]][[1]][g,]
      SPs <- out[par_ind[2]][[1]][g,]
      labg <- ifelse(g == 1, "Reference", paste0("Group ", g))
      plot(y = c(0, SEs, 1), x = c(0, 1 - SPs,1), type = "l",
           ylim = c(0,1), xlim = c(0,1), xlab = labx, ylab = laby,
           main = paste0("ROC: ", labg," (Partial Invariance)"))
      # compute AUC for the group 
      partial <- cbind(partial, auc(SEs, SPs))
    }
    colnames(partial) <- paste0("Group ", plot_group)
    AUCs[[grep("partial", plot_mods)]] <- partial
  }
  # if strict invariance is indicated and `out` does contain strict inv. results
  if("strict" %in% plot_mods && (!(is.integer(str_ind) && length(str_ind) == 0L))) {
    strict <- c()
    # for each group, extract SE and SP values, pad with 0 and 1, plot ROC
    for(g in plot_group){
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


