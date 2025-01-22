#' @title
#' Impact of deleting biased item(s) on classification accuracy indices (CAI)
#'
#' @name
#' item_deletion_h
#'
#' @description
#' \code{item_deletion_h} computes effect size indices that quantify the impact
#'  of (and changes in the impact of) measurement bias on CAI if an item is 
#'  dropped vs. retained.
#'  Comparisons are made between CAI computed for the reference group and
#'  expected CAI computed for the focal group; between CAI computed under strict
#'  factorial invariance (SFI) vs. partial factorial invariance (PFI); and 
#'  aggregate CAI computed for item subsets.
#' @param cfa_fit CFA model output from lavaan.
#' @param propsel Proportion of selection. If missing, computed using `cut_z`.
#' @param cut_z Pre-specified cutoff score on the observed composite. Ignored 
#'     when `propsel` has an input.
#' @param weights_item A vector of item weights.
#' @param weights_latent A vector of latent factor weights.
#' @param alpha A list of length `g` containing `1 x d` latent factor mean
#'     vectors where `g` is the number of groups and `d` is the number of latent
#'     dimensions. The first element is assumed to belong to the reference group.
#' @param psi A list of length `g` containing `d x d` latent factor
#'     variance-covariance matrices where `g` is the number of groups and `d` is
#'     the number of latent dimensions. The first element is assumed to belong
#'     to the reference group.
#' @param lambda A list of length `g` containing `n_i x d` factor loading matrices
#'     where `g` is the number of groups, `d` is the number of latent dimensions,
#'     and `n` is the number of items . The first element is assumed
#'     to belong to the reference group.
#' @param nu A list of length `g` containing `1 x n_i` measurement intercept
#'     vectors where `g` is the number of groups and `n_i` is the number of 
#'     items. The first element is assumed to belong to the reference group.
#' @param theta A list of length `g` containing `1 x n_i` vectors or `n_i x n_i`
#'     matrices of unique factor variances and covariances, where `g` is the
#'     number of groups and `n_i` is the number of items . The first
#'     element is assumed to belong to the reference group.
#' @param alpha_r,alpha_f,nu_r,nu_f,Theta_r,Theta_f,psi_r,psi_f,lambda_r,lambda_f,pmix_ref
#'     Deprecated; included for backward compatibility. With two groups, '_r' 
#'     and '_f' suffixes refer to the reference group and the focal group.
#' @param pmix List of length `g` containing the mixing proportions of each
#'     group. If `NULL`, defaults to `1/g` for each group (i.e., equal sizes).
#' @param plot_contour Logical; whether the contour of the two populations
#'     should be plotted; default to `TRUE`.
#' @param n_dim Number of dimensions, 1 by default.
#' @param n_i_per_dim A vector containing the number of items per dimension;
#'     `NULL` by default. If `n_dim` \eqn{> 1} and \code{n_i_per_dim = NULL}, 
#'      subscales are assumed to have an equal number of items.
#' @param delete_items A vector; default to `NULL`. If `NULL`, only items 
#'     determined to contain bias will be considered for deletion.
#' @param delete_one_cutoff User-specified cutoff to use in delete-one scenarios.
#'     `NULL` by default; if `NULL`, PS on the full item set will be used.
#' @param labels A character vector with two elements to label the reference
#'     and the focal group on the graph.
#' @param plot_contour Logical; whether the contour of the populations should be
#'     plotted; `TRUE` by default.
#' @param reference Optional argument for specifying the reference group.
#'     Currently only functional when cfa_fit is provided. If providing parameter 
#'     estimates instead, reorder estimates such that the first estimates belong to 
#'     the reference group.
#' @param digits Number of digits for rounding. 3 by default.
#' @param ... Other arguments for \code{\link[graphics]{contour}}.
#' @return An object of class `itemdeletion` containing 13 elements.
#'     \item{AI}{A data frame storing Adverse Impact (AI) values under partial 
#'      invariance by groups.}
#'     \item{ACAI}{A list storing aggregate PS, SR, SE, SP computed for the full
#'      set of items and item subsets excluding biased or user specified items
#'      under partial invariance (PFI).}
#'     \item{h_acai_p}{A list storing Cohen's h values quantifying the impact of 
#'      deleting each item considered in the `ACAI` table.}
#'     \item{h_acai_s_p}{A list storing Cohen's h values quantifying the 
#'      discrepancy between ACAI under SFI vs. ACAI under PFI.}
#'     \item{delta_h_acai_s_p}{A list storing delta h values quantifying the 
#'      impact of deleting an item on discrepancies in h_acai_s_p.}
#'     \item{h_R_EF}{A list storing Cohen's h values quantifying the discrepancy
#'      between observed CAI for the reference group and the expected CAI 
#'      computed for the focal group if it matched the distribution of the 
#'      reference group (Efocal) under PFI.}
#'     \item{delta_h_R_EF}{A list storing delta h values quantifying the impact
#'      of deleting an item on discrepancies in h_R_Ef.}
#'     \item{h_s_p}{A list containing Cohen's h values quantifying the
#'     discrepancy between CAI under SFI vs. PFI for item subsets.}
#'     \item{delta_h_s_p}{A list storing delta h values quantifying the impact
#'      of deleting an item on discrepancies in h_s_p.}
#'     \item{PartInv}{A list containing PartInv() outputs under PFI and SFI for
#'      each item deletion scenario.}
#'     \item{items}{A vector with items considered for deletion.}
#'     \item{function_call}{Function call to item_deletion_h().}
#'     \item{digits}{Number of digits utilized for rounding.}
#' @examples
#' set.seed(7) 
#' # Simulate random data to fit a multigroup CFA, invariance across languages
#' library(lavaan)
#' sim_m <-
#'   "f =~ c(1, .7, 1) * x1 + c(.8, 1.1, 1) * x2 + 1 * x3 + 1 * x4 + 1 * x5
#'    f ~~ c(1, 1.3, 1.5) * f
#'    f ~  c(0, .5, 1) * 1
#'    x1 ~ c(0, .3, 0) * 1
#'    x3 ~ c(.3, 0, -.3) * 1
#'    x1 ~~ c(1, .5, 1) * x1"
#' dat_sim <- simulateData(sim_m, sample.nobs = c(120, 90, 50))
#' dat_sim$group <- ifelse(dat_sim$group == 1, "English",
#'                  ifelse(dat_sim$group == 2, "Japanese",
#'                  ifelse(dat_sim$group == 3, "Swahili", NA)))
#' fit_sim <- lavaan::cfa(model = sim_m, data = dat_sim, group = "group")
#' del <- item_deletion_h(cfa_fit = fit_sim, propsel = .05)
#' del # formatted
#' summary(del) #formatted with additional output
#' del$AI # can access all outputs without rounding or formatting
#' del$ACAI
#' del$PartInv_outputs$`Full item set`$summary
#' # choose Japanese as the reference group and use a cutoff:
#' del <- item_deletion_h(cfa_fit = fit_sim, cut_z = 15, reference = "Japanese")
#' del
#' 
#' # Multidimensional example
#' l_mat <- matrix(0, nrow = 5, ncol = 2)
#' l_mat[1:2, 1] <- c(.322, .655); l_mat[3:5, 2] <- c(.398, .745, .543)
#' multi_dim <- item_deletion_h(propsel = .05, n_dim = 5, 
#'     weights_item = c(1/4, 1/4, 1/6, 1/6, 1/6),
#'     weights_latent = c(0.5, 0.5), alpha_r = c(0, 0), alpha_f = c(-0.3, 0.1),
#'     psi_r = matrix(c(1, 0.5, 0.5, 1), nrow = 2), lambda_r = l_mat,
#'     nu_r = c(.225, .025, .010, .240, .125),
#'     nu_f = c(.225, -.05, .240, -.025, .125), Theta_r = diag(1, 5),
#'     Theta_f = diag(c(1, .95, .80, .75, 1)), plot_contour = TRUE)
#' print(multi_dim)
#' # Single dimension example
#' single_dim <- item_deletion_h(propsel = .10, weights_item = c(1, 0.9, 1, 1),
#'     weights_latent = 0.9, alpha_r = 0.5, alpha_f = 0, psi_r = 1,
#'     lambda_r = c(.3, .5, .9, .7), nu_r = c(.225, .025, .240, .240),
#'     nu_f = c(.225, -.05, .240, -.025), Theta_r = diag(.96, 4), n_dim = 1, 
#'     plot_contour = TRUE)
#' print(single_dim)
#' # Using cfa_fit
#' HS <- HolzingerSwineford1939
#' HS$sex <- as.factor(HS$sex)
#' HS.model <- ' visual  =~ x1 + x2 + x3
#'               textual =~ x4 + x5 + x6
#'               speed   =~ x7 + x8 + x9 '
#' fit <- cfa(HS.model, data = HS, group = "sex")
#' item_deletion_h(cfa_fit = fit, propsel = .05)
#' @export
item_deletion_h <- function(cfa_fit = NULL,
                            propsel = NULL,
                            cut_z = NULL,
                            weights_item = NULL,
                            weights_latent = NULL,
                            alpha = NULL, psi = NULL, lambda = NULL, theta = NULL, nu = NULL,
                            pmix = NULL,
                            pmix_ref = 0.5,
                            plot_contour = TRUE,
                            labels = NULL, #c("Reference", "Focal"),
                            n_dim = 1,
                            n_i_per_dim = NULL,
                            delete_items = NULL,
                            delete_one_cutoff = NULL,
                            alpha_r = NULL, alpha_f = alpha_r,
                            psi_r = NULL, psi_f = psi_r,
                            lambda_r = NULL, lambda_f = lambda_r,
                            nu_r = NULL, nu_f = nu_r,
                            Theta_r = NULL, Theta_f = Theta_r, reference = NULL,
                            digits = 3,
                            ...) {
  functioncall <- match.call()            
  CAIs <- c("TP", "FP", "TN", "FN", "PS", "SR", "SE", "SP")
  ACAI_labs <- paste0(CAIs, "*")
  # make adjustments for formatting and backward compatibility
  pl <- prep_params(
    cfa_fit, propsel, cut_z, weights_item, weights_latent, alpha, psi, lambda, 
    theta, nu, pmix, pmix_ref, plot_contour, labels, n_dim = n_dim, 
    n_i_per_dim = n_i_per_dim, delete_items = delete_items, 
    delete_one_cutoff = delete_one_cutoff, alpha_r, alpha_f,
    phi_r = NULL, phi_f = NULL, psi_r, psi_f, lambda_r, lambda_f, tau_r = NULL, 
    tau_f = NULL, kappa_r = NULL, kappa_f = NULL, nu_r, nu_f, Theta_r, Theta_f, 
    reference, custom_colors = NULL, PartInv_fit = NULL)
  
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

  # determine which set of items will be returned
  biased_items <- c()
  if (is.null(delete_items)) { # default: return only the biased items.
      biased_items <- determine_biased_items(lambda, nu, theta)
      biased_items <-  setdiff(biased_items, which(weights_item == 0))
  } else {
    if (!all(delete_items == floor(delete_items))) {
      stop("'delete_items' should only contain integers corresponding to item indices.")}
    if (!all(delete_items < n_i + 1)) {
      stop("'delete_items' cannot take integers > the scale length.")}
    biased_items <- delete_items
  }

  labs <-  c(paste0("|", 1:n_i))
  if (is.null(delete_items)) {
    final_set <- paste0("|", biased_items)
  } else {
    final_set <- paste0("|", delete_items)
  }
  
  acai_p <- acai_s <- create_list_of_dfs(
    ls_len = num_g - 1, ncol = 8, nrow = n_i + 1, 
    df_cn = paste0(CAIs, "*"), df_rn = c("Full", labs), groups = labels[-1])
  h_acai_s_p <- create_list_of_dfs(
    ls_len = num_g - 1, ncol = 8, nrow = n_i + 1, 
    df_cn = paste0("h(", CAIs, "*)"), df_rn = c("Full", labs), groups = labels[-1])
  h_R_Ef <- create_list_of_dfs(
    ls_len = num_g - 1, ncol = 8, nrow = n_i + 1, 
    df_cn = paste0("h(", CAIs, ")"), df_rn = c("r_Ef", paste0("r_Ef", labs)), 
    groups = labels[-1])
  h_acai_p <- create_list_of_dfs(
    ls_len = num_g - 1, ncol = 8, nrow = n_i, df_cn = paste0("h(", CAIs, "*)"),
    df_rn = labs, groups = labels[-1])
  delta_h_R_Ef <- create_list_of_dfs(
    ls_len = num_g - 1, ncol = 8, nrow = n_i,
    df_cn = paste0("\u0394h(", CAIs, ")"), df_rn = labs, groups = labels[-1])
  delta_h_acai_s_p <- create_list_of_dfs(
    ls_len = num_g - 1, ncol = 8, nrow = n_i, 
    df_cn = paste0("\u0394h(", CAIs, "*)"), df_rn = labs, groups = labels[-1])
  delta_h_s_p <- create_list_of_dfs(
    ls_len = num_g, ncol = 8, nrow = n_i, 
    df_cn = paste0("\u0394h(", CAIs, ")"), df_rn = labs, groups = labels)
  h_s_p <- create_list_of_dfs(
    ls_len = num_g, ncol = 8, nrow = n_i + 1 , 
    df_cn = paste0("h(", CAIs, ")"), df_rn = c("Full", labs), groups = labels)
  
  store_str <- store_par <- vector(mode = "list", n_i + 1)
  names(store_str) <- names(store_par) <- c("Full item set", labs)
  AI_ratios <- as.data.frame(matrix(ncol = num_g, nrow = n_i + 1), 
                             row.names = c("Full", labs))

  colnames(AI_ratios) <- c("(SFI)", paste0(labels[-1]))
  out_par <- c("propsel", "cutpt_xi", "cutpt_z", "summary", "bivar_data", 
               "ai_ratio", "labels", "functioncall")
  out_str <- c(paste0(out_par[1:6], "_mi"), out_par[7:8])
  
  # Call PartInv with the full item set under partial and strict invariance ###
  store_par[[1]] <- 
    PartInv(cfa_fit = NULL, propsel = propsel, cut_z = cut_z, 
            weights_item = weights_item, weights_latent = weights_latent, 
            alpha = alpha, psi = psi, nu = nu, theta = theta, lambda = lambda,
            pmix = pmix, plot_contour = plot_contour, labels = labels,
            show_mi_result = TRUE)
  class(store_par[[1]]) <- "PartInv"
  store_str[[1]] <- store_par[[1]][out_str]
  class(store_str[[1]]) <- "PartInv"

  # h: strict vs. partial invariance (full item set) for all groups
  temp <- as.data.frame(Map(cohens_h, 
                            store_str[[1]]$summary_mi, 
                            store_par[[1]]$summary[1:num_g]))
  h_s_p <- update_rows_in_lists_of_dfs(h_s_p, temp, ind = 1)
 
  temp <- apply(as.data.frame(
    store_par[[1]]$summary[,(num_g + 1):(2 * num_g - 1)]), MARGIN = 2,
    FUN = function(x) cohens_h(store_par[[1]]$summary[, 1], x))
  temp <- split(temp, colnames(temp))
  h_R_Ef <- update_rows_in_lists_of_dfs(h_R_Ef, temp, ind = 1)

  # Compute aggregate CAI on the full item set
  temp <- get_aggregate_CAI(pmix, store_par[[1]]$summary, inv_cond = "partial")
  acai_p <- update_rows_in_lists_of_dfs(acai_p, temp, ind = 1)
  temp <- get_aggregate_CAI(pmix, store_str[[1]]$summary_mi, inv_cond = "strict")
  acai_s <- update_rows_in_lists_of_dfs(acai_s, temp, ind = 1)

  # compute cohen's h for the difference between the strict and partial inv. conditions
  # (on the first row of the data frames in each element of the two lists)
  temp <- Map(function(df_s, df_p) {
    cohens_h(df_s[1, ], df_p[1, ]) 
  }, acai_s, acai_p)
  # store each vector in the first row of the corresponding dataframe in h_acai_s_p
  h_acai_s_p <- update_rows_in_lists_of_dfs(h_acai_s_p, temp, ind = 1)
  
  # h: difference between strict and partial invariance for aggregate CAI
   AI_ratios[1,] <- as.vector(c(1, store_par[[1]]$ai_ratio), mode = "double")

  # If no cutoff was provided, set propsel based on PartInv output with all items
  if(is.null(delete_one_cutoff)) {
    propsel_p <- store_par[[1]]$propsel
    cut_z <- NULL
  } else {
    cut_z <- delete_one_cutoff
    propsel_p <- NULL
  }
  
  # Item deletion scenarios ####
  for (i in seq_len(length(weights_item) + 1)[-1]) {
    # Assign a weight of 0 to the item to be deleted (indexed at i - 1), and
    # redistribute the weight from this item across the retained items
    take_one_out <- redistribute_weights(weights_item, n_dim = n_dim,
                                         n_i_per_dim = n_i_per_dim, del_i = i - 1)
    
    # Call PartInv with the new weights ####
    store_par[[i]] <- 
      PartInv(cfa_fit = NULL, propsel = propsel_p, cut_z = cut_z,
              weights_item = take_one_out, weights_latent = weights_latent,
              alpha = alpha, psi = psi, nu = nu, theta = theta, lambda = lambda,
              pmix = pmix, plot_contour = plot_contour, labels = labels,
              show_mi_result = TRUE)
    class(store_par[[i]]) < "PartInv"
    store_str[[i]] <- store_par[[i]][out_str]
    class(store_str[[i]]) <- "PartInv"
    
    # Check whether improvements in ACAI may be misleading due pmix
    err_improv_acai(i = i, s_full = store_par[[1]]$summary,
                    s_del1 = store_par[[i]]$summary, num_g = num_g)
    err_improv_acai(i = i, s_full = store_str[[1]]$summary_mi,
                    s_del1 = store_str[[i]]$summary_mi, num_g = num_g)
    
    # h: strict vs. partial invariance (delete-one item set) for all groups
    temp <- as.data.frame(Map(cohens_h, 
                              store_str[[i]]$summary_mi, 
                              store_par[[i]]$summary[1:num_g]))
    h_s_p <- update_rows_in_lists_of_dfs(h_s_p, temp, ind = i)

    # delta h: comparing CAI under strict vs. partial invariance when item i is 
    # deleted (i.e. the change in h_s_p_ref and h_s_p_foc) for all groups
    delta_h_s_p <- Map(function(df, df2) {
      df[i - 1, ] <- delta_h(df2[1, ], df2[i, ])
      df
    }, delta_h_s_p, h_s_p)
    
    # h: difference in CAI under partial invariance for the ref group vs. for 
    # the expected CAI for the focal group with the full item set
    
    temp <- apply(as.data.frame(
      store_par[[i]]$summary[, (num_g + 1):(2 * num_g - 1)]), MARGIN = 2,
      FUN = function(x) cohens_h(store_par[[i]]$summary[, 1], x))
    temp <- split(temp, colnames(temp))
    h_R_Ef <- update_rows_in_lists_of_dfs(h_R_Ef, temp, ind = i)
    
    
    # change in h_R_Ef_del when item i is deleted (under partial invariance)
    delta_h_R_Ef <- Map(function(df, df2) {
      df[i - 1, ] <- delta_h(df2[1, ], df2[i, ])
      df
    }, delta_h_R_Ef, h_R_Ef)
    
    # Weight the aggregate SR, SE, SP indices under partial and strict invariance
    temp <- get_aggregate_CAI(pmix, store_par[[i]]$summary, inv_cond = "partial")
    acai_p <- update_rows_in_lists_of_dfs(acai_p, temp, ind = i)
    temp <- get_aggregate_CAI(pmix, store_str[[i]]$summary_mi, inv_cond = "strict")
    acai_s <- update_rows_in_lists_of_dfs(acai_s, temp, ind = i)
    
    # compute cohen's h for the difference between the strict and partial inv. conditions
    # (on the i-th row of the data frames in each element of the two lists)
    temp <- Map(function(df_s, df_p) {
      cohens_h(df_s[i, ], df_p[i, ]) 
    }, acai_s, acai_p)

    # store each vector in the i-th row of the corresponding dataframe in h_acai_s_p
    h_acai_s_p <- update_rows_in_lists_of_dfs(h_acai_s_p, temp, ind = i)
    
    # h: change in aggregate CAI when an item is deleted under partial invariance
    temp <- lapply(acai_p, function(df) cohens_h(df[1, ], df[i, ]))
    h_acai_p <- update_rows_in_lists_of_dfs(h_acai_p, temp, ind = i - 1) 
   
    delta_h_acai_s_p <- Map(function(df, df2) {
     df[i-1,] <- delta_h(df2[1,], df2[i,])
     df
   }, delta_h_acai_s_p, h_acai_s_p)
   AI_ratios[i,] <- as.vector(c(1, store_par[[i]]$ai_ratio), mode = "double")
  }
  
  out <- list(
  "AI" = AI_ratios,
  "ACAI" = acai_p,
  "h_acai_p" = h_acai_p, 
  "h_acai_s_p" = h_acai_s_p,
  "delta_h_acai_s_p" = delta_h_acai_s_p, 
  "h_R_Ef" = h_R_Ef, 
  "delta_h_R_Ef" = delta_h_R_Ef, 
  "h_s_p" = h_s_p, 
  "delta_h_s_p" = delta_h_s_p,
  "PartInv_outputs" = store_par,
  "items" = final_set,
  "function_call" = functioncall,
  "digits" = ifelse(is.null(digits), 3, digits)
     )
  
  class(out) <- "itemdeletion"
  return(out)
}
