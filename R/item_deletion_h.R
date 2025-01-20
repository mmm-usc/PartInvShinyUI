#' @title
#' Impact of deleting biased item(s) on classification accuracy indices (CAI)
#'
#' @name
#' item_deletion_h
#'
#' @description
#' \code{item_deletion_h} computes effect size indices that quantify the impact
#'  of (and changes in the impact of) measurement bias on classification accuracy
#'  indices (CAI) such as TP and SE if an item is dropped vs. included in analyses.
#'  Comparisons are made between CAI computed for the reference group and
#'  expected CAI computed for the focal group; between CAI computed under strict
#'  factorial invariance (SFI) vs. CAI computed under partial factorial
#'  invariance (PFI); and between aggregate CAI computed for item subsets.
#' @param cfa_fit CFA model output from lavaan.
#' @param propsel Proportion of selection. If missing, computed using `cut_z`.
#' @param cut_z Pre-specified cutoff score on the observed composite. This
#'        argument is ignored when `propsel` has an input.
#' @param weights_item A vector of item weights.
#' @param weights_latent A  vector of latent factor weights.
#' @param alpha A list of length `g` containing `1 x d` latent factor mean
#'     vectors where `g` is the number of groups and `d` is the number of latent
#'     dimensions. The first element is assumed to belong to the reference group.
#' @param psi A list of length `g` containing `d x d` latent factor
#'     variance-covariance matrices where `g` is the number of groups and `d` is
#'     the number of latent dimensions. The first element is assumed to belong
#'     to the reference group.
#' @param lambda A list of length `g` containing `n_i x d` factor loading matrices
#'     where `g` is the number of groups, `d` is the number of latent dimensions,
#'     and `n` is the number of items in the scale. The first element is assumed
#'     to belong to the reference group.
#' @param nu A list of length `g` containing `1 x n_i` measurement intercept
#'     vectors where `g` is the number of groups and `n_i` is the number of items
#'     in the scale. The first element is assumed to belong to the reference
#'     group.
#' @param theta A list of length `g` containing `1 x n_i` vectors or `n_i x n_i`
#'     matrices of unique factor variances and covariances, where `g` is the
#'     number of groups and `n_i` is the number of items in the scale. The first
#'     element is assumed to belong to the reference group.
#' @param alpha_r,alpha_f,nu_r,nu_f,Theta_r,Theta_f,psi_r,psi_f,lambda_r,lambda_f,pmix_ref
#'        Deprecated; included only for backward compatibility. When comparing two
#'        groups, parameters with the '_r' suffix refer to the reference group while
#'        parameters with the '_f' suffix refer to the focal group.
#' @param pmix List of length `g` containing the mixing proportions of each
#'     group. If `NULL`, defaults to `1/g` for each group (i.e., the populations
#'     have equal size).
#' @param plot_contour Logical; whether the contour of the two populations
#'        should be plotted; default to `TRUE`.
#' @param n_dim Number of dimensions, 1 by default. If the user does not supply
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
#' @param labels A character vector with two elements to label the reference
#'         and the focal group on the graph.
#' @param ... Other arguments passed to the \code{\link[graphics]{contour}}
#'     function.
#' @return An object of class `itemdeletion` containing 13 elements.
#'        \item{ACAI}{A matrix that stores aggregate PS, SR, SE, SP computed for
#'        the full set of items and item subsets excluding biased or user specified
#'        items under PFI.}
#'        \item{h ACAI (deletion)}{A matrix that stores Cohen's h computed for
#'        the impact of deleting each item considered in the `ACAI` table.}
#'        \item{h ACAI SFI-PFI}{A matrix that stores Cohen's h values
#'        quantifying the discrepancy between ACAI under SFI vs. ACAI under PFI.}
#'        \item{delta h ACAI SFI-PFI (deletion)}{A matrix that stores delta h
#'        values quantifying the impact of deleting an item on the discrepancy
#'        between ACAI under SFI vs. ACAI under PFI for subsets of items.}
#'        \item{AI Ratio}{A matrix storing Adverse Impact Ratio values computed
#'        for item subsets by invariance condition.}
#'        \item{h CAI Ref-EF}{A matrix that stores Cohen's h values quantifying
#'        the discrepancy between CAI computed for the reference group and the
#'        expected CAI computed for the focal group if it matched the
#'        distribution of the reference group (Efocal), under PFI for subsets of
#'        items.}
#'        \item{delta h CAI Ref-EF (deletion)}{A matrix that stores delta h
#'        values quantifying the impact of deleting an item on the discrepancy
#'        between CAI of reference vs. Efocal groups under PFI.}
#'        \item{h CAI SFI-PFI}{A list containing two items, `ref` and `foc`
#'        which are matrices storing Cohen's h values quantifying the
#'        discrepancy between CAI under SFI vs. PFI for the reference group and
#'        the focal group respectively, for subsets of items.}
#'        \item{delta h SFI-PFI (deletion)}{A list containing two items,
#'        `ref` and `foc` which are matrices storing delta h values quantifying
#'        the impact of deleting an item on the discrepancy between CAI under
#'        SFI vs. PFI for the reference group and the focal group respectively,
#'        for subsets of items.}
#'        \item{h SFI-PFI by groups}{Two lists (`reference` and `focal`). The lists
#'        contain tables for each item deletion scenario displaying raw CAI
#'        under SFI, under PFI, and the Cohen's h value associated with the
#'        difference between the invariance condition.}
#'        \item{PartInv}{Two lists (`strict` and `partial`), each containing
#'        PartInv() outputs.}
#'        \item{return_items}{A vector containing the items that will be considered
#'        for deletion.}
#' @examples
#' # Multidimensional example
#' lambda_matrix <- matrix(0, nrow = 5, ncol = 2)
#' lambda_matrix[1:2, 1] <- c(.322, .655)
#' lambda_matrix[3:5, 2] <- c(.398, .745, .543)
#'
#' multi_dim <- item_deletion_h(propsel = .05, n_dim = 5,
#'                              weights_item = c(1/4, 1/4, 1/6, 1/6, 1/6),
#'                              weights_latent = c(0.5, 0.5),
#'                              alpha_r = c(0, 0),
#'                              alpha_f = c(-0.3, 0.1),
#'                              psi_r = matrix(c(1, 0.5, 0.5, 1), nrow = 2),
#'                              lambda_r = lambda_matrix,
#'                              nu_r = c(.225, .025, .010, .240, .125),
#'                              nu_f = c(.225, -.05, .240, -.025, .125),
#'                              Theta_r = diag(1, 5),
#'                              Theta_f = diag(c(1, .95, .80, .75, 1)),
#'                              plot_contour = TRUE)
#' print(multi_dim)
#' # Single dimension example
#' single_dim <- item_deletion_h(propsel = .10,
#'                                weights_item = c(1, 0.9, 0.8, 1),
#'                                weights_latent = 0.9,
#'                                alpha_r = 0.5,
#'                                alpha_f = 0,
#'                                psi_r = 1,
#'                                lambda_r = c(.3, .5, .9, .7),
#'                                nu_r = c(.225, .025, .010, .240),
#'                                nu_f = c(.225, -.05, .240, -.025),
#'                                Theta_r = diag(.96, 4),
#'                                n_dim = 1, plot_contour = TRUE)
#' print(single_dim)
#' # Using cfa_fit
#' library(lavaan)
#' HS <- HolzingerSwineford1939
#' HS$sex <- as.factor(HS$sex)
#' HS.model <- ' visual  =~ x1 + x2 + x3
#'               textual =~ x4 + x5 + x6
#'               speed   =~ x7 + x8 + x9 '
#' 
#' fit <- cfa(HS.model, data = HS, group = "sex")
#' item_deletion_h(cfa_fit = fit, propsel = .05)
#' 
#' # Simulate random data to fit a multigroup cfa
#' 
#' set.seed(7)  
#' sim_m <-
#'   "f =~ c(1, .7, 1) * x1 + c(.8, 1.1, 1) * x2 + 1 * x3 + 1 * x4 + 1 * x5
#'    f ~~ c(1, 1.3, 1.5) * f
#'    f ~  c(0, .5, 1) * 1
#'    x1 ~ c(0, .3, 0) * 1
#'    x3 ~ c(.3, 0, -.3) * 1
#'    x1 ~~ c(1, .5, 1) * x1"
#' dat_sim <- simulateData(sim_m, sample.nobs = c(80, 100, 110))
#' fit_sim <- lavaan::cfa(model = sim_m, data = dat_sim, group = "group")
#' item_deletion_h(cfa_fit = fit_sim, propsel = .05)

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
                            labels = NULL,#c("Reference", "Focal"),
                            n_dim = 1,
                            n_i_per_dim = NULL,
                            delete_items = NULL,
                            delete_one_cutoff = NULL,
                            alpha_r = NULL, alpha_f = alpha_r,
                            psi_r = NULL, psi_f = psi_r,
                            lambda_r = NULL, lambda_f = lambda_r,
                            nu_r = NULL, nu_f = nu_r,
                            Theta_r = NULL, Theta_f = Theta_r,
                            ...) {
  functioncall <- match.call()
  CAIs <- c("TP", "FP", "TN", "FN", "PS", "SR", "SE", "SP")
  
  # for backward compatibility with different input names ####
  # pl: parameter list after adjustments
  pl <- unbiasr:::prep_params(
    cfa_fit, propsel, cut_z, weights_item, weights_latent, alpha, psi, lambda, 
    theta, nu, pmix, pmix_ref, plot_contour, labels, n_dim = n_dim, 
    n_i_per_dim = n_i_per_dim, delete_items = delete_items, 
    delete_one_cutoff = delete_one_cutoff, alpha_r, alpha_f,
    phi_r = NULL, phi_f = NULL, psi_r, psi_f, lambda_r, lambda_f, tau_r = NULL, 
    tau_f = NULL, kappa_r = NULL, kappa_f = NULL, nu_r, nu_f, Theta_r, Theta_f, 
    reference = NULL, custom_colors)
  
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

  #####
  # # Determine which set of items will be returned
  # return_items <- c()
  # if(is.null(delete_items)) { # default: return only the biased items.
  #     return_items <- determine_biased_items(lambda, nu, theta)
  #     return_items <-  setdiff(return_items, which(weights_item == 0))
  # } else {
  #   if (!all(delete_items == floor(delete_items))) {
  #     stop("'delete_items' should only contain integers corresponding
  #            to item indices.")}
  #   if (!all(delete_items < n_i + 1)) {
  #     stop("'delete_items' cannot take integers > the scale length.")}
  #   return_items <- delete_items
  # }
  # return items labels
  return_labels <- paste0("del_i", 1:n_i) # paste0("del_i", return_items)
  #####

  delta_h_str_par_h.foc.del_i <- 
    create_list_of_dfs(ls_len = num_g - 1, ncol = n_i, nrow = 8, 
                       df_cn = return_labels, df_rn = CAIs) 
  acai_p <- acai_s <- h_acai_s_p <- 
    create_list_of_dfs(ls_len = num_g - 1, ncol = 8, nrow = n_i + 1, 
                       df_cn = CAIs, df_rn = c("Full", return_labels))
  h_acai_p <- create_list_of_dfs(ls_len = num_g - 1, ncol = 8, nrow = n_i, 
                                 df_cn = CAIs, df_rn = return_labels) 
  
  store_str <- store_par <- s_p_ref <- s_p_foc <- vector(mode = "list", n_i + 1)
  delta_h_str_par_h.ref.del_i <- as.data.frame(matrix(nrow = 8, ncol = n_i), 
                                               row.names = CAIs) #delta_h_R_Ef <- 
  
  # h_R_Ef <- as.data.frame(matrix(nrow = 8, ncol = n_i + 1)) #h_s_p_ref <- h_s_p_foc <- 
  # delta_h_s_p_acai <- as.data.frame(matrix(nrow = 4, ncol = n_i))
  AI_ratios <- as.data.frame(matrix(ncol = num_g, nrow = n_i + 1))

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
  str_par_h.full <- str_par_h(store_str[[1]]$summary_mi,
                                   store_par[[1]]$summary, num_g = num_g)
  # s_p_ref is a list of length n_i dataframes, for the reference group. each df contains strict, partial, and h columns for a given item set. Here, full.
  s_p_ref[[1]] <- str_par_h.full$Reference
  # s_p_foc is a list of length n_i lists, each element is a list of dfs. the outer list is the item set, the inner list is the focal groups; e.g., s_p_foc[[1]][[1]] == s_p_foc[[1]][["Focal_1"]], which is the CAI and h for the full item set for the Focal 1 group.
  s_p_foc[[1]] <- str_par_h.full[2:num_g]
 
  
  ### h_R_Ef[1] <- apply(
  ###   store_par[[1]]$summary[,(num_g + 1):(2 * num_g - 1)],
  ###   MARGIN = 2, FUN = function(x) cohens_h(store_par[[1]]$summary[,1], x))      #cohens_h(partial$Reference, partial$`E_R(Focal)`)

  # Compute aggregate CAI on the full item set
  temp_acai_p <- get_aggregate_CAI(pmix, store_par[[1]]$summary, inv_cond = "partial")
  acai_p <- Map(function(df, vec) {
    df[1,] <- vec
    df
  }, acai_p, temp_acai_p)
  temp_acai_s <- get_aggregate_CAI(pmix, store_str[[1]]$summary_mi, inv_cond = "strict")
  acai_s <- Map(function(df, vec) {
    df[1,] <- vec
    df
  }, acai_s, temp_acai_s)
  
  # compute cohen's h for the difference between the strict and partial inv. conditions
  # (on the first row of the data frames in each element of the two lists)
  temp_h_acai_s_p <- mapply(function(df_s, df_p) {
    cohens_h(df_s[1, ], df_p[1, ]) 
  }, acai_s, acai_p, SIMPLIFY = FALSE)
  
  # store each vector in the first row of the corresponding dataframe in h_acai_s_p
  h_acai_s_p <- mapply(function(df, temp) {
    df[1, ] <- temp
    df 
  }, h_acai_s_p, temp_h_acai_s_p, SIMPLIFY = FALSE)

  # h: difference between strict and partial invariance for aggregate CAI
   AI_ratios[1,] <- as.vector(c(1, store_par[[1]]$ai_ratio), mode = "double")

  # If no cutoff was provided, set propsel based on PartInv output with all items
  if(is.null(delete_one_cutoff)) {
    propsel_p <- store_par[[1]]$propsel
    cut_z <- NULL
  } else {
    cut_z <- delete_one_cutoff
    propsel_p <- NULL; propsel_s <- NULL
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
    str_par_h.del_i <- str_par_h(store_str[[i]]$summary_mi, 
                                      store_par[[i]]$summary, num_g = num_g)
    s_p_ref[[i]] <- str_par_h.del_i$Reference
    s_p_foc[[i]] <- str_par_h.del_i[2:num_g]

    # delta h: comparing CAI under strict vs. partial invariance when item i is 
    # deleted (i.e. the change in h_s_p_ref and h_s_p_foc) for all groups
    delta_h_str_par_h.ref.del_i[,i - 1] <- delta_h(
      as.data.frame(s_p_ref[[1]])$h, 
      as.data.frame(s_p_ref[[i]])$h)
    
    focal_h.full <- lapply(s_p_foc[[1]], FUN = function(x) as.data.frame(x)$h)
    focal_h.del_i <- lapply(s_p_foc[[i]], FUN = function(x) as.data.frame(x)$h)
    
    for (k in seq_along(focal_h.del_i)) {
      delta_h_str_par_h.foc.del_i[[k]][,return_labels[i - 1]] <- 
        delta_h(focal_h.full[[k]], focal_h.del_i[[k]])
    } 
    
    # h: difference in CAI under partial invariance for the ref group vs. for 
    # the expected CAI for the focal group with the full item set
    ### h_R_Ef[i] <- cohens_h(store_par[[i]]$summary$Reference, 
    ###                       store_par[[i]]$summary$`E_R(Focal)`)

    # delta_h: change in h comparing CAI under partial invariance for the ref
    # group vs. for the expected CAI for the focal group (i.e. the change in
    # h_R_Ef_del) when item i is deleted
    ### delta_h_R_Ef[i-1] <- delta_h(h_R_Ef[1], h_R_Ef[i])

    # Weight the aggregate SR, SE, SP indices under partial and strict invariance
    temp_acai_p <- get_aggregate_CAI(pmix, store_par[[i]]$summary, inv_cond = "partial")
    acai_p <- Map(function(df, vec) {
      df[i, ] <- vec
      df
    }, acai_p, temp_acai_p)
    temp_acai_s <- get_aggregate_CAI(pmix, store_str[[i]]$summary_mi, inv_cond = "strict")
    acai_s <- Map(function(df, vec) {
      df[i, ] <- vec
      df
    }, acai_s, temp_acai_s)
    # compute cohen's h for the difference between the strict and partial inv. conditions
    # (on the i-th row of the data frames in each element of the two lists)
    temp_h_acai_s_p <- mapply(function(df_s, df_p) {
      cohens_h(df_s[i, ], df_p[i, ]) 
    }, acai_s, acai_p, SIMPLIFY = FALSE)
    
    # store each vector in the i-th row of the corresponding dataframe in h_acai_s_p
    h_acai_s_p <- mapply(function(df, temp) {
      df[i, ] <- temp
      df 
    }, h_acai_s_p, temp_h_acai_s_p, SIMPLIFY = FALSE)
    
    # h: change in aggregate CAI when an item is deleted under partial invariance
   temp_h_acai_p <- lapply(acai_p, function(df) cohens_h(df[1, ], df[i, ]))
   h_acai_p <- Map(function(df, vec) {
     df[i - 1, ] <- vec
     df
   }, h_acai_p, temp_h_acai_p)                   
    
    ### delta_h_s_p_acai[i - 1] <- delta_h(h_acai_s_p[1], h_acai_s_p[i])

      AI_ratios[i,] <- as.vector(c(1, store_par[[i]]$ai_ratio), mode = "double")
  }

  colnames(delta_h_str_par_h.ref.del_i) <- return_labels
  rownames(AI_ratios) <- c("Full", return_labels)

  out <- list(
    "AI_ratios" = AI_ratios,
    "ACAI" = acai_p,
    "h_acai_p"= h_acai_p,
    "h_acai_s_p"= h_acai_s_p,
    "delta h (impact of deleting item i)" = 
      list("referenceGroup" = delta_h_str_par_h.ref.del_i,
           "focalGroup(s)" = delta_h_str_par_h.foc.del_i)
    )
  
  ###class(out) <- "itemdeletion"
  return(out)
}

# # Format stored variables
# vars <- list("AI_ratios" = AI_ratios, "h_R_Ef" = h_R_Ef, 
#              "delta_h_str_par_h.ref.del_i" = delta_h_str_par_h.ref.del_i,
#              ###"delta_h_str_par_h.foc.del_i" = delta_h_str_par_h.foc.del_i,
#              "store_str" = store_str, "store_par" = store_par,
#              "s_p_ref" = s_p_ref, "s_p_foc" = s_p_foc, 
#              ###"acai_p" = acai_p, "h_acai_s_p" = h_acai_s_p, 
#              ### "h_acai_p" = h_acai_p, "delta_h_s_p_acai" = delta_h_s_p_acai, 
#              ###  "delta_h_R_Ef" = delta_h_R_Ef, 
#              ###  "return_items" = return_items
#              functioncall = functioncall
# )

###vlist <- format_item_del(n_i, l = vars)

# Declare classes
###  class(vlist$store_par) <- class(vlist$store_str) <-  c("PartInvList", "PartInv")
###  class(vlist$h_s_p_list_ref) <- class(vlist$h_s_p_list_foc) <-  
###    c("PartInvList", "PartInv", "PartInvGroups")

#  vlist <- list("delta_h_str_par_h.ref.del_i" = delta_h_str_par_h.ref.del_i, 
#                "delta_h_str_par_h.foc.del_i" = delta_h_str_par_h.foc.del_i)#TEMP
# list(
#   # "ACAI" = vlist$acai_p,
#   # "h ACAI (deletion)" = vlist$h_acai_p,
#   # "h ACAI SFI-PFI" = vlist$h_acai_s_p,
#   # "delta h ACAI SFI-PFI (deletion)" = vlist$delta_h_s_p_acai,
#   # "AI Ratio" = vlist$AI_ratios,
#   # "h CAI Ref-EF" = vlist$h_R_Ef,
#   # "delta h CAI Ref-EF (deletion)" = vlist$delta_h_R_Ef,
#   "delta h (impact of deleting item i)" = 
#     list("referenceGroup" = delta_h_str_par_h.ref.del_i,
#          "focalGroup(s)" = delta_h_str_par_h.foc.del_i)#,
#   # "h SFI-PFI by groups" = list("reference" = vlist$h_s_p_list_ref, 
#   #                              "focal" = vlist$h_s_p_list_foc),
#   # "PartInv" = list("strict" = vlist$store_str, "partial" = vlist$store_par),
#   # "return_items" = return_items
# )
###class(out) <- "itemdeletion"

#"h_s_p_ref" = h_s_p_ref, #"h_s_p_foc" = h_s_p_foc, 
#delta_h(h_s_p_ref[1], h_s_p_ref[i])
#delta_h(s_p_foc[[1]], x))#delta_h(h_s_p_foc[1], h_s_p_foc[i])
#h_s_p_ref <- h_s_p_foc <- list()#h_s_p_ref[[1]] <- acc$Reference#h_s_p_foc[[1]] <- acc$Focal# h_s_p_ref[1] <- acc$Reference$h; h_s_p_foc[1] <- acc$Focal$h
# h_s_p_ref[i] <- s_p_ref[[i]]$h # h_s_p_foc[i] <- s_p_foc[[i]]$h # h_s_p_ref[[i]] <- s_p_ref[[i]]#  h_s_p_foc[[i]] <- s_p_foc[[i]]
#"h CAI SFI-PFI" = list("ref"= vlist$h_s_p_ref,  #"foc" = vlist$h_s_p_foc),