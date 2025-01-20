create_list_of_dfs <- function(ls_len, ncol, nrow, df_cn, df_rn) {
  lapply(seq_len(ls_len), function(x) {
    df <- data.frame(matrix(NA, ncol = ncol, nrow = nrow))
    colnames(df) <- df_cn; rownames(df) <- df_rn
    df
  })
}

update_rows_in_lists_of_dfs <- function(l1, l2, ind) {
  Map(function(df, df2) {
    df[ind,] <- df2
    df
  }, l1, l2)
}


# Function that performs formatting for various variables to be returned in
# item_deletion_h
format_item_del <- function(n_i, l) {
  # Format stored variables
  names(l$AI_ratios) <- c("full", paste0("|", c(1:n_i)))
  rownames(l$AI_ratios) <- c("AI_SFI", "AI_PFI")
  names(l$h_R_Ef) <-  c("r-Ef", paste0("r-Ef|", c(1:n_i)))
  names(l$delta_s_p_ref) <- names(l$delta_s_p_foc) <- paste0("SFI, PFI|", 
                                                             c(1:n_i))
  names(l$store_str) <- names(l$store_par) <- names(l$s_p_ref_list) <-
    names(l$s_p_foc_list) <- c("full", paste0("|", c(1:n_i)))
  names(l$acai_p) <- c("full", paste0("|", c(1:n_i)))
  names(l$h_acai_p) <- paste0("|", c(1:n_i))
  rownames(l$h_acai_p) <- rownames(l$h_acai_s_p) <-
    c("h(PS*)", "h(SR*)", "h(SE*)", "h(SP*)")
  rownames(l$acai_p) <- c("PS*", "SR*", "SE*", "SP*")
  rownames(l$delta_h_s_p_acai) <-
    paste0("\u0394h(", c("h(PS*)", "h(SR*)", "h(SE*)", "h(SP*)"), ")")
  rownames(l$delta_s_p_ref) <- rownames(l$delta_s_p_foc) <- 
    rownames(l$delta_h_R_Ef) <-
    paste0("\u0394h(", c("TP", "FP", "TN", "FN", "PS", "SR", "SE", "SP"), ")")
  rownames(l$h_R_Ef) <- rownames(l$h_s_p_ref) <- rownames(l$h_s_p_foc) <-
    c("h(TP)", "h(FP)", "h(TN)", "h(FN)", "h(PS)", "h(SR)", "h(SE)", "h(SP)")
  
  store_par <- list(outputlist = l$store_par, condition = "partial",
                    itemset = l$return_items)
  store_str <- list(outputlist = l$store_str, condition = "strict",
                    itemset = l$return_items)
  h_s_p_list_ref <- list(outputlist = l$s_p_ref_list, condition = "ref",
                                itemset = l$return_items)
  h_s_p_list_foc <- list(outputlist = l$s_p_foc_list, condition = "foc",
                                itemset = l$return_items)
  names(l$h_s_p_ref) <- names(l$h_s_p_foc) <- c("SFI, PFI", 
                                            paste0("SFI, PFI|", c(1:N)))
  
  names(l$delta_h_R_Ef) <- paste0("r-Ef|", c(1:N))
  names(l$delta_h_s_p_acai) <- paste0("SFI, PFI|", c(1:N))
  names(l$h_acai_s_p) <- c("SFI, PFI", paste0("SFI, PFI|", c(1:N)))
  
  acai_p <- as.data.frame(cbind(l$acai_p))
  h_acai_p <- as.data.frame(l$h_acai_p)
  
  return(list("acai_p" = t(acai_p), "h_acai_p" = t(h_acai_p), 
              "h_acai_s_p" = t(l$h_acai_s_p), 
              "delta_h_s_p_acai" = t(l$delta_h_s_p_acai), 
              "AI_ratios" = t(l$AI_ratios), "h_R_Ef" = t(l$h_R_Ef), 
              "delta_h_R_Ef" = t(l$delta_h_R_Ef), "h_s_p_ref" = t(l$h_s_p_ref), 
              "h_s_p_foc" = t(l$h_s_p_foc), 
              "delta_s_p_ref" = t(l$delta_s_p_ref), 
              "delta_s_p_foc" = t(l$delta_s_p_foc), 
              "h_s_p_list_ref" = h_s_p_list_ref, 
              "h_s_p_list_foc" = h_s_p_list_foc, 
              "store_str" = store_str, "store_par" = store_par))
}



#' @title
#' Compute PS, SR, SE, SP weighted by group proportions
#'
#' @name
#' get_aggregate_CAI
#'
#' @description
#' \code{get_aggregate_CAI} computes aggregate PS, SR, SE, SP under partial or
#' strict invariance by weighting the TP, TF, TN, FP values for the reference
#' and focal groups with the group proportions.
#'
#' @param pmix Proportion of the reference group.
#' @param store_summary The summary table from [PartInv()]
#' under partial or strict invariance.
#' @param inv_cond Strict vs. partial.
#'
#' @return A vector of length 4.
#'          \item{PS}{Proportion selected, computed as \eqn{TP + FP}.}
#'          \item{SR}{Success ratio, computed as \eqn{TP/(TP + FP)}.}
#'          \item{SE}{Sensitivity, computed as \eqn{TP/(TP + FN)}.}
#'          \item{SP}{Specificity, computed as \eqn{TN/(TN + FP)}.}
#' @export
get_aggregate_CAI <- function(pmix, store_summary, inv_cond) {
  # Ensure pmix sums to 1
  if (abs(sum(pmix) - 1) > 1e-6) {
    stop("The sum of pmix must be equal to 1.")
  }
  
  num_g <- length(pmix)
  # Check for consistency between pmix and store_summary
  if ((inv_cond == "partial") && ((2 * length(pmix) - 1) != ncol(store_summary))) {
    stop("Length of pmix must match the number of groups in store_summary.")
  }  
  if ((inv_cond == "strict") && (length(pmix) != ncol(store_summary))) {
    stop("Length of pmix must match the number of groups in store_summary.")
  }
  # Compute weighted aggregates for TP, FP, TN, FN 
  # TP <- sum(pmix * store_summary[1, seq_len(num_g)]) # this gets the aggregate across groups, it should be aggregates between the reference and one focal group
  # FP <- sum(pmix * store_summary[2, seq_len(num_g)])
  # TN <- sum(pmix * store_summary[3, seq_len(num_g)])
  # FN <- sum(pmix * store_summary[4, seq_len(num_g)])
  # 
  
  store <- store_summary[, seq_len(num_g)]
  
  weighted_pair_sum <- function(vec) {
    as.numeric(
      sapply(2:num_g, function(i) vec[1] * pmix[1] + vec[i] * pmix[i])
    )
  }

  TP <- weighted_pair_sum(store[1,])
  FP <- weighted_pair_sum(store[2,])
  TN <- weighted_pair_sum(store[3,])
  FN <- weighted_pair_sum(store[4,])
  
  # Compute metrics
  PS <- TP + FP
  SR <- TP / (TP + FP)
  SE <- TP / (TP + FN)
  SP <- TN / (TN + FP)
  
  df <- rbind(TP, FP, TN, FN, PS, SR, SE, SP)
  ls <- split(as.matrix(df), col(df))
  return(ls)
}


#' @title
#' Check for misleading improvements in aggregate CAI
#'
#' @name
#' err_improv_acai
#'
#' @description
#' \code{err_improv_acai} checks if any improvement observed in aggregate CAI
#' may have resulted from the higher mixing proportion of the reference group
#' masking worsening performance for the focal group. If the effect size of any
#' change indicating worse performance for the focal group and better
#' performance for the reference group is larger than 0.1, prints a warning.
#' @param i Index of item under consideration.
#' @param s_full PartInv summary for the case where all items are retained.
#' @param s_del1 PartInv summary for the case where item i is excluded.
#' @param num_g Number of groups
err_improv_acai <- function(i, s_full, s_del1, num_g) {
  # store the focal groups' values in a list
  focals <- list(s_full[,2:num_g], s_del1[,2:num_g])
  
  # compute Cohen's h for the difference between full and drop one indices for 
  # the reference group
  h_r <- cohens_h(s_full[,1], s_del1[,1]) #cohens_h(r, r_del1)
  # compute Cohen's h for the difference between full and drop one indices for 
  # the focal groups
  h_f <- lapply(focals, FUN = function(x) cohens_h(x[[1]], x[[2]])) #h_f <- cohens_h(f_full, f_del1)
  # Check for changes (boolean)
  r_bool <- s_full[,1] < s_del1[,1] 
  f_bool_leq1 <- s_full[,2:num_g] <= s_del1[2:num_g]
  f_bool_leq <- apply(f_bool_leq1, MARGIN = 1, FUN = all)
  f_bool_geq1 <- s_full[,2:num_g] >= s_del1[2:num_g]
  f_bool_geq <- apply(f_bool_leq1, MARGIN = 1, FUN = all)
  
  # wrangle into a matrix to apply operations by row
  h_f_mat <- matrix(unlist(h_f), nrow = 8, ncol = (num_g - 1))
  # check the difference for the reference or focal groups has Cohen's h > 0.1
  h_rf.1 <- (h_r > 0.1 | apply((h_f_mat > 0.1), MARGIN = 1, FUN = all))

  vals <- c("TP", "FP", "TN", "FN")
  cat1 <- function(i, vals, val_i) {
    cat("Increases in aggregate CAI after deleting item ", i, "may be
         misleading due to the \n mixing proportion. Examine ", vals[val_i],
        "values from detailed output tables before proceeding.\n")
  }
  # TP_f decreases/remains unchanged & TP_r increases
  if(r_bool[1] && f_bool_geq[1] && h_rf.1[1]) cat1(1, vals, 1)
  # FP_r decreases and FP_f increases/remains unchanged
  if(!r_bool[2] && f_bool_leq[2] && h_rf.1[2]) cat1(2, vals, 2)
  # TN_f decreases/remains unchanged and TN_r increases
  if(r_bool[3] && f_bool_geq[3] && h_rf.1[3]) cat1(3, vals, 3)
    # FN_r decreases and FN_f increases/remains unchanged
  if(!r_bool[4] && f_bool_leq[4] && h_rf.1[4]) cat1(4, vals, 4)


  }

#' @title
#' Delete item i and redistribute its weight within subscale
#'
#' @name
#' redistribute_weights
#'
#' @description
#' \code{redistribute_weights} replaces the item weight with 0 for the item to
#' be deleted, and redistributes this item's weight across the remaining items.

#' @param weights_item A vector of item weights.
#' @param n_dim Number of dimensions, 1 by default. If the user does not supply
#'        a value, assumes that the scale is unidimensional.
#' @param n_i_per_dim A vector containing the number of items in each
#'        dimension; `NULL` by default. If the user provides a value for n_dim
#'        that is \eqn{> 1} but leaves \code{n_i_per_dim = NULL}, assumes that
#'        the subscales have an equal number of items.
#' @param del_i Index of the item to be deleted.
#'
#' @return `new_w` Weights vector with redistributed weights.
#' @examples
#' one_dim_w <- c(1:7)
#' redistribute_weights(one_dim_w, del_i = 2)
#' redistribute_weights(one_dim_w, n_dim = 1, n_i_per_dim = 7, del_i = 2)
#' sum(one_dim_w)==sum(redistribute_weights(one_dim_w, del_i = 2))
#'
#' multi_eq_w <- c(1:9)
#' redistribute_weights(multi_eq_w, n_dim = 3, del_i = 2)
#' redistribute_weights(multi_eq_w, n_dim = 3, n_i_per_dim = c(3, 3, 3), 
#' del_i = 2)
#' sum(multi_eq_w)==sum(redistribute_weights(multi_eq_w, n_dim = 3, del_i = 2))
#'
#' multi_uneq_w <- c(1:12)
#' redistribute_weights(multi_uneq_w, n_dim = 3, n_i_per_dim = c(3, 6, 3), 
#' del_i=2)
#' sum(multi_uneq_w)==sum(redistribute_weights(multi_uneq_w, n_dim = 3,
#'                                             n_i_per_dim = c(3, 6, 3),
#'                                             del_i=2))
#' @export
redistribute_weights <- function(weights_item, n_dim = 1, n_i_per_dim = NULL,
                               del_i) {
  n_items <- length(weights_item)
  new_w <- weights_item
  new_w[del_i] <- 0
  del_weight <- weights_item[del_i] # the weight to be redistributed

  # Unidimensional
  if ((n_dim == 1) && (is.null(n_i_per_dim) || length(n_i_per_dim) == 1)) {
    # Increase each non-zero item in the vector by the weight to be distributed
    # proportional to the original weighting of the items.
    new_w[new_w != 0] <- new_w[new_w != 0] +
      new_w[new_w != 0] * del_weight / sum(new_w[new_w != 0])

  # Multidimensional, equal n
  } else if ((n_dim > 1) && is.null(n_i_per_dim) && n_items %% n_dim != 0) {
    stop("Please pass a vector of subscale lengths to n_i_per_dim.")

  # Multidimensional, number of items per dimension is not specified
  } else if ((n_dim > 1) && is.null(n_i_per_dim)) {
    # Split indices into dimensions assuming dimensions have the same length.
    i_by_dim <- split(1:n_items, cut(seq_along(1:n_items), n_dim, 
                                     labels = FALSE))
    new_w <- multidim_redist(n_dim, del_i, i_by_dim, new_w, del_weight)

  # Multidimensional, unequal n
  } else if ((n_dim > 1) && !is.null(n_i_per_dim)) {
    # Split indices into dimensions
    i_by_dim <- split(1:n_items, cut(seq_along(1:n_items),
                                     breaks = cumsum(c(0, n_i_per_dim)),
                                     labels = FALSE))
    new_w <- multidim_redist(n_dim, del_i, i_by_dim, new_w, del_weight)
  } else {
    stop("Check n_dim and n_i_per_dim")
  }
  return(new_w)
}
# Helper function for redistribute_weights() that, for each dimension of the
# scale, updates the specific subscale with redistributed weights proportional
# to the original weights if the item to be deleted is in that dimension.
multidim_redist <- function(n_dim, del_i, i_by_dim, new_w, del_weight) {
  for (k in 1:n_dim) {
    if (del_i %in% i_by_dim[[k]]) { # If del_i is in dimension k
      # Create temporary vector to store the remaining indices in dimension k
      temp_i <- i_by_dim[[k]][i_by_dim[[k]] != del_i]
      new_w[temp_i][new_w[temp_i] != 0] <- new_w[temp_i][new_w[temp_i] != 0] +
        new_w[temp_i][new_w[temp_i] != 0] * del_weight / sum(new_w[temp_i])
    }
  }
  return(new_w)
}

#' @title
#' Compute Cohen's h effect size for the difference in two proportions.
#'
#' @name
#' cohens_h
#'
#' @description
#' \code{cohens_h} computes Cohen's h (Cohen, 1988) for the difference in two
#' proportions using \eqn{h = 2arcsin(\sqrt{p1}) - 2arcsin(\sqrt{p2})}.
#'
#' @param p1 The first proportion.
#' @param p2 The second proportion.
#' @return `h` The computed Cohen's h value.
#' @examples
#' cohens_h(0.7, 0.75)
#' cohens_h(0.3, 0.4)
#' @export
cohens_h <- function(p1, p2) {
  h <- 2 * asin(sqrt(p1)) - 2 * asin(sqrt(p2))
  return(h)
}

#' @title
#' Compute effect size for the impact of item deletion
#'
#' @name
#' delta_h
#'
#' @description
#' \code{delta_h} Computes the effect size of the impact of item bias by
#' comparing Cohen's h values for CAI for the full versus delete-one item sets.
#'
#' @param h_R h effect sizes for when the item is included.
#' @param h_i_del h effect sizes for when the item is deleted.
#' @return Cohen's h for the difference in the classification accuracy index
#' when the item is deleted.
#' @examples
#' delta_h(0.04, 0.01)
#' delta_h(-0.002, 0.011)
#' @export
delta_h <- function(h_R, h_i_del) {
  abs(h_R) - abs(h_i_del)
}

#' @title
#' Determine biased items
#'
#' @name
#' determine_biased_items
#'
#' @description
#' \code{determine_biased_items} takes in the factor loadings, intercepts, and
#'  uniqueness, and returns indices of noninvariant items.
#' @param nu_r,nu_f,Theta_r,Theta_f,lambda_r,lambda_f Deprecated; included only 
#' for backward compatibility.
#' @param lambda Factor loadings.
#' @param nu Measurement intercepts.
#' @param theta Uniqueness.
#' @return A vector containing the indices of the biased items.
#' @examples
#' lambda_matrix <- matrix(0, nrow = 5, ncol = 2)
#' lambda_matrix[1:2, 1] <- c(.322, .655)
#' lambda_matrix[3:5, 2] <- c(.398, .745, .543)
#' lambda_matrix2 <- lambda_matrix
#' lambda_matrix2[3,1] <- 3
#' determine_biased_items(lambda = list(lambda_matrix, lambda_matrix2),
#'                        nu = list(c(.225, .025, .010, .240, .125),
#'                                  c(.225, -.05, .240, -.025, .125)),
#'                        theta = list(diag(1, 5), diag(c(1, .95, .80, .75, 1))))
#' @export
determine_biased_items <- function(lambda, nu, theta, 
                                   lambda_r = NULL, lambda_f = lambda_r,
                                   nu_r = NULL, nu_f = nu_r,
                                   Theta_r = NULL, Theta_f = Theta_r) {
  biased_items <- c()
  
  # backward compatibility
  if (missing(nu) && !is.null(nu_r)) {
    nu <- vector(2, mode = "list")
    nu[[1]] <- nu_r; nu[[2]] <- nu_f
  }
  if (missing(lambda) && !is.null(lambda_r)) {
    lambda <- vector(2, mode = "list")
    lambda[[1]] <- lambda_r; lambda[[2]] <- lambda_f
  }
  if (missing(theta) && !is.null(Theta_r)) {
    theta <- vector(2, mode = "list")
    theta[[1]] <- Theta_r; theta[[2]] <- Theta_f
  }
  
  # compare factor loadings (lambda)
  lambda_mismatch <- find_mismatched_indices(lambda)
  if (!is.null(lambda_mismatch)) {
    biased_items <- c(biased_items, unique(lambda_mismatch[, 1]))  # row indices
  }
  # compare uniqueness (theta)
  theta_mismatch <- find_mismatched_indices(theta)
  if (!is.null(theta_mismatch)) {
    biased_items <- c(biased_items, unique(theta_mismatch[, 1]))  # row indices
  }
  # compare intercepts (nu)
  nu_mismatch <- find_mismatched_indices(nu)
  if (!is.null(nu_mismatch)) {
    biased_items <- c(biased_items, unique(nu_mismatch))  # Vector indices
  }
  
  biased <- unique(biased_items)
  if (length(biased) == 0) {
    message("Strict invariance holds for all items.")
    return(NULL)
  }
  return(sort(biased))
}




find_mismatched_indices <- function(lst) {
  if (any(sapply(lst, length) == 0)) stop("The list contains empty elements.")
  # check for non-numeric elements
  if (!all(sapply(lst, function(x) is.numeric(x) || is.matrix(x)))) {
    stop("All elements must be numeric.")
  }
  if (all(sapply(lst, is.vector))) {
    if (!all(sapply(lst, length) == length(lst[[1]]))) {
      stop("Vectors have unequal lengths.")
    }
    
    # combine vectors into a matrix
    combined_matrix <- do.call(rbind, lst)
    # check for mismatches across rows for each column (vector element index)
    mismatches <- apply(combined_matrix, 2, function(x) length(unique(x)) > 1)

    if (!any(mismatches)) return(NULL)
    
    return(which(mismatches))
  } else {
    # handle lists of matrices or mixed inputs
    lst <- lapply(lst, function(x) if (is.vector(x)) matrix(x, nrow = 1) else x)
    
    # ensure all elements have identical dimensions
    dims <- sapply(lst, dim)
    if (!all(apply(dims, 1, function(x) length(unique(x)) == 1))) {
      stop("All elements must have the same dimensions.")
    }
    
    # convert list elements into arrays for element-wise comparison
    combined_array <- array(unlist(lst), dim = c(dim(lst[[1]]), length(lst)))
    
    # check for mismatches across the third dimension
    mismatches <- apply(combined_array, c(1, 2), function(x) length(unique(x)) > 1)
    
    if (!any(mismatches)) {
      return(NULL)
    }
    # extract row and column indices of mismatches
    return(which(mismatches, arr.ind = TRUE))
  }
}




