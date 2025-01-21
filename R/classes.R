#' @importFrom methods setClass

# classes: PartInv, PartInvGroups, PartInvList, itemdeletion


# dividers
stars <-
  "***********************************************************************"
dashes <-
  "-----------------------------------------------------------------------"

summary_print <- function(x, ...) {
  rownames(x) <- c("True Positive", "False Positive", "True Negative",
                   "False Negative", "Proportion Selected",
                   "Success Ratio", "Sensitivity", "Specificity")
  print(round(x, digits = 3))
}

setClass("PartInvGroups", representation(tab = "data.frame"))

#'@export
print.PartInvGroups <- function(x, ...) {
  rownames(x) <- c("True Positive", "False Positive", "True Negative",
                        "False Negative",  "Proportion Selected",
                        "Success Ratio", "Sensitivity", "Specificity")
  print(round(x, 3))
}

setClass("PartInv",
  representation(
    propsel = "numeric",
    cutpt_xi = "numeric", cutpt_z = "numeric",
    summary = "data.frame",
    bivar_data = "list",
    ai_ratio = "numeric",
    propsel_mi = "numeric",
    cutpt_xi_mi = "numeric", cutpt_z_mi = "numeric",
    bivar_data_mi = "list",
    summary_mi = "data.frame",
    labels = "character", 
    functioncall = "character"
  )
)
#' @method print PartInv 
#' @title Print method for PartInv class
#' @description performs printing and formatting on an PartInv object.
#' @param x An object of class \code{PartInv}, the output from
#'  \code{PartInv()}.
#' @param ... Additional arguments passed to methods.
#' @return NULL
#' 
#' @examples
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
#' out <- PartInv(propsel = 0.25, pmix = c(1/4, 1/4, 1/4, 1/4),
#'   alpha = list(0, -0.70, -1.05, -1.10), psi = list(1, 1.2, 1.29, 1.3),
#'   nu = list(nu_matrix, nu_matrix1, nu_matrix2, nu_matrix3),
#'   lambda = list(lambda_matrix, lambda_matrix, lambda_matrix, lambda_matrix),
#'   theta = list(theta_matrix, theta_matrix1, theta_matrix2, theta_matrix3),
#'   plot_contour = TRUE, show_mi_result = TRUE,
#'   labels = c("Group 1", "Group 2", "Group 3", "Group 4"),
#'   custom_colors = c("salmon1", "lightgreen", "skyblue1", "pink"))
#'
#'@export
print.PartInv <- function(x, ...) {
  cat("Partial invariance results:\n\n")
  cat("Proportion selected: ", round(x$propsel, 3), "\n")
  cat("Cutpoint on the latent scale (xi): ", round(x$cutpt_xi, 3), "\n")
  cat("Cutpoint on the observed scale (Z): ", round(x$cutpt_z, 3), "\n")
  cat(paste0("Adverse impact ratio ", "(reference group: '",
             colnames(x$summary)[1], "'):\n"))
  print(as.data.frame(lapply(x$ai_ratio, round, digits = 3), row.names = ""))
  cat("\n")
  nc <- ncol(x$summary)
  if (nc > 8) {
    cat("Classification Accuracy Indices:\n")
    summary_print(x$summary[, 1:(ceiling(nc / 2))])
    cat("\n")
    cat("Expected Results if Latent Distributions Matched the Reference Group:\n")
    summary_print(x$summary[, (ceiling(nc / 2) + 2):nc])
  } else {
    cat("Classification Accuracy Indices:\n")
    summary_print(x$summary)
  }
  if (!is.null(x$summary_mi)) {
    cat("\n\nStrict invariance results:\n\n")
    cat("Proportion selected: ", round(x$propsel_mi, 3), "\n")
    cat("Cutpoint on the latent scale (xi): ", round(x$cutpt_xi_mi, 3), "\n")
    cat("Cutpoint on the observed scale (Z): ", round(x$cutpt_z_mi, 3), "\n\n")
    cat("Classification Accuracy Indices:\n")
    summary_print(x$summary_mi)
  }
}

setClass("PartInvList",
  representation(
    outputlist = "list",
    condition = "character",
    itemset = "vector"
    )
)

#'@export
print.PartInvList <- function(x, ...) {
  if (x$condition == "partial") {
    cat(paste0("\n", stars))
    cat("\nPartInv() outputs under Partial Factorial Invariance (PFI)\n")
    cat(stars)
    cat("\n\nUnder PFI, full item set:\n\n")
    print.PartInv(x$outputlist[[1]])
    for (i in x$itemset) {
      cat(dashes)
      cat("\n\nUnder PFI, if item", i, "is dropped:\n\n")
      print.PartInv(x$outputlist[[i + 1]])
    }
  }
    if (x$condition == "strict") {
    cat(paste0("\n", stars))
    cat("\nPartInv() outputs under Strict Factorial Invariance (SFI)\n")
    cat(stars)
    cat("\n\nUnder SFI, full item set:\n\n")
    print.PartInv(x$outputlist[[1]])
    for (i in x$itemset) {
      cat(dashes)
      cat("\n\nUnder SFI, if item", i, "is dropped:\n\n")
      print.PartInv(x$outputlist[[i + 1]])
    }
    }
  if (x$condition == "ref") {
    cat(paste0("\n", stars))
    cat("\nCAI under SFI vs. PFI for the reference group\n")
    cat(stars)
    cat("\n\nReference group, full item set:\n")
    print.PartInvGroups(x$outputlist[[1]])
    for (i in x$itemset) {
      cat(dashes)
      cat("\n\nReference group, if item", i, "is dropped:\n")
      print.PartInvGroups(x$outputlist[[i + 1]])
    }
  }
  if (x$condition == "foc") {
    cat(paste0("\n", stars))
    cat("\nCAI under SFI vs. PFI for the focal group\n")
    cat(stars)
    cat("\n\nFocal group, full item set:\n")
    print.PartInvGroups(x$outputlist[[1]])
    for (i in x$itemset) {
      cat(dashes)
      cat("\n\nFocal group, if item", i, "is dropped:\n")
      print.PartInvGroups(x$outputlist[[i + 1]])
    }
  }
}

setClass("itemdeletion",
  representation(
    AI_ratios = "data.frame",
    ACAI = "list",
    h_acai_p = "list",
    h_acai_s_p = "list",
    delta_h_acai_s_p =  "list",
    h_R_Ef = "list",
    delta_h_R_Ef = "list",
    delta_h_str_vs_par = "list",
    #str_vs_par = "list",
   # aggregate = "matrix",
    #h_aggregate_str_par = "matrix",
    #Ref_foc = "PartInvList",
   # PartInv = "PartInvList",
    items = "vector",
    function_call = "call",
    digits = "numeric"
    )
)

#' @method print itemdeletion 
#' @title Print method for itemdeletion class
#' @description performs printing and formatting on an itemdeletion object.
#' @param x An object of class \code{itemdeletion}, the output from
#'  \code{item_deletion_h()}.
#' @param ... Additional arguments passed to methods.
#' @return NULL
#' 
#' @examples
#' # Multidimensional example
#' lambda_matrix <- matrix(0, nrow = 5, ncol = 2)
#' lambda_matrix[1:2, 1] <- c(.322, .655)
#' lambda_matrix[3:5, 2] <- c(.398, .745, .543)
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
#'@export
print.itemdeletion <- function(x, ...) {
  digits <- x$digits
  item_set <- x$items
  # apply rounding and get the composite index columns
  x <- c(x[1],
              lapply(x[-c(1, 10:12)], function(inner_list) {
                lapply(inner_list, function(df) {
                  df <- df[, 5:8]
                  round(df, digits)
                  })
                }),
            list(x[10:12]))
  x$AI <- round(x$AI, digits)
  cat(paste0("\n", stars, "\nAdverse Impact (AI) ratios under partial invariance by group\n", stars, "\n"))
  print(x$AI[item_set, -1])
  cat("\n   (Note: AI ratios equal 1 under strict invariance by definition.)\n")
  cat(paste0(stars, "\nAGGREGATE CLASSIFICATION ACCURACY INDICES (CAI*)\n", stars))
  cat("\nAggregate CAI (CAI*) under partial invariance:\n")
  print_dfs_from_list(x$ACAI, c("Full", item_set) )
  cat("\nImpact of deleting an item on aggregate CAI (CAI*) under partial invariance:\n")
  print_dfs_from_list(x$h_acai_p, item_set)
  cat(paste0("\n", stars, "\nCOMPARING CAI FOR REFERENCE AND (EXPECTED) FOCAL GROUPS\n", stars))
  cat("\nDiscrepancy between CAI of reference vs. Efocal groups under PFI:\n")
  print_dfs_from_list(x$h_R_Ef, c("r_Ef", paste0("r_Ef", item_set))) 
  cat(paste0(dashes, "\nImpact of deleting an item on the discrepancy between observed CAI for the 
reference group and expected CAI for the focal groups (Efocal):\n"))
   print_dfs_from_list(x$delta_h_R_Ef, item_set)
   invisible(NULL)
}

# register the custom print method for the itemdeletion class
setMethod("print", "itemdeletion", print.itemdeletion)

#' @method summary itemdeletion 
#' @title Summary method for itemdeletion class
#' @description prints a detailed summary for an itemdeletion object.
#' @param object An object of class \code{itemdeletion}, the output from
#'  \code{item_deletion_h()}.
#' @param ... Additional arguments passed to methods.
#' @return NULL
#' 
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
#' summary(multi_dim)
#'@export
summary.itemdeletion <- function(object, ...) {
  item_set <- object$items
  digits <- object$digits
  # apply rounding and get the composite index columns
  object <- c(object[1],
         lapply(object[-c(1, 10:12)], function(inner_list) {
           lapply(inner_list, function(df) {
             df <- df[, 5:8]
             round(df, digits)
           })
         }),
         list(object[10:12]))
 
  cat(paste0("\n", stars, "\nAdverse Impact (AI) ratios for by invariance condition and group:\n", stars, "\n"))
  object$AI <- round(object$AI, digits)
  print(object$AI[item_set,])
  
  cat("\nAggregate CAI under PFI computed for item subsets:\n")
  lapply(object$ACAI, FUN = function(y) print(y[c("Full", item_set),]))
  cat("\nImpact of deleting an item on aggregate CAI under PFI:\n")
  lapply(object$h_acai_p, FUN = function(y) print(y[item_set,]))
  cat(paste0(dashes, "\nDiscrepancy between aggregate CAI under SFI vs. PFI: \n"))
  lapply(object$h_acai_s_p, FUN = function(y) print(y[c("Full", item_set),]))
  cat("\nImpact of deleting an item on the discrepancy between ACAI under SFI vs. PFI:\n")
  lapply(object$delta_h_acai_s_p, FUN = function(y) print(y[item_set,]))
  cat(paste0("\n", stars, "\nCOMPARING CAI FOR REFERENCE AND (EXPECTED) FOCAL GROUPS\n", stars))
  cat("\nDiscrepancy between CAI of reference vs. Efocal groups under PFI:\n")
  lapply(object$h_R_Ef, FUN = function(y) print(y[c("r_Ef", paste0("r_Ef", item_set)),]))
  cat(paste0(dashes, "\nImpact of deleting an item on the discrepancy between observed CAI for the 
reference group and expected CAI for the focal groups (Efocal):\n"))
  lapply(object$delta_h_R_Ef, FUN = function(y) print(y[item_set,]))
  cat("\nDiscrepancy between CAI under SFI vs. PFI:\n")
  lapply(object$h_s_p, FUN = function(y) print(y[c("Full", item_set),]))
  cat(paste0(dashes, "\nImpact of deleting an item on the discrepancy between
             CAI under SFI \nvs. PFI for the reference group:\n"))
  lapply(object$delta_h_s_p, FUN = function(y) print(y[item_set,]))
  
  
  # cat("\nImpact of deleting an item on the discrepancy between CAI under SFI 
  #     \nvs. PFI for the focal group:\n")
  # print(round(object$`delta h SFI-PFI (deletion)`[[2]][c(item_set), ,
  #                                                      drop = FALSE], 3))
  # print(object$`h SFI-PFI by groups`$reference)
  # print(object$`h SFI-PFI by groups`$focal)
  # print(object$PartInv$strict)
  # print(object$PartInv$partial)
}
# register the custom summary method for the itemdeletion class
setMethod("summary", "itemdeletion", summary.itemdeletion)

