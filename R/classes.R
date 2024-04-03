#' @importFrom methods setClass

# dividers
stars <-
  "***********************************************************************"
dashes <-
  "-----------------------------------------------------------------------"


summary_print <- function(x, ...) {
  rownames(x) <- c("True Positive", "False Positive", "True Negative",
                        "False Negative",  "Proportion Selected",
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
           itemset = "vector")
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
    h_aggregate_par = "matrix",
    delta_h_str_par_aggregate =  "matrix",
    AI = "matrix",
    deltaREf = "matrix",
    h_R_Ef = "matrix",
    delta_h_str_vs_par = "matrix",
    str_vs_par = "list",
    aggregate = "matrix",
    h_aggregate_str_par = "matrix",
    Ref_foc = "PartInvList",
    PartInv = "PartInvList",
    return_items = "list"
    )
)

#'@export
print.itemdeletion <- function(x, ...) {
  item_set <- x$return_items
    cat(paste0(stars, "\nAGGREGATE CLASSIFICATION ACCURACY INDICES (CAI*)\n",
               stars))
    cat("\nAggregate CAI under PFI computed for item subsets:\n")
    print(round(x$`ACAI`[c(1, item_set + 1), , drop = FALSE], 3))
    cat("\nImpact of deleting an item on aggregate CAI under PFI:\n")
    print(round(x$`h ACAI (deletion)`[c(item_set), , drop = FALSE], 3))
    cat(paste0("\n", stars, "\nAdverse Impact (AI) ratio for item subsets by 
               invariance condition:\n", stars, "\n"))
    print(round(x$AI[c(1, item_set + 1), , drop = FALSE], 3))
    cat(paste0("\n", stars, "\nCOMPARING CAI FOR REFERENCE AND (EXPECTED) FOCAL
               GROUPS\n", stars))
    cat("\nDiscrepancy between CAI of reference vs. Efocal groups under PFI:\n")
    print(round(x$`h CAI Ref-EF`[c(1, item_set + 1), 5:8, drop = FALSE], 3))
    cat(paste0(dashes, "\nImpact of deleting an item on the discrepancy between
               CAI of\nreference vs. Efocal groups under PFI:\n"))
    print(round(x$`delta h CAI Ref-EF (deletion)`[c(item_set), 5:8,
                                                  drop = FALSE], 3))
}


#'@export
summary.itemdeletion <- function(object, ...) {
  item_set <- object$return_items
  cat(paste0(stars, "\nAGGREGATE CLASSIFICATION ACCURACY INDICES (CAI*)\n",
             stars))
  cat("\nAggregate CAI computed for item subsets:\n")
  print(round(object$`ACAI`[c(1, item_set + 1), , drop = FALSE], 3))
  cat("\nImpact of deleting an item on aggregate CAI:\n")
  print(round(object$`h ACAI (deletion)`[c(item_set), , drop = FALSE], 3))
  cat(paste0(dashes, "\nDiscrepancy between aggregate CAI under SFI vs. PFI:
             \n"))
  print(round(object$`h ACAI SFI-PFI`[c(1, item_set + 1), , drop = FALSE], 3))
  cat("\nImpact of deleting an item on the discrepancy between aggregate CAI 
      under \nSFI vs. PFI:\n")
  print(round(object$`delta h ACAI SFI-PFI (deletion)`[c(item_set), ,
                                                       drop = FALSE], 3))
  cat(paste0("\n", stars, "\nAdverse Impact (AI) ratio for item subsets by 
             invariance condition:\n", stars, "\n"))
  print(round(object$AI[c(1, item_set + 1), , drop = FALSE], 3))
  cat(paste0("\n", stars, "\nCOMPARING CAI FOR REFERENCE AND (EXPECTED) FOCAL 
             GROUPS\n", stars))
  cat("\nDiscrepancy between CAI of reference vs. Efocal groups:\n")
  print(round(object$`h CAI Ref-EF`[c(1, item_set + 1), , drop = FALSE], 3))
  cat(paste0(dashes, "\nImpact of deleting an item on the discrepancy between 
             CAI of reference \nvs. Efocal groups:\n"))
  print(round(object$`delta h CAI Ref-EF (deletion)`[c(item_set), ,
                                                     drop = FALSE], 3))
  cat(paste0("\n", stars, "\nCOMPARING CAI UNDER STRICT AND PARTIAL FACTORIAL
             INVARIANCE\n", stars))
  cat("\nDiscrepancy between CAI under SFI vs. PFI for the reference group:\n")
  print(round(object$`h CAI SFI-PFI`[[1]][c(1, item_set + 1), , drop = FALSE],
              3))
  cat("\nDiscrepancy between CAI under SFI vs. PFI for the focal group:\n")
  print(round(object$`h CAI SFI-PFI`[[2]][c(1, item_set + 1), , drop = FALSE],
              3))
  cat(paste0(dashes, "\nImpact of deleting an item on the discrepancy between
             CAI under SFI \nvs. PFI for the reference group:\n"))
  print(round(object$`delta h SFI-PFI (deletion)`[[1]][c(item_set), ,
                                                       drop = FALSE], 3))
  cat("\nImpact of deleting an item on the discrepancy between CAI under SFI 
      \nvs. PFI for the focal group:\n")
  print(round(object$`delta h SFI-PFI (deletion)`[[2]][c(item_set), ,
                                                       drop = FALSE], 3))
  print(object$`h SFI-PFI by groups`$reference)
  print(object$`h SFI-PFI by groups`$focal)
  print(object$PartInv$strict)
  print(object$PartInv$partial)
}
