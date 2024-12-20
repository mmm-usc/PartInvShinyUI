#' @importFrom graphics lines
NULL

#' Plot classification accuracy indices at different proportions of selection
#'   or at different threshold (cutoff) values
#' 
#' \code{plotPropselRange} plots classification accuracy indices at different
#' proportions of selection under partial and strict invariance conditions for 
#' a given CFA fit.
#' 
#' @param cfa_fit CFA model output from lavaan.
#' @param pmix List of length `g` containing the mixing proportions of each
#'     group (where `g` is the number of groups). If `NULL`, defaults to `1/g` 
#'     for each group (i.e., the populations have equal size).
#' @param labels A character vector with `g` elements to label the reference
#'     and focal groups on the plot, where `g` is the number of groups. If not
#'     provided, groups are labeled automatically as 'Reference' (for the first
#'     group) and 'Focal_1' through 'Focal_(g-1)', where `g` is the number of
#'     groups.
#' @param cai_names A vector of strings indicating the classification accuracy
#'     indices of interest. c("PS", "SR", "SE", "SP", "AI) by default.
#' @param mod_names A vector of strings indicating the invariance conditions of
#'     interest. c("par", "str") by default.
#' @param from The lowest proportion of selection to consider. `0.01` by default.
#' @param to The largest proportion of selection to consider. `0.25` by default.
#' @param by The increment of the sequence of proportions. `0.01` by default.
#' @param cutoffs_from The lowest threshold to consider.`NULL` by default.
#' @param cutoffs_to The largest threshold to consider. `NULL` by default.
#' @return Eight plots illustrating how proportion selected (PS), success ratio 
#'     (SR), sensitivity (SE), and specificity (SP) change across different 
#'     proportions of selection under partial and strict invariance conditions.
#' @param custom_colors Optional argument for specifying colors.
#' @param ... Additional arguments.
#' @examples
#' \dontrun{
#' library(lavaan)
#' HS <- HolzingerSwineford1939
#' HS$sex <- as.factor(HS$sex)
#' HS.model <- ' visual  =~ x1 + x2 + x3
#'               textual =~ x4 + x5 + x6
#'               speed   =~ x7 + x8 + x9 '
#' 
#' fit <- cfa(HS.model, data = HS, group = "sex")
#' 
#' # plot CAI at different proportions of selection
#' plotPropselRange(fit, pmix = table(HS$sex)/sum(table(HS$sex)))
#' 
#' # plot CAI at different cutoffs
#' plotPropselRange(fit, pmix = table(HS$sex)/sum(table(HS$sex)), 
#'    cutoffs_from = 35, cutoffs_to = 50)
#'
#' # plot only SR under partial invariance for up to 10% selection.
#' plotPropselRange(fit, pmix = table(HS$sex)/sum(table(HS$sex)), 
#'     from = 0.01, to = 0.10, cai_names = "SR", mod_names = "par")
#' }
#' @export
plotPropselRange <- function(cfa_fit,
                             pmix = NULL,
                             labels = NULL,
                             cai_names = c("PS", "SR", "SE", "SP", "AI"),
                             mod_names = c("par", "str"),
                             from = 0.01,
                             to = 0.25,
                             by = 0.01,
                             cutoffs_from = NULL,
                             cutoffs_to = NULL, 
                             custom_colors = NULL, ...
                             ) {
  
  stopifnot("cai_names can only take the following values: PS, SR, SE, SP, AI." =
              (all(cai_names %in% c("PS", "SR", "SE", "SP", "AI"))))
  stopifnot("mod_names can only take the following values: par, str" =
              (all(mod_names %in% c("par", "str"))))
  plotAIs <- TRUE
  if ("AI" %in% cai_names) {
    
    cai_names <- cai_names[!cai_names %in% "AI"]
    if (length(cai_names) == 0) {
      cai_names <- NULL
    }
  } else {
    plotAIs <- FALSE
  }
  
  est <- format_cfa_partinv(cfa_fit, comp = "est")
  
  propsels <- seq(from = from, to = to, by = by)
  use <- "propsels"
  xl <- "Proportion of selection"
  rangeVals <- propsels
  
  if ((is.null(cutoffs_from) && !is.null(cutoffs_to)) ||
    (!is.null(cutoffs_from) && is.null(cutoffs_to))) {
    warning("If you would like to plot CAI at different thresholds, provide
            parameter values for both `cutoffs_to` and `cutoffs_from`.
            CAI were plotted at different proportions of selection by default.")
  }
    
  # if the user provided the max and min cutoff values, update rangeVals with 
  # a range of cutoffs
  if (!is.null(cutoffs_from) && !is.null(cutoffs_to)) {
    cutoffs <- seq(from = cutoffs_from, to = cutoffs_to, by = by)
    rangeVals <- cutoffs
    xl <- "Thresholds" # for the plots later
    use <- "cutoffs"
  }
  n_g <- cfa_fit@Data@ngroups # number of groups
  
  # if the user did not provide labels, or provided the wrong number of labels,
  if (is.null(labels) || (length(labels) != n_g)) {
    labels <- cfa_fit@Data@group.label
    labels <- paste(labels, c("(reference)", rep("(focal)", n_g - 1)))
  }
  
  ls_mat <- matrix(NA, ncol = length(rangeVals), nrow = n_g,
                   dimnames = list(labels, rangeVals))
  AIs <- matrix(NA, ncol = length(rangeVals), nrow = n_g - 1,
                dimnames = list(labels[-1], rangeVals))
  ls_names <- c(t(outer(cai_names, Y = mod_names, FUN = paste, sep = "_")))
  ls <- rep(list(ls_mat), length(ls_names))
  names(ls) <- ls_names
  
  # if pmix is missing, assume equal mixing proportions
  if (is.null(pmix)) pmix <- as.matrix(c(rep(1 / n_g, n_g)), ncol = n_g)
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
                        psi = est$psi,
                        lambda = est$lambda,
                        theta = est$theta,
                        alpha = est$alpha,
                        nu = est$nu,
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
                        psi = est$psi,
                        lambda = est$lambda,
                        theta = est$theta,
                        alpha = est$alpha,
                        nu = est$nu,
                        pmix = pmix,
                        plot_contour = FALSE,
                        labels = labels,
                        show_mi_result = TRUE)
      })
    }  
    num_comb <- length(cai_names) * length(mod_names) + 1 # for specifying the
    # index within ls
    ind <- 1

    while (ind < num_comb) {
      # for each specified CAI
      for (i in seq_along(cai_names)) {
        # swap out the acronym of the composite CAI with the full form
        cai <- lab_cai(substr(cai_names[i], 1, 2))
        # for each specified invariance condition
        for (j in seq_along(mod_names)) {
          # if the specified invariance condition is partial inv.,
          ls[[ind]][, p] <-
            ifelse(rep(mod_names[j] == "par", n_g),
              as.numeric(pinv$summary[cai, 1:n_g]),
              as.numeric(pinv$summary_mi[cai, 1:n_g])
            )

          ylabs <- c(ylabs, paste0(cai, " (", cai_names[i], ")"))

          temp <- ""
          if (mod_names[j] == "par") temp <- "partial invariance"
          if (mod_names[j] == "str") temp <- "strict invariance"

          mains <- c(mains, paste0(cai, " under ", temp))
          ind <- ind + 1
        }
      }
    }
    AIs[,p] <- as.numeric(pinv$ai_ratio)
  }
  
  mains <- mains[-1]
  ylabs <- ylabs[-1]
  
  colorlist <-  c('#000000', '#e6194b', 'lightskyblue', '#f58231',
                  '#bcf60c', '#fabebe', '#911eb4', '#46f0f0',
                  '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000',
                  '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080')
                   #https://sashamaps.net/docs/resources/20-colors/
  if (!is.null(custom_colors)) { colorlist <- custom_colors }
  
  legends <- c(rep("topright", 2), rep("bottomright", 6))
  
  if (!is.null(cai_names)) {
    # iterate over each CAI & invariance condition of interest and produce plots
    for (l in seq_along(ls_names)) {
      plot(rangeVals, ls[[ls_names[l]]][1, ], type = "l", ylim = c(0, 1),
           col = colorlist[1], lwd = 1.5, xlab = xl,
           ylab = ylabs[l],
           main = mains[l],
           cex = 1.1)
      for (i in seq_len(n_g - 1)) {
        lines(rangeVals, ls[[ls_names[l]]][i + 1, ], type = "l",
              lwd = 1.5, col = colorlist[i + 1])
      }
      legend(legends[l], legend = labels, col = colorlist[1:n_g], lty = 1,
             lwd = 1.5, cex = 0.8)  
    }
  }
  if (plotAIs) {
    colorlist <- colorlist[-1]
    plot(rangeVals, AIs[1,], type = "l", ylim = c(0, 1.5),
         col = colorlist[1], lwd = 1.5, xlab = xl,
         ylab = "AI ratio",
         main = paste0("AI ratios [reference group: ", labels[1], "]"),
         cex = 1.1)
    if (n_g > 2) {
      for (i in seq(from = 2, to = n_g - 1)) {
        lines(rangeVals, AIs[i,], type = "l",
              lwd = 1.5, col = colorlist[i])
      }
    }
    legend("bottomright", legend = labels[-1], col = colorlist[1:n_g], lty = 1,
           lwd = 1.5, cex = 0.8)  
  }
}

