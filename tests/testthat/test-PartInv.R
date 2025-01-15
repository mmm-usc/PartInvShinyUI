# Create sample output
piout <- PartInv(
  propsel = .10,
  kappa_r = 0.5,
  kappa_f = 0,
  phi_r = 1,
  lambda_r = c(.3, .5, .9, .7, .8),
  tau_r = c(.225, .025, .010, .240, .123),
  Theta_r = diag(.96, 5),
  labels = c("female", "male")
)

# Create sample output 2
piout_eq <- PartInv(
  propsel = .10,
  kappa_r = 0.5,
  phi_r = 1,
  lambda_r = c(.3, .5, .9, .7, .8),
  tau_r = c(.225, .025, .010, .240, .123),
  Theta_r = diag(.96, 5),
  labels = c("female", "male")
)

test_that("PartInv() returns a data frame", {
  pidf <- vapply(piout, FUN = inherits, what = c("data.frame"),
                 FUN.VALUE = logical(1))
  expect_true(any(pidf))
})

test_that("plot.PartInv() works successfully", {
  # case where which_result=NULL, uses piout to determine which results are 
  # available
  expect_no_error(plot(piout, labels = c("female", "male")))
  expect_no_error(plot(piout))
  # case where which_result="mi" but mi results were not previously requested
  expect_error(plot(piout, which_result = "mi"))
  piout_with_mi <- PartInv(
    propsel = .10,
    kappa_r = 0.5,
    phi_r = 1,
    lambda_r = c(.3, .5, .9, .7, .8),
    tau_r = c(.225, .025, .010, .240, .123),
    Theta_r = diag(.96, 5),
    labels = c("female", "male"),
    show_mi_result = TRUE
  )
  # case where which_result=mi, and mi results had been requested
  expect_no_error(plot(piout_with_mi, labels = c("female", "male"),
                    which_result = "mi"))
  # case where which_result=NULL, and mi results had been requested
  expect_no_error(plot(piout_with_mi, labels = c("female", "male")))
})

test_that("Identical selection with the same parameters", {
  expect_equal(piout_eq$summary[, 1],
               piout_eq$summary[, 2])
  expect_equal(piout_eq$summary["Proportion selected", 1], .10)
})

# if alpha or alpha_r or kappa_r or cfa_fit not provided
expect_error(PartInv(
  propsel = .10,
  phi_r = 1,
  lambda_r = c(.3, .5, .9, .7, .8),
  tau_r = c(.225, .025, .010, .240, .123),
  Theta_r = diag(.96, 5),
  labels = c("female", "male"),
  show_mi_result = TRUE
)
)

test_that("Duplicated results with `show_mi_result = TRUE` when inputting invariant model", {
  piout_eq2 <- PartInv(
    propsel = .10,
    kappa_r = 0.5,
    phi_r = 1,
    lambda_r = c(.3, .5, .9, .7, .8),
    tau_r = c(.225, .025, .010, .240, .123),
    Theta_r = diag(.96, 5),
    labels = c("female", "male"),
    show_mi_result = TRUE
  )
  expect_equal(piout_eq2$summary[, 1:2], piout_eq2$summary_mi)
})

test_that("PartInv() handles matrix input", {
  piout_mat <- PartInv(
    propsel = .10,
    kappa_r = matrix(0.5),
    kappa_f = matrix(0),
    phi_r = matrix(1),
    lambda_r = matrix(c(.3, .5, .9, .7, .8)),
    tau_r = matrix(c(.225, .025, .010, .240, .123)),
    Theta_r = diag(.96, 5),
    labels = c("female", "male")
  )
  expect_equal(piout[1:4], piout_mat[1:4])
})

test_that("PartInv() handles diagonal input for Theta", {
  piout_diag <- PartInv(
    propsel = .10,
    kappa_r = matrix(0.5),
    kappa_f = matrix(0),
    phi_r = matrix(1),
    lambda_r = matrix(c(.3, .5, .9, .7, .8)),
    tau_r = matrix(c(.225, .025, .010, .240, .123)),
    Theta_r = rep(.96, 5),
    labels = c("female", "male")
  )
  expect_equal(piout[1:4], piout_diag[1:4])
})

test_that("PartInv() output passes logical test", {
  psel <- .2
  # Example favoring reference group
  piout1_pstrict <- PartInv(
    propsel = psel,
    kappa_r = 0,
    kappa_f = -0.1,
    phi_r = 1.3,
    phi_f = 1.2,
    lambda_r = c(1, .5, .9, .7, .8),
    lambda_f = c(1, .8, .9, .7, .8),
    tau_r = c(.225, .025, .010, .240, .123),
    tau_f = c(.025, .025, .010, .240, .053),
    Theta_r = diag(.96, 5),
    Theta_f = diag(c(1, .65, .75, .9, .8)),
    labels = c("reference", "focal")
  )
  piout1_strict <- PartInv(
    propsel = psel,
    kappa_r = 0,
    kappa_f = -0.1,
    phi_r = 1.3,
    phi_f = 1.2,
    lambda_r = c(1, .5, .9, .7, .8),
    tau_r = c(.225, .025, .010, .240, .123),
    Theta_r = diag(.96, 5),
    labels = c("reference", "focal")
  )
  sum_ps <- piout1_pstrict$summary
  sum_s <- piout1_strict$summary
  psel_pstrict <- unlist(sum_ps["Proportion selected", 1:2])
  psel_strict <- unlist(sum_s["Proportion selected", 1:2])
  expect_equal(sum(psel_pstrict), psel * 2)
  expect_equal(sum(psel_strict), psel * 2)
  expect_gt(psel_pstrict["reference"] - psel_strict["reference"], 0)
  expect_lt(abs(
    sum_ps["Specificity", "reference"] -
      sum_ps["C (true negative)", "reference"] /
      (sum_ps["C (true negative)", "reference"] +
         sum_ps["B (false positive)", "reference"])),
    .0005)
})

test_that("`show_mi_result = TRUE` works properly", {
  psel <- .2
  lambda_r <- c(1, .5, .9, .7, .8)
  lambda_f <- c(1, .8, .9, .7, .8)
  tau_r <- c(.225, .025, .010, .240, .123)
  tau_f <- c(.025, .025, .010, .240, .053)
  Theta_r <- diag(.96, 5)
  Theta_f <- diag(c(1, .65, .75, .9, .8))
  # Example favoring reference group
  piout1_pstrict <- PartInv(
    propsel = psel,
    kappa_r = 0,
    kappa_f = -0.1,
    phi_r = 1.3,
    phi_f = 1.2,
    lambda_r = lambda_r,
    lambda_f = lambda_f,
    tau_r = tau_r,
    tau_f = tau_f,
    Theta_r = Theta_r,
    Theta_f = Theta_f,
    labels = c("reference", "focal"),
    pmix_ref = .2,
    show_mi_result = TRUE
  )
  piout1_strict <- PartInv(
    propsel = psel,
    kappa_r = 0,
    kappa_f = -0.1,
    phi_r = 1.3,
    phi_f = 1.2,
    lambda_r = lambda_r * .2 + lambda_f * .8,
    tau_r = tau_r * .2 + tau_f * .8,
    Theta_r = Theta_r * .2 + Theta_f * .8,
    pmix_ref = .2,
    labels = c("reference", "focal")
  )
  expect_equal(piout1_pstrict$summary_mi, piout1_strict$summary[, 1:2])
})

test_that("`show_mi_result = TRUE` uses same cut_z if specified", {
  cut_score <- 4
  lambda_r <- c(1, .5, .9, .7, .8)
  lambda_f <- c(1, .8, .9, .7, .8)
  tau_r <- c(.225, .025, .010, .240, .123)
  tau_f <- c(.025, .025, .010, .240, .053)
  Theta_r <- diag(.96, 5)
  Theta_f <- diag(c(1, .65, .75, .9, .8))
  # Example favoring reference group
  piout1_pstrict <- PartInv(
    cut_z = cut_score,
    kappa_r = 0,
    kappa_f = -0.1,
    phi_r = 1.3,
    phi_f = 1.2,
    lambda_r = lambda_r,
    lambda_f = lambda_f,
    tau_r = tau_r,
    tau_f = tau_f,
    Theta_r = Theta_r,
    Theta_f = Theta_f,
    labels = c("reference", "focal"),
    pmix_ref = .2,
    show_mi_result = TRUE
  )
  piout1_strict <- PartInv(
    cut_z = cut_score,
    kappa_r = 0,
    kappa_f = -0.1,
    phi_r = 1.3,
    phi_f = 1.2,
    lambda_r = lambda_r * .2 + lambda_f * .8,
    tau_r = tau_r * .2 + tau_f * .8,
    Theta_r = Theta_r * .2 + Theta_f * .8,
    pmix_ref = .2,
    labels = c("reference", "focal")
  )
  expect_equal(piout1_pstrict$summary_mi, piout1_strict$summary[, 1:2])
})

