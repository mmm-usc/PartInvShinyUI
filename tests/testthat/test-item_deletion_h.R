# Create sample output

pmix_CESD_r <- 4903/(1903+4903)
lambda_POS_r <- lambda_POS_f <- c(1.00, 1.66, 2.30, 2.29)
nu_POS_r <- c(1.54, 1.36, 1.16, 1.08)
nu_POS_f <- c(0.68, 1.36, 1.16, 1.08)
Theta_POS_r <- diag(c(1.20, 0.81, 0.32, 0.32))
Theta_POS_f <- diag(c(0.72, 0.81, 0.32, 0.32))
alpha_POS_r <- 0
alpha_POS_f <- -0.125
psi_POS_r <- 0.354^2
psi_POS_f <- 0.329^2

CESD_pos <- item_deletion_h(cut_z = 16/60 * 12, 
                            weights_item = c(rep(1, 4)), 
                            weights_latent = 1,
                            alpha_r = alpha_POS_r,
                            alpha_f = alpha_POS_f,
                            psi_r = psi_POS_r,
                            psi_f = psi_POS_f,
                            lambda_r = lambda_POS_r,
                            lambda_f = lambda_POS_f,
                            nu_r = nu_POS_r,
                            nu_f = nu_POS_f,
                            Theta_r = Theta_POS_r, 
                            Theta_f = Theta_POS_f,
                            pmix_ref = pmix_CESD_r, 
                            plot_contour = FALSE,
                            print_formatted = TRUE)

test_that("AI ratio is > 0 for all cells", {
  expect_true(all(CESD_pos$`AI Ratio` > 0))
})


test_that("item_deletion_h() handles matrix input", {
  CESD_pos_mat <- item_deletion_h(cut_z = 16/60 * 12, 
                                  weights_item = c(rep(1,4)), 
                                  weights_latent = 1,
                                  alpha_r = matrix(alpha_POS_r),
                                  alpha_f = matrix(alpha_POS_f),
                                  psi_r = matrix(psi_POS_r),
                                  psi_f = matrix(psi_POS_f),
                                  lambda_r = matrix(lambda_POS_r),
                                  lambda_f = matrix(lambda_POS_f),
                                  nu_r = matrix(nu_POS_r),
                                  nu_f = matrix(nu_POS_f),
                                  Theta_r = Theta_POS_r, 
                                  Theta_f = Theta_POS_f,
                                  pmix_ref = pmix_CESD_r, 
                                  plot_contour = FALSE,
                                  print_formatted = TRUE)
  expect_equal(CESD_pos, CESD_pos_mat)
})


### Helper functions

test_that("redistribute_weights() is working properly", {
error_ex <- c(1:12)
one_dim_w <- c(1:7)
multi_eq_w <- c(1:9)
multi_uneq_w <- c(1:12)

# checking dimension errors
expect_error(redistribute_weights(error_ex, 
                                  n_dim = -3, 
                                  n_i_per_dim = c(3, 6, 3),
                                  del_i = 2))
expect_error(redistribute_weights(error_ex, 
                                  n_dim = 5, 
                                  del_i = 2))

# single dimension
expect_equal(sum(one_dim_w), sum(redistribute_weights(one_dim_w, 
                                                  del_i = 2)))
expect_equal(sum(one_dim_w), sum(redistribute_weights(one_dim_w, 
                                                  n_dim = 1, 
                                                  n_i_per_dim = 7, 
                                                  del_i = 2)))
expect_equal(redistribute_weights(one_dim_w, 
                                  del_i = 2), 
             redistribute_weights(one_dim_w, 
                                  n_dim = 1, 
                                  n_i_per_dim = 7, 
                                  del_i = 2))
# multidimensional, equal number of items in each dim
expect_equal(redistribute_weights(multi_eq_w, 
                                  n_dim = 3, 
                                  del_i = 2),
             redistribute_weights(multi_eq_w, 
                                  n_dim = 3, 
                                  n_i_per_dim = c(3, 3, 3), 
                                  del_i = 2))
expect_equal(sum(multi_eq_w),
             sum(redistribute_weights(multi_eq_w, n_dim = 3, del_i = 2)))

# multidimensional, unequal number of items in dimensions
expect_equal(sum(multi_uneq_w),
             sum(redistribute_weights(multi_uneq_w, 
                                      n_dim = 3, 
                                       n_i_per_dim = c(3, 6, 3), 
                                      del_i=2)))
})


test_that("determine_biased_items() can identify biased items", {
lambda_matrix <- matrix(0, nrow = 5, ncol = 2)
lambda_matrix[1:2, 1] <- c(.322, .655)
lambda_matrix[3:5, 2] <- c(.398, .745, .543)
lambda_matrix2 <- lambda_matrix
lambda_matrix2[1, 1] <- c(.422)
expect_equal(c(2, 3, 4), 
             determine_biased_items(lambda_r = lambda_matrix,
                                    lambda_f = lambda_matrix, 
                                    nu_r = c(.225, .025, .010, .240, .125),
                                    nu_f = c(.225, -.05, .240, -.025, .125),
                                    Theta_r = diag(1, 5), 
                                    Theta_f = diag(c(1, .95, .80, .75, 1))))  
expect_equal(2, determine_biased_items(lambda_r = lambda_matrix,
                                       lambda_f = lambda_matrix, 
                                       nu_r = c(.225, .025, .010, .240, .125),
                                       nu_f =  c(.225, .025, .010, .240, .125),
                                       Theta_r = diag(1, 5), 
                                       Theta_f = diag(c(1, .95, 1, 1, 1))))   
expect_equal(1, determine_biased_items(lambda_r = lambda_matrix,
                                       lambda_f = lambda_matrix2, 
                                       nu_r = c(.225, .025, .010, .240, .125),
                                       nu_f =  c(.225, .025, .010, .240, .125),
                                       Theta_r = diag(1, 5), 
                                       Theta_f = diag(1, 5))) 
     
})


# more example outputs
propsel <- 0.72
cut_z <- NULL
weights_item <- c(rep(1,4)) 
weights_latent <- 1
pmix_ref <- 4903/(1903+4903)
lambda_r <- lambda_f <- c(1.00, 1.66, 2.30, 2.29)
nu_r <- c(1.54, 1.36, 1.16, 1.08)
nu_f <- c(0.68, 1.36, 1.16, 1.08)
Theta_r <- diag(c(1.20, 0.81, 0.32, 0.32))
Theta_f <- diag(c(0.72, 0.81, 0.32, 0.32))
alpha_r <- 0
alpha_f <- -0.125
psi_r <- 0.354^2
psi_f <- 0.329^2
labels <- c("Reference", "Focal")
n_dim <- 1 
n_i_per_dim <- NULL
print_formatted <- TRUE
user_specified_items <- 2
delete_one_cutoff <- NULL
plot_contour <- FALSE  

ex_strict <- PartInv(propsel = propsel, 
                     cut_z = cut_z, 
                     weights_item = weights_item, 
                     weights_latent = weights_latent,
                     alpha_r = alpha_r,
                     alpha_f = alpha_f,
                     psi_r = psi_r,
                     psi_f = psi_f,
                     lambda_r = lambda_f * (1 - pmix_ref) +
                       lambda_r * pmix_ref,
                     nu_r = nu_f * (1 - pmix_ref) + nu_r * pmix_ref,
                     Theta_r = Theta_f * (1 - pmix_ref) + 
                       Theta_r * pmix_ref,
                     pmix_ref = pmix_ref, 
                     plot_contour = plot_contour, 
                     labels = c("Reference", "Focal"), show_mi_result = TRUE)
ex_partial <- PartInv(propsel = propsel, 
                      cut_z = cut_z, 
                      weights_item = weights_item, 
                      weights_latent = weights_latent,
                      alpha_r = alpha_r,
                      alpha_f = alpha_f,
                      psi_r = psi_r,
                      psi_f = psi_f,
                      lambda_r = lambda_r,
                      lambda_f = lambda_f,
                      nu_r = nu_r,
                      nu_f = nu_f,
                      Theta_r = Theta_r,
                      Theta_f = Theta_f,
                      pmix_ref = pmix_ref, 
                      plot_contour = plot_contour,
                      labels = c("Reference", "Focal"))
test_that("cohens_h() is computed correctly for comparing reference with the
          Efocal", {
  out <- c(0.104203486, 0.166642094, -0.101844057, -0.180318247,  0.200493093,
           -0.161048764, 0.211048367, -0.373602451)
  expect_equal(out, cohens_h(ex_partial$summary$Reference, 
                             ex_partial$summary$`E_R(Focal)`))
  
})


# la <- list(
#   matrix(c(0.32, 0.65, 0, 0, 0, 3, 0, 0.3, 0.7, 0.5), nrow = 5, byrow = FALSE),
#   matrix(c(0.32, 0.65, 0, 0, 0, 0, 0, 0.3, 0.7, 0.5), nrow = 5, byrow = FALSE),
#   matrix(c(0.32, 0.65, 0, 0, 0, 0, 0, 0.3, 0.7, 0.5), nrow = 5, byrow = FALSE),
#   matrix(c(0.32, 0.65, 0, 0, 3, 0, 0, 0.3, 0.7, 0.5), nrow = 5, byrow = FALSE)
# )
# 
# lambda_matrix <- matrix(0, nrow = 5, ncol = 2)
# lambda_matrix[1:2, 1] <- c(.3, .6)
# lambda_matrix[3:5, 2] <- c(.4, .75, .53)
# 
# lambda_matrix2 <- lambda_matrix
# lambda_matrix3 <- lambda_matrix
# lambda_matrix2[1,1] <- 3
# lambda_matrix3[2,2] <- 3
# 
# print(find_mismatched_indices(
#   list(lambda_matrix, lambda_matrix2, lambda_matrix3)))
# 
# # Example list with vectors
# lv <- list(
#   c(0.3, 0, 0.3, 3),
#   c(0.3, 0, 0.3, 9),
#   c(0.3, 0, 0.3, 3),
#   c(0.3, 0, 0.3, 3)
# )
# print(find_mismatched_indices(lv))
# 
# # Example list with identical matrices
# lms <- list(
#   matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE),
#   matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE),
#   matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE),
#   matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE)
# )
# 
# print(find_mismatched_indices(lms))
# 
# nus <- list(c(.225, .025, .010, .240, .125), c(.225, -.05, .240, -.025, .125))
# print(find_mismatched_indices(nus))





