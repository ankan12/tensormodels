test_that("reconstruct as expected", {
  expect_equal(cp_reconstruct(
                list(lambda = c(1, 2),
                     mats = list(diag(2),
                                 matrix(c(1, 0,
                                          1, 2,
                                          1, 0), nrow = 3),
                                 matrix(c(1, 0,
                                          1, 2,
                                          0, 1,
                                          1, 1), nrow = 4)))),
                     array(c(1, 0, 0, 0, 1, 0,
                             0, 4, 0, 2, 0, 0,
                             1, 4, 0, 2, 1, 0,
                             2, 4, 0, 2, 2, 0), dim = c(2, 3, 4)))
})
