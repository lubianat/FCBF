library(FCBF)
library(testthat)
data("single_cell_dengue_exprs")
data("single_cell_dengue_annot")


test_that("fcbf works properly", {
  expect_error(fcbf(single_cell_dengue_exprs, single_cell_dengue_annot))
  expect_error(fcbf(t(single_cell_dengue_exprs, single_cell_dengue_annot)))
})
