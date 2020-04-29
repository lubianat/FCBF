library(FCBF)
library(testthat)
data("single_cell_dengue_exprs")
data("single_cell_dengue_annot")


test_that("fcbf works properly", {
  expect_error(fcbf(single_cell_dengue_exprs, single_cell_dengue_annot))
  expect_error(fcbf(t(single_cell_dengue_exprs, single_cell_dengue_annot)))
})



test_that("base entropy functions work", {
  expect_error(entropy(0.5))
  expect_error(entropy(as.character(c("A", "A","B", "B"))))
  expect_equivalent(entropy(as.factor(c("A", "A","B", "B"))), 
                    0.69, 
                    tolerance = 0.01)
  expect_equal(entropy(as.factor(c("A", "A","B", "B")), base =2), 
                    1,
                    tolerance = 0.01)
  
  factor_a = as.factor(c("A", "A","B", "B"))
  factor_b = as.factor(c("A", "A","A", "B"))
  
  expect_equivalent(entropy.joint (factor_a, factor_b ),
                    1.039,
                    tolerance = 0.01)
  
  expect_equivalent(entropy.cond(factor_a, factor_b ),
                    0.477,
                    tolerance = 0.01)
  
  expect_equivalent(SU(factor_a, factor_b), 
                    0.34,
                    tolerance = 0.01)
  
  expect_equivalent(IG(factor_a, factor_b), 
                    0.216,
                    tolerance = 0.01)
  
  
})





