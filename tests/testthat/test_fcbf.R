library(FCBF)
library(testthat)


####### Variables for testing:###### 

factor_a = as.factor(c("OFF", "OFF","OFF","OFF","OFF",     "ON", "ON","ON", "ON", "ON"))
factor_b = as.factor(c("OFF", "OFF","OFF","OFF","OFF",      "OFF","ON", "ON","ON", "ON"))
factor_c = as.factor(c("OFF", "OFF","OFF","OFF","OFF",      "OFF","OFF", "ON","ON", "ON"))
factor_d = as.factor(rep("ON", 10))

mock_feature_table = data.frame(A = factor_a, 
                                   B = factor_b,
                                   B2 = factor_b,
                                   B3 = factor_b,
                                   C = factor_c,
                                   C2 = factor_c,
                                   C3 = factor_c,
                                   D = factor_d)

# For FCBF, genes are in columns, and cells in rows
mock_expression_table = as.data.frame(t(mock_feature_table))


mock_cells = c()


for (n in 1:10){
  mock_cells = c(mock_cells,   paste("cell",n))
}

colnames(mock_expression_table) = mock_cells


mock_cell_class = as.factor(c(rep("monocyte", 5), rep("B cell", 5)))

####### Tests ####### 

test_that("fcbf works properly", {
  
  fcbf_output = fcbf(mock_expression_table, mock_cell_class, n_genes = 5)
  
  expect_equal(rownames(fcbf_output), c("A","B"))
  expect_equal(fcbf_output["SU"][,1][1], 1)
  
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
  
  expect_equivalent(get_joint_entropy_for_vectors (factor_a, factor_b ),
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





