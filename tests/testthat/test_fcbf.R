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
  
  fcbf_standard_output = fcbf(mock_expression_table, 
                              mock_cell_class)
  
  expect_equal(rownames(fcbf_standard_output), c("A","B"))
  expect_equal(fcbf_standard_output["SU"][,1][2], 0.618977)
  
  
  fcbf_output_for_2_genes_in_filter = fcbf(mock_expression_table, 
                                           mock_cell_class,
                                           n_genes_selected_in_first_step = 2)
  
  expect_equal(nrow(fcbf_output_for_2_genes_in_filter), 1)
  expect_equal(fcbf_output_for_2_genes_in_filter["SU"][,1][1], 1)
  
})



test_that("base entropy functions work", {
  
  expect_error(get_entropy_for_vector(0.5))
  expect_error(get_entropy_for_vector(as.character(factor_a)))
  
  
  
  expect_equal(get_entropy_for_vector(factor_a), 
                    1, 
                    tolerance = 0.01)
  expect_equivalent(get_entropy_for_vector(factor_a, base = exp(1)), 
                    0.69,
                    tolerance = 0.01)

  
  expect_equivalent(get_joint_entropy_for_vectors (factor_a, factor_b ),
                    1.36,
                    tolerance = 0.01)
  
  expect_equivalent(get_conditional_entropy_for_vectors(factor_a, factor_b ),
                    0.39,
                    tolerance = 0.01)
  expect_equivalent(get_conditional_entropy_for_vectors(factor_b, factor_a ),
                    0.36,
                    tolerance = 0.01)
  
  
  expect_equal(get_SU_for_vector_pair(factor_a, factor_a), 
                    1)
  expect_equivalent(get_SU_for_vector_pair(factor_a, factor_b), 
                    0.619,
                    tolerance = 0.01)
  
  
  expect_equal(get_IG_for_vector_pair(factor_a, factor_a), 
                    1)
  expect_equivalent(get_IG_for_vector_pair(factor_a, factor_b), 
               0.61,
               tolerance = 0.01)
  
  
})





