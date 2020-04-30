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
mock_discrete_expression_table = as.data.frame(t(mock_feature_table))


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

test_that("discretization base functions work", {
  
  twenty_numeric_vector <- c(0,1,2,3,3,3,4,5,6,7,8,8,8,8,8,8,8,8,8,9)
  
  # Binary cutoff cor n_bin = 4 should be (0+9)/4 = 2.25
  vw_expectation <- c(rep("low",3), rep("high", 17))
  expect_equal(.split_vector_in_two_varying_width(twenty_numeric_vector,
                                                  n_of_bins = 4), 
               vw_expectation)
  
  # varying_width and min_max_cutoff are related. cutoff = 1/n_of_bins  
  expect_equal(.split_vector_in_two_varying_width(twenty_numeric_vector,
                                                  n_of_bins = 4), 
               .split_vector_in_two_by_min_max_thresh(twenty_numeric_vector,
                                                      cutoff = 1/4))
  
  expect_equal(.split_vector_in_two_varying_width(twenty_numeric_vector,
                                                  n_of_bins = 6), 
               .split_vector_in_two_by_min_max_thresh(twenty_numeric_vector,
                                                      cutoff = 1/6))
  
  
  # mean is 5.75
  mean_expectation <- c(rep("low",8), rep("high", 12))
  expect_equal(.split_vector_in_two_by_mean(twenty_numeric_vector), 
               mean_expectation)
  
  
  # mean is 7.5
  median_expectation <- c(rep("low",10), rep("high", 10))
  expect_equal(.split_vector_in_two_by_median(twenty_numeric_vector), 
               median_expectation)
  
  
  # I was not able to set a test for kmwans. Results are not stable.
  expect_vector(.split_vector_by_kmeans(twenty_numeric_vector, centers = 2))
  
  # mean = 5.75, sd = 2.844663
  mean_sd_expectation <- c(rep("low",3), rep("medium",16), rep("high", 1))
  expect_vector(.split_vector_in_three_by_mean_sd(twenty_numeric_vector))
  
  
  
})


counts_a = as.factor(c(rep(0,5), rep(10,5)))
counts_b = as.factor(c(rep(0,3), rep(6,3), rep(10,4)))
counts_c = as.factor(c(rep(0,6), rep(10,4)))
counts_d = as.factor(rep(100, 10))

mock_feature_table_counts = data.frame(A = counts_a, 
                                B = counts_b,
                                B2 = counts_b,
                                B3 = counts_b,
                                C = counts_c,
                                C2 = counts_c,
                                C3 = counts_c,
                                D = counts_d)

mock_expression_table_counts = data.frame(t(mock_feature_table_counts))

test_that("discretization of expression table works",{
  
  
  discretized_exprs = discretize_exprs(mock_expression_table_counts, number_of_bins = 3,
                                       method = "varying_width")
  
  expected_a = as.factor(c(rep("low",5), rep("high", 5)))
  expect_equivalent(as.factor(t(discretized_exprs["A",])), expected_a)
  
  expected_b = as.factor(c(rep("low",3), rep("high", 7)))
  expect_equivalent(as.factor(t(discretized_exprs["B",])), expected_b)
  
  expected_c = as.factor(c(rep("low",6), rep("high", 4)))
  expect_equivalent(as.factor(t(discretized_exprs["C",])), expected_c)
  
  expected_d = as.factor(rep("high", 10))
  expect_equivalent(as.factor(t(discretized_exprs["D",])), expected_d)
  
})
