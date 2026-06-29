test_that("calculate_sparsity works correctly", {
  table <- data.frame(
    id = c("f1", "f2", "f3"),
    s1 = c(10, 0, 5),
    s2 = c(0, 0, 3),
    s3 = c(8, 4, 0)
  )

  sparsity <- calculate_sparsity(table)

  expect_type(sparsity, "double")
  # 4 zeros out of 9 cells = 0.444...
  expect_equal(sparsity, 4/9, tolerance = 0.001)
})

test_that("convert_to_relative converts correctly", {
  table <- data.frame(
    id = c("f1", "f2"),
    s1 = c(10, 90),   # Total: 100
    s2 = c(20, 80)    # Total: 100
  )

  result <- convert_to_relative(table)

  expect_s3_class(result, "data.frame")
  expect_equal(result$s1[1], 0.1, tolerance = 0.001)
  expect_equal(result$s1[2], 0.9, tolerance = 0.001)
})

test_that("calculate_feature_cv calculates coefficient of variation", {
  table <- data.frame(
    id = c("f1", "f2"),
    s1 = c(10, 20),
    s2 = c(20, 40),
    s3 = c(30, 60)
  )

  cv <- calculate_feature_cv(table)

  expect_s3_class(cv, "data.frame")
  expect_true("cv" %in% colnames(cv))
  # f2 has same CV as f1 (proportional values)
  expect_equal(cv$cv[1], cv$cv[2], tolerance = 0.001)
})
