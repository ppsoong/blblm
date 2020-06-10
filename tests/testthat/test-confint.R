test_that("confidence interval works correctly", {

  fit <- blblm(mpg ~ wt * hp, data = mtcars, m = 5, B = 100)
  cft <- confint(fit, c("wt", "hp"))

  expect_equal(mode(cft), "numeric")
  expect_equal(typeof(cft), "double")
})
