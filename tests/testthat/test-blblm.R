library(pryr)

test_that("blblm default work correctly", {

  fit = blblm(mpg ~ wt * hp, mtcars, m = 3, B = 100)

  expect_equal(otype(fit), "S3")
  expect_equal(mode(fit$formula), "call")
  expect_equal(mode(fit$estimates), "list")
})
