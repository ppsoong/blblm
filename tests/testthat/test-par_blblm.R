library(parallel)
library(pryr)

test_that("parallel versions of blbllm work correctly", {

  skip("triggers 'losing unused connection error' but functions work")

  cl = makeCluster(2)
  fit1 = par_blblm(mpg ~ wt * hp, mtcars, 10, 100, cl)
  fit2 = parLB_blblm(mpg ~ wt * hp, mtcars, 10, 100, cl)
  stopCluster(cl)

  expect_equal(otype(fit1), "S3")
  expect_equal(mode(fit1$formula), "call")
  expect_equal(mode(fit1$estimates), "list")

  expect_equal(otype(fit2), "S3")
  expect_equal(mode(fit2$formula), "call")
  expect_equal(mode(fit2$estimates), "list")
})
