test_that("model training works", {
  expect_silent(train_MagicalRsq(file = "data/toy_chr22_50k_integrated.txt.gz", nvar_use = 10000))
})
