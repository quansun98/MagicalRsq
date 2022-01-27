test_that("model training works", {
  expect_silent(train_deepRsq(file = "data/chr22_BioMe_HRC.gz", nvar_use = 10000))
})
