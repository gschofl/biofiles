context("gbFeatureList accessor checks")

nuc <- gbRecord("sequences/nucleotide.gbk")
fl <- .features(nuc)

test_that("Can extract features", {
  expect_is(fl, 'gbFeatureList')
  expect_is(getFeatures(nuc), 'gbFeatureList')
})

test_that("Can subset a gbFeatureList", {
  cds <- fl["CDS"]
  expect_equal(unique(key(cds)), "CDS")
  
  expect_is(fl[1:5], "gbFeatureList")
  expect_is(fl[[1]], "gbFeature")
  
  expect_equal(getLocus(fl[[1]]), "AF229646")  
})