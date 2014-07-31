context("gbFeatureTable accessor checks")

if (getOption('biofiles.test.parser')) {
  nuc <- gbRecord("sequences/nucleotide.gbk")
} else {
  #save(nuc, file = "sequences/nucleotide.rda")
  load("sequences/nucleotide.rda")
}
fl <- .features(nuc)

test_that("Can extract features", {
  expect_is(fl, 'gbFeatureTable')
  expect_is(getFeatures(nuc), 'gbFeatureTable')
})

test_that("Can subset a gbFeatureTable", {
  cds <- fl["CDS"]
  expect_equal(unique(key(cds)), "CDS")
  
  expect_is(fl[1:5], "gbFeatureTable")
  expect_is(fl[[1]], "gbFeature")
  
  expect_equal(getLocus(fl[[1]]), "AF229646")  
})

test_that("gbFeatureTable accessors work", {
  x <- fl
    ## index
  expect_equal(index(x), 1:17)
  
  ## key
  expect_equal(unique(key(x)), c('source', 'gene', 'CDS'))
  
  ## location
  expect_is(start(x), 'integer')
  expect_is(end(x), 'integer')
  expect_is(span(x), 'integer')
  expect_is(strand(x), 'integer')
  expect_is(location(x), "list")
  expect_equal(fuzzy(x), matrix(rep(FALSE, 2*17), ncol = 2))
  
})

