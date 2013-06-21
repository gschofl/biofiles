context("gbFeatureList getter checks")

x <- getFeatures(gbRecord("sequences/nucleotide.gbk"))

test_that("Sequence, and Seqinfo can be extracted", {
  expect_is(seqinfo(x), 'Seqinfo')
  expect_is(getSequence(x), 'DNAStringSet')
})

test_that(".dbSource and .defline work for gbFeatureLists", {
  expect_equal(.dbSource(x), '|gb|')
  expect_match(.defline(x)[1], "lcl|.+|gb|AF229646")
})

test_that("Accessors work for gbFeatureLists", {
  expect_equal(getAccession(x), "AF229646")
  expect_equal(getLength(x), 8959)
  expect_equal(getDefinition(x), "Caulobacter crescentus pilus assembly gene cluster.")
  
  expect_equal(index(x), 1:17)
  expect_true( all(key(x) %in% c("source","gene","CDS")) )
  
  x <- x["CDS"]
  expect_true( all(key(x) %in% "CDS") )
  expect_true( all(end(x) - start(x) + 1 == width(x)) )
  
  qualifier_table <- qualif(x, c('gene','protein_id',"db_xref"))
  expect_equal(dim(qualifier_table), c(8, 3))
})


