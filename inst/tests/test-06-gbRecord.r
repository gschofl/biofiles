context("gbFeature getter checks")

x <- select(gbRecord(gb="sequences/nucleotide.gbk"), key='CDS')[[3]]

test_that("Sequence, and Seqinfo can be extracted", {
  expect_is(.seqinfo(x), 'seqinfo')
  expect_is(.locus(x), 'gbLocus')
  expect_is(.header(x), 'gbHeader')
  expect_is(.sequence(x), 'DNAStringSet')
})

test_that(".dbSource and .defline work for gbFeatures", {
  expect_equal(.dbSource(x), '|gb|')
  expect_match(.defline(x), "lcl|CDS.+|gb|AF229646")
})

test_that("Accessors work for gbFeatures", {
  expect_equal(getAccession(x), "AF229646")
  expect_equal(getLength(x), 8959)
  expect_equal(getDefinition(x), "Caulobacter crescentus pilus assembly gene cluster.")
  
  expect_true(end(x) - start(x) + 1 == width(x))
  expect_equal(strand(x), 1)
  expect_equal(fuzzy(x), matrix(c(FALSE,FALSE), nrow=1))
  
  expect_equal(index(x), 7)
  expect_equal(key(x), 'CDS')
  expect_output(show(location(x)), "1521..2414")
  
  expect_equal(qualif(x, 'gene', use.names=FALSE), 'cpaB')
  expect_equal(qualif(x, 'bla', use.names=FALSE), NA_character_)
  expect_equal(length(qualif(x, c("gene","protein_id"))), 2)
  
  expect_equal(unname(dbxref(x)), "7208424")
  expect_equal(dbxref(x, 'bla'), NA_character_)
  
  expect_equal(locusTag(x), NA_character_)
  expect_equal(product(x), "CpaB")
  expect_equal(proteinID(x), "AAF40191.1")
  expect_equal(note(x), "required for pilus assembly in Caulobacter")
  expect_is(translation(x), "AAStringSet")
  
})


test_that("Ranges work for gbFeatures", {
  range <- ranges(x)
  expect_is(range, 'GRanges')
  expect_equal(start(range), 1521)
  expect_equal(end(range), 2414)
  expect_equal(IRanges::width(ranges(range)), 894)
  expect_equal(names(range), "cpaB")
})


test_that("getSequence works for gbFeatures", {
  seq <- getSequence(x)
  expect_is(seq, "DNAStringSet")
  expect_equal(names(seq), .defline(x))
  expect_equal(length(seq[[1]]), width(x))  
})


context("gbFeature setter checks")

test_that("Setters work for gbFeatures", {
  strand(x) <- -1
  expect_output(show(location(x)), "complement\\(1521..2414\\)")  
})  


