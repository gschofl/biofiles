context("gbFeature parser checks")

src <- readLines("sequences/source.gbk")
cds <- readLines("sequences/CDS.gbk")

test_that("GenBank feature tables parse correctly", {
  expect_is(gbFeature(src), 'gbFeature')
  expect_is(gbFeature(cds), 'gbFeature')
})

context("gbFeature getter checks")

test_that("gbFeature accessors work", {
  x <- gbFeature(cds, id=10)
  
  expect_warning(.seqinfo(x), "No header associated with this object")
  expect_output(.header(x), "An empty .gbHeader. instance.")
  expect_output(.locus(x), "An empty .gbLocus. instance.")
  expect_output(.sequence(x), "A BStringSet instance of length 0")
  
  expect_equal(index(x), 10)
  expect_equal(key(x), 'CDS')
  expect_equal(x[['key']], 'CDS')
  
  expect_output(show(location(x)), "338..2800")
  expect_output(show(x[['location']]), "338..2800")
  
  expect_equal(unname(qualif(x, 'locus_tag')), 'STMUK_0002')
  expect_equal(names(qualif(x, 'locus_tag')), 'locus_tag')

  expect_equal(unname(dbxref(x)), c("990282","332986953"))
  expect_equal(names(dbxref(x)), c("taxon","GI"))
  
  expect_equal(unname(dbxref(x, db="gi")), "332986953")
  expect_equal(dbxref(x, db="PFAM"), NA_character_)
  
  expect_equal(start(x), 338)
  expect_equal(end(x), 2800)
  expect_equal(width(x), 2463)
  expect_equal(strand(x), 1)
})
  
  