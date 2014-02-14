context("gbFeature parser checks")

src <- readLines("sequences/source.gbk")
cds <- readLines("sequences/CDS.gbk")

test_that("GenBank feature tables parse correctly", {
  expect_is(gbFeature(src), 'gbFeature')
  expect_is(gbFeature(cds), 'gbFeature')
})

context("gbFeature getter checks")

test_that("gbFeature accessors work", {
  x <- biofiles:::gbFeature(cds, id = 10)
  
  expect_output(biofiles:::.seqinfo(x), "Seqinfo:")
  expect_output(biofiles:::.header(x), "An empty .gbHeader. instance.")
  expect_output(biofiles:::.locus(x), "An empty .gbLocus. instance.")
  expect_output(biofiles:::.sequence(x), "A BStringSet instance of length 0")
  
  ## index
  expect_equal(index(x), 10)
  
  ## key
  expect_equal(key(x), 'CDS')
  expect_equal(x[['key']], 'CDS')
  
  ## location
  expect_equal(start(x), 338)
  expect_equal(end(x), 2800)
  expect_equal(width(x), 2463)
  expect_equal(strand(x), 1)
  expect_equal(joint_width(x), 2463)
  expect_equal(joint_range(x), c(338, 2800))
  expect_equal(fuzzy(x), matrix(c(FALSE, FALSE), nrow = 1))
  expect_output(show(location(x)), "338..2800")
  expect_output(show(x[['location']]), "338..2800")
  
  ## specific qualifiers
  expect_equal(locusTag(x), "STMUK_0002")
  expect_equal(geneID(x), "thrA")
  expect_equal(product(x), "bifunctional aspartokinase I/homoserine dehydrogenase I")
  expect_equal(proteinID(x), "AEF05936.1")
  expect_is(note(x), "character")
  expect_is(translation(x), "AAStringSet")
  
  ## db_xrefs returns named character vector or NA if a db is not present
  expect_equal(dbxref(x), c(db_xref.taxon = "990282", db_xref.GI = "332986953"))
  expect_equal(dbxref(x, "GI"), c(db_xref.GI = "332986953"))
  expect_equal(dbxref(x, "foo"), c(db_xref.foo = NA_character_))

  # qualif returns named character vector or NA if a db is not present
  expect_equal(qualif(x, 'locus_tag'), c(locus_tag = 'STMUK_0002'))
  
  # get more than one qualifier
  expect_equal(qualif(x, c("locus_tag", "gene")), c(locus_tag = "STMUK_0002", gene = "thrA"))
  expect_equal(qualif(x, "foo"), c(foo = NA_character_))
})
  
  