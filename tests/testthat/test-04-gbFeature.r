context("Feature parser checks")

src1 <- readLines("sequences/source.gbk")
cds1 <- readLines("sequences/CDS.gbk")

test_that("GenBank feature tables parse correctly", {
  expect_is(gbFeature(src1), 'gbFeature')
  expect_is(gbFeature(cds1), 'gbFeature')
})

src2 <- readLines("sequences/source.embl")
cds2 <- readLines("sequences/CDS.embl")

test_that("Embl feature tables parse correctly", {
  expect_is(gbFeature(src2), 'gbFeature')
  expect_is(gbFeature(cds2), 'gbFeature')
})

context("gbFeature getter/setter checks")

test_that("gbFeature accessors work", {
  x <- biofiles:::gbFeature(cds2, id = 10)
  
  expect_output(biofiles:::.seqinfo(x), "Seqinfo:")
  expect_output(biofiles:::.header(x), "An empty .*gbHeader.* instance.")
  expect_output(biofiles:::.locus(x), "An empty .*gbLocus.* instance.")
  expect_output(biofiles:::.sequence(x), "A BStringSet instance of length 0")
  
  ## index
  expect_equal(index(x), 10)
  
  ## key
  expect_equal(key(x), 'CDS')
  expect_equal(x[['key']], 'CDS')
  
  ## location
  expect_equal(start(x), 14)
  expect_equal(end(x), 1495)
  expect_equal(span(x), 1482)
  expect_equal(strand(x), 1)
  #expect_equal(joint_range(x), c(338, 2800))
  expect_equal(fuzzy(x), matrix(c(FALSE, FALSE), nrow = 1))
  expect_output(show(location(x)), "14..1495")
  expect_output(show(x[['location']]), "14..1495")

  x <- biofiles:::gbFeature(cds1, id = 10)
  
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

test_that("gbFeature replacement methods work", {
  x <- biofiles:::gbFeature(cds1, id = 10)
  
  ## replace start
  start(x) <- 100
  expect_equal(start(x), 100)
  
  ## replace end
  end(x) <- 200
  expect_equal(end(x), 200)
  
  ## try to make start > end
  expect_error(start(x) <- 300)
  ## try again, switching off validity checks
  start(x, check = FALSE) <- 300
  expect_output(location(x), "300..200")
  
  ## replace key
  ## Should we check for valid keys?
  key(x) <- "FOO"
  expect_equal(key(x), "FOO")
  
  ## replace qulifier
  qualif(x, "gene") <- "bar"
  expect_equal(qualif(x, "gene"), c(gene = "bar"))
  
  ## try replacing qualifier without specifying which
  expect_error(qualif(x) <- "baz")
})
  