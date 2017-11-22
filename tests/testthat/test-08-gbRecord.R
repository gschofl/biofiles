context("gbFeatureTable getter checks")

if (getOption('biofiles.test.parser')) {
  x <- getFeatures(gbRecord("sequences/nucleotide.gbk"))
} else {
  load("sequences/nucleotide.rda")
  x <- getFeatures(nuc)
}

test_that("Sequence, and Seqinfo can be extracted", {
  expect_is(.seqinfo(x), 'seqinfo')
  expect_is(.locus(x), 'gbLocus')
  expect_is(.header(x), 'gbHeader')
  expect_is(.sequence(x), 'DNAStringSet')
  expect_is(getSequence(x), 'DNAStringSet')
})

test_that(".dbSource and .defline work for gbFeatureTables", {
  expect_equal(.dbSource(x), '|gb|')
  expect_match(.defline(x)[1], "lcl|.+|gb|AF229646")
})

test_that("Global accessors work for gbFeatureTables", {
  expect_equal(getLocus(x), 'AF229646')
  expect_equal(getLength(x), 8959)
  expect_equal(getMoltype(x), 'DNA')
  expect_equal(getTopology(x), 'linear')
  expect_equal(getDivision(x), 'BCT')
  expect_equal(getDefinition(x), 'Caulobacter crescentus pilus assembly gene cluster.')
  expect_equal(getAccession(x), 'AF229646')
  expect_equal(getVersion(x), 'AF229646.1')
  expect_equal(getGeneID(x), '7208421')
  expect_equal(getDBLink(x), NA_character_)
  expect_equal(getDBSource(x), NA_character_)
  expect_equal(getSource(x), 'Caulobacter crescentus CB15')
  expect_equal(getOrganism(x), 'Caulobacter crescentus CB15')
  expect_equal(getTaxonomy(x), 'Bacteria; Proteobacteria; Alphaproteobacteria; Caulobacterales; Caulobacteraceae; Caulobacter.')
  expect_equal(getOrganism(x), 'Caulobacter crescentus CB15')
  expect_is(getReference(x), 'gbReferenceList')
  expect_equal(getKeywords(x), '.')
  expect_equal(getComment(x), NA_character_)
})

test_that("Accessors work for gbFeatureTables", {
  expect_equal(index(x), 1:17)
  expect_true( all(key(x) %in% c("source","gene","CDS")) )
  
  x <- x["CDS"]
  expect_true( all(key(x) %in% "CDS") )
  expect_true( all(end(x) - start(x) + 1 == span(x)) )
  
  qualifier_table <- qualif(x, c('gene','protein_id',"db_xref"))
  expect_equal(dim(qualifier_table), c(8, 3))
  
  ans1 <- .qual_access(x, 'product', TRUE, use.names=FALSE)
  expect_equal(ans1[[1]], 'PilA')
  expect_equal(names(ans1[[1]]), NULL)
  
  ans2 <- .qual_access(x, 'product', TRUE, use.names=TRUE)
  expect_equal(unname(ans2[[1]]), 'PilA')
  expect_equal(names(ans2[[1]]), 'product')
  
  expect_equal(.simplify(ans1),
               c("PilA","CpaA","CpaB","CpaC","CpaD","CpaE","CpaF","TadB"))
  expect_equal(.simplify(ans2, unlist=FALSE),
               data.frame(product=c("PilA","CpaA","CpaB","CpaC","CpaD","CpaE","CpaF","TadB"),
                          stringsAsFactors=FALSE))
})


