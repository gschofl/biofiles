context("gbRecord parser checks (gbk files)")

nuc <- gbRecord("sequences/nucleotide.gbk")
#nuc.efetch <- efetch('457866357', 'nuccore', 'gb')
#save(nuc.efetch, file="inst/tests/sequences/nuc.efetch.RData")
load("sequences/nuc.efetch.RData")
nuc.efetch <- gbRecord(nuc.efetch)

test_that("GenBank records parse correctly from file and efetch", {
  expect_is(nuc, 'gbRecord')
  expect_is(nuc.efetch, 'gbRecord')
  expect_is(gbRecordList(nuc, nuc.efetch), 'gbRecordList')

})

test_that("Feature list, Sequence, and Seqinfo can be extracted", {
  expect_is(seqinfo(nuc), 'Seqinfo')
  expect_is(getFeatures(nuc), 'gbFeatureList')
  expect_is(getSequence(nuc), 'DNAStringSet')
})

test_that(".dbSource and .defline work for gbRecords", {
  expect_equal(.dbSource(nuc), '|gb|')
  expect_equal(.defline(nuc),
               "gi|7208421|gb|AF229646 Caulobacter crescentus pilus assembly gene cluster.")
})

test_that("Accessors work for GenBank records", {
  expect_equal(getAccession(nuc), "AF229646")
  expect_equal(getGeneID(nuc), "7208421")
  expect_equal(getLength(nuc), 8959)
  expect_equal(getDefinition(nuc), "Caulobacter crescentus pilus assembly gene cluster.")
  
  expect_true( all(end(nuc) - start(nuc) + 1 == width(nuc)) )
  expect_equal(unique(strand(nuc)), 1)
  
  expect_is(listQualif(nuc), 'list')
})


context("gbRecord parser checks (GenPept files)")

prot <- gbRecord("sequences/protein.gp")
#prot.efetch <- efetch(c('459479542','379049216'), 'protein', 'gp')
#save(prot.efetch, file="inst/tests/sequences/prot.efetch.RData")
load("sequences/prot.efetch.RData")
prot.efetch <- gbRecord(prot.efetch)

test_that("GenPept records parse correctely from file and efetch", {
  expect_is(prot, 'gbRecord')
  expect_is(prot.efetch, 'gbRecordList')
  expect_error(gbRecordList(prot, prot.efetch),
               "All elements in '...' must be gbRecord objects")
  
})

test_that("Feature list, Sequence, and Seqinfo can be extracted", {
  expect_is(seqinfo(prot), 'Seqinfo')
  expect_is(getFeatures(prot), 'gbFeatureList')
  expect_is(getSequence(prot), 'AAStringSet')
})

test_that(".dbSource and .defline work", {
  expect_equal(.dbSource(prot), '|gb|')
  expect_equal(.defline(prot),
               "gi|404302394|gb|EJZ56356 CpaF [Pseudomonas fluorescens R124].")
  expect_equal(.dbSource(prot.efetch), c('|gb|', '|gb|'))
})

