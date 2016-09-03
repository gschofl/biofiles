context("gbRecord parser checks (GenBank files)")

if (getOption('biofiles.test.parser')) {
  #nuc.efetch <- reutils::efetch('457866357', 'nuccore', 'gb')
  #save(nuc.efetch, file = "sequences/nuc.efetch.RData")
  load("sequences/nuc.efetch.RData")
  nuc.efetch <- gbRecord(nuc.efetch)
  nuc <- gbRecord("sequences/nucleotide.gbk")
} else {
  #save(nuc.efetch, file = "sequences/nuc.efetch.rda")
  load("sequences/nuc.efetch.rda")
  load("sequences/nucleotide.rda")
}

test_that("GenBank records parse correctly from file and efetch", {
  expect_is(nuc, 'gbRecord')
  expect_is(nuc.efetch, 'gbRecord')
  expect_is(gbRecordList(nuc, nuc.efetch), 'gbRecordList')
})

test_that("Feature list, Sequence, and Seqinfo can be extracted", {
  expect_is(.locus(nuc@seqinfo), 'gbLocus')
  expect_is(.locus(nuc), 'gbLocus')
  
  expect_is(.header(nuc@seqinfo), 'gbHeader')
  expect_is(.header(nuc), 'gbHeader')
  
  expect_is(.sequence(nuc@seqinfo), 'DNAStringSet')
  expect_is(.sequence(nuc), 'DNAStringSet')
  
  expect_is(.features(nuc), 'gbFeatureTable')
  expect_identical(.contig(nuc), NULL)
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
  
  expect_true( all(end(nuc) - start(nuc) + 1 == span(nuc)) )
  expect_equal(unique(strand(nuc)), 1)
  
  expect_is(qualifList(nuc), 'list')
})


context("gbRecord parser checks (GenPept files)")


if (getOption('biofiles.test.parser')) {
  #prot.efetch <- reutils::efetch(c('459479542','379049216'), 'protein', 'gp')
  #save(prot.efetch, file = "sequences/prot.efetch.RData")
  load("sequences/prot.efetch.RData")
  prot.efetch <- gbRecord(prot.efetch)
  prot <- gbRecord("sequences/protein.gp")
} else {
  #save(prot.efetch, file = "sequences/prot.efetch.rda")
  #save(prot, file = "sequences/protein.rda")
  load("sequences/prot.efetch.rda")
  load("sequences/protein.rda")
}

test_that("GenPept records parse correctely from file and efetch", {
  expect_is(prot, 'gbRecord')
  expect_is(prot.efetch, 'gbRecordList')
  expect_is(gbRecordList(prot, prot.efetch), "gbRecordList")
})

test_that("Feature list, Sequence, and Seqinfo can be extracted", {
  expect_is(.locus(prot@seqinfo), 'gbLocus')
  expect_is(.locus(prot), 'gbLocus')
  
  expect_is(.header(prot@seqinfo), 'gbHeader')
  expect_is(.header(prot), 'gbHeader')
  
  expect_is(.sequence(prot@seqinfo), 'AAStringSet')
  expect_is(.sequence(prot), 'AAStringSet')
  
  expect_is(.features(prot), 'gbFeatureTable')
  expect_identical(.contig(prot), NULL)
})

test_that(".dbSource and .defline work", {
  expect_equal(.dbSource(prot), '|gb|')
  expect_equal(.defline(prot),
               "gi|404302394|gb|EJZ56356 CpaF [Pseudomonas fluorescens R124].")
  expect_equal(.dbSource(prot.efetch), c('|gb|', '|gb|'))
})

context("gbRecord parser checks (Embl files)")

if (getOption('biofiles.test.parser')) {
  nuc <- gbRecord("sequences/sample.embl")
  
  test_that("Emble records parse correctly from file", {
    expect_is(nuc, 'gbRecord')
  })
  
  test_that("Feature list, Sequence, and Seqinfo can be extracted", {
    expect_is(.locus(nuc@seqinfo), 'gbLocus')
    expect_is(.locus(nuc), 'gbLocus')
    
    expect_is(.header(nuc@seqinfo), 'gbHeader')
    expect_is(.header(nuc), 'gbHeader')
    
    expect_is(.sequence(nuc@seqinfo), 'DNAStringSet')
    expect_is(.sequence(nuc), 'DNAStringSet')
    
    expect_is(.features(nuc), 'gbFeatureTable')
    expect_identical(.contig(nuc), NULL)
  })
  
  test_that("Accessors work for GenBank records", {
    expect_equal(getAccession(nuc), c("X56734", "S46826"))
    expect_equal(getGeneID(nuc), NA_character_)
    expect_equal(getLength(nuc), 1859)
    expect_equal(getDefinition(nuc), "Trifolium repens mRNA for non-cyanogenic beta-glucosidase")
    
    expect_true( all(end(nuc) - start(nuc) + 1 == span(nuc)) )
    expect_equal(unique(strand(nuc)), 1)

    expect_is(qualifList(nuc), 'list')
  })
  
}