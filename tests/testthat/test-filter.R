context("Test \"filter\"")

if (getOption('biofiles.test.parser')) {
  x <- getFeatures(gbRecord("sequences/nucleotide.gbk"))
} else {
  load("sequences/nucleotide.rda")
  x <- getFeatures(nuc)
}

test_that("filter works on feature indices", {
  expect_equal(key(filter(x, idx = 1)), "source")
  expect_equal(index(filter(x, index = 2)), 2)
  expect_equal(length(filter(x, idx = 1:10)), 10)
  expect_equal(length(filter(x, idx = 0)), 0)
})

test_that("filter works on feature keys", {
  expect_true(all(index(filter(x, key = "CDS")) == index(x["CDS"])))
  expect_true(length(filter(x, key = c("source", "gene"))) == 9)
})

test_that("filter works on ranges", {
  expect_true(max(start(filter(x, range = "..1500"))) < 1500)
  expect_true(min(end(filter(x, range = "2200.."))) > 2200)
  expect_equal(geneID(filter(x, range = "1500..2400", key = "CDS")), "cpaB")
})

test_that("filter works on arbitrary qualifiers", {
  expect_equal(product(filter(x, product = "CpaF")), "CpaF")
  expect_equal(geneID(filter(x, gene = "cpaE", key = "CDS")), "cpaE")
  expect_true(all(key(filter(x, db_xref = "GI:")) == "CDS"))
})

test_that("filter used with the .cols argument", {
  df <- filter(x, key = "CDS", .cols = c("start", "end", "gene"))
  expect_is(df, "data.frame")
  expect_equal(colnames(df), c("start", "end", "gene"))
})
  

