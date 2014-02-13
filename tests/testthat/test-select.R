context("Test \"select\"")

x <- getFeatures(gbRecord("sequences/nucleotide.gbk"))
cds <- filter(x, key = "CDS")

test_that("select works with feature indices", {
  expect_equal(.select(cds, "idx"), data.frame(list(idx = c(3,5,7,9,11,13,15,17))))
})

test_that("select works with feature keys", {
  expect_equal(.select(cds, "key"), data.frame(list(key = rep("CDS", 8)), stringsAsFactors = FALSE))
})

test_that("select works with feature locations", {
  expect_equal(dim(.select(cds, "start")), c(8, 1))
  expect_equal(dim(.select(cds, "start", "end")), c(8, 2))
  expect_equal(dim(.select(cds, "start", "end", "width")), c(8, 3))
  expect_equal(dim(.select(cds, "start", "end", "width", "strand")), c(8, 4))
})

test_that("select works with feature qualifiers", {
})