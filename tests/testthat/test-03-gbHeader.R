context("Header parser checks")

gbkh <- readLines("sequences/header.gbk")

test_that("Genbank headers parse correctly", {
  expect_is(gbk_header(gbkh), class = "gbHeader")
})

emblh <- readLines("sequences/header.embl")

test_that("Embl headers parse correctly", {
  expect_is(embl_header(emblh), class = "gbHeader")
})

gbkh2 <- readLines("sequences/header_no_ref.gbk")

test_that("Genbank headers with no reference parse correctly", {
  x <- gbk_header(gbkh2)
  expect_is(x, class = "gbHeader")
  expect_equal(x@.xData$source, "")
  expect_equal(x@.xData$organism, "")
  expect_equal(x@.xData$taxonomy, ".")
})
