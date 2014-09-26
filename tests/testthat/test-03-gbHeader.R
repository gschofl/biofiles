context("Header parser checks")

gbkh <- readLines("sequences/header.gbk")

test_that("Genbank headers parse correctly", {
  expect_is(gbk_header(gbkh), class = "gbHeader")
})

emblh <- readLines("sequences/header.embl")

test_that("Embl headers parse correctly", {
  expect_is(embl_header(emblh), class = "gbHeader")
})