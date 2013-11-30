context("Parse and render LOCUS lines correctly")

l1 <- "LOCUS       AF171097                1231 bp    DNA     linear   BCT 21-AUG-2001"
l2 <- "LOCUS       CR531020                 876 bp    mRNA    linear   EST 27-DEC-2010"
l3 <- "LOCUS       FQ311868             2988332 bp    DNA     circular BCT 10-JAN-2012"
l4 <- "LOCUS       AAD51968                 143 aa            linear   BCT 21-AUG-2001"

test_that("Locus lines are parsed and rendered correctly", {
  expect_equal(gbLocus(l1)$to_string(), l1)
  expect_equal(gbLocus(l2)$to_string(), l2)
  expect_equal(gbLocus(l3)$to_string(), l3)
  expect_equal(gbLocus(l4)$to_string(), l4)
})