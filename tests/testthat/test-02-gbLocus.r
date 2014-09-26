context("LOCUS and ID line parser checks")

l1 <- "LOCUS       AF171097                1231 bp    DNA     linear   BCT 21-AUG-2001"
l2 <- "LOCUS       CR531020                 876 bp    mRNA    linear   EST 27-DEC-2010"
l3 <- "LOCUS       FQ311868             2988332 bp    DNA     circular BCT 10-JAN-2012"
l4 <- "LOCUS       AAD51968                 143 aa            linear   BCT 21-AUG-2001"

test_that("Locus lines are parsed and rendered correctly", {
  expect_is(gbk_locus(locus_line = l1), 'gbLocus')
  
  expect_equal(gbk_locus(locus_line = l1)$to_string(), l1)
  expect_equal(gbk_locus(locus_line = l2)$to_string(), l2)
  expect_equal(gbk_locus(locus_line = l3)$to_string(), l3)
  expect_equal(gbk_locus(locus_line = l4)$to_string(), l4)
})

id1 <- "ID   X56734; SV 1; linear; mRNA; STD; PLN; 1859 BP."
dt1 <- c("DT   12-SEP-1991 (Rel. 29, Created)", "DT   25-NOV-2005 (Rel. 85, Last updated, Version 11)")

test_that("Embl lines are processed correctely before parsing", {
  expect_equal(embl_line(id1, "ID", "ID"), "X56734; SV 1; linear; mRNA; STD; PLN; 1859 BP.")
  expect_equal(embl_line(dt1, "DT", "DT"), c("12-SEP-1991 (Rel. 29, Created)",
                                             "25-NOV-2005 (Rel. 85, Last updated, Version 11)"))
})

test_that("Id lines are parsed and rendered correctly", {
  ll <- embl_line(id1, "ID", "ID")
  dl <- embl_line(dt1, "DT", "DT")
  expect_is(embl_locus(locus_line = ll, date_line = dl), 'gbLocus')
})
