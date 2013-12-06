context("Test linebreak")

s1 <- "reverse complementary sequence cleaved during processing of trans-spliced tRNAs."
s2 <- "TCTCGCAGAGTTCTTTTTTGTATTAACAAACCCAAAACCCATAGAATTTAATGAACCCAAACCGCAATCGTACAAAAATT"
s3 <- "join(complement(AADE01002756.1:1..10234),gap(1206),AADE01006160.1:1..1963,gap(unk323))"

test_that("linebreak splits at the correct positions", {
  expect_equal(usplit(linebreak(s1, width=50), "\n"),
               c("reverse complementary sequence cleaved during",
                 "processing of trans-spliced tRNAs."))
  
  expect_equal(usplit(linebreak(s3, width=50, split=","), "\n"),
               c("join(complement(AADE01002756.1:1..10234),",
                 "gap(1206),AADE01006160.1:1..1963,gap(unk323))"))
  
  expect_error(linebreak(s2, width=40), "Can't break in the middle of a word.")
  
  expect_equal(nchar(usplit(linebreak(s2, width=40, FORCE=TRUE), "\n")), c(40, 40))  
})


test_that("linebreak indents to the specified position", {
  ss1 <- usplit(linebreak(s1, width=40, indent=10), "\n")
  expect_true(all(sapply(gregexpr("^\\s+", ss1), attr, "match.length") == c(10, -1, -1)))
  expect_true(all(nchar(ss1) < 40))             
  
  ss2 <- usplit(linebreak(s2, width=40, indent=10, FORCE=TRUE), "\n")
  expect_true(all(sapply(gregexpr("^\\s+", ss2), attr, "match.length") == c(10, -1, -1)))
  expect_true(all(nchar(ss2) == c(40, 40, 10))) 
})

test_that("negative indentation reduces the width of the first line", {
  ss1 <- usplit(linebreak(s1, width=50, indent=-10), "\n")
  expect_equal(ss1, c("reverse complementary sequence cleaved",
                      "during processing of trans-spliced tRNAs."))
  
  ss1 <- usplit(linebreak(s1, width=50, indent=-20), "\n")
  expect_equal(ss1, c("reverse complementary",
                      "sequence cleaved during processing of",
                      "trans-spliced tRNAs."))
  
  ss2 <- usplit(linebreak(s2, width=40, indent=-10, FORCE=TRUE), "\n")
  expect_true(all(nchar(ss2) == c(30, 40, 10)))
  
  ss2 <- usplit(linebreak(s2, width=40, indent=-20, FORCE=TRUE), "\n")
  expect_true(all(nchar(ss2) == c(20, 40, 20)))
})

test_that("linebreak offsets to the specified position", {
  ss1 <- usplit(linebreak(s1, width=40, offset=10), "\n")
  match_offset <- sapply(gregexpr("^\\s+", ss1), attr, "match.length")
  expect_true(match_offset[1] == -1)
  expect_true(all(match_offset[-1] == 10))
  expect_true(all(nchar(ss1) < 40)) 
  
  ## expect line1 (40), line2 (10)(30), line3 (10)(10)
  ss2 <- usplit(linebreak(s2, width=40, offset=10, FORCE=TRUE), "\n")
  match_offset <- sapply(gregexpr("^\\s+", ss1), attr, "match.length")
  expect_true(match_offset[1] == -1)
  expect_true(all(match_offset[-1] == 10))
  expect_true(all(nchar(ss2) == c(40, 40, 20)))
})

test_that("indent and offset work together", {
  ss1 <- usplit(linebreak(s1, width=40, indent=10, offset=10), "\n")
  match_shift <- sapply(gregexpr("^\\s+", ss1), attr, "match.length")
  expect_true(all(match_shift == 10))
  expect_true(all(nchar(ss1) < 40))
  
  ss2 <- usplit(linebreak(s2, width=40, indent=10, offset=10, FORCE=TRUE), "\n")
  match_shift <- sapply(gregexpr("^\\s+", ss2), attr, "match.length")
  expect_true(all(match_shift == 10))
  expect_true(all(nchar(ss2) == c(40, 40, 30)))
  
  ss2 <- usplit(linebreak(s2, width=40, indent=-10, offset=10, FORCE=TRUE), "\n")
  expect_true(all(nchar(ss2) == c(30, 40, 30)))
})


