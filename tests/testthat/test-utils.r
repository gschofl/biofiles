context("Test \"collapse\"")

test_that("collapse works on vectors", {
  vec <- letters[1:3]
  expect_that(collapse(vec), is_identical_to("a b c"))
  expect_that(collapse(vec, ','), is_identical_to("a,b,c"))
  
  ## automatic trimming works
  vec2 <- c("a ", " b ", " c") 
  expect_that(collapse(vec2), is_identical_to("a b c"))
  
  ## NAs preserved
  vec3 <- c("a ", NA, " c") 
  expect_that(collapse(vec3), is_identical_to("a NA c"))
})

test_that("collapse works on lists", {
  vec <- list(letters[1:2], letters[3:4], letters[5:6])
  expect_that(collapse(vec), is_identical_to(c("a b", "c d", "e f")))
  expect_that(collapse(vec, ','), is_identical_to(c("a,b", "c,d", "e,f")))
  
  vec2 <- list(c("a", "b"), NA, NULL, "c", c("d", "e")) 
  expect_that(collapse(vec2), is_identical_to(c("a b", "NA", "", "c", "d e")))
})

context("Test \"linebreak\"")

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

context("Test \".simplify\"")

test_that(".simplify deals with atomic vectors and lists of length-1 elements ", {
  x <- 1:10
  expect_equal(.simplify(x, unlist = TRUE), x)
  expect_equal(.simplify(x, unlist = FALSE), data.frame(X1 = x))
  
  x <- list(1, 2, 3, 4)
  expect_equal(.simplify(x, unlist = TRUE), unlist(x))
  expect_equal(.simplify(x, unlist = FALSE), data.frame(X1 = unlist(x)))
  
  x <- list(c(a = 1), c(a = 2), c(a = 3))
  expect_equal(.simplify(x, unlist = TRUE), unlist(x, use.names = FALSE))
  expect_equal(.simplify(x, unlist = FALSE), data.frame(a = unlist(x)))
})

test_that(".simplify deals with lists of equal-length elements greater than 1", {
  x <- list(c(1, 2), c(3, 4), c(5, 6))
  xout <- data.frame(matrix(1:6, nrow = 3, byrow = TRUE, dimnames = list(NULL, NULL)))
  expect_equal(.simplify(x), xout)
  
  x <- list(c(a = 1, b = 2), c(a = 3, b = 4), c(a = 5, b = 6))
  xout <- data.frame(matrix(1:6, nrow = 3, byrow = TRUE, dimnames = list(NULL, c("a", "b"))))
  expect_equal(.simplify(x), xout)
})

test_that(".simplify deals with lists of unequal-length elements ", {
  ##  is.null(nm)
  x <- list(c(1, 2), c(3, 4, 5), c(6, 7))
  expect_equal(.simplify(x), x)
  
  ##  length(nm) != len
  x <- list(c(a = 1), c(a = 2), c(b = 3))
  xout <- data.frame(matrix(c(1,NA,2,NA,NA,3), ncol = 2, byrow = TRUE,
                            dimnames = list(NULL, c("a", "b"))))
  expect_equal(.simplify(x), xout)
  
  x <- list(c(a = 1, b = 2), c(a = 3, c = 4), c(a = 5, b = 6))
  xout <- data.frame(matrix(c(1,2,NA,
                              3,NA,4,
                              5,6,NA), ncol = 3, byrow = TRUE,
                            dimnames = list(NULL, c("a", "b", "c"))))
  expect_equal(.simplify(x), xout)
  
  x <- list(c(a = 1, b = 2), c(a = 3, b = 4, c = 5), c(a = 6, b = 7))
  xout <- data.frame(matrix(c(1,2,NA,
                              3,4,5,
                              6,7,NA), ncol = 3, byrow = TRUE,
                            dimnames = list(NULL, c("a", "b", "c"))))
  expect_equal(.simplify(x), xout)
})
