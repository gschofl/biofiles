context("gbLocation parser checks")

test_that("Location descriptors parse correctly", {
  
  expect_output(show(gbLocation("467")), "467")
  expect_error(gbLocation("12.21"), "Cannot parse location descriptor")
  
  expect_output(show(gbLocation("340..565")), "340..565")
  expect_error(gbLocation("565..340"), "Inadmissible range")
  
  expect_output(show(gbLocation("<345..500")), "<345..500")
  expect_output( show(gbLocation("1..>888")), "1..>888")
  expect_error(gbLocation(">345..500"), "Cannot parse location descriptor")
  expect_error(gbLocation("1..<888"), "Cannot parse location descriptor")
  
  expect_output(show(gbLocation("123^124")), "123\\^124")
  expect_error(gbLocation("123^125"), "Inadmissible range")

  expect_output(show( gbLocation("<10")), "<10")
  expect_output(show( gbLocation(">10")), ">10")
  
  ## This is a stupid idiosycratic usage in IMGT/HLA
  expect_output(show(gbLocation("<1..546>")), "<1..>546")
  
  expect_output( 
    show(gbLocation("join(12..78,134..202)")),
    "join\\(12..78,134..202\\)")
  
  expect_output( 
    show(gbLocation("complement(34..126)")),
    "complement\\(34..126\\)")
  
  expect_output(
    show(gbLocation("complement(join(611..724,856..950))")),
    "complement\\(join\\(611..724,856..950\\)\\)")
  
  expect_output(
    show(gbLocation("join(complement(611..724),complement(856..950))")),
    "complement\\(join\\(611..724,856..950\\)\\)")
  
  expect_output(
    show(gbLocation("join(complement(611..724),856..950)")),
    "join\\(complement\\(611..724\\),856..950\\)")
  
  expect_output(
    show(gbLocation("join(<1..442,1524..>1983)")),
    "join\\(<1..442,1524..>1983\\)")
  
  expect_output( 
    show(gbLocation("J00194.1:100..202")),
    "J00194.1:100..202")
  
  expect_output( 
    show(gbLocation("join(1..100,J00194.1:100..202)")),
    "join\\(1..100,J00194.1:100..202\\)")
  
  expect_output(
    show(gbLocation("order(1..69,1308..1465,1524)")),
    "order\\(1..69,1308..1465,1524\\)")
  
  expect_output( 
    show(gbLocation("bond(12,63)")),
    "bond\\(12,63\\)")
})

test_that("Contig descriptors parse correctly", {
  expect_output( 
    show(gbLocation("join(AACY020610556.1:1..761,gap(51),complement(AACY020885497.1:1..846))")),
    "join\\(AACY020610556.1:1..761,gap\\(51\\),complement\\(AACY020885497.1:1..846\\)\\)")
  
  # it works but the exact output here depends on the width of the console
  #expect_output(
  #  show(gbLocation("join(complement(AADE01002756.1:1..10234),gap(1206),AADE01006160.1:1..1963,gap(unk323),AADE01002525.1:1..11915,gap(),AADE01005641.1:1..2377)")),
  #  "join\\(complement\\(AADE01002756.1:1..10234\\),gap\\(1206\\),AADE01006160.1:1..1963,gap\\(unk323\\),AADE01002525.1:1..11915,gap\\(\\),AADE01005641.1:1..2377\\)")
  
})

a <- gbLocation("join(1..100,J00194.1:100..202)")
b <- gbLocation("complement(join(2691..4571,4918..5163))")
c <- gbLocation("join(<2691..4571,4918..>5163)")

context("gbLocation getter/setter checks")

test_that("gbLocation accessors work", {
  
  expect_equal(start(a), c(1, 100))
  expect_equal(start(a, join = TRUE), 1)
  expect_equal(end(a), c(100, 202))
  expect_equal(end(a, join = TRUE), 202)
  expect_equal(span(a), c(100, 103))
  expect_equal(span(a, join = TRUE), 202)
  expect_equal(strand(a), c(1, 1))
  expect_equal(strand(b), c(-1, -1))
  expect_equal(fuzzy(c), matrix(c(TRUE,FALSE,FALSE,TRUE), nrow=2))
  expect_equal(getAccession(a), c("","J00194.1"))
  
  expect_equal(start(gbLocation("102")), 102)
  expect_equal(end(gbLocation("102")), 102)
  expect_equal(span(gbLocation("102")), 1)
  
  expect_equal(fuzzy(gbLocation("<102")), matrix(c(TRUE,FALSE), nrow=1))
  expect_equal(fuzzy(gbLocation(">102")), matrix(c(FALSE,TRUE), nrow=1))
})

test_that("gbLocation replacement methods work", {
  
  expect_error(start(a) <- 10, "contains 2 start values")
  expect_error(end(a) <- 10, "contains 2 end values")
  
  start(a) <- c(10, 110)
  expect_equal(start(a), c(10, 110))
  
  expect_error(
    end(a) <- c(9, 200),
    "One or more ranges with second endpoint before first")
  
  strand(a) = -1
  expect_equal(strand(a), c(-1, -1))

})

test_that("gbLocation shift method work", {
  a <- gbLocation("join(1..100,100..202)")
  a_shifted <- shift(a, 100)
  expect_identical(a, shift(a_shifted, -100)) 
})

