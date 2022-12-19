library(FixedPoint)
library(testthat)
library(SQUAREM)
context("Testing that previously found bugs have indeed been fixed.")

test_printreports_does_not_error = function(fun, method, startguess){
  # Previously this errored as a result of rounding of digits in the printReports option
  A = FixedPoint(fun, startguess, Method = method, Dampening = 0.5, PrintReports = TRUE)
  return(TRUE)
}

test_that("Testing that previously found errors are fixed", {
  expect_true(test_printreports_does_not_error(function(x){ cos(x) }, "Anderson", 0.3))
  expect_true(test_printreports_does_not_error(function(x){ cos(x) }, "Simple", 0.3))
  expect_true(test_printreports_does_not_error(function(x){ cos(x) }, "Aitken", 0.3))
  expect_true(test_printreports_does_not_error(function(x){ cos(x) }, "Newton", 0.3))
  expect_true(test_printreports_does_not_error(function(x){ cos(x) }, "MPE", 0.3))
  expect_true(test_printreports_does_not_error(function(x){ cos(x) }, "RRE", 0.3))
  expect_true(test_printreports_does_not_error(function(x){ cos(x) }, "VEA", 0.3))
  expect_true(test_printreports_does_not_error(function(x){ cos(x) }, "SEA", 0.3))
})
