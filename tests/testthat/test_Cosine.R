library(FixedPoint)
library(testthat)
context("Testing all convergence methods for a simple cosine function fixed point problem.")

Test_Of_Convergence = function(Function = function(x){ cos(x) }, Inputs = c(0.3, 0.4,0.5,0.6,0.7,0.8), Outputs = c(), Method = c("Newton") , ConvergenceMetric  = function(Resids){max(abs(Resids))} , ConvergenceMetricThreshold = 1e-10, MaxIter = 1e3, MaxM = 10, Dampening = 1, PrintReports = TRUE, ReportingSigFig = 5, ConditionNumberThreshold = 1e10){

  A = FixedPoint(Function = Function, Inputs = Inputs, Outputs = Outputs, Method = Method, ConvergenceMetric = ConvergenceMetric, ConvergenceMetricThreshold = ConvergenceMetricThreshold, MaxIter = MaxIter, MaxM = MaxM, Dampening = Dampening, PrintReports = PrintReports, ReportingSigFig = ReportingSigFig)

  return((A$Convergence[length(A$Convergence)] < ConvergenceMetricThreshold))
}

test_that("Testing that each method converges in the cos(x)  case to within tolerance", {
  expect_true(Test_Of_Convergence(Method = "Anderson")) # This takes 10  iterations.
  expect_true(Test_Of_Convergence(Method = "Simple"))   # This takes 58  iterations.
  expect_true(Test_Of_Convergence(Method = "Aitken"))   # This takes 11  iterations.
  expect_true(Test_Of_Convergence(Method = "Newton"))   # This takes 9   iterations.
  expect_true(Test_Of_Convergence(Method = "MPE"))      # This takes 13  iterations.
  expect_true(Test_Of_Convergence(Method = "RRE"))      # This takes 13  iterations.
  expect_true(Test_Of_Convergence(Method = "VEA"))      # This takes 13  iterations.
  expect_true(Test_Of_Convergence(Method = "SEA"))      # This takes 13  iterations.
})

Test_Of_Convergence = function(Function = function(x){ cos(x) }, Inputs = c(0.3, 0.4), Outputs = c(), Method = c("Newton") , ConvergenceMetric  = function(Resids){max(abs(Resids))} , ConvergenceMetricThreshold = 1e-10, MaxIter = 1e3, MaxM = 10, Dampening = 1, PrintReports = TRUE, ReportingSigFig = 5, ConditionNumberThreshold = 1e10){

  A = FixedPoint(Function = Function, Inputs = Inputs, Outputs = Outputs, Method = Method, ConvergenceMetric = ConvergenceMetric, ConvergenceMetricThreshold = ConvergenceMetricThreshold, MaxIter = MaxIter, MaxM = MaxM, Dampening = Dampening, PrintReports = PrintReports, ReportingSigFig = ReportingSigFig)

  return((A$Convergence[length(A$Convergence)] < ConvergenceMetricThreshold))
}

test_that("Testing that each method converges in the cos(x)  case to within tolerance", {
  expect_true(Test_Of_Convergence(Method = "Anderson")) # This takes 8   iterations.
  expect_true(Test_Of_Convergence(Method = "Simple"))   # This takes 58  iterations.
  expect_true(Test_Of_Convergence(Method = "Aitken"))   # This takes 11  iterations.
  expect_true(Test_Of_Convergence(Method = "Newton"))   # This takes 9   iterations.
  expect_true(Test_Of_Convergence(Method = "MPE"))      # This takes 17  iterations.
  expect_true(Test_Of_Convergence(Method = "RRE"))      # This takes 19  iterations.
  expect_true(Test_Of_Convergence(Method = "VEA"))      # This takes 13  iterations.
  expect_true(Test_Of_Convergence(Method = "SEA"))      # This takes 13  iterations.

})

Test_Of_Convergence = function(Function = function(x){ cos(x) }, Inputs = 0.3, Outputs = c(), Method = c("Newton") , ConvergenceMetric  = function(Resids){max(abs(Resids))} , ConvergenceMetricThreshold = 1e-10, MaxIter = 1e3, MaxM = 10, Dampening = 1, PrintReports = TRUE, ReportingSigFig = 5, ConditionNumberThreshold = 1e10){

  A = FixedPoint(Function = Function, Inputs = Inputs, Outputs = Outputs, Method = Method, ConvergenceMetric = ConvergenceMetric, ConvergenceMetricThreshold = ConvergenceMetricThreshold, MaxIter = MaxIter, MaxM = MaxM, Dampening = Dampening, PrintReports = PrintReports, ReportingSigFig = ReportingSigFig)

  return((A$Convergence[length(A$Convergence)] < ConvergenceMetricThreshold))
}

test_that("Testing that each method converges in the cos(x)  case to within tolerance", {
  expect_true(Test_Of_Convergence(Method = "Anderson")) # This takes 7   iterations.
  expect_true(Test_Of_Convergence(Method = "Simple"))   # This takes 58  iterations.
  expect_true(Test_Of_Convergence(Method = "Aitken"))   # This takes 11  iterations.
  expect_true(Test_Of_Convergence(Method = "Newton"))   # This takes 9   iterations.
  expect_true(Test_Of_Convergence(Method = "MPE"))      # This takes 19  iterations.
  expect_true(Test_Of_Convergence(Method = "RRE"))      # This takes 25  iterations.
  expect_true(Test_Of_Convergence(Method = "VEA"))      # This takes 13  iterations.
  expect_true(Test_Of_Convergence(Method = "SEA"))      # This takes 13  iterations.
})
