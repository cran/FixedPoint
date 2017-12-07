library(FixedPoint)
library(testthat)
context("Testing all convergence methods for a difficult cross vector dependent fixed point problem. This means that all entries in vector depend on all other entries.")

# The fixed point vector here is only findable through Anderson, VEA and MPE acceleration. Its approximately value is dput'ed at the bottom of this file

DifficultFunction = function(x){
  x[1] = sum(x)
  for (i in 2:length(x)){
    if (i %% 2 == 0){
      x[i] = log(i) + log(abs(x[i-1]))
    } else {
      x[i] = x[i-1]^(1+log(i))
    }
  }
  return(x)
}


set.seed(1988)
StartVector = rlnorm(12000)

Test_Of_Convergence = function(Function = DifficultFunction, Inputs = StartVector, Outputs = c(), Method = c("Newton") , ConvergenceMetric  = function(Resids){max(abs(Resids))} , ConvergenceMetricThreshold = 1e-12, MaxIter = 1e3, MaxM = 10, Dampening = 1, PrintReports = TRUE, ReportingSigFig = 5, ConditionNumberThreshold = 1e10){

  A = FixedPoint(Function = Function, Inputs = Inputs, Outputs = Outputs, Method = Method, ConvergenceMetric = ConvergenceMetric, ConvergenceMetricThreshold = ConvergenceMetricThreshold, MaxIter = MaxIter, MaxM = MaxM, Dampening = Dampening, PrintReports = PrintReports, ReportingSigFig = ReportingSigFig)

  return(A$Convergence[length(A$Convergence)] < ConvergenceMetricThreshold)
}

test_that("Testing that each method converges in the quadratic case to within tolerance", {
  expect_true(Test_Of_Convergence(Method = "Anderson")) # This takes 12  iterations.
  #expect_true(Test_Of_Convergence(Method = "Simple"))  # Does not converge.
  #expect_true(Test_Of_Convergence(Method = "Aitken"))  # Does not converge.
  #expect_true(Test_Of_Convergence(Method = "Newton"))  # Does not converge.
  expect_true(Test_Of_Convergence(Method = "VEA"))      # This takes 26  iterations.
  #expect_true(Test_Of_Convergence(Method = "SEA"))     # Does not converge.
  expect_true(Test_Of_Convergence(Method = "MPE"))      # This takes 32  iterations.
  #expect_true(Test_Of_Convergence(Method = "RRE"))     # Does not converge.
})

