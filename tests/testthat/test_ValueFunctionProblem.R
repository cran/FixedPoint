#library(FixedPoint)
#library(testthat)
#library(SQUAREM)
#library(schumaker)
#library(cubature)
#cat(paste0("Doing Value Function Test\n"))


#ValueGivenShock = function(Budget, epsilon, NextValueFunction, delta, beta, BudgetStateSpace){
#  optimize(f = function(x) epsilon*(x^delta) + beta*NextValueFunction(Budget - x + 1), lower = 0, upper = Budget, maximum = TRUE)$objective
#}

#ExpectedUtility = function(Budget, NextValueFunction, delta, beta, BudgetStateSpace){
#  if (Budget > 0.001){
#    adaptIntegrate(f = function(epsilon) ValueGivenShock(Budget, epsilon, NextValueFunction, delta, beta, BudgetStateSpace)* dlnorm(epsilon), lowerLimit = qlnorm(0.0001), upperLimit = qlnorm(0.9999))$integral
#  } else {
#    beta*NextValueFunction(1)
#  }
#}

#OneIterateBudgetValues = function(BudgetValues, delta, beta, BudgetStateSpace){
#  NextValueFunction = schumaker::Schumaker(BudgetStateSpace, BudgetValues, Extrapolation = "Linear")$Spline
#  for (i in 1:length(BudgetStateSpace)){
#    BudgetValues[i] = ExpectedUtility(BudgetStateSpace[i], NextValueFunction, delta, beta, BudgetStateSpace)
#  }
#  return(BudgetValues)
#}
#A = FixedPoint(OneIterateBudgetValues, InitialGuess, Method = c("Simple"))

#Test_Of_Convergence = function(delta, beta, BudgetStateSpace, Inputs,  Function = OneIterateBudgetValues, Outputs = c(), Method = c("Newton") , ConvergenceMetric  = function(Resids){max(abs(Resids))} , ConvergenceMetricThreshold = 1e-10, MaxIter = 2500, MaxM = 10, Dampening = 1, PrintReports = FALSE, ReportingSigFig = 5, ConditionNumberThreshold = 1e10,   ChangePerIteratexaxis = BudgetStateSpace, ConvergenceFigLags = 5, DropOldIterates = TRUE){
#  A = FixedPoint(Function = function(x) Function(x, delta, beta, BudgetStateSpace), Inputs = Inputs, Outputs = Outputs, Method = Method, ConvergenceMetric = ConvergenceMetric, ConvergenceMetricThreshold = ConvergenceMetricThreshold, MaxIter = MaxIter, MaxM = MaxM, Dampening = Dampening, PrintReports = PrintReports, ReportingSigFig = ReportingSigFig, ChangePerIteratexaxis = ChangePerIteratexaxis, ConvergenceFigLags = ConvergenceFigLags, DropOldIterates = DropOldIterates)
#  cat(paste0("The ", Method, " method took ", length(A$Convergence), " iterations and finished with ", A$Finish, "\n"))
#  return(list(Convergence = A$Convergence[length(A$Convergence)] < ConvergenceMetricThreshold, Result = A))
#}

# Note this is a slow function so we only run Anderson below.
#test_that("Testing that each method converges in the quadratic case to within tolerance", {
    #delta = 0.2
    #beta = 0.95
    #BudgetStateSpace = c(seq(0,1, 0.015), seq(1.05,3,0.05))
    #InitialGuess = sqrt(BudgetStateSpace)
    #expect_true(Test_Of_Convergence(delta, beta, BudgetStateSpace, InitialGuess, Method = "Anderson")$Convergence) # This takes 17  iterations.
    #expect_true(Test_Of_Convergence(Method = "Simple")$Convergence)   # This takes 459 iterations and so will not be normally run.
    #expect_true(Test_Of_Convergence(Method = "Aitken")$Convergence)   # Does not converge.
    #expect_true(Test_Of_Convergence(Method = "Newton")$Convergence)   # Does not converge.
    #expect_true(Test_Of_Convergence(Method = "VEA")$Convergence)      # This takes 204 iterations and so will not be normally run.
    #expect_true(Test_Of_Convergence(Method = "SEA")$Convergence)      # Does not converge.
    #expect_true(Test_Of_Convergence(Method = "MPE")$Convergence)      # This takes 61  iterations and so will not be normally run.
    #expect_true(Test_Of_Convergence(Method = "RRE")$Convergence)      # This takes 97  iterations and so will not be normally run.
#})
#Func = OneIterateBudgetValues
#Inputs = InitialGuess
#sqm = SQUAREM::squarem(Inputs,Func, control=list(tol= 1e-10*4 ))
#convergence = sum(abs(sqm$par - Func(sqm$par) ))
#acat(paste0("The squarem method took ", sqm$fpevals, " iterations and finished with convergence of ", convergence, "\n")) # This takes 93 iterations and will not normally be run.

# Note this is a slow function so we only run Anderson below.
#test_that("Testing that each method converges in the quadratic case to within tolerance", {
  #delta = 0.2
  #beta = 0.99
  #BudgetStateSpace = c(seq(0,1, 0.015), seq(1.05,3,0.05))
  #InitialGuess = sqrt(BudgetStateSpace)
  #expect_true(Test_Of_Convergence(delta, beta, BudgetStateSpace, InitialGuess, Method = "Anderson")$Convergence) # This takes 71 iterations.
  #expect_true(Test_Of_Convergence(Method = "Simple")$Convergence)   # This takes 2316 iterations and so will not be normally run.
  #expect_true(Test_Of_Convergence(Method = "Aitken")$Convergence)   # Does not converge.
  #expect_true(Test_Of_Convergence(Method = "Newton")$Convergence)   # Does not converge.
  #expect_true(Test_Of_Convergence(Method = "VEA", MaxIter = 10000, PrintReports = TRUE)$Convergence)       # Does not converge.
  #expect_true(Test_Of_Convergence(Method = "SEA")$Convergence)      # Does not converge.
  #expect_true(Test_Of_Convergence(Method = "MPE")$Convergence)      # This takes 217  iterations and so will not be normally run.
  #expect_true(Test_Of_Convergence(Method = "RRE")$Convergence)      # This takes 159  iterations and so will not be normally run.
#})
#Func = OneIterateBudgetValues
#Inputs = InitialGuess
#sqm = SQUAREM::squarem(Inputs,Func, control=list(tol= 1e-10*4 ))
#convergence = sum(abs(sqm$par - Func(sqm$par) ))
#cat(paste0("The squarem method took ", sqm$fpevals, " iterations and finished with convergence of ", convergence, "\n")) # This takes 285 iterations.



