#' A function for finding the fixed point of a contraction mapping
#'
#' This function takes in a function and an initial guess for the fixed point of that function. It then uses a fixed point method to find the
#' fixed point of that function.
#' @export
#' @param Function This is the function for which a fixed point is sought. This function must take and return a vector of the same dimension.
#' @param Inputs This can be either a vector of values that is an initial guess for a fixed point or it can be an N x A matrix of previous inputs
#' for which corresponding outputs are available. In this case N is the dimensionality of the fixed point vector you are seeking (Hence each
#' column is a matrix that is input to the "Function") and A is the number of previous Inputs/Outputs that are being provided to the fixed point.
#' Where a matrix is input, a corresponding outputs must be provided or the last column of the outputs matrix is taken as a startpoint guess and
#' the rest of the inputs and output matrices are discarded.
#' @param Outputs (Optional) This is a matrix of the Function values for each column of the input. It must be provided so that column k of the
#' outputs matrix is equal to Function(Column k of inputs matrix).
#' @param Method This is the fixed point method to be used. It can be "Anderson", "Simple", "Aitken", "Newton", "MPE", "RRE", "VEA" or "SEA". See
#' vignette and references to see explanations of these methods.
#' @param ConvergenceMetric This is a function that takes in a vector of residuals from one iterate of the function (defined as f(x) - x for
#' vector x and function f) and returns a scalar. This scalar should be low when convergence is close to being achieved. By default this is
#' the maximum residual by absolute value (the sup norm in the space of residuals).
#' @param ConvergenceMetricThreshold This is the threshold for terminating the algorithm. The algorithm will terminate when the scalar that
#' ConvergenceMetric returns is less than ConvergenceMetricThreshold. This can be set to a negative number in which case the algorithm will
#' run until MaxIter is hit or an error occurs (Note that an error is likely in trying to use any method other than "Simple" when a fixed point
#' is already found).
#' @param MaxIter This is the maximum number of iterates that will be undertaken.
#' @param MaxM This is the maximum number of saved iterates that are used in the Anderson algorithm. It has no effect if another method is chosen.
#' Note that the number of previous iterates that will actually be used is the minimum of MaxIter, the dimensionality of the function's vector and
#' the number of inputs that have been tried to previously (the width of the Outputs matrix at each given stage of the algorithm). If
#' PrintReports = TRUE, the number of previous iterates actually used is reported as the algorithm is running.
#' @param ExtrapolationPeriod This is the number of simple iterates to perform before extrapolating. This is used for the MPE, RRE, VEA and SEA
#' methods and has no effect if another method is chosen. Where an epsilon algorithm is used this should be one plus a multiple of two, ie (3,5,7,etc).
#' @param Dampening This is the dampening parameter. By default it is 1 which means no dampening takes place. It can also be less than 1
#' (indicating dampening) or more than 1 (indicating extrapolation).
#' @param PrintReports This is a boolean describing whether to print ongoing ConvergenceMetric values for each iterate.
#' @param ReportingSigFig This is the number of significant figures that will be used in printing the convergence values to the console
#' (only if PrintReports is TRUE).
#' @param ConditionNumberThreshold This is a threshold for what condition number is acceptable for solving the least squares problem for the
#' Anderson Method. If the condition number is larger than this threshold then fewer previous iterates will be used in solving the problem.
#' This has no effect unless the "Anderson" method is used.
#' @param Plot This determines whether a plot should be drawn for every iterate. It can be "NoPlot", "ConvergenceFig" or "ChangePerIterate".
#' By default it is "NoPlot" and no plot is drawn. If it is "ConvergenceFig" then a plot is shown with iterates on the x axis and convergence
#' (as defined by the ConvergenceMetric) is on the y axis. If it is "ChangePerIterate" then there is the index of the array value on the x axis
#' and the value of the array value on the y axis. The previous iterate is also shown so the change per iterate can be visualised.
#' @param ConvergenceFigLags This only affects anything if Plot == "ConvergenceFig". This gives how many previous iterates should be shown on
#' the x axis. By default it is 5. To see them all set it to a high number.
#' @param ChangePerIteratexaxis This only affects anything if Plot == "ChangePerIterate". Sometimes there is a more appropriate xaxis value
#' to use than (the default) value index for this figure. For instance in the consumption smoothing problem in the vignette every value is a
#' value function value at a given budget level. In this case the budget levels could be used for this xaxis.
#' @return A list containing the FixedPoint, the Inputs and corresponding Outputs, and convergence values (which are computed under the
#' "ConvergenceMetric"). The list will also include a "Finish" statement describing why it has finished. This is often going to be due to
#' either MaxIter or ConvergenceMetricThreshold being  reached. It may also terminate due to an error in generating a new input guess or using
#' the function with that guess. If this occurs the function will terminate early and the "Finish" statement will describe the issue.
#' In this event there will also be additional objects returned in the list "NewInputVector" and possibly "NewOutputVector" for use in debugging
#' the issue.
#' @examples
#' # For the simplest possible example we can seek the fixed point of the cos function with a scalar.
#' Inputs = 0.3
#' Function = function(x){ cos(x) }
#' A = FixedPoint(Function, Inputs, Method = "Aitken", Dampening = 0.5)
#' B = FixedPoint(Function, Inputs, Method = "Anderson", Dampening = 1.0)
#'
#' # For this next one the ConvergenceMetricThreshold is negative so the algorithm
#' # will keep running until MaxIter is met.
#' C = FixedPoint(Function, Inputs, Method = "Simple", MaxIter = 4, ConvergenceMetricThreshold = -1)
#' # This is not yet close to the fixed point but we can continue solving for this fixed point
#' # by inputting the previous inputs and outputs to the function. We can also switch methods
#' # and will switch below to the Newton Method.
#' D = FixedPoint(Function, C$Inputs, C$Outputs, Method = "Newton")
#'
#' # We can also find a 4 dimensional fixed point vector of this function.
#' Inputs = c(0.3, 98, 0, pi)
#' E = FixedPoint(Function, Inputs, Method = "Anderson")
#' F = FixedPoint(Function, Inputs, Method = "Anderson", MaxM = 4, ReportingSigFig = 13)
FixedPoint = function(Function, Inputs, Outputs = c(), Method = c("Anderson", "Simple", "Aitken", "Newton", "MPE", "RRE", "VEA", "SEA") ,
                      ConvergenceMetric  = function(Resids){max(abs(Resids))} , ConvergenceMetricThreshold = 1e-10, MaxIter = 1e3, MaxM = 10,
                      ExtrapolationPeriod = 7, Dampening = 1, PrintReports = FALSE, ReportingSigFig = 5, ConditionNumberThreshold = 1e3,
                      Plot = c("NoPlot", "ConvergenceFig", "ChangePerIterate"), ConvergenceFigLags = 5, ChangePerIteratexaxis = c()){
  # Checking inputs to the function
  Method = Method[1]
  Plot   = Plot[1]
  if (ConditionNumberThreshold < 1){stop("ConditionNumberThreshold must be at least 1.")}
  if (is.vector(Inputs)){
    Inputs  = matrix(Inputs, ncol = 1)
    Outputs = NULL
  }
  if (is.null(Outputs)){
    if (dim(Inputs)[2] > 1.5){
      warning("If you do not give outputs to the function then you can only give one vector of inputs. As you have input a matrix of input values
              everything but the last column has been discarded.\n")
      Inputs = Inputs[,dim(Inputs)[2]]
      }
  } else {
    if (sum(dim(Inputs) != dim(Outputs)) > 0.5) {
      warning("If you input a matrix of outputs as well as a matrix of inputs then inputs and outputs must be the same shape. As they differ in
               this case the last column of the inputs matrix has been taken as the starting point and everything else discarded.\n")
      Inputs  = Inputs[,dim(Inputs)[2]]
      Outputs = NULL
    }
  }
  SimpleStartIndex = dim(Inputs)[2]

  # Do an initial run if no runs have been done:
  if (is.null(Outputs)){
    Inputs          = matrix(Inputs, ncol = 1)
    NewOutputVector = try(Function(Inputs[,1]))
    # In the event of an error we return informative information on the error.
    if (class(NewOutputVector) == "try-error"){return(list(Inputs = Inputs, Outputs = NA, Convergence = NA, FixedPoint = NA,
                                                           Finish = "Could not execute function with input vector"))}
    if (sum(is.na(NewOutputVector))){return(list(Inputs = Inputs, Outputs = NA, Convergence = NA, FixedPoint = NA,
                                                 Finish = "New output vector contains NAs", NewInputVector = Inputs,
                                                 NewOutputVector = NewOutputVector))}
    # If no error we put the results into the Outputs matrix.
    Outputs = matrix(Function(Inputs[,1]), ncol = 1)
  } else {
    # This ensures that MaxIter refers to max iter excluding any previous passed in results
    MaxIter = MaxIter + dim(Outputs)[2]
    # This is to take into account previously passed in simple iterations (without jumps).
    SimpleStartIndex = SimpleStartIndex- (dim(PutTogetherIteratesWithoutJumps(Inputs, Outputs))[2] )
  }

  # First running through the last column of Inputs to test if we already have a fixed point.
  Resid = Outputs - Inputs
  iter = dim(Resid)[2]
  ConvergenceVector = sapply(1:iter, function(x) ConvergenceMetric(Resid[,x]) )
  if (ConvergenceVector[iter] < ConvergenceMetricThreshold){
    if (PrintReports){cat("The last column of Inputs matrix is already a fixed point under input convergence metric and convergence threshold")}
    return(list(FixedPoint = Outputs[,iter] , Inputs = Inputs, Outputs = Outputs, Convergence = ConvergenceVector))
  }
  # Printing a report for initial convergence
  Convergence = ConvergenceVector[iter]
  if (PrintReports){cat(paste0(format(" ", width = 49, justify = "right"), "Method: ", format(Method, width = 8, justify = "right")   ,
                               ". Iteration: ", format(iter, digits =  0, width = 5,scientific = FALSE),
                               ". Convergence: ", NicePrint(Convergence, ReportingSigFig), "\n"))}
  if (Plot == "ConvergenceFig"){ConvergenceFig(Inputs, Outputs,  Input_Convergence = ConvergenceVector, FromIterate = 1)}
  if (Plot == "ChangePerIterate"){ChangePerIterate(Inputs, Outputs, ConvergenceVector, FromIterate = dim(Inputs)[2],
                                                   ToIterate = dim(Inputs)[2], xaxis = ChangePerIteratexaxis, secondhold = -1)}
  iter = iter + 1

  while ((Convergence > ConvergenceMetricThreshold) & (iter <= MaxIter)){
    # Generating new input and output.
    NewInputVector = try({FixedPointNewInput(Inputs = Inputs, Outputs = Outputs, Method = Method, MaxM = MaxM,
                                              SimpleStartIndex = SimpleStartIndex, ExtrapolationPeriod = ExtrapolationPeriod,
                                              Dampening = Dampening, ConditionNumberThreshold = ConditionNumberThreshold,
                                              PrintReports = PrintReports)})
    if (class(NewInputVector) == "try-error"){return(list(Inputs = Inputs, Outputs = Outputs, Convergence = ConvergenceVector,
                                                          FixedPoint = NA, Finish = "Could not generate new input vector"))}
    if (sum(is.na(NewInputVector))>0.5){return(list(Inputs = Inputs, Outputs = Outputs, Convergence = ConvergenceVector,
                                                FixedPoint = NA, Finish = "New input vector contains NAs", NewInputVector = NewInputVector))}
    if (sum(is.infinite(NewInputVector))>0.5){return(list(Inputs = Inputs, Outputs = Outputs, Convergence = ConvergenceVector,
                                                    FixedPoint = NA, Finish = "New input vector contains Infs", NewInputVector = NewInputVector))}
    if (Method != "Anderson" & PrintReports){cat(paste0(format(" ", width = 49, justify = "right")))}

    NewOutputVector = try({Function(NewInputVector)})
    if (class(NewOutputVector) == "try-error"){return(list(Inputs = Inputs, Outputs = Outputs, Convergence = ConvergenceVector,
                                                           FixedPoint = NA, Finish = "Could not execute function with generated vector",
                                                           NewInputVector = NewInputVector))}
    if (sum(is.na(NewOutputVector))){return(list(Inputs = Inputs, Outputs = Outputs, Convergence = ConvergenceVector,
                                                 FixedPoint = NA, Finish = "New output vector contains NAs", NewInputVector = NewInputVector,
                                                 NewOutputVector = NewOutputVector))}

    Inputs  = matrix(c(Inputs, NewInputVector), ncol = iter, byrow = FALSE)
    Outputs = matrix(c(Outputs, NewOutputVector), ncol = iter, byrow = FALSE)
    Resid   = matrix(c(Resid, NewOutputVector - NewInputVector), ncol = iter, byrow = FALSE)
    # Checking and recording convergence
    ConvergenceVector =  c(ConvergenceVector, ConvergenceMetric(Resid[,iter]) )
    Convergence = ConvergenceVector[iter]
    # Output of report and going to next iteration.
    if (PrintReports){cat(paste0("Method: ", format(Method, width = 8, justify = "right")   , ". Iteration: ",
                                 format(iter, digits =  0, width = 5,scientific = FALSE), ". Convergence: ",
                                 NicePrint(Convergence, ReportingSigFig), "\n"))}
    if (Plot == "ConvergenceFig"){ConvergenceFig(Inputs, Outputs,  Input_Convergence =  ConvergenceVector,
                                                 FromIterate = max(1, dim(Inputs)[2]- ConvergenceFigLags))}
    if (Plot == "ChangePerIterate"){ChangePerIterate(Inputs, Outputs, ConvergenceVector, secondhold = -1,
                                                     FromIterate = dim(Inputs)[2], ToIterate = dim(Inputs)[2],
                                                     xaxis = ChangePerIteratexaxis)}
    iter  = iter + 1
  }
  FixedPoint = Outputs[,dim(Outputs)[2]]

  if (Convergence < ConvergenceMetricThreshold){Finish = "Reached Convergence Threshold"} else {Finish = "Reached MaxIter"}
  return(list(Inputs = Inputs, Outputs = Outputs, Convergence = ConvergenceVector, FixedPoint = FixedPoint, Finish = Finish))
}

#' FixedPointNewInput
#' This function takes the previous inputs and outputs from the FixedPoint function and determines what vector to try next in seeking
#' a fixed point.
#' @export
#' @param Inputs This is an N x A matrix of previous inputs for which corresponding outputs are available. In this case N is the
#' dimensionality of the fixed point vector that is being sought (Hence each column is a matrix that is input to the "Function") and A
#' is the number of previous Inputs/Outputs that are being provided to the fixed point.
#' @param Outputs This is a matrix of Function values for the each column of the "Inputs" matrix.
#' @param Method This is the fixed point method to be used. It can be "Anderson", "Simple", "Aitken", "Newton", "MPE", "RRE", "VEA", "SEA".
#' @param MaxM This is the number of saved iterates that are used in the Anderson algorithm. This has no role if another method is used.
#' @param SimpleStartIndex This is the index for what column in the input/output matrices did the algorithm start doing simple iterates
#' without jumps. This is used for all methods except the simple and Anderson methods where it has no effect.
#' @param ExtrapolationPeriod This is the number of simple iterates to perform before extrapolating. This is used for the MPE, RRE, VEA and SEA
#' methods and has no effect if another method is chosen. Where an epsilon algorithm is used this should be one plus a multiple of two, ie (3,5,7,etc).
#' @param Dampening This is the dampening parameter. By default it is 1 which means no dampening takes place. It can also be less than 1
#' (indicating dampening) or more than 1 (indicating extrapolation).
#' @param ConditionNumberThreshold This is what threshold should be chosen to drop previous iterates if the matrix is ill conditioned.
#' Only used in Anderson acceleration.
#' @param PrintReports This is a boolean describing whether to print ongoing ConvergenceMetric values for each iterate.
#' @return A vector containing the next guess in seeking a fixed point.
#' @examples
#' FPFunction = function(x){c(0.5*sqrt(abs(x[1] + x[2])), 1.5*x[1] + 0.5*x[2])}
#' A = FixedPoint( Function = FPFunction, Inputs = c(0.3,900), MaxIter = 6, Method = "Simple")
#' NewGuessAnderson = FixedPointNewInput(A$Inputs, A$Outputs, Method = "Anderson")
#' NewGuessVEA = FixedPointNewInput(A$Inputs, A$Outputs, Method = "VEA")
#' NewGuessMPE = FixedPointNewInput(A$Inputs, A$Outputs, Method = "MPE")
#' NewGuessAitken = FixedPointNewInput(A$Inputs, A$Outputs, Method = "Aitken")
FixedPointNewInput = function(Inputs, Outputs, Method = "Anderson", MaxM = 10, SimpleStartIndex = 1,
                               ExtrapolationPeriod = 7, Dampening = 1, ConditionNumberThreshold = 1e3, PrintReports = FALSE ) {

  CompletedIters = dim(Outputs)[2]

  if (Method == "Simple"){
    return( as.vector(Outputs[,CompletedIters] ))
    }
  if (Method == "Anderson"){
    if (CompletedIters < 1.5){
      if (PrintReports){cat(paste0(format(" ", width = 32, justify = "right"), "  Using",  format(0, width = 3, justify = "right")," lags. "))}
      return(as.vector( Outputs[,CompletedIters] ))
      }
    VectorLength   = dim(Outputs)[1]
    M = min(MaxM-1,CompletedIters-1,VectorLength)
    Coeffs = NA

    Outputs         = matrix(Outputs[, seq(CompletedIters-M, CompletedIters,1)], ncol = M+1)
    Inputs          = matrix(Inputs[ , seq(CompletedIters-M, CompletedIters,1)], ncol = M+1)
    Resid           = Outputs - Inputs
    DeltaOutputs    = matrix(Outputs[,2:(M+1)] - Outputs[,1:M], ncol = M)
    DeltaResids     = matrix(Resid[  ,2:(M+1)] - Resid[  ,1:M], ncol = M)

    LastResid       = Resid[,M+1]
    LastOutput      = Outputs[,M+1]

    while (sum(is.na(Coeffs))>0.5){
      ConditionNumber = rcond(DeltaResids)
      if (ConditionNumber > ConditionNumberThreshold){
        M = M-1
        DeltaOutputs= matrix(DeltaOutputs[, 2:(M+1)], ncol = M)
        DeltaResids = matrix(DeltaResids[,2:(M+1)], ncol = M)
        next
      }

      Fit             = stats::lm.fit(x = DeltaResids, y = LastResid, tol = 1e-13)
      Coeffs          = Fit$coefficients

      if (sum(is.na(Coeffs))>0.5){
      M = M-1
      if (M < 1.5){
        # This happens occasionally in test cases where the iteration is very close to a fixed point.
        if (PrintReports){cat(paste0(format(" ", width = 32, justify = "right"), "  Using",
                                     format(0, width = 3, justify = "right")," lags. "))}
        return( LastOutput )
      }
      DeltaOutputs= matrix(DeltaOutputs[, 2:(M+1)], ncol = M)
      DeltaResids = matrix(DeltaResids[,2:(M+1)], ncol = M)
      }
    }
    if (PrintReports){cat(paste0("Condition number is ", format(format(ConditionNumber, digits = 5,scientific = TRUE), width = 12, justify = "right")
                                 ,". Used:",  format(M+1, width = 3, justify = "right"), " lags. "))}
    NewGuess        = LastOutput - Dampening * Coeffs  %*% t(DeltaOutputs)
    return( as.vector(NewGuess) )
  }
  if (Method == "Aitken"){
    if ((CompletedIters + SimpleStartIndex) %% 3 == 0){
      # If we are in 3rd, 6th, 9th, 12th iterate from when we started Acceleration then we want to do a jumped Iterate,
      # First we extract the guess that started this run of 3 iterates (x), the Function applied to it (fx) and the function applied to that (ffx)
      x = Inputs[,CompletedIters -1]
      fx = Outputs[,CompletedIters -1]
      ffx = Outputs[,CompletedIters ]
      # Now using the appropriate formula to make a new guess. Note that if a vector is input here it is used elementwise.
      NewGuess = x - ((fx - x)^2/(ffx - 2*fx + x))
      # Now there is the chance that the demoninator is zero which results in an Inf.
      if (sum(is.infinite(NewGuess)) + sum(is.na(NewGuess)) > 0.5){
         NewGuess[(is.infinite(NewGuess) | is.na(NewGuess))] = Outputs[(is.infinite(NewGuess) | is.na(NewGuess)) ,CompletedIters]
        }
      return( as.vector(Dampening*NewGuess + (1-Dampening)*Outputs[,CompletedIters] ))
    } else {
      # We just do a simple iterate. We do an attempt with the latest iterate.
      return( as.vector(Outputs[,CompletedIters] ))
      }
  }
  if (Method == "Newton"){
    if ((CompletedIters + SimpleStartIndex) %% 2 == 1 & CompletedIters > 1){
      # If we are in 3rd, 6th, 9th, 12th iterate from when we started Newton Acceleration then we want to do a Newton Iterate,
      # First we extract the guess that started this run of 3 iterates (xk1), the Function applied to it (fxk1, xk) and the function applied to that (fxk)
      xk1 = Inputs[,CompletedIters -1]
      fxk1 = Outputs[,CompletedIters -1]
      gxk1 = fxk1 - xk1
      xk  = Inputs[,CompletedIters ]
      fxk = Outputs[,CompletedIters ]
      gxk = fxk - xk
      # Now using the appropriate formula to make a new guess. Note that if a vector is input here it is used elementwise.
      derivative = (gxk-gxk1)/(xk - xk1)
      NewGuess = xk - (gxk /derivative)
      # Now checking and replacing elements where the derivative was zero.
      if (sum(is.infinite(NewGuess)) + sum(is.na(NewGuess)) > 0.5){
        NewGuess[(is.infinite(derivative) | is.na(derivative))] = Outputs[(is.infinite(derivative) | is.na(derivative)) ,CompletedIters]
      }
      return(as.vector(Dampening*NewGuess + (1-Dampening)*Outputs[,CompletedIters] ))
    } else {
      # We just do a simple iterate.
      return(as.vector( Outputs[,CompletedIters] ))
    }
  }
  if (Method == "MPE" | Method == "RRE"){
    SimpleIteratesMatrix = PutTogetherIteratesWithoutJumps(Inputs, Outputs)
    if ( dim(SimpleIteratesMatrix)[2] %% ExtrapolationPeriod == 0  ){
        NewGuess = PolynomialExtrapolation(SimpleIteratesMatrix, Method = Method)
        return(as.vector(Dampening*NewGuess + (1-Dampening)*Outputs[,CompletedIters]))
    } else {
      # We just do a simple iterate.
      return( as.vector(Outputs[,CompletedIters]) )
    }
  }
  if (Method == "VEA" | Method == "SEA"){
    SimpleIteratesMatrix = PutTogetherIteratesWithoutJumps(Inputs, Outputs)
    if ( dim(SimpleIteratesMatrix)[2] %% ExtrapolationPeriod == 0  ){
      NewGuess = EpsilonExtrapolation(SimpleIteratesMatrix, Method = Method)
      return(as.vector(Dampening*NewGuess + (1-Dampening)*Outputs[,CompletedIters]))
    } else {
      # We just do a simple iterate.
      return( as.vector(Outputs[,CompletedIters] ))
    }
  }
}


#' PolynomialExtrapolation
#' This function performs Minimal Polynomial extrapolation (MPE) or Reduced Rank Extrapolation (RRE) given a matrix of previous iterates of the
#' function.
#' @export
#' @import MASS
#' @param Iterates A matrix of inputs and outputs excluding jumps. Can be pieced together from Inputs and Outputs matrices of the FixedPoint
#' function using the PutTogetherIteratesWithoutJumps function.
#' @param Method The method for polynomial extrapolation. Should be either "MPE" for minimal polynomial extrapolation or "RRE" for reduced
#' rank extrapolation.
#' @return A vector containing the extrapolated vector.
#' @examples
#' FPFunction = function(x){c(0.5*sqrt(abs(x[1] + x[2])), 1.5*x[1] + 0.5*x[2])}
#' A = FixedPoint( Function = FPFunction, Inputs = c(0.3,900), MaxIter = 6, Method = "Simple")
#' Iterates = PutTogetherIteratesWithoutJumps(A$Inputs, A$Outputs)
#' PolynomialExtrapolation(Iterates, "MPE")
#' B = FixedPoint( function(x){cos(x)}, Inputs = 1, MaxIter = 5, Method = "Simple")
#' Iterates = PutTogetherIteratesWithoutJumps(B$Inputs, B$Outputs)
#' PolynomialExtrapolation(Iterates, "RRE")
PolynomialExtrapolation = function(Iterates, Method = c("MPE", "RRE")){
  Method = Method[1]
  if (!(Method %in% c("MPE", "RRE"))){stop("Invalid method input. PolynomialExtrapolation function can only take Method as MPE or RRE.")}
  if (Method == "MPE"){
    TotalColumnsOfIterates = dim(Iterates)[2]
    Differences          = matrix(Iterates[,2:(TotalColumnsOfIterates-1)] - Iterates[,1:(TotalColumnsOfIterates-2)],
                                   ncol = TotalColumnsOfIterates-2)
    LastDifference       = matrix(  (Iterates[,TotalColumnsOfIterates]-Iterates[,(TotalColumnsOfIterates-1)])     , ncol = 1)
    InverseDifferences   = MASS::ginv(Differences)
    cVector               = -InverseDifferences %*% LastDifference
    cVector               = rbind(cVector,1)
    s                    = (Iterates[,2:TotalColumnsOfIterates]%*%cVector)/sum(cVector)
    return(s)
  }
  if (Method == "RRE"){
    TotalColumnsOfIterates = dim(Iterates)[2]
    FirstColumn          = matrix(Iterates[,1], ncol = 1)
    Differences          = matrix(Iterates[,2:(TotalColumnsOfIterates)] - Iterates[,1:(TotalColumnsOfIterates-1)],
                                   ncol = TotalColumnsOfIterates-1)
    SecondDifferences    = matrix(Differences[,2:(TotalColumnsOfIterates-1)] - Differences[,1:(TotalColumnsOfIterates-2)],
                                   ncol = TotalColumnsOfIterates-2)
    FirstDifference      = matrix(Differences[,1], ncol = 1)
    Differences          = matrix(Differences[, 1:(TotalColumnsOfIterates-2)], ncol =  TotalColumnsOfIterates-2)
    InverseSecondDifferences = MASS::ginv(SecondDifferences)
    s = FirstColumn - Differences %*% InverseSecondDifferences %*% FirstDifference
    return(s)
  }
}


#' EpsilonExtrapolation
#' This function takes a matrix with previous iterates and extrapolates the limit of the sequence.
#' @export
#' @import MASS
#' @param Iterates A matrix representing different iterates with one iterate per column. Can be pieced together from Inputs and Outputs matrices
#' of the FixedPoint function using the PutTogetherIteratesWithoutJumps function.
#' @param Method Method for epsilon extrapolation. Should be either "VEA" for the vector extrapolation algorithm or "SEA" for the scalar epsilon
#' algorithm.
#' @return A vector with the extrapolated vector.
#' @examples
#' FPFunction = function(x){c(0.5*sqrt(abs(x[1] + x[2])), 1.5*x[1] + 0.5*x[2])}
#' A = FixedPoint( Function = FPFunction, Inputs = c(0.3,900), MaxIter = 6, Method = "Simple")
#' Iterates = PutTogetherIteratesWithoutJumps(A$Inputs, A$Outputs)
#' EpsilonExtrapolation(Iterates, "VEA")
#' B = FixedPoint( function(x){cos(x)}, Inputs = 1, MaxIter = 5, Method = "Simple")
#' Iterates = PutTogetherIteratesWithoutJumps(B$Inputs, B$Outputs)
#' EpsilonExtrapolation(Iterates, "SEA")
EpsilonExtrapolation = function(Iterates, Method = c("VEA", "SEA")){
    if (dim(Iterates)[2] %% 2 == 0){Iterates = matrix(Iterates[,2:dim(Iterates)[2]], ncol = dim(Iterates)[2]-1 )}
    # The function cannot do anything to a one dimensional input so will return input unchanged.
    if (dim(Iterates)[2] == 1){return(Iterates)}
    Method = Method[1]
    if (!(Method %in% c("VEA", "SEA"))){stop("Invalid method input. EpsilonExtrapolation function can only take Method as VEA or SEA")}

    Matrix = Iterates
    RowsOfMatrix    = dim(Matrix)[1]
    TotalColumnsOfMatrix = dim(Matrix)[2]
    PreviousMatrix = matrix(rep(0, RowsOfMatrix*(TotalColumnsOfMatrix-1)), nrow = RowsOfMatrix )

    for (MatrixColumns in TotalColumnsOfMatrix:(2)){
      NewMatrix = PreviousMatrix + EpsilonExtrapolationVectorOfInverses(matrix(Matrix[,2:MatrixColumns] - Matrix[,1:(MatrixColumns-1)],
                                                                                nrow = RowsOfMatrix), Method = Method)
      PreviousMatrix = Matrix[,2:(MatrixColumns-1)]
      Matrix = NewMatrix
    }

    # The function can get NAs from the inversion (ie if differenceMatrix contains a zero for SEA then there is division by zero).
    # To avert this we try with 2 less columns.
    if (any(is.na(Matrix))){Matrix = EpsilonExtrapolation( matrix(Iterates[,3:dim(Iterates)[2]], ncol = dim(Iterates)[2]-2 )  ,Method)    }
  return(matrix(Matrix, nrow = RowsOfMatrix))
}

#' EpsilonExtrapolationVectorOfInverses
#' This is a helper function for EpsilonExtrapolation
#' @param DifferenceMatrix The matrix of the differences in elements to be inverted.
#' @param Method "SEA" or "VEA".
#' @return A vector of the result of inverting each (column) vector in a mmatrix.
EpsilonExtrapolationVectorOfInverses = function(DifferenceMatrix, Method){
  if (dim(DifferenceMatrix)[1] == 1 | Method == "SEA"){return(1/DifferenceMatrix)
    } else {
  return(apply(DifferenceMatrix, 2, function(x) MASS::ginv(x)))
    }
}



#' PutTogetherIteratesWithoutJumps
#' This function takes the previous inputs and outputs and assembles a matrix with both excluding jumps.
#' @export
#' @param Inputs This is an N x A matrix of previous inputs for which corresponding outputs are available. In this case N is the
#' dimensionality of the fixed point vector that is being sought (and each column is a matrix that is input to the "Function") and A is the
#' number of previous Inputs/Outputs that are being provided to the fixed point.
#' @param Outputs This is a matrix of "Function" values for each column of the "Inputs" matrix.
#' @param AgreementThreshold A parameter for determining when a column in Inputs and a column in Outputs match. They are deemed to match if
#' the sum of the absolute values of the difference in the columns is less than AgreementThreshold.
#' @return A matrix of inputs and outputs excluding jumps.
#' @examples
#'A = FixedPoint( function(x){c(0.5*sqrt(abs(x[1] + x[2])), 1.5*x[1] + 0.5*x[2])},
#' Inputs = c(0.3,900), MaxIter = 5, Method = "Simple")
#'A = FixedPoint( function(x){c(0.5*sqrt(abs(x[1] + x[2])), 1.5*x[1] + 0.5*x[2])},
#' Inputs = A$Inputs, Outputs = A$Outputs, MaxIter = 5, Method = "Aitken")
#'A = FixedPoint( function(x){c(0.5*sqrt(abs(x[1] + x[2])), 1.5*x[1] + 0.5*x[2])},
#' Inputs = A$Inputs, Outputs = A$Outputs, MaxIter = 5, Method = "Simple")
#'CombinedSequence = PutTogetherIteratesWithoutJumps(A$Inputs, A$Outputs)
PutTogetherIteratesWithoutJumps = function(Inputs, Outputs, AgreementThreshold = 1e-10){
  if (!(is.matrix(Inputs))){stop("Inputs is not a matrix.")}
  if (!(is.matrix(Outputs))){stop("Outputs is not a matrix.")}
  if (any(dim(Inputs) != dim(Outputs))){stop("Inputs and Outputs matrices are not comformable.")}
  Dimensions = dim(Inputs)

  if (Dimensions[2] == 1){return(matrix(cbind(matrix(Inputs, ncol = 1), matrix(Outputs, ncol = 1)), ncol =2))}
  Difference = matrix(Inputs[,2:(Dimensions[2])] - Outputs[,1:(Dimensions[2]-1)], nrow = Dimensions[1])

  Sum_Of_Differences = apply(Difference, 2, function(x) {sum(abs(x))})
  Agreements = Sum_Of_Differences < AgreementThreshold
  if (all(Agreements)){return(cbind(matrix( Inputs[,1:(Dimensions[2])], nrow=Dimensions[1]), matrix(Outputs[,Dimensions[2]],nrow=Dimensions[1])))}

  LocationsOfBreaks = which(!Agreements)
  LastBreak = LocationsOfBreaks[length(LocationsOfBreaks)]
  return(matrix(cbind(matrix(Inputs[,(LastBreak+1):(Dimensions[2])], ncol = Dimensions[2] - LastBreak),
                      matrix(Outputs[, Dimensions[2]], ncol = 1)), ncol = Dimensions[2] - LastBreak + 1))
}


#' NicePrint
#'
#' This is a small function to ensure that decimal places line up for convergence values reported by the fixedpoint function.
#' @param Number The number to be formatted.
#' @param sigfig The number of significant figures to be taken.
#' @return A nicely formatted string version of the input number for printing to the console.
NicePrint = function(Number, sigfig = 5){
  CPrint    = formatC(Number,  big.interval = 1, small.interval = sigfig, format = "e")
  return(CPrint)
}
