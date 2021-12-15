#' @title Make Prediction Using the SM-DL Model
#' @description The function can make prediction given the test data based on the SM-DNN model.
#' @param testMat A genotype matrix (N x M; N individuals, M markers) for prediction.
#' @param SMDNN_model A list contains local and global networks obtaining from SMDNN function.
#' @param subp A constant for splitting the features. It indicates how many features each subset contains after splitting the data e.g. the data originally contains 2000 features. By setting subp=500, the function split the orginal data into 4 subsets with each subset contains 500 features.
#' @param localtype (String)  This parameter indicates what networks you would like to use for local networks. The default setting for type is CNN-convolutional neural network or FNN-Feed-forward Neural Network, respectively.
#' @author Wei-Heng Huang


SMpred <- function(SMDNN_model, testMat, subp, localtype = 'CNN'){
  #extract the last hidden layer from the local networks
  hidden_test <- c()
  local_num <- length(SMDNN_model) - 1
  for(i in 1:local_num){
    if(local_num == 1){
      testMat_sub <- testMat
    }else{
      print(paste0("Retrieving the hidden layer of Local Network: ", i))
      #Split the Features
      if(i != local_num){
        testMat_sub <- testMat[,((i-1)*subp+1):(i*subp)]
      }else{
        testMat_sub <- testMat[,((i-1)*subp+1):dim(testMat)[2]]
      }
    }
    testMat_sub <- as.matrix(testMat_sub)

    localnet <- SMDNN_model[[i]]
    if(localtype == 'CNN'){
      
      dim(testMat_sub) <- c(1, ncol(testMat_sub),1,nrow(testMat_sub))
      #For the case that the local nets are CNN
      internals <- localnet$symbol$get.internals()
      internal_num <- length(internals$outputs)
      fea_symbol <- internals[[internal_num-7]]

      para_num <- length(localnet$arg.params)
      localnet$arg.params[[para_num]] <- NULL
      localnet$arg.params[[para_num-1]] <- NULL

      tmpmodel <- list(symbol = fea_symbol,
                     arg.params = localnet$arg.params,
                     aux.params = localnet$aux.params)
      class(tmpmodel) <- "MXFeedForwardModel"

      hidden <- predict(tmpmodel, testMat_sub)

      hidden_test <- rbind(hidden_test, hidden)
    } else{
      #For the case that the local nets are FNN
      internals <- localnet$symbol$get.internals()
      internal_num <- length(internals$outputs)
      fea_symbol <- internals[[internal_num-7]]

      para_num <- length(localnet$arg.params)
      localnet$arg.params[[para_num]] <- NULL
      localnet$arg.params[[para_num-1]] <- NULL

      tmpmodel <- list(symbol = fea_symbol,
                       arg.params = localnet$arg.params,
                       aux.params = localnet$aux.params)
      class(tmpmodel) <- "MXFeedForwardModel"

      hidden <- predict(tmpmodel, testMat_sub)

      hidden_test <- rbind(hidden_test, hidden)
    }

  }
  hidden_test <- data.matrix(t(hidden_test))

  pred_test <- predict(SMDNN_model[[length(SMDNN_model)]], hidden_test)

  return(pred_test)

}
