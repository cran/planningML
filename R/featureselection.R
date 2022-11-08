#' Feature selection
#'
#' This function selects important features from the dataset
#'
#' @param x a matrix of predictor variables
#' @param y a vector of binary outcome
#' @param method feature selection method, default is iHCT
#'
#' @return \code{featureselection()} returns selected features and other outcomes needed for sample size determination.
#'
#' @import tidyverse caret lubridate Matrix
#'
#' @examples
#'
#' ## load data
#' pilot.data = readRDS(system.file("extdata", "pilotdata_sub.rds", package = "planningML"))
#' x = pilot.data[,-ncol(pilot.data)]
#' y = pilot.data$DEPRESSION
#'
#' ## select important features
#' features = featureselection(x = x, y = y)
#' summary(features)
#'
#' @export
#'

featureselection = function(x=NULL, y=NULL, method = "iHCT"){

  # check contents of input dataset
  if (is.null(x) | is.null(y)){
    warning("The input dataset is not complete.")
  }

  # preprocess data
  # assume users give pilot data or historical data they have

  # pilot.data = x
  # pilot.data[,ncol(pilot.data)+1] = y
  pilot.data = data.frame(x, y)

  ## remove columns with all zeros
  pilot.data.filt<-pilot.data[vapply(pilot.data, function(x) length(unique(x)) > 1, logical(1L))]

  #Scale pilot data (has 9315 covariates)
  scaled.pilot.data<- data.frame( apply(pilot.data.filt, 2,scale))
  scaled.pilot.data[,ncol(scaled.pilot.data)]<-y

  ## remove features with no variation or whose values are almost constant for all subjects
  sd_check<-which(is.na(apply(pilot.data[,-ncol(pilot.data)],2,sd))== TRUE)
  if (length(sd_check) != 0){
    pilot.data<-pilot.data[,-sd_check]
  }

  # Calculate correlation matrix
  pilot.data.corr<-cor(as.matrix(pilot.data[,-ncol(pilot.data)]))
  pilot.data.corr<-nearPD(pilot.data.corr)
  pilot.data.corr<-as.matrix(pilot.data.corr$mat)

  # From Innovated HCT paper by Hall and Jin
  # Create cholesky decomposition
  pilot.data.corr.chol=chol(pilot.data.corr)

  #Check whether u %*%Sigma %*% u = I where u = solve(t(pilot.data.selected.cov.chol))
  A = solve(t(pilot.data.corr.chol))%*%pilot.data.corr%*%solve(pilot.data.corr.chol)
  u = solve(t(pilot.data.corr.chol))

  n = nrow(pilot.data)
  b <- log(nrow(pilot.data))

  utilde <- u
  for (k in 1:n){
    for (j in 1:n){
      if (k - b +1 <= j && j <= k){
        utilde[k,j] <- utilde[k,j]
      }
      else {
        utilde[k,j] <- 0
      }
    }
  }
  #Define function for normalizing columns
  normalize.cols<-function(mat){
    mat_squared = apply(mat,2, function(x) x^2)
    mat_squared_sum = sqrt(colSums(mat_squared))
    #mapply(function(x,y) apply(x,2, function(z) z/y),mat,mat_squared_sum)
    mat_normalized<-mat
    for(i in 1: ncol(mat))
    {
      mat_normalized[,i] = mat[,i]/mat_squared_sum[i]
    }
    return(mat_normalized)
  }

  ubar <- normalize.cols(utilde)
  v <- t(ubar)%*%u
  vx <- as.matrix(pilot.data[,-ncol(pilot.data)])%*%v #check i
  pvalue_vx <- 1 - pnorm(vx,0,1)


  #Function to calculate z-scores and obtain p-values
  z_score<-function(mat)
  {
    z_score<-rep(NA,ncol(mat)-1)
    n1 = sum(mat[,ncol(mat)])
    n2 = sum(1-mat[,ncol(mat)])
    matr = mat[,-ncol(mat)]
    sd0 = apply(matr,2,var)
    sd = sd0 + as.numeric(sd0==0)
    #sd1 = apply(mat[which(mat[,ncol(mat)] == 1),],2,var)
    #sd2 = apply(mat[which(mat[,ncol(mat)] == 0),],2,var)
    z_score<- (1/sqrt((sd/n1)+(sd/n2)))*(apply(matr[which(mat[,ncol(mat)] == 1),],2,sum)/n1 - apply(matr[which(mat[,ncol(mat)] == 0),],2,sum)/n2)
    standardized_z_score <- (z_score - mean(z_score,na.rm = TRUE)) /sd(z_score,na.rm = TRUE)
    #standardized_z_score <- z_score
    #pvalue_standardized_z_score<- 2*pnorm( abs(standardized_z_score ),0, 1, lower.tail = FALSE)
    pvalue_standardized_z_score<- 1- pnorm(standardized_z_score,0,1)
    return(data.frame(standardized_z_score,pvalue_standardized_z_score))
  }
  #Get the z-score and p-values
  vx_new<-as.matrix(data.frame(vx,pilot.data[,ncol(pilot.data)]))
  z_score_result<-z_score(vx_new)

  #convert NAN values to INF
  z_score_result$pvalue_standardized_z_score[is.nan(z_score_result$pvalue_standardized_z_score)] = Inf;

  #Sort z_score_result
  z_score_result_sort <- z_score_result[order( z_score_result$pvalue_standardized_z_score),]

  #Function to calculate hc-scores

  in_hc_score<-function(sorted_p_values){
    N<-length(sorted_p_values)
    hc_score<-vector(length = N)
    for(i in 1: N)
    {
      hc_score[i]<-sqrt(N)*((i/N - sorted_p_values[i])/ sqrt( (sorted_p_values[i] * (1 - sorted_p_values[i]))))
    }
    return(hc_score)
  }

  #Calculate the HC scores
  hc_score_pilot_data<-in_hc_score(z_score_result_sort$pvalue_standardized_z_score)
  names(hc_score_pilot_data)<-rownames(z_score_result_sort)

  # Choose alpha_0*N p-value sort
  p<- ncol(pilot.data)-1
  alpha=0.1 # Results don't change for higher values
  pilot.data.hc_score_alpha<-abs(hc_score_pilot_data)[1: floor(alpha*p)]
  ind<-(1/(sqrt((2*b) -1)))*(which.max(pilot.data.hc_score_alpha))

  #Obtain threshold value
  threshold_choice<-abs(z_score_result_sort[names(ind),]$standardized_z_score)
  #Obtain z-score at threshold choice
  zscore_selected_features<-z_score_result_sort[which(z_score_result_sort$standardized_z_score>threshold_choice),]

  selected_features<-rownames(zscore_selected_features)

  #Indices of the selected features
  index = match(selected_features,colnames(pilot.data))

  output = list(x = x,
                y = y,
                data = pilot.data,
                index =  index,
                features = selected_features)


  class(output)<-"planningML"
  return(output)


}
