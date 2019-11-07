#BRTs Model Evaluation
#This script contains functions for EcoROMS model evaluations.
#Written by S.Brodie, H.Welch, and E.Becker

library(caret)
library(mlbench)

########
#K-folds evaluation
kfolds_eval <- function(dataInput, gbm.x, gbm.y, lr=lr, tc){
  DataInput <- dataInput
  DataInput$Kset <- kfold(DataInput,10) #randomly allocate k groups
  Evaluations_kfold <- as.data.frame(matrix(data=0,nrow=10,ncol=4))
  colnames(Evaluations_kfold) <- c("k","Deviance","AUC","TSS")
  counter=1
  for (k in 1:10){
    print(k)
    DataInput_train <- DataInput[DataInput$Kset!=k,]
    DataInput_test <- DataInput[DataInput$Kset==k,]
    DataInput.kfolds <- gbm.step(data=DataInput_train, gbm.x= gbm.x, gbm.y = gbm.y, 
                                 family="bernoulli", tree.complexity=tc,
                                 learning.rate = lr, bag.fraction = 0.6)
    preds <- predict.gbm(DataInput.kfolds, DataInput_test,
                         n.trees=DataInput.kfolds$gbm.call$best.trees, type="response")
    dev <- calc.deviance(obs=DataInput_test$PresAbs, pred=preds, calc.mean=TRUE)
    d <- cbind(DataInput_test$PresAbs, preds)
    pres <- d[d[,1]==1,2]
    abs <- d[d[,1]==0,2]
    e <- evaluate(p=pres, a=abs)
    Evaluations_kfold[counter,1] <- k
    Evaluations_kfold[counter,2] <- dev
    Evaluations_kfold[counter,3] <- e@auc
    Evaluations_kfold[counter,4] <- max(e@TPR + e@TNR-1)
    counter=counter+1 
  }
  return(Evaluations_kfold)}

#DEMO
# Species <- "SWOR"
# SWOR_kfold_res2 <- kfolds_eval(DGN_ROMS_SWOR, gbm.x.res2, "PresAbs", lr=0.01, tc=3)

########
#Leave One year Out analysis
LOO_eval <- function(DataInput, gbm.x, gbm.y, lr=lr, tc){
  DataInput$Year <- format(DataInput$dt, "%Y")
  Evaluations_LOO <- as.data.frame(matrix(data=0,nrow=1,ncol=4))
  colnames(Evaluations_LOO) <- c("k","Deviance","AUC","TSS")
  counter=1
  for (y in min(DataInput$Year):max(DataInput$Year)){
    print(y)
    DataInput_train <- DataInput[DataInput$Year!=y,]
    DataInput_test <- DataInput[DataInput$Year==y,]
    DataInput.loo <- gbm.step(data=DataInput_train, gbm.x= gbm.x, gbm.y = c("PresAbs"), 
                              family="bernoulli", tree.complexity=tc,
                              learning.rate = lr, bag.fraction = 0.6)
    preds <- predict.gbm(DataInput.loo, DataInput_test,
                         n.trees=DataInput.loo$gbm.call$best.trees, type="response")
    dev <- calc.deviance(obs=DataInput_test$PresAbs, pred=preds, calc.mean=TRUE)
    d <- cbind(DataInput_test$PresAbs, preds)
    pres <- d[d[,1]==1,2]
    abs <- d[d[,1]==0,2]
    if(length(pres)>0 & length(abs)>0){
    e <- evaluate(p=pres, a=abs)
    
    Evaluations_LOO[counter,1] <- y
    Evaluations_LOO[counter,2] <- dev
    Evaluations_LOO[counter,3] <- e@auc
    Evaluations_LOO[counter,4] <- max(e@TPR + e@TNR-1)
    counter=counter+1 
  }
  }
  return(Evaluations_LOO)}

# #DEMO
# Species <- "BLSH" 
# BLSH_kfold_res3 <- LOO_eval(DGN_ROMS_BLSH, gbm.x, "PresAbs", lr=0.03, tc=3)

#######
#Make function to 75/25 split AUC test. 
#This is to 'repeat' what Elliot thinks they did for EcoCast (Kylie disagrees, see next function)
eval_7525 <- function(dataInput, gbm.x, gbm.y, lr, tc=tc, family){
  DataInput <- dataInput
  Evaluations_7525 <- as.data.frame(matrix(data=0,nrow=1,ncol=3))
  colnames(Evaluations_7525) <- c("Deviance","AUC","TSS")
  DataInput_bound <- floor((nrow(DataInput)/4)*3)         #define % of training and test set
  DataInput_train<- DataInput[sample(nrow(DataInput),DataInput_bound),]
  DataInput_test<- DataInput[sample(nrow(DataInput),nrow(DataInput)-DataInput_bound),]
  DataInput.kfolds <- gbm.step(data=DataInput_train, gbm.x= gbm.x, gbm.y = gbm.y, 
                               family="bernoulli", tree.complexity=tc,
                               learning.rate = lr, bag.fraction = 0.6)
  preds <- predict.gbm(DataInput.kfolds, DataInput_test,
                       n.trees=DataInput.kfolds$gbm.call$best.trees, type="response")
  dev <- calc.deviance(obs=DataInput_test$PresAbs, pred=preds, calc.mean=TRUE)
  d <- cbind(DataInput_test$PresAbs, preds)
  pres <- d[d[,1]==1,2]
  abs <- d[d[,1]==0,2]
  e <- evaluate(p=pres, a=abs)
  Evaluations_7525[1,1] <- dev
  Evaluations_7525[1,2] <- e@auc
  Evaluations_7525[1,3] <- max(e@TPR + e@TNR-1)
  
  return(Evaluations_7525)}

###
#Make function to 100/100  AUC test
eval_100_percent <- function(dataInput, gbm.x, gbm.y, lr=lr, tc){
  DataInput <- dataInput
  Evaluations_100_percent <- as.data.frame(matrix(data=0,nrow=1,ncol=3))
  colnames(Evaluations_100_percent) <- c("Deviance","AUC","TSS")
  DataInput_train<- DataInput
  DataInput_test<- DataInput
  DataInput.kfolds <- gbm.step(data=DataInput_train, gbm.x= gbm.x, gbm.y = gbm.y, 
                               family="bernoulli", tree.complexity=tc,
                               learning.rate = lr, bag.fraction = 0.6)
  preds <- predict.gbm(DataInput.kfolds, DataInput_test,
                       n.trees=DataInput.kfolds$gbm.call$best.trees, type="response")
  dev <- calc.deviance(obs=DataInput_test$PresAbs, pred=preds, calc.mean=TRUE)
  d <- cbind(DataInput_test$PresAbs, preds)
  pres <- d[d[,1]==1,2]
  abs <- d[d[,1]==0,2]
  e <- evaluate(p=pres, a=abs)
  Evaluations_100_percent[1,1] <- dev
  Evaluations_100_percent[1,2] <- e@auc
  Evaluations_100_percent[1,3] <- max(e@TPR + e@TNR-1)
  
  return(Evaluations_100_percent)}

########
#Percent deviance explained: SINGLE model
#((null deviance - residual deviance)/null deviance)*100
dev_eval2=function(model_object){
  null <- model_object$self.statistics$mean.null
  res <- model_object$self.statistics$mean.resid
  dev=((null - res)/null)*100 
  return(dev)
}

########
#Elizabeths function to calculate ratio of observed to predicted values for study area 
ratio <- function(dataInput,model_object){
  ## Predict on model data using the best tree for predicting
  BRTpred <- predict.gbm(model_object, dataInput, n.trees = model_object$gbm.call$best.trees, "response")
  # calculate ratio of observed to predicted values for study area
  ratio.BRTpred <- sum(dataInput$PresAbs)/sum(BRTpred)
  return(ratio.BRTpred)
}

########
#Leave CA bight out analysis
# test how well model for ROMS extent predicts in CA bight
# bight extent = -121.5,-116, 32.5, 35
#Initially built by Briana for benioff, may have some usefulness for EcoROMS
eval_LBO <- function(dataInput, gbm.x, gbm.y, lr=lr, tc){
  DataInput <- dataInput
  Evaluations_LBO <- as.data.frame(matrix(data=0,nrow=1,ncol=3))
  colnames(Evaluations_LBO) <- c("Deviance","AUC","TSS")
  DataInput_train <- DataInput[DataInput$lon < -121.5 | DataInput$lat < 32.5 | DataInput$lat > 35,]
  DataInput_test <- DataInput[DataInput$lon > -121.5 & DataInput$lat > 32.5 & DataInput$lat < 35,]
  DataInput.kfolds <- gbm.step(data=DataInput_train, gbm.x= gbm.x, gbm.y = gbm.y, 
                               family="bernoulli", tree.complexity=tc,
                               learning.rate = lr, bag.fraction = 0.6)
  preds <- predict.gbm(DataInput.kfolds, DataInput_test,
                       n.trees=DataInput.kfolds$gbm.call$best.trees, type="response")
  dev <- calc.deviance(obs=DataInput_test$PresAbs, pred=preds, calc.mean=TRUE)
  d <- cbind(DataInput_test$PresAbs, preds)
  pres <- d[d[,1]==1,2]
  abs <- d[d[,1]==0,2]
  e <- evaluate(p=pres, a=abs)
  Evaluations_LBO[1,1] <- dev
  Evaluations_LBO[1,2] <- e@auc
  Evaluations_LBO[1,3] <- max(e@TPR + e@TNR-1)
  
  return(Evaluations_LBO)
}

########
#Get predicted values for a Poisson BRT, based on 75% 25% split
eval_7525_Poisson <- function(dataInput, gbm.x, gbm.y, lr=lr, tc=tc){
DataInput <- dataInput
DataInput_bound <- floor((nrow(DataInput)/4)*3)         #define % of training and test set
DataInput_train<- DataInput[sample(nrow(DataInput),DataInput_bound),]
DataInput_test<- DataInput[sample(nrow(DataInput),nrow(DataInput)-DataInput_bound),]
DataInput.kfolds <- gbm.step(data=DataInput_train, gbm.x= gbm.x, gbm.y = gbm.y, 
                             family="poisson", tree.complexity=tc,
                             learning.rate = lr, bag.fraction = 0.6)
preds <- predict.gbm(DataInput.kfolds, DataInput_test,
                     n.trees=DataInput.kfolds$gbm.call$best.trees, type="response")
d <- cbind(preds,DataInput_test$TotCat)
colnames(d) <- c("predicted", "observed")
d <- as.data.frame(d)
return(d)
}


#####
#Randomly split to 75% 25%
split_7525 <- function(dataInput){
  DataInput <- dataInput
  DataInput_bound <- floor((nrow(DataInput)/4)*3)   #define % of training and test set
  DataInput_train<- DataInput[sample(nrow(DataInput),DataInput_bound),]
  DataInput_test<- DataInput[sample(nrow(DataInput),nrow(DataInput)-DataInput_bound),]
  l <- list()
  l[[1]] <- DataInput_train
  l[[2]] <- DataInput_test
  return(l)
}



