### model eval functions
### written by HW 2019; modified from code written by HW, SB, BA

library(caret)
library(mlbench)
library(gbm)
library(dismo)

### BRTS ####
## These functions call the gbm.fixed() function - differs from gbm.step in that is fits a predefined number of trees.
# dataInput: complete.cases input data with predictors. y term needs to be "presabs"
# gbm.x = character list of predictor names
# gbm. y = response term, must be "presabs"
# lr = learning rate
# tc = tree complexity
# bf = bag fraction
# nt = number of trees

kfolds_eval_brt <- function(dataInput, gbm.x, gbm.y, lr, tc, bf,nt){
  DataInput <- dataInput
  DataInput$Kset <- kfold(DataInput,10) #randomly allocate k groups
  Evaluations_kfold <- as.data.frame(matrix(data=0,nrow=10,ncol=3))
  colnames(Evaluations_kfold) <- c("k","AUC","TSS")
  counter=1
  for (k in 1:10){
    print(k)
    DataInput_train <- DataInput[DataInput$Kset!=k,]
    DataInput_test <- DataInput[DataInput$Kset==k,]
    DataInput.kfolds <- gbm.fixed(data=DataInput_train, gbm.x= gbm.x, gbm.y = gbm.y, 
                                 family="bernoulli", tree.complexity=tc,
                                 learning.rate = lr, bag.fraction = bf,n.trees = nt)
    preds <- predict.gbm(DataInput.kfolds, DataInput_test,
                         n.trees=nt, type="response")
    dev <- calc.deviance(obs=DataInput_test$presabs, pred=preds, calc.mean=TRUE)
    d <- cbind(DataInput_test$presabs, preds)
    pres <- d[d[,1]==1,2]
    abs <- d[d[,1]==0,2]
    e <- evaluate(p=pres, a=abs)
    Evaluations_kfold[counter,1] <- k
    Evaluations_kfold[counter,2] <- e@auc
    Evaluations_kfold[counter,3] <- max(e@TPR + e@TNR-1)
    counter=counter+1 
  }
  return(Evaluations_kfold)
}

LOO_eval_brt <- function(DataInput, gbm.x, gbm.y, lr, tc, bf,nt){
  Evaluations_LOO <- as.data.frame(matrix(data=0,nrow=1,ncol=3))
  colnames(Evaluations_LOO) <- c("k","AUC","TSS")
  counter=1
  for (y in min(DataInput$Year):max(DataInput$Year)){
    print(y)
    DataInput_train <- DataInput[DataInput$Year!=y,]
    DataInput_test <- DataInput[DataInput$Year==y,]
    DataInput.loo <- gbm.fixed(data=DataInput_train, gbm.x= gbm.x, gbm.y = gbm.y, 
                              family="bernoulli", tree.complexity=tc,
                              learning.rate = lr, bag.fraction = bf,n.trees = nt)
    preds <- predict.gbm(DataInput.loo, DataInput_test,
                         n.trees=nt, type="response")
    dev <- calc.deviance(obs=DataInput_test$presabs, pred=preds, calc.mean=TRUE)
    d <- cbind(DataInput_test$presabs, preds)
    pres <- d[d[,1]==1,2]
    abs <- d[d[,1]==0,2]
    if(length(pres)>0 & length(abs)>0){
      e <- evaluate(p=pres, a=abs)
      
      Evaluations_LOO[counter,1] <- y
      Evaluations_LOO[counter,2] <- e@auc
      Evaluations_LOO[counter,3] <- max(e@TPR + e@TNR-1)
      counter=counter+1 
    }
  }
  return(Evaluations_LOO)}

eval_7525_brt <- function(dataInput, gbm.x, gbm.y, lr, tc, bf,nt){
  DataInput <- dataInput
  Evaluations_7525 <- as.data.frame(matrix(data=0,nrow=1,ncol=2))
  colnames(Evaluations_7525) <- c("AUC","TSS")
  DataInput_bound <- floor((nrow(DataInput)/4)*3)         #define % of training and test set
  DataInput_train<- DataInput[sample(nrow(DataInput),DataInput_bound),]
  DataInput_test<- DataInput[sample(nrow(DataInput),nrow(DataInput)-DataInput_bound),]
  DataInput.kfolds <- gbm.fixed(data=DataInput_train, gbm.x= gbm.x, gbm.y = gbm.y, 
                               family="bernoulli", tree.complexity=tc,
                               learning.rate = lr, bag.fraction = bf,n.trees = nt)
  preds <- predict.gbm(DataInput.kfolds, DataInput_test,
                       n.trees=nt, type="response")
  dev <- calc.deviance(obs=DataInput_test$presabs, pred=preds, calc.mean=TRUE)
  d <- cbind(DataInput_test$presabs, preds)
  pres <- d[d[,1]==1,2]
  abs <- d[d[,1]==0,2]
  e <- evaluate(p=pres, a=abs)
  Evaluations_7525[1,1] <- e@auc
  Evaluations_7525[1,2] <- max(e@TPR + e@TNR-1)
  
  return(Evaluations_7525)}

eval_100_brt <- function(dataInput, gbm.x, gbm.y, lr, tc, bf,nt){
  DataInput <- dataInput
  Evaluations_100_percent <- as.data.frame(matrix(data=0,nrow=1,ncol=2))
  colnames(Evaluations_100_percent) <- c("AUC","TSS")
  DataInput_train<- DataInput
  DataInput_test<- DataInput
  DataInput.kfolds <- gbm.fixed(data=DataInput_train, gbm.x= gbm.x, gbm.y = gbm.y, 
                               family="bernoulli", tree.complexity=tc,
                               learning.rate = lr, bag.fraction =bf,n.trees = nt)
  preds <- predict.gbm(DataInput.kfolds, DataInput_test,
                       n.trees=nt, type="response")
  dev <- calc.deviance(obs=DataInput_test$presabs, pred=preds, calc.mean=TRUE)
  d <- cbind(DataInput_test$presabs, preds)
  pres <- d[d[,1]==1,2]
  abs <- d[d[,1]==0,2]
  e <- evaluate(p=pres, a=abs)
  Evaluations_100_percent[1,1] <- e@auc
  Evaluations_100_percent[1,2] <- max(e@TPR + e@TNR-1)
  
  return(Evaluations_100_percent)}

dev_eval_brt=function(model_object){
  null <- model_object$self.statistics$null.deviance
  res <- model_object$self.statistics$resid.deviance
  dev=((null - res)/null)*100 
  return(dev)
}

ratio_brt <- function(dataInput,model_object,nt){
  ## Predict on model data using the best tree for predicting
  dataInput$RN=1:nrow(dataInput)
  BRTpred <- predict.gbm(model_object, dataInput, n.trees = nt, "response")
  # calculate ratio of observed to predicted values for study area
  ratio.BRTpred <- sum(dataInput$presabs)/sum(BRTpred)
  return(ratio.BRTpred)
}


### GAMM and GLMM ####

eval_100_gamm_glmm <- function(dataInput, formula) {
  DataInput=dataInput
  Evaluations_100_gamm <- as.data.frame(matrix(data=0,nrow=1,ncol=2))
  colnames(Evaluations_100_gamm) <- c("AUC","TSS")
  DataInput_train<- DataInput
  DataInput_test<- DataInput
  DataInput.100 <- mgcv::gamm(formula=formula,data=DataInput_train,random=list(tag=~1), family=binomial,niterPQL=50)
  preds<-as.data.frame(mgcv::predict.gam(DataInput.100$gam, DataInput_test, se=TRUE, type="response"))
  d <- cbind(DataInput_test$presabs, preds)
  pres <- as.numeric(d[d[,1]==1,2])
  abs <- as.numeric(d[d[,1]==0,2])
  e <- evaluate(p=pres, a=abs)
  Evaluations_100_gamm[1,1] <- e@auc
  Evaluations_100_gamm[1,2] <- max(e@TPR + e@TNR-1)
  return(Evaluations_100_gamm)
} 

eval_kfold_gamm_glmm <- function(dataInput, formula){
  DataInput=dataInput
  DataInput$Kset <- kfold(DataInput,10) #randomly allocate k groups
  Evaluations_kfold_gamm <- as.data.frame(matrix(data=0,nrow=10,ncol=3))
  colnames(Evaluations_kfold_gamm) <- c("k","AUC","TSS")
  counter=1
  for (k in 1:10){
    print(k)
    DataInput_train <- DataInput[DataInput$Kset!=k,]
    DataInput_test <- DataInput[DataInput$Kset==k,]
    DataInput.kfolds <- mgcv::gamm(formula=formula,data=DataInput_train,random=list(tag=~1), family=binomial,niterPQL=50)
    preds <- as.data.frame(mgcv::predict.gam(DataInput.kfolds$gam, DataInput_test,se=TRUE, type="response"))
    d <- cbind(DataInput_test$presabs, preds)
    pres <- as.numeric(d[d[,1]==1,2])
    abs <- as.numeric(d[d[,1]==0,2])
    e <- evaluate(p=pres, a=abs)
    Evaluations_kfold_gamm[counter,1] <- k
    Evaluations_kfold_gamm[counter,2] <- e@auc
    Evaluations_kfold_gamm[counter,3] <- max(e@TPR + e@TNR-1)
    counter=counter+1 
  }
  return(Evaluations_kfold_gamm)}

eval_LOO_gamm_glmm <- function(dataInput, formula){
  DataInput=dataInput
  Evaluations_LOO_gamm <- as.data.frame(matrix(data=0,nrow=1,ncol=3))
  colnames(Evaluations_LOO_gamm) <- c("k","AUC","TSS")
  counter=1
  for (y in unique(DataInput$Year)){
    DataInput_train <- DataInput[DataInput$Year!=y,]
    DataInput_test <- DataInput[DataInput$Year==y,]
    print(y)
    DataInput.loo <- mgcv::gamm(formula=formula,data=DataInput_train,random=list(tag=~1), family=binomial,niterPQL=50)
    preds <- as.data.frame(mgcv::predict.gam(DataInput.loo$gam, DataInput_test,se=TRUE, type="response"))
    d <- cbind(DataInput_test$presabs, preds)
    pres <- as.numeric(d[d[,1]==1,2])
    abs <- as.numeric(d[d[,1]==0,2])
    e <- evaluate(p=pres, a=abs)
    Evaluations_LOO_gamm[counter,1] <- y
    Evaluations_LOO_gamm[counter,2] <- e@auc
    Evaluations_LOO_gamm[counter,3] <- max(e@TPR + e@TNR-1)
    counter=counter+1 
  }
  return(Evaluations_LOO_gamm)}

eval_7525_gamm_glmm <- function(dataInput, formula){
  Evaluations_7525 <- as.data.frame(matrix(data=0,nrow=1,ncol=2))
  colnames(Evaluations_7525) <- c("AUC","TSS")
  dataInput_bound <- floor((nrow(dataInput)/4)*3)         #define % of training and test set
  dataInput_train<- dataInput[sample(nrow(dataInput),dataInput_bound),]
  dataInput_test<- dataInput[sample(nrow(dataInput),nrow(dataInput)-dataInput_bound),]
  dataInput.7525 <- mgcv::gamm(formula=formula,data=dataInput_train,random=list(tag=~1), family=binomial,niterPQL=50)
  preds<-as.data.frame(mgcv::predict.gam(dataInput.7525$gam, dataInput_test, se=TRUE, type="response"))
  dev <- calc.deviance(obs=dataInput_test$presabs, pred=preds$fit, calc.mean=TRUE)
  d <- cbind(dataInput_test$presabs, preds$fit)
  pres <- as.numeric(d[d[,1]==1,2])
  abs <- as.numeric(d[d[,1]==0,2])
  e <- evaluate(p=pres, a=abs)
  Evaluations_7525[1,1] <- e@auc
  Evaluations_7525[1,2] <- max(e@TPR + e@TNR-1)
  
  return(Evaluations_7525)}

rsq_gamm_glmm=function(model_object){ ### still lost
   rsq=summary(glm$gam)$r.sq
  return(rsq)
}

ratio_gamm_glmm <- function(dataInput,model_object){
  ## Predict on model data using the best tree for predicting
  dataInput$RN=1:nrow(dataInput)
  GMpred <- as.data.frame(mgcv::predict.gam(dataInput.7525$gam, dataInput_test, se=TRUE, type="response"))
  # calculate ratio of observed to predicted values for study area
  ratio.GMpred <- sum(dataInput$presabs)/sum(GMpred$fit)
  return(ratio.GMpred)
}


# ####testing brts ####
# setwd("~/Dropbox/Pseudoabsence_Data")
# load("bwhalemodels.RData")
# brt=bwhaleBRT.lr005.back
# dat=readRDS("bwhaledata_buff_test.RDS") %>% .[complete.cases(.),] %>% .[1:1000,]
# 
# gbm.x <- c("Bathymetry","Oxygen_100m", "SST", "Chla_25km_monthly", 
#            "SLA", "Oxygen_Surface", "SST_SD", "Chla_4km_8day", 
#            "MLD", "Rugosity","U", "FSLE_max", "eke", "V",
#            "Theta_max")
# 
# a=kfolds_eval_brt(dataInput = dat,gbm.x = gbm.x,gbm.y = "presabs",lr=0.005,tc = 5,bf=.75,nt=nt)
# a1=LOO_eval_brt(DataInput=dat,gbm.x = gbm.x,gbm.y = "presabs",lr=0.005,tc = 5,bf=.75,nt=nt)
# a2=eval_7525_brt(dataInput=dat,gbm.x = gbm.x,gbm.y = "presabs",lr=0.005,tc = 5,bf=.75,nt=nt)
# a3=eval_100_brt(dataInput=dat,gbm.x = gbm.x,gbm.y = "presabs",lr=0.005,tc = 5,bf=.75,nt=nt)
# b=dev_eval_brt(brt)
# c=ratio_brt(dataInput = dat,model_object = brt,nt=2000)

# #### testing gamms ####
# gam=bwhalegamm_buff
# dat=readRDS("bwhaledata_buff_test.RDS") %>% .[complete.cases(.),] %>% .[1:1000,]
# formula=presabs~s(SST,bs="ts",k=3)+s(MLD,bs="ts",k=3)+s(Oxygen_100m,bs="ts",k=3)+s(Theta_max,bs="ts", k=3)+s(log10(Chla_4km_8day+0.001),k=3,bs="ts")+s(log10(Chla_25km_monthly+0.001),k=3,bs="ts")+s(SSH,k=3,bs="ts")+s(SLA,k=3,bs="ts")+s(Bathymetry,k=3,bs="ts")+s(Rugosity,k=3,bs="ts")
# a=eval_100_gamm(dataInput = dat,formula = formula)
# b=eval_kfold_gamm(dataInput = dat,formula = formula)
# c=eval_LOO_gamm(dataInput = dat,formula = formula)
# 
# ### testing glms ###
# glm=bwhaleGLMM_buff
# formula=presabs~SST + SST_SD + MLD + Oxygen_100m + FSLE_max + Chla_4km_8day + Chla_25km_monthly + SSH + SLA + Bathymetry + Rugosity
# a=eval_100_gamm_glmm(dataInput = dat,formula = formula)
# b=eval_kfold_gamm_glmm(dataInput = dat,formula = formula)
# c=eval_LOO_gamm_glmm(dataInput = dat,formula = formula)
# d=eval_7525_gamm_glmm(dataInput = dat,formula = formula)
# e=ratio_gamm_glmm(dataInput=dat,model_object = glm)
# e=rsq_gamm_glmm(model_object = glm)
