##Build and store regression models

##builds all models and stores them on synapse in rds files

## Now we can build the models

#We want to build both the LASSO and logistic regression models using the following code. This code will divide the patients into 5 batches and train/test the model holding out 1/5 of the data each time. Then we will have predictions on each of the samples. 

library(amlresistancenetworks)
source("beatAMLdata.R")

if(!exists('dataLoaded')){
  loadBeatAMLData()
  #loadUnNormPhosData()
  dataLoaded=TRUE
}


kFoldCrossValModels<-function(drugFamily=FALSE){
  all.pats <- intersect(pat.data$`AML sample`,auc.dat$`AML sample`)
  library(caret)
  folds<-createFolds(all.pats,k=5)
  ##separate out data?
  
  ###single drug or family?
  cond='Condition'
  if(drugFamily)
    cond='family'
  
  res<-lapply(folds,function(testpats){
    test.pats<-all.pats[testpats]
    train.pats<-setdiff(all.pats,test.pats)
    
    print('Now getting gene-based preds')
    
    tr.dat<-subset(pat.data,`AML sample`%in%train.pats)
    te.dat<-subset(pat.data,`AML sample`%in%test.pats)
    
    ##This list describes the features to be used  
    eval.list<-list(
      combo=list(mf=c('Gene','Gene','Gene','site'),
                 fn=c('mRNALevels','proteinLevels','binaryMutations','Phosphosite')),
      protCombo=list(mf=c("Gene","site"),
                     fn=c("proteinLevels","Phosphosite")),
      mRNA=c('Gene','mRNALevels'),
      prot=c("Gene","proteinLevels"),
      mut=c('Gene','binaryMutations'),
      phosph=c('site','Phosphosite'))
    
    message("Comparing logistic preds")
    log.reg.preds<-purrr::map_df(eval.list,
                                 ~ drugMolLogRegEval(auc.dat,
                                                     tr.dat,
                                                     mol.feature=.x[1],
                                                     mol.feature.name=.x[2],
                                                     dplyr::rename(auc.dat,Sample='AML sample'),
                                                     dplyr::rename(te.dat,Sample='AML sample'),
                                                     category=cond))
    
    
    message('Running lasso')
    ##now train model on AML and eval on depmap data
    reg.preds<-purrr::map_df(eval.list,
                             ~ amlresistancenetworks::drugMolRegressionEval(auc.dat,
                                                                            tr.dat,
                                                                            mol.feature=.x[1],
                                                                            mol.feature.name=.x[2],                                                   
                                                                            dplyr::rename(auc.dat,Sample='AML sample'),
                                                                            dplyr::rename(te.dat,Sample='AML sample'),
                                                                            category=cond))
    message("Running elastic net")
    enet.reg.preds<-purrr::map_df(eval.list,
                                  ~ drugMolRegressionEval(auc.dat,
                                                          tr.dat,
                                                          mol.feature=.x[1],
                                                          mol.feature.name=.x[2],
                                                          dplyr::rename(auc.dat,Sample='AML sample'),
                                                          dplyr::rename(te.dat,Sample='AML sample'),
                                                          category=cond,doEnet=TRUE))
    
    
    
    
    e.results<-enet.reg.preds%>%
      mutate(method='ElasticNet')
    
    full.results<-reg.preds%>%
      mutate(method='LASSO')
    
    lr.results<-log.reg.preds%>%
      mutate(method='LogisticReg')%>%
      mutate(MSE=MSE*10000)
    

    all.res<-rbind(full.results,lr.results,e.results)
    print(dim(all.res))
    return(all.res)
    
  })
  
  res<-do.call(rbind,res)
  
  fname='combinedKfoldRes.rds'
  if(drugFamily)
    fname='combinedKfoldResFamily.rds'
  saveRDS(res,fname)
  synapseStore(fname,'syn26529370')
  
}


fullyTrainedModel<-function(drugFamily=FALSE){
  
    ###next build full model
    tr.dat<- pat.data%>%
      left_join(rename(pat.phos,`AML sample`='Sample',Phosphosite='LogFoldChange'))
  
    cond='Condition'
    if(drugFamily)
      cond='family'
  
  
  
    eval.list<-list(
      combo=list(mf=c('Gene','Gene','Gene','site'),
                   fn=c('mRNALevels','proteinLevels','binaryMutations','Phosphosite')),
      protCombo=list(mf=c("Gene","site"),
                       fn=c("proteinLevels","Phosphosite")),
        
      mRNA=c('Gene','mRNALevels'),
      prot=c("Gene","proteinLevels"),
      mut=c('Gene','binaryMutations'),
      phosph=c('site','Phosphosite'))
    
    
     message('Running lasso')
    ##now train model on AML and eval on depmap data
     reg.preds<-purrr::map_df(eval.list,
                             ~ amlresistancenetworks::drugMolRegression(auc.dat,
                                                                        tr.dat,
                                                                        mol.feature=.x[1],
                                                                        mol.feature.name=.x[2],                                                   
                                                                        category=cond))
    
    
    message("Comparing logistic preds")
    
    log.reg.preds<-purrr::map_df(eval.list,
                                 ~ drugMolLogReg(auc.dat,
                                                 tr.dat,
                                                 mol.feature=.x[1],
                                                 mol.feature.name=.x[2],
                                                 category=cond))
    
    message("Running elastic net")
    
    
    enet.reg.preds<-purrr::map_df(eval.list,
                                  ~ drugMolRegression(auc.dat,
                                                      tr.dat,
                                                      mol.feature=.x[1],
                                                      mol.feature.name=.x[2],
                                                      category=cond,doEnet=TRUE))
    
    
    
    e.results<-enet.reg.preds%>%
      mutate(method='ElasticNet')
    
    full.results<-reg.preds%>%
      mutate(method='LASSO')
    
    lr.results<-log.reg.preds%>%
      mutate(method='LogisticReg')%>%
      mutate(MSE=MSE*10000)
    
    #full.results<-rbind(full.results,lr.results)
    #saveRDS(full.results,'fulllassoRegPreds.rds')
    #saveRDS(lr.results,'fulllogRegPreds.rds')
    #saveRDS(e.results,'fullenetRegPreds.rds')
    #  saveRDS(full.results,'mostlyCompletePredictions.rds')
    
    all.preds<-rbind(full.results,lr.results,e.results)
    #all.preds<-rbind(enet.preds,log.preds,reg.preds)
    fname='combinedPreds.rds'
    if(drugFamily)
      fname='combinedPredsFamily.rds'
    saveRDS(res,fname)
    synapseStore(fname,'syn26529370')
}
fullyTrainedModel(TRUE)
fullyTrainedModel(FALSE)
kFoldCrossValModels(TRUE)
kFoldCrossValModels(FALSE)