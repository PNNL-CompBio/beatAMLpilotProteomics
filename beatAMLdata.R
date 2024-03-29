#'
#'Specialized set of commands that specifically load beatAML
#'data - it's pretty big
##just looad all the data for this project in a single file

library(amlresistancenetworks)
#'summarize dataset
#'@param auc.data - AUC data and clinical
#'@param mol.data -molelcular data
#'@import pheatmap
#'@import dplyr
#'@import gridExtra
plotAllPatients<-function(auc.data,pat.data,pphos){
  library(gridExtra)
  
  numDrugs = auc.data%>%
    group_by(`AML sample`)%>%
    summarize(numDrugs = n_distinct(Condition))
  
  pat.df<-pat.data%>%
    group_by(`AML sample`)%>%
    summarize(RNA=any(mRNALevels!=0),mutations=any(geneMutations!=0),
              proteins=any(proteinLevels!=0))%>%
    right_join(numDrugs)
  
  pat.df<-pphos%>%group_by(`Sample`)%>%summarize(phosphoSites=any(LogFoldChange!=0))%>%
    dplyr::rename(`AML sample`='Sample')%>%right_join(pat.df)
    
  pdf('patientSummaryTab.pdf',height=11)
  grid.table(pat.df)
  dev.off()
  return(pat.df)
  
}



#' loadBeatAMLMolecularData
#' @import dplyr
#' @import tidyr
loadBeatAMLMolecularData<-function(){
  require(tidyr)
  require(dplyr)
  print("loading molecular data")
  orig.data<-querySynapseTable('syn22172602')

  orig.data<-orig.data%>%dplyr::rename(proteinLevels='LogFoldChange')%>%
    dplyr::rename(mRNALevels='transcriptCounts')%>%
    dplyr::rename(geneMutations='Tumor VAF')%>%
    mutate(Gene=unlist(Gene))%>%rowwise()%>%
    mutate(binaryMutations=ifelse(geneMutations==0,0,1))
  
  pat.data<<-querySynapseTable("syn22314121")%>%
    mutate(Gene=unlist(Gene))%>%
    subset(Treatment=='Vehicle')%>%
    subset(`Cell number`>=10000000)%>%
    dplyr::select(Gene,LogFoldChange,`AML sample`)%>%
    distinct()%>%
    full_join(orig.data,by=c('AML sample','Gene'))%>%
    rowwise()%>%
    mutate(proteinLevels=max(proteinLevels,LogFoldChange,na.rm=T))%>%
    dplyr::select(-LogFoldChange)%>%
    mutate(mRNALevels=tidyr::replace_na(mRNALevels,0))%>%
    mutate(geneMutations=tidyr::replace_na(geneMutations,0))%>%
        mutate(binaryMutations=tidyr::replace_na(binaryMutations,0))%>%
    dplyr::select(-countMetric)%>%
   distinct()
  
  pats.with.prot<<-pat.data%>%
    group_by(`AML sample`)%>%
    summarize(hasProt=all(proteinLevels==0))%>%
    subset(hasProt==FALSE)%>%
    dplyr::select('AML sample')
  
  pat.data<<-pat.data%>%subset(`AML sample`%in%pats.with.prot$`AML sample`)
  print(paste('Have',length(unique(pat.data$`AML sample`)),'patients with proteomic data'))
  
  # ###replace values in original table....
   print("Getting site-corrected phosphosite data")
   p.pat.phos<-querySynapseTable("syn22156830")%>%
     mutate(site=unlist(site))%>%
     mutate(Gene=unlist(Gene))
   
   extra.phos<-querySynapseTable("syn22156814")%>%
     dplyr::select(Gene,site,Peptide,LogFoldChange='value',Sample="AML sample")%>%
    mutate(site=unlist(site))%>%
     mutate(Gene=unlist(Gene))
   
   soraf.phos<-querySynapseTable("syn22314122")%>%
     subset(Treatment=='Vehicle')%>%
     subset(`Cell number`>=10000000)%>%
     dplyr::select(Gene,site,Peptide,LogFoldChange,Sample="AML sample")%>%
     distinct()%>%
     mutate(site=unlist(site))%>%
     mutate(Gene=unlist(Gene))
  
  orig.pat.phos<<-rbind(p.pat.phos,unique(extra.phos),soraf.phos)%>%
    subset(Sample%in%pats.with.prot$`AML sample`)
  # 
  # ##getting kinase
  # print('Getting kinase estimates')
  # pat.kin <<-mapPhosphoToKinase(pat.phos)
  pat.phos<<-loadUnNormPhosData()%>%subset(Sample%in%pats.with.prot$`AML sample`)
  
}


loadUnNormPhosData<-function(){
  #loads unnormalized phosphoproteomics data
    require(tidyr)
    require(dplyr)
    print("Getting un-corrected phosphosite data")
  p.pat.phos<-querySynapseTable("syn24227903")%>%
    mutate(Gene=unlist(Gene))%>%
    mutate(site=unlist(site))
  

#We're still missing the un-normalized version of this data
  extra.phos<-querySynapseTable("syn24240355")%>%
    dplyr::select(Gene,site,Peptide,LogFoldChange='value',Sample="AML sample")%>%
    mutate(Gene=unlist(Gene))%>%
    mutate(site=unlist(site))

  soraf.phos<-querySynapseTable("syn24228075")%>%
    subset(Treatment=='Vehicle')%>%
    subset(`Cell number`>=10000000)%>%
    dplyr::select(Gene,site,Peptide,LogFoldChange,Sample="AML sample")%>%distinct()%>%
      mutate(Gene=unlist(Gene))%>%
    mutate(site=unlist(site))

  pat.phos<<-rbind(p.pat.phos,unique(extra.phos),soraf.phos)
  return(pat.phos)
}

#' looadBeatAMLClinicalDrugData
#' This function looads the clinical parameters and the
#' drug dosage data for each patient
#' @param threshold is the fraction of patients with an AUC under
#' 100, required to be included
#' @import dplyr
#' @import tidyr
loadBeatAMLClinicalDrugData<-function(threshold=0.10,famThreshold=0.05){
  print("Loading patient variables and drug response")

  require(dplyr)
  require(tidyr)
  
  mut.status<<-querySynapseTable("syn23538858")%>%
    mutate(status=tolower(status))%>%
    pivot_wider(values_from='status',names_from='variant')
  
  drug.class<<-querySynapseTable("syn22156956")%>%
    dplyr::rename(Condition='inhibitor')%>%
    mutate(Condition=unlist(Condition))%>%
    mutate(family=unlist(family))
  
  #drug response
  pat.drugClin<-querySynapseTable("syn22170540")%>%
    mutate(Condition=unlist(Condition))%>%
    left_join(drug.class,by='Condition')
  
  print("Fixing Vargatef mispelling")
  pat.drugClin$Condition<-sapply(pat.drugClin$Condition,function(x){
    ifelse(x=='Vargetef','Vargatef',x)})
  pat.drugClin$Condition<-sapply(pat.drugClin$Condition,function(x){
    ifelse(x=='Doramapimod','Doramapimod (BIRB 796)',x)})
  pat.drugClin$Condition<-sapply(pat.drugClin$Condition,function(x){
    ifelse(x=='Gilteritinib','Gilteritinib (ASP-2215)',x)})
  
  print("Reformating AUC data")
  clin.dat<-pat.drugClin%>%
    dplyr::select(`AML sample`,gender,ageAtDiagnosis,vitalStatus,overallSurvival,Condition,family)%>%
    distinct()
  
  auc.dat<- subset(pat.drugClin,Metric%in%c('auc','AUC'))%>%
    dplyr::select(`AML sample`,Condition,AUC='Value')%>%
    distinct()%>%
    group_by(`AML sample`,Condition)%>%
    mutate(meanAUC=mean(AUC))%>%
    ungroup()%>%
    mutate(Condition=unlist(Condition))%>%
    group_by(Condition)%>%
    mutate(medAUC=median(AUC))%>%
    mutate(percAUC=100*AUC/medAUC)%>%
    ungroup()%>%
    left_join(clin.dat)%>%
    dplyr::select(-AUC)%>%
    dplyr::rename(AUC='meanAUC')%>%
    left_join(mut.status)
  
    numSens<-auc.dat%>%
      group_by(Condition)%>%
      subset(AUC<100)%>%summarize(numSens=n())
    
    ##counting sensitivity by family
    numSensFam <-auc.dat%>%
      group_by(family)%>%
      subset(!is.na(family))%>%
      subset(AUC<100)%>%summarize(numSensFam=n())
    
     fracSens<-auc.dat%>%group_by(Condition)%>%
      summarize(nSamps=n())%>%
      left_join(numSens)%>%mutate(fracSens=numSens/nSamps)
 
     ##grouping by family to see if we get more drugs   
    fracSensFam<-auc.dat%>%
      subset(!is.na(family))%>%
      group_by(family)%>%
      summarize(nSampsFam=n())%>%
      left_join(numSensFam)%>%mutate(fracSens=numSensFam/nSampsFam)
    
    withSens=subset(fracSens,fracSens>threshold)%>%
      subset(fracSens<(1-threshold))%>%
      subset(numSens>1)
    include<-union(withSens$Condition,'Gilteritinib (ASP-2215)')
    auc.dat<-subset(auc.dat,Condition%in%include)
  
    drug.combos<-unique(auc.dat$Condition[grep(" - |\\+",auc.dat$Condition)])
    print("Removing drug combinations")
    auc.dat<<-subset(auc.dat,!Condition%in%drug.combos)
  
    ##add in family sensitivity
    withSensFam=subset(fracSensFam,fracSens>famThreshold)%>%
      subset(fracSens<(1-famThreshold))%>%
      subset(numSensFam>1)

    auc.dat.fam<<-subset(auc.dat,family%in%withSensFam$family)
    
#    drug.combos<-unique(auc.dat$Condition[grep(" - |\\+",auc.dat$Condition)])
#    print("Removing drug combinations")
#    auc.dat<<-subset(auc.dat,!Condition%in%drug.combos)
  
}
  
  
getNetworksAndLVs<-function(){
  library(dplyr)
  ##reduce dims
  print("Getting Latent Variables")
  lv.df<-querySynapseTable('syn22274890')
  
  ##now get network dsitances
  print('Getting network distances')
  #pat.net<-querySynapseTable("syn22343177")
  pat.net <-querySynapseTable("syn23448907")
  filtered = pat.net%>% 
    filter(net2_type=='community')%>%
    filter(hyp1=='patients')%>%
    filter(hyp2=='panCan')
  
#  filtered=pat.net%>%mutate(same=(`Hypha 1`==`Hypha 2`))%>%
#    filter(same)%>%
#    subset(net2_type=='community')
#  mut.nets<<-filtered%>%subset(`Hypha 1`=='mutations')%>%
#    dplyr::select(Community='Network 2', distance,`AML sample`='Network 1')%>%distinct()
  prot.nets<<-filtered%>%dplyr::select(Community='net2',distance,`AML sample`='net1')
    #filtered%>%subset(`Hypha 1`=='proteomics')%>%
    #dplyr::select(Community='Network 2', distance,`AML sample`='Network 1')%>%distinct()
}


#'loadBeatAMLData
#'General function that calls all beatAML data into memory, required for
#'analysis code to work
#'@export
loadBeatAMLData<-function(){
  require(dplyr)
  #getNetworksAndLVs()
  loadBeatAMLMolecularData()
  loadBeatAMLClinicalDrugData()
  auc.dat<<-auc.dat%>%subset(`AML sample`%in%pats.with.prot$`AML sample`)
  res<-plotAllPatients(auc.dat,pat.data,pat.phos)
  
  full.pats<<-res%>%
    rowwise()%>%
    mutate(fullData=(as.character(RNA)=="TRUE" && as.character(proteins)=="TRUE"
                     && as.character(mutations)=="TRUE" && as.character(phosphoSites)=="TRUE"))%>%
    subset(fullData==TRUE)%>%
    dplyr::select(`AML sample`)
  
}