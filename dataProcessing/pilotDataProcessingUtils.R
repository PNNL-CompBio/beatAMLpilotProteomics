##Scripts to show howw we processed data for the pilot
##################################BEATAML PATIENT samPles##################
library(amlresistancenetworks)

getPatientTranscript<-function(patientlist){
  library(dplyr)
  library(readxl)
  syn=synapseLogin()
  #we dont need the RPKM because we have the CPM
  exp.3.megafile=syn$get('syn22130786')$path
  patient.cpm<-readxl::read_xlsx(exp.3.megafile,sheet='Table S9-Gene Counts CPM')%>%
    tidyr::pivot_longer(-c(Gene,Symbol),names_to='patient',values_to='transcriptCounts')%>%
    mutate(countMetric='CPM')%>%
    select(-Gene)%>%
    rename(Gene='Symbol',`AML sample`='patient')

  subset(patient.cpm, `AML sample`%in%patientlist)

}

getPatientVariants<-function(patientlist){
  library(dplyr)
  library(readxl)
  syn=synapseLogin()
  exp.3.megafile=syn$get('syn22130786')$path
  gene.var<-readxl::read_xlsx(exp.3.megafile,sheet='Table S7-Variants for Analysis')%>%
    subset(labId%in%patientlist)%>%
    select(labId,t_vaf,symbol)%>%distinct()%>%
    #tidyr::pivot_longer(c(t_vaf,n_vaf),names_to='Metric',values_to='Value')%>%
    rename(`AML sample`='labId',`Tumor VAF`='t_vaf',Gene='symbol')
  return(gene.var)
}

#'add in wave3/4 genetics
#'@import readxl
#'@import dplyr
addWave34GeneticsToTable<-function(){
  syn=synapseLogin()
  newtab<-readxl::read_xlsx(syn$get('syn22130786')$path,sheet=9)%>%
    dplyr::select(Gene=symbol,`AML sample`=original_id,`Tumor VAF`=t_vaf)%>%distinct()
  newtab

}


#' this grabs drug responses from three different sources, all on synapse
#' @param patientlist it will only return drugs for these patients
#' @return data frame of drug response values
#' @import dplyr
#' @import readxl
getPatientDrugResponses<-function(patientlist){
  library(dplyr)
  library(readxl)
  syn=synapseLogin()
    exp.3.megafile=syn$get('syn22130786')$path

    #this file has the primary drug responses
  dose.response<-readxl::read_xlsx(exp.3.megafile,sheet='Table S10-Drug Responses')%>%
    tidyr::pivot_longer(c(ic50,auc),names_to='Metric',values_to='Value')%>%
    dplyr::rename(`AML sample`='lab_id',Condition='inhibitor')

  #here is an additional set of data
  other.data<-readAndTidySensMetadata()%>%
    dplyr::select(-Barcode)

  #even more data
  extra.files<-c('syn22170222','syn22170223','syn22170225','syn22170226','syn22170227')
  more.data<-purrr::map_df(extra.files,function(x){
    fname=syn$get(x)$path
    if(length(grep('xlsx',fname))>0)
      ex<-readxl::read_xlsx(fname)
    else
      ex<-readxl::read_xls(fname)

    ex%>%dplyr::select(`AML sample`='Specimen: Lab ID',
                    Condition='Inhibitor Panel Definition: Drug',
                    AUC='Probit Interpretation: Area Under Curve',
                    IC50='Probit Interpretation: IC50')%>%
      tidyr::pivot_longer(c(IC50,AUC),names_to='Metric',values_to='Value')
  })
  comb.response<-rbind(dose.response,other.data,more.data)
  return(subset(comb.response,`AML sample`%in%patientlist))

}



#' getCuratedMutstatus
#' @export
#' @import dplyr
#' @import readxl
#' @import tidyr
getCuratedMutStatus<-function(){
  library(dplyr)
  syn<-synapseLogin()
  new.file<-syn$get('syn23538805')$path
#  exp.3.megafile=syn$get('syn22130786')$path
  
  patData<-readxl::read_xlsx(new.file,sheet='wv1to4')%>%
    dplyr::select(`AML sample`='labId','FLT3-ITD',NPM1,`FLT3-MUT`='FLT3')%>%distinct()%>%
      tidyr::pivot_longer(-`AML sample`,names_to='variant',values_to='status')%>%
      tidyr::replace_na(list(status='not_tested'))

  
# patData<-readxl::read_xlsx(exp.3.megafile,sheet='Clinical Summary')%>%
#  subset(labId%in%unlist(patients))%>%
#  select(`AML sample`='labId',contains("FLT3"))%>%
#  distinct()%>%
#   pivot_longer(-`AML sample`,values_to='FLT3_Status',names_to='FLT3_Assay')%>%
#   separate(FLT3_Status,into=c("FLT3_Status","FLT3_Details"),sep=' ',extra='merge')
  # mutate(Value=tidyr::replace_na(Value,0))
 synTableStore(patData,'BeatAML Pilot Mutation Status')

}
 
#' getPatientMetadata
#' @export
#' @import dplyr
#' @import readxl
getPatientMetadata<-function(synid='syn22170540'){
    require(dplyr)
#    synapser::synLogin()
  syn<-synapseLogin()
    exp.3.megafile=syn$get('syn22130786')$path

  patients<-readxl::read_xlsx(exp.3.megafile,sheet='Sample Summary')%>%
    dplyr::select('Specimen ID')%>%distinct()
  patients=c(unlist(patients),'16-01254')

  drugs<-getPatientDrugResponses(unlist(patients))

 patData<-readxl::read_xlsx(exp.3.megafile,sheet='Clinical Summary')%>%
   subset(labId%in%unlist(patients))%>%
   select(`AML sample`='labId',gender,ageAtDiagnosis, priorMalignancyType,
          vitalStatus,overallSurvival,causeOfDeath)%>%
          distinct()%>%
          left_join(drugs)%>%
   mutate(Value=tidyr::replace_na(Value,0))
#   subset(!is.na(Value))

 if(!is.null(synid)){##we need to delete rows and re-uplooad
   synTableUpdate(patData,synid)
 }else{
  synTableStore(patData,'BeatAML Pilot Drug and Clinical Data')
 }

 return(patData)

}

#' store drug class for now..
#' @import dplyr
storeDrugClassInfo<-function(){
  syn=synapseLogin()
  beat.samps<-syn$get('syn22130788')$path
  drug.class<-readxl::read_xlsx(beat.samps,sheet='Table S11-Drug Families')%>%
    distinct()
  synTableStore(drug.class,'Drug Classes')
}

# get all molecular data for patient baselines
#'@export
getPatientMolecularData<-function(synid='syn22172602'){
  syn=synapseLogin()
    exp.3.megafile=syn$get('syn22130786')$path

  patients<-readxl::read_xlsx(exp.3.megafile,sheet='Sample Summary')%>%
    dplyr::select('Specimen ID')%>%distinct()
  patients=c(unlist(patients),'16-01254')
  rna<-getPatientTranscript(unlist(patients))
  variants<-getPatientVariants(unlist(patients))
  moreVars<-addWave34GeneticsToTable()
  variants<-rbind(variants,moreVars)
  prots<-getPatientBaselines()
  patientMolecularData<-prots%>%
    full_join(rna,by=c('AML sample','Gene'),na_matches="never")%>%
    full_join(variants,by=c('AML sample','Gene'),na_matches="never")%>%
    mutate(transcriptCounts=tidyr::replace_na(transcriptCounts,0))%>%
    mutate(`Tumor VAF`=tidyr::replace_na(`Tumor VAF`,0))%>%
    mutate(LogFoldChange=tidyr::replace_na(LogFoldChange,0))

  if(!is.null(synid)){##we need to delete rows and re-uplooad
    synTableUpdate(patientMolecularData,synid)
    }else{
    synTableStore(patientMolecularData,'BeatAML Pilot Molecular Data')
  }

  return(patientMolecularData)
}

#' get proteomic data for beataml samples
#' @import dplyr
getPatientBaselines<-function(){
  library(dplyr)
  syn=synapseLogin()
  dat<-read.table(syn$get('syn22130778')$path,sep='\t',header=T)#'../../Projects/CPTAC/exp_3'/PTRC_baseline_global_std_ref_with_genes.txt,sep='\t',header=T)
  patientProtSamples<-dat%>%tidyr::pivot_longer(cols=c(5:ncol(dat)),names_to='Sample', 
                                                values_to='LogFoldChange')%>%
    dplyr::mutate(specId=stringr::str_replace(Sample,"X",""))%>%
    dplyr::mutate(Sample=stringr::str_replace(specId,stringr::fixed("."),"-"))%>%
    dplyr::select(`AML sample`='Sample',Gene, LogFoldChange)%>%
    mutate(LogFoldChange=as.numeric(LogFoldChange))%>%
    mutate(LogFoldChange=tidyr::replace_na(LogFoldChange,0))

  #this data is stored elsewhere
  metadata<-readAndTidySensMetadata()%>%
    dplyr::select(`AML sample`,Barcode)%>%
    distinct()%>%
    mutate(Barcode=as.character(Barcode))

  dat2<-read.table(syn$get('syn22130842')$path,sep='\t',header=T)
  secondData<-dat2%>%
    tidyr::pivot_longer(cols=c(3:ncol(dat2)),names_to='Sample', 
                        values_to='LogFoldChange')%>%
    dplyr::mutate(specId=stringr::str_replace(Sample,"X",""))%>%
    dplyr::mutate(Sample=stringr::str_replace(specId,stringr::fixed("."),"-"))%>%
    dplyr::select(Barcode='Sample',Gene, LogFoldChange)%>%
    mutate(LogFoldChange=as.numeric(LogFoldChange))%>%
    mutate(LogFoldChange=tidyr::replace_na(LogFoldChange,0))%>%
    left_join(metadata)%>%
    select(-Barcode)

  return(rbind(patientProtSamples,secondData))
}


#' @export
#' @import dplyr
getPatientPhosphoBaselines<-function(){
    library(dplyr)
    syn<-synapseLogin()                          
    dat<-read.table(syn$get('syn22130779')$path,sep='\t',header=T)

    patientPhosphoSamples<-dat%>%
        tidyr::pivot_longer(-c(Entry,Gene,site,Peptide,ids,Entry_name),"Sample",
                            values_to='LogFoldChange')%>%
        dplyr::mutate(specId=stringr::str_replace(Sample,"X",""))%>%
        dplyr::mutate(Sample=stringr::str_replace(specId,stringr::fixed("."),"-"))%>%
        dplyr::select(Sample,Gene, site,Peptide,LogFoldChange)%>%
      mutate(LogFoldChange=as.numeric(LogFoldChange))%>%
      subset(!is.na(LogFoldChange))
    
    #  mutate(LogFoldChange=tidyr::replace_na(LogFoldChange,0))

    synTableStore(patientPhosphoSamples,'BeatAML Pilot Phosphoproteomics')

    return(patientPhosphoSamples)
}

  ##Added 2/1/2021 --> decided we needed uncorrected phospho
getUnNormPhosphoBaselines<-function(){
  library(dplyr)
  syn<-synapseLogin()
     dat<-read.table(syn$get('syn24610481')$path,sep='\t',header=T)

    patientPhosphoSamples<-dat%>%
        tidyr::pivot_longer(-c(Entry,Gene,site,Peptide,ids,Entry_name),"Sample",
                            values_to='LogFoldChange')%>%
        dplyr::mutate(specId=stringr::str_replace(Sample,"X",""))%>%
        dplyr::mutate(Sample=stringr::str_replace(specId,stringr::fixed("."),"-"))%>%
        dplyr::select(Sample,Gene, site,Peptide,LogFoldChange)%>%
      mutate(LogFoldChange=as.numeric(LogFoldChange))%>%
      tidyr::drop_na(LogFoldChange)
      #mutate(LogFoldChange=tidyr::replace_na(LogFoldChange,0))

    synTableStore(patientPhosphoSamples,'BeatAML Pilot Phosphoproteomics Unnormalized')

    return(patientPhosphoSamples)

}


getNewSorafenibPhospho<-function(){
  ##Added 2/1/2021 --> decided we needed uncorrected phospho
  library(dplyr)
  library(tibble)
  syn<-synapseLogin()
  metadata<-readxl::read_xlsx(syn$get('syn22314059')$path)%>%
    dplyr::select(specimen='Specimen ID',`FLT ITD`, `Sorafenib sensitivity`,`Sorafenib IC50`, `Cell number`,`Treatment`)%>%
    rowwise()%>%mutate(`AML sample`=stringr::str_split_fixed(specimen,stringr::fixed('.'),2)[1])%>%
    subset(!is.na(Treatment))
   
  phdat<-read.csv2(syn$get('syn22313433')$path,sep='\t')%>%
    dplyr::select(-Entry,ids)%>%
    tidyr::pivot_longer(-c(Gene,site,Peptide),values_to='LogFoldChange',names_to='specIds')%>%
    mutate(LogFoldChange=as.numeric(as.character(LogFoldChange)))%>%
    tidyr::drop_na(LogFoldChange)%>%
   # mutate(LogFoldChange=tidyr::replace_na(LogFoldChange,0))%>%
    mutate(specimen=stringr::str_replace(specIds,'X',''))%>%
    mutate(specimen=stringr::str_replace(specimen,'\\.','-'))%>%
    left_join(metadata)%>%select(-c(specIds,specimen))%>%distinct()%>%subset(!is.na(`AML sample`))
  
    synTableStore(phdat,'Sorafenib treated Phosphoproteomics Unnormalized')

}

getNewComboPhospho<-function(){
  #Added 2/2/2020 --> uncorrected phospho
  library(dplyr)
   metadata<-amlresistancenetworks::readAndTidySensMetadata()
  syn<-synapseLogin()
  dat<-read.table(syn$get('syn24240354')$path,sep='\t',header=T)
  gilt.sens.pdata<-dat%>%
    tidyr::pivot_longer(-c(Entry,Gene,site,Peptide,ids),"Sample")%>%
    dplyr::mutate(Barcode=as.numeric(stringr::str_replace(Sample,"X","")))%>%
    dplyr::left_join(metadata,by='Barcode')%>%
  tidyr::drop_na(value)
  
   # mutate(value=tidyr::replace_na(value,0))

  synTableStore(gilt.sens.pdata,'Drug Combination Phosphoproteomic Unnormalized')
  return(gilt.sens.pdata)

}

#' getSorafenibSamples
#' Gets sorafenib trreated samples and stores to synaspe tables
#' @import dplyr
#' @import readxl
#' @import stringr
#' @export
getSorafenibSamples<-function(){
  library(dplyr)
  syn<-synapseLogin()
  metadata<-readxl::read_xlsx(syn$get('syn22314059')$path)%>%
    dplyr::select(specimen='Specimen ID',`FLT ITD`, `Sorafenib sensitivity`,`Sorafenib IC50`, `Cell number`,`Treatment`)%>%
    rowwise()%>%mutate(`AML sample`=stringr::str_split_fixed(specimen,stringr::fixed('.'),2)[1])%>%
    subset(!is.na(Treatment))
  
  pdat<-read.csv2(syn$get('syn24227680')$path,sep='\t')%>%
    dplyr::select(-ids)%>%
    tidyr::pivot_longer(-Gene,values_to='LogFoldChange',names_to='specIds')%>%
    mutate(LogFoldChange=as.numeric(as.character(LogFoldChange)))%>%
    mutate(LogFoldChange=tidyr::replace_na(LogFoldChange,0))%>%
    mutate(specimen=unlist(stringr::str_replace(specIds,'X','')))%>%
    mutate(specimen=stringr::str_replace(specimen,'\\.','-'))%>%
    left_join(metadata)%>%
    select(-c(specIds,specimen))%>%distinct()%>%
    subset(!is.na(`AML sample`))
  
  phdat<-read.csv2(syn$get('syn22313433')$path,sep='\t')%>%
    dplyr::select(-Entry,ids)%>%
    tidyr::pivot_longer(-c(Gene,site,Peptide),values_to='LogFoldChange',names_to='specIds')%>%
    mutate(LogFoldChange=as.numeric(as.character(LogFoldChange)))%>%
    mutate(LogFoldChange=tidyr::replace_na(LogFoldChange,0))%>%
    mutate(specimen=stringr::str_replace(specIds,'X',''))%>%
    mutate(specimen=stringr::str_replace(specimen,'\\.','-'))%>%
    left_join(metadata)%>%select(-c(specIds,specimen))%>%distinct()%>%subset(!is.na(`AML sample`))
  
  synTableStore(pdat,'Sorafenib treated proteomics')
  synTableStore(phdat,'Sorafenib treated Phospho-proteomics')
}
