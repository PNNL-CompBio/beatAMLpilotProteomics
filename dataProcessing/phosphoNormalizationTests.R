##Compare phospho data with normalization and without

source("beatAMLdata.R")
#get phospho data
library(dplyr)


##lets first check in AML patient data
loadBeatAMLMolecularData()
mean.orig<-orig.pat.phos%>%select(Gene,Sample,site,siteCorrected='LogFoldChange')%>%distinct()

mean.unnorm<-pat.phos%>%select(Gene,Sample,site,unCorrected='LogFoldChange')%>%distinct()

combined<-mean.orig%>%
  left_join(mean.unnorm)%>%
  left_join(select(pat.data,c('Gene',Sample='AML sample','proteinLevels'))%>%
              distinct)


merged<-combined%>%tidyr::pivot_longer(cols=c('siteCorrected','unCorrected'),
                                       names_to='Normalization',values_to='PhosphoLevels')

corvals<-merged%>%group_by(Sample,Normalization,use='pairwise.complete.obs')%>%
  summarize(corval=cor(proteinLevels,PhosphoLevels,use='pairwise.complete.obs'))

ggplot(corvals)+geom_bar(aes(x=Sample,y=corval,fill=Normalization),position='dodge',stat='identity')+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p<-ggplot(merged)+geom_point(aes(x=proteinLevels,y=PhosphoLevels,col=Normalization))+facet_grid(~Normalization)



#phosph-to-protein cor