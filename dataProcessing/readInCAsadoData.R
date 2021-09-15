##read in excel spreadsheet from casado et al.


library(readxl)
library(tidyverse)
library(reticulate)

##first read in proteomics data
tab <- readxl::read_xlsx('Copy of Data file S5 Proteins With Values.xlsx')%>%
  tidyr::pivot_longer(cols=c(6:35,37:66),names_to = 'sample',values_to='logratio')%>%
  mutate(sample=stringr::str_remove(sample,"\\.\\.\\.[0-9]*"))%>%
  dplyr::select(-c(average,ttest,fold,FDR,Significance,"Z-Score"))%>%
  mutate(symbol=stringr::str_remove(`...1`,'_HUMAN'))

##then read in chemical data
dtab<-readxl::read_xlsx('41375_2018_32_MOESM2_ESM.xlsx',skip=2,na='ND')%>%
  tidyr::pivot_longer(cols=c(10:14),names_to='compound',values_to = 'sensOverControl')%>%
  tidyr::pivot_longer(cols=c(10:14),names_to='compound2',values_to='ED50')%>%
  tidyr::pivot_longer(cols=c(26:42),names_to='Differentiation Marker',values_to='log2fc')

chemData<-dtab%>%
  dplyr::select(sample='Vial ID',compound,sensOverControl,'PAT ID',Gender,'Response to therapy',age='AGE (at diagnose)')%>%
  distinct()%>%
  mutate(compound=stringr::str_remove(compound,"\\.\\.\\.[0-9]*"))


syn <- reticulate::import('synapseclient')
sync <- syn$login()

tab1 <- syn$build_table(name='Casado 2008 Clinical Data',parent='syn22128879',values=chemData)
tab2 <-syn$build_table(name='Casado 2008 Protoemics',parent='syn22128879',values=tab)

sync$store(tab1)
sync$store(tab2)