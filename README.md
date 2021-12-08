# beatAMLproteomics

Analysis of BEAT AML patient Pilot proteomics data.
 <!-- badges: start -->
  [![Render-RMDs](https://github.com/PNNL-CompBio/beatAMLpilotProteomics/workflows/render-RMDs/badge.svg)](https://github.com/PNNL-CompBio/beatAMLpilotProteomics/actions)

  <!-- badges: end -->
  

## Repository Overview

This repository contains scripts that pull data from a [Synapse repository](http://synapse.org/ptrc) to carry out the various analysis steps needed. You will need to acquire a [synapse username](http://synapse.org/register) to access synapse, and become a [certified user](https://docs.synapse.org/articles/accounts_certified_users_and_profile_validation.html) to add data, but after that you will be set with future projects. You will then need to navigate to the [PNNL/OHSU Synapse page](http://synapse.org/ptrc) to request access.


### Before you begin...

This repository is only the basic structure of the tools needed, not the end-to-end analysis. Here are the steps you'll need to use this:

1- Read up on the tools
  - GitHub requires some basic protocols such as pull requests and commits, you should try to get a basic understanding. I found this [tutorial](https://medium.com/@jonathanmines/the-ultimate-github-collaboration-guide-df816e98fb67) that can serve as a starting point.
  - Synapse also has a bit of learning curve. To understand what Synapse is1 and isn't, check out [this document](https://docs.synapse.org/articles/getting_started.html).
2- Get [RStudio](http://rstudio.org). Basic R is essential, but RStudio will make your life a lot easier, I promise!
3- Install the [synapse python client](https://python-docs.synapse.org/build/html/index.html), create a [`.synapseConfig` file](https://python-docs.synapse.org/build/html/Credentials.html) in your home directory.
4- Clone this repository - it has all that you will need to contribute and run this analysis.

## Beat AML Processing
Here we descirbe the processing of the BeatAML Data

### Proteomics and phosphoproteomics processing

This repository contains the code for the normalization and processing.

This is derived from the P3 proteomics workflow and can be found in the [Ex10_proteomics](./Ex10_proteomics) folder. It contains the scripts required to process and normalize the data. It also contains the [study design](./proteomics/study_design) files that are required to do the processing.

Once the data is processed from DMS it is uploaded to Synapse in the [Pilot raw data folder](https://www.synapse.org/#!Synapse:syn22130775) folder.

The data was pushed from raw files to long-form tables for facile querying and viewing:

| Description | Normalization/filtering| Link |
| --- | --- | --- |
| Global Proteomics, transcriptomics and genetic data |  |[syn22172602](https://www.synapse.org/#!Synapse:syn22172602/tables/)|
| Global phosphoproteomics | Normalized |[syn22156830](https://www.synapse.org/#!Synapse:syn22156830/tables/)|
| Global phosphoproteomics | Unnormalized | [syn24227903](https://www.synapse.org/#!Synapse:syn24227903/tables/)|
| Drug response and clinical data | all of the data combined | [syn22170540](https://www.synapse.org/#!Synapse:syn22170540/tables/)


## Data Analysis
This data analysis is mostly done, but can be described below.

### Data summary

We first summarize the data we have collected in Figure 1. These panels are produced in [the Figure 1 markdown](./pilotManuscript/fig1_dataSummary.Rmd). 

### Model building

The model building takes place outside of an R markdowns, in the [buildAndStoreModels.R](./buildAndStoreModels.R). This builds the models and stores them in a [synapse folder](https://www.synapse.org/#!Synapse:syn26529370) for later use. 

### Model assessment
Using the k-fold cross validation we re-assess the error rate in the models. 

### Model selection

### Model validation