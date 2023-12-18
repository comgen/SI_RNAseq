# Stoichiometric imbalance of gene expression in bipolar disorder

This is the collection of code and analysis used for the paper 
"Bipolar patients display stoichiometric imbalance of gene expression in post-mortem brain samples" 
by Asbjørn Holmgren et al (2023)

## Description
The file DE_absolute_v08.Rmd consists of code to import and normalise the RNAseq data and do differential gene expression analysis. 
It produces the results for "absolute gene expression" analysis, i,e. figure 1. It also produces the correlation heatmaps in Figure 2.
The file SI_relative_v09.Rmd consists of code to do the expression modelling, calculating the SI score and predict diagnosis (cross validation), 
as well as control experiments (HC vs HC and model in BD). 


## Getting Started

### Dependencies

All code is run in R v4. The dependencies are listed in each Markdown document. 

### Installing

The initial data have restricted access through the PsychENCODE Consortium. Access to the data is managed by the NIMH Repository and Genomics Resource, 
and the data are distributed via Synapse under the CommonMind HBCC (syn10623034), CommonMind CMC-Pitt (syn8241760), 
BrainGVEX (syn3270015) and BipSeq (syn5844980) studies.


## Help


```

```

## Authors

Contributors names and contact info

Asbjørn Holmgren (MSc)
asbjorn.holmgren@medisin.uio.no
Timothy Hughes (PhD)
timothy.hughes@medisin.uio.no

## Version History



## License


## Acknowledgments
