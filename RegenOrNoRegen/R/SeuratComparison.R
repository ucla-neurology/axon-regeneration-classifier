
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

library(devtools)
library(roxygen2)
library(testthat)
library(garnett)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ggplot2)
library(Seurat)
library(scuttle)



SeuratComparison <-function (SEURAT,MetaData){
  RegenIndex=as.character(SEURAT$Regeneration_Index)
  SelectedMetaData=SEURAT[[MetaData]]
  Combined=cbind(RegenIndex,SelectedMetaData)
  d=table(Combined[[c("Level.1", "Regeneration_Index")]])
  ggplot(d, aes(RegenIndex, SelectedMetaData, fill= Number)) +
    geom_tile()
  return(SEURAT)
  }


