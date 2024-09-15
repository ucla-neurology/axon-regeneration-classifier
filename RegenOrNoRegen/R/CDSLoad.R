# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
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



CDSLoad <-function (CDS,GeneIDType="SYMBOL"){
  CDS <- estimateSizeFactors(CDS)
  CDS <- classify_cells(CDS, Hugo_classifier,
                        db = org.Mm.eg.db,
                        cluster_extend = TRUE,
                        cds_gene_id_type = GeneIDType)
  UMAP=CDS@reducedDimS
  Regeneration_Index=CDS@phenoData@data$cluster_ext_type
  Regeneration_Index[Regeneration_Index=="NonRegeneratingCST(Cl1)"]<- "NonRegenerating"
  Regeneration_Index[Regeneration_Index== "RegeneratingCST(Cl1)"]<-"Regenerating"
  Regeneration_Index[Regeneration_Index== "Cl1"]<- "Unknown"
  Regeneration_Index=factor(Regeneration_Index, levels=c("Unknown", "NonRegenerating","Regenerating"))
  qplot(UMAP[,1], UMAP[,2], color = Regeneration_Index,xlab="UMAP1",ylab="UMAP2") + theme_bw()+scale_color_manual(values=c("#F8766D","#00BA38","#619CFF"))

}


