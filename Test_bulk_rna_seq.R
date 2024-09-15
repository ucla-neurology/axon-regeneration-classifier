library(Seurat)
library(readxl)
library(garnett)

set.seed(123)

bulk_data <- read_excel("Mouse In Vitro RNA seq_Complete Gene List extract HL004.xlsx", sheet = 1)
bulk_data <- as.data.frame(bulk_data)

gene_means <- exp(bulk_data$logFC_HL004_vs_Vehicle)
mean(gene_means)
sd(gene_means)

n_genes <- length(gene_means)
n_cells <- 1000
standard_deviation <- 1.5

# Initialize an empty matrix to store the simulated counts
simulated_counts <- matrix(0, nrow = n_genes, ncol = n_cells)

# Loop over each gene to simulate its counts across 1000 cells
for (i in 1:n_genes) {
  # Generate random normal data for each gene with mean = gene_means[i] and sd = standard_deviation
  simulated_counts[i, ] <- rnorm(n_cells, mean = gene_means[i], sd = standard_deviation)

  # Ensure no negative counts (since gene expression counts cannot be negative)
  simulated_counts[i, simulated_counts[i, ] < 0] <- 0
}

#name the rows and columns to match gene IDs and cell numbers
rownames(simulated_counts) <- bulk_data$GeneID
colnames(simulated_counts) <- paste0("Cell", 1:n_cells)

row_means <- rowMeans(simulated_counts)
head(row_means)

# Check the dimensions
dim(simulated_counts)

Seurat_Object <- CreateSeuratObject(counts = simulated_counts)
Seurat_Object <- NormalizeData(Seurat_Object)
Seurat_Object <- FindVariableFeatures(Seurat_Object)
Seurat_Object <- ScaleData(Seurat_Object)
Seurat_Object <- RunPCA(Seurat_Object)
Seurat_Object <- RunUMAP(Seurat_Object, dims = 1:10)

# Apply the classifier to your Seurat object
source("RegenOrNoRegen/R/SeuratLoad.R")
load("RegenOrNoRegen/data/Hugo_classifier.rda")
ResultSeurat <- SeuratLoad(Seurat_Object, GeneIDType = "ENSEMBL")

# View the classifications
print(ResultSeurat)
