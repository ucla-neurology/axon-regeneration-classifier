library(Seurat)
library(readxl)
library(garnett)

set.seed(123)

bulk_data <- read_excel("Mouse In Vitro RNA seq_Complete Gene List extract HL004.xlsx", sheet = 1)
bulk_data <- as.data.frame(bulk_data)

# Set the desired parameters
# The mean and the sd are base on the gene_means from the bulk rna seq data.
desired_mean <- 1.138563
desired_sd <- 0.765353  # Variability corresponds to the standard deviation

# Generate 19405 random gene means with the specified mean and standard deviation
set.seed(123)  # Set seed for reproducibility
random_gene_means <- rnorm(19405, mean = desired_mean, sd = desired_sd)

n_genes <- length(random_gene_means)
n_cells <- 1000
standard_deviation <- 1.5

# Initialize an empty matrix to store the simulated counts
random_simulated_counts <- matrix(0, nrow = n_genes, ncol = n_cells)

# Loop over each gene to simulate its counts across 1000 cells
for (i in 1:n_genes) {
  # Generate random normal data for each gene with mean = gene_means[i] and sd = standard_deviation
  random_simulated_counts[i, ] <- rnorm(n_cells, mean = random_gene_means[i], sd = standard_deviation)

  # Ensure no negative counts (since gene expression counts cannot be negative)
  random_simulated_counts[i, random_simulated_counts[i, ] < 0] <- 0
}

#name the rows and columns to match gene IDs and cell numbers
rownames(random_simulated_counts) <- bulk_data$GeneID
colnames(random_simulated_counts) <- paste0("Cell", 1:n_cells)

row_means <- rowMeans(random_simulated_counts)
head(row_means)

# Check the dimensions
dim(random_simulated_counts)

Randomly_simulated_Seurat_Object <- CreateSeuratObject(counts = random_simulated_counts)
Randomly_simulated_Seurat_Object <- NormalizeData(Randomly_simulated_Seurat_Object)
Randomly_simulated_Seurat_Object <- FindVariableFeatures(Randomly_simulated_Seurat_Object)
Randomly_simulated_Seurat_Object <- ScaleData(Randomly_simulated_Seurat_Object)
Randomly_simulated_Seurat_Object <- RunPCA(Randomly_simulated_Seurat_Object)
Randomly_simulated_Seurat_Object <- RunUMAP(Randomly_simulated_Seurat_Object, dims = 1:10)

# Apply the classifier to your Seurat object
source("RegenOrNoRegen/R/SeuratLoad.R")
load("RegenOrNoRegen/data/Hugo_classifier.rda")
random_ResultSeurat <- SeuratLoad(Randomly_simulated_Seurat_Object, GeneIDType = "ENSEMBL")

# View the classifications
print(random_ResultSeurat)