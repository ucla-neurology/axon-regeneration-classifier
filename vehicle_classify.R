library(Seurat)
library(readxl)
library(garnett)
library(stringr)

source("RegenOrNoRegen/R/SeuratLoad.R")
source("RegenOrNoRegen/R/SeuratReturn.R")
source("SeuratLoadReturn.R")
load("RegenOrNoRegen/data/Hugo_classifier.rda")

set.seed(123)

bulk_data <- read_excel("BrkThAI_Mouse In Vitro RNA seq_FPKM Read Counts.xlsx", sheet = 3)
bulk_data <- as.data.frame(bulk_data)

vehicle_FPKM_1 <- log1p(bulk_data$`Vehicle Replicate 1`)
mean(vehicle_FPKM_1)
sd(vehicle_FPKM_1)

vehicle_FPKM_2 <- log1p(bulk_data$`Vehicle Replicate 2`)
mean(vehicle_FPKM_2)
sd(vehicle_FPKM_2)

gene_means <- apply(cbind(vehicle_FPKM_1, vehicle_FPKM_2), 1, mean)
gene_sds <- apply(cbind(vehicle_FPKM_1, vehicle_FPKM_2), 1, sd)

n_genes <- length(vehicle_FPKM_1)
n_cells <- 1000
standard_deviation <- 0.3

# Initialize an empty matrix to store the simulated counts
simulated_FPKM_vehicle <- matrix(0, nrow = n_genes, ncol = n_cells)

# Loop over each gene to simulate its counts across 1000 cells
for (i in 1:n_genes) {
  # Generate random normal data for each gene with mean = gene_means[i] and sd = standard_deviation
  simulated_FPKM_vehicle[i, ] <- rnorm(n_cells, mean = gene_means[i], sd = standard_deviation)

  # Ensure no negative counts (since gene expression counts cannot be negative)
  simulated_FPKM_vehicle[i, simulated_FPKM_vehicle[i, ] < 0] <- 0
}

#name the rows and columns to match gene IDs and cell numbers
rownames(simulated_FPKM_vehicle) <- bulk_data$GeneID
colnames(simulated_FPKM_vehicle) <- paste0("Cell", 1:n_cells)

row_means <- rowMeans(simulated_FPKM_vehicle)
head(row_means)

# Check the dimensions
dim(simulated_FPKM_vehicle)

Seurat_Object <- CreateSeuratObject(counts = simulated_FPKM_vehicle)

Seurat_Object <- NormalizeData(Seurat_Object)

Seurat_Object <- FindVariableFeatures(Seurat_Object)

Seurat_Object <- ScaleData(Seurat_Object)

Seurat_Object <- RunPCA(Seurat_Object)

Seurat_Object <- RunUMAP(Seurat_Object, dims = 1:10)

# Apply the classifier to the Seurat object
Seurat_Object <- SeuratLoadReturn(Seurat_Object, GeneIDType = "ENSEMBL")

classified_count_df = table(Seurat_Object$Regeneration_Index)
known_cells <- sum(classified_count_df)
na_cells <- 1000 - known_cells
classified_count_df <- c(classified_count_df, "NA" = na_cells)
print(classified_count_df)

proportion_df = as.data.frame(prop.table(classified_count_df))
proportion_df$Class <- rownames(proportion_df)
rownames(proportion_df) <- NULL
colnames(proportion_df) <- c("Prop", "Class")
print(proportion_df)

ggplot(proportion_df, aes(x = Class, y = Prop, fill = Class)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = scales::percent(Prop)), vjust = -0.5) +
  scale_y_continuous(labels = scales::percent) +  # Convert y-axis to percentage
  ggtitle("Regeneration Class Proportions") +
  theme_minimal() +
  ylab("Proportion") +
  xlab("Regeneration Class") +
  scale_fill_manual(values = c("Unknown" = "lightblue", "NonRegenerating" = "lightgreen", "Regenerating" = "lightcoral"))

missing_genes <- tryCatch({
  SeuratLoad(Seurat_Object, GeneIDType = "ENSEMBL")
}, warning = function(w) {
  # Extract the warning message
  warning_message <- conditionMessage(w)
  print(warning_message)

  missing_genes_list <- str_extract_all(warning_message, "ENSMUSG[0-9]+")[[1]]

  print("List of missing genes:")
  print(missing_genes_list)

  return(missing_genes_list)
})
