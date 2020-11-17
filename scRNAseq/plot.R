options(stringsAsFactors = F)

library(devtools)
install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("Seurat")
safelyLoadAPackageInCRANorBioconductor("ggplot2")

data <- list()
seurat.obj <- list()

stages <- c("E11", "E13", "E15")
figure.table <- read.delim("figures_stage_genes.txt")
# First steps
for(stage in stages){
  # Read the 10X files
  data[[stage]] <- Read10X(data.dir = paste0(stage, "_HL/"))
  # Check the number of cells and genes
  print(dim(data[[stage]]))
  # Create a Seurat object using the same filterings as in Kelly et al.
  seurat.obj[[stage]] <- CreateSeuratObject(counts = data[[stage]], project = stage,
                                            min.cells = 0, min.features = 200)
  # Compute the percentage of reads in mt genes
  seurat.obj[[stage]][["percent.mt"]] <- PercentageFeatureSet(seurat.obj[[stage]], pattern = "^mt-")
  # Filter as in Kelly et al.
  seurat.obj[[stage]] <- subset(seurat.obj[[stage]], subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
  # Check the number of cells after filtering
  print(dim(seurat.obj[[stage]]))
  # Normalize the expressions
  seurat.obj[[stage]] <- NormalizeData(seurat.obj[[stage]])
  # Find the variable features to do PCA before scaleData
  seurat.obj[[stage]] <- FindVariableFeatures(seurat.obj[[stage]], selection.method = "vst", nfeatures = 2000)
  # Scale data regressing the percent.mt
  seurat.obj[[stage]] <- ScaleData(seurat.obj[[stage]], vars.to.regress = "percent.mt")
}
# Run PCA/UMAP
for(stage in stages){
  # Run PCA using variable features
  seurat.obj[[stage]] <- RunPCA(seurat.obj[[stage]], features = VariableFeatures(object = seurat.obj[[stage]]))
  # Run UMAP using 10 PCA
  seurat.obj[[stage]] <- RunUMAP(seurat.obj[[stage]], dims = 1:10)
}
# Plot
for(i in 1:nrow(figure.table)){
  stage <- figure.table$stage[i]
  figure <- figure.table$figure[i]
  genes <- strsplit(figure.table$genes[i], ",")[[1]]
  g <- FeaturePlot(seurat.obj[[stage]], features = genes,
                   order = TRUE, reduction = "umap", ncol = 3)
  nrow.in.g <- length(genes) %/% 3
  svg(paste0("Fig", figure, ".svg"), width = 14, height = 14 / 3 * nrow.in.g)
  print(g)
  dev.off()
}
