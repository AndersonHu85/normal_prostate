library(dplyr)
library(ggplot2)
library(cowplot)
library(harmony)
library(Seurat)
library(ggplot2)
library(cowplot)
library(ggsci)
library(RColorBrewer)
options(stringsAsFactors = F)

scRNA <- readRDS("./scRNA_total.rds")

VlnPlot(scRNA, c("nUMI","nGene","percent.mito"), group.by = "sample",pt.size = 0)

scRNA <- SetIdent(scRNA, cells = NULL, 
                  value = scRNA@meta.data$celltype2)

scRNA <- NormalizeData(scRNA)

scRNA <- FindVariableFeatures(scRNA, nfeatures = 2000)
VariableFeaturePlot(scRNA)
scRNA <- ScaleData(scRNA, verbose = T)
scRNA <- RunPCA(scRNA, npcs = 50, verbose = FALSE)
PCAPlot(scRNA,cols=celltype)
ElbowPlot(scRNA, ndims = 50)
scRNA <- RunHarmony(scRNA, c("sample"),max.iter.harmony = 20)
scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims = 1:50, do.plot = T)
scRNA <- FindClusters(scRNA, resolution = 1)

scRNA  <- RunUMAP(scRNA, reduction = "harmony", dims = 1:50, n.components = 2)

scRNA <- CellCycleScoring(scRNA, g2m.features = cc.genes$g2m.genes,
                        s.features = cc.genes$s.genes)
