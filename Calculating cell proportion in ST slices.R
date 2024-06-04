library(spacexr)
library(ggpubr)
library(CARD)

spatial_count <- Read10X(exp_mat[i])
ST <- CreateSeuratObject(spatial_count, assay = "Spatial")
image <- Read10X_Image(spatial[i],filter.matrix = F)
image <- image[Cells(x = ST)]
  
spatial_location <- image@coordinates
spatial_location <- spatial_location[,3:2]
colnames(spatial_location) <- c("x","y")
spatial_location$y <- -spatial_location$y
sc_count <- spatial_reference@assays$RNA@counts
sc_meta <- spatial_reference@meta.data[,c(1,26,4)]
colnames(sc_meta) <- c("cellID","cellType","sampleInfo")
sc_meta$cellID <- rownames(sc_meta)
  
CARD_obj = createCARDObject(
    sc_count = sc_count,
    sc_meta = sc_meta,
    spatial_count = spatial_count,
    spatial_location = spatial_location,
    ct.varname = "cellType",
    ct.select = unique(sc_meta$cellType),
    sample.varname = "sampleInfo",
    minCountGene = 100,
    minCountSpot = 5) 
  
CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)

puck <- SpatialRNA(spatial_location, spatial_count)

reference = Reference(spatial_reference@assays$RNA@counts, as.factor(spatial_reference$celltype3))

myRCTD <- create.RCTD(puck, reference, max_cores = 20, test_mode = FALSE,CELL_MIN_INSTANCE = 20) # here puck is the SpatialRNA object, and reference is the Reference object.
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')

metadata <- myRCTD@results$weights %>% normalize_weights() %>% data.frame()
norm_weights <- normalize_weights(weights)
metadata2 <- CARD_obj@Proportion_CARD %>% data.frame()
all(rownames(metadata)==rownames(metadata2))
all(colnames(metadata2)==colnames(metadata))

metadata <- (metadata[,colnames(metadata2)]+metadata2)/2
colnames(metadata) <- gsub(".","_",colnames(metadata),fixed = T)
ST <- AddMetaData(ST,metadata)
