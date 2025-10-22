sceObj <- scRNAseq::BaronPancreasData('human')
sceObj <- scuttle::logNormCounts(sceObj)

seuratObj <- suppressWarnings(Seurat::as.Seurat(sceObj))
seuratObj <- Seurat::FindVariableFeatures(seuratObj)
seuratObj <- Seurat::ScaleData(seuratObj)
seuratObj <- Seurat::RunPCA(seuratObj)
seuratObj <- suppressWarnings(Seurat::RunUMAP(seuratObj, dims=1:15))

sceObj <- scater::runPCA(sceObj)
sceObj <- hammers::computeSilhouette(sceObj, 'label')
sceObj <- scater::runUMAP(sceObj)
