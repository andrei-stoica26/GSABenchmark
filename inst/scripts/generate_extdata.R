# A self-contained script to generate the data in extdata
#
createExtData <- function(){
    if (requireNamespace(c('GSABenchmark', 'hammers', 'qs2', 'scRNAseq',
                           'scuttle', 'Seurat'))){

        scObj <- scRNAseq::BaronPancreasData('human')
        scObj <- scuttle::logNormCounts(scObj)
        scObj <- Seurat::as.Seurat(scObj)
        scObj <- subset(scObj, label %in% c('acinar', 'ductal'))
        scObj <- scObj[, seq(251, 350)]
        scObj <- hammers::removeRareGenes(scObj)
        scObj <- Seurat::FindVariableFeatures(scObj)
        scObj <- Seurat::ScaleData(scObj)
        scObj <- Seurat::RunPCA(scObj)
        scObj <- Seurat::RunUMAP(scObj, dims=seq(10))

        acinarMarkers <- c('PRSS1', 'KLK1', 'CTRC', 'PNLIP', 'AKR1C3', 'CTRB1',
                           'DUOXA2', 'ALDOB', 'REG3A', 'SERPINA3', 'PRSS3',
                           'REG1B','CFB',	'GDF15', 'MUC1','ANPEP', 'ANGPTL4',
                           'OLFM4','GSTA1', 'LGALS2', 'PDZK1IP1', 'RARRES2',
                           'CXCL17','UBD', 'GSTA2', 'LYZ', 'RBPJL', 'PTF1A',
                           'CELA3A', 'SPINK1', 'ZG16', 'CEL', 'CELA2A',
                           'CPB1', 'CELA1','PNLIPRP1', 'RNASE1', 'AMY2B',
                           'CPA2','CPA1', 'CELA3B', 'CTRB2', 'PLA2G1B',
                           'PRSS2', 'CLPS', 'REG1A', 'SYCN')

        ductalMarkers <-  c('CFTR', 'SERPINA5', 'SLPI', 'TFF1', 'CFB', 'LGALS4',
                            'CTSH',	'PERP', 'PDLIM3', 'WFDC2', 'SLC3A1', 'AQP1',
                            'ALDH1A3', 'VTCN1',	'KRT19', 'TFF2', 'KRT7',
                            'CLDN4', 'LAMB3', 'TACSTD2', 'CCL2', 'DCDC2',
                            'CXCL2', 'CLDN10', 'HNF1B', 'KRT20', 'MUC1',
                            'ONECUT1', 'AMBP', 'HHEX', 'ANXA4', 'SPP1', 'PDX1',
                            'SERPINA3', 'GDF15', 'AKR1C3', 'MMP7', 'DEFB1',
                            'SERPING1', 'TSPAN8', 'CLDN1', 'S100A10',
                            'PIGR')

        acinarMarkers <- intersect(acinarMarkers, rownames(scObj))
        ductalMarkers <- intersect(ductalMarkers, rownames(scObj))
        geneSets <- setNames(list(acinarMarkers, ductalMarkers),
                             c('acinar', 'ductal'))
        qs2::qs_save(geneSets, 'inst/extdata/geneSets.qs2')

        gsaMethods <- GSABenchmark::supportedMethods()
        scObj <- GSABenchmark::runGSAMethods(scObj, 'label', geneSets,
                                             gsaMethods)
        qs2::qs_save(scObj, 'inst/extdata/scObj.qs2')

        smr <- GSABenchmark::runBenchmark(scObj, 'label', geneSets,
                                          gsaMethods)
        qs2::qs_save(smr, 'inst/extdata/smr.qs2')
    }
}

createExtData()
