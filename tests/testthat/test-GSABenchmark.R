test_that("runMethods works", {
    gsaMethods <- supportedMethods()
    scObj <- runGSAMethods(scObj, 'label', geneSets, gsaMethods)
    nCols <- length(colnames(scObj[[]]))
    nRuns <- length(geneSets) * length(gsaMethods)
    expect_equal(nCols, 6 + nRuns)
    runsDF <- scObj[[]][, 7:nCols]
    expect_equal(max(runsDF), 1)
    expect_equal(min(runsDF), 0)
})

test_that("runMethodShuffle works", {
    scObj <- runMethodShuffle(scObj, 'label', geneSets, 'CSOA', 0.2, 0.2)
    runsDF <- scObj[[]][, 41:42]
    expect_equal(max(runsDF), 1)
    expect_equal(min(runsDF), 0)
})

test_that("runBenchmark works", {
    smr <- runBenchmark(scObj, 'label', geneSets, supportedMethods(), FALSE)
    expect_equal(length(smr), 4)
    expect_equal(length(smr[[1]]), 8)
    expect_equal(length(smr[[2]]), 2)
    expect_equal(length(smr[[3]]), 9)
    nMethods <- length(supportedMethods())
    for (i in seq(3)){
        temp <- unlist(smr[[i]])
        expect_equal(nrow(smr[[i]][[1]]), nMethods)
        expect_equal(ncol(smr[[i]][[1]]), 3)
        expect_gte(min(temp), 0)
        expect_lte(max(temp), 1)
    }
})

test_that("efficiencyBenchmark works", {
    smr <- efficiencyBenchmark(scObj, 'label', geneSets, c('CSOA', 'Zscore'),
                               verbose=FALSE)
    expect_equal(length(smr), 2)
    for (i in seq(2)){
        expect_equal(nrow(smr[[i]]), 2)
        expect_equal(ncol(smr[[i]]), 3)
        expect_gte(max(smr[[i]]), 0)
    }
})

test_that("simple visualization functions work", {
    p <- scorePlot(smr[[1]][[1]])
    expect_equal(length(intersect(is(p), c('gg', 'ggplot2::ggplot'))), 1)

    p <- timePlot(smr[[5]])
    expect_equal(length(intersect(is(p), c('gg', 'ggplot2::ggplot'))), 1)

    p <- memoryPlot(smr[[5]])
    expect_equal(length(intersect(is(p), c('gg', 'ggplot2::ggplot'))), 1)

    p <- aggregateRankPlot(smr)
    expect_equal(length(intersect(is(p), c('gg', 'ggplot2::ggplot'))), 1)

    p <- ratioPlot(smr)
    expect_equal(length(intersect(is(p), c('gg', 'ggplot2::ggplot'))), 1)
})

test_that("complex visualization functions work", {
    titleSuffix <- 'Test'

    plots <- benchmarkPlots(smr[[1]], titleSuffix=titleSuffix)
    expect_equal(length(intersect(is(plots[[1]]),
                                  c('gg', 'ggplot2::ggplot'))), 1)

    plots <- allBenchmarkPlots(smr, titleSuffix=titleSuffix)
    expect_equal(length(intersect(is(plots[[1]]),
                                  c('gg', 'ggplot2::ggplot'))), 1)

    plots <- geneSetRankPlots(smr, titleSuffix=titleSuffix)
    expect_equal(length(intersect(is(plots[[1]]),
                                  c('gg', 'ggplot2::ggplot'))), 1)

    plots <- metricRankPlots(smr, titleSuffix=titleSuffix)
    expect_equal(length(intersect(is(plots[[1]]),
                                  c('gg', 'ggplot2::ggplot'))), 1)

    plots <- mdsPlots(scObj, smr, titleSuffix=titleSuffix)
    expect_equal(length(intersect(is(plots[[1]]),
                                  c('gg', 'ggplot2::ggplot'))), 1)

    plots <- corrPlots(scObj, smr, titleSuffix=titleSuffix)
    expect_equal(length(intersect(is(plots[[1]]),
                                  c('gg', 'ggplot2::ggplot'))), 1)

    plots <- predJaccardPlots(smr$predictions, titleSuffix=titleSuffix)
    expect_equal(length(intersect(is(plots[[1]]),
                                  c('gg', 'ggplot2::ggplot'))), 1)
})

test_that("rankAlignment score works", {
    v <- c(2, 3, 6, 7, 8, 4, 12, 9, 10)
    w <- c(3, 4, 5, 6, 2, 7, 8, 13, 3)
    expect_equal(rankAlignmentScore(v, w), 0.7936118, tolerance=0.0001)
})

test_that("nearestNeighbors works", {
    df <- data.frame(v = c(1, 2, 4, 5, 6),
                     w = c(2, 3, 1, 5, 8),
                     x = c(2, 8, 7, 1, 1),
                     y = c(2, 3, 2, 2, 4),
                     z = c(1, 9, 9, 7, 6))
    distMat <- as.matrix(stats::dist(df))
    rownames(distMat) <- c('v', 'w', 'x', 'y', 'z')
    colnames(distMat) <- c('v', 'w', 'x', 'y', 'z')
    res <- nearestNeighbors(distMat)
    expected <- setNames(c('y', 'x', 'w', 'z', 'y'), rownames(distMat))
    expect_equal(res, expected)
})

test_that("proximity works", {
    expect_equal(proximity(2, 3, 6), 0.8333333, tolerance=0.0001)
})

test_that("shuffleGenes works", {
    genes <- c('Gene_0226', 'Gene_0210', 'Gene_0280', 'Gene_0202',
               'Gene_0313', 'Gene_0101', 'Gene_0195')
    newGenes <- shuffleGenes(scObj, genes, 0.3, 0.9)
    expect_equal(length(intersect(genes, newGenes)), 5)
    expect_equal(length(newGenes), 50)
})

test_that("tabulateVector works", {
    v <- c(2, 3, 4, 19, 15, 25, 32, 8)
    res <- tabulateVector(v, paste0('r', seq(4)), paste0('c', seq(2)))
    df <- data.frame(c1 = c(2, 3, 4, 19),
                     c2 = c(15, 25, 32, 8),
                     row.names = paste0('r', seq(4)))
    expect_equal(res, df)
})

test_that("numCosine works", {
    res <- numCosine(c(2, 3, 6), c(4, 3, 2))
    expect_equal(res, 0.7693093, tolerance=0.0001)
})
