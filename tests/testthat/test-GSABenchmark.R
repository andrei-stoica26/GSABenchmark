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
    smr <- runBenchmark(scObj, 'label', geneSets, supportedMethods(),
                        runEFBenchmark=FALSE)
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
    smr <- efficiencyBenchmark(scObj, 'label', geneSets[1], 'CSOA',
                               verbose=FALSE)
    expect_equal(length(smr), 2)
    for (i in seq(2)){
        expect_equal(nrow(smr[[i]]), 1)
        expect_equal(ncol(smr[[i]]), 2)
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
