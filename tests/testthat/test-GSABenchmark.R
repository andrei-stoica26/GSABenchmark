test_that("runBenchmark works", {
    smr <- runBenchmark(scObj, 'label', geneSets, supportedMethods(), FALSE)
    expect_equal(length(smr), 3)
    expect_equal(length(smr[[1]]), 8)
    expect_equal(length(smr[[2]]), 2)
    expect_equal(length(smr[[3]]), 7)
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

test_that("visualization functions work", {
    p <- scorePlot(smr[[1]][[1]])
    expect_equal(length(intersect(is(p), c('gg', 'ggplot2::ggplot'))), 1)
    p <- timePlot(smr[[4]])
    expect_equal(length(intersect(is(p), c('gg', 'ggplot2::ggplot'))), 1)
    p <- memoryPlot(smr[[4]])
    expect_equal(length(intersect(is(p), c('gg', 'ggplot2::ggplot'))), 1)
    plots <- benchmarkPlots(smr[[1]])
    expect_equal(length(intersect(is(plots[[1]]),
                                  c('gg', 'ggplot2::ggplot'))), 1)
    plots <- allBenchmarkPlots(smr)
    expect_equal(length(intersect(is(plots[[1]]),
                                  c('gg', 'ggplot2::ggplot'))), 1)
    plots <- geneSetRankPlots(smr)
    expect_equal(length(intersect(is(plots[[1]]),
                                  c('gg', 'ggplot2::ggplot'))), 1)
    plots <- metricRankPlots(smr)
    expect_equal(length(intersect(is(plots[[1]]),
                                  c('gg', 'ggplot2::ggplot'))), 1)
    p <- aggregateRankPlot(smr)
    expect_equal(length(intersect(is(p), c('gg', 'ggplot2::ggplot'))), 1)
    p <- ratioPlot(smr)
    expect_equal(length(intersect(is(p), c('gg', 'ggplot2::ggplot'))), 1)
    plots <- mdsPlots(scObj, smr, colorScheme='orichalc')
    expect_equal(length(intersect(is(plots[[1]]),
                                  c('gg', 'ggplot2::ggplot'))), 1)
    plots <- corrPlots(scObj, smr)
    expect_equal(length(intersect(is(plots[[1]]),
                                  c('gg', 'ggplot2::ggplot'))), 1)
})

test_that("rankAlignment score works", {
    v <- c(2, 3, 6, 7, 8, 4, 12, 9, 10)
    w <- c(3, 4, 5, 6, 2, 7, 8, 13, 3)
    expect_equal(rankAlignmentScore(v, w), 0.7936118, tolerance=0.0001)
})
