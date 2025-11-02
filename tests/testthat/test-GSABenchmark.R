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

