context('Intron Retention')

test_that('construction',
    {

        i2t <- data.frame(
            intron = c('i1', 'i1', 'i2', 'i3', 'i3'),
            target_id = c('t1', 't2', 't1', 't4', 't1'),
            stringsAsFactors = F)

        targExp <- data.frame(
            samp1 = c(1, 0.5, 4, 1, 0, 10, 2),
            samp2 = c(1, 0.25, 0, 2, 8, 3, 3),
            target_id = c('i1', 'i2', 't1', 't2', 't3', 't4', 'i3'),
            stringsAsFactors = F)

        # print(targExp)
        # print(i2t)

        ir <- newIntronRetention(targExp, i2t, factor(c('c1', 'c2')))

        expect_equal(ir@numerator$samp1, c(1, 0.5, 2), tolerance = 0.01)
        expect_equal(ir@numerator$samp2, c(1, 0.25, 3), tolerance = 0.01)
        expect_equal(ir@retention$samp1, c(1 / (4 + 1), 0.5 / 4, 2 / (10 + 4)),
            tolerance = 0.01)
        expect_equal(ir@retention$samp2, c(1 / (0 + 2), 0.25 / 0, 3 / (3 + 0)),
            tolerance = 0.01)

        # debugonce(newIntronRetention)
        ir <- newIntronRetention(targExp, i2t, factor(c('c1', 'c2')), TRUE)
        expect_equal(ir@retention$samp1, c(1 / (4 + 1 + 1),
                0.5 / (4 + 0.5), 2 / (10 + 4 + 2)), tolerance = 0.0001)
        expect_equal(ir@retention$samp2, c(1 / (0 + 2 + 1),
                0.25 / (0 + 0.25), 3 / (3 + 0 + 3)), tolerance = 0.0001)

    })

test_that('lowExpressionFilter',
    {
        i2t <- data.frame(
            intron = c('i1', 'i2', 'i3', 'i4'),
            target_id = c('t1', 't2', 't3', 't4'),
            stringsAsFactors = F)

        targExp <- data.frame(
            samp1 = c(rep(0, 4), 1:4),
            samp2 = c(rep(0, 4), c(1, 0, 2, 3)),
            samp3 = c(rep(0, 4), 1:4),
            target_id = c('i1', 'i2', 'i3', 'i4', 't1', 't2', 't3', 't4'),
            stringsAsFactors = F)

        ir <- newIntronRetention(targExp, i2t, factor(c('c1', 'c1', 'c2')))
        ir <- lowExpressionFilter(ir, 0.25)

        expect_equal(ir@validIntrons$c1, c(F, F, T, T))
        expect_equal(ir@validIntrons$c2, c(F, T, T, T))
    })

test_that('retentionTest',
    {
        i2t <- data.frame(
            intron = c('i1', 'i2', 'i3', 'i4'),
            target_id = c('t1', 't2', 't3', 't4'),
            stringsAsFactors = F)

        targExp <- data.frame(
            samp1 = c(rep(0, 4), 1:4),
            samp2 = c(rep(0, 4), c(1, 0, 2, 3)),
            samp3 = c(rep(0, 4), 1:4),
            target_id = c('i1', 'i2', 'i3', 'i4', 't1', 't2', 't3', 't4'),
            stringsAsFactors = F)

        ir <- newIntronRetention(targExp, i2t, factor(c('c1', 'c1', 'c2')))
        ir <- lowExpressionFilter(ir, 0.25)

    })

test_that('retentionTestSingle',
    {
        set.seed(42)
        data.frame(
            samp1 = runif(10),
            samp2 = runif(10),
            samp3 = runif(10))
        # TODO: write me
    })
