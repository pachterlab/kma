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

        ir <- newIntronRetention(targExp, i2t, factor(c('c1', 'c2')))

    })
