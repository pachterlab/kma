context('Intron Retention')

test_that('construction',
    {

        i2t <- data.frame(
            intron = c('i1', 'i1', 'i2', 'i3', 'i3'),
            target_id = c('t1', 't2', 't1', 't4', 't1'),
            stringsAsFactors = T)

        targExp <- data.frame(
            samp1 = c(1, 0.5, 4, 1, 0, 10, 2),
            samp2 = c(1, 0.25, 0, 2, 8, 3, 3),
            target_id = c('i1', 'i2', 't1', 't2', 't3', 't4', 'i3'),
            stringsAsFactors = T)

        # targExp <- data.table(targExp)
        # i2t <- data.table(i2t)
        # setkey(targExp, target_id)
        # setkey(i2t, target_id)

        # debug(newIntronRetention)
        ir <- newIntronRetention(targExp, i2t, factor(c('c1', 'c2')))


    })

# i2t <- data.frame(
#     intron = c('i1', 'i1', 'i2', 'i3', 'i3'),
#     target_id = c('t1', 't2', 't1', 't4', 't1'),
#     stringsAsFactors = T)
# 
# targExp <- data.frame(
#     samp1 = c(1, 0.5, 4, 1, 0, 10, 2),
#     samp2 = c(1, 0.25, 0, 2, 8, 3, 3),
#     target_id = c('i1', 'i2', 't1', 't2', 't3', 't4', 'i3'),
#     stringsAsFactors = T)
# 
# targExp[i2t][,lapply(.SD,sum),by = iCol]
# 
# 
# testFunction <- function(df1, df2)
# {
#     df1 <- data.table(df1, key = 'target_id')
#     df2 <- data.table(df2, key = 'target_id')
#     df1[df2][,lapply(.SD, sum), by = iCol]
# }
# 
# testFunction(as.data.frame(targExp), as.data.frame(i2t))
# 
# hi <- newIntronRetention(as.data.frame(targExp), as.data.frame(i2t),
#     factor(c('s1', 's2')))
