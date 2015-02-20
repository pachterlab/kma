# context('Read/write')
#
# test_that('construction',
#     {
#
#         i2t <- data.frame(
#             intron = c('i1', 'i2', 'i3', 'i3', 'i1'),
#             target_id = c('t1', 't1', 't4', 't1', 't2'),
#             gene = c('g1', 'g2', 'g1', 'g1', 'g1'),
#             intron_extension = c('i1_ext', 'i2_ext', 'i3_ext', 'i3_ext',
#                 'i1_ext'),
#             stringsAsFactors = F)
#
#         targExp <- data.frame(
#             samp1 = c(1, 0.5, 4, 1, 0, 10, 2),
#             samp2 = c(1, 0.25, 0, 2, 8, 3, 3),
#             target_id = c('i1_ext', 'i2_ext', 't1', 't2', 't3', 't4', 'i3_ext'),
#             stringsAsFactors = F)
#
#         ir <- newIntronRetention(targExp, i2t, factor(c('c1', 'c2')), FALSE)
#
#         expect_equal(rownames(ir@numerator), rownames(ir@denominator))
#
#         expect_equal(ir@numerator$samp1, c(1, 0.5, 2), tolerance = 0.01)
#         expect_equal(ir@numerator$samp2, c(1, 0.25, 3), tolerance = 0.01)
#         expect_equal(ir@retention$samp1, c(1 / (4 + 1), 0.5 / 4, 2 / (10 + 4)),
#             tolerance = 0.01)
#         expect_equal(ir@retention$samp2, c(1 / (0 + 2), 0.25 / 0, 3 / (3 + 0)),
#             tolerance = 0.01)
#
#         ir <- newIntronRetention(targExp, i2t, factor(c('c1', 'c2')), TRUE)
#         expect_equal(ir@numerator$samp1, c(1, 0.5, 2), tolerance = 0.01)
#         expect_equal(ir@numerator$samp2, c(1, 0.25, 3), tolerance = 0.01)
#         expect_equal(ir@retention$samp1, c(1 / (4 + 1 + 1),
#                 0.5 / (4 + 0.5), 2 / (10 + 4 + 2)), tolerance = 0.0001)
#         expect_equal(ir@retention$samp2, c(1 / (0 + 2 + 1),
#                 0.25 / (0 + 0.25), 3 / (3 + 0 + 3)), tolerance = 0.0001)
#
#     })
#
