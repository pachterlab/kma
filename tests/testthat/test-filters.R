context('Filters')

test_that('aggregate',
    {
        tmp_ir <- list()
        tmp_ir$flat <- data.frame(
            intron = c("i1", "i1", "i1", "i1", "i2", "i2", "i2", "i2"),
            sample = c("s1", "s2", "s3", "s4", "s1", "s2", "s3", "s4"),
            condition = c(rep("c1", 3), "c2", rep("c1", 3), "c2"),
            f_1 = c(TRUE, TRUE, FALSE, TRUE,
                TRUE, TRUE, TRUE, FALSE),
            f_2 = c(TRUE, TRUE, TRUE, TRUE,
                TRUE, TRUE, TRUE, TRUE),
            retention = rnorm(8))
        expect_equal(aggregate_filters(tmp_ir)$f_all, c(FALSE, TRUE, TRUE, FALSE))
    })

