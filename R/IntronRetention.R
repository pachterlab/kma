#' Class IntronRetention
#'
#' Class \code{IntronRentention}.
#'
#' @name IntronRetention-class
#' @rdname IntronRetention-class
#' @aliases IntronRetention-class
#' @export
setClass("IntronRetention",
         slots = list(retention = "data.frame",
                      numerator = "data.frame",
                      denominator = "data.frame",
                      labels = "character",
                      groups = "factor",
                      features = "data.frame",
                      validIntrons = "data.frame"))

#' Compute intron retention
#'
#' Create an IntronRetention object.
#'
#' @param targetExpression a matrix of nothing but expression values and one
#' column containing 'target_id's
#' @param intronToUnion a table with one column labeled 'intron' and another
#' labeled 'target_id' that gives all the target_ids which map to the intron
#' @param groups a vector with the grouping
#' @param psi if TRUE, compute the psi value. otherwise, compute a rate
#' @return a IntronRetention object
#' @export
newIntronRetention <- function(targetExpression, intronToUnion, groups, psi = FALSE)
{

    # TODO: verify all 'introns' are in targetExpression and all target_ids in
    # targetExpression
    labs <- setdiff(colnames(targetExpression), 'target_id')

    if (length(groups) != length(labs))
        stop("length(groups) must be the same as the number of experiments included (and also in the same order)")

    if (psi)
    {
        # XXX: this is to add introns in the denom
        print("in here! (psi)")
        repIntrons <- data_frame(intron = unique(intronToUnion$intron),
            target_id = unique(intronToUnion$intron))
        intronToUnion <- rbind_list(intronToUnion, repIntrons)
    }

    denomExp <- left_join(intronToUnion, targetExpression, by = "target_id") %>%
        group_by(intron) %>%
        select(-(target_id)) %>%
        summarise_each(funs(sum)) %>%
        arrange(intron)

    numExp <- select(denomExp, intron) %>%
        inner_join(targetExpression, by = c("intron" = "target_id")) %>%
        arrange(intron)

    print("all.equal: ")
    print(all.equal(numExp$intron, denomExp$intron))

    denomExp <- as.data.frame(denomExp)
    # denomExp <- denomExp[order(denomExp$intron),]
    # rownames(denomExp) <- denomExp$intron
    denomExp$intron <- NULL

    numExp <- as.data.frame(numExp)
    # numExp <- numExp[order(numExp$intron),]
    # rownames(numExp) <- numExp$intron
    numExp$intron <- NULL

    retentionExp <- numExp / denomExp

    rownames(targetExpression) <- targetExpression$target_id
    targetExpression$target_id <- NULL

    nGroups <- length(unique(groups))

    validIntrons <- data.frame(
        matrix(rep(TRUE, nrow(retentionExp) * nGroups),
            ncol = nGroups))
    colnames(validIntrons) <- unique(groups)

    new("IntronRetention",
        retention = retentionExp,
        numerator = numExp,
        denominator = denomExp,
        labels = labs,
        groups = groups,
        features = targetExpression,
        validIntrons = validIntrons
        )
}

#' Filter low expression
#'
#' Filter out low expression features.
#'
#' @name lowExpressionFilter
#' @param obj the type
#' @param ... additional arguments
#' @rdname lowExpressionFilter-methods
#' @export lowExpressionFilter
#' @examples
#' lowExpressionFilter()
setGeneric(
    name = "lowExpressionFilter",
    def = function(obj, ...)
        standardGeneric("lowExpressionFilter"))

setMethod("lowExpressionFilter", signature("IntronRetention"),
    function(obj, lower = 0.25, ...)
    {
        curFilt <- complete.cases(obj@retention)

        gtq <- apply(obj@denominator, 2,
            function(col)
            {
                colGt0 <- col[col > 0]
                # removing duplicates gives you more
                colGt0 <- colGt0[!duplicated(colGt0)]
                q <- quantile(colGt0, probs = lower)
                col >= q
            })

        repGt <- lapply(unique(obj@groups), function(grp)
            {
                curGroups <- obj@groups %in% grp
                apply(data.frame(gtq[,curGroups]), 1, all)
            })
        repGt <- do.call(cbind, repGt)
        colnames(repGt) <- unique(obj@groups)

        obj@validIntrons <- as.data.frame(repGt)

        obj
    })


#' Single condition intron-retention test
#'
#' For a single condition, test whether or not the intron was retained
#'
#' @name retentionTest
#' @param obj the object
#' @param ... additional arguments
#' @rdname retentionTest-methods
#' @export retentionTest
setGeneric(
    name = "retentionTest",
    def = function(obj, ...)
        standardGeneric("retentionTest"))

setMethod("retentionTest", signature("IntronRetention"),
    function(obj, level = 0.5, conditions = c("all"), adjust = T)
    {
        if( sum(unique(obj@groups) %in% conditions) != length(conditions) &&
            conditions != c('all') )
            stop("Not all conditions were found in the original group
                definition.")

        if (conditions == 'all')
            conditions <- unique(obj@groups)

        # bs <- bootstrapMean(obj)
        retentionRes <- lapply(conditions, ret_test, obj, level, adjust)
        rbind_all(retentionRes)
    })


#' @export
ret_test <- function(cond, obj, level, adjust)
{
    print(adjust)
    whichSamp <- obj@groups %in% cond
    valid <- obj@validIntrons[,cond]
    bs <- bootstrapMean(obj)
    if (adjust)
    {
        samp_bs_mean <- apply(bs[,whichSamp], 2, mean, na.rm = T)
        bs_diff <- apply(obj@retention[valid,whichSamp], 2, mean) - samp_bs_mean
        # obj@retention[valid,whichSamp] <- t(t(matrix(obj@retention[valid,whichSamp])) - samp_bs_mean)
        # hello <- t(t(obj@retention[valid,whichSamp]) - bs_diff)
        obj@retention[valid,whichSamp] <- sweep(obj@retention[valid,whichSamp], 2, bs_diff, "-")
        # obj@retention[valid,whichSamp] <- t(apply(obj@retention[valid,whichSamp],
        #         2, pmax, 0))
        obj@retention[valid,whichSamp] <- t(apply(obj@retention[valid,whichSamp],
                2, pmax, 0))
    }

    retRes <- retentionTestSingleCond(obj@retention[, whichSamp], level)

    testStat <- rep(NA, length(valid))
    testStat[valid] <- retRes$testStat[valid]

    retention <- retRes$avg
    variance <- retRes$variance

    # lretRes <- retentionTestSingleCond(log(obj@retention[, whichSamp]), level)

    # log_retention <- lretRes$avg
    # log_variance <- lretRes$variance

    # data.frame(intron = rownames(obj@retention), condition = cond,
    #     testStat = testStat, retention = retention,
    #     variance = variance, log_retention = log_retention,
    #     log_variance = log_variance, stringsAsFactors = F)

    data.frame(intron = rownames(obj@retention), condition = cond,
        testStat = testStat, retention = retention,
        variance = variance, stringsAsFactors = F)
}

#' @export
retentionTestSingleCond <- function(retentionMat, level = 0.0, offset = 0.05)
{
    stopifnot(ncol(retentionMat) > 1)
    m <- apply(retentionMat, 1, mean)
    v <- apply(retentionMat, 1, var)
    v <- v + offset
    testStat <- (m - level) / sqrt(v)
    list(avg = m, variance = v, testStat = testStat)
}

#' Get summary statistics for an IntronRetention object
#'
#' Get summary statistics for an IntronRetention object.
#'
#' @name summaryStats
#' @param obj the object
#' @param ... additional arguments
#' @rdname summaryStats-methods
#' @export summaryStats
setGeneric(
    name = "summaryStats",
    def = function(obj, ...)
        standardGeneric("summaryStats"))

setMethod("summaryStats", signature("IntronRetention"),
    function(obj, conditions = c("all"))
    {
        if( sum(unique(obj@groups) %in% conditions) != length(conditions) &&
            conditions != c('all') )
            stop("Not all conditions were found in the original group
                definition.")

        if (length(conditions) == 1 && conditions == 'all')
            conditions <- unique(obj@groups)

        allDenom <- getExpCondition(obj, "denominator", conditions)
        allSummary <- mapply(function(cond, expDf)
            {
                m <- apply(expDf, 1, mean)
                v <- apply(expDf, 1, var)
                data.frame(intron = rownames(expDf), condition = cond, avg = m,
                    variance = v, stringsAsFactors = F, row.names = NULL)
            }, conditions, allDenom, SIMPLIFY = F)
        print(head(allSummary))
        as.data.frame(rbindlist(allSummary))
    })


#' Get experimental condition
#'
#' Get a filtered data.frame back when specifying a condition and slot
#'
#' @name getExpCondition
#' @param obj the object
#' @param ... additional arguments
#' @rdname getExpCondition-methods
#' @export getExpCondition
setGeneric(
    name = "getExpCondition",
    def = function(obj, ...)
        standardGeneric("getExpCondition"))


setMethod("getExpCondition", signature("IntronRetention"),
    function(obj, expType = "retention", conditions = c("all"))
    {
        stopifnot(expType %in% c("retention", "numerator", "denominator"))

        if( sum(unique(obj@groups) %in% conditions) != length(conditions) &&
            conditions != c('all') )
            stop("Not all conditions were found in the original group
                definition.")

        if (length(conditions) == 1 && conditions == 'all')
            conditions <- unique(obj@groups)

        allExp <- slot(obj, expType)
        lapply(conditions,
            function(cond)
            {
                valid <- obj@validIntrons[,cond]
                whichSamp <- obj@groups %in% cond
                allExp[valid, whichSamp]
            })
    })


#' @export
setGeneric(
    name = "bootstrapMean",
    def = function(obj, ...)
        standardGeneric("bootstrapMean"))

setMethod("bootstrapMean", signature("IntronRetention"),
    function(obj, nsamp = 2000)
    {
        denom <- obj@denominator - obj@numerator
        valid <- obj@validIntrons[,1]
        denom <- denom[valid,]
        num <- obj@numerator[valid,]

        denom_samp <- sample.int(nrow(denom), nsamp, replace = T)
        num_samp <- sample.int(nrow(num), nsamp, replace = T)

        num_mat <- num[num_samp,]
        denom_mat <- denom[denom_samp,] + num_mat

        num_mat / denom_mat
    })
