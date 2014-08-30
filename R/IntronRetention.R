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
#' @return a IntronRetention object
#' @export
newIntronRetention <- function(targetExpression, intronToUnion, groups)
{

    # TODO: verify all 'introns' are in targetExpression and all target_ids in
    # targetExpression
    labs <- setdiff(colnames(targetExpression), 'target_id')

    if (length(groups) != length(labs))
        stop("length(groups) must be the same as the number of experiments included (and also in the same order)")

    te <- data.table(targetExpression, key = 'target_id')
    i2u <- data.table(intronToUnion, key = 'target_id')

    denomExp <- te[i2u][,!'target_id', with = F][,lapply(.SD, sum), by = intron]
    setkey(denomExp, intron)

    # get just the expression of the introns
    numExp <- te[ denomExp[,intron,] ]
    setnames(numExp, c('target_id'), c('intron'))
    setkey(numExp, intron)

    denomExp <- as.data.frame(denomExp)
    rownames(denomExp) <- denomExp$intron
    denomExp$intron <- NULL

    numExp <- as.data.frame(numExp)
    rownames(numExp) <- numExp$intron
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
    function(obj, level = 0.5, conditions = c("all"))
    {
        if( sum(unique(obj@groups) %in% conditions) != length(conditions) &&
            conditions != c('all') )
            stop("Not all conditions were found in the original group
                definition.")

        if (conditions == 'all')
            conditions <- unique(obj@groups)

        retentionRes <- lapply(conditions,
            function(cond)
            {
                whichSamp <- obj@groups %in% cond
                valid <- obj@validIntrons[,cond]
                retRes <- retentionTestSingleCond(obj@retention[valid, whichSamp], level)

                testStat <- rep(NA, length(valid))
                testStat[valid] <- retRes$testStat

                retention <- rep(NA, length(valid))
                retention[valid] <- retRes$avg

                variance <- rep(NA, length(valid))
                variance[valid] <- retRes$variance

                data.frame(intron = rownames(obj@retention), condition = cond,
                    testStat = testStat, retention = retention,
                    variance = variance)

            })
        as.data.frame(rbindlist(retentionRes))
    })


#' @export
retentionTestSingleCond <- function(retentionMat, level = 0.0)
{
    stopifnot(ncol(retentionMat) > 1)
    m <- apply(retentionMat, 1, mean)
    v <- apply(retentionMat, 1, var)
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
