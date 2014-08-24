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

    # TODO: refactor data.table -> data.frame
    retentionExp <- numExp[,!c('intron'),with = F] / denomExp[,!c('intron'),
        with = F]
    retentionExp[,intron := denomExp[,intron,],]
    setkey(retentionExp, intron)

    # TODO: write check ensuring in same ordering
    retentionExp <- as.data.frame(retentionExp)
    rownames(retentionExp) <- retentionExp$intron
    retentionExp$intron <- NULL

    numExp <- as.data.frame(numExp)
    rownames(numExp) <- numExp$intron
    numExp$intron <- NULL

    denomExp <- as.data.frame(denomExp)
    rownames(denomExp) <- denomExp$intron
    denomExp$intron <- NULL

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
