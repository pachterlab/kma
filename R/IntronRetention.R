setClass("IntronRetention",
         slots = list(retention = "data.frame",
                      numerator = "data.frame",
                      denominator = "data.frame",
                      labels = "character",
                      groups = "factor",
                      features = "data.frame",
                      validIntrons = "logical"))

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

    allDenomExp <- merge(intronToUnion, targetExpression,
                         by = 'target_id', all.x = T)
    expressionCols <- setdiff(colnames(targetExpression), 'target_id')
    denomExp <- ddply(allDenomExp, .(intron),
                      function(df) {
                          apply(df[,expressionCols], 2, sum)
                      }, .parallel = TRUE)
    rownames(denomExp) <- denomExp$intron
    denomExp$intron <- NULL
    rownames(targetExpression) <- targetExpression$target_id
    targetExpression$target_id <- NULL
    numExp <- targetExpression[rownames(denomExp),]

    retentionExp <- numExp / denomExp
    groups <- factor(groups)

    new("IntronRetention",
        retention = retentionExp,
        numerator = numExp[rownames(retentionExp),],
        denominator = denomExp[rownames(retentionExp),],
        labels = labs,
        groups = groups,
        features = targetExpression,
        validIntrons = rep(TRUE, nrow(retentionExp))
        )
}
