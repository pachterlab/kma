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

    te <- data.table(targetExpression, key = 'target_id')
    i2u <- data.table(intronToUnion, key = 'target_id')
    # setkey(targetExpression, target_id)
    # setkey(intronToUnion, target_id)

    # print('bye')
    print(class(i2u))
    print(intronToUnion)
    print(targetExpression)
    denomExp <- te[i2u,,][,lapply(.SD, sum), by = intron]
    denomExp <- as.data.frame(denomExp)
    denomExp$target_id <- NULL
    # denomExp <- targetExpression[intronToUnion,,]
    # allDenomExp <- merge(targetExpression, intronToUnion)

    # TODO: change all the merges to joins using data.table
    # allDenomExp <- merge(intronToUnion, targetExpression,
    #                      by = 'target_id', all.x = T)

    # denomExp <- ddply(allDenomExp, .(intron),
    #                   function(df) {
    #                       apply(df[,expressionCols], 2, sum)
    #                   }, .parallel = TRUE)

    # eek

    expressionCols <- setdiff(colnames(targetExpression), 'target_id')

    rownames(denomExp) <- as.character(denomExp$intron)
    denomExp$intron <- NULL
    rownames(targetExpression) <- as.character(targetExpression$target_id)
    targetExpression$target_id <- NULL
    numExp <- targetExpression[rownames(denomExp),]

    retentionExp <- numExp / denomExp
    # groups <- factor(groups)

    new("IntronRetention",
        retention = retentionExp,
        numerator = numExp[rownames(retentionExp),],
        denominator = denomExp[rownames(retentionExp),],
        labels = labs,
        groups = groups,
        features = targetExpression,
        validIntrons = rep(TRUE, nrow(retentionExp))
        )
    # denomExp
}
