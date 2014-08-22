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

    print('join')
    allDenomExp <- left_join(intronToUnion, targetExpression, by = c('target_id'))
    print('summarize')
    denomExp <- allDenomExp %>% group_by(intron) %>% select(-target_id) %>% 
        summarise_each(c('sum'))

    print('rest')
    expressionCols <- setdiff(colnames(targetExpression), 'target_id')

    rownames(denomExp) <- as.character(denomExp$intron)
    denomExp$intron <- NULL
    rownames(targetExpression) <- as.character(targetExpression$target_id)
    targetExpression$target_id <- NULL
    numExp <- targetExpression[rownames(denomExp),]

    retentionExp <- numExp / denomExp
    # groups <- factor(groups)

    new("IntronRetention",
        retention = as.data.frame(retentionExp),
        numerator = as.data.frame(numExp[rownames(retentionExp),]),
        denominator = as.data.frame(denomExp[rownames(retentionExp),]),
        labels = labs,
        groups = groups,
        features = as.data.frame(targetExpression),
        validIntrons = rep(TRUE, nrow(retentionExp))
        )
    # denomExp
}
