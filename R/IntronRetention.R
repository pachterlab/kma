setClass("IntronRetention",
         slots = list(retention = "data.table",
                      numerator = "data.table",
                      denominator = "data.table",
                      labels = "character",
                      groups = "factor",
                      features = "data.table",
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

    print('join')
    denomExp <- te[i2u][,!'target_id', with = F][,lapply(.SD, sum), by = intron]
    setkey(denomExp, intron)

    print('rest')

    # get just the expression of the introns
    numExp <- te[ denomExp[,intron,] ]
    setnames(numExp, c('target_id'), c('intron'))
    setkey(numExp, intron)

    retentionExp <- numExp[,!c('intron'),with = F] / denomExp[,!c('intron'),with = F]
    retentionExp[,intron := denomExp[,intron,],]
    setkey(retentionExp, intron)

    new("IntronRetention",
        retention = retentionExp,
        numerator = numExp,
        denominator = denomExp,
        labels = labs,
        groups = groups,
        features = te,
        validIntrons = rep(TRUE, nrow(retentionExp))
        )
}
