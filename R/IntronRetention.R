# TODO: update docs below...

#' Compute intron retention
#'
#' Create an IntronRetention object.
#'
#' @param targetExpression a matrix of nothing but expression values and one
#' column containing 'target_id's
#' @param intronToUnion a table with columns (all character vectors):
#' \describe{
#'  \item{intron}{the identifier of the intron}
#'  \item{target_id}{the target_id of transcripts compatible with the intron}
#'  \item{gene}{the gene name that this intron belongs to}
#'  \item{intron_extension}{the actual coordinates of the intron quantified
#'  (including the region that overlaps the intronic region)}
#' }
#' @param groups a vector with the grouping
#' @param psi if TRUE, compute the psi value. otherwise, compute a rate
#' @return an IntronRetention object
#' @export
newIntronRetention <- function(targetExpression, intronToUnion, groups, psi = TRUE)
{
    targetExpression <- as.data.frame(targetExpression, stringsAsFactors = F)
    intronToUnion <- as.data.frame(intronToUnion, stringsAsFactors = F)

    # TODO: verify all 'introns' are in targetExpression and all target_ids in
    # targetExpression
    labs <- setdiff(colnames(targetExpression), 'target_id')

    if (length(groups) != length(labs))
        stop("length(groups) must be the same as the number of experiments included (and also in the same order)")

    if (psi)
    {
        # XXX: this is to add introns in the denom
        # repIntrons <- data_frame(intron = unique(intronToUnion$intron),
        #     target_id = unique(intronToUnion$intron))
        repIntrons <- intronToUnion %>%
            select(intron, gene, intron_extension) %>%
            distinct() %>%
            mutate(target_id = intron_extension)
        intronToUnion <- rbind_list(intronToUnion, repIntrons)
    }

    cat("computing denominator\n")
    denomExp <- left_join(intronToUnion, targetExpression, by = "target_id") %>%
        group_by(intron) %>%
        select(-(target_id)) %>%
        summarise_each(funs(sum), -matches("gene"),
            -matches("intron_extension")) %>%
        arrange(intron) %>%
        left_join(
            select(intronToUnion, intron, intron_extension) %>%
                distinct(),
            by = c("intron"))

    cat("computing numerator\n")
    numExp <- select(denomExp, intron, intron_extension) %>%
        inner_join(targetExpression, by = c("intron_extension" = "target_id")) %>%
        arrange(intron_extension)

    denomExp <- as.data.frame(denomExp) %>%
        arrange(intron)
    rownames(denomExp) <- denomExp$intron
    denomExp <- select(denomExp, -c(intron, intron_extension))

    numExp <- as.data.frame(numExp) %>%
        arrange(intron)
    rownames(numExp) <- numExp$intron
    numExp <- select(numExp, -c(intron, intron_extension))

    cat("computing retention\n")
    retentionExp <- numExp / denomExp

    rownames(targetExpression) <- targetExpression$target_id
    targetExpression$target_id <- NULL

    cat("'melting' expression\n")
    flat = melt_retention(retentionExp, numExp, denomExp, groups)

    cat("sorting and grouping by (intron, condition)\n")
    flat <- flat %>%
        arrange(intron, condition) %>%
        group_by(intron, condition)

    structure(list(retention = retentionExp,
            numerator = numExp,
            denominator = denomExp,
            labels = labs,
            groups = groups,
            features = targetExpression,
            flat = flat
            ),
        class = "IntronRetention")
}

#' @export
melt_retention <- function(ret, num, denom, groupings)
{
    samp_to_condition <- data.frame(sample = colnames(ret),
        condition = groupings, stringsAsFactors = FALSE)

    ret <- ret %>% mutate(intron = rownames(ret))

    ret <- melt(ret, id.vars = "intron",
        variable.name = "sample",
        value.name = "retention")


    denom <- denom %>% mutate(intron = rownames(denom))
    denom <- melt(denom, id.vars = "intron",
        variable.name = "sample",
        value.name = "denominator")

    num <- num %>% mutate(intron = rownames(num))
    num <- melt(num, id.vars = "intron",
        variable.name = "sample",
        value.name = "numerator")

    m_res <- inner_join(num, denom, by = c("intron", "sample")) %>%
        inner_join(ret, by = c("intron", "sample"))

    left_join(m_res, samp_to_condition, by = "sample")
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

#' @export
lowExpressionFilter_upd <- function(obj, lower = 0.25, ...)
{
    curFilt <- complete.cases(obj@retention)

    denom_minus <- obj@denominator - obj@numerator

    gtq <- apply(denom_minus, 2,
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
}


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
retentionTestSingleCond <- function(retentionMat, level = 0.0, offset = 0.00)
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

#' @export
lowTpmFilter <- function(obj, min_tpm)
{
    rep_filt <- lapply(unique(obj@groups), function(cond)
        {
            whichCond <- obj@groups %in% cond
            apply(obj@denominator[,whichCond], 1, function(row) all(row >= min_tpm))
        })
    rep_filt <- do.call(cbind, rep_filt)
    colnames(rep_filt) <- unique(obj@groups)
    obj@validIntrons <- as.data.frame(obj@validIntrons & rep_filt)

    obj
}

#' @export
lowTpmFilter_upd <- function(obj, min_tpm)
{
    rep_filt <- lapply(unique(obj@groups), function(cond)
        {
            whichCond <- obj@groups %in% cond
            apply(obj@denominator[,whichCond] - obj@numerator[,whichCond], 1, function(row) all(row >= min_tpm))
        })
    rep_filt <- do.call(cbind, rep_filt)
    colnames(rep_filt) <- unique(obj@groups)
    obj@validIntrons <- as.data.frame(obj@validIntrons & rep_filt)

    obj
}

#' @export
lowCountFilter <- function(obj, counts, n)
{
    counts <- data_frame(target_id = rownames(obj@retention)) %>%
        inner_join(counts, by = "target_id") %>%
        arrange(target_id) %>% select(-(target_id))

    rep_filt <- lapply(unique(obj@groups), function(cond)
        {
            whichCond <- obj@groups %in% cond
            apply(counts[,whichCond], 1, function(row) all(row > n));
        })
    rep_filt <- do.call(cbind, rep_filt)
    colnames(rep_filt) <- unique(obj@groups)
    obj@validIntrons <- as.data.frame(obj@validIntrons & rep_filt)

    obj
}

#' @export
perfectPsiFilter <- function(obj)
{
    rep_mean <- lapply(unique(obj@groups), function(grp)
        {
            whichCond <- obj@groups %in% grp

            mean_ret <- apply(obj@retention[,whichCond], 1, mean)
            round_mean <- round(mean_ret, 6)

            tmp <- round_mean != 1.0 & round_mean != 0.0 & obj@validIntrons[,grp]
            tmp[is.na(tmp)] <- FALSE
            tmp
        })
    tmp_names <- colnames(obj@validIntrons)
    obj@validIntrons <- data.frame(do.call(cbind, rep_mean))
    colnames(obj@validIntrons) <- tmp_names

    obj
}

#' Get intron lengths from an identifier
#'
#' Transform an intron identifier into lengths.
#'
#' @param intron_names the names of the introns in format 'chrom:start-stop'
#' @return an integer vector of the lengths
#' @export
intron_length <- function(intron_names)
{
    unlist(lapply(strsplit(intron_names, ":"), function(x)
        {
            coords <- as.integer(strsplit(x[2], '-')[[1]])
            coords[2] - coords[1]
        }))
}

#' Generate the null distribution
#'
#'
#' @param grouped_df
#' @return a numeric with means
#' @export
intron_null_dist <- function(grouped_df, n_samp = 10000, test_stat = mean)
{
    all_dat <- dcast(grouped_df, intron ~ sample, value.var = "retention")
    all_dat <- select(all_dat, -c(intron))
    all_dat <- all_dat[complete.cases(all_dat),]

    samps <- as.matrix(as.data.frame(lapply(all_dat, sample,
                n_samp,
                replace = T)))

    data <- apply(samps, 1, test_stat)

    list(data = data, ecdf = ecdf(data))
}

intron_pval <- function(mean_val, null_dist)
{
    1 - null_dist$ecdf(meal_val)
}

#' @export
low_tpm_filter <- function(obj, tpm, filter_name = paste0("low_tpm_", round(tpm, 2)))
{
    # TODO: check if grouping exists
    obj$flat <- obj$flat %>%
        mutate_(.dots = setNames(list(~all((denominator - numerator) >= tpm)),
                c(filter_name)))
    obj
}

check_groupings <- function(dat, valid_groups = c("intron", "condition"))
{
    grouping_valid <- identical(sort(as.character(groups(dat))),
        sort(valid_groups))

    if (!grouping)
    {
        # TODO: group them!
    }

    dat
}

perfect_zero_filter <- function(obj)
{

}
