#' Read zero coverage data
#'
#' Read the zero coverage data otuput into a \code{data.frame}.
#'
#' @param fname a character vector of length one containing the file name
#' @param exp_name a character containing the name of the experiment
#' @param regex a regular expression to filter from. For intron retention, you
#' usually want to only get the introns. If you want to get everything back,
#' use the regex \code{"*"}.
#' @param adjust_method the method to use for multiple testing
#' @param rename_target_id if TRUE, rename 'target_id' to 'intron_extension'
#' @return a data.frame with columns:
#' \describe{
#'  \item{\code{intron_extension}}{the intron_extension coordinates}
#'  \item{\code{len}}{the length of the target}
#'  \item{\code{zc_len}}{the length of the region with zero coverage}
#'  \item{zc_pval}{the raw p-value under the Poisson model}
#'  \item{zc_qval}{the corrected p-value using method \code{adjust_method}}
#' }
#' @export
get_intron_zc <- function(fname, sample_name, regex = "^chr",
    adjust_method = "bonferroni", rename_target_id = TRUE)
{
    zc_data <- fread(fname, header = TRUE, data.table = FALSE)
    zc_max <- filter(zc_data, grepl(regex, reference)) %>%
        mutate(len = end - start) %>%
        group_by(reference) %>%
        filter(len == max(len)) %>%
        distinct() %>%
        ungroup()

    zc_max <- zc_max %>%
        mutate(pval = exp(pvalue),
            qval = p.adjust(pval, method = adjust_method)) %>%
        select(reference, len, pval, qval) %>%
        rename("target_id" = reference, "zc_len" = len, "zc_pval" = pval,
            "zc_qval" = qval) %>%
        mutate(sample = sample_name)

    if (rename_target_id)
        zc_max <- rename(zc_max, "intron_extension" = target_id)

    zc_max
}


#' Get several zero coverage results and aggregate them
#'
#' Takes a list of zero coverage results and aggregates them into a data.frame
#' that can be used by an IntronRetention class as a filter.
#'
#' @param fnames a character vector containing the path to files
#' @param sample_names a character vector containing names for each sample. Note: this should correspond to names that are in the IntronRetention table
#' @param condition_names this should correspond to the coditions used in newIntronRetention
#' @return a \code{data.frame} containing all of the samples
#' @export
get_batch_intron_zc <- function(fnames,
    sample_names,
    condition_names,
    regex = "^chr")
{
    if (!is(sample_names, "character")) {
        stop("'sample_names' must be a character vector")
    }
    if (!is(condition_names, "character")) {
        stop("'condition_names' must be a character vector")
    }
    if (length(sample_names) != length(condition_names)) {
        stop("sample_names length must be equal to condition_names")
    }

    sample_to_conditions <- data.frame(sample = sample_names,
        condition = condition_names)

    # call get_intron_zc on all the samples and dump them into one table
    all_zc <- Map(get_intron_zc, fnames, sample_names, regex) %>%
        rbind_all()

    # add in the condition names
    all_zc <- all_zc %>%
        left_join(sample_to_conditions, by = "sample")

    all_zc
}

#' Summarize zero coverage results per condition
#'
#' Check whether all zero coverage results agree per intron and condition
#'
#' @param ir an \code{IntronRetention} object
#' @param zc a zero coverage \code{data.frame}
#' @param prop_length the maximum proportion of the intron length the zero
#' coverage region should span before returning FALSE
#' @param level if the adjust p-value is less than \code{level}, then the
#' intron will fail the filter
#' @return an \code{IntronRetention} object with the data.frame updated to
#' include the zero coverage filter
#' @export
summarize_zero_coverage <- function(ir, zc, prop_length = 0.20, level = 0.05)
{
    stopifnot(is(ir, "IntronRetention"))

    # TODO: think about what the requirements are (e.g. conditions and stuff
    # matching w/ IR

    # TODO: accept prop_length as a vector and return several different
    # filter in one call
    zc <- inner_join(zc, ir$intron_to_extension, by = "intron_extension")

    filter_name <- paste0("f_zc_", prop_length)

    zc <- zc %>%
        group_by(intron, condition) %>%
        summarise_(.dots = setNames(
                list(~all(zc_qval > level |
                        zc_len < (prop_length * extension_len))),
                    c(filter_name)))

    ir$flat <- ir$flat %>%
        left_join(zc, by = c("intron", "condition")) %>%
        ungroup()

    ir$flat[filter_name] <- replace(ir$flat[filter_name],
        is.na(ir$flat[filter_name]), TRUE)

    ir
}
