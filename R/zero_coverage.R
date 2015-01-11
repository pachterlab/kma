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
#'  \item{\code{target_id}}{the target name}
#'  \item{\code{len}}{the length of the target}
#'  \item{\code{zc_len}}{the length of the region with zero coverage}
#'  \item{zc_pval}{the raw p-value under the Poisson model}
#'  \item{zc_qval}{the corrected p-value using method \code{adjust_method}}
#' }
#' @export
get_intron_zc <- function(fname, exp_name, regex = "^chr",
    adjust_method = "bonferroni", rename_target_id = TRUE)
{
    zc_data <- as.data.frame(fread(fname, header = T))
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
        mutate(sample = exp_name)

    if (rename_target_id)
        zc_max <- rename(zc_max, "intron_extension" = target_id)

    zc_max
}


