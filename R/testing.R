# Keep Me Around: Intron Retention Detection
# Copyright (C) 2015  Harold Pimentel <haroldpimentel@gmail.com>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#' Retention test
#'
#' Compute the retention test using resampling.
#'
#' @param obj an IntronRetention object
#' @param filtered_df a data.frame with columns intron, condition, and f_all.
#' If NULL, use the filters already computed and merge all of them.
#' @param n_samp number of samples to generate the null distribution
#' @return a \code{data.frame} with columns intron, condition, mean_retention,
#' var_retention, pvalue, and qvalue
#' @export
retention_test <- function(obj, filtered_df = NULL, n_samp = 10000)
{
    if (is.null(filtered_df)) {
        cat("aggregating filters\n")
        filtered_df <- aggregate_filters(obj)
    }
    else {
        # TODO: check whether or not all columns here
        stopifnot(inherits(filtered_df, "data.frame"))

        if (!(colnames(filtered_df) %in% c("intron", "condition", "f_all"))) {
            stop("require at least the following columns: 'intron', 'condition', 'f_all'")
        }
    }

    cat("joining filtered data set\n")
    tmp_df <- inner_join(obj$flat,
        select(filtered_df, intron, condition, f_all),
            by = c("intron", "condition"))

    stopifnot(nrow(tmp_df) == nrow(obj$flat))

    tmp_df <- tmp_df %>%
        filter(f_all)

    null_dists <- tmp_df %>%
        group_by(condition) %>%
        do({
            res <- intron_null_dist(., n_samp = n_samp)
            data_frame(
                samples = list(res$data),
                ecdf = list(res$ecdf)
                )
        }) %>%
        select(-c(samples))

    cat("computing p-values\n")

    tmp_df %>%
        group_by(intron, condition) %>%
        summarise(mean_retention = mean(retention),
            var_retention = var(retention)) %>%
        inner_join(null_dists, by = "condition") %>%
        mutate(pvalue = intron_pval(mean_retention, ecdf[[1]])) %>%
        select(-c(ecdf)) %>%
        ungroup() %>%
        group_by(condition) %>%
        mutate(qvalue = p.adjust(pvalue, method = "BH"))
}
