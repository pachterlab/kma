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

#' Low TPM filter
#'
#'
#' @export
filter_low_tpm <- function(obj, tpm, filter_name = paste0("f_low_tpm_", round(tpm, 2)))
{
    # TODO: check if grouping exists
    obj$flat <- check_groupings(obj$flat, c("intron", "condition"))

    obj$flat <- obj$flat %>%
        mutate_(.dots = setNames(list(~all((denominator - numerator) >= tpm)),
                c(filter_name)))
    obj
}

#' Perfect psi filter
#'
#' Remove things that have "perfect" psi scores. A perfect score is when all
#' samples contain exactly 0 or 1 retention after rounding.
#'
#'
#' @export
filter_perfect_psi <- function(obj, digits = 5, filter_name = "f_perf_psi")
{
    # TODO: check if psi is actually calculated
    obj$flat <- check_groupings(obj$flat, c("intron", "condition"))

    obj$flat <- obj$flat %>%
        mutate(round_ret = round(retention, digits)) %>%
        mutate_(.dots = setNames(list(~(!(
                            all(round_ret == 0.0) ||
                            all(round_ret == 1.0) ||
                            any(is.nan(retention))
                            ))),
                c(filter_name))) %>%
        select(-c(round_ret))


    obj
}


#' Filter introns with too few unique fragments.
#'
#' Filter introns with too few unique fragments.
#'
#' @param obj IntronRetention object
#' @param min_frags the minimum number of fragments required to pass the filter
#' in every sample
#' @param filter_name a character string denoting the column name in the final
#' table
#' @return an IntronRetention object with \code{flat} containing a new column
#' named \code{filter_name}
#' @export
filter_low_frags <- function(obj, min_frags,
    filter_name = paste0("f_low_count_", min_frags))
{
    # TODO: can refactor this into an independent calculation
    if (is.null(obj$unique_counts)) {
        stop("Please recreate the retention object with unique counts")
    }


    # require grouping by intron, condition
    obj$flat <- check_groupings(obj$flat, c("intron", "condition"))

    obj$flat <- obj$flat %>%
        mutate_(.dots = setNames(list(~all(unique_counts >= min_frags)),
                c(filter_name)))

    obj
}

#' @export
aggregate_filters <- function(obj)
{
    tmp_summarise <- obj$flat %>%
        group_by(intron, condition) %>%
        select(starts_with("f_")) %>%
        summarise_each(funs(all)) %>%
        ungroup()

    all_filt <- tmp_summarise %>%
        select(starts_with("f_")) %>%
        t() %>%
        apply(2, all)

    tmp_summarise <- tmp_summarise %>%
        mutate(f_all = all_filt)

        # ungroup() %>%
        # rowwise() %>%
        # do({
        #     filt_vals <- .[dplyr:::starts_with(names(.), "f_")]
        #     f_all <- all(unlist(filt_vals))
        #     data.frame(., f_all = f_all, stringsAsFactors = FALSE)
        # }) %>%
        # ungroup()
    tmp_summarise
}
