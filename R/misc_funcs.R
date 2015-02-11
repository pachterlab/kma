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

#' Convert intron coords string into data.frame
#'
#' Convert intron coords string into data.frame.
#'
#' @param introns a character vector from the \code{intron_to_transcript} table automatically generated (usually the "intron" column)
#' @return a \code{data.frame} with columns:
#' \describe{
#'  \item{chrom}{The reference name (usually chromosome).}
#'  \item{start}{The start position}
#'  \item{stop}{The stop posision}
#' }
#' @export
melt_intron_coords <- function(introns)
{
    if ( !is(introns, "character") )
        stop("Requires that 'introns' be a character vector")

    coords_list <- lapply(strsplit(introns, ":"), function(s)
        {
            ref <- as.character(s[1])
            coords <- as.integer(strsplit(s[2], "-")[[1]])
            data.frame(chrom = ref, start = coords[1],
                stop = coords[2], stringsAsFactors = F)
        })

    res <- rbind_all(coords_list)
    mutate(res, intron = introns)
}
