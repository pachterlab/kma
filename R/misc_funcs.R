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
