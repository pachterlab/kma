#' Get the columns from eXpress output
#'
#' Given a list of data.frames (all containing the same type of layout),
#'
#' @param aList a list of xprs data.frames
#' @param col the column you wish to extract
#' @param subStringReplace replace this substring from the target_id
#' @return a data.frame with the particular column from express output
#' @export
getCol <- function(aList, col, subStringReplace = NULL)
{
    # TODO: check that col is valid
    targets <- sort(aList[[1]]$target_id)

    df <- lapply(aList,
        function(x)
        {
            rownames(x) <- x$target_id
            x[targets, col]
        })
    class(df) <- 'data.frame'
    attr(df, 'row.names') <- 1:length(targets)
    attr(df, 'names') <- 1:ncol(df)
    df$targets <- targets

    if (!is.null(subStringReplace))
        df$target_id <- sub(subStringReplace, '',
                            sort(aList[[1]]$target_id))
    df
}
