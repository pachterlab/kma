#' Get the columns from eXpress output
#'
#' Given a list of data.frames (all containing the same type of layout),
#'
#' @param aList a list of xprs data.frames
#' @param col the column you wish to extract
#' @param substring_replace replace this substring from the target_id
#' @return a data.frame with the particular column from express output
#' @export
getCol <- function(aList, col, column_ids, substring_replace = NULL)
{
    df <- data.frame(lapply(aList, function(x) x[col]))

    setnames(df, colnames(df), column_ids)

    df$target_id <- aList[[1]]$target_id

    if (!is.null(substring_replace))
        df$target_id <- sub(substring_replace, '', aList[[1]]$target_id)

    df
}

#' @export
read_express <- function(file_names, sample_names, condition_names,
    cols = c("tpm", "uniq_counts"), substring_replace = NULL)
{
    cat("Reading all eXpress data...\n")
    all_xprs <- lapply(file_names, fread, header = TRUE, stringsAsFactors = FALSE,
        data.table = FALSE)

    cat("Sorting...\n")
    all_xprs <- lapply(all_xprs, arrange, target_id)

    res <- lapply(cols, function(col){
        cat("Extracting column: ", col, "\n")
        getCol(all_xprs, col, sample_names, substring_replace)
        })

    names(res) <- cols
    res$all_data <- all_xprs
    res$condition <- condition_names
    res$sample <- sample_names

    res
}
