#' @importFrom magrittr %>%

SpatialWeight <- function(df, shape, snap, queen){

    w <- NULL
    nb <- unclass(spdep::poly2nb(shape, snap, queen))
    original_id <- shape$id
    nb <- lapply(nb, function(data){
        return(original_id[data])
    })
    names(nb) <- original_id

    did <- unlist(nb)
    oid <- substr(names(did), 1, nchar(names(did))-1)
    nb_frame <- data.frame(oid = oid)
    nb_frame$did <- did
    nb_frame$w <- 1
    result <- df %>%
        dplyr::left_join(nb_frame, by = c("oid" = "oid", "did" = "did")) %>%
        dplyr::mutate(w = ifelse(oid == did, 1, w)) %>%
        dplyr::arrange(oid, did)

    return(result)

}
