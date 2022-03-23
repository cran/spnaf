#' @importFrom magrittr %>%
#' @importFrom dplyr bind_rows

Gstat <- function(SpecialWeights, method, star = TRUE){

    t <- nrow(SpecialWeights)
    r_bar <- sum(SpecialWeights$n)/t
    s_sq <- sum((SpecialWeights$n - r_bar)**2)/(t-1)
    s <- sqrt(s_sq)

    SpecialWeightsl <- SpecialWeights %>%
        dplyr::mutate(seq = paste(oid, did, sep="-"))
    SWL <- split(SpecialWeightsl,
                 f = SpecialWeightsl$seq)

    oid <- did <- w <- NULL

    subframe <- function(l, ref = SpecialWeights, m = method){
        o <- l$oid
        d <- l$did

        if(m == 't'){
            origins <- ref %>%
                dplyr::filter(oid == o, w == 1) %>%
                dplyr::select(did) %>% unlist()
            destinations <- ref %>%
                dplyr::filter(oid == d, w == 1) %>%
                dplyr::select(did) %>% unlist()
        }else if(m == 'o'){
            origins <- o
            destinations <- ref %>%
                dplyr::filter(oid == d, w == 1) %>%
                dplyr::select(did) %>% unlist()
        }else if(m == 'd'){
            origins <- ref %>%
                dplyr::filter(oid == o, w == 1) %>%
                dplyr::select(did) %>% unlist()
            destinations <- d
        }else{
            stop("method must be one of c('t', 'o', 'd') \n")
        }
        set <- ref %>%
            dplyr::filter(oid %in% origins, did %in% destinations)

        ## Wij*
        Wij_star <- length(origins) * length(destinations)
        ## Wij*^2
        Wij_star_sq <- Wij_star**2
        ## S1
        S1 <- Wij_star
        ## Gij calculation
        sigma <- sum(set$n, na.rm = T)
        ## numerator
        numerator <- sigma - Wij_star*r_bar
        ## denominator
        denominator <- s * sqrt((t*S1 - Wij_star_sq)/(t-1))

        l$Gij <- numerator/denominator

        return(l)
    }
    G <- do.call("bind_rows", lapply(SWL, subframe))

    return(G)
}
