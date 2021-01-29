#'
#' Randomly generated roll call voting data set from a one-dimensional spherical latent factor model, which contains 20 rows(legislators) and 10 columns(bills).
#' The data set is binary and 1 corresponds to 'yea' and 0 corresponds to 'nay'
#' @docType data
#'
#' @usage data(ymat)
#'
#' @format An object of class \code{"matrix"}
#' \describe{
#' \item{Bill 1}{Generated voting record for bill 1}
#' \item{Bill 2}{Generated voting record for bill 2}
#' \item{Bill 3}{Generated voting record for bill 3}
#' \item{Bill 4}{Generated voting record for bill 4}
#' \item{Bill 5}{Generated voting record for bill 5}
#' \item{Bill 6}{Generated voting record for bill 6}
#' \item{Bill 7}{Generated voting record for bill 7}
#' \item{Bill 8}{Generated voting record for bill 8}
#' \item{Bill 9}{Generated voting record for bill 9}
#' \item{Bill 10}{Generated voting record for bill 10}
#' }
#' @references This data set was randomly generated for this package.
#' @keywords datasets
#' @examples
#'
#' data(ymat)
#' head(ymat)
"ymat"
