#' Sample migration data by counties in California.
#'
#' A dataframe containing migration of CA counties with origins and destinations
#' from US Census
#'
#' @source
#' Census.gov > Population > Migration/Geographic Mobility > Guidance for Data Users > County-to-County Migration Flows
#' https://www.census.gov/topics/population/migration/guidance/county-to-county-migration-flows.html
#'
#'
#' @format A \code{data.frame} object with 2580 rows and 12 variables
#' \describe{
#'   \item{State.Code.of.Geography.A}{Destinations' State code}
#'   \item{FIPS.County.Code.of.Geography.A}{Destinations' FIPS County code}
#'   \item{State.U.S..Island.Area.Foreign.Region.Code.of.Geography.B}{Destinations' State U.S. Island Area Foreign Region Code}
#'   \item{FIPS.County.Code.of.Geography.B}{Origins' FIPS County code}
#'   \item{State.Name.of.Geography.A}{Destinations' State name}
#'   \item{County.Name.of.Geography.A}{Destinations' County name}
#'   \item{State.U.S..Island.Area.Foreign.Region.of.Geography.B}{Origins' State U.S. Island Area Foreign Region Code}
#'   \item{County.Name.of.Geography.B}{Origins' County name}
#'   \item{Flow.from.Geography.B.to.Geography.A}{Flow count from the origin to the destination}
#'   \item{Counterflow.from.Geography.A.to.Geography.B}{Counterflow count from the destination to the origin}
#'   \item{Net.Migration.from.Geography.B.to.Geography.A}{Net migration count from the origin to the destination}
#'   \item{Gross.Migration.between.Geography.A.and.Geography.B}{Gross migration count between counties}
#' }
#'
#' @examples
#' CA
"CA"
