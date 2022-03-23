## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo = FALSE------------------------------------------------------
library(spnaf)
knitr::opts_chunk$set(warning = FALSE, message = FALSE) 

## -----------------------------------------------------------------------------
dim(CA)
head(CA)

## -----------------------------------------------------------------------------
library(sf)
plot(CA_polygon, col = 'white', main = 'CA polygon')

## -----------------------------------------------------------------------------
args(Gij.polygon)

## ---- warnings = FALSE--------------------------------------------------------
# Data manipulation
CA <- spnaf::CA
OD <- cbind(CA$FIPS.County.Code.of.Geography.B, CA$FIPS.County.Code.of.Geography.A)
OD <- cbind(OD, CA$Flow.from.Geography.B.to.Geography.A)
OD <- data.frame(OD)
names(OD) <- c("oid", "did", "n")
OD$n <- as.numeric(OD$n)
OD <- OD[order(OD[,1], OD[,2]),]
head(OD) # check the input df's format

# Load sf polygon
CA_polygon <- spnaf::CA_polygon
head(CA_polygon) # it has geometry column

# Execution of Gij.polygon with data above and given parameters
result <- Gij.polygon(df = OD, shape = CA_polygon, queen = TRUE, snap = 1,
method = 't', n = 1000)

## ---- eval = TRUE-------------------------------------------------------------
# positive clusters at the significance level of 0.05
head(result[[1]][result[[1]]$pval < 0.05,])
# positive clusters at the significance level of 0.05 in lines class
head(result[[2]][result[[2]]$pval < 0.05,])

## ---- warning = FALSE, fig.show = "hold", out.width = "45%"-------------------
library(tmap)
# plot all flows with the polygon (left)
tm_shape(CA_polygon) +
  tm_polygons()+
  tm_shape(result[[2]]) +
  tm_lines()
# plot significant flows only with the polygon (right)
tm_shape(CA_polygon) +
  tm_polygons()+
  tm_shape(result[[2]][result[[2]]$pval < 0.05,]) +
  tm_lines(col='pval')


