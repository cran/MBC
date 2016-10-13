\name{cccma}
\alias{cccma}
\title{
Sample CanESM2 and CanRCM4 data
}
\description{
Sample CanESM2 (T63 grid) and CanRCM4 (0.22-deg grid) data
(122.5 deg W, 50 deg N).
\preformatted{pr: precipitation (mm day-1) 
tas: average surface temperature (deg. C)
dtr: diurnal temperature range (deg. C)
sfcWind: surface wind speed (m s-1)
ps: surface pressure (ps)
huss: surface specific humidity (kg kg-1)
rsds: surface downwelling shortwave radiation (W m-2)
rlds: surface downwelling longwave radiation (W m-2)}
}
\value{
     a list of with elements consisting of:
      \item{gcm.c}{matrix of CanESM2 variables for the calibration period.}
      \item{gcm.p}{matrix of CanESM2 variables for the validation period.}
      \item{rcm.c}{matrix of CanRCM4 variables for the calibration period.}
      \item{rcm.p}{matrix of CanRCM4 variables for the validation period.}
      \item{ratio.seq}{vector of logical values indicating if samples are of a ratio quantity.}
      \item{trace}{numeric values indicating trace thresholds for each ratio quantity.}
}

