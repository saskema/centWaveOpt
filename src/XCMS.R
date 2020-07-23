## This program is free software: you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation, either version 3 of the License, or (at your option) any later
## version.

## This program is distributed in the hope that it will be useful, but WITHOUT ANY
## WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
## PARTICULAR PURPOSE.  See the GNU General Public License for more details.

## You should have received a copy of the GNU General Public License along with
## this program.  If not, see https://www.gnu.org/licenses/.


## Enter your parameters for the peak picking algorithm as received by
## centWaveOpt.R
parameter <- c(5, 20, 5, 4, 0.01, 3, 100, 5)  # XCMS-Online preset

## Load necessary packages
library("xcms")
library("CAMERA")

## Create directory for output files
dir.create("XCMSResults", showWarnings = FALSE)

## Load file paths
files <- list.files(path = getwd(), recursive = TRUE, full.names = TRUE, pattern = ".mzXML")

## Preprocessing
set <- xcmsSet(files, method = "centWave", peakwidth = c(parameter[1], parameter[2]),
               ppm = parameter[3], snthresh = parameter[4], mzdiff = parameter[5],
               prefilter = c(parameter[6], parameter[7]))
set <- group.density(set, bw = parameter[8])
set.retcor <- retcor(set, method = "obiwarp", plottype = "none")
set.retcor <- group.density(set.retcor, bw = parameter[8])
set.retcor.2 <- retcor(set.retcor, method = "obiwarp", plottype = "none")
set.retcor.2 <- group.density(set.retcor.2, bw = parameter[8])
set.filled <- fillPeaks(set.retcor.2)

## Annotation of isotopes and adducts
set.CAMERA <- CAMERA::annotate(set.filled, polarity = "positive")

## Export feature table for identification of analytes
features <- getPeaklist(set.CAMERA)
rownames(features) <- groupnames(set.filled)
write.csv(features, "XCMSResults/features.csv")
