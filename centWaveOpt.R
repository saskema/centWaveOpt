## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see https://www.gnu.org/licenses/.


## Place XCMS readable files in your working directory and enter your file
## extension!
format <- ".mzXML"

## Load necessary packages
library("tidyverse")
library("xcms")

## Load file paths and directory names
files <- list.files(path = "QC", recursive = TRUE, full.names = TRUE, pattern = format)

## Create directories for output files
dir.create("centWaveOptResults", showWarnings = FALSE)

## Generate matrix with preprocessing parameters
preset <- c(5, 20, 5, 4, 0.01, 3, 100, 5)

peakwidth.min <- seq(from = 1, to = 10, by = 0.1)
peakwidth.max <- seq(from = 10, to = 100, by = 1)
ppm <- seq(from = 1, to = 2.5, by = 0.1)
snthresh <- seq(from = 1, to = 100, by = 1)
mzdiff <- seq(from = -0.1, to = 0.1, by = 0.002)
prefilter.scan.number <- seq(from = 1, to = 100, by = 1)
prefilter.scan.abundance <- seq(from = 100, to = 10000, by = 100)
bandwidth <- seq(from = 0.1, to = 1, by = 0.1)

## Create one matrices with preset parameters for every paramter to optimise
parameters.peakwidth.min <- matrix(data = preset, ncol = 8, nrow = length(peakwidth.min), 
    byrow = TRUE, dimnames = list(NULL, c("Minimum Peakwidth", "Maximum Peakwidth", 
        "ppm", "S/N Threshold", "m/z Difference", "Prefilter 1", "Prefilter 2", "bandwidth")))

parameters.peakwidth.max <- matrix(data = preset, ncol = 8, nrow = length(peakwidth.max), 
    byrow = TRUE, dimnames = list(NULL, c("Minimum Peakwidth", "Maximum Peakwidth", 
        "ppm", "S/N Threshold", "m/z Difference", "Prefilter 1", "Prefilter 2", "bandwidth")))

parameters.ppm <- matrix(data = preset, ncol = 8, nrow = length(ppm), byrow = TRUE, 
    dimnames = list(NULL, c("Minimum Peakwidth", "Maximum Peakwidth", "ppm", "S/N Threshold", 
        "m/z Difference", "Prefilter 1", "Prefilter 2", "bandwidth")))

parameters.snthresh <- matrix(data = preset, ncol = 8, nrow = length(snthresh), byrow = TRUE, 
    dimnames = list(NULL, c("Minimum Peakwidth", "Maximum Peakwidth", "ppm", "S/N Threshold", 
        "m/z Difference", "Prefilter 1", "Prefilter 2", "bandwidth")))

parameters.mzdiff <- matrix(data = preset, ncol = 8, nrow = length(mzdiff), byrow = TRUE, 
    dimnames = list(NULL, c("Minimum Peakwidth", "Maximum Peakwidth", "ppm", "S/N Threshold", 
        "m/z Difference", "Prefilter 1", "Prefilter 2", "bandwidth")))

parameters.prefilter.scan.number <- matrix(data = preset, ncol = 8, nrow = length(prefilter.scan.number), 
    byrow = TRUE, dimnames = list(NULL, c("Minimum Peakwidth", "Maximum Peakwidth", 
        "ppm", "S/N Threshold", "m/z Difference", "Prefilter 1", "Prefilter 2", "bandwidth")))

parameters.prefilter.scan.abundance <- matrix(data = preset, ncol = 8, nrow = length(prefilter.scan.abundance), 
    byrow = TRUE, dimnames = list(NULL, c("Minimum Peakwidth", "Maximum Peakwidth", 
        "ppm", "S/N Threshold", "m/z Difference", "Prefilter 1", "Prefilter 2", "bandwidth")))

parameters.bandwidth <- matrix(data = preset, ncol = 8, nrow = length(bandwidth), 
    byrow = TRUE, dimnames = list(NULL, c("Minimum Peakwidth", "Maximum Peakwidth", 
        "ppm", "S/N Threshold", "m/z Difference", "Prefilter 1", "Prefilter 2", "bandwidth")))

## Replace preset parameters in matrix with parameter sweeping parameters
for (i in 1:length(peakwidth.min)) {
    parameters.peakwidth.min[i, 1] <- peakwidth.min[i]
}

for (i in 1:length(peakwidth.max)) {
    parameters.peakwidth.max[i, 2] <- peakwidth.max[i]
}

for (i in 1:length(ppm)) {
    parameters.ppm[i, 3] <- ppm[i]
}

for (i in 1:length(snthresh)) {
    parameters.snthresh[i, 4] <- snthresh[i]
}

for (i in 1:length(mzdiff)) {
    parameters.mzdiff[i, 5] <- mzdiff[i]
}

for (i in 1:length(prefilter.scan.number)) {
    parameters.prefilter.scan.number[i, 6] <- prefilter.scan.number[i]
}

for (i in 1:length(prefilter.scan.abundance)) {
    parameters.prefilter.scan.abundance[i, 7] <- prefilter.scan.abundance[i]
}

for (i in 1:length(bandwidth)) {
    parameters.bandwidth[i, 8] <- bandwidth[i]
}
## Initialize matrix for description of density distribution
results.peakwidth.min <- cbind(parameters.peakwidth.min, matrix(nrow = length(peakwidth.min), 
    ncol = 7, dimnames = list(NULL, c("cvmin", "cv0.25", "cvmed", "cv0.75", "cvmax", 
        "cviqr", "time"))))
results.peakwidth.max <- cbind(parameters.peakwidth.max, matrix(nrow = length(peakwidth.max), 
    ncol = 7, dimnames = list(NULL, c("cvmin", "cv0.25", "cvmed", "cv0.75", "cvmax", 
        "cviqr", "time"))))
results.ppm <- cbind(parameters.ppm, matrix(nrow = length(ppm), ncol = 7, dimnames = list(NULL, 
    c("cvmin", "cv0.25", "cvmed", "cv0.75", "cvmax", "cviqr", "time"))))
results.snthresh <- cbind(parameters.snthresh, matrix(nrow = length(snthresh), ncol = 7, 
    dimnames = list(NULL, c("cvmin", "cv0.25", "cvmed", "cv0.75", "cvmax", "cviqr", 
        "time"))))
results.mzdiff <- cbind(parameters.mzdiff, matrix(nrow = length(mzdiff), ncol = 7, 
    dimnames = list(NULL, c("cvmin", "cv0.25", "cvmed", "cv0.75", "cvmax", "cviqr", 
        "time"))))
results.prefilter.scan.number <- cbind(parameters.prefilter.scan.number, matrix(nrow = length(prefilter.scan.number), 
    ncol = 7, dimnames = list(NULL, c("cvmin", "cv0.25", "cvmed", "cv0.75", "cvmax", 
        "cviqr", "time"))))
results.prefilter.scan.abundance <- cbind(parameters.prefilter.scan.abundance, matrix(nrow = length(prefilter.scan.abundance), 
    ncol = 7, dimnames = list(NULL, c("cvmin", "cv0.25", "cvmed", "cv0.75", "cvmax", 
        "cviqr", "time"))))
results.bandwidth <- cbind(parameters.bandwidth, matrix(nrow = length(bandwidth), 
    ncol = 7, dimnames = list(NULL, c("cvmin", "cv0.25", "cvmed", "cv0.75", "cvmax", 
        "cviqr", "time"))))

## Initialize loop for minimum peakwidth parameter sweeping
for (i in 1:length(peakwidth.min)) {
    ## Record time at start of loop
    start.time <- Sys.time()
    
    ## Preprocessing of raw files
    set <- xcmsSet(files, method = "centWave", peakwidth = c(peakwidth.min[i], preset[2]), 
        ppm = preset[3], snthresh = preset[4], mzdiff = preset[5], prefilter = c(preset[6], 
            preset[7]))
    
    set <- group.density(set, bw = preset[8])
    set.retcor <- retcor(set, method = "obiwarp", plottype = "none")
    set.retcor <- group.density(set.retcor, bw = preset[8])
    set.retcor.2 <- retcor(set.retcor, method = "obiwarp", plottype = "none")
    set.retcor.2 <- group.density(set.retcor.2, bw = preset[8])
    set.filled <- fillPeaks(set.retcor.2)
    
    ## Extraction of feature list
    peaklist <- peakTable(set.filled)
    
    ## Statistical calculations
    number.of.replicates <- length(files)
    peaklist.modified <- peaklist
    peaklist.modified$abundance.mean <- apply(peaklist.modified[, 9:(8 + number.of.replicates)], 
        1, FUN = mean, na.rm = T)
    peaklist.modified$abundance.sd <- apply(peaklist.modified[, 9:(8 + number.of.replicates)], 
        1, FUN = sd, na.rm = T)
    peaklist.modified$abundance.cv <- (peaklist.modified$abundance.sd * 100)/peaklist.modified$abundance.mean
    peaklist.modified <- arrange(peaklist.modified, rt, mz)
    
    ## Filter features by retention time (only features between 2 and 12 min) and
    ## abundance (only feature higher E6)
    peaklist.modified.filtered <- filter(peaklist.modified, rt > 60, rt < 720, peaklist.modified[, 
        8] == number.of.replicates)
    
    ## Summarising cv distribution in results table
    results.peakwidth.min[i, 9] <- round(min(peaklist.modified.filtered$abundance.cv), 
        0)
    results.peakwidth.min[i, 10] <- round(quantile(peaklist.modified.filtered$abundance.cv, 
        probs = 0.25), 0)
    results.peakwidth.min[i, 11] <- round(median(peaklist.modified.filtered$abundance.cv), 
        0)
    results.peakwidth.min[i, 12] <- round(quantile(peaklist.modified.filtered$abundance.cv, 
        probs = 0.75), 0)
    results.peakwidth.min[i, 13] <- round(max(peaklist.modified.filtered$abundance.cv), 
        0)
    results.peakwidth.min[i, 14] <- round(results.peakwidth.min[i, 11] - results.peakwidth.min[i, 
        9], 0)
    
    ## Recording time at end of loop and add difference to results table
    end.time <- Sys.time()
    results.peakwidth.min[i, 15] <- round(as.numeric(end.time) - as.numeric(start.time), 
        0)
}
results.peakwidth.min <- arrange(as.data.frame(results.peakwidth.min), cvmed, cviqr, 
    time)
write.csv(results.peakwidth.min, file = "centWaveOptResults/peakwidthMinResults.csv")

## Initialize loop for maximum peakwidth parameter sweeping
for (i in 1:length(peakwidth.max)) {
    ## Record time at start of loop
    start.time <- Sys.time()
    
    ## Preprocessing of raw files
    set <- xcmsSet(files, method = "centWave", peakwidth = c(preset[1], peakwidth.max[i]), 
        ppm = preset[3], snthresh = preset[4], mzdiff = preset[5], prefilter = c(preset[6], 
            preset[7]))
    
    set <- group.density(set, bw = preset[8])
    set.retcor <- retcor(set, method = "obiwarp", plottype = "none")
    set.retcor <- group.density(set.retcor, bw = preset[8])
    set.retcor.2 <- retcor(set.retcor, method = "obiwarp", plottype = "none")
    set.retcor.2 <- group.density(set.retcor.2, bw = preset[8])
    set.filled <- fillPeaks(set.retcor.2)
    
    ## Extraction of feature list and saving to hard disc
    peaklist <- peakTable(set.filled)
    
    ## Statistical calculations
    number.of.replicates <- length(files)
    peaklist.modified <- peaklist
    peaklist.modified$abundance.mean <- apply(peaklist.modified[, 9:(8 + number.of.replicates)], 
        1, FUN = mean, na.rm = T)
    peaklist.modified$abundance.sd <- apply(peaklist.modified[, 9:(8 + number.of.replicates)], 
        1, FUN = sd, na.rm = T)
    peaklist.modified$abundance.cv <- (peaklist.modified$abundance.sd * 100)/peaklist.modified$abundance.mean
    peaklist.modified <- arrange(peaklist.modified, rt, mz)
    
    ## Filter features by retention time (only features between 2 and 12 min) and
    ## abundance (only feature higher E6)
    peaklist.modified.filtered <- filter(peaklist.modified, rt > 60, rt < 720, peaklist.modified[, 
        8] == number.of.replicates)
    
    ## Summarising cv distribution in results table
    results.peakwidth.max[i, 9] <- round(min(peaklist.modified.filtered$abundance.cv), 
        0)
    results.peakwidth.max[i, 10] <- round(quantile(peaklist.modified.filtered$abundance.cv, 
        probs = 0.25), 0)
    results.peakwidth.max[i, 11] <- round(median(peaklist.modified.filtered$abundance.cv), 
        0)
    results.peakwidth.max[i, 12] <- round(quantile(peaklist.modified.filtered$abundance.cv, 
        probs = 0.75), 0)
    results.peakwidth.max[i, 13] <- round(max(peaklist.modified.filtered$abundance.cv), 
        0)
    results.peakwidth.max[i, 14] <- round(results.peakwidth.max[i, 11] - results.peakwidth.max[i, 
        9], 0)
    
    ## Recording time at end of loop and add difference to results table
    end.time <- Sys.time()
    results.peakwidth.max[i, 15] <- round(as.numeric(end.time) - as.numeric(start.time), 
        0)
}
results.peakwidth.max <- arrange(as.data.frame(results.peakwidth.max), cvmed, cviqr, 
    time)
write.csv(results.peakwidth.max, file = "centWaveOptResults/peakwidthMaxResults.csv")
save.image(file = "centWaveOpt.RData")

## Initialize loop for ppm parameter sweeping
for (i in 1:length(ppm)) {
    ## Record time at start of loop
    start.time <- Sys.time()
    
    ## Preprocessing of raw files
    set <- xcmsSet(files, method = "centWave", peakwidth = c(preset[1], preset[2]), 
        ppm = ppm[i], snthresh = preset[4], mzdiff = preset[5], prefilter = c(preset[6], 
            preset[7]))
    
    set <- group.density(set, bw = preset[8])
    set.retcor <- retcor(set, method = "obiwarp", plottype = "none")
    set.retcor <- group.density(set.retcor, bw = preset[8])
    set.retcor.2 <- retcor(set.retcor, method = "obiwarp", plottype = "none")
    set.retcor.2 <- group.density(set.retcor.2, bw = preset[8])
    set.filled <- fillPeaks(set.retcor.2)
    
    ## Extraction of feature list and saving to hard disc
    peaklist <- peakTable(set.filled)
    
    ## Statistical calculations
    number.of.replicates <- length(files)
    peaklist.modified <- peaklist
    peaklist.modified$abundance.mean <- apply(peaklist.modified[, 9:(8 + number.of.replicates)], 
        1, FUN = mean, na.rm = T)
    peaklist.modified$abundance.sd <- apply(peaklist.modified[, 9:(8 + number.of.replicates)], 
        1, FUN = sd, na.rm = T)
    peaklist.modified$abundance.cv <- (peaklist.modified$abundance.sd * 100)/peaklist.modified$abundance.mean
    peaklist.modified <- arrange(peaklist.modified, rt, mz)
    
    ## Filter features by retention time (only features between 2 and 12 min) and
    ## abundance (only feature higher E6)
    peaklist.modified.filtered <- filter(peaklist.modified, rt > 60, rt < 720, peaklist.modified[, 
        8] == number.of.replicates)
    
    ## Summarising cv distribution in results table
    results.ppm[i, 9] <- round(min(peaklist.modified.filtered$abundance.cv), 0)
    results.ppm[i, 10] <- round(quantile(peaklist.modified.filtered$abundance.cv, 
        probs = 0.25), 0)
    results.ppm[i, 11] <- round(median(peaklist.modified.filtered$abundance.cv), 
        0)
    results.ppm[i, 12] <- round(quantile(peaklist.modified.filtered$abundance.cv, 
        probs = 0.75), 0)
    results.ppm[i, 13] <- round(max(peaklist.modified.filtered$abundance.cv), 0)
    results.ppm[i, 14] <- round(results.ppm[i, 11] - results.ppm[i, 9], 0)
    
    ## Recording time at end of loop and add difference to results table
    end.time <- Sys.time()
    results.ppm[i, 15] <- round(as.numeric(end.time) - as.numeric(start.time), 0)
}
results.ppm <- arrange(as.data.frame(results.ppm), cvmed, cviqr, time)
write.csv(results.ppm, file = "centWaveOptResults/ppmResults.csv")
save.image(file = "centWaveOpt.RData")


## Initialize loop for snthresh parameter sweeping
for (i in 1:length(snthresh)) {
    ## Record time at start of loop
    start.time <- Sys.time()
    
    ## Preprocessing of raw files
    set <- xcmsSet(files, method = "centWave", peakwidth = c(preset[1], preset[2]), 
        ppm = preset[3], snthresh = snthresh[i], mzdiff = preset[5], prefilter = c(preset[6], 
            preset[7]))
    
    set <- group.density(set, bw = preset[8])
    set.retcor <- retcor(set, method = "obiwarp", plottype = "none")
    set.retcor <- group.density(set.retcor, bw = preset[8])
    set.retcor.2 <- retcor(set.retcor, method = "obiwarp", plottype = "none")
    set.retcor.2 <- group.density(set.retcor.2, bw = preset[8])
    set.filled <- fillPeaks(set.retcor.2)
    
    ## Extraction of feature list and saving to hard disc
    peaklist <- peakTable(set.filled)
    
    ## Statistical calculations
    number.of.replicates <- length(files)
    peaklist.modified <- peaklist
    peaklist.modified$abundance.mean <- apply(peaklist.modified[, 9:(8 + number.of.replicates)], 
        1, FUN = mean, na.rm = T)
    peaklist.modified$abundance.sd <- apply(peaklist.modified[, 9:(8 + number.of.replicates)], 
        1, FUN = sd, na.rm = T)
    peaklist.modified$abundance.cv <- (peaklist.modified$abundance.sd * 100)/peaklist.modified$abundance.mean
    peaklist.modified <- arrange(peaklist.modified, rt, mz)
    
    ## Filter features by retention time (only features between 2 and 12 min) and
    ## abundance (only feature higher E6)
    peaklist.modified.filtered <- filter(peaklist.modified, rt > 60, rt < 720, peaklist.modified[, 
        8] == number.of.replicates)
    
    ## Summarising cv distribution in results table
    results.snthresh[i, 9] <- round(min(peaklist.modified.filtered$abundance.cv), 
        0)
    results.snthresh[i, 10] <- round(quantile(peaklist.modified.filtered$abundance.cv, 
        probs = 0.25), 0)
    results.snthresh[i, 11] <- round(median(peaklist.modified.filtered$abundance.cv), 
        0)
    results.snthresh[i, 12] <- round(quantile(peaklist.modified.filtered$abundance.cv, 
        probs = 0.75), 0)
    results.snthresh[i, 13] <- round(max(peaklist.modified.filtered$abundance.cv), 
        0)
    results.snthresh[i, 14] <- round(results.snthresh[i, 11] - results.snthresh[i, 
        9], 0)
    
    ## Recording time at end of loop and add difference to results table
    end.time <- Sys.time()
    results.snthresh[i, 15] <- round(as.numeric(end.time) - as.numeric(start.time), 
        0)
}
results.snthresh <- arrange(as.data.frame(results.snthresh), cvmed, cviqr, time)
write.csv(results.snthresh, file = "centWaveOptResults/snthreshResults.csv")

## Initialize loop for mzdiff parameter sweeping
for (i in 1:length(mzdiff)) {
    ## Record time at start of loop
    start.time <- Sys.time()
    
    ## Preprocessing of raw files
    set <- xcmsSet(files, method = "centWave", peakwidth = c(preset[1], preset[2]), 
        ppm = preset[3], snthresh = preset[4], mzdiff = mzdiff[i], prefilter = c(preset[6], 
            preset[7]))
    
    set <- group.density(set, bw = preset[8])
    set.retcor <- retcor(set, method = "obiwarp", plottype = "none")
    set.retcor <- group.density(set.retcor, bw = preset[8])
    set.retcor.2 <- retcor(set.retcor, method = "obiwarp", plottype = "none")
    set.retcor.2 <- group.density(set.retcor.2, bw = preset[8])
    set.filled <- fillPeaks(set.retcor.2)
    
    ## Extraction of feature list and saving to hard disc
    peaklist <- peakTable(set.filled)
    
    ## Statistical calculations
    number.of.replicates <- length(files)
    peaklist.modified <- peaklist
    peaklist.modified$abundance.mean <- apply(peaklist.modified[, 9:(8 + number.of.replicates)], 
        1, FUN = mean, na.rm = T)
    peaklist.modified$abundance.sd <- apply(peaklist.modified[, 9:(8 + number.of.replicates)], 
        1, FUN = sd, na.rm = T)
    peaklist.modified$abundance.cv <- (peaklist.modified$abundance.sd * 100)/peaklist.modified$abundance.mean
    peaklist.modified <- arrange(peaklist.modified, rt, mz)
    
    ## Filter features by retention time (only features between 2 and 12 min) and
    ## abundance (only feature higher E6)
    peaklist.modified.filtered <- filter(peaklist.modified, rt > 60, rt < 720, peaklist.modified[, 
        8] == number.of.replicates)
    
    ## Summarising cv distribution in results table
    results.mzdiff[i, 9] <- round(min(peaklist.modified.filtered$abundance.cv), 0)
    results.mzdiff[i, 10] <- round(quantile(peaklist.modified.filtered$abundance.cv, 
        probs = 0.25), 0)
    results.mzdiff[i, 11] <- round(median(peaklist.modified.filtered$abundance.cv), 
        0)
    results.mzdiff[i, 12] <- round(quantile(peaklist.modified.filtered$abundance.cv, 
        probs = 0.75), 0)
    results.mzdiff[i, 13] <- round(max(peaklist.modified.filtered$abundance.cv), 
        0)
    results.mzdiff[i, 14] <- round(results.mzdiff[i, 11] - results.mzdiff[i, 9], 
        0)
    
    ## Recording time at end of loop and add difference to results table
    end.time <- Sys.time()
    results.mzdiff[i, 15] <- round(as.numeric(end.time) - as.numeric(start.time), 
        0)
}
results.mzdiff <- arrange(as.data.frame(results.mzdiff), cvmed, cviqr, time)
write.csv(results.mzdiff, file = "centWaveOptResults/mzdiffResults.csv")
save.image(file = "centWaveOpt.RData")


## Initialize loop for prefilter scan number parameter sweeping
for (i in 1:length(prefilter.scan.number)) {
    ## Record time at start of loop
    start.time <- Sys.time()
    
    ## Preprocessing of raw files
    set <- xcmsSet(files, method = "centWave", peakwidth = c(preset[1], preset[2]), 
        ppm = preset[3], snthresh = preset[4], mzdiff = preset[5], prefilter = c(prefilter.scan.number[i], 
            preset[7]))
    
    set <- group.density(set, bw = preset[8])
    set.retcor <- retcor(set, method = "obiwarp", plottype = "none")
    set.retcor <- group.density(set.retcor, bw = preset[8])
    set.retcor.2 <- retcor(set.retcor, method = "obiwarp", plottype = "none")
    set.retcor.2 <- group.density(set.retcor.2, bw = preset[8])
    set.filled <- fillPeaks(set.retcor.2)
    
    ## Extraction of feature list and saving to hard disc
    peaklist <- peakTable(set.filled)
    
    ## Statistical calculations
    number.of.replicates <- length(files)
    peaklist.modified <- peaklist
    peaklist.modified$abundance.mean <- apply(peaklist.modified[, 9:(8 + number.of.replicates)], 
        1, FUN = mean, na.rm = T)
    peaklist.modified$abundance.sd <- apply(peaklist.modified[, 9:(8 + number.of.replicates)], 
        1, FUN = sd, na.rm = T)
    peaklist.modified$abundance.cv <- (peaklist.modified$abundance.sd * 100)/peaklist.modified$abundance.mean
    peaklist.modified <- arrange(peaklist.modified, rt, mz)
    
    ## Filter features by retention time (only features between 2 and 12 min) and
    ## abundance (only feature higher E6)
    peaklist.modified.filtered <- filter(peaklist.modified, rt > 60, rt < 720, peaklist.modified[, 
        8] == number.of.replicates)
    
    ## Summarising cv distribution in results table
    results.prefilter.scan.number[i, 9] <- round(min(peaklist.modified.filtered$abundance.cv), 
        0)
    results.prefilter.scan.number[i, 10] <- round(quantile(peaklist.modified.filtered$abundance.cv, 
        probs = 0.25), 0)
    results.prefilter.scan.number[i, 11] <- round(median(peaklist.modified.filtered$abundance.cv), 
        0)
    results.prefilter.scan.number[i, 12] <- round(quantile(peaklist.modified.filtered$abundance.cv, 
        probs = 0.75), 0)
    results.prefilter.scan.number[i, 13] <- round(max(peaklist.modified.filtered$abundance.cv), 
        0)
    results.prefilter.scan.number[i, 14] <- round(results.prefilter.scan.number[i, 
        11] - results.prefilter.scan.number[i, 9], 0)
    
    ## Recording time at end of loop and add difference to results table
    end.time <- Sys.time()
    results.prefilter.scan.number[i, 15] <- round(as.numeric(end.time) - as.numeric(start.time), 
        0)
}
results.prefilter.scan.number <- arrange(as.data.frame(results.prefilter.scan.number), 
    cvmed, cviqr, time)
write.csv(results.prefilter.scan.number, file = "centWaveOptResults/prefilterScanNumberResults.csv")
save.image(file = "centWaveOpt.RData")


## Initialize loop for prefilter scan abundance parameter sweeping
for (i in 1:length(prefilter.scan.abundance)) {
    ## Record time at start of loop
    start.time <- Sys.time()
    
    ## Preprocessing of raw files
    set <- xcmsSet(files, method = "centWave", peakwidth = c(preset[1], preset[2]), 
        ppm = preset[3], snthresh = preset[4], mzdiff = preset[5], prefilter = c(preset[6], 
            prefilter.scan.abundance[i]))
    
    set <- group.density(set, bw = preset[8])
    set.retcor <- retcor(set, method = "obiwarp", plottype = "none")
    set.retcor <- group.density(set.retcor, bw = preset[8])
    set.retcor.2 <- retcor(set.retcor, method = "obiwarp", plottype = "none")
    set.retcor.2 <- group.density(set.retcor.2, bw = preset[8])
    set.filled <- fillPeaks(set.retcor.2)
    
    ## Extraction of feature list and saving to hard disc
    peaklist <- peakTable(set.filled)
    
    ## Statistical calculations
    number.of.replicates <- length(files)
    peaklist.modified <- peaklist
    peaklist.modified$abundance.mean <- apply(peaklist.modified[, 9:(8 + number.of.replicates)], 
        1, FUN = mean, na.rm = T)
    peaklist.modified$abundance.sd <- apply(peaklist.modified[, 9:(8 + number.of.replicates)], 
        1, FUN = sd, na.rm = T)
    peaklist.modified$abundance.cv <- (peaklist.modified$abundance.sd * 100)/peaklist.modified$abundance.mean
    peaklist.modified <- arrange(peaklist.modified, rt, mz)
    
    ## Filter features by retention time (only features between 2 and 12 min) and
    ## abundance (only feature higher E6)
    peaklist.modified.filtered <- filter(peaklist.modified, rt > 60, rt < 720, peaklist.modified[, 
        8] == number.of.replicates)
    
    ## Summarising cv distribution in results table
    results.prefilter.scan.abundance[i, 9] <- round(min(peaklist.modified.filtered$abundance.cv), 
        0)
    results.prefilter.scan.abundance[i, 10] <- round(quantile(peaklist.modified.filtered$abundance.cv, 
        probs = 0.25), 0)
    results.prefilter.scan.abundance[i, 11] <- round(median(peaklist.modified.filtered$abundance.cv), 
        0)
    results.prefilter.scan.abundance[i, 12] <- round(quantile(peaklist.modified.filtered$abundance.cv, 
        probs = 0.75), 0)
    results.prefilter.scan.abundance[i, 13] <- round(max(peaklist.modified.filtered$abundance.cv), 
        0)
    results.prefilter.scan.abundance[i, 14] <- round(results.prefilter.scan.abundance[i, 
        11] - results.prefilter.scan.abundance[i, 9], 0)
    
    ## Recording time at end of loop and add difference to results table
    end.time <- Sys.time()
    results.prefilter.scan.abundance[i, 15] <- round(as.numeric(end.time) - as.numeric(start.time), 
        0)
}
results.prefilter.scan.abundance <- arrange(as.data.frame(results.prefilter.scan.abundance), 
    cvmed, cviqr, time)
write.csv(results.prefilter.scan.abundance, file = "centWaveOptResults/prefilterScanAbundanceResults.csv")
save.image(file = "centWaveOpt.RData")


## Initialize loop for bandwidth parameter sweeping
for (i in 1:length(bandwidth)) {
    ## Record time at start of loop
    start.time <- Sys.time()
    
    ## Preprocessing of raw files
    set <- xcmsSet(files, method = "centWave", peakwidth = c(preset[1], preset[2]), 
        ppm = preset[3], snthresh = preset[4], mzdiff = preset[5], prefilter = c(preset[6], 
            prefilter.scan.abundance[7]))
    
    set <- group.density(set, bw = bandwidth[i])
    set.retcor <- retcor(set, method = "obiwarp", plottype = "none")
    set.retcor <- group.density(set.retcor, bw = bandwidth[i])
    set.retcor.2 <- retcor(set.retcor, method = "obiwarp", plottype = "none")
    set.retcor.2 <- group.density(set.retcor.2, bw = bandwidth[i])
    set.filled <- fillPeaks(set.retcor.2)
    
    ## Extraction of feature list and saving to hard disc
    peaklist <- peakTable(set.filled)
    
    ## Statistical calculations
    number.of.replicates <- length(files)
    peaklist.modified <- peaklist
    peaklist.modified$abundance.mean <- apply(peaklist.modified[, 9:(8 + number.of.replicates)], 
        1, FUN = mean, na.rm = T)
    peaklist.modified$abundance.sd <- apply(peaklist.modified[, 9:(8 + number.of.replicates)], 
        1, FUN = sd, na.rm = T)
    peaklist.modified$abundance.cv <- (peaklist.modified$abundance.sd * 100)/peaklist.modified$abundance.mean
    peaklist.modified <- arrange(peaklist.modified, rt, mz)
    
    ## Filter features by retention time (only features between 2 and 12 min) and
    ## abundance (only feature higher E6)
    peaklist.modified.filtered <- filter(peaklist.modified, rt > 60, rt < 720, peaklist.modified[, 
        8] == number.of.replicates)
    
    ## Summarising cv distribution in results table
    results.bandwidth[i, 9] <- round(min(peaklist.modified.filtered$abundance.cv), 
        0)
    results.bandwidth[i, 10] <- round(quantile(peaklist.modified.filtered$abundance.cv, 
        probs = 0.25), 0)
    results.bandwidth[i, 11] <- round(median(peaklist.modified.filtered$abundance.cv), 
        0)
    results.bandwidth[i, 12] <- round(quantile(peaklist.modified.filtered$abundance.cv, 
        probs = 0.75), 0)
    results.bandwidth[i, 13] <- round(max(peaklist.modified.filtered$abundance.cv), 
        0)
    results.bandwidth[i, 14] <- round(results.bandwidth[i, 11] - results.bandwidth[i, 
        9], 0)
    
    ## Recording time at end of loop and add difference to results table
    end.time <- Sys.time()
    results.bandwidth[i, 15] <- round(as.numeric(end.time) - as.numeric(start.time), 
        0)
}
results.bandwidth <- arrange(as.data.frame(results.bandwidth), cvmed, cviqr, time)
write.csv(results.bandwidth, file = "centWaveOptResults/bandwidthResults.csv")
save.image(file = "centWaveOpt.RData")


## Evaluate results tables and summarise parameter values
parameter <- NULL
results.peakwidth.min.cropped <- results.peakwidth.min[results.peakwidth.min$cvmed != min(results.peakwidth.min$cvmed),]
parameter[1] <- max(results.peakwidth.min.cropped$Minimum.Peakwidth[results.peakwidth.min.cropped$cvmed == min(results.peakwidth.min.cropped$cvmed)])

results.peakwidth.max.cropped <- results.peakwidth.max[results.peakwidth.max$cvmed != min(results.peakwidth.max$cvmed),]
parameter[2] <- max(results.peakwidth.max.cropped$Maximum.Peakwidth[results.peakwidth.max.cropped$cvmed == min(results.peakwidth.max.cropped$cvmed)])

results.ppm.cropped <- results.ppm[results.ppm$cvmed != min(results.ppm$cvmed),]
parameter[3] <- max(results.ppm.cropped$ppm[results.ppm.cropped$cvmed == min(results.ppm.cropped$cvmed)])

results.snthresh.cropped <- results.snthresh[results.snthresh$cvmed != min(results.snthresh$cvmed),]
parameter[4] <- max(results.snthresh.cropped$S.N.Threshold[results.snthresh.cropped$cvmed == min(results.snthresh.cropped$cvmed)])

results.mzdiff.cropped <- results.mzdiff[results.mzdiff$cvmed != min(results.mzdiff$cvmed),]
parameter[5] <- max(results.mzdiff.cropped$m.z.Diff[results.mzdiff.cropped$cvmed == min(results.mzdiff.cropped$cvmed)])

results.prefilter.scan.number.cropped <- results.prefilter.scan.number[results.prefilter.scan.number$cvmed != min(results.prefilter.scan.number$cvmed),]
parameter[6] <- max(results.prefilter.scan.number.cropped$Prefilter.1[results.prefilter.scan.number.cropped$cvmed == min(results.prefilter.scan.number.cropped$cvmed)])

results.prefilter.scan.abundance.cropped <- results.prefilter.scan.abundance[results.prefilter.scan.abundance$cvmed != min(results.prefilter.scan.abundance$cvmed),]
parameter[7] <- max(results.prefilter.scan.abundance.cropped$Prefilter.2[results.prefilter.scan.abundance.cropped$cvmed == min(results.prefilter.scan.abundance.cropped$cvmed)])

results.bandwidth.cropped <- results.bandwidth[results.bandwidth$cvmed != min(results.bandwidth$cvmed),]
parameter[8] <- max(results.bandwidth.cropped$bandwidth[results.bandwidth.cropped$cvmed == min(results.bandwidth.cropped$cvmed)])

## Export sorted parameter file
write.csv(t(parameter), "centWaveOptResults/parameter.csv", row.names = FALSE)

## Save data set for hard drive
save.image(file = "centWaveOpt.RData")
