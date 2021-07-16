#load necessary packages
library(foreach)
library(doParallel)
library(Cardinal)

#VARIABLES

name = "S1-1"

#preprocessing and processing on decisions
normalize_on = TRUE
pick_on = TRUE
align_on = TRUE
filter_on = TRUE

#Option to save separate objects after each (pre)processing method
save_intermediates = TRUE

#turn off queueing for (pre)processing functions, make new objects with each operation
options(Cardinal.delay=FALSE)

#Import data
data_mse <- readImzML(name,resolution = 5, units = "ppm", mass.range = c(70,600),as="MSImagingExperiment",attach.only = TRUE)

#select two pixels on entire mass range, reduces file size for non-Slurm analysis
#COMMENT OUT BELOW LINE FOR FULL DATASET
data_mse <- data_mse[1:214844, 5000:5001]

#coerce into MSContinuousImagingExperiment due to documented and apparently uncorrected Cardinal error
data <- as(data_mse, "MSContinuousImagingExperiment")
#extra copy
data_backup <- data

#initial plot of data
png(paste(name, "first_2_pixels.png", sep = "_"))
plot(data,pixel=1:2,main="spectra for first 2 pixels")
dev.off()

if (save_intermediates == TRUE) {
  save(data, file=paste(name, "initial_data.RData", sep ="_"))
}

#conditionally normalized data
if (normalize_on == TRUE) {
data_n <- normalize(data, method="rms")
rm(data)
data <- data_n
if (save_intermediates == TRUE) {
  save(data_n, file=paste(name, "normalized_data.RData", sep ="_"))
}
rm(data_n)
}

#conditionally picked data, must not occur before normalizing
if (pick_on == TRUE) {
data_p <- peakPick(data, method="adaptive", SNR=6)
rm(data)
data <- data_p
if (save_intermediates == TRUE) {
  save(data_p, file=paste(name, "picked_data.RData", sep ="_"))
}
rm(data_p)
}

#conditionally aligned data, must not occur before preprocessing
if (align_on == TRUE) {
data_a <- peakAlign(data, method="diff",diff.max=5,units="ppm")
rm(data)
data <- data_a
if (save_intermediates == TRUE) {
  save(data_a, file=paste(name, "aligned_data.RData", sep ="_"))
}
rm(data_a)
}

#plot processed pixels in image, plot, and spectral variations
source("/Genomics/grid/users/nrpark/cardinalelucidata/Cardinal/R/utils2.R")
png(paste(name, "data_norm_pick_align_image.png", sep = "_"))
.setup.layout(c(1,1),bg="black")
image(data,mz=283.26,plusminus=5*283.26/1e6,main="parallel")
dev.off()

source("/Genomics/grid/users/nrpark/cardinalelucidata/Cardinal/R/utils2.R")
png(paste(name, "data_norm_pick_align_plot.png", sep = "_"))
.setup.layout(c(1,1),bg="white")
plot(data,pixel=1:ncol(data),main="parallel")
dev.off()

png(paste(name, "data_norm_pick_align_spectra.png", sep = "_"))
plot(data,pixel=seq(1,ncol(data),10),main="peak aligned_spectra")
dev.off()

#select peaks w/ frequency greater than 0.5, removes approx. 2000 peaks, maintains 1600
if (filter_on ==TRUE){
data_f <- peakFilter(data, freq.min= 0.5)
rm(data)
data <- data_f
if (save_intermediates == TRUE) {
  save(data_f, file=paste(name, "filtered_data.RData", sep ="_"))
}
rm(data_f)
}

#image and plot of picked peaks
source("/Genomics/grid/users/nrpark/cardinalelucidata/Cardinal/R/utils2.R")
png(paste(name, "data_norm_pick_align_filter_plot.png", sep = "_"))
.setup.layout(c(1,1),bg="white")
plot(data_norm_pick_align_filter,pixel=1:ncol(data_norm_pick_align_filter),main="parallel")
dev.off()

source("/Genomics/grid/users/nrpark/cardinalelucidata/Cardinal/R/utils2.R")
png(paste(name, "data_norm_pick_align_filter_image.png", sep = "_"))
.setup.layout(c(1,1),bg="black")
image(data_norm_pick_align_filter,mz=283.26,plusminus=5*283.26/1e6,main="parallel")
dev.off()

#cleanup
 gc()
 ls()