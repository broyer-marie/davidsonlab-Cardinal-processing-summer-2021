#load necessary packages
library(foreach)
library(doParallel)
library(Cardinal)

#turn off queueing for (pre)processing functions, make new objects with each operation
options(Cardinal.delay=FALSE)

#Import data
data_mse <- readImzML(name="S1-1",resolution = 5, units = "ppm",mass.range = c(70,600),as="MSImagingExperiment",attach.only = TRUE)

#select two pixels on entire mass range, reduces file size for non-Slurm analysis
data_reduced <- data_mse[1:214844, 5000:5001]

#coerce into MSContinuousImagingExperiment due to documented and apparently uncorrected Cardinal error
data <- as(data_reduced, "MSContinuousImagingExperiment")

#initial plot of data
png("S1-1_first_2_pixels.png")
plot(data,pixel=1:2,main="spectra for first 2 pixels")
dev.off()

#normalized data
data_norm <- normalize(data, method="rms")

#peak picked data serially after normalizing
data_norm_pick <- peakPick(data_norm, method="adaptive", SNR=6)

#aligned data serially after peak picking
data_norm_pick_align <- peakAlign(data_norm_pick, method="diff",diff.max=5,units="ppm")

#plot processed pixels in image, plot, and spectral variations
source("/Genomics/grid/users/nrpark/cardinalelucidata/Cardinal/R/utils2.R")
png("S1-1_data_norm_pick_align_image.png")
.setup.layout(c(1,1),bg="black")
image(data_norm_pick_align,mz=283.26,plusminus=5*283.26/1e6,main="parallel")
dev.off()

source("/Genomics/grid/users/nrpark/cardinalelucidata/Cardinal/R/utils2.R")
png("S1-1_data_norm_pick_align_plot.png")
.setup.layout(c(1,1),bg="white")
plot(data_norm_pick_align,pixel=1:ncol(data_norm_pick_align),main="parallel")
dev.off()

png("S1-1_data_norm_pick_align_spectra.png")
plot(data_norm_pick_align,pixel=seq(1,ncol(data_norm_pick_align),10),main="peak aligned_spectra")
dev.off()

#select peaks w/ frequency greater than 0.5, removes approx. 2000 peaks, maintains 1600
data_norm_pick_align_filter <- peakFilter(data_norm_pick_align, freq.min= 0.5)

#image and plot of picked peaks
source("/Genomics/grid/users/nrpark/cardinalelucidata/Cardinal/R/utils2.R")
png("S1-1_data_norm_pick_align_filter_plot.png")
.setup.layout(c(1,1),bg="white")
plot(data_norm_pick_align_filter,pixel=1:ncol(data_norm_pick_align_filter),main="parallel")
dev.off()

source("/Genomics/grid/users/nrpark/cardinalelucidata/Cardinal/R/utils2.R")
png("S1-1_data_norm_pick_align_filter_image.png")
.setup.layout(c(1,1),bg="black")
image(data_norm_pick_align_filter,mz=283.26,plusminus=5*283.26/1e6,main="parallel")
dev.off()

#cleanup
 gc()
 ls()