#load settings file into data frame file_data
file = "settings_local.txt"
file_data <- read.delim(file, header = FALSE, sep = "\t", dec = ".")

#load necessary packages
library(foreach)
library(doParallel)
library(Cardinal)

#VARIABLES

#Consider maintaining a global variable to increment some file suffix which indicates
#which process run or variable change occured for multi-run Slurm projects

#name of imzML file to read in
name = as.character(file_data[1, "V2"])

#preprocessing and processing on decisions
#check if this coercion is valid
normalize_on = as.logical(file_data[2, "V2"])
pick_on = as.logical(file_data[3, "V2"])
align_on = as.logical(file_data[4, "V2"])
filter_on = as.logical(file_data[5, "V2"])

SNR = as.integer(file_data[6, "V2"])

#Option to save separate objects after each (pre)processing method
#Should make way to save these to their own directory
save_intermediates = as.logical(file_data[7, "V2"])

#directory created for output
write_dir = as.character(file_data[8, "V2"])

#turn off queueing for (pre)processing functions, make new objects with each operation
options(Cardinal.delay=FALSE)

#Import data
data_mse <- readImzML(name,resolution = 5, units = "ppm", mass.range = c(70,600),as="MSImagingExperiment",attach.only = TRUE)

#coerce into MSContinuousImagingExperiment due to documented and apparently uncorrected Cardinal error
data <- as(data_mse, "MSContinuousImagingExperiment")

#extra copy
data_backup <- data

dir.create(write_dir)

#initial plot of data
png(paste(write_dir, name, "_first_100_pixels.png", sep = ""))
plot(data,pixel=1:100,main="spectra for first 100 pixels")
dev.off()

if (save_intermediates == TRUE) {
  save(data, file=paste(write_dir, name, "_initial_data.RData", sep =""))
}

#conditionally normalized data
if (normalize_on == TRUE) {
data_n <- normalize(data, method="rms")
rm(data)
data <- data_n
if (save_intermediates == TRUE) {
  save(data_n, file=paste(write_dir, name, "_normalized_data.RData", sep =""))
}
rm(data_n)
}

#conditionally picked data, must not occur before normalizing
if (pick_on == TRUE) {
data_p <- peakPick(data, method="adaptive", SNR)
rm(data)
data <- data_p
if (save_intermediates == TRUE) {
  save(data_p, file=paste(write_dir, name, "_picked_data.RData", sep =""))
}
rm(data_p)
}

#conditionally aligned data, must not occur before preprocessing
if (align_on == TRUE) {
data_a <- peakAlign(data, method="diff",diff.max=5,units="ppm")
rm(data)
data <- data_a
if (save_intermediates == TRUE) {
  save(data_a, file=paste(write_dir, name, "_aligned_data.RData", sep =""))
}
rm(data_a)
}

#plot processed pixels in image, plot, and spectral variations
source("utils2.R")
png(paste(write_dir, name, "_data_norm_pick_align_image.png", sep = ""))
.setup.layout(c(1,1),bg="black")
image(data,mz=283.26,plusminus=5*283.26/1e6,main="parallel")
dev.off()

source("utils2.R")
png(paste(write_dir, name, "_data_norm_pick_align_plot.png", sep = ""))
.setup.layout(c(1,1),bg="white")
plot(data,pixel=1:ncol(data),main="parallel")
dev.off()

png(paste(write_dir, name, "_data_norm_pick_align_spectra.png", sep = ""))
plot(data,pixel=seq(1,ncol(data),10),main="peak aligned_spectra")
dev.off()

#select peaks w/ frequency greater than 0.5, removes approx. 2000 peaks, maintains 1600
if (filter_on ==TRUE){
data_f <- peakFilter(data, freq.min= 0.5)
rm(data)
data <- data_f
if (save_intermediates == TRUE) {
  save(data_f, file=paste(write_dir, name, "_filtered_data.RData", sep =""))
}
rm(data_f)
}

#image and plot of picked peaks
source("utils2.R")
png(paste(write_dir, name, "_data_norm_pick_align_filter_plot.png", sep = ""))
.setup.layout(c(1,1),bg="white")
plot(data,pixel=1:ncol(data),main="parallel")
dev.off()

source("utils2.R")
png(paste(write_dir, name, "_data_norm_pick_align_filter_image.png", sep = ""))
.setup.layout(c(1,1),bg="black")
image(data, mz=283.26, plusminus=5*283.26/1e6, main="parallel")
dev.off()

#cleanup
 gc()
 ls()