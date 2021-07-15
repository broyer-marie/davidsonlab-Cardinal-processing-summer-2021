library(foreach)
library(doParallel)
library(Cardinal)

options(Cardinal.delay=FALSE)

data_mse <- readImzML(name="S1-1",resolution = 5, units = "ppm",mass.range = c(70,600),as="MSImagingExperiment",attach.only = TRUE)

data <- as(data_mse, "MSContinuousImagingExperiment")

png("S1-1_first_100_pixels.png")
plot(data,pixel=1:100,main="spectra for first 100 pixels")
dev.off()

data_norm <- normalize(data, method="rms")

data_norm_pick <- peakPick(data_norm, method="adaptive", SNR=6)

data_norm_pick_align <- peakAlign(data_norm_pick, method="diff",diff.max=5,units="ppm")

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

data_norm_pick_align_filter <- peakFilter(data_norm_pick_align, method="freq",freq.min=0.01)

source("/Genomics/grid/users/nrpark/cardinalelucidata/Cardinal/R/utils2.R")

df <- read.csv("/Genomics/grid/users/nrpark/P1_neg_ionization.csv")
df <- data.frame(df$Name,df$mz_neg)

dir.create("/Genomics/davidsonlab/Brianna/Cardinal_Test_Analysis_1")

for(i in 1:nrow(df)){
  print(df$df.Name[i])
  png(paste("/Genomics/davidsonlab/Brianna/Cardinal_Test_Analysis_1/",df$df.Name[i],".png",sep = ""), height = 750, width = 750)
  .setup.layout(c(1,1),bg = "black")
  col=bw.colors
  image(subset,mz = df$df.mz_neg[i], plusminus = 5*df$df.mz_neg[i]/1e6,contrast.enhance="histogram", col = bw.colors)
  dev.off()
}