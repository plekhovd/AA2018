library(raster)
library(rgdal)
library(dplyr)
library(ggpubr)
library(RStoolbox)


projection = "+init=EPSG:32618"   # set projection (UTM 18N)
projection = "+init=EPSG:32619"   # set projection (UTM 19N)
crs = crs(projection) # create crs objection

#import LiDAR data
lidar = raster("")

e = extent(lidar)
e = as(e, "SpatialPolygons")
proj4string(e) = CRS(projection)

#process Planet data
#stack time series Planet data. Must have same projection and extent.
planet.data = "/Users/danplekhov/Desktop/SAA2018/GIS/Gay City/merged"
planet.dir = dir(planet.data, pattern = ".tif$", full.names = T)
planet = lapply(planet.dir, stack)

#NDVI
ndvi = list()
for (i in 1:length(planet)){
  image = stack(planet[[i]])
  names(image) = paste0("B", 1:length(image@layers))
  ndvi.image = (image$B4-image$B3)/(image$B4+image$B3)
  ndvi[i] = ndvi.image
}

for (i in 1:length(ndvi)){
  if (ndvi[[i]]@ncols != ndvi[[1]]@ncols){
    ndvi[[i]] = projectRaster(ndvi[[i]], ndvi[[1]])
    }
  }

ndvi = stack(ndvi)

###Analysis
dates = unique(substr(basename(dir(planet.data)), 1, 8))
dates = paste0(substr(dates, 5,6), "/", substr(dates, 7,8))

#import stone wall polygons and use to mask NDVI images
boundary = readOGR("polygon.shp")
boundary.mask = rasterize(boundary, ndvi[[1]])
boundary.mask = reclassify(boundary.mask, c(-Inf, 0,NA))
ndvi.in = mask(ndvi, mask = boundary.mask, inverse = F)
ndvi.out = mask(ndvi, mask = boundary.mask, inverse = T)

#count how many cells are within polygons
sum(!is.na(getValues(ndvi.in[[1]])))

#set sampling parameters
iterations = 50
n = 100
seq.i = seq(from = 1, to = iterations*n, by = n)

#create matrices to store sampled values for "in" and "out" spectra
in.spectra = matrix(NA, nrow = n*iterations, ncol = length(ndvi@layers))
out.spectra = matrix(NA, nrow = n*iterations, ncol = length(ndvi@layers))

for (i in 1:iterations){
  in.values = sampleRandom(ndvi.in, size = n)
  in.spectra[seq.i[i]:(seq.i[i]+(n-1)),] = in.values
  
  out.values = sampleRandom(ndvi.out, size = n)
  out.spectra[seq.i[i]:(seq.i[i]+(n-1)),] = out.values
}

spectra = matrix(c(in.spectra,out.spectra), nrow = (n*iterations*2*length(ndvi@layers)), ncol = 3, byrow = FALSE)
spectra[,2] = rep(1:length(ndvi@layers), time=2, each=n*iterations)
spectra[,3] = rep(1:2, time=1, each=n*iterations*length(ndvi@layers))
spectra = as.data.frame(spectra)
colnames(spectra)  = c("NDVI", "Date", "Location")
spectra$Date = factor(spectra$Date, 
                      levels = seq(1,length(ndvi@layers),1),
                      labels = dates)
spectra$Date = as.Date(spectra$Date, "%m/%d")

spectra$Location = factor(spectra$Location, 
                          levels = c(1, 2),
                          labels = c("Enclosed", "Non-Enclosed"))

table(spectra$Date, spectra$Location)

#Box plot figure
ggboxplot(spectra, x = "Date", y = "NDVI", fill = "Location", palette = c("grey","white"))

#Line plot figure
ggline(spectra, x = "Date", y = "NDVI", color = "Location", palette = c("black", "grey"), plot_type = "l",
                    add = "mean_ci")


#ANCOVA Statistic
summary(lm(spectra$NDVI~spectra$Date*spectra$Location))
anova(lm(spectra$NDVI~spectra$Date*spectra$Location))


