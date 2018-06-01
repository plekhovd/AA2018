library(raster)
library(rgdal)
library(dplyr)
library(ggpubr)
library(RStoolbox)

wd = "/Users/danplekhov/Desktop/SAA2018/GIS/Phillipston"
setwd(wd)

projection = "+init=EPSG:32618"   # set projection (UTM 18N)
projection = "+init=EPSG:32619"   # set projection (UTM 19N)
crs = crs(projection) # create crs objection

#process image tiles
#image.dir = dir("image", full.names = T)
#image.raster = lapply(image.dir, stack)
#image.merge = do.call("merge", image.raster)
#image = projectRaster(image.merge, crs=crs)

#process LIDAR tiles
lidar.dir = dir("lidar", pattern = ".img$", full.names = T)
lidar.raster = lapply(lidar.dir, raster)
lidar.merge = do.call("merge",lidar.raster)
#crs(lidar.merge) = crs("+proj=lcc +lat_1=41.2 +lat_2=41.86666666666667 +lat_0=40.83333333333334 +lon_0=-72.75 +x_0=304800.6096 +y_0=152400.3048 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs")
#lidar.merge = projectRaster(lidar.merge, crs = crs)
writeRaster(lidar.merge, "lidar/Purgatory_lidar.tif", overwrite = T)

e = extent(lidar.merge)
e = as(e, "SpatialPolygons")
proj4string(e) = CRS(projection)
shapefile(e, "Pottapaug_extent", overwrite=T)
e.planet = spTransform(e, CRS("+init=EPSG:4326"))
shapefile(e.planet, "Purgatory_extent_planet", overwrite=T)

#process Planet data
p.dir = dir("planet_order_186916", full.names = T, pattern = "2017")
names = unique(substr(basename(p.dir), 1, 8))
for (i in 1:length(names)){
  p.dir2 = p.dir[grepl(names[i], p.dir)]
  p.raster = list()
  
  for (j in 1:length(p.dir2)){
    p.raster[j] = stack(dir(p.dir2[j], pattern = "MS_SR.tif$", full.names = T))#/10000
    if (!compareCRS(p.raster[[j]],crs)){
      p.raster[[j]] = projectRaster(p.raster[[j]], crs=crs, res = c(3,3))
    }
  }
  
  if (length(p.raster) > 1){
    p.merge = do.call("merge", p.raster)
    p.merge = crop(p.merge, e)
    p.raster = p.merge
  }  else{
    p.raster = stack(crop(p.raster[[1]], e))
  }
  
  writeRaster(p.raster, paste0("merged/",names[i],".tif"), overwrite = T)
}

#NDVI
ndvi.dir = dir("merged", full.names = T, pattern = ".tif$")
ndvi = list()
evi = list()
for (i in 1:length(ndvi.dir)){
  image = stack(ndvi.dir[i])
  names(image) = paste0("B", 1:length(image@layers))
  ndvi.image = (image$B4-image$B3)/(image$B4+image$B3)
  ndvi[i] = ndvi.image
  evi.image = (2.5*(image$B4-image$B3))/((image$B2+(6*image$B3)-(7.5*image$B1)+1))
  evi[i] = evi.image
}

for (i in 1:length(ndvi)){
  if (ndvi[[i]]@ncols != ndvi[[1]]@ncols){
    ndvi[[i]] = projectRaster(ndvi[[i]], ndvi[[1]])
    }
  }

ndvi = stack(ndvi)
writeRaster(ndvi, "ndvi_stack.tif", overwrite = T)
evi = stack(evi)
writeRaster(evi, "evi_stack.tif", overwrite = T)

#Analysis
dates = unique(substr(basename(dir("merged/")), 1, 8))
dates = paste0(substr(dates, 5,6), "/", substr(dates, 7,8))

lidar.merge = raster("/Users/danplekhov/Desktop/SAA2018/GIS/Phillipston/lidar/Phillipston_lidar.tif")
e = extent(lidar.merge)
e = as(e, "SpatialPolygons")
proj4string(e) = CRS(projection)

ndvi = stack("ndvi_stack.tif")

boundary = readOGR("polygon/polygon.shp")
boundary.mask = rasterize(boundary, ndvi[[1]])
boundary.mask = reclassify(boundary.mask, c(-Inf, 0,NA))
ndvi.in = mask(ndvi, mask = boundary.mask, inverse = F)
ndvi.out = mask(ndvi, mask = boundary.mask, inverse = T)

sum(!is.na(getValues(ndvi.in[[1]])))

iterations = 50
n = 100
seq.i = seq(from = 1, to = iterations*n, by = n)

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

ggboxplot(spectra, x = "Date", y = "NDVI", color = "Location", 
          palette = c("pink", "light blue"))

ggline(spectra, x = "Date", y = "NDVI", color = "Location", palette = c("pink", "light blue"), plot_type = "l",
                    add = "mean_ci")


summary(lm(spectra$NDVI~spectra$Date*spectra$Location))
anova(lm(spectra$NDVI~spectra$Date*spectra$Location))

#check for proper statistic test
rquery.t.test(in.spectra[,1], out.spectra[,1])

#check normality
shapiro.test(spectra$NDVI[spectra$Location=="Non-Enclosed" & spectra$Date ==dates[1]])
shapiro.test(spectra$logNDVI[spectra$Location=="Non-Enclosed" & spectra$Date ==dates[1]])
shapiro.test(spectra$sqrtNDVI[spectra$Location=="Non-Enclosed" & spectra$Date ==dates[1]])

#Wilcoxon Test
wilcox.list = list()
for (i in 1:ncol(in.spectra)){
  m = spectra[spectra$Date==dates[i],]
  group_by(m, Location) %>%
    summarise(
      count = n(),
      median = median(NDVI, na.rm = TRUE),
      IQR = IQR(NDVI, na.rm = TRUE)
    )
  ggboxplot(m, x = "Location", y = "NDVI", 
            color = "Location", palette = c("#00AFBB", "#E7B800"),
            ylab = "NDVI", xlab = "Location")
  res = wilcox.test(m$NDVI[m$Location=="Enclosed"], m$NDVI[m$Location=="Non-Enclosed"])
  wilcox.list[[i]] = res
}

dunn.test::dunn.test(m$NDVI[m$Location=="Enclosed"], m$NDVI[m$Location=="Non-Enclosed"], method = "bonferroni")

#Two Way ANOVA
res.aov2 = aov(NDVI ~ Location + Date, data = spectra)
summary(res.aov2)
model.tables(res.aov2, type="means", se = TRUE)
TukeyHSD(res.aov2, which = "Date") 
plot(res.aov2, 1) #Homogeneity of Variance
leveneTest(NDVI ~ Location * Date, data = spectra) #Homogeneity of Variance
plot(res.aov2, 2) #Normality
