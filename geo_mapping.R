#######################################
###      create sampling frame      ###
#######################################

library(dplyr)
library(lubridate)
library(rgdal)  # read GeoJSON data sets
library(rgeos)  # spatial manipulation

#######################################
###          get data sets          ###
#######################################

# ward data set
url.ward <- "http://opendata.dc.gov/datasets/0ef47379cbae44e88267c01eaec2ff6e_31.geojson"
ward.spdf <- rgdal::readOGR(url.ward, layer = "OGRGeoJSON")

# station entrance data set
url.sta <- "http://opendata.dc.gov/datasets/ab5661e1a4d74a338ee51cd9533ac787_50.geojson"
sta.spdf <- rgdal::readOGR(url.sta, layer = "OGRGeoJSON")

# station centroids, needed for report plotting
url.sta.ctrd <- "http://opendata.dc.gov/datasets/54018b7f06b943f2af278bbe415df1de_52.geojson"
sta.ctrd.spdf <- rgdal::readOGR(url.sta.ctrd, layer = "OGRGeoJSON")

# Metro lines data set, needed for report plotting
url.lines <- "http://opendata.dc.gov/datasets/a29b9dbb2f00459db2b0c3c56faca297_106.geojson"
lines.spdf <- rgdal::readOGR(url.lines, layer = "OGRGeoJSON")

#######################################
###      cluster the stations       ###
#######################################

wgs84 <- sta.spdf@proj4string
nad83md <- CRS("+init=esri:102285")

source("helper_functions.R")

# for various numbers of clusters, how large a circle can we draw around
# each cluster centroid?
# k.dist <- max_radius(sta.spdf, max.dist = 500, min.dist = 450)

# DC blocks are roughly 150 m in length, so we're going to choose a 1-block
# radius. to draw a circle of radius r around a point such that the circles do
# not overlap for points A and B, we need A and be to be at least r * 2 apart.
min.dist <- 300

# we might not cluster anything
ncl <- length(sta.spdf)
dif <- ncl
# do an initial clustering
x <- cluster_points(sta.spdf, min.dist)

while(dif > 0) {
  x <- sp::SpatialPointsDataFrame(x, data = data.frame(OBJECTID = seq_along(x)))
  x <- cluster_points(x, min.dist, verbose = TRUE)  
  dif <- ncl - length(x)
  ncl <- length(x)
}

# create an arbitrary (unique) centroid ID
x <- sp::SpatialPointsDataFrame(x, 
                                data = data.frame(CTRD_ID = seq_along(x)))

# now create the buffers using the centroids. we need to transform the 
# coordinates to a planar projection before we can use gBuffer; may as
# well use MD state plane NAD 83.
ctrd.buffer <- rgeos::gBuffer(spgeom = sp::spTransform(x, nad83md), 
                              byid = T, width = min.dist / 2)

ctrd.buffer <- sp::spTransform(sp::disaggregate(ctrd.buffer), wgs84)

# add centroid ID back, same as we assigned for t above
ctrd.buffer$CTRD_ID <- sp::over(ctrd.buffer, x)$CTRD_ID

# for each centroid, add the ward number
ctrd.spdf <- sp::SpatialPoints(ctrd.buffer, proj4string = ctrd.buffer@proj4string)
ctrd.buffer$WARD <- sp::over(ctrd.spdf, ward.spdf)$WARD

# add centroid IDs of clusted stations
ctrd.id <- get_ctrd_id(sta.spdf, x)
sta.spdf <- merge(sta.spdf, ctrd.id)

# a little spatial subsetting, i.e., station entrances by ward
# sta.spdf[ward.spdf[which(ward.spdf@data$WARD == 8, arr.ind = T),],]

#######################################
###        add the crime data       ###
#######################################

# load crime data, obtained from MPD database
crime <- dplyr::tbl_df(read.csv("alldc_2014-01-01_2015-07-02.csv"))

# get the dates into a sensible format
crime$REPORT_DAT <- as.POSIXct(crime$REPORT_DAT, format = "%m/%d/%Y %I:%M:%S %p")
crime$START_DATE <- as.POSIXct(crime$REPORT_DAT, format = "%m/%d/%Y %I:%M:%S %p")
crime$END_DATE <- as.POSIXct(crime$REPORT_DAT, format = "%m/%d/%Y %I:%M:%S %p")

# find the day of the week on which the crime was reported
crime$REPORT_DAY <- lubridate::wday(crime$REPORT_DAT, label = T)

# drop non-weekend crimes
crime.wknd <- crime[crime$REPORT_DAY %in% c("Sat", "Sun"),]
crime.wknd$OFFENSE <- as.character(crime.wknd$OFFENSE)

# add weekend start and end
crime.wknd$WKND_START <- as.POSIXct(trunc(crime.wknd$REPORT_DAT, units = "days") -
                                      lubridate::ddays(ifelse(crime.wknd$REPORT_DAY == "Sat", 0, 1)))

crime.wknd$WKND_END <- as.POSIXct(trunc(crime.wknd$REPORT_DAT, units = "days") +
                                    lubridate::ddays(ifelse(crime.wknd$REPORT_DAY == "Sun", 0, 1)))

# the data appear to be pretty clean from 2014 onward, but still
# check for (and remove) NA values; we can't do anything with them
crime.wknd <- dplyr::filter(crime.wknd, 
                            !is.na(XBLOCK), !is.na(YBLOCK), !is.na(OFFENSE))

# map in crime category; it will be nice to use dplyr to join tables, 
# but we'll have to do a little work to get there
crime.types <- unique(crime.wknd$OFFENSE)
violent <- c("HOMICIDE", "SEX ABUSE", "ROBBERY", "ASSAULT W/DANGEROUS WEAPON")
property <- crime.types[!(crime.types %in% violent)]
crime.wknd <-
  dplyr::bind_rows(data.frame(OFFENSE = violent, CATEGORY = "violent", 
                              stringsAsFactors = F),
                   data.frame(OFFENSE = property, CATEGORY = "property", 
                              stringsAsFactors = F)) %>%
  dplyr::inner_join(crime.wknd, by = "OFFENSE")

# all crimes should be categorized
any(is.na(crime.wknd$CATEGORY))

# now we need to transform the projection to WGS 84 to align with the 
# ward and station entrance data. from http://crimemap.dc.gov/Download.aspx:
# "Location is approximated to the center of the street block. Values are 
# in the Maryland State Plane meters NAD 83 map projection."

# thanks to https://github.com/eshilts/dc_crime_data, and to SO, where
# he got his answer:
# http://stackoverflow.com/questions/11649580/how-to-convert-nad-83-coordinates-to-latitude-and-longitude-with-rgdal-package
crime.spdf <- 
  sp::SpatialPointsDataFrame(as.data.frame(dplyr::select(crime.wknd, XBLOCK, YBLOCK)),
                             data = as.data.frame(crime.wknd),
                             proj4string = nad83md)
# now build the SpatialPointsDataFrame
crime.spdf <- sp::spTransform(crime.spdf, wgs84)

# now, map a centroid ID for each crime. NA means the crime doesn't fall in 
# any of the polygons in ctrd.buffer.
crime.spdf$CTRD_ID <- sp::over(crime.spdf, ctrd.buffer)$CTRD_ID
crime.spdf$WARD <- sp::over(crime.spdf, ctrd.buffer)$WARD

# now we need to 
# create a data set that consists of (centroid, weekend) tuples with the
# number of crimes by crime category.
df.wknd.sta <- as.data.frame(crime.spdf[which(!is.na(crime.spdf$CTRD_ID)), ])
df.wknd.sta <- dplyr::tbl_df(df.wknd.sta)

# generate a list of all weekends
int <- lubridate::interval(ymd(20140101, tz = ""), ymd(20150701, tz = ""))
# all the days in the period
alldays <- int@start + 
  lubridate::days(0:as.numeric(int_end(int) - int@start, unit = "days"))

wknd.all <- data.frame(WKND_START = alldays[lubridate::wday(alldays) == 7],
                       WKND_END = alldays[lubridate::wday(alldays) == 7] + 
                         lubridate::days(1))
wknd.all$ID <- seq_len(nrow(wknd.all))
wknd.all <- dplyr::tbl_df(wknd.all)

# add weekend ID
df.wknd.sta <- 
  df.wknd.sta %>% 
  dplyr::inner_join(wknd.all, by = c("WKND_START", "WKND_END")) %>% 
  dplyr::rename(WKND_ID = ID)

# write.csv(df.wknd.sta, file = "weekend_crimes.csv", row.names = F)

# now do the Cartesian join
ctrd.wknd.all <- merge(wknd.all, as.data.frame(ctrd.buffer), by = NULL)
ctrd.wknd.all <- merge(ctrd.wknd.all, 
                       data.frame(CATEGORY = c("violent","property")), by = NULL)
ctrd.wknd.all <- dplyr::tbl_df(ctrd.wknd.all)
ctrd.wknd.all$CATEGORY <- as.character(ctrd.wknd.all$CATEGORY)

# add latitude and longitude
ctrd.wknd.all <- 
  dplyr::tbl_df(as.data.frame(x)) %>%
  dplyr::select(CTRD_ID, CTRD_LAT = coords.x2, CTRD_LONG = coords.x1) %>%
  dplyr::inner_join(ctrd.wknd.all, by = "CTRD_ID")

# now join the crime data to the (centroid, weekend) data
grp.vars <- c("CTRD_ID", "WKND_START", "WKND_END", "CATEGORY")

data.final <- 
  df.wknd.sta %>%
  dplyr::group_by_(.dots = grp.vars) %>%
  dplyr::tally() %>%
  dplyr::right_join(ctrd.wknd.all, by = grp.vars) %>%
  dplyr::select(WKND_ID = ID, WKND_START, WKND_END, CTRD_ID, WARD,
                CTRD_LAT, CTRD_LONG, CATEGORY, NUM_CRIMES = n) %>%
  dplyr::ungroup()

# replace NA by zero
data.final$NUM_CRIMES[is.na(data.final$NUM_CRIMES)] <- 0

# and write out the data file
# save(data.final, file = "station_crimes.RData")

#######################################
###  various statistics for report  ###
#######################################

# how many violent crimes near stations on the weekend?
df.wknd.sta %>% dplyr::filter(CATEGORY == "violent") %>% dplyr::count()

# and how many property crimes near stations on the weekend?
df.wknd.sta %>% dplyr::filter(CATEGORY == "property") %>% dplyr::count()

# how many crimes near stations overall?
crime %>%
  dplyr::mutate(CATEGORY = dplyr::if_else(OFFENSE %in% violent, 
                                          "violent", "property")) %>%
  dplyr::group_by(CATEGORY) %>%
  dplyr::tally()

# how many stations have more than 1 exit?
as.data.frame(sta.spdf) %>%
  dplyr::group_by(GIS_ID) %>%
  dplyr::tally() %>%
  dplyr::filter(n > 1)

# how many stations have exits in different wards?
sta.spdf$STA_WARD <- sp::over(sta.spdf, ward.spdf)$WARD

as.data.frame(sta.spdf) %>%
  dplyr::distinct(GIS_ID, STA_WARD) %>%
  dplyr::group_by(GIS_ID) %>%
  dplyr::tally() %>%
  dplyr::filter(n > 1)

# what is the average length of a block?
# url.block <- "http://opendata.dc.gov/datasets/ba2539327dcf448789dc65a55ebe3d16_5.geojson"
# block.spdf <- rgdal::readOGR(url.block, layer = "OGRGeoJSON")

# table of stations with centroid ID for report
# sta.tbl <-
#   as.data.frame(sta.spdf) %>%
#   dplyr::inner_join(tbl_df(as.data.frame(ctrd.buffer))) %>%
#   dplyr::mutate(LONGITUDE = formatC(coords.x1, digits = 6, format = "f"),
#                 LATITUDE = formatC(coords.x2, digits = 6, format = "f")) %>%
#   dplyr::select(WARD, CTRD_ID, NAME, LONGITUDE, LATITUDE) %>%
#   dplyr::arrange(WARD, CTRD_ID, NAME)

# write.csv(sta.tbl, file = "sta_tbl.csv", row.names = F)

#######################################
###          test plotting          ###
#######################################

# plot the centroids and label them by ward

# plot(ward.spdf)
# add the ward numbers
# text(as.data.frame(rgeos::gCentroid(ward.spdf, byid = T)),
#      labels = ward.spdf$WARD, cex = 0.75, col = "red", font = 2)
# now plot the centroids, checking that they are in the correct wards
# plot(ctrd.buffer, add = T)
# text(ctrd.spdf@coords, labels = ctrd.buffer$WARD, cex = 0.5, col = "blue")

# plot crimes to make sure they are near the stations
# plot(ward.spdf)
# plot(ctrd.buffer, add = T)
# plot(crime.spdf[which(!is.na(crime.spdf$CTRD_ID)), ], add = T)

# plot crimes in a single ward
# plot(ward.spdf)
# text(as.data.frame(gCentroid(ward.spdf, byid = T)),
#      labels = ward.spdf$WARD, cex = 0.75, col = "red", font = 2)
# plot(crime.spdf[which(!is.na(crime.spdf$CTRD_ID) & 
#                         crime.spdf$WARD == 6), ], add = T)

# plot crimes by centroid
# plot(ward.spdf)
# plot(ctrd.buffer[ctrd.buffer$CTRD_ID == 7,], add = T)
# plot(crime.spdf[which(!is.na(crime.spdf$CTRD_ID) & 
#                         crime.spdf$CTRD_ID == 7), ], add = T)
