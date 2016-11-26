#######################################
###         report plotting         ###
#######################################

library(dplyr)
library(ggplot2)
library(maptools)
library(rgeos)

# as calculated in the sampling script
ybars <- data.frame(y_bar = ybars_SRS, sample_type = "SRS")
ybars <- rbind(ybars, 
               data.frame(y_bar = ybars_Stratified, sample_type = "Stratified"))

#######################################
###   density estimate for \bar{y}  ###
#######################################

# as sampled in the sampling script
hist_brks <- c(seq(0.20, 0.40, 0.05), mean(surveyData$PROPERTY_CRIMES))
hist_brks <- hist_brks[order(hist_brks)]
hist_labs <- expression("0.20", "0.25", bar(y)[U], "0.30", "0.35", "0.40")

g <- ggplot2::ggplot(ybars, aes(x = y_bar, y = ..density..))
g <- g + ggplot2::geom_density(data = ybars[ybars$sample_type == "SRS",],
                               alpha = 0.5, aes(fill = "SRS"))
g <- g + ggplot2::geom_density(data = ybars[ybars$sample_type == "Stratified",],
                      alpha = 0.5, aes(fill = "Stratified"))
g <- g + ggplot2::geom_vline(xintercept = mean(surveyData$PROPERTY_CRIMES))
g <- g + ggplot2::scale_x_continuous(breaks = hist_brks, labels = hist_labs)
g <- g + ggplot2::labs(x = expression(bar(y)))
g <- g + ggplot2::scale_fill_manual(values = c("SRS" = "blue", 
                                               "Stratified" = "orange"),
                                    name = "")
print(g)

#######################################
###           metro lines           ###
#######################################

ward.df <- ggplot2::fortify(ward.spdf, region = "WARD")
ward.df$id <- as.numeric(ward.df$id)
ward.df <- dplyr::rename(ward.df, WARD = id)

lines.df <- ggplot2::fortify(lines.spdf, region = "GIS_ID")
lines.df <- dplyr::rename(lines.df, mline = id)
lines.df$mline <- as.numeric(lines.df$mline) + 1

lines.df <- 
  as.data.frame(lines.spdf) %>%
  dplyr::select(OBJECTID, GIS_ID, NAME) %>%
  dplyr::inner_join(lines.df, by = c("OBJECTID" = "mline"))

plot(ward.spdf)
text(as.data.frame(rgeos::gCentroid(ward.spdf, byid = T)), 
     labels = ward.spdf$WARD, cex = 0.75, col = "black", font = 2)
plot(lines.spdf[lines.spdf$GIS_ID == "Metro_003",], add = T, col = "orange",
     lty = 1, lwd = 2)
plot(lines.spdf[lines.spdf$GIS_ID == "Metro_001",], add = T, col = "blue",
     lty = 1, lwd = 2)
plot(lines.spdf[lines.spdf$GIS_ID == "Metro_004",], add = T, col = "red",
     lty = 1, lwd = 2)
plot(lines.spdf[lines.spdf$GIS_ID == "Metro_005",], add = T, col = "yellow",
     lty = 1, lwd = 2)
plot(lines.spdf[lines.spdf$GIS_ID == "Metro_002",], add = T, col = "green",
     lty = 1, lwd = 2)
plot(sta.ctrd.spdf, add = T, pch = 18)

#######################################
###  property crimes near stations  ###
#######################################

ctrd.ward <- rgeos::gCentroid(ward.spdf, byid = T)
ctrd.ward$WARD <- sp::over(ctrd.ward, ward.spdf)$WARD
ctrd.ward <- as.data.frame(ctrd.ward)
ctrd.ward <- dplyr::select(ctrd.ward, WARD, long = x, lat = y)

ctrd.plot <- as.data.frame(x)

crime.df <- 
  as.data.frame(crime.spdf[which(!is.na(crime.spdf$CTRD_ID) & 
                                   crime.spdf$CATEGORY == "property" & 
                                   crime.spdf$WKND_START >= lubridate::ymd(20140701)),])

g <- ggplot2::ggplot(ward.df)
g <- g + ggplot2::geom_polygon(aes(x = long, y = lat, group = id), 
                               fill = "white", color = "black")
g <- g + ggplot2::geom_point(data = crime.df, 
                             aes(x = XBLOCK.1, y = YBLOCK.1, color = "Crime"))
g <- g + ggplot2::geom_point(data = ctrd.plot, aes(x = coords.x1, 
                                                   y = coords.x2, 
                                                   color = "Station centroid"), 
                             pch = 3, size = 4)
g <- g + ggplot2::geom_density2d(data = crime.df, aes(x = XBLOCK.1, 
                                                      y = YBLOCK.1))
g <- g + ggplot2::geom_text(data = ctrd.ward, aes(x = long, 
                                                  y = lat, label = WARD), 
                            color = "red", fontface = "bold")
g <- g + ggplot2::scale_color_manual(name = "", 
                                     values = c("Crime" = "black", 
                                                "Station centroid" = "red"))
g <- g + ggplot2::guides(color = guide_legend(override.aes = list(linetype = c(0,0), 
                                                                  shape = c(16,3), 
                                                                  size = c(2,4))))
g <- g + ggplot2::coord_equal()
g <- g + ggplot2::theme(axis.title.x = element_blank(),
               axis.title.y = element_blank(),
               axis.text.x = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.x = element_blank(),
               axis.ticks.y = element_blank(),
               legend.position = "bottom")
print(g)
