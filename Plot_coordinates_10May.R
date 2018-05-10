#Estimating the coordinates of all plots from the coordinates of the corner plots
#"Waypoints_01-MAY-18.gpx" from day trip to knersvlakte 1 May
#Honours 2018
#10 May

#----------------------------------------------------

rm(list = ls())
setwd("/Users/larawootton/Documents/Honours/Data")

if(!require(pacman)){install.packages("pacman", dependencies=TRUE); library(pacman)}
p_load(dplyr, plotKML, sp)

#read in gps data
waypts <- readGPX("Waypoints_01-MAY-18.gpx")
glimpse(waypts)

#Extract the waypoints

my_points <- waypts$waypoints
glimpse(my_points)


#Generate the lon and lat of all plots along the right and LH sides of the three sites

site2_rh_lon <- seq(my_points$lon[my_points$name == "006"], my_points$lon[my_points$name == "007"], length.out = 11)
site2_rh_lat <- seq(my_points$lat[my_points$name == "006"], my_points$lat[my_points$name == "007"], length.out = 11)
site2_lh_lon <- seq(my_points$lon[my_points$name == "005"], my_points$lon[my_points$name == "004"], length.out = 11)
site2_lh_lat <- seq(my_points$lat[my_points$name == "005"], my_points$lat[my_points$name == "004"], length.out = 11)

site1_rh_lon <- seq(my_points$lon[my_points$name == "012"], my_points$lon[my_points$name == "013"], length.out = 11)
site1_rh_lat <- seq(my_points$lat[my_points$name == "012"], my_points$lat[my_points$name == "013"], length.out = 11)
site1_lh_lon <- seq(my_points$lon[my_points$name == "011"], my_points$lon[my_points$name == "010"], length.out = 11)
site1_lh_lat <- seq(my_points$lat[my_points$name == "011"], my_points$lat[my_points$name == "010"], length.out = 11)


site3_rh_lon <- seq(my_points$lon[my_points$name == "015"], my_points$lon[my_points$name == "014"], length.out = 11)
site3_rh_lat <- seq(my_points$lat[my_points$name == "015"], my_points$lat[my_points$name == "014"], length.out = 11)
site3_lh_lon <- seq(my_points$lon[my_points$name == "016"], my_points$lon[my_points$name == "017"], length.out = 11)
site3_lh_lat <- seq(my_points$lat[my_points$name == "016"], my_points$lat[my_points$name == "017"], length.out = 11)




plot(site2_rh_lon, site2_rh_lat, xlim= c(18.5464, 18.549 ), ylim = c(-31.267, -31.2649))
points(site2_lh_lon, site2_lh_lat)
points(site1_lh_lon, site1_lh_lat)
points(site1_rh_lon, site1_rh_lat)
points(site3_lh_lon, site3_lh_lat)
points(site3_rh_lon, site3_rh_lat)

