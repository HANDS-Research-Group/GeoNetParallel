
library(rgdal)
## Update the path to the folder
file_path <- "C:/Users/rohit/OneDrive - Syracuse University/GeoNet/Repo/GeoNet2022/"

shape = rgdal::readOGR(dsn=paste0(file_path, "data/Flowline_R_Package"), layer="Flowline_R_Package")
save(shape, file = paste0(file_path, "shape.RData"))
