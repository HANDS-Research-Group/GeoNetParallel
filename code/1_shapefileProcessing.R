
library(rgdal)
## Update the path to the folder
file_path <- "C:/GeoNet/GeoNet_2021_packageDataset/"

shape = rgdal::readOGR(dsn=paste0(file_path, "data/Flowline_R_Package"), layer="Flowline_R_Package")
save(shape, file = paste0(file_path, "shape.RData"))
