##
## Creatd on 2017-SEP-27 
## Created by Chel Hee Lee <chelhee.lee@ucalgary>
## Updated on 2017-OCT-11
## 
## Description:  Practice mapping techniques using R
## 
## Notes:  Sorry, the codes are not cleaned. 
## 

rm(list=ls(all.names=TRUE))

##
## Read the Calgary Police Crime Stats 
## 

cstats <- gdata::read.xls(xls="http://www.calgary.ca/cps/Documents/statistical-reports/2017%20Community%20Crime%20Statistics.xls", sheet="By Category", header=TRUE, stringsAsFactors=FALSE, skip=1)

##
## structuring data for my interest 
## 

str(cstats)
names(cstats)
tmp <- subset(cstats, select=names(cstats)[2:14], subset=(X == "Assault (Non-domestic)"))
head(tmp)
tmp$sum12 <- rowSums(tmp[,-1], na.rm=TRUE)
tmp <- tmp[!(names(tmp) %in% toupper(month.abb))]
names(tmp)[1] <- "community"


##
## Defining a safety index
## 

pcut <- c(0, 0.6, 0.9, 0.95, 0.99, 0.994, 1)
breaks <- quantile(tmp$sum12, probs=pcut)
k <- length(pcut)-1
idx <- cut(tmp$sum12, breaks=breaks, include.lowest=TRUE, labels=1:k)
cc <- RColorBrewer::brewer.pal((length(breaks)-1), "Greens") # cc: color code
tmp$cc <- cc[idx]


##
## Reading the GIS info of Community Boundary from the City of Calgary
## 

download.file(url="https://data.calgary.ca/api/geospatial/ab7m-fwn6?method=export&format=Shapefile", destfile="~/Downloads/community_boundary.zip", method="wget")
unzip(zipfil="~/Downloads/community_boundary.zip", exdir="~/Downloads/")
flist <- list.files("~/Downloads")
fname <- flist[grep(".shp", flist)]
cb <- maptools::readShapePoly(file.path("~/Downloads", fname))

##
## mapping
## 

# pdf(file="./map2012_assault_calgary.pdf")
png(file="./map2012_assault_calgary.png")
sp::plot(cb, axes=TRUE)
for(i in 1:k){
    cb1 <- cb[which(cb$name %in% tmp$community[tmp$cc == cc[i]]),]
    sp::plot(cb1, add=TRUE, border=TRUE, col=cc[i])
}
legend(x=-114.3046, y=51.00231, legend=maptools::leglabs(floor(breaks)), fill=cc, bty="n", cex=0.8)
points(-114.0890, 51.0413, col="red", cex=1.5, pch=19)  # beltline
points(-114.0691, 51.0460, col="blue", cex=1.5, pch=19) # downtown commercial core 

## Annotation for Calgary Police Crime Data 
mtext("Naive number of assaults in 2012 from Calgary Police Crime Statistics \n (http://www.calgary.ca/cps/Documents/statistical-reports/2017%20Community%20Crime%20Statistics.xls)", side=3, line=0, cex=0.7) 

## Annotation for Open data license from the city of Calgary 
text(x=-114.3, y=50.84, labels=paste("Contains information licensed under the Open Government Licence", "City of Calgary", sep=" - "), adj=0, cex=0.6)
text(x=-114.3, y=50.83, labels=paste("License from ", "https://data.calgary.ca/stories/s/u45n-7awa/"), adj=0, cex=0.6)
text(x=-114.3, y=50.82, labels=paste("Community Boundary Shape File from ", "https://data.calgary.ca/Base-Maps/Community-Boundaries/ab7m-fwn6/data"), adj=0, cex=0.6)

text(x=-114.32, y=51.21, labels=c("This example is presented by Chel on 2017.09.27 at CalgaryR for the purpose of practicing mapping techniques.\n 'maptools', 'RColorBrewer', 'sp' packages are used.  \n TODO, FIX, WISHLIST are in future meetings."), col="red", cex=0.6, adj=c(0,0))
dev.off()

##
## TODO::2017-OCT-11
## + Cannot be listed here since too many things to do (sorry)
## 