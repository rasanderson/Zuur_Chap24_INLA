# This file was taken from the support material from Jagger et al. (2015).

# Jagger TH, Elsner JB, Widen HM (2015). A statistical model for regional
# tornado climate studies. PLoS ONE 10(8): e0131876. doi:10.1371/
# journal.pone.0131876


#This contains support functions for the countywide statistics project.
print(paste("Sourcing CountyStatisticsSupport.R","on",date(),sep=" "))
print("loading packages")
print("Version 1.2 Fri Dec  5 13:19:14 2014")
slibrary(sp)
slibrary(rgeos)
slibrary(rgdal)
slibrary(plyr)
slibrary(maps)
slibrary(moments)
slibrary(spdep)
slibrary(maptools)
slibrary(mapproj)
slibrary(raster) 
slibrary(stringr)
slibrary(scales)
slibrary(munsell)
slibrary(ggplot2)
slibrary(grid)
slibrary(INLA)
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#  Author: winston@stdout.org
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols), byrow=TRUE)
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
make.precip.df<-function(x=precip.raw){
  precip.split = lapply(x, strsplit, "\t")
  precip.header = precip.split[[1]][[1]]
  precip.header[c(1, 2)] = c("county", "Year")
  lengthprecip = sapply(precip.split, function(x) length(x[[1]]))
  precip.df = do.call("rbind", lapply(precip.split[-c(1, which(lengthprecip == 14))],
                                      function(x) tolower(str_trim(x[[1]]))))
  cn = rle(precip.df[, 1])
  precip.df[,1] = rep(cn$values[cn$values != ""], each = sum(cn$lengths[1:2]))
  colnames(precip.df) = str_trim(precip.header)
  precip.df      
}
getBordersAndCenters<-function(x,CRS0){  
  centers.sp = rgeos::gCentroid(x, byid=TRUE)
  centers.sp = spTransform(centers.sp, CRS0)
  centers.df = as.data.frame(coordinates(centers.sp)) 
  names(centers.df) = c("long", "lat")
  centers.df$county = rownames(centers.df)
  borders.df = fortify(spTransform(x,CRS0),region="county")
  borders.df$county = borders.df$id
  return(list(borders=borders.df,centers=centers.df))
}
plotTorn<-function(x,tracks){
  CRS0=CRS("+proj=longlat +datum=WGS84")
  #x is the spatial polygons object
  #tracks is the track file for each polygon.
  bc<-getBordersAndCenters(x,CRS0)
  bc$centers$nT=x$nT
  bc$centers$rate = round(x$rate, 2)
  borders.df = merge(bc$centers[c("county","nT","rate")], bc$borders, by = "county")

  #Use this data to make the plot of intensity
  mapT = ggplot(borders.df, aes(x = long, y = lat, group = county, fill = nT)) +
    geom_polygon(colour = "gray", size = .2) +
    coord_map(project="polyconic") + 
    theme(panel.grid.minor=element_blank(), 
          panel.grid.major=element_blank(), 
          panel.background=element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.position="none") +
    #  labs(title="Number of Tornadoes\n(1954-2012)") + 
    xlab("") + ylab("") +
    scale_fill_gradient() +
    geom_text(data=bc$centers, aes(x = long, y = lat, group=NULL, label = nT), 
              vjust = .5, size = 5, color = "white")
  #generate the tracks data from the track.df file 
  tracks = tracks[!duplicated(tracks$OM), ]
  tracks = subset(tracks, ELAT != 0)
  mapT +
    geom_segment(aes(x = SLON, xend = ELON, y = SLAT, yend = ELAT, 
                     group = NULL, fill = NULL),
                 data = tracks, lineend = "round",
                 color = "orangered", size = 1.5, alpha = .6) +
    geom_text(data=bc$centers, aes(x = long, y = lat, group = NULL, label = nT), 
              vjust = .5, size = 5, color="white")
}
read.df<-function(x,values="t"){
  for(i in values) x[x==i]<-NA
  if(is.matrix(x)) x<-data.frame(x,stringsAsFactors=FALSE)
  out<-data.frame(lapply(x,function(x) if(is.character(x)) type.convert(x,as.is=T) else x ) ,stringsAsFactors=FALSE) 
  if(!is.null(colnames(x))) colnames(out)<-colnames(x)
  if(!is.null(rownames(x))) rownames(out)<-rownames(x)                  
  out
}
processCity<-function(pop,field){
  ctys = sapply(strsplit(pop$GEO.display.label, ","), function(x) x[1])
  ctys2 = tolower(sapply(strsplit(ctys, " "), function(x) x[1]))
  pop = pop[order(ctys2), ]
  if(!all(tolower(substring(pop$GEO.display.label, 1, nchar(spdf$county))) == spdf$county )) stop("processCity failed mismatch countynames")
  pop[[field]]
}
processPop<-function(x){
  #Process state county yearly population data
  #Similar to melt() but shows details
years<-as.numeric(colnames(x)[-1]) #year values
data<-as.matrix(x[,-1])
nc = nrow(data)  #number of counties nc
ny = ncol(data)  #number of years    ny should equal length(years)
county<-tolower(gsub(" Co.","",x[[1]],fixed=T)) #get lowercase county names
data1<-data.frame(county=rep(county,ny),Year=rep(years,each=nc),pop=as.vector(data),stringsAsFactors=F)
return(data1[with(data1,order(county,Year)),])
}
changeTornDataTypes<-function(x=TornL){
  #converts character data to strings, length from miles to meters and width from yards to meters.
  #adds a field called sequence id, which can be used to identify similar events split over multiple counties
  x$Year = as.integer(x$YR)
  x$Month = as.integer(x$MO)
  x$EF = as.integer(x$MAG) 
  x$Date = as.Date(x$DATE)
  x$Length = as.numeric(x$LEN) * 1609.34
  x$Width = as.numeric(x$WID) * 0.9144
  x$FAT = as.integer(x$FAT)
  x$SLON = as.numeric(x$SLON)
  x$SLAT = as.numeric(x$SLAT)
  x$ELON = as.numeric(x$ELON)
  x$ELAT = as.numeric(x$ELAT)
  x$INJ = as.numeric(x$INJ)
  x$SeqID = 1:nrow(x)
  return(x)
}
make.equal.positions<-function(nx,ny,by.row=T,topdown=T){
  xval<-seq(0,1,length=nx+1)
  yval<-seq(0,1,length=ny+1)
  if(topdown) yval<-rev(yval)
  positions<-expand.grid(x=xval,y=yval)
  loc<-matrix(1:((nx+1)*(ny+1)),nx+1,ny+1)
  nx1<-1:nx
  ny1<-1:ny
  d=topdown*1  
  if(!by.row)
  out<-do.call("rbind",lapply(nx1,function(ix) do.call("rbind",lapply(ny1, function(iy) cbind(positions[loc[ix,iy+d],],positions[loc[ix+1,iy+1-d],])))))
  else
    out<-do.call("rbind",lapply(ny1,function(iy) do.call("rbind",lapply(nx1, function(ix) cbind(positions[loc[ix,iy+d],],positions[loc[ix+1,iy+1-d],])))))
  out<-(t((out)))
  colnames(out)<-paste("fig",(1:(nx*ny)),sep=".")
  rownames(out)<-c("xstart","ystart","xend","yend")
  data.frame(out)
}
multipleplot<-function(tlist,nx,ny,by.row=T){
  if(length(tlist) != nx*ny) warning("Missing Corners")
  positions<-make.equal.positions(nx,ny,by.row)
  for(i in 1:(length(tlist)))
  {
    more = (i != length(tlist))
    print(tlist[[i]],position=positions[[i]],more=more)
  }        
}    
ggplotSpdata<-function(x,field, tfun=function(x) I(x), color=muted("yellow",l=50,c=80)){
  CRS0=CRS("+proj=longlat +datum=WGS84")
  bc = getBordersAndCenters(x,CRS0)
  bc$centers$newfield = tfun(x[[field]])
  borders.df = merge(bc$borders, bc$centers[c("county","newfield")], by="county")
  mapO = ggplot(borders.df, aes(x=long, y=lat, group=group, fill = newfield)) +
    geom_polygon(colour = "gray", size = .2) +
    coord_map(project="polyconic") + 
    theme(panel.grid.minor=element_blank(), 
          panel.grid.major=element_blank(), 
          panel.background=element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.position="none") +
    xlab("") + ylab("") +
    scale_fill_gradient2() +
    geom_text(data=bc$centers, aes(x = long, y = lat, group=NULL, label = newfield), 
              vjust = .5, size = 5, color = color)
    mapO
}  
ggplotmargin<-function(x,type,effect,xlab,ylab=expression(Posterior~~Density),int.value=c(value=0,2.5,97.5),color=c("red","green","blue")){
xx<-as.data.frame(inla.smarginal(x[[paste("marginals",type,sep=".")]][[effect]]))
  out<-ggplot(xx, aes(x, y)) + geom_line() +  ylab(ylab) + xlab(xlab)    
if(length(int.value)==0) int.value=0
int.value <-lapply(int.value,function(x) if(is.character(x)) type.convert(x,as.is=T) else x)
int.value <-lapply(int.value, function(x) if(is.character(x)) 
            lapply(strsplit(x,"=")[[1]], type.convert, as.is=T) else x)
nx<-names(int.value)
if(!is.null(nx))
   for(i in which(nx !=""))  int.value[[i]] <- list(nx[i],int.value[[i]])
    int.value<-sapply(int.value, function(x) {
                      if(is.numeric(x)) xx$x[which.max(cumsum(xx$y)/sum(xx$y) >= as.numeric(x/100))]
                      else switch(x[[1]], 
                      mean= sum(xx$y*xx$x)/sum(xx$y), 
                      median=xx$x[which.max(cumsum(xx$y)/sum(xx$y) >=.5)],
                      mode=xx$x[which.max(xx$y)],
                      value=x[[2]],
                      zero=0)})

if(length(color) <= length(int.value)) color<-rep(color,length=length(int.value))
for(i in 1:length(int.value)) out<- out + geom_vline(xintercept = int.value[i], color = color[i]) 
out
}

printinla <-  function (x, call=F,time=F, notes=F, digits=3,...) {
    if(class(x)=="inla") x<-summary(x,digits=digits)
    digits = digits
    if(call){    cat("\nCall:\n", inla.formula2character(x$call), "\n\n", sep = "")}
    if(time) {if (inla.is.element("cpu.used", x)) {
      cat("Time used:\n")
      print(x$cpu.used)
      cat("\n")
    }
    }
    if (inla.is.element("fixed", x)) {
      cat("Fixed effects:\n")
      print.default(x$fixed)
      cat("\n")
    }
    else {
      cat("The model has no fixed effects\n\n")
    }
    if (inla.is.element("lincomb.derived", x)) {
      cat("Linear combinations (derived):\n")
      print.default(x$lincomb.derived)
      cat("\n")
    }
    if (inla.is.element("lincomb", x)) {
      cat("Linear combinations:\n")
      print.default(x$lincomb)
      cat("\n")
    }
    if (inla.is.element("random.names", x)) {
      cat("Random effects:\n")
      cat("Name\t ", "Model\n ")
      for (i in 1:length(x$random.names)) cat(paste(inla.nameunfix(x$random.names[i]), 
                                                    " ", x$random.model[i], "\n"))
      cat("\n")
    }
    else {
      cat("The model has no random effects\n\n")
    }
    if (inla.is.element("hyperpar", x)) {
      cat("Model hyperparameters:\n")
      print(format(x$hyperpar, digits = digits, nsmall = 2), 
            quote = FALSE)
      cat("\n")
    }
    else {
      cat("The model has no hyperparameters\n\n")
    }
    if (inla.is.element("neffp", x)) {
      cat("Expected number of effective parameters(std dev): ", 
          format(x$neffp[1], digits = digits, nsmall = 2), 
          "(", format(x$neffp[2], digits = digits, nsmall = 2), 
          ")\n", sep = "")
      cat("Number of equivalent replicates :", format(x$neffp[3], 
                                                      digits = digits, nsmall = 2), "\n")
      cat("\n")
    }
    else {
      cat("Expected number of effective parameters and Number of equivalent replicates not computed\n")
      cat("\n")
    }
    if (inla.is.element("dic", x)) {
      cat(paste("Deviance Information Criterion: ", format(x$dic$dic, 
                                                           digits = digits, nsmall = 2), "\n", "Effective number of parameters: ", 
                format(x$dic$p.eff, digits = digits, nsmall = 2), 
                "\n", sep = ""))
      cat("\n")
    }
    if(notes){
      if (inla.is.element("mlik", x)) {
      cat(paste("Marginal Likelihood: ", format(x$mlik[2], 
                                                digits = digits, nsmall = 2), "\n"))
    }
    if (inla.is.element("cpo", x)) {
      cat("CPO and PIT are computed\n")
      cat("\n")
    }
    if (inla.is.element("linear.predictor", x)) {
      cat("Posterior marginals for linear predictor and fitted values computed\n\n")
    }
    }
  }
environment(printinla)<-environment(INLA:::print.summary.inla)

#The following functions below allow one to aggregate over large raster files, as long as small sections can be read.
mytabulate<-function(x,levels=NULL){
    if(!is.null(levels)) {
      values<- match(x,levels,nomatch=NA)
      if(any(is.na(values))) stop("Additional Value found in data:",unique(x[is.na(values)])) 
    }
    else stop("Need levels for mytabulate")
    out<-tabulate(values,nbins=length(levels))
    names(out)<-levels
    out
  }
aggRasterOverSpatialPolygon<-function(r,poly,myfun=mytabulate,...,printinfo=T){
  #r is a raster (layer)
  #values r@data@min:r@data@max
  #poly is a spatialPolygons (or dataframe)
  #Returns a dataframe, but could be a Spatial Polygons dataframe
  #function returns a vector of named values for each polygons object in the SpatialPolygons object
  #poly has polygons which has Polygons (so we have functions of polygon and Polygon)
  #We assume one layer here.... More to follow ...
  if(printinfo) print(r@data@attributes[[1]][r@data@attributes[[1]]$COUNT>0,])
  polygons<-lapply(poly@polygons,getValues,r=r,myfun=myfun,...)
  names(polygons)<-sapply(poly@polygons, slot, name="ID")
  do.call("rbind",polygons)
}
getExtentMulti<-function(x,y){ 
  #Find cell locations in x from extents (bbox) in spatial polygons y
  bboxes<-lapply(y@polygons,function(xx) getextent(x,xx))
  names(bboxes) <-sapply(y@polygons, function(xx) xx@ID)
  bboxes
}
getExtent<-function(x,y) {
  z<-getRowColLoc(xx=bbox(x),yy=bbox(y),dimrc=dim(x))
  z[1,]<-z[1,2:1]
  return(cbind(min=trunc(z[,1]+1),max=ceiling((z[,2]))))
}
getValues<-function(r,polygon,myfun=mytabulate,...){
  bbox<-getExtent(r,polygon)
  values<-getValuesBlock1(r,bbox)
  rowcol<-getCoordsExtentCenters(bbox) 
  cells <- rowSums(sapply(polygon@Polygons,function(x) getCells(rowcol,getPoints(r,x@coords)) *(1-2*x@hole)))
  chosen<-values[cells>0]
  myfun(x=chosen,...)
}  
getValuesBlock1<-function(r,bbox){
  getValuesBlock(r,format="",bbox[1,1],bbox[1,2]-bbox[1,1]+1,bbox[2,1],bbox[2,2]-bbox[2,1]+1)
}
getRowColLoc<-function(xx,yy,dimrc) { #xx is the xy extent of raster, yy are points im t(xy),  returns coordinates one column for each pt by 2 rows}
  z<-(yy-xx[,1]) /(xx[,2]-xx[,1])
  rbind(row=(1-z[2,]),col=z[1,])*dimrc[1:2]
}
getPoints<-function(r,pts){
  #if the row is the x coordinate and column is y coordiante, we do not need to rotate poly for direction as
  #the transformation from x, y to row,col is a 90 degree counter clockwise rotation within the 0,1 by 0,1 scaled 
  #location withinraster extent box times the dimensions of the raster. 
  t(getRowColLoc(xx=bbox(r),yy=t(pts),dimrc=dim(r)))  
}
getCoordsExtentCenters<-function(x,offset=c(-.5,-.5)){
  #if the coordinates are the row column centers relative to the plot extent of [0,Nrow]x [0,Ncol]
  #and the upper left is row and column (1,1) and lower right is (Nrow,Ncol) then the centers are given here
  x<-round(x)
  if(x[1,1]>x[1,2] || x[2,1] > x[2,2]) {warning("Empty extent") ; return(NULL)}
  rows<-x[1,1]:(x[1,2])+offset[1]
  cols<-x[2,1]:(x[2,2])+offset[2]
  lr<-length(rows)
  lc<-length(cols)
  return(data.frame(row=rep.int(rows,rep.int(lc,lr)),col=rep.int(cols,lr)))
}
getCells<-function(pts,poly){  
  sp::point.in.polygon(pts$row,pts$col,poly[,1],poly[,2])  
}


#use point.in.polygon to determine if a point is in polygon. 

#How to extract values

#1. Get cells in extent of polygon
#2. Get cell locations in extent
#3. Use over to get index 1 in 0 out.
#4. Calculate statistic.

samplePlot<-function(){
f = "~/nlcd_2011_landcover_2011_edition_2014_03_31/nlcd_2011_landcover_2011_edition_2014_03_31.img"
r = raster(f)
rr<-GDAL.open(r@file@name)
tmpv<-getValuesBlock(r,row=50000,nrow=4000,col=40000,ncol=50000,format="matrix")
tmpvv<-t(tmpv[3000:1,16000:20000])
whichcols<-do.call("seq",as.list(range(tmpvv)))+1 #for colors
x11()
rr<-GDAL.open(r@file@name)
cols<-getColorTable(rr, band = 1)
image(tmpvv,col=cols[whichcols])
return(list(tmpvv=tmpvv,cols=cols, whichcols=whichcols))
}

generatedayofweekdata<-function(x,lon=-100, EF0=1, Year0=1954)
{
  t2=c("Sun","Mon","Tue","Wed","Thu","Fri","Sat")
  x<-subset(x,SLON>=lon & Year >= Year0 & EF >= EF0)@data
  xt<-table(x$Year,  factor(weekdays(x$Date,abbrev=T),levels=t2), x$ST)
  xm<-melt(xt)
  colnames(xm)<-c("Year","Day","State","Count")
  xm$IDyear<-xm$Year-Year0+1
  xm$IDday<-as.numeric(xm$Day)-1
  xm
}

inlaboxplot<-function(x,type="random",effect="IDday",levels=c("Sun","Mon","Tue","Wed","Thur","Fri","Sat"),p=c(.025,.25,.5,.75,.975),xlab="Day",ylab="effect"){
   marginals<-data.frame(t(sapply(x[[paste("marginals",type,sep=".")]][[effect]],function(x) inla.qmarginal(p=p,marginal=x) )))
   marginals$x=factor(levels,levels=levels)
   ggplot(data=marginals)+geom_boxplot(aes(ymin=X1,lower=X2,middle=X3,upper=X4,ymax=X5,x=x),stat="identity")+xlab(xlab)+ylab(ylab)                        
}

testday<-function(x,lon=-100,EF0=1,Year0=1954,formula=Count~f(IDyear,model="rw2")+State+f(IDday,model="seasonal",season.length=7),family="nbinomial"){
testdata<-generatedayofweekdata(x=x,lon=lon,EF0=EF0,Year0=Year0)
testinla<-inla(family=family,data=testdata,formula=formula)
plot<-inlaboxplot(x=testinla)
return(list(data=testdata,inla=testinla,plot=plot))
}

#End May 24, 2104 Thomas Jagger. 