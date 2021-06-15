#    Highland Statistics Ltd.
#    www.highstat.com
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.



# Section 24.2
library(raster)
library(rgdal)
library(sp)
library(spdep)
library(INLA)
library(mgcv)
library(ggplot2)
source("HighstatLibV11.R")


#Illinois map and data
Illi.shp <- readOGR("Shapefiles/Tornado_Illinois.shp")
Torn <- read.csv(file = "Illinois_popNtor2.csv",
                 header = TRUE)

names(Illi.shp@data)
names(Torn)
head(Torn)
##########################################



# Figure 24.1
plot(Illi.shp)  
# text(coordinates(Illi.shp),
#           labels(Illi.shp$ID),
#     cex =0.7)




# Figure 24.3
Illi.nb <- poly2nb(Illi.shp)
Illi.nb
Coords <- coordinates(Illi.shp)


plot(Illi.shp, border=grey(0.5))  
plot(Illi.nb, 
     coords = Coords, 
     add = TRUE,
     pch = 16,
     lwd = 2)


nb2INLA("Illi.graph", Illi.nb)
Illi.adj <- paste(getwd(), "Illi.graph", sep ="/")

Illi.Inla.nb <- inla.read.graph(filename = "Illi.graph")
inla.debug.graph("Illi.graph")


#Not in the book:
# plot(Illi.Inla.nb, par(cex = 0.75))



##############
# 24.3 Tornado data


#Figure 24.4
par(mfrow = c(2,2), mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = Torn$nT,
     y = 1:nrow(Torn),
     xlab = "Tornado counts",
     ylab = "Order of the data")

plot(x = log(Torn$density),
     y = 1:nrow(Torn),
     xlab = "Population density",
     ylab = "Order of the data")

plot(x = Torn$elevS,
     y = 1:nrow(Torn),
     xlab = "Variation in elevation",
     ylab = "Order of the data")

plot(x = Torn$area,
     y = 1:nrow(Torn),
     xlab = "Area",
     ylab = "Order of the data")




cor( Torn$elevS, Torn$area )
cor( Torn$elevS, log(Torn$density) )

Torn$LogDensity <- log(Torn$density)

# Source HighstatLibV11.R
MyX <- c("LogDensity", "elevS", "area", "Year")
MyMultipanel.ggp2(Z = Torn,
                  varx = MyX,
                  vary = "nT",
                  ylab = "Number of tornados",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = FALSE) 

  


# Percentages of zeros
sum(Torn$nT == 0) / nrow(Torn)





# 24.5 CAR correlation
par(mfrow=c(1,1))
# iCAR correlation
XY <- coordinates(Illi.shp)
range(XY[,1])
range(XY[,2])
plot(Illi.shp,
     xlim = c(-48912.2,  57413.6),
     ylim = c(-157811.7,  100197.6))  

Poly3 = subset(Illi.shp,
               ID == 3 |
                 ID == 5 |
                 ID == 37 | 
                 ID == 42 |
                 ID == 60 |
                 ID == 91 |
                 ID == 98)
plot(Poly3, add=T, col = grey(0.7), lwd=1)

MyID <- Illi.shp$ID

Neigbours1 <- MyID == 3  | MyID== 5 | MyID == 37 | 
  MyID == 42 | MyID == 60 | MyID == 91 |
  MyID == 98

text(x = XY[Neigbours1,1],
     y =  XY[Neigbours1,2],
     MyID[Neigbours1],
     cex = 1)

text(x = XY[1,1],
     y =  XY[1,2],
     MyID[1],
     cex = 2)






# 24.6 GAM with iCAR in R-INLA for the tornado data

# 24.6.2 Running the Poisson GLM with iCAR in R-INLA
# Not in book
image(inla.graph2matrix(Illi.nb))
# plot(Illi.nb) # Uses Rgraphviz
# inla.debug.graph("Illi.graph")




# Standardize all parametric covariates
MyStd <- function(x) {(x - mean(x)) / sd(x)}

Torn$ElevS.std <- MyStd(Torn$elevS)
Torn$LogDensity.std <- MyStd(Torn$LogDensity)
Torn$Area.std <- MyStd(Torn$area)


library(mgcv)
sm <- smoothCon(s(Year, bs = "cr", k = 6, fx = TRUE), 
                data = Torn,  
                absorb.cons = TRUE)[[1]]
Year.cr <- sm$X
Year.cr.df <- data.frame(Year.cr)
names(Year.cr.df) <- paste0("Year.cr", 1:ncol(Year.cr))
head(Year.cr.df)
lcs.Year <- inla.make.lincombs(Year.cr.df)


f1 <- nT ~ ElevS.std + LogDensity.std + Area.std +
           Year.cr 
  
M1 <- inla(formula = f1, 
               family = "poisson", 
               data = Torn,
               control.compute = list(dic = TRUE, waic = TRUE))


f2 <- nT ~ ElevS.std + LogDensity.std + Area.std +
           Year.cr +
          f(ID, 
            model = "besag", 
            graph = Illi.adj,
            scale.model = TRUE) 

M2 <- inla(formula = f2, 
           lincomb = lcs.Year,
                family = "poisson", 
                data = Torn,
                control.compute = list(dic = TRUE, waic = TRUE))

mu2 <- M2$summary.fitted.values[,"mean"]
head(cbind(mu2, Torn$nT))
cor(mu2, Torn$nT)


summary(M2)



# And compare the models with DICs and WAICs
dic  <- c(M1$dic$dic, M2$dic$dic)   
waic <- c(M1$waic$waic, M2$waic$waic)
Z.out     <- cbind(dic, waic)
rownames(Z.out) <- c("Poisson GAM",  
                     "Poisson GAM + Besag")
Z.out



# 24.6.3 Covariate effects

round(M2$summary.fixed[, c("mean", "0.025quant", "0.975quant")], 3)


# Figure 24.7
Ns <- nrow(Torn)
f.Year    <- M2$summary.lincomb.derived[1:Ns + 0 * Ns, "mean"] 
SeLo.Year <- M2$summary.lincomb.derived[1:Ns + 0 * Ns,"0.025quant"] 
SeUp.Year <- M2$summary.lincomb.derived[1:Ns + 0 * Ns,"0.975quant"]

IYear <- order(Torn$Year)

MyData <- data.frame(
  mu   = c( f.Year[IYear]), 
  SeUp = c(SeUp.Year[IYear]), 
  SeLo = c(SeLo.Year[IYear]), 
  Xaxis = c(sort(Torn$Year)),
  ID    = factor(rep(c("Year smoother"), each = nrow(Torn))))


library(ggplot2)
p <- ggplot()
p <- p + xlab("Year") + ylab("Smoother")
p <- p + theme(text = element_text(size = 15))
p <- p + geom_line(data = MyData, 
                   aes(x = Xaxis, y = mu))

p <- p + geom_ribbon(data = MyData, 
                     aes(x = Xaxis, 
                         ymax = SeUp, 
                         ymin = SeLo),
                     alpha = 0.6)

my.ggp.yrange <- c(0, 0)
XPos <- c(Torn$Year)
XID  <- rep(c("Year smoother"), each = nrow(Torn) )

MyData2 <- data.frame(Y     = rep(my.ggp.yrange, each = nrow(Torn)),
                      Xaxis = XPos,
                      ID    = factor(XID))
p <- p + geom_text(data = MyData2,
                   aes(y = Y,
                       x = Xaxis,
                       label = "|"),
                   size = 1)
p





# 24.6.4 Spatial random effects
Illi.shp$u.pm <- M2$summary.random$ID$mean

u.pm <- M2$summary.random$ID$mean
Id  <- 1:102





# Plot the spatial random effect
Illi.shp$RandomMean <- M2$summary.random$ID$mean
Illi.shp$RandomQ50  <- M2$summary.random$ID$'0.5quant'
Illi.shp$sd         <- M2$summary.random$ID$sd


# Figure 24.8
library(RColorBrewer)
crq = brewer.pal(5, "GnBu")
range(Illi.shp$RandomMean)
rng = seq(-4, 3, 1)
spplot(Illi.shp, "RandomMean", col = "white", at = rng, 
       col.regions = crq,
       colorkey = list(space = "bottom", labels = paste(rng)),
       par.settings = list(axis.line = list(col = NA)),
       sub = "Spatial random effects")



df <- M2$summary.random$ID
names(df) = c("ID", "mean", "sd", "QL", "QM", "QH", "mode", "kld")
df$QL = (exp(df$QL) - 1) * 100
df$QM = (exp(df$QM) - 1) * 100
df$QH = (exp(df$QH) - 1) * 100
df$Sig = sign(df$QL) == sign(df$QH)
table(df$Sig)
Illi.shp$Sig = df$Sig







ggplotSpdata <- function(x, field, tfun = function(x) I(x), 
                         color = muted("yellow", l = 50, c = 80)){
  CRS0 = CRS("+proj=longlat +datum=WGS84")
  bc = getBordersAndCenters(x, CRS0)
  bc$centers$newfield = tfun(x[[field]])
  borders.df = merge(bc$borders, bc$centers[c("county", "newfield")], by = "county")
  mapO = ggplot(borders.df, aes(x = long, y = lat, group = group, fill = newfield)) +
    geom_polygon(colour = "gray", size = .2) +
    coord_map(project = "polyconic") + 
    theme(panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.position = "bottom") +
    xlab("") + ylab("") +
    geom_text(data = bc$centers, 
              aes(x = long, y = lat, group = NULL, label = newfield), 
              vjust = .5, size = 3, color = color)
  mapO
}


slibrary <- function(...) suppressMessages(library(...))

# This file was taken from the support material from Jagger et al. (2015).
source("CountyStatisticsSupport.R")

MygetBordersAndCenters <- function(x,CRS0){  
  centers.sp = rgeos::gCentroid(x, byid=TRUE)
  centers.sp = spTransform(centers.sp, CRS0)
  centers.df = as.data.frame(coordinates(centers.sp)) 
  names(centers.df) = c("long", "lat")
  centers.df$county = x$county  #rownames(centers.df)
  borders.df = fortify(spTransform(x,CRS0),region="county")
  borders.df$county = borders.df$id
  return(list(borders=borders.df,centers=centers.df))
}

CRS0 = CRS("+proj=longlat +datum=WGS84")
bc   = MygetBordersAndCenters(Illi.shp, CRS0)





df2  = bc$centers
df3  = cbind(df, df2)
names(df3)[5] = 'newfield'
borders.df = merge(bc$borders, 
                   df3[c("county", "newfield", "Sig")], 
                   by = "county")
borders2.df = borders.df[borders.df$Sig, ]



MyggplotSpdata <- function(x,field, tfun=function(x) I(x), color=muted("yellow",l=50,c=80)){
  CRS0=CRS("+proj=longlat +datum=WGS84")
  bc = MygetBordersAndCenters(x,CRS0)
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



p <- MyggplotSpdata(Illi.shp, "RandomMean", tfun = function(x) round((exp(x) - 1) * 100))
p <- p + scale_fill_gradient2( 
  low = muted("blue"), 
  mid = "white", 
  high = muted("red"),
  midpoint = 0, 
  space = "Lab", 
  na.value = "grey50", 
  guide = "colourbar",
  aesthetics = "fill") 

p
p <- p + geom_polygon(data = borders2.df[borders2.df$newfield > 0, ],
                      aes(x = long, y = lat, group = group, fill = newfield),
                      color = muted("red"), size = 2) 
p

p <- p +  geom_polygon(data = borders2.df[borders2.df$newfield < 0, ],
                       aes(x = long, y = lat, group = group, fill = newfield),
                       color = muted("blue"), size = 2) 
p <- p+
  geom_text(data = df3[df3$Sig,], 
            aes(x = long, y = lat, group = NULL, label = round(newfield)), 
            vjust = .5, size = 3, color = 'black')
p




# Figure 24.9
Ui <- M2$summary.random$ID$mean

u.breaks <- seq(min(Ui) * 1.01, max(Ui) * 1.01, length = 7)
u.group <- cut(Ui, 
               breaks = u.breaks, 
               include.lowest =TRUE)
table(u.group)


NCounties <- length(Ui)

ui <- M2$marginals.random$ID[1:NCounties]
MyFun0 <- function(x){x}
MyFun1 <- function(x){exp(x)}
MyFun2 <- function(x){(exp(x) - 1) * 100}
u <- lapply(ui, function(x) inla.emarginal(MyFun2, x))
u <- lapply(ui, function(x) inla.emarginal(MyFun0, x))

cbind(unlist(u), Ui)
u <- unlist(u)
zeta.cutoff <- seq(min(u) * 1.01, max(u) * 1.01, length = 7)
fzeta <- cut(unlist(u), 
             breaks = zeta.cutoff, 
             include.lowest =TRUE)
table(fzeta)


map.u <- data.frame(ID      = Illi.shp@data$ID, 
                    u.group = u.group)
Data2 <- attr(Illi.shp, "data")
attr(Illi.shp, "data") <- merge(Data2, map.u, by = "ID")
Illi.shp$u <- u


spplot(obj = Illi.shp, 
       zcol = "u.group",
       colorkey = TRUE,
       col.regions = heat.colors(6, alpha = 0.5),
       asp = 1,
       panel = function(x,y,z,subscripts,...) {
         panel.polygonsplot(x,y,z,subscripts,...)
         #sp.text(coordinates(Illi.shp), 
         #Illi.shp$ID[subscripts]
         #round(Illi.shp$u[subscripts])
         #)
       }
)







######################
# 24.6.5 Model validation

# Fitted values Poisson GAM + Besag model
mu2      <- M2$summary.fitted.values[,"mean"]
E2       <- (Torn$nT - mu2) / sqrt(mu2)

# Figure 24.10
par(mfrow = c(2,2), mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = mu2, 
     y = E2,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)

plot(x = mu2, 
     y = Torn$nT,
     xlab = "Fitted values",
     ylab = "Observed nT ",
     xlim = c(0, 12),
     ylim = c(0, 12))


plot(x = Torn$Area.std, 
     y = E2,
     xlab = "Area",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)


plot(x = Torn$LogDensity, 
     y = E2,
     xlab = "Log density",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)
############################




# Simulation study
M2.sim <- inla(formula = f2, 
           family = "poisson", 
           control.compute = list(config = TRUE),
           data = Torn)

NSim <- 1000
SimData <- inla.posterior.sample(n = NSim, result = M2.sim)

MyParams <- rownames(M2.sim$summary.fixed)
MyID <- function(x){ which(rownames(SimData[[1]]$latent) == x) }
RowNum.Betas <- lapply(MyParams, MyID)
RowNum.Betas <- as.numeric(RowNum.Betas)
# The MyID function does not work as the parameters have :1 suffix but
# they always weem to be rows 4489 to 4497 (same as page 662)
RowNum.Betas <- 4489:4497
RowNum.Betas

N1 <- 102
MyParamsID <- paste("ID:", 1:N1, sep = "")
RowNum.ID <- lapply(MyParamsID, MyID)
RowNum.ID <- unlist(RowNum.ID)
RowNum.ID

ui <- M2$summary.random$ID$mean

df.102 <- data.frame(ID = Illi.shp$ID,
                     ui = Ui)
df.4386 <- data.frame(ID = Torn$ID)
merge(df.102, df.4386, by = "ID", sort = FALSE)[,"ui"]


library(plyr)
Join1 <- join(df.102, df.4386, by = "ID", type = "inner")["ui"]
Join2 <- df.102[df.4386[,"ID"],"ui"]

Join1-Join2
# Start a loop to extract betas and ws for each component, 
# calculate the fitted values mu and pi,
# and simulate count data from the model.
N  <- nrow(Torn)
Ysim <- matrix(nrow = N, ncol = NSim)
SSsim <- vector(length = N)

X <- model.matrix(~ 1 + ElevS.std + LogDensity.std + Area.std + Year.cr,
                     data = Torn)
X <- as.matrix(X)


for (i in 1: NSim){
  Betas <- SimData[[i]]$latent[RowNum.Betas]
  ui    <- SimData[[i]]$latent[RowNum.ID] 
  
  df.102 <- data.frame(ID = Illi.shp$ID, ui = ui)
  ui.Big <- join(df.102, df.4386, by = "ID", type = "inner")[,'ui']

  mu        <- exp(X %*% Betas + ui.Big)
  Ysim[,i]  <- rpois(N, lambda = mu)
  E         <- (Ysim[,i] - mu) / sqrt(mu) 
  SSsim[i]  <- sum(E^2) / N
  }


# Now we have 1000 simulated data sets from the model.
# What shall we do with these simulated data sets?
# We could calculate the number of zeros in each of the 1,000
# data sets.
zeros <- vector(length = NSim)
for(i in 1:NSim){
  zeros[i] <- sum(Ysim[,i] == 0)
}


#Let's plot this as a table

# Figure 24.11
par(mfrow = c(2,1), mar = c(5,5,2,2), cex.lab = 1.5)

Z <- table(zeros)
Range <- range(as.numeric(names(Z)))
SumZeros <- sum(Torn$nT == 0)
x1 <- min(Range[1], SumZeros)
x2 <- max(Range[2], SumZeros)

plot(table(zeros), 
     xlab = "How often do we have 0, 1, 2, 3, etc. number of zeros",
     ylab = "Frequency",
     xlim = c(0.98 * x1, 1.02 * x2),
     main = "")
points(x = sum(Torn$nT == 0), 
       y = 0, 
       pch = 16, 
       cex = 3, 
       col = 2)
#The red dot is the number of zeros in the original data set.
#The data simulated from the Poisson model
#         does not contain enough zeros.


Range <- range(SSsim)
SSE <- sum(E2^2) / N
x1 <- min(Range[1], SSE)
x2 <- max(Range[2], SSE)

hist(SSsim, 
     xlim = c(0.98 * x1, 1.02 * x2),
     main = "",
     xlab = "Sum of squared Pearson residuals / N")

points(x = sum(E2^2) / N, 
       y = 0, 
       pch = 16, 
       cex = 3, 
       col = 2)

#####################################################






# 24.7 Zero-inflated GAM with iCAR correlation
M3 <- inla(formula = f2, 
           lincomb = lcs.Year,
           family = "zeroinflatedpoisson1", 
           data = Torn,
           control.compute = list(dic = TRUE, waic = TRUE))



# And compare the models with DICs and WAICs
dic  <- c(M1$dic$dic, M2$dic$dic,  M3$dic$dic)   
waic <- c(M1$waic$waic, M2$waic$waic, M3$waic$waic)
Z.out     <- cbind(dic, waic)
rownames(Z.out) <- c("Poisson GAM",  
                     "Poisson GAM + Besag",
                     "ZIP GAM + Besag")
Z.out



# Not in the book:
M4 <- inla(formula = f2, 
           family = "nbinomial", 
           data = Torn,
           quantiles = c(.05, .5, .95),
           control.compute = list(dic = TRUE, waic = TRUE))


# And compare the models with DICs and WAICs
dic  <- c(M1$dic$dic, M2$dic$dic,  M3$dic$dic, M4$dic$dic)   
waic <- c(M1$waic$waic, M2$waic$waic, M3$waic$waic, M4$waic$waic)
Z.out     <- cbind(dic, waic)
rownames(Z.out) <- c("Poisson GAM",  
                     "Poisson GAM + Besag",
                     "ZIP GAM + Besag",
                     "NB GAM + Besag")
Z.out






# 24.8 Changing the priors and the BYM2 model
f2 <- nT ~ ElevS.std + LogDensity.std + Area.std +
  Year.cr +
  f(ID, 
    model = "besag", 
    graph = Illi.adj,
    scale.model=TRUE) 

M2 <- inla(formula = f2, 
           lincomb = lcs.Year,
           family = "poisson", 
           data = Torn,
           control.compute = list(dic = TRUE, waic = TRUE))

summary(M2)
M2$summary.hyperpar

tau_CAR <- M2$marginals.hyperpar$`Precision for ID`
sigma_CAR <- inla.emarginal(function(x) (1/sqrt(x)), tau_CAR)
sigma_CAR



############################################
# 24.8.2 Penalised complexity prior

HyperPC = list(prec = list(prior = "pc.prec",
                           fixed = FALSE,
                           param = c(0.1, 0.0001))) # originally commented out; output matches book
                           #param = c(1, 0.5))) # original in code sigma.pc_CAR = 0.368
                           #param = c(0.1, 0.001))) # text in book but sigma.pc_CAR = 0.226

f2.pc <- nT ~ ElevS.std + LogDensity.std + Area.std + Year.cr +
         f(ID, 
           model = "besag", 
           graph = Illi.adj,
           hyper = HyperPC, 
           adjust.for.con.comp = FALSE,
           constr = TRUE,
           scale.model=TRUE)

M2.pc <- inla(formula = f2.pc, 
           lincomb = lcs.Year,
           family = "poisson", 
           data = Torn,
           control.compute = list(dic = TRUE, waic = TRUE))

tau.pc_CAR <- M2.pc$marginals.hyperpar$`Precision for ID`
sigma.pc_CAR <- inla.emarginal(function(x) (1/sqrt(x)), tau.pc_CAR)
sigma.pc_CAR
#################################################





# 24.8.3 The BYM model for spatial correlation
#HyperBYM = list(prec = list(prior = "pc.prec",
#                            param = c(1, 0.5)))

HyperBYM <- list(
  prec.unstruct=list(prior = "pc.prec",param = c(1, 0.5)),
  prec.spatial=list(prior = "pc.prec",param = c(1, 0.5)))



f4 <- nT ~ ElevS.std + LogDensity.std + Area.std + Year.cr +
  f(ID, 
    model = "bym", 
    graph = Illi.adj,
    hyper = HyperBYM, 
    adjust.for.con.comp = FALSE,
    constr = TRUE,
    scale.model=TRUE)

M4 <- inla(formula = f4, 
              lincomb = lcs.Year,
              family = "poisson", 
              data = Torn,
              control.compute = list(dic = TRUE, waic = TRUE))
summary(M4)
tau.pc_CAR <- M4$marginals.hyperpar$`Precision for ID (spatial component)`
sigma.pc_CAR <- inla.emarginal(function(x) (1/sqrt(x)), tau.pc_CAR)
sigma.pc_CAR


tau.iid <- M4$marginals.hyperpar$`Precision for ID (iid component)`
sigma.iid <- inla.emarginal(function(x) (1/sqrt(x)), tau.iid)
sigma.iid



# And compare the models with DICs and WAICs
dic  <- c(M1$dic$dic, M2$dic$dic,  M4$dic$dic)   
waic <- c(M1$waic$waic, M2$waic$waic,  M4$dic$dic)
Z.out     <- cbind(dic, waic)
rownames(Z.out) <- c("Poisson GAM",  
                     "Poisson GAM + Besag",
                     "Poisson GAM + BYM")
Z.out

M4$summary.random$ID


Torn$County  <- Torn$ID
Torn$County2 <- Torn$County
Torn$ID2     <- Torn$ID



# 28.8.4 The BYM2 model for spatial correlation

# Can't remember what this one was. # This comment is written by Alain Zuur!!
# Roy has commented f5 and M5 out as not in book
# f5 <- nT ~ ElevS.std + LogDensity.std + Area.std + Year.cr +
#   f(County, 
#     model = "besag", 
#     graph = Illi.adj,
#     adjust.for.con.comp = FALSE,
#     constr = TRUE,
#     scale.model=TRUE) +
#   f(County2, model = "iid") 
#   
# M5 <- inla(formula = f5, 
#            lincomb = lcs.Year,
#            family = "poisson", 
#            data = Torn,
#            control.compute = list(dic = TRUE, waic = TRUE))



# BYM2 model
HyperBYM2 <- list(
  prec = list(prior = "pc.prec", param = c(1 , 0.01)),
  phi  = list(prior = "pc", param = c(0.5, 0.5)))


f6 <- nT ~ ElevS.std + LogDensity.std + Area.std + Year.cr +
  f(ID, 
    model = "bym2", 
    graph = Illi.adj,
    hyper = HyperBYM2, 
    adjust.for.con.comp = TRUE,
    constr = TRUE,
    scale.model=TRUE)

M6 <- inla(formula = f6,      # Very slow; about 5 minutes
           lincomb = lcs.Year,
           family = "poisson", 
           data = Torn,
           control.compute = list(dic = TRUE, waic = TRUE))

summary(M6)
tau <- M6$marginals.hyperpar$`Precision for ID`
sigma <- inla.emarginal(function(x) (1/sqrt(x)), tau)
sigma




#######################################################
# 24.9 Spatial-temporal correlation for areal data
Torn$County <- Torn$ID
Torn$County1 <- Torn$ID
Torn$Year0 <- Torn$Year - min(Torn$Year) + 1
Torn$Year1 <- Torn$Year - min(Torn$Year) + 1

Hyper <- list(theta = list(prior = "pc.prec", 
                           param = c(0.01, 0.0001)))

f7 <- nT ~ ElevS.std + LogDensity.std + Area.std + 
           f(County, 
             model = "bym", 
             graph = Illi.adj,
             adjust.for.con.comp = FALSE,
             constr = TRUE,
             scale.model = TRUE) +
           f(Year, 
             model = "rw2") +
           f(Year1,
             model = "iid")

lcs <- inla.make.lincombs(Year = diag(43),  Year1 = diag(43))


M7 <- inla(formula = f7, 
           family = "poisson", 
           lincomb=lcs,
           data = Torn,
           control.compute = list(dic = TRUE, waic = TRUE))

f.Year    <- M7$summary.random$Year[, "mean"] 
SeLo.Year <- M7$summary.random$Year[,"0.025quant"] 
SeUp.Year <- M7$summary.random$Year[,"0.975quant"]

f.Yeariid    <- M7$summary.random$Year1[, "mean"] 
SeLo.Yeariid <- M7$summary.random$Year1[,"0.025quant"] 
SeUp.Yeariid <- M7$summary.random$Year1[,"0.975quant"]


# Combine the smoother, and 95% CIs
MyData <- data.frame(mu   = c(f.Year),
                     SeUp = c(SeUp.Year),
                     SeLo = c(SeLo.Year),
                     Xaxis = c(M7$summary.random$Year[,1]))

MyData.iid <- data.frame(mu   = c(f.Yeariid),
                        SeUp = c(SeUp.Yeariid),
                        SeLo = c(SeLo.Yeariid),
                        Xaxis = c(M7$summary.random$Year1[,1]+1969))

# Figure 24.13
p <- ggplot()
p <- p + xlab("Time") + ylab("Temporal trend and iid noise")
p <- p + theme(text = element_text(size = 15))
p <- p + geom_line(data = MyData, 
                   aes(x = Xaxis, y = mu),
                   lwd = 3)

p <- p + geom_ribbon(data = MyData, 
                     aes(x = Xaxis, 
                         ymax = SeUp, 
                         ymin = SeLo),
                     alpha = 0.6)

p <- p + geom_line(data = MyData.iid, 
                   aes(x = Xaxis, y = mu))
p





################################################
# 24.9.3 Type I interaction
Torn$County.Year <- as.numeric(as.factor(paste(Torn$County, Torn$Year, sep = ".")))

f8 <- nT ~ ElevS.std + LogDensity.std + Area.std + 
  f(County, 
    model = "bym", 
    graph = Illi.adj,
    adjust.for.con.comp = FALSE,
    constr = TRUE,
    scale.model = TRUE) +
  f(Year, 
    model = "rw2") +
  f(Year1,
    model = "iid") +
  f(County.Year, model = "iid")



M8 <- inla(formula = f8, 
           family = "poisson", 
           data = Torn,
           control.compute = list(dic = TRUE, waic = TRUE))
f.Year    <- M8$summary.random$Year1[, "mean"] 
SeLo.Year <- M8$summary.random$Year1[,"0.025quant"] 
SeUp.Year <- M8$summary.random$Year1[,"0.975quant"]


# Not ion the book
# Combine the smoother, and 95% CIs
MyData <- data.frame(mu   = c(f.Year),
                     SeUp = c(SeUp.Year),
                     SeLo = c(SeLo.Year),
                     Xaxis = c(M7$summary.random$Year1[,1]))

p <- ggplot()
p <- p + xlab("Time") + 
  ylab("Smoother")
p <- p + theme(text = element_text(size = 15))
p <- p + geom_line(data = MyData, 
                   aes(x = Xaxis, y = mu))

p <- p + geom_ribbon(data = MyData, 
                     aes(x = Xaxis, 
                         ymax = SeUp, 
                         ymin = SeLo),
                     alpha = 0.6)
p
summary(M8)



################################################
# 24.9.4 Type II interaction
f9 <- nT ~ ElevS.std + LogDensity.std + Area.std + 
  f(County, 
    model = "bym", 
    graph = Illi.adj,
    adjust.for.con.comp = FALSE,
    constr = TRUE,
    scale.model = TRUE) +
  f(Year0, 
    model = "rw2") +
  f(Year1,
    model = "iid") +
  f(County1, 
    model = "iid",
    group = Year1,
    control.group = list(model = "rw2"))

M9 <- inla(formula = f9, 
           family = "poisson", 
           data = Torn,
           control.compute = list(dic = TRUE, waic = TRUE))

round(M9$summary.hyperpar[,c("mean", "mode")], 2)

tau <- M9$marginals.hyperpar$`Precision for County1`
sigma <- inla.emarginal(function(x) (1/sqrt(x)), tau)
sigma



f.Year    <- M9$summary.random$Year1[, "mean"] 
SeLo.Year <- M9$summary.random$Year1[,"0.025quant"] 
SeUp.Year <- M9$summary.random$Year1[,"0.975quant"]


names(M9$summary.random)
nrow(M9$summary.random$Year0) 

Time <- M9$summary.random$Year0[,"ID"] + 1969
MainTr <- M9$summary.random$Year0[,"mean"]
TimeComplete <- rep(1970:2013, each = 102)

ST.II <- M9$summary.random$County1
ST.II$TimeComplete <- TimeComplete

ST.IIa <- subset(ST.II, TimeComplete != 2007)

ST.IIa$Time <- rep(Time, each = 102)
ST.IIa$MainTrend <- rep(MainTr, each = 102)
RE.year <- M9$summary.random$Year0[,"mean"]
ST.IIa$RE.year <- rep(RE.year, each = 102)

head(ST.IIa)

# Interaction trend area 1

library(ggplot2)
# Figure 24.16 and 24.17 (drop the exp for Figure 24.16)
p <- ggplot()
p <- p + xlab("Time") + ylab("Time + space-time interaction for each county")
p <- p + theme(text = element_text(size = 15), legend.position="none")
p <- p + geom_line(data = ST.IIa, col =1,
                   aes(x = Time, y = exp(mean + MainTrend + RE.year), group = ID))
p


# Not in the book:
delta.intII <- data.frame(delta = ST.IIa$mean,
                           year  = Torn$Year,
                           ID.area = Torn$County)
head(delta.intII)

delta.intII.matrix <- matrix(delta.intII[,1], 
                              nrow =102,
                              ncol = 43,
                              byrow = TRUE)
delta.intII.matrix[1:5, 1:5]
rownames(delta.intII.matrix)<- 1:102


cutoff.interaction <- c(-2,-1, 1,2)

delta.intII.factor <- data.frame(NAME=1:102)
for(i in 1:43){
  delta.factor.temp <- cut(delta.intII.matrix[,i],
                           breaks = cutoff.interaction,
                           include.lowest=TRUE) 
  delta.intII.factor <- cbind(delta.intII.factor,delta.factor.temp)
}
colnames(delta.intII.factor)<- c("NAME",Time)


attr(Illi.shp, "data") <- data.frame(Torn, intII = delta.intII.factor)
library(lattice)
trellis.par.set(axis.line=list(col=NA))


spplot(obj=Illi.shp, zcol=c("intII.1973", "intII.1974", "intII.1975",
                            "intII.1976", "intII.1977", "intII.1978",
                            "intII.1979", "intII.1980", "intII.1981",
                            "intII.1982", "intII.1983", "intII.1984",
                            "intII.1985", "intII.1986", "intII.1987",
                            "intII.1988", "intII.1989", "intII.1990",
                            "intII.1991", "intII.1992", "intII.1993",
                            "intII.1994", "intII.1995", "intII.1996",
                            "intII.1997", "intII.1998", "intII.1999",
                            "intII.2000", "intII.2001", "intII.2002",
                            "intII.2003", "intII.2004", "intII.2005",
                            "intII.2006", "intII.2008", "intII.2009",
                            "intII.2010", "intII.2011", "intII.2012",
                            "intII.2013"), 
       col.regions=c("blue", "white", "red"), #gray(2.5:0.5/3),
       names.attr=c(seq(1973, 2006),seq(2008, 2013)),#Time,
       main="")












###############
# 24.9.5 Type III interaction
Torn$Year2 <- Torn$Year - min(Torn$Year) + 1

f10 <- nT ~ ElevS.std + LogDensity.std + Area.std + 
  f(County, 
    model = "bym", 
    graph = Illi.adj,
    adjust.for.con.comp = FALSE,
    constr = TRUE,
    scale.model = TRUE) +
  f(Year0, 
    model = "rw2") +
  f(Year1,
    model = "iid") +
  f(Year2, 
    model = "iid",
    group = County1,
    control.group = list(model = "besag", graph = Illi.adj))

M10 <- inla(formula = f10, # Very slow to converge (10 to 15 minutes)
           family = "poisson", 
           data = Torn,
           control.compute = list(dic = TRUE, waic = TRUE))




delta.intIII <- data.frame(delta = M10$summary.random$Year2[,2],
                           year  = Torn$Year,
                           ID.area = Torn$County)

head(delta.intIII)

delta.intIII.matrix <- matrix(delta.intIII[,1], 
                              nrow =102,
                              ncol = 43,
                              byrow = TRUE)
delta.intIII.matrix[1:5, 1:5]
rownames(delta.intIII.matrix)<- 1:102


cutoff.interaction <- c(-2,-1, 1,2)

delta.intIII.factor <- data.frame(NAME=1:102)
for(i in 1:43){
  delta.factor.temp <- cut(delta.intIII.matrix[,i],
                          breaks = cutoff.interaction,
                          include.lowest=TRUE) 
  delta.intIII.factor <- cbind(delta.intIII.factor,delta.factor.temp)
}
colnames(delta.intIII.factor)<- c("NAME",Time)


attr(Illi.shp, "data") <- data.frame(Torn, intIII = delta.intIII.factor)
library(lattice)
trellis.par.set(axis.line=list(col=NA))

#Figure 24.19
spplot(obj=Illi.shp, zcol=c("intIII.1973", "intIII.1974", "intIII.1975",
                            "intIII.1976", "intIII.1977", "intIII.1978",
                            "intIII.1979", "intIII.1980", "intIII.1981",
                            "intIII.1982", "intIII.1983", "intIII.1984",
                            "intIII.1985", "intIII.1986", "intIII.1987",
                            "intIII.1988", "intIII.1989", "intIII.1990",
                            "intIII.1991", "intIII.1992", "intIII.1993",
                            "intIII.1994", "intIII.1995", "intIII.1996",
                            "intIII.1997", "intIII.1998", "intIII.1999",
                            "intIII.2000", "intIII.2001", "intIII.2002",
                            "intIII.2003", "intIII.2004", "intIII.2005",
                            "intIII.2006", "intIII.2008", "intIII.2009",
                            "intIII.2010", "intIII.2011", "intIII.2012",
                            "intIII.2013"), 
       col.regions=c("blue", "white", "red"), #gray(2.5:0.5/3),
       names.attr=c(seq(1973, 2006),seq(2008, 2013)),#Time,
       main="")
######################




####################
# 24.9.6 Type IV interaction
f11 <- nT ~ ElevS.std + LogDensity.std + Area.std + 
  f(County, 
    model = "bym", 
    graph = Illi.adj,
    adjust.for.con.comp = FALSE,
    constr = TRUE,
    scale.model = TRUE) +
  f(Year0, 
    model = "rw2") +
  f(Year1,
    model = "iid") +
  f(County1, 
    model = "besag", 
    graph = Illi.adj,
    group = Year1,
    control.group = list(model = "rw2"))

M11 <- inla(formula = f11,  # Very slow. 10 to 15 mins CPU
            family = "poisson", 
            data = Torn,
            control.compute = list(dic = TRUE, waic = TRUE))








############################################
# Figure a la Figure 1 in Knorr-Held (1999)
#Figures 24.12, 24.14, 24.15, 24.18, 
par(mfrow=c(1,1)) # Too cramp to read in RStudio if on one plot
plot(x = 0, y = 0, axes = FALSE, type = "n",
     xlim = c(0, 11),
     ylim = c(0,11),
     xlab = "",
     ylab = "")

lines(x=c(3,3),y=c(1,10))
lines(x=c(1,10),y=c(3,3))
xy.Ind <- expand.grid(4:9, seq(4,9, by = 0.5))
points(xy.Ind, pch = 16)
text(x=10, y = 3.5, "Time", cex = 1.5)
text(x=3, y = 10.5, "Space", cex = 1.5)

# At the bottom:
points(x = 4:9, y = rep(1,6), pch = 1, cex = 2.5)
text(10, 1, expression(italic(phi[s])), cex = 1.5)
text(10, 2, expression(italic(gamma[s])), cex = 1.5)

rect(xleft = 3.75, ybottom = 1.75, xright = 9.25, ytop = 2.25,
     col = grey(0.5))

# left
points(y = seq(4,9, by = 0.5), x = rep(1,11), pch = 1, cex = 2)
text(1, 10, expression(italic(v[i])), cex = 1.5, font = 3)
text(2, 10, expression(italic(u[i])), cex = 1.5)

rect(xleft = 1.75, ybottom = 3.75, ytop = 9.25, xright = 2.25,
     col = grey(0.5))


# Add this for type I interaction
# points(xy.Ind, pch = 1, cex = 2)

# Add this for type II interaction
#for (i in 1:11){ 
#  rect(xleft = 3.5, ybottom = 3.9 +0.5*(i-1), ytop = 4.1+0.5*(i-1), xright = 9.5,
#       col = grey(0.7))
#}
#points(xy.Ind, pch = 16)



# Add this for type III interaction
#for (i in 1:6){ 
#  rect(xleft = 3.9 +(i-1), 
#       ybottom = 3.9, 
#       ytop = 9.5, 
#       xright = 4.1 + (i-1),
#       col = grey(0.7))
#}
points(xy.Ind, pch = 16)


