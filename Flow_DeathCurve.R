# the script is designed to analyse real-time changes in PI-positiveness
#  of population
# Library loading ---------------------------------------------------------

library("flowViz")
library("flowStats")

# Cleaning enviroment -----------------------------------------------------

rm(list=ls())

# List of transformation --------------------------------------------------

log.Tr <- logTransform(transformationId="log10-transformation", logbase=10, r=1, d=1)
tr.Tr<-truncateTransform("myTrTr", a=1)
time.Tr <- linearTransform(transformationId = "Timetransformation", a = 1/10/60, b = 0)

# Define sample names, based on input data --------------------------------

FSC.list<-list.files()[grep("fcs",list.files())]
Quant<-c(1:length(FSC.list))
ArRaY<-read.flowSet(FSC.list, transformation=FALSE, alter.names=T)
SampleDef <- data.frame(Quant, Quant)
for (i in Quant) SampleDef[i,1] <- description(ArRaY[[i]])$`FIL`
for (i in Quant) SampleDef[i,2] <- description(ArRaY[[i]])$`#SAMPLE`
print(SampleDef)
# remove flowSet, based on old FSC.list and write the 
# definition of sample to the file with following reading it
remove(ArRaY)
colnames(SampleDef)<-c("File","Treatment")
write.table(x = SampleDef, file = paste(getwd(), "/SampleList.csv", sep=""),row.names = F)
SampleDef<-read.table(file = paste(getwd(), "/SampleList.csv", sep=""), header = T,
                      stringsAsFactors = FALSE)

# Creating flowSet object and its transformation --------------------------

Quant<-c(1:length(FSC.list))
# define the matrix layout for plotting all sample in one device
Quant2<-Quant

ArRaY<-read.flowSet(FSC.list, transformation=FALSE, alter.names=T)
ArRaY.tr1<-transform(ArRaY, `FL3.A`=tr.Tr(`FL3.A`))
ARray<-transform(ArRaY.tr1, `FL3.A`=log.Tr(`FL3.A`))
remove(ArRaY.tr1)
remove(ArRaY)
LiveDead<-ARray[seq(along=ARray)]

#---- Naming of samples----

# Names - names of files with gates,
# Samples - names of samples in the plots

Names<-as(phenoData(ARray), "data.frame")
Names<-Names[,1]
Names <- gsub(pattern = " ", replacement = "_", x = Names)
NamesP1 <- gsub(pattern = ".fcs", replacement = "P1", x = Names)
NamesP2 <- gsub(pattern = ".fcs", replacement = "P2", x = Names)
NamesP3 <- gsub(pattern = ".fcs", replacement = "P3", x = Names)

Samples<-SampleDef[,2]

# Introduction of "internal" values ---------------------------------------

PercentageP1<-Quant
PercentageP2<-Quant
PercentageP3<-Quant
W <- 6 # width of device
H <- 6 # height of device
NRs <- 2 # numbers of rows in the matrix layout


# Defining the P1 gate: "Cells and cell derived debris" -------------------

# FSC vs SSC plots to define limits
dev.new(width = W, heigth = H)
layout(matrix(Quant, nrow = NRs, byrow = TRUE), respect = TRUE)
for (i in Quant) { 
  plot(ARray[[i]], c("FSC.A", "SSC.A"),
       colramp = colorRampPalette(c("white", "gray90", "gray70","gray50","gray30","gray20", "gray10")), 
       bandwidth=1048, nbin = 512, main = paste(Samples[i]))
}
locator () # delay to think and place limits

FSCLim <- c(0,9000000)
SSCLim <- c(0,3500000)

# show all FSC-SSC plot with limits to oversee and choose one for gating

dev.new(width = W, heigth = H)
layout(matrix(Quant2, nrow = NRs, byrow = TRUE), respect = TRUE)
for (i in Quant) { 
  plot(ARray[[i]], c("FSC.A", "SSC.A"), xlim = FSCLim, ylim = SSCLim,
       colramp = colorRampPalette(c("white", "gray90", "gray70","gray50","gray30","gray20", "gray10")), 
       bandwidth=1048, nbin = 512, main = paste(Samples[i]))
}

# Define P1 gate: Live + Dead cells

# the cycle looks if there is gate description for this sample in the working directory
# and if "yes" then apply it, if not cycles looks for gate information for previous
# sample and if "yes" then apply it, if both are "no" then plot shows random define gate


for (i in Quant) {
  if (i>1) {      
    if (system(paste("ls | grep ", NamesP1[i], ".txt", sep=""))==0) {
      EST.P1.points<-read.table(file=paste(NamesP1[i], ".txt", sep=""), sep="\t", header=T)
    } else {      
      if (system(paste("ls | grep ", NamesP1[i-1], ".txt", sep=""))==0) {
        EST.P1.points<-read.table(file=paste(NamesP1[i-1], ".txt", sep=""), sep="\t", header=T)
      } else {
        EST.P1.points <- matrix(c(1469297,  1538984,
                                  233817,   186899,
                                  674293, 19530,
                                  5283172, 12966,
                                  5594727,  193463,
                                  5519524, 1516012), ncol = 2, byrow = T)
      }
    }
  } else {
    if (system(paste("ls | grep ", NamesP1[i], ".txt", sep=""))==0) {
      EST.P1.points<-read.table(file=paste(NamesP1[i], ".txt", sep=""), sep="\t", header=T)
    } else {
      EST.P1.points <- matrix(c(1469297,  1538984,
                                233817,   186899,
                                674293, 19530,
                                5283172, 12966,
                                5594727,  193463,
                                5519524, 1516012), ncol = 2, byrow = T)
    }
  }
  colnames(EST.P1.points) <- c("FSC.A","SSC.A")
  EST.P1.points<-tr.Tr(EST.P1.points)
  EST.P1.filter <- polygonGate( .gate = EST.P1.points)
  EST.P1 <- filter(ARray[[i]], EST.P1.filter)
  
  dev.new(title = paste(Samples[i]))
  flowPlot(ARray[[i]], EST.P1, plotParameters = c("FSC.A", "SSC.A"),
           xlim = FSCLim, ylim = SSCLim, colParent = "gray40",
           colChild = "gray40",pch = 15, cex = 0.2, main = "choose the gate (P1)")
  polygon(EST.P1.points[,1], EST.P1.points[,2], border = "red", lwd = 2)
  a <- locator(type = "l")
  ifelse (is.null(a), {
    write.table(file = paste(NamesP1[i], ".txt", sep=""), EST.P1.points, col.names=T,  sep = "\t")
    P1.points <- read.table(file = paste(NamesP1[i], ".txt", sep=""), sep="\t")
  }, {
    write.table(file = paste(NamesP1[i], ".txt", sep=""), a, sep = "\t")
    P1.points <- read.table(file = paste(NamesP1[i], ".txt", sep=""), sep="\t")
  })
  
#   after confirmation of the gate, the corresponding gate information is written in 
#   working directory in file Sample_XXXPX.txt. Then it used and for "final" plotting,
#   showing percentage, and subsetting the flowFrame by this gate. The separate flowSet 
#   is created.
   
  colnames(P1.points) <- c("FSC.A","SSC.A")
  P1.points <- tr.Tr(P1.points)
  P1.filter <- polygonGate(.gate = P1.points)
  P1 <- filter(ARray[[i]], P1.filter)
  PercentageP1[i]<-round(summary(P1)$p*100,2)
  flowPlot(ARray[[i]], P1, plotParameters = c("FSC.A", "SSC.A"),
           xlim = FSCLim, ylim = SSCLim, colParent = "gray80",
           colChild = "gray40",pch = 15, cex = 0.2, main=Samples[i])
  text(x=350000, y=1800000, labels = paste(PercentageP1[i],"%", sep = ""), cex = 1.2, col = "red")
  polygon(P1.points[,1], P1.points[,2], border = "red", lwd = 2)
  LiveDead[[i]]<-Subset(ARray[[i]], P1)
  locator()
}

# Block to close open devices ---------------------------------------------
source("~/Dropbox/R-dev/DeviseCloser.R")

# Defining gate P2: removing the dublets ----------------------------------

# Set up the limits for width and check how it looks

WidthLim <- c(0,400)

dev.new(width = W, heigth = H)
layout(matrix(Quant2, nrow = NRs, byrow = TRUE), respect = TRUE)
for (i in Quant) { 
  plot(LiveDead[[i]], c("FSC.A", "Width"), bandwidth=0.01, ylim = WidthLim, xlim = FSCLim,
       colramp = colorRampPalette(c("white", "gray90", "gray70","gray50","gray30","gray20", "gray10")), 
       nbin = 512, main = paste(Samples[i]))
}

# Defining P2 gate on FSCvsWidth: removing dublets 

NoDublets<-LiveDead[seq(along=LiveDead)]

# the cycle looks if there is gate description for this sample in the working directory
# and if "yes" then apply it, if not cycles looks for gate information for previous
# sample and if "yes" then apply it, if both are "no" then plot shows random define gate

for (i in Quant) {
  if (i>1) {      
    if (system(paste("ls | grep ", NamesP2[i], ".txt", sep=""))==0) {
      EST.P2.points<-read.table(file=paste(NamesP2[i], ".txt", sep=""), sep="\t", header=T)
    } else {      
      if (system(paste("ls | grep ", NamesP2[i-1], ".txt", sep=""))==0) {
        EST.P2.points<-read.table(file=paste(NamesP2[i-1], ".txt", sep=""), sep="\t", header=T)
      } else {
        EST.P2.points <- matrix(c(1469297,  1538984,
                                  233817,   186899,
                                  674293, 19530,
                                  5283172, 12966,
                                  5594727,  193463,
                                  5519524, 1516012), ncol = 2, byrow = T)
      }
    }
  } else {
    if (system(paste("ls | grep ", NamesP2[i], ".txt", sep=""))==0) {
      EST.P2.points<-read.table(file=paste(NamesP2[i], ".txt", sep=""), sep="\t", header=T)
    } else {
      EST.P2.points <- matrix(c(1469297,  1538984,
                                233817,   186899,
                                674293, 19530,
                                5283172, 12966,
                                5594727,  193463,
                                5519524, 1516012), ncol = 2, byrow = T)
    }
  }
  colnames(EST.P2.points) <- c("FSC.A","Width")
  EST.P2.points<-tr.Tr(EST.P2.points)
  EST.P2.filter <- polygonGate( .gate = EST.P2.points)
  EST.P2 <- filter(LiveDead[[i]], EST.P2.filter)
  
#   after confirmation of the gate, the corresponding gate information is written in 
#   working directory in file Sample_XXXPX.txt. Then it used and for "final" plotting,
#   showing percentage, and subsetting the flowFrame by this gate. No separate flowSet 
#   is created.
  
  dev.new(title = paste(Samples[i]))
  flowPlot(LiveDead[[i]], EST.P2, plotParameters = c("FSC.A", "Width"),
           xlim = FSCLim, ylim = WidthLim, colParent = "gray40",
           colChild = "gray40",pch = 15, cex = 0.2, main = "choose the gate (P1)")
  polygon(EST.P2.points[,1], EST.P2.points[,2], border = "red", lwd = 2)
  a <- locator(type = "l")
  ifelse (is.null(a), {
    write.table(file = paste(NamesP2[i], ".txt", sep=""), EST.P2.points, col.names=T,  sep = "\t")
    P2.points <- read.table(file = paste(NamesP2[i], ".txt", sep=""), sep="\t")
  }, {
    write.table(file = paste(NamesP2[i], ".txt", sep=""), a, sep = "\t")
    P2.points <- read.table(file = paste(NamesP2[i], ".txt", sep=""), sep="\t")
  })
  colnames(P2.points) <- c("FSC.A","Width")
  P2.points <- tr.Tr(P2.points)
  P2.filter <- polygonGate(.gate = P2.points)
  P2 <- filter(LiveDead[[i]], P2.filter)
  PercentageP2[i]<-round(summary(P2)$p*100,2)
  flowPlot(LiveDead[[i]], P2, plotParameters = c("FSC.A", "Width"),
           xlim = FSCLim, ylim = WidthLim, colParent = "gray80",
           colChild = "gray40",pch = 15, cex = 0.2, main=Samples[i])
  text(x=350000, y=180, labels = paste(PercentageP2[i],"%", sep = ""), cex = 1.2, col = "red")
  polygon(P2.points[,1], P2.points[,2], border = "red", lwd = 2)
  NoDublets[[i]]<-Subset(LiveDead[[i]], P2)
  locator()
}


# Block to close open devices ---------------------------------------------
source("~/Dropbox/R-dev/DeviseCloser.R")

# Defining P3 gate: removing "crap" ---------------------

# Creating the flowSet for this gate

NoCrap<-NoDublets[seq(along=NoDublets)]

# calculating crude limits for FL3.A

TempVec<-c(1:6)
for (i in Quant) {
  temps<-summary(NoDublets[[i, "FL3.A"]])[,]
  TempVec<-rbind(TempVec, temps)
}
FL3.SumUp<-TempVec[2:(nrow(TempVec)),]
rownames(FL3.SumUp)<-Names
FL3.SumUp<-as.data.frame(FL3.SumUp)
FL3.Max <- ceiling(max(FL3.SumUp$Max.))

TempVec<-Quant
for (i in Quant) {
  TempVec[i]<-quantile(exprs(NoDublets[[i, "FL3.A"]]), 0.01)
}
FL3.Min <- trunc(min(TempVec))

# Defining of ylim in FSC.A vs FL3.A plot

FL3Lim <- c(FL3.Min, FL3.Max)

dev.new(width = W, heigth = H)
layout(matrix(Quant2, nrow = NRs, byrow = TRUE), respect = TRUE)
for (i in Quant) { 
  plot(NoDublets[[i]], c("FSC.A", "FL3.A"), bandwidth=0.01, ylim = FL3Lim, xlim = FSCLim,
       colramp = colorRampPalette(c("white", "gray90", "gray70","gray50","gray30","gray20", "gray10")), 
       nbin = 512, main = paste(Samples[i]))
}

# Defining P3 gate on FSCvsFL3: removing "crap"

# the cycle looks if there is gate description for this sample in the working directory
# and if "yes" then apply it, if not cycles looks for gate information for previous
# sample and if "yes" then apply it, if both are "no" then plot shows random define gate

for (i in Quant) {
  if (i>1) {      
    if (system(paste("ls | grep ", NamesP3[i], ".txt", sep=""))==0) {
      EST.P3.points<-read.table(file=paste(NamesP3[i], ".txt", sep=""), sep="\t", header=T)
    } else {      
      if (system(paste("ls | grep ", NamesP3[i-1], ".txt", sep=""))==0) {
        EST.P3.points<-read.table(file=paste(NamesP3[i-1], ".txt", sep=""), sep="\t", header=T)
      } else {
        EST.P3.points <- matrix(c(1469297,  1538984,
                                  233817,   186899,
                                  674293, 19530,
                                  5283172, 12966,
                                  5594727,  193463,
                                  5519524, 1516012), ncol = 2, byrow = T)
      }
    }
  } else {
    if (system(paste("ls | grep ", NamesP3[i], ".txt", sep=""))==0) {
      EST.P3.points<-read.table(file=paste(NamesP3[i], ".txt", sep=""), sep="\t", header=T)
    } else {
      EST.P3.points <- matrix(c(1469297,  1538984,
                                233817,   186899,
                                674293, 19530,
                                5283172, 12966,
                                5594727,  193463,
                                5519524, 1516012), ncol = 2, byrow = T)
    }
  }
  colnames(EST.P3.points) <- c("FSC.A","FL3.A")
  EST.P3.points<-tr.Tr(EST.P3.points)
  EST.P3.filter <- polygonGate( .gate = EST.P3.points)
  EST.P3 <- filter(NoDublets[[i]], EST.P3.filter)
  
  #   after confirmation of the gate, the corresponding gate information is written in 
  #   working directory in file Sample_XXXPX.txt. Then it used and for "final" plotting,
  #   showing percentage, and subsetting the flowFrame by this gate. No separate flowSet 
  #   is created.
  
  dev.new(title = paste(Samples[i]))
  par(xaxs='i',yaxs='i')
  flowPlot(NoDublets[[i]], EST.P3, plotParameters = c("FSC.A", "FL3.A"),
           xlim = FSCLim, ylim = FL3Lim, colParent = "gray40",
           colChild = "gray40",pch = 15, cex = 0.2, main = "choose the gate (P1)")
  polygon(EST.P3.points[,1], EST.P3.points[,2], border = "red", lwd = 2)
  a <- locator(type = "l")
  ifelse (is.null(a), {
    write.table(file = paste(NamesP3[i], ".txt", sep=""), EST.P3.points, col.names=T,  sep = "\t")
    P3.points <- read.table(file = paste(NamesP3[i], ".txt", sep=""), sep="\t")
  }, {
    write.table(file = paste(NamesP3[i], ".txt", sep=""), a, sep = "\t")
    P3.points <- read.table(file = paste(NamesP3[i], ".txt", sep=""), sep="\t")
  })
  colnames(P3.points) <- c("FSC.A","FL3.A")
  P3.points <- tr.Tr(P3.points)
  P3.filter <- polygonGate(.gate = P3.points)
  P3 <- filter(NoDublets[[i]], P3.filter)
  PercentageP3[i]<-round(summary(P3)$p*100,2)
  flowPlot(NoDublets[[i]], P3, plotParameters = c("FSC.A", "FL3.A"),
           xlim = FSCLim, ylim = FL3Lim, colParent = "gray80",
           colChild = "gray40",pch = 15, cex = 0.2, main=Samples[i])
  text(x=350000, y=3.5, labels = paste(PercentageP3[i],"%", sep = ""), cex = 1.2, col = "red")
  polygon(P3.points[,1], P3.points[,2], border = "red", lwd = 2)
  NoCrap[[i]]<-Subset(NoDublets[[i]], P3)
  locator()
}


# Block to close open devices ---------------------------------------------
source("~/Dropbox/R-dev/DeviseCloser.R")

# Plotting death curve: FL3 vs Time ---------------------------------------

NoCrap<-transform(NoCrap, `Time` = time.Tr(`Time`))

# calculating crude limits for Time

TempVec<-c(1:6)
for (i in Quant) {
  temps<-summary(NoCrap[[i, "Time"]])[,]
  TempVec<-rbind(TempVec, temps)
}
Time.Sum<-TempVec[2:(nrow(TempVec)),]
rownames(Time.Sum)<-Names
Time.Sum<-as.data.frame(Time.Sum)
Time.Sum$Min. <- trunc(Time.Sum$Min.)
Time.Sum$Max. <- ceiling(Time.Sum$Max.)
Time.Sum <- Time.Sum[,c(1,6)]

# Creating the data frame for ticks and labels at plots

# this allows to plot all data uniformly starting with 2 min on x-axis
# that reflect the approximate time of starting aquisition
# all ticks are placed at 5-min intervals

a<-matrix(data = 1, nrow = length(Quant), ncol = 3)
a[,3]<-5
a<-as.data.frame(a)
for (i in Quant) a[i, 1:2] <- Time.Sum[i,]

Files <- gsub(pattern = ".fcs", replacement = "", x = Names)

for (i in Quant)
{ 
  
# plotting the Time vs FL3 plots: creating them and then writing 
# jpeg file on drive with showing it with simple viewer - ristretto
  
  png(file=paste( Files[i], ".png", sep=""), width = 5, height = 5, units = "in", res = 300)
  par(xaxs='i',yaxs='i')
  plot(NoCrap[[i]], c("Time", "FL3.A"), bandwidth=c(0.005, 0.0025), nbin = 768,
     xlab = "", ylab = "", xaxt = "n", yaxt = "n",
     colramp = colorRampPalette(c("white", "gray90","gray80","gray70","gray20", "gray15", "gray10")),
     xlim = as.numeric(Time.Sum[i,]), ylim = FL3Lim,  main = Samples[i], cex.main = 1.5)
  axis(side = 1, at = seq(a[i,1], a[i,2], a[i,3]), 
       labels = seq(2, (length(seq(a[i,1], a[i,2], a[i,3]))-1)*5+2, 5), lwd = 2,  cex.axis = 1.2)
  axis(side = 2, at = seq(1, 7, by = 1), lwd = 2, cex.axis = 1.2)
  title(xlab = "Time, min", ylab = expression("FL3-A," ~ log[10]), cex.lab = 1.5, line = 2.5)
  dev.off()
  system(paste("ristretto ", Files[i], ".png", sep=""))
}

