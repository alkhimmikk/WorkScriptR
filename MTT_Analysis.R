### - Libraries section - ####
library(ggplot2)
library(plyr)
library(reshape2)
library(plotrix)
library(grid)
library(extrafont)
library(dplyr)

### - Presets section - ####
rm(list=ls())


# input data - temporarily manually changed in script
Samples<-c("UT","ZAP","CT\nDNA","poly\nIC","DNA+\nDNAse")
SamOrder<-c("UT","ZAP","CT\nDNA","poly\nIC","DNA+\nDNAse")
Label<-c("30-06-14 MTT Assay on chicken BMM")
GenLabel<-c("30-06-14 MTT Assay on chicken BMM")
File<-"15-07-01_ChBMM_NA"
Hours<-c(1,3,6) # time points
Plicates <- 4 # number of well for one sample-treatment

#defining plot options
TxtSize<-element_text(size = rel(2.2), face=1)
LegTxtSize<-element_text(size = rel(28), angle=0, face=1)
LegTxtSize2<-element_text(size = rel(20), angle=0, face=1)
GridType<-element_blank()
AxLines<-element_line(size=1.5)
AR=2.75/length(SamOrder)

# creating derivative object by input data
RawSamples<-rep(Samples, each=2)
FilesIn <- paste(rep(File, times = length(Hours)), "_", Hours, "h.pda.txt", sep="")
Labels<-paste(rep(Label, times = length(Hours)), " ", Hours, "h", sep="")
FilesOut <- paste(rep(File, times = length(Hours)), "_Cuvettes_", Hours, "h.ps", sep="")
FilesDuplOut <- paste(rep(File, times = length(Hours)), "_", Hours, "h.ps", sep="")
TableTPOut <- paste(rep(File, times = length(Hours)), "_", Hours, "h.csv", sep="")
SummaryResults <- vector("list", length = length(Hours))
FileOut <- paste(File, ".ps", sep = "")

for (i in c(1:length(Hours)))
{

# Read the file of time point and cut off everything above actual data
InputData <- head(read.table(file=FilesIn[i], sep="\t", fill=T), -2)
Indice<-which(InputData=="Wells", arr.ind=T)
InputData<-head(read.table(file=FilesIn[i], sep="\t", fill=T, skip=Indice[1]), -2)
Indice<-which(InputData=="Wells", arr.ind=T)
InputData <- InputData[-1*c(1:Indice[1]),]

# Write meaningful columns
InputData[,1] <- rep(rep(Samples, each = Plicates), each = 2)
InputData[,2] <- rep(rep(c(" 1st  ", " 2nd  "), each = Plicates), times = length(Samples))
head(InputData)
CleanData <- InputData[,c(1,2,4)]
colnames(CleanData) <- c("Sample", "Cuvette", "OD")

PlotDataSingle <- ddply(melt(CleanData), .(Sample, Cuvette), summarise,
                        mean = mean(value), min = min(value), max = max(value))
SingleErrBars <- aes(ymin = PlotDataSingle$min, ymax = PlotDataSingle$max)
YLim <- ceiling(max(PlotDataSingle$max)*10/2)*2/10

g<-ggplot(PlotDataSingle, aes(x=Sample, y=mean, fill=Cuvette, width=0.7 ))
g1<-g+geom_bar(position = position_dodge(width = 0.75), stat="identity")+
  scale_fill_manual(values=c("gray90","gray20"))
g2<-g1+scale_x_discrete(limits=SamOrder)
g3<-g2+geom_errorbar(SingleErrBars, position=position_dodge(width=0.75), width=0.2, size=1.2)
g4<-g3+scale_y_continuous(breaks=c(seq(0, YLim, by=0.2)),limits=c(0,YLim),expand=c(0,0))
g5<-g4+xlab("")+ylab("Relative viability")+ggtitle(Labels[i])
g6<-g5+theme_classic()+theme(axis.title = TxtSize, axis.text = TxtSize,aspect.ratio=AR,
                             plot.margin=unit(c(1,1,1,3), "line"), plot.title=TxtSize,
                             axis.title.y=element_text(vjust=2), 
                             axis.line = AxLines, axis.ticks = AxLines,
                             panel.grid.major = GridType, panel.grid.minor = GridType)
g7<-g6+ guides(fill = guide_legend(title="Cuvette: ", title.position="left", 
                                   keywidth=1.5, keyheight=1, 
                                   title.hjust=0.5, title.theme=LegTxtSize, 
                                   label.theme=LegTxtSize, label.position = "right",
                                   ncol=2))+
  theme(legend.position="top", legend.key=element_rect(colour="black", size=1))
g8<-g7+geom_bar(position = position_dodge(width = 0.75),
                colour="black", stat="identity", size=1.5, show_guide=FALSE)+
  scale_fill_manual(values=c("gray90","gray20"))+
  geom_errorbar(SingleErrBars, position=position_dodge(width=0.75), width=0.2, size=1.2)

# print(g8)
postscript(file = FilesOut[i],family = "AvantGarde")
print(g8)
dev.off()
system(paste("gv ", FilesOut[i], " &", sep = ""))

TimePointTable <- left_join(unique(CleanData[,1:2]), PlotDataSingle, by = colnames(CleanData[,1:2]))
TimePointTable$Sample <- gsub(pattern = "\n", replacement = " ", x = TimePointTable$Sample)
write.csv(file = TableTPOut[i], x = TimePointTable, row.names = T)

PlotDataDupl <- ddply(melt(PlotDataSingle[,c(1,3)], id.vars = "Sample"), .(Sample), summarise,
                      mean = mean(value), min = min(value), max = max(value))

# Check whether more than one time-point is analysed.
# if only one the Plot for mean time-point value is the same as general one
if (length(Hours) > 1) {
DuplErrBars <- aes(ymin = PlotDataDupl$min, ymax = PlotDataDupl$max)
YLimDupl <- ceiling(max(PlotDataDupl$max)*10/2)*2/10


g<-ggplot(PlotDataDupl, aes(x=Sample, y=mean, width=0.7))
g1<-g+geom_bar(position = position_dodge(width = 0.75), stat="identity", fill = "gray70", 
               colour = "black", size = 2)
g2<-g1+scale_x_discrete(limits=SamOrder)
g3<-g2+geom_errorbar(DuplErrBars, position=position_dodge(width=0.75), width=0.2, size=2)
g4<-g3+scale_y_continuous(breaks=c(seq(0, YLim, by=0.2)),limits=c(0,YLim),expand=c(0,0))
g5<-g4+xlab("")+ylab("Relative viability")+ggtitle(Labels[i])
g6<-g5+theme_classic()+theme(axis.title = TxtSize, axis.text = TxtSize,aspect.ratio=AR*1.5,
                             plot.margin=unit(c(1,1,1,3), "line"), plot.title=TxtSize,
                             axis.title.y=element_text(vjust=2), 
                             axis.line = AxLines, axis.ticks = AxLines,
                             panel.grid.major = GridType, panel.grid.minor = GridType)
# print(g6)

postscript(file = FilesDuplOut[i],family = "AvantGarde")
print(g6)
dev.off()
system(paste("gv ", FilesDuplOut[i], " &", sep = ""))
}
# Place the table of means for time point into the list of results

ListData <- data.frame("Sample" = PlotDataDupl[,1], 
                       "TP" = rep(paste(Hours[i], "hour(s)", sep = ""), times = length(Samples)),
                       PlotDataDupl[,2:4])
SummaryResults[[i]] <- ListData
}

CombineData <- rbind_all(SummaryResults)
CombErrBars <- aes(ymin = CombineData$min, ymax = CombineData$max)
YLimComb <- ceiling(max(CombineData$max)*10/2)*2/10

g<-ggplot(CombineData, aes(x=Sample, y=mean, fill=TP, width=0.7 ))
g1<-g+geom_bar(position = position_dodge(width = 0.75), stat="identity")+
  scale_fill_manual(values=c("gray90","gray50", "gray20"))
g2<-g1+scale_x_discrete(limits=SamOrder)
g3<-g2+geom_errorbar(CombErrBars, position=position_dodge(width=0.75), width=0.2, size=1.2)
g4<-g3+scale_y_continuous(breaks=c(seq(0, YLimComb, by=0.2)),limits=c(0,YLimComb),expand=c(0,0))
g5<-g4+xlab("")+ylab("Relative viability")+ggtitle(GenLabel)
g6<-g5+theme_classic()+theme(axis.title = TxtSize, axis.text = TxtSize,aspect.ratio=AR*1.35,
                             plot.margin=unit(c(1,1,1,3), "line"), plot.title=TxtSize,
                             axis.title.y=element_text(vjust=2), title = element_text(vjust = 2),
                             axis.line = AxLines, axis.ticks = AxLines,
                             panel.grid.major = GridType, panel.grid.minor = GridType)
g7<-g6+ guides(fill = guide_legend(title="Time: ", title.position="left",
                                   keywidth = 1, keyheight=1.5, 
                                   title.hjust=0.5, title.theme=LegTxtSize, 
                                   label.theme=LegTxtSize2, label.position = "right",
                                   ncol=1))+
  theme(legend.position="bottom", legend.key=element_rect(colour="black", size=1))
g8<-g7+geom_bar(position = position_dodge(width = 0.75),
                colour="black", stat="identity", size=1.5, show_guide=FALSE)+
  scale_fill_manual(values=c("gray90","gray50", "gray20"))+
  geom_errorbar(CombErrBars, position=position_dodge(width=0.75), width=0.2, size=1.2)

#print(g8)
postscript(file = FileOut,family = "AvantGarde")
print(g8)
dev.off()
system(paste("gv ", FileOut, " &", sep = ""))

