

#Name is first column, 
do_boxplot_basic <- function(suffix, prefix, main, ylab, use31p=TRUE, allneg=FALSE, transform = identity){
  
  vb_F<-read.table(paste("../results/Vb_functional_10_Nov_2016",prefix,suffix,sep=""),header=F,sep="\t",stringsAsFactors=FALSE)
  
  vb_31<-read.table(paste("../results/Vb31_10_Nov_2016",prefix,suffix,sep=""),header=F,sep="\t",stringsAsFactors=FALSE)
  
  vb <- vb_F
  if(use31p){
    vb<-rbind(vb,vb_31)  
  }
  
  colnames(vb)<-c("name","x")
  
  cRSS_pos_nonz<-read.table(paste("../results/cRSS_nonz_pos_10_Nov_2016",prefix,suffix,sep=""),header=F,sep="\t",stringsAsFactors=FALSE)  
  cRSS_pos_z<-read.table(paste("../results/cRSS_z_pos_10_Nov_2016",prefix,suffix,sep=""),header=F,sep="\t",stringsAsFactors=FALSE)  
  cRSS_pos<-rbind(cRSS_pos_nonz,cRSS_pos_z)
  colnames(cRSS_pos)<-c("name","x")
  
  vb_pseudo_pos<-read.table(paste("../results/Vb_pseudo_pos_10_Nov_2016",prefix,suffix,sep=""),header=F,sep="\t",stringsAsFactors=FALSE)
  colnames(vb_pseudo_pos)<-c("name","x")
  
  if(allneg){
    vb_pseudo_neg<-read.table(paste("../results/Vb_pseudo_neg_10_Nov_2016.regions.CA",prefix,suffix,sep=""),header=F,sep="\t",stringsAsFactors=FALSE)  
  } else {
    vb_pseudo_neg<-read.table(paste("../results/Vb_pseudo_neg_10_Nov_2016.regions.CA.max_REC",prefix,suffix,sep=""),header=F,sep="\t",stringsAsFactors=FALSE)  
  }
  colnames(vb_pseudo_neg)<-c("name","x")
  
  if(allneg){
    cRSS_neg_nonz<-read.table(paste("../results/cRSS_nonz_neg_10_Nov_2016.regions.CA",prefix,suffix,sep=""),header=F,sep="\t",stringsAsFactors=FALSE)  
    cRSS_neg_z<-read.table(paste("../results/cRSS_z_neg_10_Nov_2016.regions.CA",prefix,suffix,sep=""),header=F,sep="\t",stringsAsFactors=FALSE)  
  } else {
    cRSS_neg_nonz<-read.table(paste("../results/cRSS_nonz_neg_10_Nov_2016.regions.CA.max_REC",prefix,suffix,sep=""),header=F,sep="\t",stringsAsFactors=FALSE)  
    cRSS_neg_z<-read.table(paste("../results/cRSS_z_neg_10_Nov_2016.regions.CA.max_REC",prefix,suffix,sep=""),header=F,sep="\t",stringsAsFactors=FALSE)  
  }
  cRSS_neg<-rbind(cRSS_neg_nonz,cRSS_neg_z)
  colnames(cRSS_neg)<-c("name","x")
  
  vb$x=transform(vb$x)
  cRSS_pos$x = transform(cRSS_pos$x)
  vb_pseudo_pos$x = transform(vb_pseudo_pos$x)
  cRSS_neg$x = transform(cRSS_neg$x)
  vb_pseudo_neg$x = transform(vb_pseudo_neg$x)
  
  
  pos<-rbind(cRSS_pos,vb_pseudo_pos)
  neg<-rbind(cRSS_neg,vb_pseudo_neg)

  vb$type<-"vRSS"
  pos$type<-"cRSS"
  neg$type<-"nonRSS"
  
  all<-rbind(vb,pos,neg)
  
  return( all )
  
}


do_boxplot <- function(suffix,main,ylab,use31p=TRUE,allneg=FALSE,prefix=".500bp_around.50bp_windows.txt.", ylims=c(-500,500), transform = identity, aggfun = mean){
  
  vb_F<-read.table(paste("../results/Vb_functional_10_Nov_2016",prefix,suffix,sep=""),header=F,sep="\t",stringsAsFactors=FALSE)
  vb_31<-read.table(paste("../results/Vb31_10_Nov_2016",prefix,suffix,sep=""),header=F,sep="\t",stringsAsFactors=FALSE)
  
  vb_F$window<-vb_F$V2 - vb_F$V4
  
  #Need to reverse vb_31 (carefull may need to remove border windows)
  vb_31$window<-vb_31$V4 - vb_31$V2  
  
  vb <- vb_F
  if(use31p){
    vb<-rbind(vb,vb_31)  
  }
  
  #all the rest are forward...
  cRSS_pos_nonz<-read.table(paste("../results/cRSS_nonz_pos_10_Nov_2016",prefix,suffix,sep=""),header=F,sep="\t",stringsAsFactors=FALSE)  
  cRSS_pos_z<-read.table(paste("../results/cRSS_z_pos_10_Nov_2016",prefix,suffix,sep=""),header=F,sep="\t",stringsAsFactors=FALSE)  
  cRSS_pos<-rbind(cRSS_pos_nonz,cRSS_pos_z)
  
  vb_pseudo_pos<-read.table(paste("../results/Vb_pseudo_pos_10_Nov_2016",prefix,suffix,sep=""),header=F,sep="\t",stringsAsFactors=FALSE)
  
  if(allneg){
    vb_pseudo_neg<-read.table(paste("../results/Vb_pseudo_neg_10_Nov_2016.regions.CA",prefix,suffix,sep=""),header=F,sep="\t",stringsAsFactors=FALSE)
  } else {
    vb_pseudo_neg<-read.table(paste("../results/Vb_pseudo_neg_10_Nov_2016.regions.CA.max_REC",prefix,suffix,sep=""),header=F,sep="\t",stringsAsFactors=FALSE)
  }
  
  if(allneg){
    cRSS_neg_nonz<-read.table(paste("../results/cRSS_nonz_neg_10_Nov_2016.regions.CA",prefix,suffix,sep=""),header=F,sep="\t",stringsAsFactors=FALSE)
    cRSS_neg_z<-read.table(paste("../results/cRSS_z_neg_10_Nov_2016.regions.CA",prefix,suffix,sep=""),header=F,sep="\t",stringsAsFactors=FALSE)
  } else {
    cRSS_neg_nonz<-read.table(paste("../results/cRSS_nonz_neg_10_Nov_2016.regions.CA.max_REC",prefix,suffix,sep=""),header=F,sep="\t",stringsAsFactors=FALSE)
    cRSS_neg_z<-read.table(paste("../results/cRSS_z_neg_10_Nov_2016.regions.CA.max_REC",prefix,suffix,sep=""),header=F,sep="\t",stringsAsFactors=FALSE)
  }
  
  cRSS_neg<-rbind(cRSS_neg_z, cRSS_neg_nonz)
  
  cRSS_pos$window<-cRSS_pos$V2 - cRSS_pos$V4  
  cRSS_neg$window<-cRSS_neg$V2 - cRSS_neg$V4
  vb_pseudo_pos$window<-vb_pseudo_pos$V2 - vb_pseudo_pos$V4
  vb_pseudo_neg$window<-vb_pseudo_neg$V2 - vb_pseudo_neg$V4
  
  #aggregate the filtered data on RSS $V4 (RSS name)
  data_mean_vb<-aggregate(vb$V8[vb$window>=ylims[1] & vb$window<=ylims[2]],list(name=vb$V4[vb$window>=ylims[1] & vb$window<=ylims[2]]),aggfun)  
  data_mean_cRSS_pos<-aggregate(cRSS_pos$V8[cRSS_pos$window>=ylims[1] & cRSS_pos$window<=ylims[2]],list(name=cRSS_pos$V4[cRSS_pos$window>=ylims[1] & cRSS_pos$window<=ylims[2]]),aggfun)  
  data_mean_cRSS_neg<-aggregate(cRSS_neg$V8[cRSS_neg$window>=ylims[1] & cRSS_neg$window<=ylims[2]],list(name=cRSS_neg$V4[cRSS_neg$window>=ylims[1] & cRSS_neg$window<=ylims[2]]),aggfun)  
  data_mean_vb_pseudo_pos<-aggregate(vb_pseudo_pos$V8[vb_pseudo_pos$window>=ylims[1] & vb_pseudo_pos$window<=ylims[2]],list(name=vb_pseudo_pos$V4[vb_pseudo_pos$window>=ylims[1] & vb_pseudo_pos$window<=ylims[2]]),aggfun)    
  data_mean_vb_pseudo_neg<-aggregate(vb_pseudo_neg$V8[vb_pseudo_neg$window>=ylims[1] & vb_pseudo_neg$window<=ylims[2]],list(name=vb_pseudo_neg$V4[vb_pseudo_neg$window>=ylims[1] & vb_pseudo_neg$window<=ylims[2]]),aggfun)   
  
  data_mean_vb$x<-transform(data_mean_vb$x)
  data_mean_cRSS_pos$x<-transform(data_mean_cRSS_pos$x)
  data_mean_vb_pseudo_pos$x<-transform(data_mean_vb_pseudo_pos$x)  
  data_mean_cRSS_neg$x<-transform(data_mean_cRSS_neg$x)
  data_mean_vb_pseudo_neg$x<-transform(data_mean_vb_pseudo_neg$x)
  
  data_mean_pos<-rbind(data_mean_cRSS_pos,data_mean_vb_pseudo_pos)
  data_mean_neg<-rbind(data_mean_cRSS_neg,data_mean_vb_pseudo_neg)
  
  data_mean_vb$type<-"vRSS"
  data_mean_pos$type<-"cRSS"
  data_mean_neg$type<-"nonRSS"
  
  data_mean_all<-rbind(data_mean_vb,data_mean_pos,data_mean_neg)
  
  return( data_mean_all )
  
}

do_merge<-function(input,regionids,var){
  merged<-merge(input,regionids, by.x="name", by.y="name", all.x=TRUE)
  merged$name[!is.na(merged$region)]<-merged$region[!is.na(merged$region)]
  merged<-aggregate(merged$x,list(merged$name,merged$type),mean)
  colnames(merged)<-c("name","type",var)
  return(merged)
}


regionids_pseudo<-read.table("../results/Vb_pseudo_neg_10_Nov_2016.regions.CA.ref.txt",sep="\t",header=F)
regionids_crss_nonz<-read.table("../results/cRSS_nonz_neg_10_Nov_2016.regions.CA.ref.txt",sep="\t",header=F)
regionids_crss_z<-read.table("../results/cRSS_z_neg_10_Nov_2016.regions.CA.ref.txt",sep="\t",header=F)
regionids<-rbind(regionids_pseudo,regionids_crss_nonz,regionids_crss_z)
colnames(regionids)<-c("name","region")

eb_maxrec<-do_boxplot(suffix="Eb-Tcell.batch2.s1.bed",main="Eb (s1) mean of Counts in 50bp windows along 500bp around RSS",ylab="3C Counts",use31p=TRUE,allneg=FALSE)
colnames(eb_maxrec)[2]<-"Eb"

faire_maxrec<-do_boxplot(suffix="FAIRE_Rag_unpublished.gff.processed.tab",main="FAIRE mean Chip score in 50bp windows 250bp downstream RSS",ylab="Chip score",use31p=TRUE,allneg = FALSE,ylims=c(0,250))
colnames(faire_maxrec)[2]<-"FAIRE"

h3k27Ac_maxrec<-do_boxplot(suffix="GSM1051313-27034_H3K27Ac.txt.processed.NCBIM37.bed.chr6.bed",main="H3K27Ac mean Chip score in 50bp windows along 500bp around RSS",ylab="Chip score",use31p=TRUE,allneg = FALSE)
colnames(h3k27Ac_maxrec)[2]<-"H3K27ac"

p300_maxrec<-do_boxplot(suffix="GSM1051314-27036_P300.txt.processed.NCBIM37.bed.chr6.bed",main="P300 mean Chip score in 50bp windows along 500bp around RSS","Chip Score",use31p=TRUE,allneg = FALSE)
colnames(p300_maxrec)[2]<-"P300"

h3ac_maxrec<-do_boxplot(suffix="GSM1185222_Ach3_DN.wig.NCBIM37.bed.chr6.bed",main="Ach3 mean Log Counts in 50bp windows 500bp upstream RSS",ylab="Log Counts",use31p=TRUE,ylims=c(-500,0),transform=log,allneg = FALSE)
h3ac_maxrec$x[is.infinite(h3ac_maxrec$x)]<-NA
colnames(h3ac_maxrec)[2]<-"H3Ac"

h3k4me2_maxrec<-do_boxplot(suffix="H3k4me2-KO-RAG_published.gff.processed.tab",main="H3k4me2 mean Chip Score in 50bp windows along 500bp around RSS",ylab="Chip Score",use31p=TRUE,allneg = FALSE)
colnames(h3k4me2_maxrec)[2]<-"H3K4me2"

h3k9me2_maxrec<-do_boxplot(suffix="H3k9-KO-RAG_published.gff.processed.tab",main="H3k9 mean Chip Score in 50bp windows along 500bp around RSS",ylab="Chip Score",use31p=TRUE,allneg = FALSE)
colnames(h3k9me2_maxrec)[2]<-"H3K9me2"

h3k27me3_maxrec<-do_boxplot(suffix="H3k27-KO-RAG_published.gff.processed.tab",main="H3k27 mean Chip Score in 50bp windows along 500bp around RSS",ylab="Chip Score",use31p=TRUE,allneg = FALSE)
colnames(h3k27me3_maxrec)[2]<-"H3K27me3"

h3k36me3_maxrec<-do_boxplot(suffix="H3k36-KO-RAG_published.gff.processed.tab",main="H3k36 mean Chip Score in 50bp windows along 500bp around RSS",ylab="Chip Score",use31p=TRUE,allneg = FALSE)
colnames(h3k36me3_maxrec)[2]<-"H3K36me3"

pol2_seq_r1_maxrec<-do_boxplot(suffix="GSM1360720_Pol2_Rag_N20_rep1.NCBIM37.chr6.bed",main="Pol_II Rep1 mean counts in 50bp windows along 250bp around RSS",ylab="Log Counts",use31p=TRUE,ylims=c(-250,250),allneg = FALSE)
colnames(pol2_seq_r1_maxrec)[2]<-"PolII_Rep1"
pol2_seq_r2_maxrec<-do_boxplot(suffix="GSM1360721_Pol2_Rag_N20_rep2.NCBIM37.chr6.bed",main="Pol_II Rep2 mean counts in 50bp windows along 250bp around RSS",ylab="Log Counts",use31p=TRUE,ylims=c(-250,250),allneg = FALSE)
colnames(pol2_seq_r2_maxrec)[2]<-"PolII_Rep2"
pol2_maxrec<-merge(pol2_seq_r1_maxrec,pol2_seq_r2_maxrec,by.x=c('name','type'),by.y=c('name','type'),all=T)
pol2_maxrec[,5]<-log((pol2_maxrec$PolII_Rep1+pol2_maxrec$PolII_Rep2)/2)
colnames(pol2_maxrec)[5]<-"PolII"
pol2_maxrec$PolII_Rep1<-NULL
pol2_maxrec$PolII_Rep2<-NULL
pol2_maxrec[is.infinite(pol2_maxrec[,3]),3]<-NA

h3k4me1_maxrec<-do_boxplot(suffix="Rag_CIML_H3K4me1_146_146_merged__no_artefacts.wig.NCBIM37.bed.chr6.bed",main="H3K4me1 Chip Seq mean Log Counts in 50bp windows along 500bp around RSS",ylab="Log Counts",use31p=TRUE,transform=log,allneg = FALSE)
h3k4me1_maxrec$x[is.infinite(h3k4me1_maxrec$x)]<-NA
colnames(h3k4me1_maxrec)[2]<-"H3K4me1"

h3k4me3_maxrec<-do_boxplot(suffix="Rag_CIML_H3K4me3_136_136_merged__no_artefacts.wig.NCBIM37.bed.chr6.bed",main="H3K4me3 Chip Seq mean Log Counts in 50bp windows along 500bp around RSS",ylab="Log Counts",use31p=TRUE,transform=log,allneg = FALSE)
h3k4me3_maxrec$x[is.infinite(h3k4me3_maxrec$x)]<-NA
colnames(h3k4me3_maxrec)[2]<-"H3K4me3"

gltm_maxrec<-do_boxplot(suffix="WIGfs_DN_Total_negativeStrand_el0_bins50_(ArtefactParam_NA)_Illumina.wig.NCBIM37.bed.chr6.bed",main="Total RNA (- strand) mean log counts in 50bp windows along 500bp around RSS",ylab="Log Counts",use31p=TRUE,allneg = FALSE)
colnames(gltm_maxrec)[2]<-"GLT-"

gltp_maxrec<-do_boxplot(suffix="WIGfs_DN_Total_positiveStrand_el0_bins50_(ArtefactParam_NA)_Illumina.wig.NCBIM37.bed.chr6.bed",main="Total RNA (+ strand) mean log counts in 50bp windows along 250bp around RSS",ylab="Log Counts",use31p=TRUE, ylims=c(-250,250),allneg = FALSE)
colnames(gltp_maxrec)[2]<-"GLT+"

ctcf_alt_maxrec<-do_boxplot_basic(suffix="GSM1023416_CTCF_DN_peaks_track.NCBIM37.chr6.bed",main="Closest CTCF from RSS",ylab=" Log Distance (bp)",use31p=TRUE,prefix=".point.txt.", transform=function(x){log(abs(x+1))},allneg = FALSE)
colnames(ctcf_alt_maxrec)[2]<-"CTCFdist"

rec_maxrec<-do_boxplot_basic(suffix="",main="REC Score",ylab="rec",use31p=TRUE,prefix='.REC.txt',allneg = FALSE)
colnames(rec_maxrec)[2]<-"REC"

nonamer_maxrec<-do_boxplot_basic(suffix="",main="Nonamer Scores (all CAs in negative regions)",ylab="nonamer score",use31p=TRUE,prefix='.nonamer.scores',allneg = FALSE)
colnames(nonamer_maxrec)[2]<-"nonamer"

heptamer_maxrec<-do_boxplot_basic(suffix="",main="heptamer Scores (all CAs in negative regions)",ylab="heptamer score",use31p=TRUE,prefix='.heptamer.scores',allneg = FALSE)
colnames(heptamer_maxrec)[2]<-"heptamer"

map75_maxrec<-do_boxplot(suffix="wgEncodeCrgMapabilityAlign75mer.bedGraph.NCBIM37.chr6.bed",main="Mean 75bp Mappability in 50bp windows along 500bp around RSS",ylab="Mappability",use31p=TRUE,allneg = FALSE)
colnames(map75_maxrec)[2]<-"map75"

# Problem: on meanCA, identifier is beginning of region, while for other variables, it is the locus of the best REC...
list.of.data.frames.maxREC <- list(eb_maxrec, faire_maxrec, h3k27Ac_maxrec, p300_maxrec, h3ac_maxrec, h3k4me2_maxrec, h3k9me2_maxrec, h3k27me3_maxrec, h3k36me3_maxrec, pol2_maxrec, h3k4me1_maxrec, h3k4me3_maxrec, gltm_maxrec, gltp_maxrec, ctcf_alt_maxrec, rec_maxrec, heptamer_maxrec, nonamer_maxrec, map75_maxrec)
merged.data.maxREC = Reduce(function(...) merge(..., by=c("name","type"),all=T), list.of.data.frames.maxREC)
merged.data.maxREC$type<-as.factor(merged.data.maxREC$type)

#Now I don't need that info on
merged.data.maxREC$map75 <- NULL
rownames(merged.data.maxREC)<-merged.data.maxREC$name
merged.data.maxREC$name<-NULL
colnames(merged.data.maxREC)[2]<-"3C-Eb"


library(reshape2)
library(ggplot2)
library(Rmisc)

do_summary_plot<-function(input,names){
  merged.data<-input
  merged.data[,-c(1)]<-apply(merged.data[,-c(1)],2,scale)
  means<-aggregate(merged.data[,-c(1,2)],by=list(merged.data$type),mean,na.rm=TRUE)
  means.long<-melt(means,id.vars="Group.1")
  melted.data<-melt(merged.data,id="type")
  cbPalette <- c("#555555", "#888888", "#BBBBBB")
  melted.data$type<-ordered(melted.data$type, levels = c("vRSS", "cRSS", "nonRSS"))
  summaries<-summarySE(melted.data,measurevar = "value", groupvars = c("type","variable"),na.rm = TRUE)
  print(summaries)
  
  ggplot(summaries, aes(x=variable, y=value, fill=type)) + geom_bar(position=position_dodge(width=0.7), stat="identity",width=0.8) + geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0.5, position=position_dodge(width=0.7))  + scale_fill_manual("RSS Type",values=cbPalette) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=14), axis.text.y = element_text(size=18),axis.title.x=element_text(size=18), axis.title.y=element_text(size=18),legend.text = element_text(size = 16),legend.title=element_text(size = 18)) + ylab("Normalized Value") + scale_x_discrete(labels=names)
}

list.of.data.frames.maxREC <- list(eb_maxrec, faire_maxrec, h3k27Ac_maxrec, p300_maxrec, h3ac_maxrec, h3k4me2_maxrec, h3k9me2_maxrec, h3k27me3_maxrec, h3k36me3_maxrec, pol2_maxrec, h3k4me1_maxrec, h3k4me3_maxrec, gltm_maxrec, gltp_maxrec, ctcf_alt_maxrec, rec_maxrec, heptamer_maxrec, nonamer_maxrec, map75_maxrec)
merged.data.maxREC = Reduce(function(...) merge(..., by=c("name","type"),all=T), list.of.data.frames.maxREC)
merged.data.maxREC$type<-as.factor(merged.data.maxREC$type)
merged.data.maxREC$name<-NULL

names(merged.data.maxREC)[14]<-"GLTm"
names(merged.data.maxREC)[15]<-"GLTp"

#imputation
library(mice)
merged.data.maxREC<-complete(mice(merged.data.maxREC))


colnames(merged.data.maxREC)[2]<-"3C-Eb"
colnames(merged.data.maxREC)[14]<-"GLT-"
colnames(merged.data.maxREC)[15]<-"GLT+"

#Fig5A
merged.data.maxREC.selection<-merged.data.maxREC[,c("type","REC", "CTCFdist", "3C-Eb","H3K9me2", "FAIRE","H3Ac","H3K4me3","PolII","GLT+")]
names<-colnames(merged.data.maxREC.selection)
names[4]<-expression(paste("3C-E",beta))
names[1]<-NULL
png(filename="Fig5a.png",width=1200, height=600,units="px",type="cairo")
do_summary_plot(merged.data.maxREC.selection, names)
dev.off()
#End Fig5A


#Fig S7
merged.data.maxREC.numeric<-merged.data.maxREC
merged.data.maxREC.numeric$type<-NULL
merged.data.maxREC.numeric$map75<-NULL
library(gplots)
colpal<-colorRampPalette(c("blue","white","red"))
cnames<-colnames(merged.data.maxREC.numeric)
cnames[1]<-expression(paste("3C-E",beta))
png(filename="FigS7.png",width=600, height=600,units="px",type="cairo")
heatmap.2(cor(as.matrix(scale(merged.data.maxREC.numeric)),method="spearman",use="pairwise.complete.obs"),scale='none',symm=T,margins=c(8,8),key.title = "",key.xlab = "Spearman Correlation",key.ylab = "",trace = 'none',col=colpal(50),labRow = cnames, labCol=cnames)
dev.off()
#End Fig S7


read_tcr_data<-function(input,logp=FALSE, add=0){
  tcrdata<-read.table(input,sep="\t")
  if(add>0){
    tcrdata$V3<-tcrdata$V3+add
  }
  tcrdata<-tcrdata[tcrdata$V2>40840000 & tcrdata$V3<41510000,]
  if(logp){
    tcrdata$V4<-log(tcrdata$V4)
  }
  tcrdata<-tcrdata[is.finite(tcrdata$V4),]
  
  ##Eliminate outliers from the scale
  box<-boxplot(tcrdata$V4,plot=FALSE)
  adjmin<-box$stats[1,1]
  adjmax<-box$stats[5,1]
  tcrdata$V4[tcrdata$V4 > adjmax]<-adjmax
  tcrdata$V4[tcrdata$V4 < adjmin]<-adjmin
  
  #The full scale
  tcrdata$perc<-round(100*(tcrdata$V4-min(tcrdata$V4))/(max(tcrdata$V4)-min(tcrdata$V4)))
  
  return(tcrdata)
}

threecdata<-read_tcr_data("../data/Eb-Tcell.batch2.s1.bed.tcrb.tab",logp=TRUE)
fairedata<-read_tcr_data("../data/FAIRE_Rag_unpublished.gff.processed.tab.tcrb.tab",logp=FALSE)
h3k9data<-read_tcr_data("../data/H3k9-KO-RAG_published.gff.processed.tab.tcrb.tab",logp=FALSE)
h3acdata<-read_tcr_data("../data/GSM1185222_Ach3_DN.wig.NCBIM37.bed.chr6.bed.tcrb.tab",logp=TRUE)
h3k4data<-read_tcr_data("../data/Rag_CIML_H3K4me3_136_136_merged__no_artefacts.wig.NCBIM37.bed.chr6.bed.tcrb.tab",logp=TRUE,add=50)
poliidata<-read_tcr_data("../data/PolII_N20_Rag_bins_100_elongation_245_published.wig.NCBIM37.bed.chr6.bed.tcrb.tab", logp=TRUE, add=100)
gltpdata<-read_tcr_data("../data/WIGfs_DN_Total_positiveStrand_el0_bins50_(ArtefactParam_NA)_Illumina.wig.NCBIM37.bed.chr6.bed.tcrb.tab", logp=TRUE,add=50)

vRSSf<-read.table("../data/Vb_functional_10_Nov_2016.txt")
v31<-read.table("../data/Vb31_10_Nov_2016.txt")
vRSS<-rbind(vRSSf,v31)
pRSS_pos<-read.table("../data/Vb_pseudo_pos_10_Nov_2016.txt")
pRSS_neg<-read.table("../data/Vb_pseudo_neg_10_Nov_2016.txt")
pRSS<-rbind(pRSS_pos,pRSS_neg)
tested<-read.table("../data/6_Tables_SUPP_150822_S3_coord_tested_sites.tab",header=T,sep="\t")
primers<-read.table("../data/6_Tables_SUPP_150822_S1_primers.csv.NCBIM38.bowtie.sam.coords.tcrb",header=FALSE,row.names=1)
primers<-primers[-grep("second",rownames(primers)),]
cRSS_nonz_pos<-read.table("../data/cRSS_nonz_pos_10_Nov_2016.txt")
cRSS_nonz_neg<-read.table("../data/cRSS_nonz_neg_10_Nov_2016.txt")
cRSS_z_pos<-read.table("../data/cRSS_z_pos_10_Nov_2016.txt")
cRSS_z_neg<-read.table("../data/cRSS_z_neg_10_Nov_2016.txt")
cRSS_pos<-rbind(cRSS_nonz_pos,cRSS_z_pos,pRSS_pos)
cRSS_neg<-rbind(cRSS_nonz_neg,cRSS_z_neg,pRSS_neg)

#One of the Z gave reverse strand. Correct the strand
primers['Z10',1]<-0

nonzprimers<-primers[-grep("Z",rownames(primers)),c(1,2)]
zdata<-read.table("../data/tcrbeta_ncbim37.bed.fasta.abcc_zdna.tab",header=T,sep="\t")
zdata<-zdata[zdata$Start>10000 & zdata$Stop<680000,]
recdata<-read.table("../data/tcrbeta_ncbim37.bed.fasta.CA.fasta.REC.tab",sep="\t")
recdata<-recdata[recdata$V2>10000 & recdata$V2<680000,]
recdata_stringent<-recdata[recdata$V8==1,]
recdata_sensitive<-recdata[recdata$V9==1,]
ctcfdata<-read.table("../data/GSM1023416_CTCF_DN_peaks_track.NCBIM37.chr6.bed.tcrb.tab",sep="\t")

plot_data<-function(xstart,xend, extra=0, space=30000, title="TCRbeta locus", legendp=TRUE, legendcex=0.8, legendadj=1, plotaxp=TRUE){
  
  plot(c(40800000,41500000),c(0,11),type='n',axes=FALSE,ylab="",xlab="",main=title,ylim=c(3.5,11),xlim=c(xstart-space,xend))
  if(plotaxp){
    axis(side=3,at=c(xstart,xend))
  }
  
  rect(vRSS$V1,9.5,vRSS$V1,10,border ="red",lwd=2)
  rect(pRSS$V1,9.5,pRSS$V1,10,border="orange",lwd=2)

  rect(primers[grep("Z",rownames(primers)),2],9,primers[grep("Z",rownames(primers)),2]+2000,9.5,density=90, angle=30,border = "transparent",col="grey")
  rect(nonzprimers[nonzprimers$V2==0,2],9,nonzprimers[nonzprimers$V2==0,2]+2000,9.5,density=90, angle=30,border = "transparent",col="grey")
  rect(nonzprimers[nonzprimers$V2==16,2]-2000,9,nonzprimers[nonzprimers$V2==16,2],9.5,density=90, angle=30,border = "transparent",col="darkgrey")

  rect(tested$Position,9,tested$Position,9.5,border ="black",lwd=1)
  
  rect(pRSS_pos$V1,8.5,pRSS_pos$V1,9,border ="red",lwd=1)
  rect(cRSS_nonz_pos$V1,8.5,cRSS_nonz_pos$V1,9,border ="red",lwd=1)
  rect(cRSS_z_pos$V1,8.5,cRSS_z_pos$V1,9,border ="red",lwd=1)

  rect(recdata_stringent$V1+recdata_stringent$V2,8,recdata_stringent$V1+recdata_stringent$V2+1,8.5, col ="black",border="black")

  rect(zdata$Sequence_name+zdata$Start,7.5,zdata$Sequence_name+zdata$Stop,8, col ="darkgrey",border="darkgrey")

  rect(ctcfdata$V2,7,ctcfdata$V3,7.5, col="darkgrey", border="darkgrey")

  colpalette<-colorRampPalette(c("blue", "tomato"))
  cols<-colpalette(101)

  datalist<-list(threecdata,h3k9data,fairedata,h3acdata,h3k4data,poliidata,gltpdata)
  legendlist<-c("Authentic","Tested","Detected","REC+","ZDNA+","CTCF",expression(paste("3C-E", beta)),"H3K9me2","FAIRE","H3AC","H3K4me3","PolII","GLT+")

  i<-6.5
  j<-7
  for(data in datalist){
    
    rect(data$V2,i,data$V3,j,col=cols[data$perc+1],border=cols[data$perc+1])
    i<-i-0.5
    j<-j-0.5
  }

  if(legendp){
    rect(xstart-1000000,0,xstart+space,10,col="white",border="white")
    p<-9.75
    for(t in legendlist){ 
      text(xstart+extra,p,t,cex=legendcex,adj = legendadj)
      #text(41530000,p,t,cex=0.8,adj = 0)
      p<-p-0.5 
    }
  }

}

#Fig4
png(filename="Fig4.png",width=1200, height=600,units="px",type="cairo")
par(cex.axis=1.3, cex.main=1.5)
plot_data(40810000,41520000, extra=15000, space=30000,title=expression(paste("TCR",beta," locus")),legendcex=1.3, plotaxp = TRUE)
text(40841885,10.5,expression(paste("V",beta,"1")),cex=1)
text(40997995,10.5,expression(paste("V",beta,"2")),cex=1)
text(41231990,10.5,expression(paste("V",beta,"30")),cex=1)
text(41507653,10.5,expression(paste("V",beta,"31")),cex=1)
dev.off()

#FigS6
png(filename="FigS6.png",width=1200, height=600,units="px",type="cairo")
par(cex.axis=1.3, cex.main=1.5)
par(mfrow=c(1,2))
par(mar=c(3,3,3,3))
plot_data(40842000,40846000, title=expression("Z2 / A"), extra=200, space=800, legendp=TRUE,legendcex=1.3)
plot_data(41036800,41040000, title=expression(paste("V",beta,"9")), extra=200, space=800, legendp=TRUE, legendcex=1.3)
dev.off()
#END FigS6

#Scale
#png("FigS6_scale.png",width=6, height=2)
#colpalette<-colorRampPalette(c("blue", "tomato"))
#cols<-colpalette(101)
#image(0:100,y=0,as.matrix(1:100), col=cols, ylab="", xlab="Relative Signal (%)",tck = 1, col.axis="white")
#axis(side=1)
#dev.off()


#Fig6
bcl11start = 109136613
bcl11stop = 109251812

bcl11exons<-read.table("../data/bcl11b_NCBI37.exons.bed",sep="\t")
bcl11introns<-read.table("../data/bcl11b_NCBI37.introns.bed",sep="\t")

primersbcl11<-read.table("../data/6_Tables_SUPP_150822_S1_primers.csv.NCBIM38.bowtie.sam.coords.bcl11",header=FALSE,row.names=1)
primersbcl11<-primersbcl11[-grep("second",rownames(primersbcl11)),]

testedbcl11<-read.table("../data/6_Tables_SUPP_150822_S5_coord_tested_sites.tab",sep="\t")
posbcl11<-read.table("../data/4_Tables_Table2_positive_sites.pos",sep="\t")
zdatabcl11<-read.table("../data/bcl11b_NCBI38.abcc.tsv",header=T,sep="\t")

recdatabcl11<-read.table("../data/bcl11b_NCBI38.revcomp.line.CA.fa.REC.tab",sep="\t",header=T)
recdatabcl11_stringent<-recdatabcl11[recdatabcl11$PASS.FAIL.REC.stringent==1,]
recdatabcl11_sensitive<-recdatabcl11[recdatabcl11$PASS.FAIL.REC.sensitive==1,]

ctcfdatabcl11<-read.table("../data/GSM1023416_CTCF_DN_peaks_track.NCBIM37.bcl11b_NCBI37.bed",sep="\t")

png(filename="Fig6.png",width=1200, height=600,units="px",type="cairo")

par(cex.axis=1.2)
plot(c(bcl11start,bcl11stop),c(0,10),type='n',axes=FALSE,ylab="",xlab="",main="Bcl11b locus",ylim=c(3.5,10))
axis(side=3,at=c(bcl11start,bcl11stop))

#Bcl11 gene annotation
rect(bcl11exons$V2,9.5,bcl11exons$V3,10,border ="red",lwd=2)
rect(bcl11introns$V2,9.74,bcl11introns$V3,9.76,border ="black",lwd=2,lty=2)

#Tested regions
rect(primersbcl11[,2]-2000,9,primersbcl11[,2],9.5,density=90, angle=30,border = "transparent",col="grey")
rect(testedbcl11[,2],9,testedbcl11[,2]+1,9.5,density=90, angle=30,border = "black",col="black")

#Detected
rect(posbcl11$V1,8.5,posbcl11$V1+1,9,border="red",lwd=2)

#REC+
rect(bcl11stop-recdatabcl11_stringent$Sequence.ID,8,bcl11stop-recdatabcl11_stringent$Sequence.ID+1,8.5, col ="black",border="black")

#ZDNA
rect(bcl11start+zdatabcl11$Start,7.5,bcl11start+zdatabcl11$Stop,8, col ="darkgrey",border="darkgrey")

#CTCF
rect(ctcfdatabcl11$V2,7,ctcfdatabcl11$V3,7.5, col="darkgrey", border="darkgrey")

#adjust to be the same scale as the tcrdata
read_bcl_data<-function(input, tcrdata, logp=FALSE, add=0){
  inputdata<-read.table(input,sep="\t")
  
  if(add>0){
    inputdata$V3<-inputdata$V3+add
  }
  if(logp){
    inputdata$V4<-log(inputdata$V4)
  }
  inputdata<-inputdata[is.finite(inputdata$V4),]
  
  #assume tcrdata is already adjusted
  inputdata$perc<-round(100*(inputdata$V4-min(tcrdata$V4))/(max(tcrdata$V4)-min(tcrdata$V4)))
  inputdata$perc[inputdata$perc>100]<-100
  inputdata$perc[inputdata$perc<0]<-0
  
  return(inputdata)
}

#Epigenetic data
h3acdatabcl11<-read_bcl_data("../data/GSM1185222_Ach3_DN.wig.NCBIM37.bcl11b_NCBI37.bed", h3acdata, logp=TRUE)
h3k4databcl11<-read_bcl_data("../data/Rag_CIML_H3K4me3_136_136_merged__no_artefacts.wig.NCBIM37.bcl11b_NCBI37.bed",h3k4data,logp=TRUE,add=50)
poliidatabcl11<-read_bcl_data("../data/PolII_N20_Rag_bins_100_elongation_245_published.wig.NCBIM37.bcl11b_NCBI37.bed", poliidata, logp=TRUE, add=100)

colpalette<-colorRampPalette(c("blue", "tomato"))
cols<-colpalette(101)

datalist<-list(h3acdatabcl11,h3k4databcl11,poliidatabcl11)
legendlist<-c("Bcl11","Tested","Detected","REC+","ZDNA+","CTCF","H3AC","H3K4me3","PolII")

i<-6.5
j<-7
for(data in datalist){
  
  #carefull: index 0 does not exist and makes vectors uneven!!
  rect(data$V2,i,data$V3,j,col=cols[data$perc+1],border=cols[data$perc+1])
  i<-i-0.5
  j<-j-0.5
}

rect(bcl11start-5000,0,bcl11start+9000,9,col="white",border="white")
p<-9.75
for(t in legendlist){ 
  text(bcl11start+7000,p,t,cex=1.3,adj = 1)
  p<-p-0.5 
}

dev.off()
#END Fig6



#Fig5b 
library(randomForest)
#Need to change column names or this thing will break...
colnames(merged.data.maxREC)[c(2,14,15)]<-c("Eb3C","GLTm","GLTp")
set.seed(111)
repeats = 100
empty<-rep(x = NA, times = repeats)
data_fvd = data.frame(cRSS_cRSS=empty,cRSS_nonRSS=empty,cRSS_vRSS=empty,nonRSS_cRSS=empty,nonRSS_nonRSS=empty,nonRSS_vRSS=empty,vRSS_cRSS=empty,vRSS_nonRSS=empty,vRSS_vRSS=empty)
for (i in 1:repeats) { 
  rf<-randomForest(type ~ . , data=merged.data.maxREC, proximity=TRUE,importance=TRUE)
  data_fvd[i,]<-as.vector(t(rf$confusion[,c(1,2,3)]/rowSums(rf$confusion[,c(1,2,3)])))
}
melted.data<-melt(data_fvd)
melted.data$variable<-ordered(melted.data$variable, levels = c("vRSS_vRSS","vRSS_cRSS","vRSS_nonRSS","cRSS_vRSS","cRSS_cRSS","cRSS_nonRSS","nonRSS_vRSS","nonRSS_cRSS", "nonRSS_nonRSS"))
#Change order to be consistent with Fig5A
grayPallete<-c("#555555","#555555","#555555","#888888","#888888","#888888","#BBBBBB","#BBBBBB","#BBBBBB")
summaries<-summarySE(melted.data,measurevar = "value", groupvars = c("variable"),na.rm = TRUE)

png(filename="Fig5b.png",width=1000, height=600,units="px",type="cairo")
ggplot(summaries, aes(x=variable, y=value, fill=variable)) + geom_bar(position=position_dodge(), stat="identity") + geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=.2, position=position_dodge(.9),colour="#000000")  + scale_fill_manual("Prediction_Real",values=grayPallete, guide=FALSE) + theme_bw() + theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"), axis.text.x = element_text(angle = 45, hjust = 1,size=14),axis.text.y=element_text(size=14) ) +xlab("")+ylab("Frequency of Prediction")
dev.off()
#END Fig5b


#Fig S1
repert_reduced<-read.table("../data/repertoires_correl_reduced.csv",header=T,row.names=1,sep="\t")
colnames(repert_reduced)<-c("Ndifon Spleen HTS","Ndifon Spleen FACS","Bonnet Lymph Nodes FACS","Matsutani DN3 MHA","Matsutani DN4 MHA","Wilson DN3 FACS","Wilson DN3L FACS", "Wilson DN4 FACS", "Li DP HTS","Gopalakrishnan DN3 qPCR")
library(gplots)
png(filename="FigS1.png",width=600, height=600,units="px",type="cairo")
heatmap.2(cor(repert_reduced,use = 'pairwise.complete.obs',method="spearman"),scale='none',symm=T,margins=c(10,10),key.title = "",key.xlab = "Spearman Correlation",key.ylab = "",trace = 'none', cexRow=1,cexCol=1,main=expression(paste("Correlation of V",beta," usage")))
dev.off()
#END Fig S1


