#This section contains all the parameters that the user should change
#Group order for plotting
#Each group MUST be named EXACTLY as in the Group folder
#The same MUST be written between ""
#Example: in a case where there are 3 groups (so 3 folders)
#Group_order<- c("Control", "Treatment1", "Treatment2")
#ALERT! IMPORTANT! --- > The number of group names MUTS
#correspond to the number of group folders

Group_order<- c("Control", "Treatment1", "Treatment2")


#FACTOR_STALL times the average step defines the variance threshold that needs to exist between to consecutive steps to be considered a moving
#Values close to 0 will constrain the variance threshold, therefore allow for more steps to be considered movement
#Values close to 1 will increase the variance threshold, therefore allow for more steps to be considered stall
FACTOR_STALL<-0.5 #recommended values between 0 and 1



#------DON'T CHANGE ANYTHING BELLOW THIS LINE-------------------



#Verification that all packages are installed and running

if(!require(tidyverse)){
  install.packages("tidyverse")
}
if(!require(viridis)){
  install.packages("viridis")
}
if(!require(viridisLite)){
  install.packages("viridisLite")
}
if(!require(corrplot)){
  install.packages("corrplot")
}

library(tidyverse)
library(viridis)
library(viridisLite)
library(corrplot)

#Load the FIJI/ImageJ data files into R/R-Studios

Cyt_folder<-getwd()
folder.name<-dir()
CytonemeMAXlength<-c()
for (w in 1:length(folder.name)){
  setwd(paste(Cyt_folder,"/",folder.name[w],sep = ""))
  files<-dir()
  files.name<-gsub(".txt", "", files)
  for (i in 1:length(files)) {
    assign(files.name[i], 
           read.csv(files[i], header = T, sep = "\t",
           ))
  }
  for(g in 1:length(files.name)){
    f<-get(files.name[g])
    CytonemeMAXlength<-max(c(max(f$Length, na.rm = T),CytonemeMAXlength))
  }
}

for (w in 1:length(folder.name)){
setwd(paste(Cyt_folder,"/",folder.name[w],sep = ""))
files<-dir()
files.name<-gsub(".txt", "", files)
for (i in 1:length(files)) {
  assign(files.name[i], 
         read.csv(files[i], header = T, sep = "\t",
         ))
}

Cyt_parameter_combined <- as.data.frame(matrix(data=NA,nrow = 0, ncol = 18))
colnames(Cyt_parameter_combined)<-c("Genotype",      #1
                           "Sample",                 #2
                           "CytonemeID",             #3
                           "Average_stepsize",       #4  
                           "Ave_FWD_stepsize",       #5
                           "Ave_REV_stepsize",       #6
                           "Ave_STALL_stepsize",     #7
                           "Total_Displacement",     #8
                           "Displa_FWD",             #9
                           "Displa_REV",             #10
                           "Displa_STALL",           #11
                           "Max_length",             #12
                           "Stall_time",             #13
                           "Average_stall_speed",    #14
                           "Average_speed",          #15
                           "FWD_Speed",              #16
                           "REV_Speed",              #17
                           "Lifetime")               #18

Cytoneme_combined_byID<-list()

for(g in 1:length(files.name)){
  f<-get(files.name[g])
  f$Step_displacement<-NA
  f$Length<- gsub(NaN,NA,f$Length)
  f[,"Length"]<-as.numeric(f[,"Length"])
  vec<-rep(NA, nrow(f))
  f$TimeRelative<-vec
  
  #Counting cytonemes
  n<-1
  for(c in 2:nrow(f)){
    if (f[c, 1] != f[(c-1),1]){ n<-n+1}
  }
  
  #Split data by cytoneme ID
  
  Cytoneme_byID<-split(f,f$ID) ## split data based on Track number collum 
  
  Cytoneme_Count_Accuracy<- length(Cytoneme_byID) == n
  
  #Calculate step displacement
  
  for(s in 1:length(Cytoneme_byID)){
    Cytoneme_byID[[s]][[1,8]]<-0
    for(d in 2:nrow(Cytoneme_byID[[s]])){
      Cytoneme_byID[[s]][[d,8]]<-as.numeric(Cytoneme_byID[[s]][[d,7]])-as.numeric(Cytoneme_byID[[s]][[(d-1),7]])
    }
  }
  
  teta<-c()
  for(s in 1:length(Cytoneme_byID)){
  if(any(duplicated(Cytoneme_byID[[s]]$Time))){  
    teta<-c(teta, Cytoneme_byID[[s]]$ID[1])
  }}  
  if (length(teta)>0){
  for(t in 1:length(teta)){
  Cytoneme_byID<- Cytoneme_byID[names(Cytoneme_byID) != teta[t]]
  }}
  
  
  for(s in 1:length(Cytoneme_byID)){
    timegap<-Cytoneme_byID[[s]]$Time[2]-Cytoneme_byID[[s]]$Time[1]
    maxrelative<-max(Cytoneme_byID[[s]]$Time, na.rm = T)-min(Cytoneme_byID[[s]]$Time, na.rm = T)
    Cytoneme_byID[[s]]$TimeRelative<-seq(0,maxrelative,timegap)
  }
  

  Cytoneme_combined_byID<-c(Cytoneme_combined_byID, Cytoneme_byID)

  
  #Calculate cytonemes parameters
  

  Cyt_parameter <- as.data.frame(matrix(data=NA,nrow = length(Cytoneme_byID), ncol = 18))
  colnames(Cyt_parameter)<-c("Genotype",               #1
                             "Sample",                 #2
                             "CytonemeID",             #3
                             "Average_stepsize",       #4  
                             "Ave_FWD_stepsize",       #5
                             "Ave_REV_stepsize",       #6
                             "Ave_STALL_stepsize",     #7
                             "Total_Displacement",     #8
                             "Displa_FWD",             #9
                             "Displa_REV",             #10
                             "Displa_STALL",           #11
                             "Max_length",             #12
                             "Stall_time",             #13
                             "Average_stall_speed",    #14
                             "Average_speed",          #15
                             "FWD_Speed",              #16
                             "REV_Speed",              #17
                             "Lifetime")               #18
  
  for(p in 1:length(Cytoneme_byID)){
    Cyt_parameter[p,1]<-folder.name[w] #Cytoneme Genotype
    
    Cyt_parameter[p,2]<-files.name[g] #Sample
    
    Cyt_parameter[p,3]<-Cytoneme_byID[[p]][[1,1]]#Cytoneme ID
    
    Cyt_parameter[p,4]<-mean(sqrt((Cytoneme_byID[[p]]$Step_displacement[2:nrow(Cytoneme_byID[[p]])])^2),na.rm = T) #Average_stepsize
    
    FWD<-Cytoneme_byID[[p]]$Step_displacement[2:nrow(Cytoneme_byID[[p]])]> FACTOR_STALL*Cyt_parameter[p,4] #Forward step
    if(all(FWD == F)){
      Cyt_parameter[p,5]<-NA
      Cyt_parameter[p,9]<-NA
      Cyt_parameter[p,16]<-NA
    }
    else{
      Cyt_parameter[p,5]<-mean(Cytoneme_byID[[p]]$Step_displacement[2:nrow(Cytoneme_byID[[p]])][FWD],na.rm = T) #Ave_FWD_stepsize
      Cyt_parameter[p,9]<-sum(sqrt((Cytoneme_byID[[p]]$Step_displacement[2:nrow(Cytoneme_byID[[p]])][FWD])^2),na.rm = T) #Displa_FWD
      Cyt_parameter[p,16]<-Cyt_parameter[p,9]/((sum(FWD,na.rm = T)+1)*timegap) #FWD_Speed
    }
    REV<-Cytoneme_byID[[p]]$Step_displacement[2:nrow(Cytoneme_byID[[p]])]<(-FACTOR_STALL*(Cyt_parameter[p,4]))
    if(all(REV == F)){
      Cyt_parameter[p,6]<-NA
      Cyt_parameter[p,10]<-NA
      Cyt_parameter[p,17]<-NA
    }
    else{
      Cyt_parameter[p,6]<-mean(Cytoneme_byID[[p]]$Step_displacement[2:nrow(Cytoneme_byID[[p]])][REV],na.rm = T) #Ave_REV_stepsize
      Cyt_parameter[p,10]<-sum(sqrt((Cytoneme_byID[[p]]$Step_displacement[2:nrow(Cytoneme_byID[[p]])][REV])^2),na.rm = T) #Displa_REV
      Cyt_parameter[p,17]<-Cyt_parameter[p,10]/((sum(REV,na.rm = T)+1)*timegap) #REV_Speed
    }
    STALL<- FWD == FALSE & REV == FALSE
    if(all(STALL == F)){
      Cyt_parameter[p,7]<-NA
      Cyt_parameter[p,11]<-NA
      Cyt_parameter[p,14]<-NA
    }
    else{
      Cyt_parameter[p,7]<-mean(sqrt((Cytoneme_byID[[p]]$Step_displacement[2:nrow(Cytoneme_byID[[p]])][STALL])^2),na.rm = T) #Ave_STALL_stepsize
      Cyt_parameter[p,11]<-sum(sqrt((Cytoneme_byID[[p]]$Step_displacement[2:nrow(Cytoneme_byID[[p]])][STALL])^2),na.rm = T) #Displa_STALL
      Cyt_parameter[p,14]<-Cyt_parameter[p,11]/((sum(STALL,na.rm = T)+1)*timegap) #Average_stall_speed
    }
    
    Cyt_parameter[p,8]<-sum(sqrt((Cytoneme_byID[[p]]$Step_displacement[2:nrow(Cytoneme_byID[[p]])])^2),na.rm = T) #Total_Displacement
    Cyt_parameter[p,12]<-max(Cytoneme_byID[[p]]$Length,na.rm = T)   #Max_length
    Cyt_parameter[p,13]<-sum(STALL,na.rm = T,na.rm = T)*timegap  #Stall_time
    Cyt_parameter[p,15]<-Cyt_parameter[p,8]/max(Cytoneme_byID[[p]]$Time,na.rm = T) #Average_speed
    Cyt_parameter[p,18]<-max(Cytoneme_byID[[p]]$Time,na.rm = T) #"Lifetime"
  }
  
  Cyt_parameter[,3:18] <- lapply(Cyt_parameter[,3:18],as.numeric)
  res <- cor(Cyt_parameter[,4:18], use = "complete.obs")
  png(height=1800, width=1800, file=paste("CorPlot_",files.name[g],".png"), type = "cairo")
  corrplot(res, type = "lower",tl.col = "black", tl.srt = 45,  tl.cex = 3, cl.cex = 3)
  dev.off()
  
  Data.stat<-data.frame(1:15,4,NA)
  colnames(Data.stat)<-c("Mean","Median", "s.d.")

  for(m in 4:18){Data.stat[m-3,1] <- colnames(Cyt_parameter[m])}
  for(m in 4:18){Data.stat[m-3,2] <- mean(Cyt_parameter[1:nrow(Cyt_parameter),m], na.rm = T )}
  for(m in 4:18){Data.stat[m-3,3] <- median(Cyt_parameter[1:nrow(Cyt_parameter),m], na.rm = T )}
  for(m in 4:18){Data.stat[m-3,4] <- sd(Cyt_parameter[1:nrow(Cyt_parameter),m], na.rm = T )}
  
  write_csv(Cyt_parameter , file= paste("Cyt_Para_",files.name[g]))
  write_csv(Data.stat , file= paste("Cyt_Stats_",files.name[g]))
  
  Cyt_parameter_combined<-rbind(Cyt_parameter_combined,Cyt_parameter)
  write_csv(Cyt_parameter_combined , file= paste("Cyt_Para_combined_",folder.name[w]))
  
  Cyt_parameter_combined[,3:18] <- lapply(Cyt_parameter_combined[,3:18],as.numeric)
  res <- cor(Cyt_parameter_combined[,4:18], use = "complete.obs")
  png(height=1800, width=1800, file=paste("CorPlot_combined_",folder.name[w],".png"), type = "cairo")
  corrplot(res, type = "lower",tl.col = "black", tl.srt = 45,  tl.cex = 3, cl.cex = 3)
  dev.off()
  
  Data.stat<-data.frame(1:15,3,NA)
  colnames(Data.stat)<-c("Mean","Median", "s.d.")
  
  for(m in 4:18){Data.stat[m-3,1] <- colnames(Cyt_parameter[m])}
  for(m in 4:18){Data.stat[m-3,2] <- mean(Cyt_parameter_combined[1:nrow(Cyt_parameter_combined),m], na.rm = T )}
  for(m in 4:18){Data.stat[m-3,3] <- median(Cyt_parameter_combined[1:nrow(Cyt_parameter_combined),m], na.rm = T )}
  for(m in 4:18){Data.stat[m-3,4] <- sd(Cyt_parameter_combined[1:nrow(Cyt_parameter_combined),m], na.rm = T )}
  
  write_csv(Data.stat , file= paste("Cyt_Stats_combined",folder.name[w]))
}

n.time<-c()
for(s in 1:length(Cytoneme_combined_byID)){
  n.time<-c(n.time, nrow(Cytoneme_combined_byID[[s]]))
}

max(n.time)

lengthMATRIX<-matrix(0, nrow = length(Cytoneme_combined_byID), ncol = max(n.time))

for(s in 1:length(Cytoneme_combined_byID)){
  lengthMATRIX[s,1:nrow(Cytoneme_combined_byID[[s]])]<-Cytoneme_combined_byID[[s]]$Length
}

rownames(lengthMATRIX)<- Cyt_parameter_combined[,3]

Cyt_dynamics<-data.frame(matrix(NA, nrow = 1,ncol = 7))
colnames(Cyt_dynamics)<- c("Genotype","ID","Time","Length","Step_displacement", "TimeRelative", "LogSeq")

for(s in 1:length(Cytoneme_combined_byID)){
  LogSeq<-rep(NA,nrow(Cytoneme_combined_byID[[s]]))
  preM<- Cytoneme_combined_byID[[s]][,c(1,6,7,8,9)]
  FWD<-Cytoneme_combined_byID[[s]]$Step_displacement[2:nrow(Cytoneme_combined_byID[[s]])]> FACTOR_STALL*Cyt_parameter_combined[s,4]
  REV<-Cytoneme_combined_byID[[s]]$Step_displacement[2:nrow(Cytoneme_combined_byID[[s]])]<(-FACTOR_STALL*(Cyt_parameter_combined[s,4]))
  STALL<- FWD == FALSE & REV == FALSE
  FWD[is.na(FWD)]<-FALSE
  REV[is.na(REV)]<-FALSE
  STALL[is.na(STALL)]<-FALSE
  for(y in 1:length(FWD)){
    if(FWD[y]==T){LogSeq[y+1]<-"FWD"}
    if(REV[y]==T){LogSeq[y+1]<-"REV"}
    if(STALL[y]==T){LogSeq[y+1]<-"STALL"}
  }
  Genotype<-rep(folder.name[w],nrow(Cytoneme_combined_byID[[s]]))
  preM<-cbind(Genotype,preM,LogSeq)
  Cyt_dynamics<-rbind(Cyt_dynamics,preM)
}

Cyt_dynamics<-Cyt_dynamics[2:nrow(Cyt_dynamics),]
Cyt_dynamics[,2:5] <- lapply(Cyt_dynamics[,2:5],as.numeric)


dist <- dist(lengthMATRIX, diag=TRUE)
hc <- hclust(dist)
cell_order<-rownames(lengthMATRIX)[hc$order]

#plotting Length change by cytoneme over time

max_limit_length <- ceiling(max(Cyt_dynamics[,4],na.rm = T))
tilew<- (Cyt_dynamics[2,3]-Cyt_dynamics[1,3])*1.1
Legth_Time_Cytoneme<-ggplot(Cyt_dynamics, aes(x = Cyt_dynamics[,3], 
                                              y = factor(Cyt_dynamics[,2], level = cell_order), 
                                              width=tilew, height=1)) + 
  geom_tile(aes(fill = Cyt_dynamics[,4]))+
  scale_fill_gradient(low="purple", high = "yellow", limits = c(0,CytonemeMAXlength))+
  ggtitle("Length change by cytoneme over time")+
  ylab("Cytoneme")+
  xlab("Time (min)")+
  labs(fill = "Length")+
  scale_y_discrete (labels = NULL, expand = c(0,0))+
  scale_x_continuous(limits = c(-1,max(Cyt_dynamics[,3], na.rm = T)+1), expand = c(0,0))+
  theme_classic()

ggsave(paste("Legth_Time_Cytoneme_",folder.name[w],".png"), device = "png", width = 15, height = 15, units = "cm", dpi=300,limitsize = F)

#plot Motion status by cytoneme over time
Status_Time_Cytoneme<-ggplot(Cyt_dynamics, aes(x = Cyt_dynamics[,3], 
                                               y = factor(Cyt_dynamics[,2], level = cell_order), 
                                               width=tilew, height=1)) + 
  geom_tile(aes(fill = Cyt_dynamics[,7]))+
  scale_fill_brewer(palette = "Accent")+
  ggtitle("Motion status by cytoneme over time")+
  ylab("Cytoneme")+
  xlab("Time (min)")+
  labs(fill = "Motion Status")+
  scale_y_discrete (labels = NULL, expand = c(0,0))+
  scale_x_continuous(limits = c(-1,max(Cyt_dynamics[,3], na.rm = T)+1), expand = c(0,0))+
  theme_classic()

ggsave(paste("Status_Time_Cytoneme_",folder.name[w],".png"), device = "png", width = 15, height = 15, units = "cm", dpi=300,limitsize = F)


#plot Steps length by cytonemes over time
Steps_Time_Cytoneme<-ggplot(Cyt_dynamics, aes(x = Cyt_dynamics[,3], 
                                              y = factor(Cyt_dynamics[,2], level = cell_order),
                                              width=tilew, height=1)) + 
  ggtitle("Steps length by cytonemes over time")+
  geom_tile(aes(fill = Cyt_dynamics[,5]))+
  scale_fill_gradient2(low="magenta", mid = "gray" ,high = "yellow", midpoint = 0, 
                       limits = c(floor(min(Cyt_dynamics[,5])),ceiling(max(Cyt_dynamics[,5]))))+
  ylab("Cytoneme")+
  xlab("Time (min)")+
  labs(fill = "Step length")+
  scale_y_discrete (labels = NULL, expand = c(0,0))+
  scale_x_continuous(limits = c(-1,max(Cyt_dynamics[,3], na.rm = T)+1), expand = c(0,0))+
  theme_classic()

ggsave(paste("Steps_Time_Cytoneme_",folder.name[w],".png"), device = "png", width = 15, height = 15, units = "cm", dpi=300,limitsize = F)


#plotting Length change by cytoneme over relative time


Legth_Time_Cytoneme<-ggplot(Cyt_dynamics, aes(x = Cyt_dynamics[,6], 
                                              y = factor(Cyt_dynamics[,2], level = cell_order),
                                              width=tilew, height=1)) + 
  geom_tile(aes(fill = Cyt_dynamics[,4]))+
  scale_fill_gradient(low="purple", high = "yellow", limits = c(0,CytonemeMAXlength))+
  ggtitle("Length change by cytoneme over time")+
  ylab("Cytoneme")+
  xlab("Time (min)")+
  labs(fill = "Length")+
  scale_y_discrete (labels = NULL, expand = c(0,0))+
  scale_x_continuous(limits = c(-1,max(Cyt_dynamics[,3], na.rm = T)+1), expand = c(0,0))+
  theme_classic()

ggsave(paste("Legth_RelativeTime_Cytoneme_",folder.name[w],".png"), device = "png", width = 15, height = 15, units = "cm", dpi=300,limitsize = F)

#plot Motion status by cytoneme over relative time
Status_Time_Cytoneme<-ggplot(Cyt_dynamics, aes(x = Cyt_dynamics[,6], 
                                               y = factor(Cyt_dynamics[,2], level = cell_order),
                                               width=tilew, height=1)) + 
  geom_tile(aes(fill = Cyt_dynamics[,7]))+
  scale_fill_brewer(palette = "Accent")+
  ggtitle("Motion status by cytoneme over time")+
  ylab("Cytoneme")+
  xlab("Time (min)")+
  labs(fill = "Motion Status")+
  scale_y_discrete (labels = NULL, expand = c(0,0))+
  scale_x_continuous(limits = c(-1,max(Cyt_dynamics[,3], na.rm = T)+1), expand = c(0,0))+
  theme_classic()

ggsave(paste("Status_RelativeTime_Cytoneme_",folder.name[w],".png"), device = "png", width = 15, height = 15, units = "cm", dpi=300,limitsize = F)

#plot Steps length by cytonemes over relative time
Steps_Time_Cytoneme<-ggplot(Cyt_dynamics, aes(x = Cyt_dynamics[,6], 
                                              y = factor(Cyt_dynamics[,2], level = cell_order), 
                                              width=tilew, height=1)) + 
  ggtitle("Steps length by cytonemes over time")+
  geom_tile(aes(fill = Cyt_dynamics[,5]))+
  scale_fill_gradient2(low="magenta", mid = "gray" ,high = "yellow", midpoint = 0, 
                       limits = c(floor(min(Cyt_dynamics[,5])),ceiling(max(Cyt_dynamics[,5]))))+
  ylab("Cytoneme")+
  xlab("Time (min)")+
  labs(fill = "Step length")+
  scale_y_discrete (labels = NULL, expand = c(0,0))+
  scale_x_continuous(limits = c(-1,max(Cyt_dynamics[,3], na.rm = T)+1), expand = c(0,0))+
  theme_classic()

ggsave(paste("Steps_RelativeTime_Cytoneme_",folder.name[w],".png"), device = "png", width = 15, height = 15, units = "cm", dpi=300,limitsize = F)


setwd(Cyt_folder)
write_csv(Cyt_parameter_combined , file= paste("Cyt_Para_combined_",folder.name[w]))

}

setwd(Cyt_folder)

#Load load Cyt_Para_combined files and combine then

files<-list.files(pattern = "Cyt_Para_combined_")
files.name<-gsub("Cyt_Para_combined_ ", "", files)
for (i in 1:length(files)) {
  assign(files.name[i], 
         read.csv(files[i], header = T, sep = ",",
         ))
}

# Combine the parameters files from all groups
Combine_ALLgroups<-data.frame()
for(g in 1:length(files.name)){
  f<-get(files.name[g])
  Combine_ALLgroups<-rbind(Combine_ALLgroups,f)
}

write_csv(Combine_ALLgroups , file= "Combine_ALLgroups")

dir.create(paste(Cyt_folder,"/Graphs", sep = ""))
setwd(paste(Cyt_folder,"/Graphs", sep = ""))

for (g in 4:17){
  plot_compare <- ggplot(Combine_ALLgroups, aes(x=Combine_ALLgroups[,18], 
                                                y=Combine_ALLgroups[,g]))+
    ylab(colnames(Combine_ALLgroups)[g])+
    xlab(expression(paste("Lifetime (min)" )))+
    # ylim(0,30)+
    geom_point(aes(color=Genotype),size=3, alpha=0.75, show.legend = T)+
    geom_smooth(method=loess,linetype=2, aes(color=Genotype), alpha=0.0,size=0.5)+
    scale_color_viridis_d()+
    theme_classic()
  
  ggsave(paste(colnames(Combine_ALLgroups)[g]," vs. Lifetime",".png"), 
         device = "png", width = 15, height = 15, units = "cm", dpi=300,limitsize = F)
}

for (g in c(4,5,6,7,9,10,11,12,13,14,15,16,17,18)){
  plot_compare <- ggplot(Combine_ALLgroups, aes(x=Combine_ALLgroups[,8], 
                                                y=Combine_ALLgroups[,g]))+
    ylab(colnames(Combine_ALLgroups)[g])+
    xlab(expression(paste("Total displacement (",mu,"m)" )))+
    # ylim(0,30)+
    geom_point(aes(color=Genotype),size=3, alpha=0.75, show.legend = T)+
    geom_smooth(method=loess,linetype=2, aes(color=Genotype), alpha=0.0,size=0.5)+
    scale_color_viridis_d()+
    theme_classic()
  
  ggsave(paste(colnames(Combine_ALLgroups)[g]," vs. Total displacement",".png"), 
         device = "png", width = 15, height = 15, units = "cm", dpi=300,limitsize = F)
}

for (g in 4:18){
  plot_compare <- ggplot(Combine_ALLgroups, aes(y=Combine_ALLgroups[,g], 
                                                x=factor(Genotype, levels = Group_order), fill= Genotype, na.rm = TRUE))+
    xlab("Groups")+
    ylab(colnames(Combine_ALLgroups)[g])+
    geom_violin(size=0.5, alpha=0.4,show.legend = F)+
    geom_point(shape = 21,size=2,show.legend = F, position = position_jitterdodge(), color="black",alpha=0.8)+
    scale_fill_viridis_d()+
    theme_classic()
  
  ggsave(paste(colnames(Combine_ALLgroups)[g],".png"), 
         device = "png", width = 15, height = 15, units = "cm", dpi=300,limitsize = T)
}

rm(list = ls())