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
library(tidyverse)
library(viridis)
library(viridisLite)

#Load load Cyt_Para_combined files and combine then

folder.name<-basename(getwd())
files<-dir()
files.name<-gsub("Cyt_Para_combined_ ", "", files)
for (i in 1:length(files)) {
  assign(files.name[i], 
         read.csv(files[i], header = T, sep = ",",
         ))
}

#Define order of the samples
files.name<-c("Wild type", "HhGal4-ttvi") #list the name of the groups inside c(), as writen in the folder of origin, between "" (i.e.: "Control"), separare the diferent groups by commaa (i.e.: "Control", "Treatment") 

Combine_ALLgroups<-data.frame()
for(g in 1:length(files.name)){
  f<-get(files.name[g])
  Combine_ALLgroups<-rbind(Combine_ALLgroups,f)
}


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
                         x=factor(Genotype, levels = files.name), fill= Genotype, na.rm = TRUE))+
    xlab("Groups")+
    ylab(colnames(Combine_ALLgroups)[g])+
    geom_violin(size=0.5, alpha=0.4,show.legend = F)+
    geom_point(shape = 21,size=2,show.legend = F, position = position_jitterdodge(), color="black",alpha=0.8)+
    scale_fill_viridis_d()+
    theme_classic()
  
  ggsave(paste(colnames(Combine_ALLgroups)[g],".png"), 
         device = "png", width = 15, height = 15, units = "cm", dpi=300,limitsize = T)
}