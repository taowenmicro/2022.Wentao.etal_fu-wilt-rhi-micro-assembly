library(picante)
library(ape)
library(vegan)
library(FSA)
library(eulerr)
library(grid)
library(gridExtra)
require(minpack.lm)
require(Hmisc)
require(stats4)
library(parallel)


ps16s = readRDS("./data/Sample_Fwilt_4/16s/ps_16s.rds")
ps16s
map = sample_data(ps16s)

map[8,2] = "WD"

sample_data(ps16s) = map

tree = read.tree("./data/Sample_Fwilt_4/16s/97_otus.tree")
phy_tree(ps16s) = tree



ps16s1 = filter_taxa(ps16s, function(x) sum(x ) > 500, TRUE)

pssec <- subset_samples(ps16s1,!Group %in% c("BH","BD"));pssec
pssec <- subset_samples(pssec,Group2 %in% c("D"));pssec
pssec = filter_taxa(pssec, function(x) sum(x ) > 0 , TRUE) 


# # 建立第二级目录
result_path <- paste("./result16s_4",sep = "")
dir.create(res1path)



# allid = c("Anaerolinea","Sphingomonas")
allid = c("Anaerolinea")
allid0 = c("Anaerolinea", "Anaeromyxobacter", "Bdellovibrio", "Conexibacter", "Flavobacterium", "Gemmatimonas","Janibacter")

allid =allid0[7]

phypath = paste(result_path,"/result_rm_",allid,sep = "")
dir.create(phypath)

#--首先去除根瘤菌做全套系统发育#------
library(tidyverse)
library(ggClusterNet)
ps16s <- pssec %>%
  subset_taxa(
    !Genus  == allid
  )
ps16s





# #--中性模型#-----
# source("G:\\Shared_Folder\\Function_local\\R_function\\micro\\phylo_Micro\\neutralModel.R")
# result = neutralModel(ps = ps16s,group  = "Group",ncol =5,nrow  = 1)
# #--合并图表
# p1 =  result[[1]]
# p1
# # #--分开的图表--存储在一个list中
# # plist = result[[2]]
# # #-提取单个的图表
# # plist[[1]]
# 
# 
# FileName <- paste(phypath,"1_neutral_modelCul", ".pdf", sep = "")
# ggsave(FileName, p1,width = 18,height = 6)




#--BNTI#-
source("G:\\Shared_Folder\\Function_local\\R_function\\micro\\phylo_Micro\\bNTICul.R")
result = bNTICul(ps = ps16s ,group  = "Group",num = 10,thread = 1)
bNTI = result[[1]]
head(bNTI)
filename = paste(phypath,"/4_bNTI.csv",sep = "")
write.csv(bNTI, filename)
#--计算RCbray#------------
source("G:\\Shared_Folder\\Function_local\\R_function\\micro\\phylo_Micro\\RCbary.R")
result = RCbary(ps = ps16s,group  = "Group",num = 10,thread = 1)
RCbary = result[[1]]
head(RCbary)
filename = paste(phypath,"/5_RCb.csv",sep = "")
write.csv(RCbary,filename)

#--BetaNTI和RCbray联合出图#---------
source("G:\\Shared_Folder\\Function_local\\R_function\\micro\\phylo_Micro\\bNTIRCPlot.R")
bNTI = read.csv(paste(phypath,"/4_bNTI.csv",sep = ""),row.names = 1)
head(bNTI)
# RCbray 数据读入，修改列名
RCb = read.csv(paste(phypath,"/5_RCb.csv",sep = ""),row.names = 1) %>%
  dplyr::mutate(Sample_1 = Site2, Sample_2 = Site1)
head(RCb)

head(map)


result = bNTIRCPlot(ps = ps16s ,RCb  = RCb,bNTI = bNTI,group  = "Group")

#--bNTI出图片
p3 <- result[[1]] 
p3

#RCbary可视化
p4 <- result[[2]] 
p4
#组合图片BNTI，RCbray
p5 <- result[[3]]
p5
plotdata = result[[4]]
head(plotdata)

filename = paste(phypath,"/6_bNTI_RCbray.csv",sep = "")
write.csv(plotdata,filename)

FileName <- paste(phypath,"6_bNTI", ".pdf", sep = "")
ggsave(FileName, p3,width =8,height = 6)

FileName <- paste(phypath,"6_RCbary", ".pdf", sep = "")
ggsave(FileName, p4,width = 6,height = 6)

FileName <- paste(phypath,"6_BNTI_RCbray", ".pdf", sep = "")
ggsave(FileName, p5,width = 12,height = 8)






#----随机去除两个属#----------

tax = vegan_tax(pssec) %>% as.data.frame()
id <- tax$Genus %>% unique()
id
i = 7

ps16sr <- pssec %>%
  subset_taxa(
    #Kingdom == "Bacteria" &
    !Genus  == id[c(i)]
    # Species %in%c("Fusarium_oxysporum","Fusarium_keratoplasticum")
    # row.names(tax_table(ps1_rela ))%in%c("SH010924.07FU_KF986690_reps_singleton","SH020983.07FU_JN235282_refs")
  )
ps16sr





phpath = paste(result_path,"/result_rn_g__anyone/",sep = "")
dir.create(phpath)
phypath = paste(phpath,"/",id[c(i)],sep = "")
dir.create(phypath)

#--BNTI#-
source("G:\\Shared_Folder\\Function_local\\R_function\\micro\\phylo_Micro\\bNTICul.R")
result = bNTICul(ps = ps16sr,group  = "Group",num = 10,thread = 1)
bNTI = result[[1]]
head(bNTI)
filename = paste(phypath,"/4_bNTI.csv",sep = "")
write.csv(bNTI, filename)
#--计算RCbray#------------
source("G:\\Shared_Folder\\Function_local\\R_function\\micro\\phylo_Micro\\RCbary.R")
result = RCbary(ps = ps16sr,group  = "Group",num = 10,thread = 1)
RCbary = result[[1]]
head(RCbary)
filename = paste(phypath,"/5_RCb.csv",sep = "")
write.csv(RCbary,filename)

#--BetaNTI和RCbray联合出图#---------
source("G:\\Shared_Folder\\Function_local\\R_function\\micro\\phylo_Micro\\bNTIRCPlot.R")
bNTI = read.csv(paste(phypath,"/4_bNTI.csv",sep = ""),row.names = 1)
head(bNTI)
# RCbray 数据读入，修改列名
RCb = read.csv(paste(phypath,"/5_RCb.csv",sep = ""),row.names = 1) %>%
  dplyr::mutate(Sample_1 = Site2, Sample_2 = Site1)
head(RCb)
head(map)
result = bNTIRCPlot(ps = ps16s ,RCb  = RCb,bNTI = bNTI,group  = "Group")
#--bNTI出图片
p3 <- result[[1]]
p3

#RCbary可视化
p4 <- result[[2]] 
p4
#组合图片BNTI，RCbray
p5 <- result[[3]]
p5
plotdata = result[[4]]
head(plotdata)

filename = paste(phypath,"/6_bNTI_RCbray.csv",sep = "")
write.csv(plotdata,filename)

FileName <- paste(phypath,"/6_bNTI", ".pdf", sep = "")
ggsave(FileName, p3,width =8,height = 6)

FileName <- paste(phypath,"/6_RCbary", ".pdf", sep = "")
ggsave(FileName, p4,width = 6,height = 6)

FileName <- paste(phypath,"/6_BNTI_RCbray", ".pdf", sep = "")
ggsave(FileName, p5,width = 12,height = 8)

#----

