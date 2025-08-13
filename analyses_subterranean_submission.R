
# R script to replicate analyses for "Capturing the unseen: a novel, low-cost method for stratified subterranean sampling of soil invertebrates in drylands"
# Last updated: July 06 2025

#start by clearing everything in your local R memory
rm(list = ls())

#load necessary packages

library(lubridate)
library(mgcv)
library(dplyr)
library(ggplot2)
library(viridis)

library(glmmTMB) 

### LOAD AND PREPARE NECESSARY DATA

# set working directory

setwd("YOUR_LOCAL_PATH")

# macrohabitats: 
hab=read.csv("subterraneanHabitat_202502.csv")

hab$site=substr(hab$Site, 2, 4) # get site info from reference name

unique(hab$site) 
length(unique(hab$site))# should be 10 

hab$line=substr(hab$Site, 5, 7) # get bucket reference info from reference name

unique(hab$line) # should be "L05" and "L01"


df=read.csv("TableS1_forAnalysis.csv")

df=df[,-1]# remove the "ref" column

# Change data format
df$date=as.Date(df$Timestamp,format="%d/%m/%y %H:%M")

head(df$date) # make sure your see yyyy-mm-dd format

#order data dy date
df=df[order(df$date),]

#remove data where we have no date
df=df[!is.na(df$date),]

df$month=month(df$date) # get month from date

df$year=year(df$date) # get year from date

df$site=substr(df$Ref, 2, 4) # get site info from reference name

unique(df$site)

length(unique(df$site))# should be 10 

df$line=substr(df$Ref, 5, 7) # get bucket reference info from reference name

unique(df$line) # should be "L05" and "L01"

df$layer=substr(df$Ref, 8, 10) # get subterranean layer from Ref

unique(df$layer) # should be "L01", "L02", and "L03"

unique(df$Method) # should be 2 things

### Clean diversity data:

# Remove unidentified Classes

df=droplevels(df[!df$Class=="Unidentified",])

# Create different taxa
df$ID=paste(df$Class,	df$Order,	df$Family,	df$Subfamily,	df$Tribe,	df$Genera,	df$Species)

table(df$ID)


# Do analyses from Oct 2024 onward, where traps were placed at 10 sites

# 1. Do richness = number of species caught

df_sub=df[!is.na(df$Number_caught)&df$date>"2024-10-01",]

unique(df_sub$line[df_sub$Method=="double-stratified subterranean"]) # should be L01

unique(df_sub$line[df_sub$Method=="three-stratified subterranean" ]) # should be L05

# create sampling occasions (sampling month):

df_sub$occ=NA

df_sub$occ[df_sub$month%in%c(10,11)]=1
df_sub$occ[df_sub$month%in%c(12)]=2
df_sub$occ[df_sub$month%in%c(1)]=3
df_sub$occ[df_sub$month%in%c(2)]=4
df_sub$occ[df_sub$month%in%c(3)]=5


# If we split by trap type and layer, we get lots of 0 inflation (uncomment to check)

# rich=aggregate(ID~layer+line+site+occ,function(x) length(x),drop=F,data=df_sub)

# 1.1. So, we first get number of different species caught do by line + site, month, and land-use

rich=aggregate(ID~line+site+occ,function(x) length(unique(x)),drop=F,data=df_sub)

rich$land_use=left_join(rich,hab[!duplicated(hab$site),],by="site")$land_use

rich$hab=left_join(rich,hab,by=c("site","line"))$Sand

rich$ID[is.na(rich$ID)]=0

rich$occ=factor(rich$occ)

### Do some visuals

ggplot(rich,aes(x=ID,group=occ))+
  geom_density(aes(colour=occ))

ggplot(rich,aes(x=ID,group=site))+
  geom_density(aes(colour=site))

ggplot(rich,aes(x=ID,group=land_use))+
  geom_density(aes(colour=land_use))

ggplot(rich,aes(x=ID,group=hab))+
  geom_density(aes(colour=hab))

# Summary: no clear trends

### MODEL RICHNESS

# Poisson model

### Test for zero inflation

poisson0a <- glmmTMB(ID~1+(1|site), 
                    data=rich,
                    ziformula=~0,
                    family=poisson)

poisson0b <- glmmTMB(ID~1+(1|site), 
                    data=rich,
                    ziformula=~1,
                    family=poisson)

AIC(poisson0a,poisson0b) # with 0 inflation

# Predictors on richness

poisson1 <- glmmTMB(ID~land_use+(1|site), 
                    data=rich,
                    ziformula=~1,
                    family=poisson)

summary(poisson1)

poisson2 <- glmmTMB(ID~occ+(1|site), 
                    data=rich,
                    ziformula=~1,
                    family=poisson)

summary(poisson2)

poisson3 <- glmmTMB(ID~hab+(1|site), 
                    data=rich,
                    ziformula=~1,
                    family=poisson)


summary(poisson3)

AIC(poisson0b,poisson1,poisson2,poisson3) # differences between occasions

poisson2a <- glmmTMB(ID~occ+(1|site), 
                    data=rich,
                    ziformula=~occ,
                    family=poisson)

poisson2b <- glmmTMB(ID~occ+(1|site), 
                     data=rich,
                     ziformula=~hab,
                     family=poisson)

poisson2c <- glmmTMB(ID~occ+(1|site), 
                     data=rich,
                     ziformula=~land_use,
                     family=poisson)


AIC(poisson0b,poisson2,poisson2a,poisson2b,poisson2c)

summary(poisson2c)
summary(poisson2)

# Best model: Seasonal differences 

# Negative Binomial model

negbin0a <- glmmTMB(ID~1+(1|site), 
                   data=rich,
                   ziformula=~0,
                   family=nbinom1)

negbin0b <- glmmTMB(ID~1+(1|site), 
                    data=rich,
                    ziformula=~1,
                    family=nbinom1)

AIC(negbin0a,poisson0b) # Poisson model better

summary(negbin0a)

negbin1 <- glmmTMB(ID~occ+(1|site), 
                    data=rich,
                    ziformula=~0,
                    family=nbinom1)

summary(negbin1)

negbin2 <- glmmTMB(ID~occ+(1|site), 
                   data=rich,
                   ziformula=~1,
                   family=nbinom2)

AIC(poisson2,negbin0a,negbin1,negbin2)

summary(poisson2)

### Summary: poisson2 is the best and show increasing richness in time

rich$samp=rich$occ
levels(rich$samp)=c("Nov", "Dec", "Jan", "Feb", "Mar")
# Plot: 

ggplot(rich,aes(x=samp,y=ID))+
  geom_boxplot(aes(colour=samp),outliers = F)+
  scale_colour_viridis(discrete = TRUE, alpha=0.8) +
  geom_jitter(color="grey", size=2, alpha=0.9) +
  scale_y_continuous(limits = c(0, 5), oob = scales::squish)+
  xlab("Month of sampling")+
  ylab("Richness (# family & genus taxa)")+
  theme_bw(base_size = 20)+
  theme(panel.grid = element_blank())+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+
  theme(strip.background =element_blank(),
        strip.text =element_text(size=20),
        legend.position = "none")

ggsave(filename ="richness.pdf",width = 7,height = 5)

#### 2. Differences among trap types and layers in richness

rich=aggregate(ID~layer+line+site+occ,function(x) length(unique(x)),drop=F,data=df_sub)

rich=rich[!(rich$line=="L01"&rich$layer=="L03"),] # this combination doesn't exist

rich$land_use=left_join(rich,hab[!duplicated(hab$site),],by="site")$land_use

rich$hab=left_join(rich,hab,by=c("site","line"))$Sand

rich$ID[is.na(rich$ID)]=0

rich$occ=factor(rich$occ)

##### Subset data to only show instances where at least 1 capture occurred

rich$analyze="no"

# for each occasion, site, and line, if there was at least one cpature, analyze = "yes"

for(a in 1:length(unique(rich$occ))){
  
  for(b in 1:length(unique(rich$site))){
    
    for(c in 1:length(unique(rich$line))){
      
      if(any(rich$ID[rich$occ==unique(rich$occ)[a]&rich$site==unique(rich$site)[b]&rich$line==unique(rich$line)[c]]>0)){
        
        
        rich$analyze[rich$occ==unique(rich$occ)[a]&rich$site==unique(rich$site)[b]&rich$line==unique(rich$line)[c]] <- "yes"
        
      }
      
    }
  }
  

}

### Do some models

### DOUBLE STRATIFIED

# Poisson model

poisson0a <- glmmTMB(ID~1+(1|site), 
                    data=rich[rich$line=="L01"&rich$analyze=="yes",],
                    ziformula=~0,
                    family=poisson)

poisson0b <- glmmTMB(ID~1+(1|site), 
                     data=rich[rich$line=="L01"&rich$analyze=="yes",],
                     ziformula=~0,
                     family=poisson)

AIC(poisson0a,poisson0b)


poisson1 <- glmmTMB(ID~layer+(1|site), 
                    data=rich[rich$line=="L01"&rich$analyze=="yes",],
                    ziformula=~0,
                    family=poisson)

summary(poisson1)


poisson2 <- glmmTMB(ID~layer+(1|site), 
                    data=rich[rich$line=="L01"&rich$analyze=="yes",],
                    ziformula=~1,
                    family=poisson)

summary(poisson2)

poisson3 <- glmmTMB(ID~1+(1|site), 
                    data=rich[rich$line=="L01"&rich$analyze=="yes",],
                    ziformula=~layer,
                    family=poisson)

summary(poisson3)

AIC(poisson0a,poisson1,poisson2,poisson3)

### THREE STRATIFIED

poisson4 <- glmmTMB(ID~layer+(1|site), 
                    data=rich[rich$line=="L05"&rich$analyze=="yes",],
                    ziformula=~0,
                    family=poisson)

summary(poisson4)


poisson5 <- glmmTMB(ID~layer+(1|site), 
                    data=rich[rich$line=="L05"&rich$analyze=="yes",],
                    ziformula=~1,
                    family=poisson)

summary(poisson5)

poisson6 <- glmmTMB(ID~1+(1|site), 
                    data=rich[rich$line=="L05"&rich$analyze=="yes",],
                    ziformula=~layer,
                    family=poisson)

summary(poisson6)

AIC(poisson4,poisson5,poisson6)


# Negative Binomial model
negbin1 <- glmmTMB(ID~layer+(1|site), 
                   data=rich[rich$line=="L01"&rich$analyze=="yes",],
                   ziformula=~0,
                   family=nbinom1)


summary(negbin1)

negbin3 <- glmmTMB(ID~layer+(1|site), 
                   data=rich[rich$line=="L05"&rich$analyze=="yes",],
                   ziformula=~0,
                   family=nbinom1)


summary(negbin3)

negbin4 <- glmmTMB(ID~layer+(1|site), 
                   data=rich[rich$line=="L05"&rich$analyze=="yes",],
                   ziformula=~1,
                   family=nbinom1)

summary(negbin4)

# Summary: no differences between layers in richness

# Plot:

# double statified: 

a=ggplot(rich[rich$line=="L01"&rich$analyze=="yes",],aes(x=layer,y=ID))+
  geom_boxplot(aes(colour=layer),outliers = F)+
  scale_colour_viridis(discrete = TRUE, alpha=0.8) +
  geom_jitter(color="grey", size=2, alpha=0.9) +
  scale_y_continuous(limits = c(0, 5), oob = scales::squish)+
  xlab("Sampling layer")+
  ylab("Richness (# family & genus taxa)")+
  ggtitle("(a) Double-stratified design")+
  theme_bw(base_size = 20)+
  theme(panel.grid = element_blank())+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+
  theme(strip.background =element_blank(),
        strip.text =element_text(size=20),
        legend.position = "none")

b=ggplot(rich[rich$line=="L05"&rich$analyze=="yes",],aes(x=layer,y=ID))+
  geom_boxplot(aes(colour=layer), outliers = F)+
  scale_colour_viridis(discrete = TRUE, alpha=0.8) +
  geom_jitter(color="grey", size=2, alpha=0.9) +
  scale_y_continuous(limits = c(0, 5), oob = scales::squish)+
  xlab("Sampling layer")+
  ylab("Richness (# family & genus taxa)")+
  ggtitle("(b) Three-stratified design")+
  theme_bw(base_size = 20)+
  theme(panel.grid = element_blank())+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+
  theme(strip.background =element_blank(),
        strip.text =element_text(size=20),
        legend.position = "none")

library(patchwork)

c=a+b

ggsave(c,filename ="richness_layer.pdf",width = 10,height = 6)


##### 3. BIOMASS ANALYSES


library(readxl)

path="YOUR_PATH/biomass_data"

nfiles=list.files(path=path)

df.m=NULL

### Join all files

# NOVEMBER (and end of October)
temp=read_excel(paste0(path,nfiles[grep("08_11_2024",nfiles)]),sheet = 1)
unique(temp$Ref.) # should be 25
temp$occ=1

df.m=rbind(df.m,temp)

temp=read_excel(paste0(path,nfiles[grep("25_10_2024",nfiles)]),sheet = 1)
unique(temp$Ref.) # should be 25

temp$occ=1

df.m=rbind(df.m,temp)

# DECEMBER
temp=read_excel(paste0(path,nfiles[grep("21_12_2024",nfiles)]),sheet = 1)
unique(temp$Ref.) # should be 25

temp$occ=2

df.m=rbind(df.m,temp)

temp=read_excel(paste0(path,nfiles[grep("27_12_2024",nfiles)]),sheet = 1)
unique(temp$Ref.) # should be 25

temp$occ=2

df.m=rbind(df.m,temp)

# JANUARY
temp=read_excel(paste0(path,nfiles[grep("30_01_2025",nfiles)]),sheet = 1)
unique(temp$Ref.) # should be 25

temp$occ=3

df.m=rbind(df.m,temp[,c(1:4,6)])

temp=read_excel(paste0(path,nfiles[grep("28_01_2025",nfiles)]),sheet = 1)
unique(temp$Ref.) # should be 25

temp$occ=3

df.m=rbind(df.m,temp)

# FEBRUARY
temp=read_excel(paste0(path,nfiles[grep("20_02_2025",nfiles)]),sheet = 1)
unique(temp$Ref.) 

temp$occ=4

df.m=rbind(df.m,temp)

temp=read_excel(paste0(path,nfiles[grep("27_02_2025",nfiles)]),sheet = 1)
unique(temp$Ref.) 

temp$occ=4

df.m=rbind(df.m,temp)

# MARCH
temp=read_excel(paste0(path,nfiles[grep("26_03_2025",nfiles)]),sheet = 1)
unique(temp$Ref.) 

temp$occ=5

df.m=rbind(df.m,temp)

temp=read_excel(paste0(path,nfiles[grep("31_03_2025",nfiles)]),sheet = 1)
unique(temp$Ref.) 

temp$occ=5

df.m=rbind(df.m,temp)

colnames(df.m)=c("date","Ref","Weight","Number","occ")

#order data dy date
df.m=df.m[order(df.m$date),]

#remove data where we have no date
df.m=df.m[!is.na(df.m$date),]

df.m$month=month(df.m$date) # get month from date

df.m$year=year(df.m$date) # get year from date

df.m$site=substr(df.m$Ref, 2, 4) # get site info from reference name

length(unique(df.m$site)) # should be 10 

df.m$line=substr(df.m$Ref, 5, 7) # get bucket reference info from reference name

unique(df.m$line) # should be "L05" and "L01"


df.m$layer=substr(df.m$Ref, 8, 10) # get subterranean layer from Ref

unique(df.m$layer) # should be "L01", "L02", and "L03"

bio=aggregate(cbind(Weight,Number)~line+site+occ,function(x) sum(x,na.rm=T),data=df.m)

bio$land_use=left_join(bio,hab[!duplicated(hab$site),],by="site")$land_use

bio$hab=left_join(bio,hab,by=c("site","line"))$Sand

bio$Weight[is.na(bio$Weight)]=0

bio$occ=factor(bio$occ)

### Do some visuals

ggplot(bio,aes(x=Weight,group=occ))+
  geom_density(aes(colour=occ))

ggplot(bio,aes(x=Weight,group=site))+
  geom_density(aes(colour=site))

ggplot(bio,aes(x=Weight,group=land_use))+
  geom_density(aes(colour=land_use))

ggplot(bio,aes(x=Weight,group=hab))+
  geom_density(aes(colour=hab))

# Summary: no clear trend

### Do some models

hist(bio$Weight)

hist(bio$Number)


# BIOMASSS

m0a <- glmmTMB(Weight~1+(1|site), 
                data=bio,
                ziformula=~0,
                family=tweedie)

m0b <- glmmTMB(Weight~1+(1|site), 
                 data=bio,
                 ziformula=~1,
                 family=tweedie)

AIC(m0a,m0b)

m1 <- glmmTMB(Weight~land_use+(1|site), 
                data=bio,
                ziformula=~0,
                family=tweedie)


AIC(m0a,m1)

summary(m1)

m2 <- glmmTMB(Weight~occ+(1|site), 
                data=bio,
                ziformula=~0,
                family=tweedie)


AIC(m0a,m2)

summary(m2)

m3 <- glmmTMB(Weight~hab+(1|site), 
                data=bio,
                ziformula=~0,
                family=tweedie)


AIC(m0a,m3)

summary(m3)

m3b <- glmmTMB(Weight~hab+(1|site), 
              data=bio,
              ziformula=~1,
              family=tweedie)

AIC(m0a,m3)

m4 <- glmmTMB(Weight~line+(1|site), 
              data=bio,
              ziformula=~0,
              family=tweedie)


AIC(m0a,m4)

summary(m4)


# Plot: 

bio$samp=bio$occ
levels(bio$samp)=c("Nov", "Dec", "Jan", "Feb", "Mar")

ggplot(bio,aes(x=samp,y=Weight))+
  geom_boxplot(aes(colour=samp),outliers = F)+
  scale_colour_viridis(discrete = TRUE, alpha=0.8) +
  geom_jitter(color="grey", size=2, alpha=0.9) +
  scale_y_continuous(limits = c(0, 0.25), oob = scales::squish)+
  xlab("Month of sampling")+
  ylab("Total biomass (g)")+
  theme_bw(base_size = 20)+
  theme(panel.grid = element_blank())+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+
  theme(strip.background =element_blank(),
        strip.text =element_text(size=20),
        legend.position = "none")

ggsave(filename ="biomass.pdf",width = 7,height = 5)


#### 2. Differences among trap types and layers in richness

bio=aggregate(cbind(Weight,Number)~layer+line+site+occ,function(x) sum(x,na.rm=T),data=df.m)

bio$Weight[is.na(bio$Weight)]=0

bio$occ=factor(bio$occ)

bio=bio[!(bio$line=="L01"&bio$layer=="L03"),] # this combination doesn't exist

bio$land_use=left_join(bio,hab[!duplicated(hab$site),],by="site")$land_use

bio$hab=left_join(bio,hab,by=c("site","line"))$Sand

bio$occ=factor(bio$occ)

##### Subset data to only show instances where at least 1 capture occurred

bio$analyze="no"

# for each occasion, site, and line, if there was at least one cpature, analyze = "yes"

for(a in 1:length(unique(bio$occ))){
  
  for(b in 1:length(unique(bio$site))){
    
    for(c in 1:length(unique(bio$line))){
      
      if(any(bio$Weight[bio$occ==unique(bio$occ)[a]&bio$site==unique(bio$site)[b]&bio$line==unique(bio$line)[c]]>0)){
        
        
        bio$analyze[bio$occ==unique(bio$occ)[a]&bio$site==unique(bio$site)[b]&bio$line==unique(bio$line)[c]] <- "yes"
        
      }
      
    }
  }
  
  
}

### Do some models

### DOUBLE STRATIFIED

tweedie0a <- glmmTMB(Weight~1+(1|site), 
                     data=bio[bio$line=="L01"&bio$analyze=="yes",],
                     ziformula=~0,
                     family=tweedie)

tweedie0b <- glmmTMB(Weight~1+(1|site), 
                     data=bio[bio$line=="L01"&bio$analyze=="yes",],
                     ziformula=~1,
                     family=tweedie)

AIC(tweedie0a,tweedie0b)


tweedie1 <- glmmTMB(Weight~layer+(1|site), 
                    data=bio[bio$line=="L01"&bio$analyze=="yes",],
                    ziformula=~0,
                    family=tweedie)

summary(tweedie1)


tweedie2 <- glmmTMB(Weight~layer+(1|site), 
                    data=bio[bio$line=="L01"&bio$analyze=="yes",],
                    ziformula=~1,
                    family=tweedie)

summary(tweedie2)

tweedie3 <- glmmTMB(Weight~1+(1|site), 
                    data=bio[bio$line=="L01"&bio$analyze=="yes",],
                    ziformula=~layer,
                    family=tweedie)

summary(tweedie3)

AIC(tweedie0a,tweedie1,tweedie2,tweedie3) # tweedie3 has lowest AIC, but estimates suggests that it is overfitted

### THREE STRATIFIED

tweedie4 <- glmmTMB(Weight~layer+(1|site), 
                    data=bio[bio$line=="L05"&bio$analyze=="yes",],
                    ziformula=~0,
                    family=tweedie)

summary(tweedie4)


tweedie5 <- glmmTMB(Weight~layer+(1|site), 
                    data=bio[bio$line=="L05"&bio$analyze=="yes",],
                    ziformula=~1,
                    family=tweedie)

summary(tweedie5)

tweedie6 <- glmmTMB(Weight~1+(1|site), 
                    data=bio[bio$line=="L05"&bio$analyze=="yes",],
                    ziformula=~layer,
                    family=tweedie)

summary(tweedie6)

AIC(tweedie4,tweedie5)


# Summary: no differences between layers in biomass

# Plot:

# double stratified: 

a=ggplot(bio[bio$line=="L01"&bio$analyze=="yes",],aes(x=layer,y=Weight))+
  geom_boxplot(aes(colour=layer),outliers=F)+
  scale_colour_viridis(discrete = TRUE, alpha=0.8) +
  geom_jitter(color="grey", size=2, alpha=0.9) +
  scale_y_continuous(limits = c(0, 1), oob = scales::squish)+
  xlab("Sampling layer")+
  ylab("Biomass (g)")+
  ggtitle("(a) Double-stratified design")+
  theme_bw(base_size = 20)+
  theme(panel.grid = element_blank())+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+
  theme(strip.background =element_blank(),
        strip.text =element_text(size=20),
        legend.position = "none")

b=ggplot(bio[bio$line=="L05"&bio$analyze=="yes",],aes(x=layer,y=Weight))+
  geom_boxplot(aes(colour=layer),outliers=F)+
  scale_colour_viridis(discrete = TRUE, alpha=0.8) +
  geom_jitter(color="grey", size=2, alpha=0.9) +
  scale_y_continuous(limits = c(0, 1), oob = scales::squish)+
  xlab("Sampling layer")+
  ylab("Biomass (g)")+
  ggtitle("(b) Three-stratified design")+
  theme_bw(base_size = 20)+
  theme(panel.grid = element_blank())+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+
  theme(strip.background =element_blank(),
        strip.text =element_text(size=20),
        legend.position = "none")

library(patchwork)

c=a+b

ggsave(c,filename ="biomass_layer.pdf",width = 11,height = 6)


##### 4. ABUNDANCE ANALYSES

library(MCMCglmm)
library(tidyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(hdrcde)
library(coda)

###########################   MCMC analyses

# See the distribution of abundances 

table(df_sub$Family,df_sub$Number_caught) # most of this is just one individual caught 1 time

abund=aggregate(Number_caught~line+site+occ+Family,function(x) sum(x,na.rm=T),drop=F,data=df_sub[df_sub$Class=="Insecta",])

abund$Number_caught[is.na(abund$Number_caught)]=0

abund$land_use=left_join(abund,hab[!duplicated(hab$site),],by="site")$land_use

abund$hab=left_join(abund,hab,by=c("site","line"))$Sand

abund$occ=factor(abund$occ)


table(abund$occ,abund$Number_caught,abund$Family)

# Take everything that was caught at least 6 times

ab_sub=abund[abund$Family%in%c("Scarabaeidae","Tenebrionidae","Termitidae"),]

table(ab_sub$Family,ab_sub$Number_caught)

ab.w=reshape(ab_sub, timevar="Family", idvar = c("line","site", "occ"),v.names = "Number_caught",direction="wide")

colnames(ab.w)[6:8]=rownames(table(ab_sub$Family,ab_sub$Number_caught))

### NULL MODEL
prior = list(R = list(V = diag(3)/4, n = 3, nu=0.002),
             G = list(G1 = list(V = diag(3)/4, n = 3, nu=0.002)))

m0.1=MCMCglmm(cbind( Scarabaeidae, Tenebrionidae, Termitidae)~trait , # + trait:Latitude,
            random = ~ us(trait):site,rcov = ~us(trait):units,prior = prior, family = rep("poisson", 3), nitt = 700000, burnin = 110000,
            pr=F,thin=250, data = ab.w)

m0.2=MCMCglmm(cbind(Scarabaeidae, Tenebrionidae, Termitidae)~trait , # + trait:Latitude,
            random = ~ us(trait):site,rcov = ~us(trait):units,prior = prior, family = rep("poisson", 3), nitt = 700000, burnin = 110000,
            pr=F,thin=250, data = ab.w)

m0.3=MCMCglmm(cbind(Scarabaeidae, Tenebrionidae, Termitidae)~trait , # + trait:Latitude,
            random = ~ us(trait):site,rcov = ~us(trait):units,prior = prior, family = rep("poisson", 3), nitt = 700000, burnin = 110000,
            pr=F,thin=250, data = ab.w)

param.coda.fixed0=mcmc.list(list(mcmc(m0.1$Sol),mcmc(m0.2$Sol),mcmc(m0.3$Sol)))

summary(param.coda.fixed0)
gelman.diag(param.coda.fixed0,multivariate=F)

par(mar=c(2,2,2,2))
plot(param.coda.fixed0,smooth=F) # The different colors indicate different chains


# Termitidae are pretty unstable

### Remove Termitidae and add habitat

prior = list(R = list(V = diag(2)/3, n = 2, nu=0.002),
             G = list(G1 = list(V = diag(2)/3, n = 4, nu=0.002)))

m1=MCMCglmm(cbind(Scarabaeidae, Tenebrionidae)~trait+trait:hab, # + trait:Latitude,
            random = ~ us(trait):site,rcov = ~us(trait):units,prior = prior, family = rep("poisson", 2), nitt = 600000, burnin = 100000,
            pr=F,thin=250, data = ab.w)

m2=MCMCglmm(cbind(Scarabaeidae, Tenebrionidae)~trait+trait:hab, # + trait:Latitude,
            random = ~ us(trait):site,rcov = ~us(trait):units,prior = prior, family = rep("poisson", 2), nitt = 600000, burnin = 100000,
            pr=F,thin=250, data = ab.w)

m3=MCMCglmm(cbind(Scarabaeidae, Tenebrionidae)~trait+trait:hab , # + trait:Latitude,
            random = ~ us(trait):site,rcov = ~us(trait):units,prior = prior, family = rep("poisson", 2), nitt = 600000, burnin = 100000,
            pr=F,thin=250, data = ab.w)

param.coda.fixed.Hab=mcmc.list(list(mcmc(m1$Sol),mcmc(m2$Sol),mcmc(m3$Sol)))

summary(param.coda.fixed.Hab) # 95 CI indicate a macrohabitat effect
gelman.diag(param.coda.fixed.Hab,multivariate=F)

par(mar=c(2,2,2,2))
plot(param.coda.fixed.Hab,smooth=F) # The different colors indicate different chains

### Add OCCASION

prior = list(R = list(V = diag(2)/3, n = 2, nu=0.002),
             G = list(G1 = list(V = diag(2)/3, n = 2, nu=0.002)))

m1.oc=MCMCglmm(cbind(Scarabaeidae, Tenebrionidae)~trait+trait:hab+trait:occ, # + trait:Latitude,
            random = ~ us(trait):site,rcov = ~us(trait):units,prior = prior, family = rep("poisson", 2), nitt = 600000, burnin = 100000,
            pr=F,thin=250, data = ab.w)

m2.oc=MCMCglmm(cbind(Scarabaeidae, Tenebrionidae)~trait+trait:hab+trait:occ, # + trait:Latitude,
            random = ~ us(trait):site,rcov = ~us(trait):units,prior = prior, family = rep("poisson", 2), nitt = 600000, burnin = 100000,
            pr=F,thin=250, data = ab.w)

m3.oc=MCMCglmm(cbind(Scarabaeidae, Tenebrionidae)~trait+trait:hab +trait:occ, # + trait:Latitude,
            random = ~ us(trait):site,rcov = ~us(trait):units,prior = prior, family = rep("poisson", 2), nitt = 600000, burnin = 100000,
            pr=F,thin=250, data = ab.w)

param.coda.fixed.oc=mcmc.list(list(mcmc(m1.oc$Sol),mcmc(m2.oc$Sol),mcmc(m3.oc$Sol)))

summary(param.coda.fixed.oc)
gelman.diag(param.coda.fixed.oc,multivariate=F)

par(mar=c(2,2,2,2))
plot(param.coda.fixed.oc,smooth=F) # The different colors indicate different chains

library(MCMCvis)

pdf("abund_param_fix_V2.pdf",width=6,height=7)
MCMCplot(param.coda.fixed.oc,ref_ovl = T,xlab="Fixed effects")
dev.off()

# Site-level covariance

m1.sub.site=m1.oc$VCV[,c(1:2,4)]
m2.sub.site=m1.oc$VCV[,c(1:2,4)]
m3.sub.site=m1.oc$VCV[,c(1:2,4)]

colnames(m1.sub.site)=colnames(m2.sub.site)=colnames(m3.sub.site)=c("S-S","S-T",
                                                                    "T-T")

param.coda.vcv=mcmc.list(list(mcmc(m1.sub.site),mcmc(m2.sub.site),mcmc(m3.sub.site)))

pdf("abund_param_site_V2.pdf",width=7,height=5)
MCMCplot(param.coda.vcv,ref_ovl = T,xlab="Random covariance")
dev.off()

# Site-level covariance

m1.sub.res=m1.oc$VCV[,c(5,6,8)]
m2.sub.res=m1.oc$VCV[,c(5,6,8)]
m3.sub.res=m1.oc$VCV[,c(5,6,8)]

colnames(m1.sub.res)=colnames(m2.sub.res)=colnames(m3.sub.res)=c("S-S","S-T",
                                                                 "T-T")

param.coda.vcv=mcmc.list(list(mcmc(m1.sub.res),mcmc(m2.sub.res),mcmc(m3.sub.res)))

pdf("abund_param_resid_V2.pdf",width=7,height=5)
MCMCplot(param.coda.vcv,ref_ovl = T,xlab="Residual covariance")
dev.off()

### ADD LINE

prior = list(R = list(V = diag(2)/3, n = 2, nu=0.002),
             G = list(G1 = list(V = diag(2)/3, n = 2, nu=0.002)))

m1.l=MCMCglmm(cbind(Scarabaeidae, Tenebrionidae)~trait+trait:hab+trait:occ+line, # + trait:Latitude,
               random = ~ us(trait):site,rcov = ~us(trait):units,prior = prior, family = rep("poisson", 2), nitt = 600000, burnin = 100000,
               pr=F,thin=250, data = ab.w)

m2.l=MCMCglmm(cbind(Scarabaeidae, Tenebrionidae)~trait+trait:hab+trait:occ+line, # + trait:Latitude,
               random = ~ us(trait):site,rcov = ~us(trait):units,prior = prior, family = rep("poisson", 2), nitt = 600000, burnin = 100000,
               pr=F,thin=250, data = ab.w)

m3.l=MCMCglmm(cbind(Scarabaeidae, Tenebrionidae)~trait+trait:hab +trait:occ+line, # + trait:Latitude,
               random = ~ us(trait):site,rcov = ~us(trait):units,prior = prior, family = rep("poisson", 2), nitt = 600000, burnin = 100000,
               pr=F,thin=250, data = ab.w)

param.coda.fixed.l=mcmc.list(list(mcmc(m1.l$Sol),mcmc(m2.l$Sol),mcmc(m3.l$Sol)))

summary(param.coda.fixed.l)
gelman.diag(param.coda.fixed.l,multivariate=F)

par(mar=c(2,2,2,2))
plot(param.coda.fixed.l,smooth=F) 

### No differences between lines 

########### ORTHER FAMILIES POOLED BUT NOT OMITTED

ab_sub=abund[!abund$Family%in%c("Unidentified"),]

table(ab_sub$Family,ab_sub$Number_caught)


ab.w=reshape(ab_sub, timevar="Family", idvar = c("line","site", "occ"),v.names = "Number_caught",direction="wide")

ab.w$Others=rowSums(ab.w[,c(6:15,17,19)])
colnames(ab.w)[6:19]=rownames(table(ab_sub$Family,ab_sub$Number_caught))

### NULL MODEL
prior = list(R = list(V = diag(3)/4, n = 3, nu=0.002),
             G = list(G1 = list(V = diag(3)/4, n = 3, nu=0.002)))

m0.1=MCMCglmm(cbind( Scarabaeidae, Tenebrionidae, Others)~trait , # + trait:Latitude,
              random = ~ us(trait):site,rcov = ~us(trait):units,prior = prior, family = rep("poisson", 3), nitt = 600000, burnin = 100000,
              pr=F,thin=250, data = ab.w)

m0.2=MCMCglmm(cbind(Scarabaeidae, Tenebrionidae, Others)~trait , # + trait:Latitude,
              random = ~ us(trait):site,rcov = ~us(trait):units,prior = prior, family = rep("poisson", 3), nitt = 600000, burnin = 100000,
              pr=F,thin=250, data = ab.w)

m0.3=MCMCglmm(cbind(Scarabaeidae, Tenebrionidae, Others)~trait , # + trait:Latitude,
              random = ~ us(trait):site,rcov = ~us(trait):units,prior = prior, family = rep("poisson", 3), nitt = 600000, burnin = 100000,
              pr=F,thin=250, data = ab.w)

param.coda.fixed0=mcmc.list(list(mcmc(m0.1$Sol),mcmc(m0.2$Sol),mcmc(m0.3$Sol)))

summary(param.coda.fixed0)
gelman.diag(param.coda.fixed0,multivariate=F)

par(mar=c(2,2,2,2))
plot(param.coda.fixed0,smooth=F) # The different colors indicate different chains

### ADD HABITAT
prior = list(R = list(V = diag(3)/4, n = 3, nu=0.002),
             G = list(G1 = list(V = diag(3)/4, n = 3, nu=0.002)))

m1.h=MCMCglmm(cbind( Scarabaeidae, Tenebrionidae, Others)~trait +trait:hab, # + trait:Latitude,
              random = ~ us(trait):site,rcov = ~us(trait):units,prior = prior, family = rep("poisson", 3), nitt = 600000, burnin = 100000,
              pr=F,thin=250, data = ab.w)

m2.h=MCMCglmm(cbind(Scarabaeidae, Tenebrionidae, Others)~trait +trait:hab, # + trait:Latitude,
              random = ~ us(trait):site,rcov = ~us(trait):units,prior = prior, family = rep("poisson", 3), nitt = 600000, burnin = 100000,
              pr=F,thin=250, data = ab.w)

m3.h=MCMCglmm(cbind(Scarabaeidae, Tenebrionidae, Others)~trait +trait:hab, # + trait:Latitude,
              random = ~ us(trait):site,rcov = ~us(trait):units,prior = prior, family = rep("poisson", 3), nitt = 600000, burnin = 100000,
              pr=F,thin=250, data = ab.w)

param.coda.fixedH=mcmc.list(list(mcmc(m1.h$Sol),mcmc(m2.h$Sol),mcmc(m3.h$Sol)))

summary(param.coda.fixedH)
gelman.diag(param.coda.fixedH,multivariate=F)

par(mar=c(2,2,2,2))
plot(param.coda.fixedH,smooth=F) # The different colors indicate different chains


### ADD OCCASION
prior = list(R = list(V = diag(3)/4, n = 3, nu=0.002),
             G = list(G1 = list(V = diag(3)/4, n = 3, nu=0.002)))

m1.occ=MCMCglmm(cbind( Scarabaeidae, Tenebrionidae, Others)~trait +trait:hab+trait:occ, # + trait:Latitude,
              random = ~ us(trait):site,rcov = ~us(trait):units,prior = prior, family = rep("poisson", 3), nitt = 600000, burnin = 100000,
              pr=F,thin=250, data = ab.w)

m2.occ=MCMCglmm(cbind(Scarabaeidae, Tenebrionidae, Others)~trait +trait:hab+trait:occ, # + trait:Latitude,
              random = ~ us(trait):site,rcov = ~us(trait):units,prior = prior, family = rep("poisson", 3), nitt = 600000, burnin = 100000,
              pr=F,thin=250, data = ab.w)

m3.occ=MCMCglmm(cbind(Scarabaeidae, Tenebrionidae, Others)~trait +trait:hab+trait:occ, # + trait:Latitude,
              random = ~ us(trait):site,rcov = ~us(trait):units,prior = prior, family = rep("poisson", 3), nitt = 600000, burnin = 100000,
              pr=F,thin=250, data = ab.w)

param.coda.fixedOcc=mcmc.list(list(mcmc(m1.occ$Sol),mcmc(m2.occ$Sol),mcmc(m3.occ$Sol)))

summary(param.coda.fixedOcc)
gelman.diag(param.coda.fixedOcc,multivariate=F)

par(mar=c(2,2,2,2))
plot(param.coda.fixedOcc,smooth=F) # The different colors indicate different chains

library(MCMCvis)

pdf("abund_param_fix.pdf",width=6,height=7)
MCMCplot(param.coda.fixedOcc,ref_ovl = T,xlab="Fixed effects")
dev.off()

# Site-level covariance

m1.sub.site=m1.occ$VCV[,c(1:3,5,6,9)]
m2.sub.site=m2.occ$VCV[,c(1:3,5,6,9)]
m3.sub.site=m3.occ$VCV[,c(1:3,5,6,9)]

colnames(m1.sub.site)=colnames(m2.sub.site)=colnames(m3.sub.site)=c("S-S","S-T","S-O",
                                                              "T-T","T-O","O-O")

param.coda.vcv=mcmc.list(list(mcmc(m1.sub.site),mcmc(m2.sub.site),mcmc(m3.sub.site)))

pdf("abund_param_site.pdf",width=7,height=5)
MCMCplot(param.coda.vcv,ref_ovl = T,xlab="Random covariance")
dev.off()

# Site-level covariance

m1.sub.res=m1.occ$VCV[,c(10:12,14,15,18)]
m2.sub.res=m2.occ$VCV[,c(10:12,14,15,18)]
m3.sub.res=m3.occ$VCV[,c(10:12,14,15,18)]

colnames(m1.sub.res)=colnames(m2.sub.res)=colnames(m3.sub.res)=c("S-S","S-T","S-O",
                                                                    "T-T","T-O","O-O")

param.coda.vcv=mcmc.list(list(mcmc(m1.sub.res),mcmc(m2.sub.res),mcmc(m3.sub.res)))

pdf("abund_param_resid.pdf",width=7,height=5)
MCMCplot(param.coda.vcv,ref_ovl = T,xlab="Residual covariance")
dev.off()

# Plot

ab_sub$Family[!ab_sub$Family%in%c("Scarabaeidae","Tenebrionidae")]="Others"

ab.mu=aggregate(Number_caught~line+site+occ+Family+hab,data=ab_sub,function(x) sum(x))

levels(ab.mu$occ)=c("Nov","Dec","Jan","Feb","Mar")
ab.mu$Family=factor(ab.mu$Family,levels=c("Scarabaeidae","Tenebrionidae","Others"))

ggplot(ab.mu,aes(x=occ,y=Number_caught,colour =Family))+
  geom_boxplot(outliers=F)+
  facet_grid(hab~.)+
  scale_colour_manual(values = c("darkred","orange","grey")) +
  geom_point(alpha = 0.4, position = position_jitterdodge(jitter.width = 0.8),
             size = 2, aes(group = Family))+
  xlab("Month of sampling")+
  ylab("Abundance (# caught)")+
  coord_cartesian(ylim=c(0,11))+
  theme_bw(base_size = 20)+
  theme(panel.grid = element_blank())+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+
  theme(strip.background =element_blank(),
        strip.text =element_text(size=20),
        legend.title = element_blank(),legend.position = "top")

ggsave(filename ="abund.pdf",width = 8,height = 7)


### Check out effect of layer

abundL=aggregate(Number_caught~layer+line+site+occ+Family,function(x) sum(x,na.rm=T),drop=F,data=df_sub[df_sub$Class=="Insecta",])

abundL=abundL[!(abundL$line=="L01"&abundL$layer=="L03"),] # this combination doesn't exist

abundL$Number_caught[is.na(abundL$Number_caught)]=0

abundL$land_use=left_join(abundL,hab[!duplicated(hab$site),],by="site")$land_use

abundL$hab=left_join(abundL,hab,by=c("site","line"))$Sand

abundL$occ=factor(abundL$occ)

ab_sub=abundL[!abundL$Family%in%c("Unidentified"),]

table(ab_sub$Family,ab_sub$Number_caught)

ab.w=reshape(ab_sub, timevar="Family", idvar = c("layer","line","site", "occ"),v.names = "Number_caught",direction="wide")

ab.w$Others=rowSums(ab.w[,c(7:16,18,20)])
colnames(ab.w)[7:20]=rownames(table(ab_sub$Family,ab_sub$Number_caught))

### ADD LAYER
prior = list(R = list(V = diag(3)/4, n = 3, nu=0.002),
             G = list(G1 = list(V = diag(3)/4, n = 3, nu=0.002)))

### Change between L01 and L05
m1.L=MCMCglmm(cbind( Scarabaeidae, Tenebrionidae, Others)~trait +trait:layer, # + trait:Latitude,
              random = ~ us(trait):site,rcov = ~us(trait):units,prior = prior, family = rep("poisson", 3), nitt = 600000, burnin = 100000,
              pr=F,thin=250, data = ab.w[ab.w$line=="L05",])

m2.L=MCMCglmm(cbind(Scarabaeidae, Tenebrionidae, Others)~trait +trait:layer, # + trait:Latitude,
              random = ~ us(trait):site,rcov = ~us(trait):units,prior = prior, family = rep("poisson", 3), nitt = 600000, burnin = 100000,
              pr=F,thin=250, data = ab.w[ab.w$line=="L05",])

m3.L=MCMCglmm(cbind(Scarabaeidae, Tenebrionidae, Others)~trait +trait:layer, # + trait:Latitude,
              random = ~ us(trait):site,rcov = ~us(trait):units,prior = prior, family = rep("poisson", 3), nitt = 600000, burnin = 100000,
              pr=F,thin=250, data = ab.w[ab.w$line=="L05",])

param.coda.fixedL=mcmc.list(list(mcmc(m1.L$Sol),mcmc(m2.L$Sol),mcmc(m3.L$Sol)))

summary(param.coda.fixedL)
gelman.diag(param.coda.fixedL,multivariate=F) # models do not converge

par(mar=c(2,2,2,2))
plot(param.coda.fixedL,smooth=F) # The different colors indicate different chains

### Summary: Layer does not affect abundance 

### Check LAYER effect removing 0s 

abundL$analyze="no"

# for each occasion, site, and line, if there was at least one cpature, analyze = "yes"

for(a in 1:length(unique(abundL$occ))){
  
  for(b in 1:length(unique(abundL$site))){
    
    for(c in 1:length(unique(abundL$line))){
      
      if(any(abundL$Number_caught[abundL$occ==unique(abundL$occ)[a]&abundL$site==unique(abundL$site)[b]&abundL$line==unique(abundL$line)[c]]>0)){
        
        
        abundL$analyze[abundL$occ==unique(abundL$occ)[a]&abundL$site==unique(abundL$site)[b]&abundL$line==unique(abundL$line)[c]] <- "yes"
        
      }
      
    }
  }
  
  
}

ab_sub=abundL[!abundL$Family%in%c("Unidentified"),]

ab.w=reshape(ab_sub, timevar="Family", idvar = c("layer","line","site", "occ","analyze"),v.names = "Number_caught",direction="wide")

ab.w$Others=rowSums(ab.w[,c(8:17,19,21)])
colnames(ab.w)[8:21]=rownames(table(ab_sub$Family,ab_sub$Number_caught))

prior = list(R = list(V = diag(3)/4, n = 3, nu=0.002),
             G = list(G1 = list(V = diag(3)/4, n = 3, nu=0.002)))

### Change between L01 and L05
m1.L=MCMCglmm(cbind( Scarabaeidae, Tenebrionidae, Others)~trait +trait:layer, # + trait:Latitude,
              random = ~ us(trait):site,rcov = ~us(trait):units,prior = prior, family = rep("poisson", 3), nitt = 600000, burnin = 100000,
              pr=F,thin=250, data = ab.w[ab.w$analyze=="yes"&ab.w$line=="L01",])

m2.L=MCMCglmm(cbind(Scarabaeidae, Tenebrionidae, Others)~trait +trait:layer, # + trait:Latitude,
              random = ~ us(trait):site,rcov = ~us(trait):units,prior = prior, family = rep("poisson", 3), nitt = 600000, burnin = 100000,
              pr=F,thin=250, data = ab.w[ab.w$analyze=="yes"&ab.w$line=="L01",])

m3.L=MCMCglmm(cbind(Scarabaeidae, Tenebrionidae, Others)~trait +trait:layer, # + trait:Latitude,
              random = ~ us(trait):site,rcov = ~us(trait):units,prior = prior, family = rep("poisson", 3), nitt = 600000, burnin = 100000,
              pr=F,thin=250, data = ab.w[ab.w$analyze=="yes"&ab.w$line=="L01",])

param.coda.fixedL=mcmc.list(list(mcmc(m1.L$Sol),mcmc(m2.L$Sol),mcmc(m3.L$Sol)))

summary(param.coda.fixedL)
gelman.diag(param.coda.fixedL,multivariate=F)

par(mar=c(2,2,2,2))
plot(param.coda.fixedL,smooth=F) # The different colors indicate different chains

### Summary: Layer does not affect abundance 

# Plot:

ab_sub$Family[!ab_sub$Family%in%c("Scarabaeidae", "Tenebrionidae")]="Others"

unique(ab_sub$Family)

# double stratified: 

a=ggplot(ab_sub[ab_sub$line=="L01"&ab_sub$analyze=="yes",],aes(x=layer,y=Number_caught))+
  geom_boxplot(aes(colour=layer),outliers=F)+
  facet_grid(Family~.)+
  scale_colour_viridis(discrete = TRUE, alpha=0.8) +
  geom_jitter(color="grey", size=2, alpha=0.9) +
  scale_y_continuous(limits = c(0, 11), oob = scales::squish)+
  xlab("Sampling layer")+
  ylab("Abundance (# caught)")+
  ggtitle("(a) Double-stratified design")+
  theme_bw(base_size = 20)+
  theme(panel.grid = element_blank())+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+
  theme(strip.background =element_blank(),
        strip.text =element_text(size=20),
        legend.position = "none")

b=ggplot(ab_sub[ab_sub$line=="L05"&ab_sub$analyze=="yes",],aes(x=layer,y=Number_caught))+
  geom_boxplot(aes(colour=layer),outliers=F)+
  facet_grid(Family~.)+
  scale_colour_viridis(discrete = TRUE, alpha=0.8) +
  geom_jitter(color="grey", size=2, alpha=0.9) +
  scale_y_continuous(limits = c(0, 11), oob = scales::squish)+
  xlab("Sampling layer")+
  ylab("Abundance (# caught)")+
  ggtitle("(b) Three-stratified design")+
  theme_bw(base_size = 20)+
  theme(panel.grid = element_blank())+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+
  theme(strip.background =element_blank(),
        strip.text =element_text(size=20),
        legend.position = "none")

library(patchwork)

c=a+b

ggsave(c,filename ="abundance_layer.pdf",width = 11,height = 10)

