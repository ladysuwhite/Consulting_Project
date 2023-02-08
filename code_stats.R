height <- read.csv("height_2021.csv")
wp_2021 <- read.csv("wp_2021.csv")




### imputation code
for (j in 183:192){
  d=100000
  d1=10000000
  index=0
  print(j)
  for(i in 1:182){
    if( height$Treatment[j]==height$Treatment[i]){
      if( height$Water[j]==height$Water[i]){
        d1=sqrt(abs(height$height.pre[j]^2-height$height.pre[i]^2))
      }}
    if(d1<d)
    {d=d1
    index=i}
  }
  print(height$PlantID[index])
}

##################
### imputation results are shred in all three cases
##################

missing_id=c(104,115,128,302,326,407,417,420,523)
imputation_id=c(123,501,301,116,110,423,425,214,301)


for (j in 1:9){
  miss_index=which(missing_id[j]==wp_2021$PlantID)
  inpute_index=which(imputation_id[j]==wp_2021$PlantID)
  wp_2021$WaterPot[miss_index]=wp_2021$WaterPot[inpute_index]
}

###############
## 1 water potential
################

###
# make ANCOVA data for wp
###
preheight=rep(1,192)
for(i in 1:192){
  j=1
  while(height$PlantID[j]!=wp_2021$PlantID[i]){
    j=j+1
  }
  preheight[i]=height$height.pre[j]
}

wp_2021_red1=data.frame(wp_2021,preheight)

### 
# ANOVA impotation data
###

wp=aov(WaterPot ~ Genotype+Water+Nitrogen+Genotype*Water+
         Water*Nitrogen+Genotype*Nitrogen+Genotype*Water*Nitrogen, 
      data = wp_2021)
summary(wp)
###
# ANOVA original data
###

wp_2021_red=na.omit(wp_2021)
wp_red=aov(WaterPot ~ Genotype+Water+Nitrogen+Genotype*Water+
         Water*Nitrogen+Genotype*Nitrogen+Genotype*Water*Nitrogen, 
       data = wp_2021_red)
summary(wp_red)


wp_2021_red1=na.omit(data.frame(wp_2021,preheight))
water_potential=wp_2021_red1$WaterPot
wp_2021_red1=data.frame(wp_2021_red1,water_potential)

###
# ANCOVA original data
###
ggscatter(
  wp_2021_red1, x = "preheight", y = "water_potential",
  facet.by  = c("Water","Genotype"), 
  short.panel.labs = FALSE
)+
stat_smooth(method = "loess", span = 0.8)

wp_2021_red1_out=aov( water_potential ~ Genotype+Water+Nitrogen+Genotype*Water+
                               Water*Nitrogen+Genotype*Nitrogen+Genotype*Water*Nitrogen+preheight, 
                   data = wp_2021_red1) 
summary(wp_2021_red1_out)


###
# plot tucky HSD
###
tukey_wp<-TukeyHSD(wp)
tukey_wp<-TukeyHSD(wp_red)

rownames(tukey_wp$Genotype)
library(ggplot2)
plot_TukeyHSD<-function(aa=aa){
ggplot(data=aa,
       aes(y = Variables,x =difference, xmin = LowerLimit, xmax = UpperLimit))+
  geom_point() + 
  geom_errorbarh(height=.1) + 
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  xlab('TukeyHSD')+ ylab("Groups")+
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.x=element_text(size=8,face="bold"),
        axis.title.y=element_text(size=8,face="bold"),
        legend.key.size = unit(0.1, 'cm'),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold")
  )
}
aa=data.frame(difference=tukey_wp$Genotype[,1],Variables=rownames(tukey_wp$Genotype),
              LowerLimit=tukey_wp$Genotype[,2], UpperLimit=tukey_wp$Genotype[,3])
plot_TukeyHSD(aa=aa)

aa=data.frame(difference=tukey_wp$Water[,1],Variables=rownames(tukey_wp$Water),
              LowerLimit=tukey_wp$Water[,2], UpperLimit=tukey_wp$Water[,3])
plot_TukeyHSD(aa=aa)

aa=data.frame(difference=tukey_wp$Nitrogen[,1],Variables=rownames(tukey_wp$Nitrogen),
              LowerLimit=tukey_wp$Nitrogen[,2], UpperLimit=tukey_wp$Nitrogen[,3])
plot_TukeyHSD(aa=aa)

aa=data.frame(difference=tukey_wp$`Genotype:Water`[,1],Variables=rownames(tukey_wp$`Genotype:Water`),
              LowerLimit=tukey_wp$`Genotype:Water`[,2], UpperLimit=tukey_wp$`Genotype:Water`[,3])
plot_TukeyHSD(aa=aa)



###############
## 2 herbivore abundance
################
b2_2021 <- read.csv("b2_2021.csv")
b3_2021 <- read.csv("b3_2021.csv")
b4_2021 <- read.csv("b4_2021.csv")

###function to summary bet data
data_bet<-function(newdata,olddata){
    out=rep(9999,192)
    for(i in 1:192){
        for(j in 1:192){
            if(newdata$PlantID[j]==olddata$PlantID[i]){
                
                out[i]=sum(newdata[j,6:9])
            }}}
    return(out)
}

b2=data_bet(b2_2021,wp_2021)
b3=data_bet(b3_2021,wp_2021)
b4=data_bet(b4_2021,wp_2021)

setwd("")

PlantID=rbind(as.matrix(wp_2021$PlantID,ncol=1),as.matrix(wp_2021$PlantID,ncol=1),
as.matrix(wp_2021$PlantID,ncol=1))
Nitrogen=rbind(as.matrix(wp_2021$Nitrogen,ncol=1),as.matrix(wp_2021$Nitrogen,ncol=1),
as.matrix(wp_2021$Nitrogen,ncol=1))
Genotype=rbind(as.matrix(wp_2021$Genotype,ncol=1),as.matrix(wp_2021$Genotype,ncol=1),
as.matrix(wp_2021$Genotype,ncol=1))
Water=rbind(as.matrix(wp_2021$Water,ncol=1),as.matrix(wp_2021$Water,ncol=1),
as.matrix(wp_2021$Water,ncol=1))
Bet=rbind(as.matrix(b2,ncol=1),as.matrix(b3,ncol=1),
as.matrix(b4,ncol=1))
Round=rbind(as.matrix(rep("b2",192),ncol=1),as.matrix(rep("b3",192),ncol=1),
as.matrix(rep("b4",192),ncol=1))

Bet_data=data.frame(PlantID,Nitrogen,Genotype,Water,
Bet,Round)
Bet_data_red=na.omit(Bet_data)
save(Bet_data_red,file = "Bet_data_red.Rdata")

# full data
Bet_data_full=Bet_data
for (j in 1:10){
    miss_index=which(missing_id[j]==Bet_data$PlantID)
    inpute_index=which(imputation_id[j]==Bet_data$PlantID)
    Bet_data_full$Bet[miss_index]=Bet_data$Bet[inpute_index]
}
save(Bet_data_full,file = "Bet_data_full.Rdata")

library(lme4)
library(MASS)



### data without imputation, with random factor
mod11 = glmer.nb(Bet ~ Genotype*Water*Nitrogen+ (1|Round),
data=Bet_data_red)
car::Anova(mod11, type='II')

Bet_data_red$WN <- paste0(Bet_data_red$Water, Bet_data_red$Nitrogen)
Bet_data_red$GWN <- paste0(Bet_data_red$Genotype,Bet_data_red$Water, Bet_data_red$Nitrogen)
mod11 = glmer.nb(Bet ~ Genotype+Water+Nitrogen+ factor(WN) + (1|Round),
data=Bet_data_red)
summary(mod11)
mod11.geno <- confint(glht(mod11, mcp(Genotype="Tukey")))  ##Confidence interval
mod11.ni <- confint(glht(mod11, mcp(Nitrogen="Tukey")))
#mod11.geni <- confint(glht(mod11,linfct= mcp(WaterNitrogen="Tukey")))
mod11.geni <- confint(lsmeans(mod11, pairwise ~ WN)[[2]])

aa=data.frame(difference=mod11.geno$confint[,1],Variables=rownames(tukey_wp$Genotype),
LowerLimit=mod11.geno$confint[,2], UpperLimit=mod11.geno$confint[,3])
plot_TukeyHSD(aa=aa)

aa=data.frame(difference=mod11.ni$confint[,1],Variables=rownames(tukey_wp$Nitrogen),
LowerLimit=mod11.ni$confint[,2], UpperLimit=mod11.ni$confint[,3])
plot_TukeyHSD(aa=aa)

aa=data.frame(difference=-mod11.geni[,2],Variables=rownames(tukey_wp$`Water:Nitrogen`),
LowerLimit=-mod11.geni[,6], UpperLimit=-mod11.geni[,5])
plot_TukeyHSD(aa=aa)

### data with imputation, with random factor
mod11_f = glmer.nb(Bet ~ Genotype*Water*Nitrogen+ (1|Round),
data=Bet_data_full)
car::Anova(mod11_f, type='II')

Bet_data_full$WN <- paste0(Bet_data_full$Water, Bet_data_full$Nitrogen)
Bet_data_full$GWN <- paste0(Bet_data_full$Genotype,Bet_data_full$Water, Bet_data_full$Nitrogen)
mod11_f = glmer.nb(Bet ~ Genotype+Water+Nitrogen+ WN + GWN+ (1|Round),
data=Bet_data_full)
summary(mod11_f)
mod11.geno.f <- confint(glht(mod11_f, mcp(Genotype="Tukey")))
mod11.ni.f <- confint(glht(mod11_f, mcp(Nitrogen="Tukey")))
mod11.geni.f <- confint(lsmeans(mod11_f, pairwise ~ WN)[[2]])
mod11.geniwa.f <- confint(lsmeans(mod11_f, pairwise ~ GWN)[[2]])

aa=data.frame(difference=mod11.geno.f$confint[,1],Variables=rownames(tukey_wp$Genotype),
LowerLimit=mod11.geno.f$confint[,2], UpperLimit=mod11.geno.f$confint[,3])
plot_TukeyHSD(aa=aa)

aa=data.frame(difference=mod11.ni.f$confint[,1],Variables=rownames(tukey_wp$Nitrogen),
LowerLimit=mod11.ni.f$confint[,2], UpperLimit=mod11.ni.f$confint[,3])
plot_TukeyHSD(aa=aa)

aa=data.frame(difference=-mod11.geni.f[,2],Variables=rownames(tukey_wp$`Water:Nitrogen`),
LowerLimit=-mod11.geni.f[,6], UpperLimit=-mod11.geni.f[,5])
plot_TukeyHSD(aa=aa)


### data without imputation, without random factor
mod12 = glm.nb(Bet ~ Genotype*Water*Nitrogen,  data=Bet_data_red)
car::Anova(mod12, type='II')

mod12 = glm.nb(Bet ~ Genotype+Water+Nitrogen+WN, data=Bet_data_red)
summary(mod12)
mod12.geno <-  confint(lsmeans(mod12, pairwise ~ Genotype)[[2]])
mod12.ni <- confint(lsmeans(mod12, pairwise ~ Nitrogen)[[2]])
mod12.geni <- confint(lsmeans(mod12, pairwise ~ WN)[[2]])

aa=data.frame(difference=-mod12.geno[,2],Variables=rownames(tukey_wp$Genotype),
LowerLimit=-mod12.geno[,6], UpperLimit=-mod12.geno[,5])
plot_TukeyHSD(aa=aa)

aa=data.frame(difference=-mod12.ni[,2],Variables=rownames(tukey_wp$Nitrogen),
LowerLimit=-mod12.ni[,6], UpperLimit=-mod12.ni[,5])
plot_TukeyHSD(aa=aa)

aa=data.frame(difference=-mod12.geni[,2],Variables=rownames(tukey_wp$`Water:Nitrogen`),
LowerLimit=-mod12.geni[,6], UpperLimit=-mod12.geni[,5])
plot_TukeyHSD(aa=aa)


### data with imputation, without random factor
mod12_f = glm.nb(Bet ~ Genotype*Water*Nitrogen, data=Bet_data_full)
car::Anova(mod12_f, type='II')

mod12_f = glm.nb(Bet ~ Genotype+Water+Nitrogen+ WN, data=Bet_data_full)
summary(mod12_f)
summary(mod12)
mod12.geno.f <-  confint(lsmeans(mod12_f, pairwise ~ Genotype)[[2]])
mod12.ni.f <- confint(lsmeans(mod12_f, pairwise ~ Nitrogen)[[2]])
mod12.geni.f <- confint(lsmeans(mod12_f, pairwise ~ WN)[[2]])

aa=data.frame(difference=-mod12.geno.f[,2],Variables=rownames(tukey_wp$Genotype),
LowerLimit=-mod12.geno.f[,6], UpperLimit=-mod12.geno.f[,5])
plot_TukeyHSD(aa=aa)

aa=data.frame(difference=-mod12.ni.f[,2],Variables=rownames(tukey_wp$Nitrogen),
LowerLimit=-mod12.ni.f[,6], UpperLimit=-mod12.ni.f[,5])
plot_TukeyHSD(aa=aa)

aa=data.frame(difference=-mod12.geni.f[,2],Variables=rownames(tukey_wp$`Water:Nitrogen`),
LowerLimit=-mod12.geni.f[,6], UpperLimit=-mod12.geni.f[,5])
plot_TukeyHSD(aa=aa)

#############################
### 3 potato leafhopper damage
##############################

plh_data <- read.csv("plh_2021.csv")
##### split the data into two rounds
plh1 <- plh_data[1:192,]
plh2 <- plh_data[193:384,]

pre_height <- c()
for (i in 1:192){
    pre_height[i] <- height$height.pre[which(height$PlantID == plh1$PlantID[i])]
}

##### add preheight column
plh1["height"] <- pre_height
plh2["height"] <- pre_height

plh <- rbind(plh1, plh2)
##### original data
plh_red <- na.omit(plh)

plh1_full <- plh1
for (j in 1:10){
    miss_index=which(missing_id[j]==plh1$PlantID)
    inpute_index=which(imputation_id[j]==plh1$PlantID)
    plh1_full$DamagedShoots[miss_index]=plh1$DamagedShoots[inpute_index]
}

plh2_full <- plh2
for (j in 1:10){
    miss_index=which(missing_id[j]==plh2$PlantID)
    inpute_index=which(imputation_id[j]==plh2$PlantID)
    plh2_full$DamagedShoots[miss_index]=plh2$DamagedShoots[inpute_index]
}

##### imputation data
plh_full <- rbind(plh1_full, plh2_full)

######### ANOVA for original data
model1_original <- lmer(DamagedShoots/10 ~ Genotype*Water*Nitrogen+ (1|Observer), data = plh_red)
anova(model1_original)

######### ANOVA for imputed data
model1_imputed <- lmer(DamagedShoots/10 ~ Genotype*Water*Nitrogen+ (1|Observer), data = plh_full)
anova(model1_imputed)

##### plot Tukey HSD
library(multcomp)
tukey_N <- summary(glht(model1_imputed, linfct=mcp(Nitrogen = "Tukey")))
tukey_N_interval <- confint(tukey_N)
aa=data.frame(difference=tukey_N_interval$confint[,1],Variables=rownames(tukey_N_interval$confint),
LowerLimit=tukey_N_interval$confint[,2], UpperLimit=tukey_N_interval$confint[,3])
plot_TukeyHSD(aa=aa)

tukey_W <- summary(glht(model1_imputed, linfct=mcp(Water = "Tukey")))
tukey_W_interval <- confint(tukey_W)
aa=data.frame(difference=tukey_W_interval$confint[,1],Variables=rownames(tukey_W_interval$confint),
LowerLimit=tukey_W_interval$confint[,2], UpperLimit=tukey_W_interval$confint[,3])
plot_TukeyHSD(aa=aa)

tukey_G <- summary(glht(model1_imputed, linfct=mcp(Genotype = "Tukey")))
tukey_G_interval <- confint(tukey_G)
aa=data.frame(difference=tukey_G_interval$confint[,1],Variables=rownames(tukey_G_interval$confint),
LowerLimit=tukey_G_interval$confint[,2], UpperLimit=tukey_G_interval$confint[,3])
plot_TukeyHSD(aa=aa)

plh_full$GN <- interaction(plh_full$Genotype, plh_full$Nitrogen)
mod <- lmer(DamagedShoots/10 ~ GN + Water + Genotype:Water + Water:Nitrogen + Genotype:Water:Nitrogen
+ (1|Observer), data = plh_full )
tukey_GN <- summary(glht(mod, linfct=mcp(GN="Tukey")))
tukey_GN_interval <- confint(tukey_GN)
aa=data.frame(difference=tukey_GN_interval$confint[,1],Variables=rownames(tukey_GN_interval$confint),
LowerLimit=tukey_GN_interval$confint[,2], UpperLimit=tukey_GN_interval$confint[,3])
plot_TukeyHSD(aa=aa)

plh_full$GW <- interaction(plh_full$Genotype, plh_full$Water)
mod <- lmer(DamagedShoots/10 ~ GW + Nitrogen + Genotype:Nitrogen + Water:Nitrogen + Genotype:Water:Nitrogen
+ (1|Observer), data = plh_full )
tukey_GW <- summary(glht(mod, linfct=mcp(GW="Tukey")))
tukey_GW_interval <- confint(tukey_GW)
aa=data.frame(difference=tukey_GW_interval$confint[,1],Variables=rownames(tukey_GW_interval$confint),
LowerLimit=tukey_GW_interval$confint[,2], UpperLimit=tukey_GW_interval$confint[,3])
plot_TukeyHSD(aa=aa)


######  ANCOVA for original
plh_red_scale <- plh_red
plh_red_scale$DamagedShoots <- plh_red$DamagedShoots / 10
###### check assumptions
library(ggpubr)
ggscatter(
plh_red_scale, x = "height", y = "DamagedShoots",
facet.by  = c("Nitrogen","Water"),
short.panel.labs = FALSE
)+
stat_smooth(method = "loess", span = 0.9)

ggscatter(
plh_red_scale, x = "height", y = "DamagedShoots",
facet.by  = c("Nitrogen","Genotype"),
short.panel.labs = FALSE
)+
stat_smooth(method = "loess", span = 0.9)

ggscatter(
plh_red_scale, x = "height", y = "DamagedShoots",
facet.by  = c("Water","Genotype"),
short.panel.labs = FALSE
)+
stat_smooth(method = "loess", span = 0.9)

model2=lmer(DamagedShoots/10 ~ height + Genotype*Water*Nitrogen + (1|Observer), data = plh_red)
anova(model2)

##### plot Tukey HSD
tukey_N <- summary(glht(model2, linfct=mcp(Nitrogen = "Tukey")))
tukey_N_interval <- confint(tukey_N)
aa=data.frame(difference=tukey_N_interval$confint[,1],Variables=rownames(tukey_N_interval$confint),
LowerLimit=tukey_N_interval$confint[,2], UpperLimit=tukey_N_interval$confint[,3])
plot_TukeyHSD(aa=aa)

tukey_W <- summary(glht(model2, linfct=mcp(Water = "Tukey")))
tukey_W_interval <- confint(tukey_W)
aa=data.frame(difference=tukey_W_interval$confint[,1],Variables=rownames(tukey_W_interval$confint),
LowerLimit=tukey_W_interval$confint[,2], UpperLimit=tukey_W_interval$confint[,3])
plot_TukeyHSD(aa=aa)

tukey_G <- summary(glht(model2, linfct=mcp(Genotype = "Tukey")))
tukey_G_interval <- confint(tukey_G)
aa=data.frame(difference=tukey_G_interval$confint[,1],Variables=rownames(tukey_G_interval$confint),
LowerLimit=tukey_G_interval$confint[,2], UpperLimit=tukey_G_interval$confint[,3])
plot_TukeyHSD(aa=aa)





