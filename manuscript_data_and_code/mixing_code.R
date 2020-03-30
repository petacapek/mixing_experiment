################################################################################################
#####Stoichiometric and microbial controls over nutrients mineralization and immobilization#####
################################################################################################
################################################################################################
#Libraries
library(dplyr)
library(reshape)
library(ggplot2)
library(gridExtra)
library(vegan)
library(deSolve)
library(FME)
library(DEoptim)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Defining ggplot theme
theme_min<-theme(axis.text.x=element_text(vjust=0.2, size=18, colour="black"),
                 axis.text.y=element_text(hjust=0.2, size=18, colour="black"),
                 axis.title=element_text(size=18, colour="black"),
                 axis.line=element_line(size=0.5, colour="black"),
                 strip.text=element_text(size=18, face="bold"),
                 axis.ticks=element_line(size=1, colour="black"),
                 axis.ticks.length=unit(-0.05, "cm"),
                 panel.background=element_rect(colour="black", fill="white"),
                 panel.grid=element_line(linetype=0),
                 legend.text=element_text(size=14, colour="black"),
                 legend.title=element_text(size=14, colour="black"),
                 legend.position=c("right"),
                 legend.key.size=unit(1, "cm"),
                 strip.background=element_rect(fill="grey98", colour="black"),
                 legend.key=element_rect(fill="white", size=1.2),
                 legend.spacing=unit(0.5, "cm"),
                 plot.title=element_text(size=18, face="bold", hjust=-0.05))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Reading data
##Experiment data
mix<-read.csv("./manuscript_data_and_code/mixing_data.csv")
summary(mix)


#~~~~~~~~~~~~~~~~~~~~~~~Not included in manuscript~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ##Communities
# ###Bacteria
# bac<-t(read.csv("./manuscript_data_and_code/mixing_bacteria.csv", header = F)[-1, -c(1, 42:49)])
# rownames(bac)<-read.csv("./manuscript_data_and_code/mixing_bacteria.csv", header = F)[1,-c(1, 42:49)]
# #Fungi
# fungi<-t(read.csv("./manuscript_data_and_code/mixing_fungi.csv", header = F)[-1, -c(1, 42:49)])
# rownames(fungi)<-read.csv("./manuscript_data_and_code/mixing_fungi.csv", header = F)[1,-c(1, 42:49)]
# 
# ###Labels
# lb<-read.csv("./manuscript_data_and_code/community_labels.csv")
# lb$SampleID<-as.factor(lb$SampleID)
# lb$SampleID<-factor(lb$SampleID, levels = rownames(bac))
# lb<-lb[order(lb$SampleID), ]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
###############################################################################################
######################################Statistics###############################################
###############################################################################################
#Initial conditions
##POOLS
###Microbial biomass carbon
Cmic_i<-glm(Cmic~Horizon+Plesne/Horizon, data = mix[(mix$Labelling=="NO" &
                                                         mix$TIME==0), ], family = Gamma)
summary(Cmic_i)
anova(Cmic_i, test="F")

###Microbial biomass nitrogen
Nmic_i<-glm(Nmic~Horizon+Plesne/Horizon, data = mix[(mix$Labelling=="NO" &
                                                       mix$TIME==0), ], family = Gamma)
summary(Nmic_i)
anova(Nmic_i, test="F")

###Microbial biomass phosphorus
Pmic_i<-glm(Pmic~Horizon+Plesne/Horizon, data = mix[(mix$Labelling=="NO" &
                                                       mix$TIME==0), ], family = Gamma)
summary(Pmic_i)
anova(Pmic_i, test="F")

###Water extractable organic carbon
DOC_i<-glm(DOC~Horizon+Plesne/Horizon, data = mix[(mix$Labelling=="NO" &
                                                       mix$TIME==0), ], family = Gamma)
summary(DOC_i)
anova(DOC_i, test="F")

###Water extractable organic nitrogen
DON_i<-glm(DON~Horizon+Plesne/Horizon, data = mix[(mix$Labelling=="NO" &
                                                     mix$TIME==0), ], family = Gamma)
summary(DON_i)
anova(DON_i, test="F")

###Water extractable organic phosphorus
DOP_i<-glm(DOP~Horizon+Plesne/Horizon, data = mix[(mix$Labelling=="NO" &
                                                     mix$TIME==0), ], family = Gamma)
summary(DOP_i)
anova(DOP_i, test="F")

###Water extractable amonia
NH_i<-glm(NH4~Horizon+Plesne/Horizon, data = mix[(mix$Labelling=="NO" &
                                                     mix$TIME==0), ], family = Gamma)
summary(NH_i)
anova(NH_i, test="F")

###Water extractable nitrates
NO_i<-glm(NO3~Horizon+Plesne/Horizon, data = mix[(mix$Labelling=="NO" &
                                                    mix$TIME==0), ], family = Gamma)
summary(NO_i)
anova(NO_i, test="F")

###Water extractable phosphates
PO_i<-glm(PO4~Horizon+Plesne/Horizon, data = mix[(mix$Labelling=="NO" &
                                                    mix$TIME==0), ], family = Gamma)
summary(PO_i)
anova(PO_i, test="F")

##STOICHIOMETRY
###Microbial biomass carbon to nitrogen ratio
CNmic_i<-glm(Cmic/Nmic~Horizon+Plesne/Horizon, data = mix[(mix$Labelling=="NO" &
                                                    mix$TIME==0), ], family = Gamma)
summary(CNmic_i)
anova(CNmic_i, test="F")

###Microbial biomass carbon to phosphorus ratio
CPmic_i<-glm(Cmic/Pmic~Horizon+Plesne/Horizon, data = mix[(mix$Labelling=="NO" &
                                                             mix$TIME==0), ], family = Gamma)
summary(CPmic_i)
anova(CPmic_i, test="F")

###Water extractable organic carbon to nitrogen ratio
DOCN_i<-glm(DOC/DON~Horizon+Plesne/Horizon, data = mix[(mix$Labelling=="NO" &
                                                             mix$TIME==0), ], family = Gamma)
summary(DOCN_i)
anova(DOCN_i, test="F")

###Water extractable organic carbon to phosphorus ratio
DOCP_i<-glm(DOC/DOP~Horizon+Plesne/Horizon, data = mix[(mix$Labelling=="NO" &
                                                          mix$TIME==0), ], family = Gamma)
summary(DOCP_i)
anova(DOCP_i, test="F")

#~~~~~~~~~~~~~~~~~~~~~~~Not included in manuscript~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ##MICROBIAL COMMUNITY
# ###Bacteria
# bnorm<-decostand(bac, method=c("normalize"))
# bdist<-vegdist(bnorm, method = "jaccard")
# 
# bad<-adonis(bdist~horizon+PL/horizon, lb)
# bad
# 
# ###Fungi
# fnorm<-decostand(fungi, method=c("normalize"))
# fdist<-vegdist(fnorm, method = "jaccard")
# 
# fad<-adonis(fdist~horizon+PL/horizon, lb)
# fad
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##TEMPORAL CHANGES
###Microbial biomass carbon
Cmic_t<-glm(Cmic~Horizon+Plesne/Horizon+TIME/Plesne/Horizon, data = mix[(mix$Labelling=="NO"), ], family = Gamma)

summary(Cmic_t)
anova(Cmic_t, test="F")

###Microbial biomass nitrogen
Nmic_t<-glm(Nmic~Horizon+Plesne/Horizon+TIME/Plesne/Horizon, data = mix[(mix$Labelling=="NO"), ], family = Gamma)

summary(Nmic_t)
anova(Nmic_t, test="F")

###Microbial biomass phosphorus
Pmic_t<-glm(Pmic~Horizon+Plesne/Horizon+TIME/Plesne/Horizon, data = mix[(mix$Labelling=="NO"), ], family = Gamma)

summary(Pmic_t)
anova(Pmic_t, test="F")

###Water extractable mineral nitrogen
Nm_t<-glm((NH4+NO3)~Horizon+Plesne/Horizon+TIME/Plesne/Horizon, data = mix[(mix$Labelling=="NO"), ], family = Gamma)

summary(Nm_t)
anova(Nm_t, test="F")

###Water extractable mineral phosphorus
Pm_t<-glm(PO4~Horizon+Plesne/Horizon+TIME/Plesne/Horizon, data = mix[(mix$Labelling=="NO"), ], family = Gamma)

summary(Pm_t)
anova(Pm_t, test="F")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
###########################################################################################
###############################Figures and calculations####################################
###########################################################################################
#MAIN TEXT
##Figure 1: Initial values of microbial biomass carbon to (A) nitrogen 
##(MBC:MBN) and (B) phosphorus (MBC:MBP) ratios, and water extractable 
##organic carbon to (C) nitrogen (DOC:DON) and (D) phosphorus (DOC:DOP) 
##ratios in two spruce forest soils (PL - Plešné and CT - Čertovo catchments) 
##mixed at five different ratios (i.e., 0:1, 0.25:0.75, 0.5:0.5, 0.75:0.25, 
##and 1:0 in respect to PL). Grey circles denote the litter horizon and 
##empty circles the topsoil organic horizon. Symbols show mean values and 
##error bars standard error of the mean (n = 4). Notice that plots have 
##different y-axes scales.

grid.arrange(
  mix %>% filter(Labelling=="NO" & TIME==0) %>% 
    group_by(Plesne, Horizon) %>% summarize(y.sd=sd(Cmic/Nmic)/sqrt(4), y=mean(Cmic/Nmic)) %>%
    ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Horizon), 
                                              position = position_dodge(width = 0.6))+
    geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd, colour=Horizon), width=0.1, lwd=0.5, 
                  position = position_dodge(width = 0.6))+
    theme_min+theme(legend.position = c(0.6, 0.7), legend.title = element_blank())+
    ylab(expression(paste("MBC:MBN (mol:mol)" )))+
    xlab("PL : CT mixing ratio")+scale_color_manual(values = c("black", "black"))+
    scale_fill_manual(values = c("grey", "white"))+
    ggtitle("A)")+scale_y_continuous(limits = c(10, 20),
                                     breaks = c(10, 12, 14, 16, 18, 20)),
  
  mix %>% filter(Labelling=="NO" & TIME==0) %>% 
    group_by(Plesne, Horizon) %>% summarize(y.sd=sd(Cmic/Pmic)/sqrt(4), y=mean(Cmic/Pmic)) %>%
    ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Horizon), 
                                              position = position_dodge(width = 0.6), show.legend = F)+
    geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd, colour=Horizon), width=0.1, lwd=0.5, 
                  position = position_dodge(width = 0.6), show.legend = F)+
    theme_min+theme(legend.position = c(0.8, 0.9), legend.title = element_blank(),
                    plot.margin = unit(c(0.05, 0.05, 0.05, 0.35), "in"))+
    ylab(expression(paste("MBC:MBP (mol:mol)" )))+
    xlab("PL : CT mixing ratio")+scale_color_manual(values = c("black", "black"))+
    scale_fill_manual(values = c("grey", "white"))+
    ggtitle("B)")+scale_y_continuous(limits = c(0, 80)),
  
  mix %>% filter(Labelling=="NO" & TIME==0) %>% 
    group_by(Plesne, Horizon) %>% summarize(y.sd=sd(DOC/DON)/sqrt(4), y=mean(DOC/DON)) %>%
    ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Horizon), 
                                              position = position_dodge(width = 0.6), show.legend = F)+
    geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd, colour=Horizon), width=0.1, lwd=0.5, 
                  position = position_dodge(width = 0.6), show.legend = F)+
    theme_min+theme(legend.position = c(0.8, 0.9), legend.title = element_blank())+
    ylab(expression(paste("DOC:DON (mol:mol)" )))+
    xlab("PL : CT mixing ratio")+scale_color_manual(values = c("black", "black"))+
    scale_fill_manual(values = c("grey", "white"))+
    ggtitle("C)")+scale_y_continuous(limits = c(10, 30)),
  
  mix %>% filter(Labelling=="NO" & TIME==0) %>% 
    group_by(Plesne, Horizon) %>% summarize(y.sd=sd(DOC/DOP)/sqrt(4), y=mean(DOC/DOP)) %>%
    ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Horizon), 
                                              position = position_dodge(width = 0.6), show.legend = F)+
    geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd, colour=Horizon), width=0.1, lwd=0.5, 
                  position = position_dodge(width = 0.6), show.legend = F)+
    theme_min+theme(legend.position = c(0.8, 0.9), legend.title = element_blank())+
    ylab(expression(paste("DOC:DOP (mol:mol)" )))+
    xlab("PL : CT mixing ratio")+scale_color_manual(values = c("black", "black"))+
    scale_fill_manual(values = c("grey", "white"))+
    ggtitle("D)")+scale_y_continuous(limits = c(0, 2000)), ncol=2)

##Figure 2: Net changes of soil microbial biomass (A) carbon (MBC), (B) 
##nitrogen (MBN) and (C) phosphorus (MBP) after 48-hour incubation of two 
##spruce forest soils (PL - Plešné and CT - Čertovo catchments) mixed at 
##five different ratios (i.e., 0:1, 0.25:0.75, 0.5:0.5, 0.75:0.25, and 1:0 
##in respect to PL). Grey circles denote the litter horizon and empty circles 
##the topsoil organic horizon. Symbols show mean values and error bars 
##standard error of the mean (n = 4). Notice that plots have different 
##y-axes scales.

mix_diff<-mix[c(1:80), c("Plesne", "Certovo", "Horizon", "Labelling")]
mix_diff$dCmic<-mix[c(81:160), c("Cmic")]-mix[c(1:80), c("Cmic")]
mix_diff$dNmic<-mix[c(81:160), c("Nmic")]-mix[c(1:80), c("Nmic")]
mix_diff$dPmic<-mix[c(81:160), c("Pmic")]-mix[c(1:80), c("Pmic")]

grid.arrange(
  mix_diff %>% filter(Labelling=="NO") %>% group_by(Plesne, Horizon) %>% 
    summarize(y.sd=sd(dCmic)/sqrt(5), y=mean(dCmic)) %>%
    ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=Horizon),
                                              position = position_dodge(width = 0.6))+
    theme_min+
    geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd, colour=Horizon), width=0.1, lwd=0.5,
                  position = position_dodge(width = 0.6))+
    ylab(expression(paste(Delta~MBC, " (", mu, "mol ",g^{-1}, ")" )))+
    xlab("PL : CT mixing ratio")+scale_fill_manual(values = c("grey", "white"))+
    scale_color_manual(values = c("black", "black"))+scale_y_continuous(breaks = c(-150, -100, -50, 0, 50, 100))+
    geom_hline(yintercept = 0, lwd=1)+theme(legend.position = c(0.65, 0.88),
                                            legend.key.size = unit(0.3, "in"),
                                            legend.title = element_blank())+
    ggtitle("A)"),
  mix_diff %>% filter(Labelling=="NO") %>% group_by(Plesne, Horizon) %>% 
    summarize(y.sd=sd(dNmic)/sqrt(5), y=mean(dNmic)) %>%
    ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=Horizon), show.legend = F,
                                              position = position_dodge(width = 0.6))+
    theme_min+theme(legend.position = c(0.2, 0.2))+
    geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd, colour=Horizon), width=0.1, lwd=0.5,
                  position = position_dodge(width = 0.6), show.legend = F)+
    ylab(expression(paste(Delta~MBN, " (", mu, "mol ",g^{-1}, ")" )))+
    xlab("PL : CT mixing ratio")+scale_fill_manual(values = c("grey", "white"))+
    scale_color_manual(values = c("black", "black"))+#scale_y_continuous(breaks = c(-150, -100, -50, 0, 50, 100))+
    geom_hline(yintercept = 0, lwd=1)+theme(legend.position = c(0.8, 0.8),
                                            legend.key.size = unit(0.3, "in"),
                                            legend.title = element_blank())+
    ggtitle("B)"),
  mix_diff %>% filter(Labelling=="NO") %>% group_by(Plesne, Horizon) %>% 
    summarize(y.sd=sd(dPmic)/sqrt(5), y=mean(dPmic)) %>%
    ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=Horizon), show.legend = F,
                                              position = position_dodge(width = 0.6))+
    theme_min+theme(legend.position = c(0.2, 0.2))+
    geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd, colour=Horizon), width=0.1, lwd=0.5,
                  position = position_dodge(width = 0.6), show.legend = F)+
    ylab(expression(paste(Delta~MBP, " (", mu, "mol ",g^{-1}, ")" )))+
    xlab("PL : CT mixing ratio")+scale_fill_manual(values = c("grey", "white"))+
    scale_color_manual(values = c("black", "black"))+#scale_y_continuous(breaks = c(-150, -100, -50, 0, 50, 100))+
    geom_hline(yintercept = 0, lwd=1)+theme(legend.position = c(0.8, 0.8),
                                            legend.key.size = unit(0.3, "in"),
                                            legend.title = element_blank())+
    ggtitle("C)"), ncol=3)

##Figure 3: Isotopic signals of respired CO2 (R, black symbols), K2SO4 
##extractable organic carbon (K2SO4 - EC, grey symbols), and microbial 
##biomass carbon (MBC, empty symbols) in the litter and topsoil organic 
##horizons of two spruce forest soils (PL - Plešné and CT - Čertovo catchments) 
##mixed at five different ratios (i.e., 0:1, 0.25:0.75, 0.5:0.5, 0.75:0.25, 
##and 1:0 in respect to PL). Symbols show mean values and error bars standard 
##error of the mean (n = 4). The horizontal dashed line represents an 
##approximation of the isotope signal of organic compounds consumed by the 
##microbial biomass. The solid arrow shows the change in the isotopic signal 
##of respired CO2 across soil mixtures that was used to calculate fMBC 
##(see section 2.5. for details).

izo<-mix[c(81:160), c("Plesne", "Certovo", "Horizon", "Labelling", "CO2atm",
                   "DOCatm", "Cmicatm")]
Izo<-melt(izo, id.vars = c("Plesne", "Certovo", "Horizon", "Labelling"))
###Calculate delta values
Izo$delta<-with(Izo, (value/(1-value)/0.011237-1)*1000)

St<-as.numeric(Izo %>% filter(Plesne==0 & Horizon=="Litter" & variable=="CO2atm" & Labelling=="NO") %>%
  summarize(mean(delta, na.rm=T)))
St2<-as.numeric(Izo %>% filter(Plesne==1 & Horizon=="Litter" & variable=="CO2atm" & Labelling=="NO") %>%
                 summarize(mean(delta, na.rm=T)))
Izo$St<-NA
Izo[(Izo$Horizon=="Litter"), "St"]<-St

Izo$St2<-NA
Izo[(Izo$Horizon=="Litter"), "St2"]<-St2

Izo %>% filter(Labelling=="NO") %>%
  group_by(Plesne, Horizon, variable) %>% 
  summarize(y.sd=sd(delta, na.rm = T)/sqrt(5), y=mean(delta, na.rm = T)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=variable))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(.~Horizon)+theme_min+
  theme(legend.position = c(0.2, 0.85),
        legend.title = element_blank(),
        legend.text.align = 0)+
  geom_hline(data=Izo, aes(yintercept=St), lwd=1, lty=2, color="grey30")+
  geom_segment(data=Izo, aes(x=1, xend=5, y=St, yend=St2),
               arrow = arrow(length = unit(0.5, "cm")),
               color="grey30")+
  ylab(expression(paste(delta^{13},C)))+
  xlab("PL : CT mixing ratio")+
  scale_fill_manual(values = c("black", "grey", "white"),
                    name = '',
                    labels = expression(R, K[2]~SO[4]-EC, MBC))

sd(Izo[c(Izo$variable=="CO2atm" & Izo$Horizon=="Organic topsoil" &
             Izo$Labelling=="NO"), "delta"]-
       Izo[c(Izo$variable=="DOCatm" & Izo$Horizon=="Organic topsoil" &
               Izo$Labelling=="NO"), "delta"])

##Figure 4: Net changes of water extractable mineral nitrogen (ΔMN) in the 
##litter and topsoil organic horizons of two spruce forest soils 
##(PL - Plešné and CT - Čertovo catchments) mixed at five different ratios 
##(i.e., 0:1, 0.25:0.75, 0.5:0.5, 0.75:0.25, and 1:0 in respect to PL). 
##White bars represent measured ΔMN, black and grey bars represent ΔMN 
##calculated using eq. 1 with CUE estimated using isotopic approach 
##(denoted as CUE1) or mass balance approach (denoted as CUE2), 
##respectively. ΔMN was calculated with (A) and without (B) fMBC 
##(see section 2.5. for details). Error bars represent standard error of 
##the mean (n = 4). Correspondence between the predictions and observations 
##are reported as log likelihood (LL, eq. 10). 
##The higher the LL (less negative), the better the correspondence.

###Assuming no contribution of decaying microbial biomass
mix_diff$CUE<-mix[c(81:160), c("CUE")]
mix_diff$CO2<-mix[c(81:160), c("CCO2")]
mix_diff$DOC<-mix[c(1:80), c("DOC")]
mix_diff$DON<-mix[c(1:80), c("DON")]
mix_diff$Cmic<-mix[c(1:80), c("Cmic")]
mix_diff$Nmic<-mix[c(1:80), c("Nmic")]
mix_diff$dNH4<-mix[c(81:160), c("NH4")]-mix[c(1:80), c("NH4")]
mix_diff$dNO3<-mix[c(81:160), c("NO3")]-mix[c(1:80), c("NO3")]

###~~~~~~~~~~~~~~~~MASS BALANCE ESTIMATE OF CUE~~~~~~~~~~~~~~~~~~~~~~~~~~###
source("./manuscript_data_and_code/Cuest.R")

cl<-makeCluster(4)
registerDoParallel(cl)

Cuest_out<-Cuest(data = mix)

stopImplicitCluster()

for(i in 1:40){
  mix_diff$CUEbalance[i]<-Cuest_out[[i]]$pars[["CUE"]]
}
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

####Amount of organic carbon utilized
mix_diff$U1<-with(mix_diff, CO2/(1-CUE))#labelled glucose
mix_diff$U2<-with(mix_diff, CO2/(1-CUEbalance))#mass balance
####Expected amount of mineral nitrogen exchanged
mix_diff$dNm_preda1<-with(mix_diff, U1*(DON/(DOC+Gl)-Nmic*CUE/Cmic))
mix_diff$dNm_preda2<-with(mix_diff, U2*(DON/(DOC+Gl)-Nmic*CUEbalance/Cmic))

###Assuming contribution of decaying microbial biomass
####Calculate relative contribution of decaying microbial biomass to food source
####of the growing microbial biomass
mix_diff$deltaCO2<-Izo[Izo$variable=="CO2atm", "delta"]
mix_diff$deltaCmic<-Izo[Izo$variable=="Cmicatm", "delta"]
mix_diff$deltaDOC<-Izo[Izo$variable=="DOCatm", "delta"]
mix_diff$deltaBase<-mean(Izo[(Izo$variable=="CO2atm" & Izo$Horizon=="Litter" &
                                Izo$Plesne==0 & Izo$Labelling=="NO"), "delta"], na.rm=T)
mix_diff$Cmic.prop<-with(mix_diff, (deltaCO2-deltaBase)/(deltaCmic-deltaBase))
mix_diff$Cmic.prop2<-with(mix_diff, (deltaCO2+0.83-deltaDOC)/(deltaCmic-deltaDOC))


####Organic topsoil Cmic.prop is 0
mix_diff[(mix_diff$Horizon=="Organic topsoil"), "Cmic.prop"]<-0
####All Cmic.prop higher than 1 and lower than 0 are set to be 1 and 0 respectively
mix_diff[(mix_diff$Cmic.prop>1), "Cmic.prop"]<-1
mix_diff[(mix_diff$Cmic.prop<0), "Cmic.prop"]<-0
mix_diff[c(41:80), "Cmic.prop"]<-mix_diff[c(1:40), "Cmic.prop"]
mix_diff[(mix_diff$Cmic.prop2>1), "Cmic.prop2"]<-1
mix_diff[(mix_diff$Cmic.prop2<0), "Cmic.prop2"]<-0
mix_diff[c(41:80), "Cmic.prop2"]<-mix_diff[c(1:40), "Cmic.prop2"]


mix_diff$dNm_predb1<-with(mix_diff, U1*(((1-Cmic.prop)*(DON/(DOC))+Cmic.prop*(Nmic/(Cmic)))-Nmic*CUE/(Cmic)))
mix_diff$dNm_predb2<-with(mix_diff, U2*(((1-Cmic.prop)*(DON/(DOC))+Cmic.prop*(Nmic/(Cmic)))-Nmic*CUEbalance/(Cmic)))


#Log likelihoods
##prediction a
with(subset(mix_diff, Labelling=="NO"),
     -2/length(dNH4)*log(2*pi*sd(dNH4+dNO3)^2) - 
       sum((dNH4+dNO3-dNm_preda1)^2/2/sd(dNH4+dNO3)^2))
with(subset(mix_diff, Labelling=="NO"),
     -2/length(dNH4)*log(2*pi*sd(dNH4+dNO3)^2) - 
       sum((dNH4+dNO3-dNm_preda2)^2/2/sd(dNH4+dNO3)^2))
##prediction b
with(subset(mix_diff, Labelling=="NO"),
     -2/length(dNH4)*log(2*pi*sd(dNH4+dNO3)^2) - 
       sum((dNH4+dNO3-dNm_predb1)^2/2/sd(dNH4+dNO3)^2))
with(subset(mix_diff, Labelling=="NO"),
     -2/length(dNH4)*log(2*pi*sd(dNH4+dNO3)^2) - 
       sum((dNH4+dNO3-dNm_predb2)^2/2/sd(dNH4+dNO3)^2))

#likelihood ratio test
-2*(with(subset(mix_diff, Labelling=="NO"),
         -2/length(dNH4)*log(2*pi*sd(dNH4+dNO3)^2) - 
           sum((dNH4+dNO3-dNm_preda1)^2/2/sd(dNH4+dNO3)^2)) - 
      with(subset(mix_diff, Labelling=="NO"),
           -2/length(dNH4)*log(2*pi*sd(dNH4+dNO3)^2) - 
             sum((dNH4+dNO3-dNm_predb1)^2/2/sd(dNH4+dNO3)^2)))

pchisq(-2*(with(subset(mix_diff, Labelling=="NO"),
                -2/length(dNH4)*log(2*pi*sd(dNH4+dNO3)^2) - 
                  sum((dNH4+dNO3-dNm_preda1)^2/2/sd(dNH4+dNO3)^2)) - 
             with(subset(mix_diff, Labelling=="NO"),
                  -2/length(dNH4)*log(2*pi*sd(dNH4+dNO3)^2) - 
                    sum((dNH4+dNO3-dNm_predb1)^2/2/sd(dNH4+dNO3)^2))), df=1, lower.tail=FALSE)

#likelihood ratio test
-2*(with(subset(mix_diff, Labelling=="NO"),
         -2/length(dNH4)*log(2*pi*sd(dNH4+dNO3)^2) - 
           sum((dNH4+dNO3-dNm_preda2)^2/2/sd(dNH4+dNO3)^2)) - 
      with(subset(mix_diff, Labelling=="NO"),
           -2/length(dNH4)*log(2*pi*sd(dNH4+dNO3)^2) - 
             sum((dNH4+dNO3-dNm_predb2)^2/2/sd(dNH4+dNO3)^2)))

pchisq(-2*(with(subset(mix_diff, Labelling=="NO"),
                -2/length(dNH4)*log(2*pi*sd(dNH4+dNO3)^2) - 
                  sum((dNH4+dNO3-dNm_preda2)^2/2/sd(dNH4+dNO3)^2)) - 
             with(subset(mix_diff, Labelling=="NO"),
                  -2/length(dNH4)*log(2*pi*sd(dNH4+dNO3)^2) - 
                    sum((dNH4+dNO3-dNm_predb2)^2/2/sd(dNH4+dNO3)^2))), df=1, lower.tail=FALSE)

ann_texta1 <- data.frame(Plesne = 0.75, y = -2, lab = "LL[(CUE[1])]==-254",
                        Horizon = factor("Organic topsoil",
                                         levels = c("Organic topsoil", "Litter")))
ann_texta2 <- data.frame(Plesne = 0.75, y = -3, lab = "LL[(CUE[2])]==-163",
                         Horizon = factor("Organic topsoil",
                                          levels = c("Organic topsoil", "Litter")))


ann_textb1 <- data.frame(Plesne = 0.75, y = -2, lab = "LL[(CUE[1])]==-69",
                        Horizon = factor("Organic topsoil",
                                         levels = c("Organic topsoil", "Litter")))
ann_textb2 <- data.frame(Plesne = 0.75, y = -3, lab = "LL[(CUE[2])]==-81",
                         Horizon = factor("Organic topsoil",
                                          levels = c("Organic topsoil", "Litter")))


mix_diff$Mn<-mix_diff$dNH4+mix_diff$dNO3

Nplot_data1<-melt(mix_diff[mix_diff$Labelling=="NO", 
                          c("Plesne", "Horizon", "Mn", "dNm_preda1", "dNm_preda2")],
                  id.vars = c("Plesne", "Horizon"))
Nplot_data1_labs<-c(expression(Measured), 
                    expression(Predicted~-~CUE[1]),
                    expression(Predicted~-~CUE[2]))

(Nplot1<-Nplot_data1[-5, ] %>% group_by(Plesne, Horizon, variable) %>% 
  summarize(y.sd=sd(value)/sqrt(5), y=mean(value)) %>%
  ggplot(aes(factor(Plesne), y)) + 
  geom_bar(aes(fill=variable), stat = "identity", position = "dodge", color="black", show.legend=F) +
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd, color=variable), position="dodge", show.legend=F)+
  scale_fill_manual(labels = Nplot_data1_labs, values = c("white", "grey", "black"))+
  scale_color_manual(values = rep("black", 3))+
  facet_grid(~Horizon) + theme_min +
  theme(legend.position = c(0.8, 0.15),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.text.align = c(0))+
  ylab(expression(paste(italic(Delta~M[N]), " (", mu, "mol ",g^{-1}, ")" )))+
  xlab("PL : CT mixing ratio") + ggtitle("A)") +
  geom_text(data = ann_texta1, label=ann_texta1$lab,
            fontface="italic", size=7, parse = T)+
  geom_text(data = ann_texta2, label=ann_texta2$lab,
            fontface="italic", size=7, color="grey30", parse = T))

Nplot_data2<-melt(mix_diff[mix_diff$Labelling=="NO", 
                           c("Plesne", "Horizon", "Mn", "dNm_predb1", "dNm_predb2")],
                  id.vars = c("Plesne", "Horizon"))

(Nplot2<-Nplot_data2[-5, ] %>% group_by(Plesne, Horizon, variable) %>% 
  summarize(y.sd=sd(value)/sqrt(5), y=mean(value)) %>%
  ggplot(aes(factor(Plesne), y)) + 
  geom_bar(aes(fill=variable), stat = "identity", position = "dodge", color="black") +
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd, color=variable), position="dodge", show.legend=F)+
  scale_fill_manual(labels = Nplot_data1_labs, values = c("white", "grey", "black"))+
  scale_color_manual(values = rep("black", 3))+
  facet_grid(~Horizon) + theme_min +
  theme(legend.position = c(0.15, 0.25),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank(),
        legend.text.align = c(0))+
  ylab(expression(paste(italic(Delta~M[N]), " (", mu, "mol ",g^{-1}, ")" )))+
  xlab("PL : CT mixing ratio") + ggtitle("B)") +
  geom_text(data = ann_textb1, label=ann_textb1$lab,
            fontface="italic", size=7, parse = T)+
  geom_text(data = ann_textb2, label=ann_textb2$lab,
            fontface="italic", size=7, color="grey30", parse = T)+
  scale_y_continuous(limits = c(-5, 1), breaks = c(-5, -4, -3, -2, -1, 0, 1)))




grid.arrange(Nplot1, Nplot2, nrow=2)

##Figure 5: Net changes of water extractable soluble reactive phosphorus 
##(ΔMP) in the litter and topsoil organic horizons of two spruce forest soils
##(PL - Plešné and CT - Čertovo catchments) mixed at five different ratios 
##(i.e., 0:1, 0.25:0.75, 0.5:0.5, 0.75:0.25, and 1:0 in respect to PL). 
##White bars represent measured ΔMP, black and grey bars represent ΔMP 
##calculated using eq. 1 with CUE estimated using isotopic approach 
##(denoted as CUE1) or mass balance approach (denoted as CUE2), respectively.
##ΔMP was calculated with (A) and without (B) fMBC (see section 2.5. for 
##details). Error bars represent standard error of the mean (n = 4). 
##Correspondences between the predictions and observations are reported as 
##log likelihood (LL, eq. 10). The higher the LL (less negative), 
##the better the correspondence.

###Assuming no contribution of decaying microbial biomass
mix_diff$DOP<-mix[c(1:80), c("DOP")]
mix_diff$Pmic<-mix[c(1:80), c("Pmic")]
mix_diff$PO4<-mix[c(1:80), c("PO4")]
mix_diff$dPO4<-mix[c(81:160), c("PO4")]-mix[c(1:80), c("PO4")]

mix_diff$dPm_preda1<-with(mix_diff, U1*(DOP/DOC-Pmic*CUE/Cmic))
mix_diff$dPm_preda2<-with(mix_diff, U2*(DOP/DOC-Pmic*CUEbalance/Cmic))

###Assuming contribution of decaying microbial biomass
mix_diff$dPm_predb1<-with(mix_diff, U1*(((1-Cmic.prop)*(DOP/DOC)+Cmic.prop*(Pmic/Cmic))-Pmic*CUE/Cmic))
mix_diff$dPm_predb2<-with(mix_diff, U2*(((1-Cmic.prop)*(DOP/DOC)+Cmic.prop*(Pmic/Cmic))-Pmic*CUEbalance/Cmic))


#Log likelihoods
##prediction a
with(subset(mix_diff, Labelling=="NO"),
     -2/length(dPO4)*log(2*pi*sd(dPO4)^2) - 
       sum((dPO4-dPm_preda1)^2/2/sd(dPO4)^2))
with(subset(mix_diff, Labelling=="NO"),
     -2/length(dPO4)*log(2*pi*sd(dPO4)^2) - 
       sum((dPO4-dPm_preda2)^2/2/sd(dPO4)^2))
##prediction b
with(subset(mix_diff, Labelling=="NO"),
     -2/length(dPO4)*log(2*pi*sd(dPO4)^2) - 
       sum((dPO4-dPm_predb1)^2/2/sd(dPO4)^2))
with(subset(mix_diff, Labelling=="NO"),
     -2/length(dPO4)*log(2*pi*sd(dPO4)^2) - 
       sum((dPO4-dPm_predb2)^2/2/sd(dPO4)^2))
# ##prediction c
# with(subset(mix_diff, Labelling=="NO"),
#      -2/length(dPO4)*log(2*pi*sd(dPO4)^2) - 
#        sum((dPO4-PO4sorbed)^2/2/sd(dPO4)^2))

#likelihood ratio test
-2*(with(subset(mix_diff, Labelling=="NO"),
         -2/length(dPO4)*log(2*pi*sd(dPO4)^2) - 
           sum((dPO4-dPm_preda1)^2/2/sd(dPO4)^2)) - 
      with(subset(mix_diff, Labelling=="NO"),
           -2/length(dPO4)*log(2*pi*sd(dPO4)^2) - 
             sum((dPO4-dPm_predb1)^2/2/sd(dPO4)^2)))

pchisq(-2*(with(subset(mix_diff, Labelling=="NO"),
                -2/length(dPO4)*log(2*pi*sd(dPO4)^2) - 
                  sum((dPO4-dPm_preda1)^2/2/sd(dPO4)^2)) - 
             with(subset(mix_diff, Labelling=="NO"),
                  -2/length(dPO4)*log(2*pi*sd(dPO4)^2) - 
                    sum((dPO4-dPm_predb1)^2/2/sd(dPO4)^2))), df=1, lower.tail=FALSE)

-2*(with(subset(mix_diff, Labelling=="NO"),
         -2/length(dPO4)*log(2*pi*sd(dPO4)^2) - 
           sum((dPO4-dPm_preda2)^2/2/sd(dPO4)^2)) - 
      with(subset(mix_diff, Labelling=="NO"),
           -2/length(dPO4)*log(2*pi*sd(dPO4)^2) - 
             sum((dPO4-dPm_predb2)^2/2/sd(dPO4)^2)))

pchisq(-2*(with(subset(mix_diff, Labelling=="NO"),
                -2/length(dPO4)*log(2*pi*sd(dPO4)^2) - 
                  sum((dPO4-dPm_preda2)^2/2/sd(dPO4)^2)) - 
             with(subset(mix_diff, Labelling=="NO"),
                  -2/length(dPO4)*log(2*pi*sd(dPO4)^2) - 
                    sum((dPO4-dPm_predb2)^2/2/sd(dPO4)^2))), df=1, lower.tail=FALSE)

# -2*(with(subset(mix_diff, Labelling=="NO"),
#          -2/length(dPO4)*log(2*pi*sd(dPO4)^2) - 
#            sum((dPO4-dPm_preda)^2/2/sd(dPO4)^2)) - 
#       with(subset(mix_diff, Labelling=="NO"),
#            -2/length(dPO4)*log(2*pi*sd(dPO4)^2) - 
#              sum((dPO4-PO4sorbed)^2/2/sd(dPO4)^2)))
# 
# pchisq(-2*(with(subset(mix_diff, Labelling=="NO"),
#                 -2/length(dPO4)*log(2*pi*sd(dPO4)^2) - 
#                   sum((dPO4-dPm_preda)^2/2/sd(dPO4)^2)) - 
#              with(subset(mix_diff, Labelling=="NO"),
#                   -2/length(dPO4)*log(2*pi*sd(dPO4)^2) - 
#                     sum((dPO4-PO4sorbed)^2/2/sd(dPO4)^2))), df=1, lower.tail=FALSE)
# 
# -2*(with(subset(mix_diff, Labelling=="NO"),
#          -2/length(dPO4)*log(2*pi*sd(dPO4)^2) - 
#            sum((dPO4-dPm_predb)^2/2/sd(dPO4)^2)) - 
#       with(subset(mix_diff, Labelling=="NO"),
#            -2/length(dPO4)*log(2*pi*sd(dPO4)^2) - 
#              sum((dPO4-PO4sorbed)^2/2/sd(dPO4)^2)))
# 
# pchisq(-2*(with(subset(mix_diff, Labelling=="NO"),
#                 -2/length(dPO4)*log(2*pi*sd(dPO4)^2) - 
#                   sum((dPO4-dPm_predb)^2/2/sd(dPO4)^2)) - 
#              with(subset(mix_diff, Labelling=="NO"),
#                   -2/length(dPO4)*log(2*pi*sd(dPO4)^2) - 
#                     sum((dPO4-PO4sorbed)^2/2/sd(dPO4)^2))), df=1, lower.tail=FALSE)


p_texta1 <- data.frame(Plesne = 0.5, y = -3, lab = "LL[(CUE[1])]==-67~.~10^{4}",
                        Horizon = factor("Organic topsoil",
                                         levels = c("Organic topsoil", "Litter")))
p_texta2 <- data.frame(Plesne = 0.5, y = -3.5, lab = "LL[(CUE[2])]==-130~.~10**4",
                       Horizon = factor("Organic topsoil",
                                        levels = c("Organic topsoil", "Litter")))
p_textb1 <- data.frame(Plesne = 0.5, y = -3, lab = "LL[(CUE[1])]==-37~.~10^{4}",
                        Horizon = factor("Organic topsoil",
                                         levels = c("Organic topsoil", "Litter")))
p_textb2 <- data.frame(Plesne = 0.5, y = -3.5, lab = "LL[(CUE[2])]==-110~.~10^{4}",
                      Horizon = factor("Organic topsoil",
                                       levels = c("Organic topsoil", "Litter")))

Pplot_data1<-melt(mix_diff[mix_diff$Labelling=="NO", 
                           c("Plesne", "Horizon", "dPO4", "dPm_preda1", "dPm_preda2")],
                  id.vars = c("Plesne", "Horizon"))

(Pplot1<-Pplot_data1[-5, ] %>% group_by(Plesne, Horizon, variable) %>% 
    summarize(y.sd=sd(value)/sqrt(5), y=mean(value)) %>%
    ggplot(aes(factor(Plesne), y)) + 
    geom_bar(aes(fill=variable), stat = "identity", position = "dodge", color="black", show.legend=F) +
    geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd, color=variable), position="dodge", show.legend=F)+
    scale_fill_manual(labels = Nplot_data1_labs, values = c("white", "grey", "black"))+
    scale_color_manual(values = rep("black", 3))+
    facet_grid(~Horizon) + theme_min +
    theme(legend.position = c(0.8, 0.15),legend.key.size = unit(0.3, "in"),
          legend.title = element_blank(),
          axis.title.x = element_blank(),
          legend.text.align = c(0))+
    ylab(expression(paste(italic(Delta~M[P]), " (", mu, "mol ",g^{-1}, ")" )))+
    xlab("PL : CT mixing ratio") + ggtitle("A)") +
    geom_text(data = p_texta1, label=p_texta1$lab,
              fontface="italic", size=7, parse = T)+
    geom_text(data = p_texta2, label=p_texta2$lab,
              fontface="italic", size=7, color="grey30", parse = T))

Pplot_data2<-melt(mix_diff[mix_diff$Labelling=="NO", 
                           c("Plesne", "Horizon", "dPO4", "dPm_predb1", "dPm_predb2")],
                  id.vars = c("Plesne", "Horizon"))

(Pplot2<-Pplot_data2[-5, ] %>% group_by(Plesne, Horizon, variable) %>% 
    summarize(y.sd=sd(value)/sqrt(5), y=mean(value)) %>%
    ggplot(aes(factor(Plesne), y)) + 
    geom_bar(aes(fill=variable), stat = "identity", position = "dodge", color="black") +
    geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd, color=variable), position="dodge", show.legend=F)+
    scale_fill_manual(labels = Nplot_data1_labs, values = c("white", "grey", "black"))+
    scale_color_manual(values = rep("black", 3))+
    facet_grid(~Horizon) + theme_min +
    theme(legend.position = c(0.15, 0.25),legend.key.size = unit(0.3, "in"),
          legend.title = element_blank(),
          legend.text.align = c(0))+
    ylab(expression(paste(italic(Delta~M[P]), " (", mu, "mol ",g^{-1}, ")" )))+
    xlab("PL : CT mixing ratio") + ggtitle("B)") +
    geom_text(data = p_textb1, label=p_textb1$lab,
              fontface="italic", size=7, parse = T)+
    geom_text(data = p_textb2, label=p_textb2$lab,
              fontface="italic", size=7, color="grey30", parse = T))




grid.arrange(Pplot1, Pplot2, nrow=2)

##Figure 6: A) Expected growth rates of the soil microbial community in the 
##litter (grey circles) and organic topsoil (empty circles) horizons without 
##an external source of phosphorus (i.e., phosphors demand is covered 
##exclusively from internal resources). Symbols show mean values and error 
##bars standard error of the mean (n = 4). The expected growth rate (solid 
##black line) was calculated from measured MBC:MBP ratios using a previously
##derived relationship between MBC:MBP ratios and growth rate in the absence
##of external sources of phosphorus (light grey points; Čapek et al., 2016).
##The horizontal solid line denotes a zero growth rate. B) Net changes of 
##water extractable soluble reactive phosphorus (ΔMP) in the litter and 
##topsoil organic horizons of two spruce forest soils (PL - Plešné and 
##CT - Čertovo catchments) mixed at five different ratios 
##(i.e., 0:1, 0.25:0.75, 0.5:0.5, 0.75:0.25, and 1:0 in respect to PL). 
##White bars represent measured ΔMP, black bars represent theoretical 
##amount of chemically adsorbed SRP. Error bars represent standard error of 
##the mean (n = 4).

####Read data
gr<-read.csv("./manuscript_data_and_code/capek2016_data_SI.csv")

####Estimate the growth rate for all sample
cpred<-coef(lm(u0~CP, gr))
mix_diff$grate_pred<-with(mix_diff, cpred[1]+(Cmic/Pmic)*cpred[2])

gr_pred<-mix_diff %>% filter(Labelling=="NO") %>%
  group_by(Plesne, Horizon) %>% summarize(CP=mean(Cmic/Pmic, na.rm = T),
                                          CP.sd=sd(Cmic/Pmic, na.rm = T),
                                          growth=mean(grate_pred, na.rm = T),
                                          growth.sd=sd(grate_pred, na.rm = T))
(Pplot3<-ggplot(data=gr_pred, aes(CP, growth))+geom_point(cex=6, pch=21, aes(fill=Horizon))+
  theme_min+
  xlab("MBC:MBP (mol:mol)")+
  ylab(expression(paste(mu[0])))+
  geom_hline(yintercept = 0, color='grey')+
  geom_errorbarh(aes(xmin=CP-CP.sd, xmax=CP+CP.sd))+
  geom_errorbar(aes(ymin=growth-growth.sd, ymax=growth+growth.sd))+
  geom_point(data=gr, cex=6, color='grey', alpha=0.5, aes(CP, u0))+
  stat_smooth(data=gr, aes(CP, u0), method = lm, se=F, color='black', lwd=0.5, fullrange = T)+
  scale_fill_manual(values = c("grey", "white"))+
  theme(legend.position = c(0.8, 0.8),
        legend.title = element_blank())+
  scale_x_continuous(breaks = c(20, 40, 60, 80, 100))+
  scale_y_continuous(breaks=c(-1, -0.5, 0, 0.5, 1, 1.5, 2))+ ggtitle("A)"))

###Assuming independence of microbial growth on soil derived phosphorus
###water extractable phosphate is chemicaly sorbed instead 
mix_diff$Sorption<-mix[c(1:80), c("Sorption")]
####All values greater than 1 is set to 1
mix_diff[mix_diff$Sorption>1, "Sorption"]<-1

####Calculating the amount of chemicaly sorbed phosphate
mix_diff$PO4sorbed<-with(mix_diff, -PO4*(1-exp(-(1-Sorption)/0.75*48)))

Pplot_data3<-melt(mix_diff[mix_diff$Labelling=="NO", 
                           c("Plesne", "Horizon", "dPO4", "PO4sorbed")],
                  id.vars = c("Plesne", "Horizon"))
Pplot_data3_labs<-c(expression(Measured~M[P]), 
                    expression(Adsorbed~P-PO[4]))

(Pplot4<-Pplot_data3 %>% group_by(Plesne, Horizon, variable) %>% 
    summarize(y.sd=sd(value)/sqrt(5), y=mean(value)) %>%
    ggplot(aes(factor(Plesne), y)) + 
    geom_bar(aes(fill=variable), stat = "identity", position = "dodge", color="black") +
    geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd, color=variable), position="dodge", show.legend=F)+
    scale_fill_manual(labels = Pplot_data3_labs, values = c("white", "black"))+
    scale_color_manual(values = rep("black", 3))+
    facet_grid(~Horizon) + theme_min +
    theme(legend.position = c(0.25, 0.18),legend.key.size = unit(0.3, "in"),
          legend.title = element_blank(),
          legend.text.align = c(0))+
    ylab(expression(paste(italic(Delta~M[P]), " (", mu, "mol ",g^{-1}, ")" )))+
    xlab("PL : CT mixing ratio") + ggtitle("B)"))

grid.arrange(Pplot3, Pplot4, nrow=1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Supplementary information
##Stoichiometry
##Figure S1: Molar (A) C/N and (B) C/P ratios of plant litterfall in Plešné (red circles) and 
##Čertovo (black circles) catchments. The chemical composition of litterfall was measured at 
##the localities with aboveground vegetation unaffected by bark beetle infestation. Symbols 
##represent medians and error bars lower and upper quartiles. The data shown in this figure 
##can be downloaded at 
##https://link.springer.com/article/10.1007/s10533-018-0470-x#SupplementaryMaterial. 

#Please specify the path of the downloaded file
fluxes<-read.xlsx(xlsxFile = c("./10533_2018_470_MOESM7_ESM.xlsx"), 2)

grid.arrange(
  fluxes %>% filter(Locality == "Lower") %>% group_by(Catchment, Year) %>%
    summarize(y=median(Clitter/Nlitter, na.rm = T),
              yl=quantile(Clitter/Nlitter, na.rm = T, probs=0.25),
              yu=quantile(Clitter/Nlitter, na.rm = T, probs=0.75)) %>%
    ggplot(aes(Year, y))+geom_point(cex=6, pch=21, aes(fill=Catchment))+
    geom_errorbar(aes(ymin=yl, ymax=yu, color=Catchment), width=0.1)+
    theme_min+xlim(2010, 2016)+
    scale_fill_manual(values = c("black", "red"))+
    scale_color_manual(values = c("black", "red"))+
    theme(legend.title = element_blank(),
          legend.position = c(0.15, 0.8))+
    ylab("C/N of plant litter input (mol/mol)")+
    ggtitle("A"),
  fluxes %>% filter(Locality == "Lower") %>% group_by(Catchment, Year) %>%
    summarize(y=median(Clitter/Plitter, na.rm = T), 
              yl=quantile(Clitter/Plitter, na.rm = T, probs=0.25),
              yu=quantile(Clitter/Plitter, na.rm = T, probs=0.75)) %>%
    ggplot(aes(Year, y))+geom_point(cex=6, pch=21, aes(fill=Catchment), show.legend = F)+
    geom_errorbar(aes(ymin=yl, ymax=yu, color=Catchment), width=0.1, show.legend = F)+
    theme_min+xlim(2010, 2016)+
    scale_fill_manual(values = c("black", "red"))+
    scale_color_manual(values = c("black", "red"))+
    ylab("C/P of plant litter input (mol/mol)")+
    ggtitle("B"), ncol=1)

##Figure S2: Molar (A) DOC/DON (water extractable organic carbon to nitrogen), (B) DOC/DOP 
##(water extractable organic carbon to phosphorus), (C) MBC/MBN ( microbial biomass carbon to 
##nitrogen), (D) MBC/MBP (microbial biomass carbon to phosphorus) and (E) MBN/MBN (microbial 
##biomass nitrogen to phosphorus) ratios in the litter (empty symbols) and topsoil organic 
##horizons (grey symbols) of two spruce forest soils (Plešné and Čertovo catchments) mixed 
##at five different ratios (i.e., 0:1, 0.25:0.75, 0.5:0.5, 0.75:0.25, and 1:0 in respect to Plešné). 
##Triangles and circles denote the values measured at the beginning and at the end of the incubation, 
##respectively. Symbols represent mean values and error bars standard deviation of the mean (n = 4). 
##The effect of soil horizon (Horizon), the proportion between Plešné and Čertovo in soil mixture 
##(Plesne : Certovo) and time on respective ratios is reported. 
##Levels of significance: ***, p < 0.001; **, p < 0.01; *, p < 0.05; n.s., not significant.

mix$Legend<-ifelse(mix$TIME==0, "Before incubation", "After incubation")
anova(glm(DOC/DON~Horizon+Plesne/Horizon+TIME/Plesne/Horizon, 
          data = subset(mix, Labelling=="NO"),
          family = Gamma), test="F")
mix %>% filter(Labelling=="NO") %>% 
  group_by(Plesne, Horizon, Legend) %>% summarize(y.sd=sd(DOC/DON), y=mean(DOC/DON)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, aes(fill=Horizon, shape=Legend), show.legend = F)+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme_min+theme(legend.position = "right",
                  legend.title = element_blank())+
  #facet_wrap(~Horizon, scales="free")+
  ylab(expression(paste("DOC/DON (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  annotate("text", label="Horizon***", 0.5, 45, hjust="left", size=6, fontface="italic")+
  annotate("text", label="Plesne : Certovo n.s.", 0.5, 42, hjust="left", size=6, fontface="italic")+
  annotate("text", label="Time n.s.", 0.5, 39, hjust="left", size=6, fontface="italic")+
  ggtitle("A")+scale_shape_manual(values = c(21, 24))
anova(glm(DOC/DOP~Horizon+Plesne/Horizon+TIME/Plesne/Horizon, 
          data = subset(mix, Labelling=="NO"),
          family = Gamma), test="F")
mix %>% filter(Labelling=="NO") %>% 
  group_by(Plesne, Horizon, Legend) %>% summarize(y.sd=sd(DOC/DOP), y=mean(DOC/DOP)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, aes(fill=Horizon, shape=Legend), show.legend = F)+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme_min+theme(legend.position = "right",
                  legend.title = element_blank())+
  #facet_wrap(~Horizon, scales="free")+
  ylab(expression(paste("DOC/DOP (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  annotate("text", label="Horizon**", 0.5, 2400, hjust="left", size=6, fontface="italic")+
  annotate("text", label="Plesne : Certovo n.s.", 0.5, 2250, hjust="left", size=6, fontface="italic")+
  annotate("text", label="Time n.s.", 0.5, 2100, hjust="left", size=6, fontface="italic")+
  ggtitle("B")+scale_shape_manual(values = c(21, 24))
anova(glm(Cmic/Nmic~Horizon+Plesne/Horizon+TIME/Plesne/Horizon, 
          data = subset(mix, Labelling=="NO"),
          family = Gamma), test="F")
mix %>% filter(Labelling=="NO") %>% 
  group_by(Plesne, Horizon, Legend) %>% summarize(y.sd=sd(Cmic/Nmic), 
                                                  y=mean(Cmic/Nmic)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, aes(fill=Horizon, shape=Legend), show.legend = F)+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme_min+theme(legend.position = "right",
                  legend.title = element_blank())+
  #facet_wrap(~Horizon, scales="free")+
  ylab(expression(paste("MBC/MBN (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  annotate("text", label="Horizon n.s.", 3.5, 16, hjust="left", size=6, fontface="italic")+
  annotate("text", label="Plesne : Certovo n.s.", 3.5, 15, hjust="left", size=6, fontface="italic")+
  annotate("text", label="Time n.s.", 3.5, 14, hjust="left", size=6, fontface="italic")+
  ggtitle("C")+scale_shape_manual(values = c(21, 24))
anova(glm(Cmic/Pmic~Horizon+Plesne/Horizon+TIME/Plesne/Horizon, 
          data = subset(mix, Labelling=="NO"),
          family = Gamma), test="F")
mix %>% filter(Labelling=="NO") %>% 
  group_by(Plesne, Horizon, Legend) %>% summarize(y.sd=sd(Cmic/Pmic), 
                                                  y=mean(Cmic/Pmic)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, aes(fill=Horizon, shape=Legend), show.legend = F)+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme_min+theme(legend.position = "right",
                  legend.title = element_blank())+
  #facet_wrap(~Horizon, scales="free")+
  ylab(expression(paste("MBC/MBP (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  annotate("text", label="Horizon***", 0.5, 80, hjust="left", size=6, fontface="italic")+
  annotate("text", label="Plesne : Certovo n.s.", 0.5, 75, hjust="left", size=6, fontface="italic")+
  annotate("text", label="Time n.s.", 0.5, 70, hjust="left", size=6, fontface="italic")+
  ggtitle("D")+scale_shape_manual(values = c(21, 24))
anova(glm(Nmic/Pmic~Horizon+Plesne/Horizon+TIME/Plesne/Horizon, 
          data = subset(mix, Labelling=="NO"),
          family = Gamma), test="F")
mix %>% filter(Labelling=="NO") %>% 
  group_by(Plesne, Horizon, Legend) %>% summarize(y.sd=sd(Nmic/Pmic), 
                                                  y=mean(Nmic/Pmic)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, aes(fill=Horizon, shape=Legend), show.legend = F)+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme_min+theme(legend.position = "right",
                  legend.title = element_blank())+
  #facet_wrap(~Horizon, scales="free")+
  ylab(expression(paste("MBN/MBP (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  annotate("text", label="Horizon***", 0.5, 8, hjust="left", size=6, fontface="italic")+
  annotate("text", label="Plesne : Certovo n.s.", 0.5, 7.5, hjust="left", size=6, fontface="italic")+
  annotate("text", label="Time n.s.", 0.5, 7, hjust="left", size=6, fontface="italic")+
  ggtitle("E")+scale_shape_manual(values = c(21, 24))


##Initial conditions
##Figure S3: Initial values of (A) pH, (B) amount of water extractable organic carbon, (C) ammonia,
##(D) nitrates, (E) organic nitrogen, (F) organic phosphorus, (G) phosphate, (H) microbial biomass 
##carbon, (I) nitrogen and (J) phosphorus in the litter (empty circles) and topsoil organic horizons
##(grey circles) of two spruce forest soils (Plešné and Čertovo catchments) mixed at five different 
##ratios (i.e., 0:1, 0.25:0.75, 0.5:0.5, 0.75:0.25, and 1:0 in respect to Plešné). Symbols represent
##mean values and error bars standard deviation of the mean (n = 4). The effect of soil horizon 
##(Horizon) and the proportion between Plešné and Čertovo in soil mixture (Plešné : Čertovo) on 
##respective soil characteristic is reported.
##Levels of significance: ***, p < 0.001; **, p < 0.01; *, p < 0.05; n.s., not significant.

anova(glm(pH~Horizon+Plesne/Horizon, mix, family = Gamma), test="F")
mix %>% group_by(Plesne, Horizon) %>% summarize(y.sd=sd(pH), y=mean(pH)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Horizon))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme_min+theme(legend.position = c(0.8, 0.15),
                  legend.key.size = unit(0.3, "in"),
                  legend.title = element_blank())+
  ylab(expression(paste("pH")))+ylim(2.3, 3.35)+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  annotate("text", label="Horizon***", 0.5, 3.35, hjust="left", size=6, fontface="italic")+
  annotate("text", label="Plesne : Certovo***", 0.5, 3.28, hjust="left", size=6, fontface="italic")+
  ggtitle("A")
mix %>% filter(Labelling=="NO" & TIME==0) %>% 
  group_by(Plesne, Horizon) %>% summarize(y.sd=sd(DOC), y=mean(DOC)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Horizon), show.legend = F)+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme_min+theme(legend.position = c(0.8, 0.15),
                  legend.title = element_blank())+
  ylim(0, 80)+
  ylab(expression(paste("DOC ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  annotate("text", label="Horizon***", 0.5, 80, hjust="left", size=6, fontface="italic")+
  annotate("text", label="Plesne : Certovo n.s.", 0.5, 73, hjust="left", size=6, fontface="italic")+
  ggtitle("B")
mix %>% filter(Labelling=="NO" & TIME==0) %>% 
  group_by(Plesne, Horizon) %>% summarize(y.sd=sd(NH4), y=mean(NH4)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Horizon), show.legend = F)+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme_min+theme(legend.position = c(0.8, 0.15),
                  legend.title = element_blank())+
  ylim(0, 8)+
  ylab(expression(paste(NH[4]," ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  annotate("text", label="Horizon***", 3.5, 7.5, hjust="left", size=6, fontface="italic")+
  annotate("text", label="Plesne : Certovo***", 3.5, 6.5, hjust="left", size=6, fontface="italic")+
  ggtitle("C")
mix %>% filter(Labelling=="NO" & TIME==0) %>% 
  group_by(Plesne, Horizon) %>% summarize(y.sd=sd(NO3), y=mean(NO3)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Horizon), show.legend = F)+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme_min+theme(legend.position = c(0.8, 0.15),
                  legend.title = element_blank())+
  ylim(0, 0.5)+
  ylab(expression(paste(NO[3]," ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  annotate("text", label="Horizon n.s.", 3.5, 0.5, hjust="left", size=6, fontface="italic")+
  annotate("text", label="Plesne : Certovo***", 3.5, 0.45, hjust="left", size=6, fontface="italic")+
  ggtitle("D")
mix %>% filter(Labelling=="NO" & TIME==0) %>% 
  group_by(Plesne, Horizon) %>% summarize(y.sd=sd(DON), y=mean(DON)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Horizon), show.legend = F)+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme_min+theme(legend.position = c(0.8, 0.15),
                  legend.title = element_blank())+
  ylim(0, 3.6)+
  ylab(expression(paste("DON ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  annotate("text", label="Horizon***", 3.5, 3.5, hjust="left", size=6, fontface="italic")+
  annotate("text", label="Plesne : Certovo n.s.", 3.5, 3.1, hjust="left", size=6, fontface="italic")+
  ggtitle("E")
mix %>% filter(Labelling=="NO" & TIME==0) %>% 
  group_by(Plesne, Horizon) %>% summarize(y.sd=sd(DOP), y=mean(DOP)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Horizon), show.legend = F)+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme_min+theme(legend.position = c(0.8, 0.15),
                  legend.title = element_blank())+
  #ylim(0, 3.6)+
  ylab(expression(paste("DOP ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  annotate("text", label="Horizon***", 0.15, 0.085, hjust="left", size=6, fontface="italic")+
  annotate("text", label="Plesne : Certovo*", 0.15, 0.075, hjust="left", size=6, fontface="italic")+
  ggtitle("F")
mix %>% filter(Labelling=="NO" & TIME==0) %>% 
  group_by(Plesne, Horizon) %>% summarize(y.sd=sd(PO4), y=mean(PO4)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Horizon), show.legend = F)+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme_min+theme(legend.position = c(0.8, 0.15),
                  legend.title = element_blank())+
  ylim(0, 0.06)+
  ylab(expression(paste(PO[4], " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  annotate("text", label="Horizon n.s.", 0.5, 0.01, hjust="left", size=6, fontface="italic")+
  annotate("text", label="Plesne : Certovo n.s.", 0.5, 0.005, hjust="left", size=6, fontface="italic")+
  ggtitle("G")
mix %>% filter(Labelling=="NO" & TIME==0) %>% 
  group_by(Plesne, Horizon) %>% summarize(y.sd=sd(Cmic), y=mean(Cmic)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Horizon), show.legend = F)+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme_min+theme(legend.position = c(0.8, 0.15),
                  legend.title = element_blank())+
  ylim(0, 900)+
  ylab(expression(paste("MBC ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  annotate("text", label="Horizon***", 0.5, 850, hjust="left", size=6, fontface="italic")+
  annotate("text", label="Plesne : Certovo***", 0.5, 740, hjust="left", size=6, fontface="italic")+
  ggtitle("H")
mix %>% filter(Labelling=="NO" & TIME==0) %>% 
  group_by(Plesne, Horizon) %>% summarize(y.sd=sd(Nmic), y=mean(Nmic)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Horizon), show.legend = F)+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme_min+theme(legend.position = c(0.8, 0.15),
                  legend.title = element_blank())+
  ylim(0, 90)+
  ylab(expression(paste("MBN ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  annotate("text", label="Horizon***", 0.5, 85, hjust="left", size=6, fontface="italic")+
  annotate("text", label="Plesne : Certovo***", 0.5, 74, hjust="left", size=6, fontface="italic")+
  ggtitle("I")
mix %>% filter(Labelling=="NO" & TIME==0) %>% 
  group_by(Plesne, Horizon) %>% summarize(y.sd=sd(Pmic), y=mean(Pmic)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Horizon), show.legend = F)+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme_min+theme(legend.position = c(0.8, 0.15),
                  legend.title = element_blank())+
  ylim(0, 30)+
  ylab(expression(paste("MBP ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  annotate("text", label="Horizon n.s.", 0.5, 26, hjust="left", size=6, fontface="italic")+
  annotate("text", label="Plesne : Certovo***", 0.5, 24, hjust="left", size=6, fontface="italic")+
  ggtitle("J")

##Microbial physiology
##Figure S5: CUE (carbon use efficiency) of microbial biomass in the litter 
##(empty symbols) and topsoil organic horizon (grey symbols) of two spruce 
##forest soils (Plešné and Čertovo catchments) that have been mixed at five 
##different ratios (i.e. 0:1, 0.25:0.75, 0.5:0.5, 0.75:0.25 and, 1:0 in 
##respect to Plešné). CUE was estimated by two independent approaches – A) 
##using 13C enriched glucose; B) using mass balance. For details see the 
##Material and methods section in the main text. Symbols represent mean 
##values and error bars standard deviation of the mean (n = 4). 
##The effect of soil horizon (Horizon) and the proportion between Plešné and 
##Čertovo in soil mixture (Plešné : Čertovo) on CUE is reported. 
##Levels of significance: ***, p < 0.001; **, p < 0.01; *, p < 0.05; n.s., 
##not significant.

anova(glm(CUE~Horizon+Plesne/Horizon, mix_diff, subset = Labelling =="NO",
          family=Gamma), test="F")
anova(glm(CUEbalance~Horizon+Plesne/Horizon, mix_diff, subset = Labelling =="NO",
          family=Gamma), test="F")

grid.arrange(mix_diff %>% filter(Labelling=="NO") %>% group_by(Plesne, Horizon) %>% 
               summarize(y.sd=sd(CUE), y=mean(CUE)) %>%
               ggplot(aes(factor(Plesne),y))+geom_point(cex=6, pch=21, aes(fill=Horizon), show.legend = F)+
               theme_min+geom_errorbar(width=0.1, aes(ymin=y-y.sd,ymax=y+y.sd))+
               scale_fill_manual(values = c("white", "grey"))+
               theme(legend.title = element_blank(),
                     legend.position = c(0.8, 0.15))+
               ylim(0.4, 0.9)+
               ylab("CUE")+xlab("Plesne : Certovo mixing ratio")+
               annotate("text", label="Horizon***", 0.5, 0.5, hjust="left", size=6, fontface="italic")+
               annotate("text", label="Plesne : Certovo**", 0.5, 0.45, hjust="left", size=6, fontface="italic") +
               ggtitle("A)"), 
             mix_diff %>% filter(Labelling=="NO") %>% group_by(Plesne, Horizon) %>% 
               summarize(y.sd=sd(CUEbalance), y=mean(CUEbalance)) %>%
               ggplot(aes(factor(Plesne),y))+geom_point(cex=6, pch=21, aes(fill=Horizon))+
               theme_min+geom_errorbar(width=0.1, aes(ymin=y-y.sd,ymax=y+y.sd))+
               scale_fill_manual(values = c("white", "grey"))+
               theme(legend.title = element_blank(),
                     legend.position = c(0.8, 0.2))+
               annotate("text", label="Horizon n.s.", 0.5, 0.5, hjust="left", size=6, fontface="italic")+
               annotate("text", label="Plesne : Certovo n.s.", 0.5, 0.45, hjust="left", size=6, fontface="italic") +
               ylim(0.4, 0.9)+
               ylab("CUE")+xlab("Plesne : Certovo mixing ratio")+
               ggtitle("B)"), nrow=1)



##Figure S6: Cumulative CO2 loss from the litter (empty symbols) and topsoil organic horizon 
##(grey symbols) of two spruce forest soils (Plešné and Čertovo catchments) that have been mixed 
##at five different ratios (i.e., 0:1, 0.25:0.75, 0.5:0.5, 0.75:0.25, and 1:0 in respect to Plešné). 
##Symbols represent mean values and error bars standard deviation of the mean (n = 4). The effect of 
##soil horizon (Horizon) and the proportion between Plešné and Čertovo in soil mixture 
##(Plešné : Čertovo) on cumulative CO2 loss is reported.
##Levels of significance: ***, p < 0.001; **, p < 0.01; *, p < 0.05; n.s., not significant.

anova(glm(CO2~Horizon+Plesne/Horizon, mix_diff, subset = Labelling=="NO",
          family=Gamma), test="F")
mix_diff %>% filter(Labelling=="NO") %>% group_by(Plesne, Horizon) %>% 
  summarize(y.sd=sd(CO2, na.rm = T), y=mean(CO2, na.rm = T)) %>%
  ggplot(aes(factor(Plesne),y))+geom_point(cex=6, pch=21, aes(fill=Horizon), show.legend = T)+
  theme_min+geom_errorbar(width=0.1, aes(ymin=y-y.sd,ymax=y+y.sd))+
  scale_fill_manual(values = c("white", "grey"))+
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.15))+
  ylim(0, 25)+
  ylab(expression(paste("Cumulative ",CO[2]," loss  (", mu, "mol ", g^{-1}, ")")))+
  xlab("Plesne : Certovo mixing ratio")+
  annotate("text", label="Horizon***", 0.5, 24, hjust="left", size=6, fontface="italic")+
  annotate("text", label="Plesne : Certovo**", 0.5, 22, hjust="left", size=6, fontface="italic")

##Figure S7: Amount of lost carbon in form of CO2 (black circles), microbial biomass carbon 
##(MBC, grey circles) and water extractable organic carbon (DOC, empty circles) in the litter 
##and topsoil organic horizon of two spruce forest soils (Plešné and Čertovo catchments) that 
##have been mixed at five different ratios (i.e., 0:1, 0.25:0.75, 0.5:0.5, 0.75:0.25, and 1:0 in 
##respect to Plešné). Symbols represent mean values and error bars standard deviation of the mean 
##(n = 4). Note that graphs have different y-axis scales.

mix_diff$dDOC<-mix[c(81:160), c("DOC")]-mix[c(1:80), c("DOC")]
cbalance<-mix_diff[mix_diff$Labelling=="NO", 
                   c("Plesne", "Horizon", "CO2", "dCmic", "dDOC")]
cbalance$dCmic<-(-cbalance$dCmic)
cbalance$dDOC<-(-cbalance$dDOC)
Cbalance<-melt(cbalance, id.vars = c("Plesne", "Horizon"))

Cbalance %>% group_by(Plesne, Horizon, variable) %>% 
  summarize(y.sd=sd(value), y=mean(value)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=variable))+
  facet_wrap(~Horizon, scales="free")+theme_min+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.65, 0.8),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank(),
        legend.text.align = 0)+
  ylab(expression(paste("C losses (", mu, "mol ", g^{-1}, ")")))+
  xlab("Plesne : Certovo mixing ratio")+
  theme(legend.title = element_blank())+
  geom_hline(yintercept = 0)+
  scale_fill_manual(values = c("black", "grey", "white"),
                    name = '',
                    labels = expression(C-CO[2], MBC, DOC))

##Figure S8: Measured versus predicted concentration of soil microbial 
##biomass carbon (MBC), dissolved organic carbon (DOC), and headspace 
##carbon dioxide (CO2) at the beginning and at the end of 48-hour incubation. 
##Different colors of symbols represent different Plešné to Čertovo mixing 
##ratios. Symbols represent mean values and error bars standard error of the
##mean (n = 4). Solid black lines represent 1:1 relationship, thus the 
##absolute correspondence between the predictions and observations. 
##The change in carbon pools over time was simulated by microbially-explicit
##model depicted on Fig. S4 and described in Material and methods section 
##in the main text. Note that graphs have different y-axis scales.

Preds<-Cuest_out[[1]]$fit

for(i in 2:40){
  Preds<-rbind(Preds, Cuest_out[[i]]$fit)
}

Preds$outliers<-"NO"
Preds[(Preds$variable=="CO2" & Preds$value>100), "outliers"]<-"YES"
Preds$Legend<-rep(c(rep("PL:CT=1:0", times=6),
                rep("PL:CT=0.75:0.25", times=6),
                rep("PL:CT=0.5:0.5", times=6),
                rep("PL:CT=0.25:0.75", times=6),
                rep("PL:CT=0:1", times=6)), 8)
Preds$variable2<-Preds$variable
levels(Preds$variable2)<-c("MBC", "DOC", "CO[2]")

Preds %>% filter(outliers=="NO") %>% group_by(time, variable2, Legend, Horizon) %>%
  summarise(y=mean(value, na.rm=T), x=mean(obs, na.rm=T),
            y.se=sd(value, na.rm=T)/sqrt(4), x.se=sd(obs, na.rm=T)/sqrt(4)) %>%
  ggplot(aes(x, y)) + geom_point(cex=6, pch=21, aes(fill=Legend)) + 
  geom_errorbar(aes(ymin=y-y.se, ymax=y+y.se))+
  geom_errorbarh(aes(xmin=x-x.se, xmax=x+x.se))+
  facet_wrap(.~variable2, scales="free", labeller = label_parsed) + 
  geom_abline(intercept = 0, slope=1, lwd=1.5) + 
  theme_min +
  xlab(expression(paste("Predicted value (", mu, "mol(C)", g^{-1}, ")"))) +
  ylab(expression(paste("Measured value (", mu, "mol(C)", g^{-1}, ")")))
