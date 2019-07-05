################################################################################################
#Biological and stoichiometric controls over nutrients mineralization and immobilization in two#
###################################spruce forest soils##########################################
################################################################################################
#Libraries
library(dplyr)
library(reshape)
library(ggplot2)
library(gridExtra)
library(vegan)
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

##Communities
###Bacteria
bac<-t(read.csv("./manuscript_data_and_code/mixing_bacteria.csv", header = F)[-1, -c(1, 42:49)])
rownames(bac)<-read.csv("./manuscript_data_and_code/mixing_bacteria.csv", header = F)[1,-c(1, 42:49)]
#Fungi
fungi<-t(read.csv("./manuscript_data_and_code/mixing_fungi.csv", header = F)[-1, -c(1, 42:49)])
rownames(fungi)<-read.csv("./manuscript_data_and_code/mixing_fungi.csv", header = F)[1,-c(1, 42:49)]

###Labels
lb<-read.csv("./manuscript_data_and_code/community_labels.csv")
lb$SampleID<-as.factor(lb$SampleID)
lb$SampleID<-factor(lb$SampleID, levels = rownames(bac))
lb<-lb[order(lb$SampleID), ]

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

##MICROBIAL COMMUNITY
###Bacteria
bnorm<-decostand(bac, method=c("normalize"))
bdist<-vegdist(bnorm, method = "jaccard")

bad<-adonis(bdist~horizon+PL/horizon, lb)
bad

###Fungi
fnorm<-decostand(fungi, method=c("normalize"))
fdist<-vegdist(fnorm, method = "jaccard")

fad<-adonis(fdist~horizon+PL/horizon, lb)
fad

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
##Figure 2: Initial values of microbial biomass carbon to nitrogen (MBC:MBN; A) and phosphorus 
##(MBC:MBP; C) ratio, and water extractable organic carbon to nitrogen (DOC:DON; B) and phosphorus (DOC:DOP; D)  
##ratio in two spruce forest soils (Plešné and Čertovo) mixed at three different ratios (i.e. ¼, ½, and ¾). 
##Grey circles denote litter layer and empty circles topsoil organic layer. 
##Symbols represent mean values (n = 4) and error bars standard deviation of the mean. 
##Notice that all plots have different y-axes scales.

grid.arrange(
  mix %>% filter(Labelling=="NO" & TIME==0) %>% 
    group_by(Plesne, Horizon) %>% summarize(y.sd=sd(Cmic/Nmic), y=mean(Cmic/Nmic)) %>%
    ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Horizon), 
                                              position = position_dodge(width = 0.6))+
    geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd, colour=Horizon), width=0.1, lwd=0.5, 
                  position = position_dodge(width = 0.6))+
    theme_min+theme(legend.position = c(0.6, 0.8), legend.title = element_blank())+
    ylab(expression(paste("MBC:MBN (mol:mol)" )))+
    xlab("Plesne : Certovo mixing ratio")+scale_color_manual(values = c("black", "black"))+
    scale_fill_manual(values = c("grey", "white"))+
    ggtitle("A)")+scale_y_continuous(limits = c(10, 20)),
  mix %>% filter(Labelling=="NO" & TIME==0) %>% 
    group_by(Plesne, Horizon) %>% summarize(y.sd=sd(DOC/DON), y=mean(DOC/DON)) %>%
    ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Horizon), 
                                              position = position_dodge(width = 0.6), show.legend = F)+
    geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd, colour=Horizon), width=0.1, lwd=0.5, 
                  position = position_dodge(width = 0.6), show.legend = F)+
    theme_min+theme(legend.position = c(0.8, 0.9), legend.title = element_blank())+
    ylab(expression(paste("DOC:DON (mol:mol)" )))+
    xlab("Plesne : Certovo mixing ratio")+scale_color_manual(values = c("black", "black"))+
    scale_fill_manual(values = c("grey", "white"))+
    ggtitle("B)")+scale_y_continuous(limits = c(10, 30)),
  mix %>% filter(Labelling=="NO" & TIME==0) %>% 
    group_by(Plesne, Horizon) %>% summarize(y.sd=sd(Cmic/Pmic), y=mean(Cmic/Pmic)) %>%
    ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Horizon), 
                                              position = position_dodge(width = 0.6), show.legend = F)+
    geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd, colour=Horizon), width=0.1, lwd=0.5, 
                  position = position_dodge(width = 0.6), show.legend = F)+
    theme_min+theme(legend.position = c(0.8, 0.9), legend.title = element_blank())+
    ylab(expression(paste("MBC:MBP (mol:mol)" )))+
    xlab("Plesne : Certovo mixing ratio")+scale_color_manual(values = c("black", "black"))+
    scale_fill_manual(values = c("grey", "white"))+
    ggtitle("C)")+scale_y_continuous(limits = c(0, 80)),
  mix %>% filter(Labelling=="NO" & TIME==0) %>% 
    group_by(Plesne, Horizon) %>% summarize(y.sd=sd(DOC/DOP), y=mean(DOC/DOP)) %>%
    ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Horizon), 
                                              position = position_dodge(width = 0.6), show.legend = F)+
    geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd, colour=Horizon), width=0.1, lwd=0.5, 
                  position = position_dodge(width = 0.6), show.legend = F)+
    theme_min+theme(legend.position = c(0.8, 0.9), legend.title = element_blank())+
    ylab(expression(paste("DOC:DOP (mol:mol)" )))+
    xlab("Plesne : Certovo mixing ratio")+scale_color_manual(values = c("black", "black"))+
    scale_fill_manual(values = c("grey", "white"))+
    ggtitle("D)")+scale_y_continuous(limits = c(0, 2000)), ncol=2)

##Figure 3: Net change of soil microbial biomass carbon (MBC; A), nitrogen (MBN; B) 
##and phosphorus (MBP; C) 
##after 48 hours incubation of two spruce forest soils (Plešné and Čertovo) mixed at three 
##different ratios (i.e. ¼, ½, and ¾). Grey circles denote litter layer and empty circles 
##topsoil organic layer. Symbols represent mean values (n = 4) and error bars standard 
##deviation of the mean. Notice that both plots have different y-axes scales.
mix_diff<-mix[c(1:80), c("Plesne", "Certovo", "Horizon", "Labelling")]
mix_diff$dCmic<-mix[c(81:160), c("Cmic")]-mix[c(1:80), c("Cmic")]
mix_diff$dNmic<-mix[c(81:160), c("Nmic")]-mix[c(1:80), c("Nmic")]
mix_diff$dPmic<-mix[c(81:160), c("Pmic")]-mix[c(1:80), c("Pmic")]

grid.arrange(
  mix_diff %>% filter(Labelling=="NO") %>% group_by(Plesne, Horizon) %>% 
    summarize(y.sd=sd(dCmic), y=mean(dCmic)) %>%
    ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=Horizon),
                                              position = position_dodge(width = 0.6))+
    theme_min+
    geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd, colour=Horizon), width=0.1, lwd=0.5,
                  position = position_dodge(width = 0.6))+
    ylab(expression(paste(Delta~MBC, " (", mu, "mol ",g^{-1}, ")" )))+
    xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("grey", "white"))+
    scale_color_manual(values = c("black", "black"))+scale_y_continuous(breaks = c(-150, -100, -50, 0, 50, 100))+
    geom_hline(yintercept = 0, lwd=1)+theme(legend.position = c(0.65, 0.8),
                                            legend.key.size = unit(0.3, "in"),
                                            legend.title = element_blank())+
    ggtitle("A)"),
  mix_diff %>% filter(Labelling=="NO") %>% group_by(Plesne, Horizon) %>% 
    summarize(y.sd=sd(dNmic), y=mean(dNmic)) %>%
    ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=Horizon), show.legend = F,
                                              position = position_dodge(width = 0.6))+
    theme_min+theme(legend.position = c(0.2, 0.2))+
    geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd, colour=Horizon), width=0.1, lwd=0.5,
                  position = position_dodge(width = 0.6), show.legend = F)+
    ylab(expression(paste(Delta~MBN, " (", mu, "mol ",g^{-1}, ")" )))+
    xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("grey", "white"))+
    scale_color_manual(values = c("black", "black"))+#scale_y_continuous(breaks = c(-150, -100, -50, 0, 50, 100))+
    geom_hline(yintercept = 0, lwd=1)+theme(legend.position = c(0.8, 0.8),
                                            legend.key.size = unit(0.3, "in"),
                                            legend.title = element_blank())+
    ggtitle("B)"),
  mix_diff %>% filter(Labelling=="NO") %>% group_by(Plesne, Horizon) %>% 
    summarize(y.sd=sd(dPmic), y=mean(dPmic)) %>%
    ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=Horizon), show.legend = F,
                                              position = position_dodge(width = 0.6))+
    theme_min+theme(legend.position = c(0.2, 0.2))+
    geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd, colour=Horizon), width=0.1, lwd=0.5,
                  position = position_dodge(width = 0.6), show.legend = F)+
    ylab(expression(paste(Delta~MBP, " (", mu, "mol ",g^{-1}, ")" )))+
    xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("grey", "white"))+
    scale_color_manual(values = c("black", "black"))+#scale_y_continuous(breaks = c(-150, -100, -50, 0, 50, 100))+
    geom_hline(yintercept = 0, lwd=1)+theme(legend.position = c(0.8, 0.8),
                                            legend.key.size = unit(0.3, "in"),
                                            legend.title = element_blank())+
    ggtitle("C)"), ncol=3)

##Figure 4: Isotopic signal of respired CO2 (black symbols), K2SO4 extractable organic carbon 
##(K2SO4 - EC, grey symbols) and microbial biomass carbon (MBC, empty symbols) in litter and 
##topsoil organic layer of two spruce forest soils (Plešné and Čertovo) mixed at three 
##different ratios (i.e. ¼, ½, and ¾). Symbols represent mean values (n = 4) and error bars 
##standard deviation of the mean. 
izo<-mix[c(81:160), c("Plesne", "Certovo", "Horizon", "Labelling", "CO2atm",
                   "DOCatm", "Cmicatm")]
Izo<-melt(izo, id.vars = c("Plesne", "Certovo", "Horizon", "Labelling"))
###Calculate delta values
Izo$delta<-with(Izo, (value/(1-value)/0.011237-1)*1000)

Izo %>% filter(Labelling=="NO") %>%
  group_by(Plesne, Horizon, variable) %>% 
  summarize(y.sd=sd(delta, na.rm = T), y=mean(delta, na.rm = T)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=variable))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(.~Horizon)+theme_min+
  theme(legend.position = c(0.1, 0.85),
        legend.title = element_blank(),
        legend.text.align = 0)+
  ylab(expression(paste(delta^{13},C)))+
  xlab("Plesne : Certovo mixing ratio")+
  scale_fill_manual(values = c("black", "grey", "white"),
                    name = '',
                    labels = expression(C-CO[2], K[2]~SO[4]-EC, MBC))

##Figure 6: Net change of water extractable mineral nitrogen (ΔMN) in litter (grey symbols) 
##and topsoil organic layer (empty circles) of two spruce forest soils (Plešné and Čertovo) 
##mixed at three different ratios (i.e. ¼, ½, and ¾). Blue circles represent ΔMN calculated 
##using eq. 1 with (A) and without (B) fMBC (for details see Materials and methods section). 
##Symbols represent mean values (n = 4) and error bars standard deviation of the mean.

###Assuming no contribution of decaying microbial biomass
mix_diff$CUE<-mix[c(1:80), c("CUE")]
mix_diff$CO2<-mix[c(1:80), c("CCO2")]
mix_diff$DOC<-mix[c(1:80), c("DOC")]
mix_diff$DON<-mix[c(1:80), c("DON")]
mix_diff$Cmic<-mix[c(1:80), c("Cmic")]
mix_diff$Nmic<-mix[c(1:80), c("Nmic")]
mix_diff$dNH4<-mix[c(81:160), c("NH4")]-mix[c(1:80), c("NH4")]
mix_diff$dNO3<-mix[c(81:160), c("NO3")]-mix[c(1:80), c("NO3")]
####Amount of organic carbon utilized
mix_diff$U<-with(mix_diff, CO2+CUE*CO2/(1-CUE))
####Expected amount of mineral nitrogen exchanged
mix_diff$dNm_preda<-with(mix_diff, U*(DON/DOC-Nmic*CUE/Cmic))

###Assuming contribution of decaying microbial biomass
####Calculate relative contribution of decaying microbial biomass to food source
####of the growing microbial biomass
mix_diff$deltaCO2<-Izo[Izo$variable=="CO2atm", "delta"]
mix_diff$deltaCmic<-Izo[Izo$variable=="Cmicatm", "delta"]
mix_diff$deltaBase<-mean(Izo[(Izo$variable=="CO2atm" & Izo$Horizon=="Litter" &
                                Izo$Plesne==0 & Izo$Labelling=="NO"), "delta"], na.rm=T)
mix_diff$Cmic.prop<-with(mix_diff, (deltaCO2-deltaBase)/(deltaCmic-deltaBase))
####Organic topsoil Cmic.prop is 0
mix_diff[(mix_diff$Horizon=="Organic topsoil"), "Cmic.prop"]<-0
####All Cmic.prop higher than 1 and lower than 0 are set to be 1 and 0 respectively
mix_diff[(mix_diff$Cmic.prop>1), "Cmic.prop"]<-1
mix_diff[(mix_diff$Cmic.prop<0), "Cmic.prop"]<-0


mix_diff$dNm_predb<-with(mix_diff, U*(((1-Cmic.prop)*(DON/DOC)+Cmic.prop*(Nmic/Cmic))-Nmic*CUE/Cmic))

grid.arrange(
  mix_diff %>% filter(Labelling=="NO") %>% group_by(Plesne, Horizon) %>% 
    summarize(y.sd=sd(dNm_preda), y=mean(dNm_preda),
              y2.sd=sd(dNH4+dNO3), y2=mean(dNH4+dNO3)) %>%
    ggplot(aes(Plesne, y))+geom_point(cex=6, pch=21, fill="dodgerblue3")+
    facet_grid(.~Horizon, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
    geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
    theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
          legend.title = element_blank(),
          axis.title.x = element_blank())+
    ylab(expression(paste(italic(Delta~M[N]), " (", mu, "mol ",g^{-1}, ")" )))+
    xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("grey", "white"))+
    geom_hline(yintercept = 0, lwd=1)+theme(legend.title = element_blank())+
    geom_point(cex=6, pch=21, aes(Plesne, y2, fill=Horizon), show.legend = F)+
    geom_errorbar(aes(ymin=y2-y2.sd, ymax=y2+y2.sd), width=0.1, lwd=0.5)+ylim(-6.5, 1.5)+
    ggtitle("A)"),
  mix_diff %>% filter(Labelling=="NO") %>% group_by(Plesne, Horizon) %>% 
    summarize(y.sd=sd(dNm_predb), y=mean(dNm_predb),
              y2.sd=sd(dNH4+dNO3), y2=mean(dNH4+dNO3)) %>%
    ggplot(aes(Plesne, y))+geom_point(cex=6, pch=21, fill="dodgerblue3")+
    facet_grid(.~Horizon, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
    geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
    theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
          legend.title = element_blank())+
    ylab(expression(paste(italic(Delta~M[N]), " (", mu, "mol ",g^{-1}, ")" )))+
    xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("grey", "white"))+
    geom_hline(yintercept = 0, lwd=1)+theme(legend.title = element_blank())+
    geom_point(cex=6, pch=21, aes(Plesne, y2, fill=Horizon), show.legend = F)+
    geom_errorbar(aes(ymin=y2-y2.sd, ymax=y2+y2.sd), width=0.1, lwd=0.5)+ylim(-6.5, 1.5)+
    ggtitle("B)"), nrow=2)

##Figure 7: Net change of water extractable mineral phosphorus in litter (grey symbols) and 
##topsoil organic layer (empty circles) of two spruce forest soils (Plešné and Čertovo) mixed 
##at three different ratios (i.e. ¼, ½, and ¾). Blue circles represent expected change of water 
##extractable mineral phosphorus assuming that (A): the entire source of organic carbon and 
##phosphorus for microbial growth is soil derived, (B): decaying microbial biomass partly 
##contribute to the source of organic carbon and phosphorus (for details see Materials and 
##methods section) or (C): microbial biomass growth is independent of soil phosphorus and the 
##sorption characteristics of the soil affects mineral phosphorus concentration. Symbols 
##represent mean values (n = 4) and error bars standard deviation of the mean.

###Assuming no contribution of decaying microbial biomass
mix_diff$DOP<-mix[c(1:80), c("DOP")]
mix_diff$Pmic<-mix[c(1:80), c("Pmic")]
mix_diff$PO4<-mix[c(1:80), c("PO4")]
mix_diff$dPO4<-mix[c(81:160), c("PO4")]-mix[c(1:80), c("PO4")]

mix_diff$dPm_preda<-with(mix_diff, U*(DOP/DOC-Pmic*CUE/Cmic))

###Assuming contribution of decaying microbial biomass
mix_diff$dPm_predb<-with(mix_diff, U*(((1-Cmic.prop)*(DOP/DOC)+Cmic.prop*(Pmic/Cmic))-Pmic*CUE/Cmic))

###Assuming independence of microbial growth on soil derived phosphorus
###water extractable phosphate is chemicaly sorbed instead 
mix_diff$Sorption<-mix[c(1:80), c("Sorption")]
####All values greater than 1 is set to 1
mix_diff[mix_diff$Sorption>1, "Sorption"]<-1

####Calculating the amount of chemicaly sorbed phosphate
mix_diff$PO4sorbed<-with(mix_diff, PO4*(1-exp(-(1-Sorption)/0.75*48)))

grid.arrange(
  mix_diff %>% filter(Labelling=="NO") %>% group_by(Plesne, Horizon) %>% 
    summarize(y.sd=sd(dPm_preda), y=mean(dPm_preda),
              y2.sd=sd(dPO4), y2=mean(dPO4)) %>%
    ggplot(aes(Plesne, y))+geom_point(cex=6, pch=21, fill="dodgerblue3")+
    facet_grid(.~Horizon, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
    geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
    theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
          legend.title = element_blank(),
          axis.title.x = element_blank())+
    ylab(expression(paste(italic(Delta~M[P]), " (", mu, "mol ",g^{-1}, ")" )))+
    xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("grey", "white"))+
    geom_hline(yintercept = 0, lwd=1)+theme(legend.title = element_blank())+
    geom_point(cex=6, pch=21, aes(Plesne, y2, fill=Horizon), show.legend = F)+
    geom_errorbar(aes(ymin=y2-y2.sd, ymax=y2+y2.sd), width=0.1, lwd=0.5)+
    ggtitle("A)"),
  mix_diff %>% filter(Labelling=="NO") %>% group_by(Plesne, Horizon) %>% 
    summarize(y.sd=sd(dPm_predb), y=mean(dPm_predb),
              y2.sd=sd(dPO4), y2=mean(dPO4)) %>%
    ggplot(aes(Plesne, y))+geom_point(cex=6, pch=21, fill="dodgerblue3")+
    facet_grid(.~Horizon, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
    geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
    theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
          legend.title = element_blank(),
          axis.title.x = element_blank())+
    ylab(expression(paste(italic(Delta~M[P]), " (", mu, "mol ",g^{-1}, ")" )))+
    xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("grey", "white"))+
    geom_hline(yintercept = 0, lwd=1)+theme(legend.title = element_blank())+
    geom_point(cex=6, pch=21, aes(Plesne, y2, fill=Horizon), show.legend = F)+
    geom_errorbar(aes(ymin=y2-y2.sd, ymax=y2+y2.sd), width=0.1, lwd=0.5)+
    ggtitle("B)"),
  mix_diff %>% filter(Labelling=="NO") %>% group_by(Plesne, Horizon) %>% 
    summarize(y.sd=sd(-PO4sorbed), y=mean(-PO4sorbed),
              y2.sd=sd(dPO4), y2=mean(dPO4)) %>%
    ggplot(aes(Plesne, y))+geom_point(cex=6, pch=21, fill="dodgerblue3")+
    facet_grid(.~Horizon, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
    geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
    theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"))+
    ylab(expression(paste(italic(Delta~M[P]), " (", mu, "mol ",g^{-1}, ")" )))+
    xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("grey", "white"))+
    geom_hline(yintercept = 0, lwd=1)+theme(legend.title = element_blank())+
    geom_point(cex=6, pch=21, aes(Plesne, y2, fill=Horizon), show.legend = F)+
    geom_errorbar(aes(ymin=y2-y2.sd, ymax=y2+y2.sd), width=0.1, lwd=0.5)+
    ggtitle("C)"), ncol=1)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Supplementary information
##Stoichiometry
###Figure S1: Molar C/N (A) and C/P ratios of plant litter inputs into low elevation localities of Plesne 
###(red circles) and Certovo (black circles) catchments with tree vegetation mostly unaffected by bark beetle 
###infestation. Symbols represent medians and error bars lower and upper quartiles. The data for shown in this
###figure can be downloaded at https://link.springer.com/article/10.1007/s10533-018-0470-x#SupplementaryMaterial. 

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

###Figure S2: Molar DOC/DON (A, water extractable organic carbon to nitrogen), 
###DOC/DOP (B, water extractable organic carbon to phosphorus), 
###MBC/MBN (C, microbial biomass carbon to nitrogen), 
###MBC/MBP (D, microbial biomass carbon to phosphorus) and 
###MBN/MBN (E, microbial biomass nitrogen to phosphorus) ratio in litter (empty symbols) and 
###topsoil organic layer (grey symbols) of two spruce forest soils (Plesne and Certovo) 
###mixed at three different ratios (i.e. ¼, ½, and ¾). Triangles denote the values of 
###respective ratio at the beginning of the incubation and circles the values at the end of 
###the incubation. The effect of soil layer (Horizon), the proportion between Plesne and 
###Certovo in soil mixture (Plesne : Certovo) and time on respective ratios is reported. 
###Levels of significance: ***, p < 0.001; **, p < 0.01; *, p < 0.05; n.s., not significant.
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
###Figure S3: Initial pH (A), amount of water extractable organic carbon (B), ammonia (C), 
###nitrates (D), organic nitrogen (E), organic phosphorus (F), phosphate (G), microbial biomass 
###carbon (H), nitrogen (I) and phosphorus (J) in litter (empty circles) and topsoil organic 
###layer (grey circles) of two spruce forest soils (Plesne and Certovo) mixed 
###at three different ratios (i.e. ¼, ½, and ¾). The effect of soil layer (Horizon) and the 
###proportion between Plesne and Certovo in soil mixture (Plesne : Certovo) on respective soil 
###characteristic is reported. 
###Levels of significance: ***, p < 0.001; **, p < 0.01; *, p < 0.05; n.s., not significant.

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
###Figure S4: CUE (A, carbon use efficiency) and microbial biomass carbon turnover rate (B) 
###in litter (empty symbols) and topsoil organic layer (grey symbols) of two spruce forest 
###soils (Plesne and Certovo) mixed at three different ratios (i.e. ¼, ½, and ¾). 
###The effect of soil layer (Horizon) and the proportion between Plesne and Certovo in soil 
###mixture (Plesne : Certovo) on respective microbial characteristics is reported. 
###Levels of significance: ***, p < 0.001; **, p < 0.01; *, p < 0.05; n.s., not significant.

anova(glm(CUE~Horizon+Plesne/Horizon, mix_diff, subset = Labelling =="NO",
          family=Gamma), test="F")
mix_diff %>% filter(Labelling=="NO") %>% group_by(Plesne, Horizon) %>% 
  summarize(y.sd=sd(CUE), y=mean(CUE)) %>%
  ggplot(aes(factor(Plesne),y))+geom_point(cex=6, pch=21, aes(fill=Horizon))+
  theme_min+geom_errorbar(width=0.1, aes(ymin=y-y.sd,ymax=y+y.sd))+
  scale_fill_manual(values = c("white", "grey"))+
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.15))+
  ylim(0.4, 0.9)+
  ylab("CUE")+xlab("Plesne : Certovo mixing ratio")+
  annotate("text", label="Horizon***", 0.5, 0.5, hjust="left", size=6, fontface="italic")+
  annotate("text", label="Plesne : Certovo**", 0.5, 0.45, hjust="left", size=6, fontface="italic")+
  ggtitle("A")

####Calculate the turnover rate
mix_diff$growth<-with(mix_diff, CO2*CUE/(1-CUE))
mix_diff$turnover<-(-log(mix_diff$dCmic+mix_diff$growth)+log(mix_diff$Cmic))/48

anova(glm(turnover~Horizon+Plesne/Horizon, mix_diff, subset = Labelling =="NO",
          family=Gamma), test="F")
mix_diff %>% filter(Labelling=="NO") %>% group_by(Plesne, Horizon) %>% 
  summarize(y.sd=sd(turnover, na.rm = T), y=mean(turnover, na.rm = T)) %>%
  ggplot(aes(factor(Plesne),y))+geom_point(cex=6, pch=21, aes(fill=Horizon), show.legend = F)+
  theme_min+geom_errorbar(width=0.1, aes(ymin=y-y.sd,ymax=y+y.sd))+
  scale_fill_manual(values = c("white", "grey"))+
  theme(legend.title = element_blank(),
        legend.position = c(0.9, 0.15))+
  #ylim(0, 0.09)+
  ylab(expression(paste("turnover rate (", h^{-1}, ")")))+
  xlab("Plesne : Certovo mixing ratio")+
  annotate("text", label="Horizon n.s.", 0.5, 0.005, hjust="left", size=6, fontface="italic")+
  annotate("text", label="Plesne : Certovo n.s", 0.5, 0, hjust="left", size=6, fontface="italic")+
  ggtitle("B")


###Figure S5: Cumulative CO2 loss from litter (empty symbols) and topsoil organic layer 
###(grey symbols) of two spruce forest soils (Plesne and Certovo) that have been mixed at 
###three different ratios (i.e. ¼, ½, and ¾). The effect of soil layer (Horizon) and the 
###proportion between Plesne and Certovo in soil mixture (Plesne : Certovo) on cumulative CO2 
###loss is reported. 
###Levels of significance: ***, p < 0.001; **, p < 0.01; *, p < 0.05; n.s., not significant.

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

###Figure S6: Amount of lost carbon in form of CO2 (black circles), microbial biomass carbon 
###(MBC, grey circles) and water extractable organic carbon (DOC, empty circles) in litter and
###topsoil organic layer of two spruce forest soils (Plesne and Certovo) mixed at three 
###different ratios (i.e. ¼, ½, and ¾). Symbols represent mean values (n = 4) and error bars 
###standard deviation of the mean. Note that both graphs have different y-axis scales.

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

###Figure S7: Expected growth rate of microbial community in litter (grey circles) and organic
###topsoil (empty circles) from Plešné and Čertovo catchments without the external source of 
###phosphorus (i.e. phosphors demand is covered exclusively from internal resources). Symbols 
###represent mean values (n = 4) and error bars standard deviation of the mean. Light grey 
###points represent original dataset of topsoil horizons published in Čapek et al. (2016). 
###Solid black line represents linear relationship between microbial biomass carbon to phosphorus 
###ratio and growth rate in absence of external source of phosphorus. 
###Horizontal grey line denotes growth rate zero.

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
ggplot(data=gr_pred, aes(CP, growth))+geom_point(cex=6, pch=21, aes(fill=Horizon))+
  theme_min+
  xlab("MBC/MBP (mol/mol)")+
  ylab(expression(paste(mu[0])))+
  geom_hline(yintercept = 0, color='grey')+
  geom_errorbarh(aes(xmin=CP-CP.sd, xmax=CP+CP.sd))+
  geom_errorbar(aes(ymin=growth-growth.sd, ymax=growth+growth.sd))+
  geom_point(data=gr, cex=6, color='grey', alpha=0.5, aes(CP, u0))+
  stat_smooth(data=gr, aes(CP, u0), method = lm, se=F, color='black', lwd=0.5, fullrange = T)+
  scale_fill_manual(values = c("grey", "white"))+
  theme(legend.position = c(0.8, 0.85),
        legend.title = element_blank())+
  scale_x_continuous(breaks = c(20, 40, 60, 80, 100))+
  scale_y_continuous(breaks=c(-1, -0.5, 0, 0.5, 1, 1.5, 2))
