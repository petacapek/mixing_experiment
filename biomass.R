############################################Biomass data##########################################
##Cmic
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(Cmic), y=mean(Cmic)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.88))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste("Microbial biomass carbon ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##Nmic
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(Nmic), y=mean(Nmic)) %>%
ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.88))+
  ylab(expression(paste("Microbial biomass nitrogen ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##Pmic
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(Pmic), y=mean(Pmic)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.17, 0.88))+
  ylab(expression(paste("Microbial biomass phosphorus ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##C/N mic
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(Cmic/Nmic), y=mean(Cmic/Nmic)) %>%
ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.33, 0.9),
                                                              legend.key.size = unit(0.3, "in"))+
  ylab(expression(paste("Microbial C/N (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##C/P mic
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(Cmic/Pmic), y=mean(Cmic/Pmic)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont)+theme_min+theme(legend.position = c(0.17, 0.88),
                                                              legend.key.size = unit(0.3, "in"))+
  ylab(expression(paste("Microbial C/P (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##N/P mic
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(Nmic/Pmic), y=mean(Nmic/Pmic)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.17, 0.88),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste("Microbial N/P (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

####zmeny
mb.zmeny<-mich2[1:80, c("Plesne", "Certovo", "horizont", "znaceni")]
mb.zmeny$dCmic<-mich2[81:160, "Cmic"]-mich2[1:80, "Cmic"]
mb.zmeny$dNmic<-mich2[81:160, "Nmic"]-mich2[1:80, "Nmic"]
mb.zmeny$dPmic<-mich2[81:160, "Pmic"]-mich2[1:80, "Pmic"]

##Cmic
mb.zmeny %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(dCmic), y=mean(dCmic)) %>%
ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(Delta, "Microbial biomass carbon ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.position = c(0.15, 0.2),
                                            legend.key.size = unit(0.3, "in"),
                                            legend.title = element_blank())

##Nmic
mb.zmeny %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(dNmic), y=mean(dNmic)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(Delta, "Microbial biomass nitrogen ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.position = c(0.15, 0.2),
                                            legend.key.size = unit(0.3, "in"),
                                            legend.title = element_blank())
##Pmic
mb.zmeny %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(dPmic), y=mean(dPmic)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(Delta, "Microbial biomass phosphorus ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.position = c(0.15, 0.2),
                                            legend.key.size = unit(0.3, "in"),
                                            legend.title = element_blank())
##C/N mic
mb.zmeny[-59, ] %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(dCmic/dNmic), y=mean(dCmic/dNmic)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(Delta, "C/N biomass ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.position = c(0.7, 0.2),
                                            legend.key.size = unit(0.3, "in"),
                                            legend.title = element_blank())

##C/P mic
mb.zmeny[-c(73, 31, 60, 78), ] %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(dCmic/dPmic), y=mean(dCmic/dPmic)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(Delta, "C/P biomass ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.position = c(0.2, 0.8),
                                            legend.key.size = unit(0.3, "in"),
                                            legend.title = element_blank())

##N/P mic
mb.zmeny[-c(73, 60, 78), ] %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(dNmic/dPmic), y=mean(dNmic/dPmic)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(Delta, "N/P biomass ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.position = c(0.2, 0.8),
                                            legend.key.size = unit(0.3, "in"),
                                            legend.title = element_blank())

###########################################################################################

