##########################################enzymy#################################################
##beta-glucosidase
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(BG), y=mean(BG)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.67, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(beta,"-glucosidase ( ", h^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##cellobiosidase
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(CELL), y=mean(CELL)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.67, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste("cellobiosidase ( ", h^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##alanin-aminopeptidase
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(ALA), y=mean(ALA)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste("alanin aminopeptidase ( ", h^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##chitinase
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(CHIT), y=mean(CHIT)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.8, 0.93),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste("chitinase ( ", h^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##phosphatase
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(PH), y=mean(PH)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste("phosphatase ( ", h^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##phenoloxidase
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(PHO), y=mean(PHO)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.7, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste("phenoloxidase ( ", h^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##peroxidase
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(PER), y=mean(PER)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.7, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste("peroxidase ( ", h^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##################################pomery####################################################
##Ec/En
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd((BG+CELL)/(ALA+CHIT)), y=mean((BG+CELL)/(ALA+CHIT))) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.66, 0.95),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(E[C]/E[N])))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##Ec/Ep
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd((BG+CELL)/(PH)), y=mean((BG+CELL)/(PH))) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.13, 0.93),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(E[C]/E[P])))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##En/Ep
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd((ALA+CHIT)/(PH)), y=mean((ALA+CHIT)/(PH))) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.67, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(E[N]/E[P])))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##############################zmeny##########################################################
####zmeny
e.zmeny<-mich2[1:80, c("Plesne", "Certovo", "horizont", "znaceni")]
e.zmeny$dBG<-mich2[81:160, "BG"]-mich2[1:80, "BG"]
e.zmeny$dCELL<-mich2[81:160, "CELL"]-mich2[1:80, "CELL"]
e.zmeny$dALA<-mich2[81:160, "ALA"]-mich2[1:80, "ALA"]
e.zmeny$dCHIT<-mich2[81:160, "CHIT"]-mich2[1:80, "CHIT"]
e.zmeny$dPH<-mich2[81:160, "PH"]-mich2[1:80, "PH"]
e.zmeny$dPHO<-mich2[81:160, "PHO"]-mich2[1:80, "PHO"]
e.zmeny$dPER<-mich2[81:160, "PER"]-mich2[1:80, "PER"]

##BG
e.zmeny %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(dBG), y=mean(dBG)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
       legend.title = element_blank())+
  ylab(expression(paste(Delta~beta, "-glucosidase ( ", h^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)

##CELL
e.zmeny %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(dCELL), y=mean(dCELL)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste(Delta, " cellobiosidase ( ", h^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)

##ALA
e.zmeny %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(dALA), y=mean(dALA)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste(Delta, " alaninaminopeptidase ( ", h^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())

##CHIT
e.zmeny %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(dCHIT), y=mean(dCHIT)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste(Delta, " leucinaminopeptidase ( ", h^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())

##PH
e.zmeny %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(dPH), y=mean(dPH)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste(Delta, " phosphatase ( ", h^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())

######################################pomery s biomasou###################################
##Ec/Cmic
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd((BG+CELL)/Cmic), y=mean((BG+CELL)/Cmic)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.7, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(E[C]/C[MB])))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##En/Nmic
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd((ALA+CHIT)/Nmic), y=mean((ALA+CHIT)/Nmic)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.7, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(E[N]/N[MB])))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##Ep/Pmic
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd((PH)/Pmic), y=mean((PH)/Pmic)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  ylab(expression(paste(E[P]/P[MB])))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))
