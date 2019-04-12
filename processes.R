##########################################procesy#################################################
##CO2
mich2[mich2$cas==1, ] %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(C.CO2c), y=mean(C.CO2c)) %>%
ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6)+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.88),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste("Cumulative ",CO[2]," ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

# ##obrat
# model_res$pars %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(Vmax), y=mean(Vmax)) %>%
#   ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6)+
#   geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
#   facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.88),
#                                                               legend.key.size = unit(0.3, "in"),
#                                                               legend.title = element_blank())+
#  ylab(expression(paste(C[MB], " turnover rate ( ",day^{-1}, ")" )))+
#  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##CO2 glukoza
mich2[(mich2$cas==1 & mich2$znaceni=="Labelled"), ] %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(C.CO2.g), y=mean(C.CO2.g)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6)+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.88),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste("Cumulative glucose  ",CO[2]," ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##Cmic glukoza
mich2[(mich2$cas==1 & mich2$znaceni=="Labelled"), ] %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(Cmic.g), y=mean(Cmic.g)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6)+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.88),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(C[MB], " glucose ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##Priming 
mich2$Priming<-numeric(length = nrow(mich2))
mich2[(mich2$cas==1 & mich2$znaceni=="Labelled"),"Priming" ]<-mich2[(mich2$cas==1 & mich2$znaceni=="Labelled"),"C.CO2.s" ]-mich2[(mich2$cas==1 & mich2$znaceni=="Unlabelled"),"C.CO2.s" ]

mich2[(mich2$cas==1 & mich2$znaceni=="Labelled"), ] %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(Priming), y=mean(Priming)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6)+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.88),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(CO[2]," soil surplus ( ", mu, "mol ",g^{-1}, ")")))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##CUE 
mich2[(mich2$cas==1 & mich2$znaceni=="Labelled"), ] %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(CUE), y=mean(CUE)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6)+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.88),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(CUE)))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))


mich2[(mich2$cas==1 & mich2$znaceni=="Labelled"), ] %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(CUE), y=mean(CUE),
                                                                                                          x.sd=sd(pH), x=mean(pH)) %>%
  ggplot(aes(x,y))+geom_point(cex=6, pch=21)+theme_min+geom_errorbar(width=0.1, aes(ymin=y-y.sd,
                                                                                    ymax=y+y.sd))+
  geom_errorbarh(width=0.1, aes(xmin=x-x.sd, xmax=x+x.sd))+ylab("CUE")+xlab("pH")

##C residualni 
mich2[(mich2$cas==1 & mich2$znaceni=="Labelled"), ] %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(100-Gl.found), y=mean(100-Gl.found)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6)+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(~horizont)+theme_min+theme(legend.position = c(0.2, 0.88),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste("Residual C (%)")))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

ziv.zmeny$Cres<-(100-mich2[(mich2$cas==1), "Gl.found"])*mich2[(mich2$cas==1), "Gl"]/100

mich2[(mich2$cas==1 & mich2$znaceni=="Labelled"), ] %>% group_by(Plesne, horizont) %>% 
  summarize(y.sd=sd((100-Gl.found)), y=mean((100-Gl.found)), x.sd=sd(bg), x=mean(bg)) %>%
  ggplot(aes(x, y))+geom_point(pch=21, cex=6, aes(fill=horizont))+#facet_wrap(~horizont, scales="free")+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  geom_errorbarh(aes(xmin=x-x.sd, xmax=x+x.sd), width=0.1, lwd=0.5)+
  theme_min+theme(legend.position = c(0.2, 0.88),
                                        legend.key.size = unit(0.3, "in"),
                                        legend.title = element_blank())+
  ylab(expression(paste("Residual C (%)")))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))
