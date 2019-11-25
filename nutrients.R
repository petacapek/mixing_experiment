############################################Nutrients data##########################################
##DOC
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(DOC), y=mean(DOC)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.9),
                                                              legend.key.size = unit(0.3, "in"))+
  ylab(expression(paste("DOC ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

mich2 %>% filter(znaceni=="Unlabelled" & Legend=="Before incubation") %>% 
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

anova(glm(DOC~horizont+Plesne/horizont, 
          data = subset(mich2, znaceni=="Unlabelled" & Legend=="Before incubation"),
          family = Gamma), test="F")

##DOC2
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(DOC2), y=mean(DOC2)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.9),
                                                              legend.key.size = unit(0.3, "in"))+
  ylab(expression(paste(K[2],SO[4], " DOC ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##DOC2/DOC
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd((DOC2-DOC)/DOC2), y=mean((DOC2-DOC)/DOC2)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.9),
                                                              legend.key.size = unit(0.3, "in"))+
  ylab(expression(paste(K[2],SO[4], " DOC ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##NH4
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(NH4), y=mean(NH4)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.85, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(NH[4], " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

mich2 %>% filter(znaceni=="Unlabelled" & Legend=="Before incubation") %>% 
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

anova(glm(NH4~horizont+Plesne/horizont, 
          data = subset(mich2, znaceni=="Unlabelled" & Legend=="Before incubation"),
          family = Gamma), test="F")

##NH42
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(NH42), y=mean(NH42)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.85, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(NH[4], " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##NH42/NH4
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd((NH42-NH4)/NH42), y=mean((NH42-NH4)/NH42)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.85, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(NH[4], " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##NO3
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(NO3), y=mean(NO3)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.85, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(NO[3], " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

mich2 %>% filter(znaceni=="Unlabelled" & Legend=="Before incubation") %>% 
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

anova(glm(NO3~horizont+Plesne/horizont, 
          data = subset(mich2, znaceni=="Unlabelled" & Legend=="Before incubation"),
          family = Gamma), test="F")

##NO32
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(NO32), y=mean(NO32)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.85, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(NO[3], " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##NO32-NO3
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd((NO32-NO3)/NO32), y=mean((NO32-NO3)/NO32)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.85, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(NO[3], " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##DON
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(DON), y=mean(DON)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.18, 0.92),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(DON, " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

mich2 %>% filter(znaceni=="Unlabelled" & Legend=="Before incubation") %>% 
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

anova(glm(DON~horizont+Plesne/horizont, 
          data = subset(mich2, znaceni=="Unlabelled" & Legend=="Before incubation"),
          family = Gamma), test="F")

##DON2
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(DON2), y=mean(DON2)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.18, 0.92),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(K[2], SO[4]~DON, " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))
##DON2/DON
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd((DON2-DON)/DON2), y=mean((DON2-DON)/DON2)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.18, 0.92),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(K[2], SO[4]~DON, " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))



##PO4
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(PO4), y=mean(PO4)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.17, 0.4),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(PO[4], " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

mich2 %>% filter(znaceni=="Unlabelled" & Legend=="Before incubation") %>% 
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

anova(glm(PO4~horizont+Plesne/horizont, 
          data = subset(mich2, znaceni=="Unlabelled" & Legend=="Before incubation"),
          family = Gamma), test="F")

##DOP
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(DOP), y=mean(DOP)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.17, 0.4),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(DOP, " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

mich2 %>% filter(znaceni=="Unlabelled" & Legend=="Before incubation") %>% 
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

anova(glm(DOP~horizont+Plesne/horizont, 
          data = subset(mich2, znaceni=="Unlabelled" & Legend=="Before incubation"),
          family = Gamma), test="F")

##PO42
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(PO42), y=mean(PO42)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.17, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(NaHCO[3]~PO[4], " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##PO42/PO4
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd((PO42-PO4)/PO42), y=mean((PO42-PO4)/PO42)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.17, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(NaHCO[3]~PO[4], " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##DOC/NH4 
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(DOC/NH4), y=mean(DOC/NH4)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.7, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(DOC/NH[4], " (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##DOC2/NH42
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(DOC2/NH42), y=mean(DOC2/NH42)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.7, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(DOC/NH[4], " (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##DOC/NO3 
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(DOC/NO3), y=mean(DOC/NO3)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.67, 0.93),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(DOC/NO[3], " (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##DOC2/NO32 
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(DOC2/NO32), y=mean(DOC2/NO32)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.67, 0.93),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(DOC/NO[3], " (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##DOC/DON 
mich2[-151, ] %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(DOC/DON), y=mean(DOC/DON)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.67, 0.32),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(DOC/DON, " (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

mich2 %>% filter(znaceni=="Unlabelled") %>% 
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

anova(glm(DOC/DON~horizont+Plesne/horizont+cas/Plesne/horizont, 
          data = subset(mich2, znaceni=="Unlabelled" & Legend=="Before incubation"),
          family = Gamma), test="F")

mich2 %>% filter(znaceni=="Unlabelled") %>% 
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

anova(glm(DOC/DOP~horizont+Plesne/horizont+cas/Plesne/horizont, 
          data = subset(mich2, znaceni=="Unlabelled" & Legend=="Before incubation"),
          family = Gamma), test="F")

mich2 %>% filter(znaceni=="Unlabelled" & Legend=="Before incubation") %>% group_by(Plesne, horizont) %>% 
  summarize(y.sd=sd(DOC/DON), y=mean(DOC/DON)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, fill="grey")+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(~horizont, scales="free")+theme_min+theme(legend.position = c(0.67, 0.32),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(DOC/DON, " (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##DOC2/DON2
mich2[-138, ] %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(DOC2/DON2), y=mean(DOC2/DON2)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(K[2],SO[4]~DOC/DON, " (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##DOC/PO4
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(DOC/PO4), y=mean(DOC2/PO4)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(DOC/PO[4], " (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##DOC/DOP
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(DOC/DOP), y=mean(DOC/DOP)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(DOC/DOP, " (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##DOC/DOP
mich2 %>% filter(znaceni=="Unlabelled" & Legend=="Before incubation") %>% 
  group_by(Plesne, Horizon) %>% summarize(y.sd=sd(DOC/DOP), y=mean(DOC/DOP)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Horizon))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme_min+theme(legend.position = c(0.2, 0.9),
                  legend.key.size = unit(0.3, "in"),
                  legend.title = element_blank())+
  ylab(expression(paste(DOC/DOP, " (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##DOC2/PO42
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(DOC2/PO42), y=mean(DOC2/PO42)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.7, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(K[2],SO[4]~DOC/NaHCO[3]~PO[4], " (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##NH4/PO4
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(NH4/PO4), y=mean(NH4/PO4)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.32, 0.32),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(NH[4]/PO[4], " (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##NH42/PO42
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(NH42/PO42), y=mean(NH42/PO42)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.32, 0.32),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(NH[4]/PO[4], " (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##NO3/PO4
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(NO3/PO4), y=mean(NO3/PO4)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.32, 0.32),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(NO[3]/PO[4], " (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##NO32/PO42
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(NO32/PO42), y=mean(NO32/PO42)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.32, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(K[2],SO[4]~NO[3]/NaHCO[3]~PO[4], " (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##DON/PO4
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(DON/PO4), y=mean(DON/PO4)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.32, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(DON/PO[4], " (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##DON2/PO42
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(DON2/PO42), y=mean(DON2/PO42)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.8, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(K[2],SO[4]~DON/NaHCO[3]~PO[4], " (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

####zmeny
ziv.zmeny<-mich2[1:80, c("Plesne", "Certovo", "horizont", "znaceni")]
ziv.zmeny$dDOC<-mich2[81:160, "DOC"]-mich2[1:80, "DOC"]
ziv.zmeny$dDOC2<-mich2[81:160, "DOC2"]-mich2[1:80, "DOC2"]
ziv.zmeny$dNH4<-mich2[81:160, "NH4"]-mich2[1:80, "NH4"]
ziv.zmeny$dNH42<-mich2[81:160, "NH42"]-mich2[1:80, "NH42"]
ziv.zmeny$dNO3<-mich2[81:160, "NO3"]-mich2[1:80, "NO3"]
ziv.zmeny$dNO32<-mich2[81:160, "NO32"]-mich2[1:80, "NO32"]
ziv.zmeny$dDON<-mich2[81:160, "DON"]-mich2[1:80, "DON"]
ziv.zmeny$dDON2<-mich2[81:160, "DON2"]-mich2[1:80, "DON2"]
ziv.zmeny$dPO4<-mich2[81:160, "PO4"]-mich2[1:80, "PO4"]
ziv.zmeny$dPO42<-mich2[81:160, "PO42"]-mich2[1:80, "PO42"]
ziv.zmeny$CO2<-mich2[81:160, "CCO2c"]
ziv.zmeny$CO2s<-mich2[81:160, "CCO2.s"]/mich2[81:160, "CCO2c"]
ziv.zmeny$CUE<-mich2[81:160, "CUE"]
ziv.zmeny$Sorption<-mich0$Sorption
##DOC
ziv.zmeny %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(dDOC), y=mean(dDOC)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste(Delta~DOC, " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())

ggplot(ziv.zmeny, aes(1-Sorption, CO2s))+geom_point(cex=6, pch=21, aes(fill=znaceni))+theme_min+facet_grid(.~znaceni, scales="free")

##DOC2
ziv.zmeny %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(dDOC2), y=mean(dDOC2)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste(Delta~K[2],SO[4]~DOC, " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())


##NH4
ziv.zmeny %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(dNH4), y=mean(dNH4)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste(Delta~NH[4], " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())

ziv.zmeny %>% filter(znaceni=="Unlabelled") %>% group_by(Plesne, horizont) %>% 
  summarize(y.sd=sd(dNH4), y=mean(dNH4)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, fill="grey")+
  facet_wrap(~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste(Delta~NH[4], " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())

##NH4
ziv.zmeny %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(dNH42), y=mean(dNH42)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste(Delta~NH[4], " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())

##NO3
ziv.zmeny %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(dNO3), y=mean(dNO3)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.2))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.15, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste(Delta~NO[3], " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())

##NO32
ziv.zmeny %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(dNO32), y=mean(dNO32)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.4, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste(Delta~K[2],SO[4]~NO[3], " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())

##DON
ziv.zmeny %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(dDON), y=mean(dDON)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.8, 0.9),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste(Delta~DON, " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())

##DON2
ziv.zmeny %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(dDON2), y=mean(dDON2)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.17, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste(Delta~K[2],SO[4]~DON, " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())

##PO4
ziv.zmeny %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(dPO4), y=mean(dPO4)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.17, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste(Delta~PO[4], " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())

##PO42
ziv.zmeny %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(dPO42), y=mean(dPO42)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(cex=6, pch=21, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme(legend.position = c(0.17, 0.2),legend.key.size = unit(0.3, "in"),
        legend.title = element_blank())+
  ylab(expression(paste(Delta~NaHCO[3]~PO[4], " ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())

#NH4/PO4
ggplot(ziv.zmeny, aes(factor(Plesne), dNH4/dPO4))+geom_boxplot(lwd=1.2, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.8))+
  ylab(expression(paste(Delta~NH[4]/Delta~PO[4],  " (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())#+
  ylim(-5000, 100)
  
##NO3/PO4
ggplot(ziv.zmeny, aes(factor(Plesne), dNO3/dPO4))+geom_boxplot(lwd=1.2, aes(fill=znaceni))+
    facet_grid(.~horizont, scales="free")+theme_min+theme(legend.position = c(0.3, 0.2))+
    ylab(expression(paste(Delta~NO[3]/Delta~PO[4],  " (mol/mol)" )))+
    xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
    geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())+
    ylim(-30, 15)

##DOC/PO4
ggplot(ziv.zmeny, aes(factor(Plesne), dDON/dPO4))+geom_boxplot(lwd=1.2, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+theme(legend.position = c(0.3, 0.2))+
  ylab(expression(paste(Delta~DON/Delta~PO[4],  " (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())+
  ylim(-200, 150)

##NH4/NO3
ggplot(ziv.zmeny, aes(factor(Plesne), dNH4/dNO3))+geom_boxplot(lwd=1.2, aes(fill=znaceni))+
  facet_wrap(~horizont, scales="free")+theme_min+theme(legend.position = c(0.3, 0.2))+
  ylab(expression(paste(Delta~NH[4]/Delta~NO[3],  " (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())

##NH4/DON
ggplot(ziv.zmeny, aes(factor(Plesne), dNH4/dDON))+geom_boxplot(lwd=1.2, aes(fill=znaceni))+
  facet_grid(.~horizont, scales="free")+theme_min+theme(legend.position = c(0.3, 0.2))+
  ylab(expression(paste(Delta~NH[4]/Delta~DON,  " (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())+
  ylim(-20, 10)

##NO3/DON
ggplot(ziv.zmeny, aes(factor(Plesne), dNO3/dDON))+geom_boxplot(lwd=1.2, aes(fill=znaceni))+
  facet_wrap(~horizont, scales="free")+theme_min+theme(legend.position = c(0.3, 0.2))+
  ylab(expression(paste(Delta~NO[3]/Delta~DON,  " (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))+
  geom_hline(yintercept = 0, lwd=1.5)+theme(legend.title = element_blank())

###################################ziviny vuci biomase#######################################
##Cmic/DOC
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(Cmic/DOC), y=mean(Cmic/DOC)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.17, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(C[MB], "/DOC (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##Cmic/DOC2
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(Cmic/DOC2), y=mean(Cmic/DOC2)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.17, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(C[MB], "/DOC (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##Nmic/NH4
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(Nmic/NH4), y=mean(Nmic/NH4)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.15, 0.32),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(N[MB], "/", NH[4] ," (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##Nmic/NH42
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(Nmic/NH42), y=mean(Nmic/NH42)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.15, 0.32),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(N[MB], "/", NH[4] ," (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##Nmic/NO3
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(Nmic/NO3), y=mean(Nmic/NO3)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.15, 0.32),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(N[MB], "/", NO[3] ," (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##Nmic/NO32
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(Nmic/NO32), y=mean(Nmic/NO32)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.15, 0.32),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(N[MB], "/", NO[3] ," (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##Nmic/DON
mich2[-151, ] %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(Nmic/DON), y=mean(Nmic/DON)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.15, 0.9),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(N[MB], "/", DON ," (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##Nmic/DON2
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(Nmic/DON2), y=mean(Nmic/DON2)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.67, 0.92),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(N[MB], "/", K[2],SO[4]~DON ," (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##Nmic/TN
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(Nmic/TN), y=mean(Nmic/TN)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.15, 0.88),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(N[MB], "/", TN ," (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##Nmic/TN2
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(Nmic/TN2), y=mean(Nmic/TN2)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.15, 0.88),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(N[MB], "/", TN ," (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##Pmic/PO4
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(Pmic/PO4), y=mean(Pmic/PO4)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.15, 0.92),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(P[MB], "/", PO[4] ," (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##Pmic/PO42
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(Pmic/PO42), y=mean(Pmic/PO42)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.15, 0.33),
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(P[MB], "/", PO[4] ," (mol/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

###################################ziviny vuci enzymum#######################################
##BG+CELL/DOC
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd((BG+CELL)/DOC), y=mean((BG+CELL)/DOC)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = "bottom",
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(E[C], "/", DOC ," (", h^{-1}, "/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##BG+CELL/DOC2
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd((BG+CELL)/DOC2), y=mean((BG+CELL)/DOC2)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = "bottom",
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(E[C], "/", DOC ," (", h^{-1}, "/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##ALA+CHIT/NH4
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd((ALA+CHIT)/NH4), y=mean((ALA+CHIT)/NH4)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = "bottom",
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(E[N], "/", NH[4] ," (", h^{-1}, "/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##ALA+CHIT/NH42
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd((ALA+CHIT)/NH42), y=mean((ALA+CHIT)/NH42)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = "bottom",
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(E[N], "/", NH[4] ," (", h^{-1}, "/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##ALA+CHIT/NO3
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd((ALA+CHIT)/NO3), y=mean((ALA+CHIT)/NO3)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = "bottom",
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(E[N], "/", NO[3] ," (", h^{-1}, "/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##ALA+CHIT/NO32
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd((ALA+CHIT)/NO32), y=mean((ALA+CHIT)/NO32)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = "bottom",
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(E[N], "/", NO[3] ," (", h^{-1}, "/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##ALA+CHIT/DON
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd((ALA+CHIT)/DON), y=mean((ALA+CHIT)/DON)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = "bottom",
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(E[N], "/", DON ," (", h^{-1}, "/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##ALA+CHIT/DON2
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd((ALA+CHIT)/DON2), y=mean((ALA+CHIT)/DON2)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = "bottom",
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(E[N], "/", DON ," (", h^{-1}, "/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##ALA+CHIT/TN
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd((ALA+CHIT)/TN), y=mean((ALA+CHIT)/TN)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = "bottom",
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(E[N], "/", TN ," (", h^{-1}, "/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##ALA+CHIT/TN2
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd((ALA+CHIT)/TN2), y=mean((ALA+CHIT)/TN2)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = "bottom",
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(E[N], "/", TN ," (", h^{-1}, "/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##PH/PO4
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd((PH)/PO4), y=mean((PH)/PO4)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = "bottom",
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(E[P], "/", PO[4] ," (", h^{-1}, "/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##PH/PO42
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd((PH)/PO42), y=mean((PH)/PO42)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(znaceni~horizont, scales="free")+theme_min+theme(legend.position = "bottom",
                                                              legend.key.size = unit(0.3, "in"),
                                                              legend.title = element_blank())+
  ylab(expression(paste(E[P], "/", NaHCO[3]~PO[4] ," (", h^{-1}, "/mol)" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))


##pH
mich2 %>% group_by(Plesne, znaceni, Horizon, Legend) %>% summarize(y.sd=sd(pH), y=mean(pH)) %>%
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

anova(glm(pH~horizont+Plesne/horizont, mich2, family = Gamma), test="F")
  
##W
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(DWovlhceni), y=mean(DWovlhceni)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=horizont))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  theme_min+theme(legend.position = c(0.1, 0.8),
                  legend.key.size = unit(0.3, "in"),
                  legend.title = element_blank())+
  ylab(expression(paste("DW")))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

##Sorpce
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(1-Sorption), y=mean(1-Sorption)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=horizont), show.legend = F)+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(.~horizont)+
  theme_min+theme(legend.position = c(0.1, 0.8),
                  legend.key.size = unit(0.3, "in"),
                  legend.title = element_blank())+
  ylab(expression(paste("Relative sorption capacity")))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

zmeny_vse<-cbind(mb.zmeny[, -c(1:4)], ziv.zmeny)
zmeny_vse$Cmic<-mich0$Cmic
zmeny_vse$Nmic<-mich0$Nmic
zmeny_vse$Pmic<-mich0$Pmic

ggplot(zmeny_vse[zmeny_vse$znaceni=="Unlabelled", ], aes(Cmic/Nmic, dDOC2/(dDON2+dNH4+dNO3)))+
  geom_point(cex=6, aes(colour=horizont))+ylim(0,20)+
  geom_abline(intercept = 0, slope=1)+facet_wrap(~horizont, scales="free")

#Spocitam kriticky pomery a porovnam je s pomerama pudnich zivin a treba mi neco 
#z toho vysvetli ztratu biomasy

zmeny_vse$CNcrit<-mich0$Cmic/mich0$Nmic/
  c(mich1[mich1$znaceni=="Labelled", "CUE"],mich1[mich1$znaceni=="Labelled", "CUE"])
zmeny_vse$CPcrit<-mich0$Cmic/mich0$Pmic/
  c(mich1[mich1$znaceni=="Labelled", "CUE"],mich1[mich1$znaceni=="Labelled", "CUE"])


zmeny_vse %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(CNcrit), y=mean(CNcrit)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, fill="grey", show.legend = F)+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(~horizont, scales = "free")+
  theme_min+theme(legend.position = c(0.1, 0.8),
                  legend.key.size = unit(0.3, "in"),
                  legend.title = element_blank())+
  ylab(expression(paste("C:N critical")))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

zmeny_vse %>% group_by(Plesne, znaceni, horizont) %>% summarize(y.sd=sd(CPcrit), y=mean(CPcrit)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, fill="grey", show.legend = F)+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_wrap(~horizont, scales = "free")+
  theme_min+theme(legend.position = c(0.1, 0.8),
                  legend.key.size = unit(0.3, "in"),
                  legend.title = element_blank())+
  ylab(expression(paste("C:P critical")))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

#Ziviny
zmeny_vse$CO2un<-mich1$CCO2
zmeny_vse$CO2gl<-c(mich1[mich1$znaceni=="Labelled", "CCO2"]*((mich1[mich1$znaceni=="Labelled", "CO2atm"]-mich1[mich1$znaceni=="Unlabelled", "CO2atm"])/
                                                             (0.06296353-mich1[mich1$znaceni=="Unlabelled", "CO2atm"])),
                   mich1[mich1$znaceni=="Labelled", "CCO2"]*((mich1[mich1$znaceni=="Labelled", "CO2atm"]-mich1[mich1$znaceni=="Unlabelled", "CO2atm"])/
                                                               (0.06296353-mich1[mich1$znaceni=="Unlabelled", "CO2atm"])))
zmeny_vse$Cmicgl<-c(mich1[mich1$znaceni=="Labelled", "Cmicg"],
                    mich1[mich1$znaceni=="Labelled", "Cmicg"])
zmeny_vse$CUE2<-with(zmeny_vse, Cmicgl/0.38/(CO2gl+Cmicgl/0.38))
zmeny_vse$dTN<-mich2[81:160, "TN"]-mich2[1:80, "TN"]
zmeny_vse$dTN2<-mich2[81:160, "TN2"]-mich2[1:80, "TN2"]
zmeny_vse$NH4<-mich0$NH4
zmeny_vse$NH42<-mich0$NH42
zmeny_vse$DON<-mich0$DON
zmeny_vse$DONb<-mich1$DON
zmeny_vse$DON2<-mich0$DON2
zmeny_vse$PO4<-mich0$PO4
zmeny_vse$PO42<-mich0$PO42
zmeny_vse$NO3<-mich0$NO3
zmeny_vse$NO32<-mich0$NO32
zmeny_vse$DOC<-mich0$DOC
zmeny_vse$DOCb<-mich1$DOC
zmeny_vse$DOC2<-mich0$DOC2
zmeny_vse$gluVmax<-mich2e[mich2e$cas==0, "gluVmax"]
zmeny_vse$celVmax<-mich2e[mich2e$cas==0, "celVmax"]
zmeny_vse$chitVmax<-mich2e[mich2e$cas==0, "chitVmax"]
zmeny_vse$leuVmax<-mich2e[mich2e$cas==0, "leuVmax"]
zmeny_vse$phosVmax<-mich2e[mich2e$cas==0, "phosVmax"]
zmeny_vse$phosKm<-mich2e[mich2e$cas==0, "phosKm"]
zmeny_vse$PO4sorbed<-with(zmeny_vse, PO4*(1-exp(-Sorption/0.75*48)))
zmeny_vse$PO42sorbed<-with(zmeny_vse, PO42*(1-exp(-Sorption/0.75*48)))
zmeny_vse$DOP<-mich2[c(1:80), "DOP"]
zmeny_vse$dDOP<-mich2[c(81:160), "DOP"]-mich2[c(1:80), "DOP"]
zmeny_vse$Ndemand<-with(zmeny_vse, Dekompozice*(Nmic*0.38*CUE2/Cmic/0.54))
zmeny_vse$Pdemand<-with(zmeny_vse, Dekompozice*(Pmic*0.38*CUE2/Cmic/0.4))
zmeny_vse$Pdemand2<-with(zmeny_vse, Dekompozice*(CUE2/57.163))
zmeny_vse$pH<-mich0$pH
zmeny_vse$CNcrit<-with(zmeny_vse, Cmic*0.54/Nmic/0.38/CUE2)
zmeny_vse$CPcrit<-with(zmeny_vse, Cmic*0.4/Pmic/0.38/CUE2)
zmeny_vse$Bdist<-merge(zmeny_vse, bdiff, by.x = c("Plesne", "horizont"),
                       by.y = c("Plesne", "Horizon"))[, 65]
zmeny_vse$Bdist.sd<-merge(zmeny_vse, bdiff, by.x = c("Plesne", "horizont"),
                       by.y = c("Plesne", "Horizon"))[, 67]
zmeny_vse$Fdist<-merge(zmeny_vse, fdiff, by.x = c("Plesne", "horizont"),
                       by.y = c("Plesne", "Horizon"))[, 67]
zmeny_vse$Fdist.sd<-merge(zmeny_vse, fdiff, by.x = c("Plesne", "horizont"),
                       by.y = c("Plesne", "Horizon"))[, 69]
zmeny_vse$Ptot<-vector("numeric", length = 80)

for(i in 1:80){
  if(zmeny_vse$horizont[i]=="O"){
    zmeny_vse$Ptot[i]<-81*zmeny_vse$Plesne[i]+76*zmeny_vse$Certovo[i]
  }else{
    zmeny_vse$Ptot[i]<-54*zmeny_vse$Plesne[i]+51*zmeny_vse$Certovo[i]
  }
}

zmeny_vse$PON<-with(zmeny_vse, Ptot-PO4)
zmeny_vse$obrat<-Obrat$obrat
zmeny_vse$obrat2<-Obrat$obrat2
zmeny_vse$BNMDS1a<-c(t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 1, cols = c(2))),
                    t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 1, cols = c(2))))
zmeny_vse$BNMDS2a<-c(t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 1, cols = c(3))),
                     t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 1, cols = c(3))))
zmeny_vse$BNMDS3a<-c(t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 1, cols = c(4))),
                     t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 1, cols = c(4))))
zmeny_vse$FNMDS1a<-c(t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 4, cols = c(2))),
                     t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 4, cols = c(2))))
zmeny_vse$FNMDS2a<-c(t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 4, cols = c(3))),
                     t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 4, cols = c(3))))
zmeny_vse$FNMDS3a<-c(t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 4, cols = c(4))),
                     t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 4, cols = c(4))))

zmeny_vse$BNMDS1b<-c(t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 2, cols = c(2))),
                     t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 3, cols = c(2))),
                     t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 2, cols = c(2))),
                     t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 3, cols = c(2))))
zmeny_vse$BNMDS2b<-c(t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 2, cols = c(3))),
                     t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 3, cols = c(3))),
                     t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 2, cols = c(3))),
                     t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 3, cols = c(3))))
zmeny_vse$BNMDS3b<-c(t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 2, cols = c(4))),
                     t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 3, cols = c(4))),
                     t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 2, cols = c(4))),
                     t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 3, cols = c(4))))
zmeny_vse$FNMDS1b<-c(t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 5, cols = c(2))),
                     t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 6, cols = c(2))),
                     t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 5, cols = c(2))),
                     t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 6, cols = c(2))))
zmeny_vse$FNMDS2b<-c(t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 5, cols = c(3))),
                     t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 6, cols = c(3))),
                     t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 5, cols = c(3))),
                     t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 6, cols = c(3))))
zmeny_vse$FNMDS3b<-c(t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 5, cols = c(4))),
                     t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 6, cols = c(4))),
                     t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 5, cols = c(4))),
                     t(read.xlsx("./Chemical_and_microbial_data/mich_exp_nmds_scores.xlsx", 6, cols = c(4))))
                     

summary(lm(CUE2~CUE, zmeny_vse))

ggplot(zmeny_vse[zmeny_vse$znaceni=="Unlabelled", ], aes(NH4/PO4, CNcrit))+
  geom_point(cex=6, pch=21, aes(fill=horizont))+
  theme_min

ggplot(zmeny_vse[zmeny_vse$znaceni=="Unlabelled", ], aes(Ndemand, dTN2))+
  geom_point(cex=6, pch=21, aes(fill=horizont))+
  theme_min+geom_abline(intercept = 0, slope = -1)

#Co tvori rozdil mezi predpokladanou a pocitanou hodnotu?
zmeny_vse$Nrozdil<-with(zmeny_vse, Ndemand+dTN2)
ggplot(zmeny_vse[zmeny_vse$znaceni=="Unlabelled", ], aes(deltaCO2-deltaCmic, Nrozdil))+
  geom_point(cex=6, pch=21, aes(fill=horizont))



ggplot(zmeny_vse[(zmeny_vse$znaceni=="Unlabelled" ), ], aes(FNMDS3b, obrat2))+
  geom_point(cex=6, pch=21, aes(fill=horizont))+facet_wrap(~horizont, scales="free")+
  theme_min+stat_smooth(method = lm, se=F, aes(colour=horizont))

summary(lm(obrat2~BNMDS1b, zmeny_vse[(zmeny_vse$znaceni=="Unlabelled" & zmeny_vse$obrat2<0.06), ], subset=horizont=="O"))
summary(lm(obrat2~FNMDS3b, zmeny_vse[(zmeny_vse$znaceni=="Unlabelled" ), ], subset=horizont=="A"))


zmeny_vse$Horizon<-zmeny_vse$horizont
zmeny_vse$Horizon<-ifelse(zmeny_vse$horizont=="O", "Litter", "Organic topsoil")

(BO1<-ggplot(zmeny_vse[(zmeny_vse$znaceni=="Unlabelled" & zmeny_vse$horizont=="O"), ],
             aes(BNMDS1b, obrat2))+geom_point(cex=6, pch=21, fill="grey")+
    theme_min+stat_smooth(method = lm, se=F, color="black")+facet_grid(.~Horizon)+
    xlab("16S NMDS axis 1")+ylab(expression(paste("turnover rate (", h^{-1}, ")")))+
    annotate("text", 0.12, 0.06, label = "p = 0.002", size=6, fontface="bold.italic")
  )
(BA1<-ggplot(zmeny_vse[(zmeny_vse$znaceni=="Unlabelled" &  zmeny_vse$horizont=="A"), ],
             aes(FNMDS3b, obrat2))+geom_point(cex=6, pch=21, fill="white")+
    theme_min+stat_smooth(method = lm, se=F, color="black")+facet_grid(.~Horizon, drop = T)+
    xlab("ITS NMDS axis 3")+ylab(expression(paste("turnover rate (", h^{-1}, ")")))+
    annotate("text", 0.15, 0.07, label = "p = 0.026", size=6, fontface="bold.italic")
)



grid.arrange(BO1, BA1, ncol=2)

summary(nls(CPcrit~(NH4/PO4)*max/(NH4/PO4+K), zmeny_vse[zmeny_vse$znaceni=="Unlabelled", ],
            start = list(max=60, K=50)))

tab1<-as.data.frame(mich2 %>% filter(znaceni=="Unlabelled" & Legend=="Before incubation") %>% 
                      group_by(horizont, Plesne) %>%
                      summarize(pH=mean(pH, na.rm=T), pHsd=sd(pH, na.rm=T),
            MBC = mean(Cmic/0.38, na.rm=T), MBCsd = sd(Cmic/0.38, na.rm=T),
            MBN = mean(Nmic/0.54, na.rm=T), MBNsd = sd(Nmic/0.54, na.rm=T),
            MBP = mean(Pmic/0.4, na.rm=T), MBPsd = sd(Pmic/0.4, na.rm=T),
            WEC = mean(DOC, na.rm=T), WECsd = sd(DOC, na.rm=T),
            DON = mean(DON, na.rm=T), DONsd = sd(DON, na.rm=T),
            NH4 = mean(NH4, na.rm=T), NH4sd = sd(NH4, na.rm=T),
            NO3 = mean(NO3, na.rm=T), NO3sd = sd(NO3, na.rm=T),
            SRP = mean(PO4, na.rm=T), PO4sd = sd(PO4, na.rm=T)))


tab1$pHsd<-as.data.frame(mich2 %>% filter(znaceni=="Unlabelled" & Legend=="Before incubation") %>% 
                           group_by(horizont, Plesne) %>%
                           summarize(pHsd=sd(pH, na.rm=T)))[,3]
tab1$DONsd<-as.data.frame(mich2 %>% filter(znaceni=="Unlabelled" & Legend=="Before incubation") %>% 
                           group_by(horizont, Plesne) %>%
                           summarize(pHsd=sd(DON, na.rm=T)))[,3]
tab1$NH4sd<-as.data.frame(mich2 %>% filter(znaceni=="Unlabelled" & Legend=="Before incubation") %>% 
                            group_by(horizont, Plesne) %>%
                            summarize(pHsd=sd(NH4, na.rm=T)))[,3]
tab1$NO3sd<-as.data.frame(mich2 %>% filter(znaceni=="Unlabelled" & Legend=="Before incubation") %>% 
                            group_by(horizont, Plesne) %>%
                            summarize(pHsd=sd(NO3, na.rm=T)))[,3]

tab1x<-tab1[order(tab1$horizont, tab1$Plesne, decreasing = T), ]

write.xlsx(tab1x, file = c("C:/Users/cape159/Documents/pracovni/publikace/Michani/Figures/tab1.xlsx"))


as.data.frame(mich2 %>% filter(znaceni=="Unlabelled" & Legend=="Before incubation") %>% 
  group_by(horizont, Plesne) %>%
  summarize(av=mean(DOP*1000, na.rm=T), s = sd(DOP*1000, na.rm=T)))
