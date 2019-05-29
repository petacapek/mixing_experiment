############################################Nutrients data##########################################
##DOC
mich2 %>% group_by(Plesne, znaceni, horizont, Legend) %>% summarize(y.sd=sd(DOC), y=mean(DOC)) %>%
  ggplot(aes(factor(Plesne), y))+geom_point(pch=21, cex=6, aes(fill=Legend))+
  geom_errorbar(aes(ymin=y-y.sd, ymax=y+y.sd), width=0.1, lwd=0.5)+
  facet_grid(znaceni~horizont, scales="free")+theme_min+theme(legend.position = c(0.2, 0.9),
                                                              legend.key.size = unit(0.3, "in"))+
  ylab(expression(paste("DOC ( ", mu, "mol ",g^{-1}, ")" )))+
  xlab("Plesne : Certovo mixing ratio")+scale_fill_manual(values = c("white", "grey"))

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
ziv.zmeny$CO2s<-mich2[81:160, "CCO2s"]
ziv.zmeny$CUE<-mich2[81:160, "CUE"]

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



  
