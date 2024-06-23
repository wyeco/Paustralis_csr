
setwd("./Intraspecific_csr")

# download the data and import 
phenotypic.group.final <- read.csv("./data_Pa_CSR.csv", header = T, sep = ';')
str(phenotypic.group.final)

## here the df is ready
library(ggplot2)
library(ggtern)
library(viridis)

## Fig. 1  ##
ggtern(data= phenotypic.group.mean.final , aes(R., C., S., color = group,  shape = group )) +
  geom_point(size=3) +
   scale_color_manual(values=c("#1b9e77", "#d95f02", "#66a61e", "#7570b3", "#e7298a", "#e6ab02"))+
  scale_shape_manual(values = c(16, 20, 6, 7, 18, 10))+
  theme_bw() +
  theme_showarrows()  +
  labs(x="R (%)",y="C (%)",z="S (%)") +
  #theme(legend.position = c(0.8,0.9))+
  
   tern_limits(T=.6,L=0.36,R=0.8) +## subset scale
   theme_rgbw() +  theme_clockwise() +
   theme(legend.position=c(0.1,0.9), legend.justification=c(0,1))
 
 pdf("Fig.CSR_ggtern.pdf", useDingbats=FALSE, width=8, height=8)
 
 dev.off()
 
 
 ## regression along latitude  ######
 s_latitude_model <- lm(S. ~ abs(Latitude2), phenotypic.group.mean.final)
 plot(s_latitude_model)
 summary(s_latitude_model)
 
p_latitude_S <-  ggplot(phenotypic.group.mean.final, aes(x = abs(Latitude2), y = S.)) +
   geom_point(size =2, aes(color = group, shape = group)) +
   scale_color_manual(values=c("#1b9e77", "#d95f02", "#66a61e", "#7570b3", "#e7298a", "#e6ab02"))+
   scale_shape_manual(values = c(16, 20, 6, 7, 18, 10))+
   geom_smooth(method = "lm", se = T, color = "red", formula = my.formula)+
   stat_cor(label.x = 40, label.y = 71) +
   stat_regline_equation(
     aes(label =  paste(..eq.label.., ..adj.rr.label..,  sep = "~~~~")),
     formula = my.formula,label.x = 40, label.y = 73
   )+
   ggtitle("(b)") +
   theme_bw() + 
   ylab("S-score (%)") +
   xlab("Latitude") + 
   theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         strip.background = element_rect(colour="white", fill="white"),
         text = element_text(size = 16),
         legend.position = "none",
         axis.text.x = element_text(face="bold", size=14),
         axis.text.y = element_text(face="bold", size=14)) +
   theme(axis.title = element_text())
 
 C_latitude_model <- lm(C. ~ abs(Latitude2), phenotypic.group.mean.final)
 plot(C_latitude_model)
 summary(C_latitude_model)
 
 library(ggpubr)
 my.formula <- y ~ x
 
 p_latitude_C <-  ggplot(phenotypic.group.mean.final, aes(x = abs(Latitude2), y = C.)) +
   geom_point(size =2, aes(color = group, shape = group)) +
   scale_color_manual(values=c("#1b9e77", "#d95f02", "#66a61e", "#7570b3", "#e7298a", "#e6ab02"))+
   scale_shape_manual(values = c(16, 20, 6, 7, 18, 10))+
   geom_smooth(method = "lm", se = T, color = "red", formula = my.formula)+
   stat_cor(label.x = 40, label.y = 28) +
   # stat_regline_equation(
   #  # aes(label =  paste(..eq.label.., ..adj.rr.label..,  sep = "~~~~")),
   #   formula = my.formula,label.x = 40, label.y = 30
   # )+
   ggtitle("(a)") +
   theme_bw() + 
   ylab("C-score (%)") +
   xlab("Latitude") + 
   theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         strip.background = element_rect(colour="white", fill="white"),
         text = element_text(size = 16),
         legend.position = c(0.85,0.8),
         axis.text.x = element_text(face="bold", size=14),
         axis.text.y = element_text(face="bold", size=14)) +
   theme(axis.title = element_text())
 
 R_latitude_model <- lm(R. ~ Latitude2, phenotypic.group.mean.final)
 plot(R_latitude_model)
 summary(R_latitude_model)
 
 
 p_latitude_R <-  ggplot(phenotypic.group.mean.final, aes(x = abs(Latitude2), y = R.)) +
   geom_point(size =2, aes(color = group, shape = group)) +
   scale_color_manual(values=c("#1b9e77", "#d95f02", "#66a61e", "#7570b3", "#e7298a", "#e6ab02"))+
   scale_shape_manual(values = c(16, 20, 6, 7, 18, 10))+
   geom_smooth(method = "lm", se = T, color = "red", formula = my.formula)+
   stat_cor(label.x = 36, label.y = 8) +
   stat_regline_equation(
     aes(label =  paste(..eq.label.., ..adj.rr.label..,  sep = "~~~~")),
     formula = my.formula,label.x = 36, label.y = 9
   )+
  # ggtitle("(c)") +
   theme_bw() + 
   ylab("R-score (%)") +
   xlab("Latitude") + 
   theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         strip.background = element_rect(colour="white", fill="white"),
         text = element_text(size = 16),
         legend.position = "none",
         axis.text.x = element_text(face="bold", size=14),
         axis.text.y = element_text(face="bold", size=14)) +
   theme(axis.title = element_text())
 
 library(cowplot)
csr_latitude = plot_grid(p_latitude_C, p_latitude_S, ncol = 1)

pdf("Fig.2 CS vs latitude.pdf", useDingbats=FALSE, width=6, height=8)
 csr_latitude
 dev.off()
 
 pdf("Fig.S1 R vs latitude.pdf", useDingbats=FALSE, width=6, height=4)
 p_latitude_R
 dev.off()

 ##  ends here 
 
 ### plot the Fig 3 and table 1 ########
### remove the NAs in genome size data (all are SA samples)
library(tidyr)
phenotypic.group.mean.noNA.final <- phenotypic.group.mean.final %>% drop_na(Cx)

# C model
c_model <- lm(C. ~ Cx +biocru_z_s1 + biocru_z_s2 + biocru_z_s3, phenotypic.group.mean.final) ## Cx has several NAs, 
## thus the final data size is 84.
vif(c_model)
plot(c_model)
summary(c_model)

install.packages("relaimpo", dep=c("Depends"))  
library(relaimpo) 
library(car)
calc.relimp( c_model, type = c("first", "betasq", "car", "lmg", "genizi", "last"), rela = T )

## plot the partial residuals for each significant factor
## Load effects
library(effects)
## Plot effect of each term
allEffects.c.full <- allEffects(c_model )
plot(allEffects.c.full)# all four predictors
plot(allEffects.c.full, selection=1) ## plot Cx only
# Cx
cx.c.full <- effect("Cx",c_model, partial.residuals =T)
plot(cx.c.full, smooth.residuals=F)

closest <- function(x, x0) apply(outer(x, x0, FUN=function(x, x0) abs(x - x0)), 1, which.min)

x.fit <- unlist(cx.c.full$x.all)
trans <- I
x <- data.frame(lower = cx.c.full$lower, upper = cx.c.full$upper, fit = cx.c.full$fit, Cx = cx.c.full$x$Cx)
xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$Cx)] + cx.c.full$residuals)
xy$group <-phenotypic.group.mean.noNA.final$group

c.cx <- ggplot(x, aes(x = Cx, y = fit)) +
  geom_line(size = 1, col = "red") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
 # geom_point(data = xy, aes(x = x, y = y, color = group, shape = group),  size = 2) +
  geom_jitter(data = xy, aes(x = x, y = y, color = group, shape = group),  size = 2, width = 0.0008) +
    scale_color_manual(values=c("#1b9e77", "#d95f02", "#66a61e", "#7570b3", "#e7298a", "#e6ab02"))+
  scale_shape_manual(values = c(16, 20, 6, 7, 18, 10))+
    ggtitle("(a)") +
  theme_bw() + 
 # coord_cartesian(ylim = c(20,50))+ 
  ylab("C-score (%)") +
  xlab("Monoploid genome size (pg)") + 
  theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      strip.background = element_rect(colour="white", fill="white"),
                      text = element_text(size = 16),
                      legend.position = "none",
                      axis.text.x = element_text(face="bold", size=14),
                      axis.text.y = element_text(face="bold", size=14)) +
  theme(axis.title = element_text())


# biocru_z_s2
c.biocru_z_s2.full <- effect("biocru_z_s2",c_model, partial.residuals =T)
plot(c.biocru_z_s2.full, smooth.residuals=F)
library(ggplot2)
library(gridExtra)
x.C.biocru_z_s2.fit <- unlist(c.biocru_z_s2.full$x.all)
x.c.biocru_z_s2 <- data.frame(lower = c.biocru_z_s2.full$lower, upper = c.biocru_z_s2.full$upper, fit = c.biocru_z_s2.full$fit, biocru_z_s2 = c.biocru_z_s2.full$x$biocru_z_s2)
xy.c.biocru_z_s2 <- data.frame(x = x.C.biocru_z_s2.fit, y = x.c.biocru_z_s2$fit[closest(trans(x.C.biocru_z_s2.fit), x.c.biocru_z_s2$biocru_z_s2)] + c.biocru_z_s2.full$residuals)
xy.c.biocru_z_s2$group <-phenotypic.group.mean.noNA.final$group
c.pc2 <- ggplot(x.c.biocru_z_s2, aes(x = biocru_z_s2, y = fit)) +
  geom_line(size = 1, col = "red") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_point(data = xy.c.biocru_z_s2, aes(x = x, y = y, color = group, shape = group), size = 2) +
  scale_color_manual(values=c("#1b9e77", "#d95f02", "#66a61e", "#7570b3", "#e7298a", "#e6ab02"))+
  scale_shape_manual(values = c(16, 20, 6, 7, 18, 10))+
    ggtitle("(b)") +
 # coord_cartesian(ylim = c(20,50))+  
  theme_bw() + 
  ylab("C-score (%)") +
  xlab("Climate PC2") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        text = element_text(size = 16),
        legend.position = "none",
        axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14)) +
  theme(axis.title = element_text())


pdf("Fig.C-partialRedisual_regression.pdf", useDingbats=FALSE, width=8, height=4)
dev.off()


###  S model
s_model <- lm(S. ~ Cx +biocru_z_s1 + biocru_z_s2 + biocru_z_s3, phenotypic.group.mean.final)
plot(s_model)
vif(s_model)
summary(s_model)
allEffects.s.full <- allEffects(s_model)
plot(allEffects.s.full)
## relative importance
calc.relimp( s_model, type = c("first", "betasq", "car", "lmg", "genizi", "last"), rela = T )

# Cx
cx.s.full <- effect("Cx",s_model, partial.residuals =T)
plot(cx.s.full, smooth.residuals=F)

x.s.cx.fit <- unlist(cx.s.full$x.all)
x.s.cx <- data.frame(lower = cx.s.full$lower, upper = cx.s.full$upper, fit = cx.s.full$fit, Cx = cx.s.full$x$Cx)
xy.s.cx <- data.frame(x = x.s.cx.fit, y = x.s.cx$fit[closest(trans(x.s.cx.fit), x.s.cx$Cx)] + cx.s.full$residuals)
xy.s.cx$group <-phenotypic.group.mean.noNA.final$group


s.cx <- ggplot(x.s.cx, aes(x = Cx, y = fit)) +
  geom_line(size = 1, col = "red") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  # geom_point(data = xy, aes(x = x, y = y, color = group, shape = group),  size = 2) +
  geom_jitter(data = xy.s.cx, aes(x = x, y = y, color = group, shape = group),  size = 2, width = 0.0008) +
  scale_color_manual(values=c("#1b9e77", "#d95f02", "#66a61e", "#7570b3", "#e7298a", "#e6ab02"))+
  scale_shape_manual(values = c(16, 20, 6, 7, 18, 10))+
  ggtitle("(c)") +
  theme_bw() + 
  ylab("S-score (%)") +
  xlab("Monoploid genome size (pg)") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        text = element_text(size = 16),
        legend.position = "none",
        axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14)) +
  theme(axis.title = element_text())

### S.biocru_z_s1
s.biocru_z_s1.full <- effect("biocru_z_s1",s_model, partial.residuals =T)
plot(s.biocru_z_s1.full, smooth.residuals=F)

x.s.biocru_z_s1.fit <- unlist(s.biocru_z_s1.full$x.all)
x.s.biocru_z_s1 <- data.frame(lower = s.biocru_z_s1.full$lower, upper = s.biocru_z_s1.full$upper, fit = s.biocru_z_s1.full$fit, biocru_z_s1 = s.biocru_z_s1.full$x$biocru_z_s1)
xy.s.biocru_z_s1 <- data.frame(x = x.s.biocru_z_s1.fit, y = x.s.biocru_z_s1$fit[closest(trans(x.s.biocru_z_s1.fit), x.s.biocru_z_s1$biocru_z_s1)] + s.biocru_z_s1.full$residuals)
xy.s.biocru_z_s1$group <-phenotypic.group.mean.noNA.final$group

s.pc1 <- ggplot(x.s.biocru_z_s1, aes(x = biocru_z_s1, y = fit)) +
  geom_line(size = 1, col = "red") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_point(data = xy.s.biocru_z_s1, aes(x = x, y = y, color = group, shape = group), size = 2) +
  scale_color_manual(values=c("#1b9e77", "#d95f02", "#66a61e", "#7570b3", "#e7298a", "#e6ab02"))+
  scale_shape_manual(values = c(16, 20, 6, 7, 18, 10))+
  ggtitle("(d)") +
  theme_bw() + 
  ylab("S-score (%)") +
  xlab("Climate PC1") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        text = element_text(size = 16),
        legend.position = "none",
        axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14)) +
  theme(axis.title = element_text())


### S.biocru_z_s2
s.biocru_z_s2.full <- effect("biocru_z_s2",s_model, partial.residuals =T)
plot(s.biocru_z_s2.full, smooth.residuals=F)

x.s.biocru_z_s2.fit <- unlist(s.biocru_z_s2.full$x.all)
x.s.biocru_z_s2 <- data.frame(lower = s.biocru_z_s2.full$lower, upper = s.biocru_z_s2.full$upper, fit = s.biocru_z_s2.full$fit, biocru_z_s2 = s.biocru_z_s2.full$x$biocru_z_s2)
xy.s.biocru_z_s2 <- data.frame(x = x.s.biocru_z_s2.fit, y = x.s.biocru_z_s2$fit[closest(trans(x.s.biocru_z_s2.fit), x.s.biocru_z_s2$biocru_z_s2)] + s.biocru_z_s2.full$residuals)
xy.s.biocru_z_s2$group <-phenotypic.group.mean.noNA.final$group

s.pc2 <- ggplot(x.s.biocru_z_s2, aes(x = biocru_z_s2, y = fit)) +
  geom_line(size = 1, col = "red") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_point(data = xy.s.biocru_z_s2, aes(x = x, y = y, color = group, shape = group),  size = 2) +
  scale_color_manual(values=c("#1b9e77", "#d95f02", "#66a61e", "#7570b3", "#e7298a", "#e6ab02"))+
  scale_shape_manual(values = c(16, 20, 6, 7, 18, 10))+
  ggtitle("(e)") +
  theme_bw() + 
  ylab("S-score (%)") +
  xlab("Climate PC2") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        text = element_text(size = 16),
        legend.position = "none",
        axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14)) +
  theme(axis.title = element_text())

grid.arrange(s.cx, s.pc1, s.pc2, ncol = 2)
pdf("Fig.S-partialRedisual_regression.pdf", useDingbats=FALSE, width=8, height=8)
dev.off()



# R model
r_model <- lm(R. ~ Cx +biocru_z_s1 + biocru_z_s2 + biocru_z_s3, phenotypic.group.mean.final)
plot(r_model)
vif(r_model)

summary(r_model)
### plot the regressions
allEffects.r.full <- allEffects(r_model)
plot(allEffects.r.full)
## relative importance
calc.relimp( r_model, type = c("first", "betasq", "car", "lmg", "genizi", "last"), rela = T )

# Cx
cx.r.full <- effect("Cx",r_model, partial.residuals =T)
plot(cx.r.full, smooth.residuals=F)

x.r.cx.fit <- unlist(cx.r.full$x.all)
x.r.cx <- data.frame(lower = cx.r.full$lower, upper = cx.r.full$upper, fit = cx.r.full$fit, Cx = cx.r.full$x$Cx)
xy.r.cx <- data.frame(x = x.r.cx.fit, y = x.r.cx$fit[closest(trans(x.r.cx.fit), x.r.cx$Cx)] + cx.r.full$residuals)
xy.r.cx$group <-phenotypic.group.mean.noNA.final$group



r.cx <- ggplot(x.r.cx, aes(x = Cx, y = fit)) +
  geom_line(size = 1, col = "red") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  # geom_point(data = xy, aes(x = x, y = y, color = group, shape = group),  size = 2) +
  geom_jitter(data = xy.r.cx, aes(x = x, y = y, color = group, shape = group),  size = 2, width = 0.0008) +
  scale_color_manual(values=c("#1b9e77", "#d95f02", "#66a61e", "#7570b3", "#e7298a", "#e6ab02"))+
  scale_shape_manual(values = c(16, 20, 6, 7, 18, 10))+
  ggtitle("(a)") +
  theme_bw() + 
  ylab("R-score (%)") +
  xlab("Monoploid genome size (pg)") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        text = element_text(size = 16),
        legend.position = "none",
        axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14)) +
  theme(axis.title = element_text())

### r.biocru_z_s1
r.biocru_z_s1.full <- effect("biocru_z_s1",r_model, partial.residuals =T)
plot(r.biocru_z_s1.full, smooth.residuals=F)

x.r.biocru_z_s1.fit <- unlist(r.biocru_z_s1.full$x.all)
x.r.biocru_z_s1 <- data.frame(lower = r.biocru_z_s1.full$lower, upper = r.biocru_z_s1.full$upper, fit = r.biocru_z_s1.full$fit, biocru_z_s1 = r.biocru_z_s1.full$x$biocru_z_s1)
xy.r.biocru_z_s1 <- data.frame(x = x.r.biocru_z_s1.fit, y = x.r.biocru_z_s1$fit[closest(trans(x.r.biocru_z_s1.fit), x.r.biocru_z_s1$biocru_z_s1)] + r.biocru_z_s1.full$residuals)
xy.r.biocru_z_s1$group <-phenotypic.group.mean.noNA.final$group


r.pc1 <- ggplot(x.r.biocru_z_s1, aes(x = biocru_z_s1, y = fit)) +
  geom_line(size = 1, col = "red") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_point(data = xy.r.biocru_z_s1, aes(x = x, y = y, color = group, shape = group),  size = 2) +
  scale_color_manual(values=c("#1b9e77", "#d95f02", "#66a61e", "#7570b3", "#e7298a", "#e6ab02"))+
  scale_shape_manual(values = c(16, 20, 6, 7, 18, 10))+
    ggtitle("(b)") +
  theme_bw() + 
  ylab("R-score (%)") +
  xlab("Climate PC1") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        text = element_text(size = 16),
        legend.position = "none",
        axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14)) +
  theme(axis.title = element_text())


### r.biocru_z_s3
r.biocru_z_s3.full <- effect("biocru_z_s3",r_model, partial.residuals =T)
plot(r.biocru_z_s3.full, smooth.residuals=F)

x.r.biocru_z_s3.fit <- unlist(r.biocru_z_s3.full$x.all)
x.r.biocru_z_s3 <- data.frame(lower = r.biocru_z_s3.full$lower, upper = r.biocru_z_s3.full$upper, fit = r.biocru_z_s3.full$fit, biocru_z_s3 = r.biocru_z_s3.full$x$biocru_z_s3)
xy.r.biocru_z_s3 <- data.frame(x = x.r.biocru_z_s3.fit, y = x.r.biocru_z_s3$fit[closest(trans(x.r.biocru_z_s3.fit), x.r.biocru_z_s3$biocru_z_s3)] + r.biocru_z_s3.full$residuals)
xy.r.biocru_z_s3$group <-phenotypic.group.mean.noNA.final$group

r.pc3 <- ggplot(x.r.biocru_z_s3, aes(x = biocru_z_s3, y = fit)) +
  geom_line(size = 1, col = "red") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_point(data = xy.r.biocru_z_s3, aes(x = x, y = y, color = group, shape = group),  size = 2) +
  scale_color_manual(values=c("#1b9e77", "#d95f02", "#66a61e", "#7570b3", "#e7298a", "#e6ab02"))+
  scale_shape_manual(values = c(16, 20, 6, 7, 18, 10))+
  ggtitle("(c)") +
  theme_bw() + 
  ylab("R-score (%)") +
  xlab("Climate PC3") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        text = element_text(size = 16),
        legend.position = "none",
        axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14)) +
  theme(axis.title = element_text())

## plot them together

library(cowplot)
first_col = plot_grid(c.cx, NULL, c.pc2, ncol = 3)
second_col = plot_grid(s.cx, s.pc1, s.pc2, ncol = 3)


final = plot_grid(first_col, second_col,  ncol = 1)
pdf("Fig.3 CS driver.pdf", useDingbats=FALSE, width=13, height=9)
final
dev.off()

## plot R as the supporting figure
pdf("Fig.S2 R driver.pdf", useDingbats=FALSE, width=13, height=4)
plot_grid(r.cx, r.pc1,r.pc3, ncol = 3)
dev.off()


####### run the multiple regression ####

summary(lm(C. ~ Cx +biocru_z_s1 + biocru_z_s2 + biocru_z_s3, phenotypic.group.mean.final))
summary(lm(S. ~ Cx +biocru_z_s1 + biocru_z_s2 + biocru_z_s3, phenotypic.group.mean.final))


### select the three haplotype M groups
target <- c("Europe", "NAmerica_in", "NAmerica_na")

haplotypeM <- phenotypic.group.mean.final %>%  filter(group %in% target) %>% droplevels()


library(lmerTest)
###C
c_lm_3groups <-lm(C. ~ group,haplotypeM)
summary(c_lm_3groups)
plot(c_lm_3groups)
anova(c_lm_3groups)  ##fixed factors

library(emmeans)
c.mean <- emmeans(c_lm_3groups, ~ group)
pairs(c.mean)
c.mean.table <- as.data.frame(summary(c.mean))[c('group','emmean', 'SE')]
is.data.frame(c.mean.table)
## S
s_lm_3groups <-lm(S. ~ group,  haplotypeM)
summary(s_lm_3groups)
plot(s_lm_3groups)
anova(s_lm_3groups) 

s.mean <- emmeans(s_lm_3groups, ~ group)
pairs(s.mean)
s.mean.table <- as.data.frame(summary(s.mean))[c('group','emmean', 'SE')]
is.data.frame(s.mean.table)


## R

r_lm_3groups <-lm(R. ~ group,  haplotypeM)
summary(r_lm_3groups)
anova(r_lm_3groups) 
plot(r_lm_3groups)

r.mean <- emmeans(r_lm_3groups, ~ group)
pairs(r.mean)
r.mean.table <- as.data.frame(summary(r.mean))[c('group','emmean', 'SE')]
is.data.frame(r.mean.table)
# ####plot the mean and SE does not use any more 01May2022 ######



library(cowplot)

s1 <- ggplot(s.mean.table, aes(x=group, y= emmean)) + 
  geom_point( size=3.5) +
  geom_errorbar(aes(ymin= emmean - SE, ymax= emmean + SE),
                width=.1) +
  #geom_hline(yintercept=0) +
  ggtitle("(c) S-score") +
  #scale_color_manual(values=c("#c2a5cf", "#a6dba0", "#008837")) +    ###add the dots
  # scale_shape_manual(values = c(15,16,17)) +
  theme_bw() +
  theme(plot.title = element_text(size=15, face="bold"),
        axis.title.x = element_text( size=14, face="bold"),
        axis.title.y = element_text( size=14, face="bold"),
        axis.text.x = element_text(face="bold", size=12),
        axis.text.y = element_text(face="bold", size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        #  legend.position = c(0.8, 0.2),legend.box = "horizontal",
        legend.text = element_text(size=15, face="bold"))

c1 <- ggplot(c.mean.table, aes(x=group, y= emmean)) + 
  geom_point( size=3.5) +
  geom_errorbar(aes(ymin= emmean - SE, ymax= emmean + SE),
                width=.1) +
  #geom_hline(yintercept=0) +
  ggtitle("(a) C-score") +
  #scale_color_manual(values=c("#c2a5cf", "#a6dba0", "#008837")) +    ###add the dots
  # scale_shape_manual(values = c(15,16,17)) +
  theme_bw() +
  theme(plot.title = element_text(size=15, face="bold"),
        axis.title.x = element_text( size=14, face="bold"),
        axis.title.y = element_text( size=14, face="bold"),
        axis.text.x = element_text(face="bold", size=12),
        axis.text.y = element_text(face="bold", size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        #  legend.position = c(0.8, 0.2),legend.box = "horizontal",
        legend.text = element_text(size=15, face="bold"))
r1 <- ggplot(r.mean.table, aes(x=group, y= emmean)) + 
  geom_point( size=3.5) +
  geom_errorbar(aes(ymin= emmean - SE, ymax= emmean + SE),
                width=.1) +
  #geom_hline(yintercept=0) +
  ggtitle("R-score") +
  #scale_color_manual(values=c("#c2a5cf", "#a6dba0", "#008837")) +    ###add the dots
  # scale_shape_manual(values = c(15,16,17)) +
  theme_bw() +
  theme(plot.title = element_text(size=15, face="bold"),
        axis.title.x = element_text( size=14, face="bold"),
        axis.title.y = element_text( size=14, face="bold"),
        axis.text.x = element_text(face="bold", size=12),
        axis.text.y = element_text(face="bold", size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        #  legend.position = c(0.8, 0.2),legend.box = "horizontal",
        legend.text = element_text(size=15, face="bold"))

pdf("c_S_3groups.pdf",  useDingbats=FALSE, width=5, height=9)

plot_grid(c1, s1, nrow=2, labels=c('(a)', '(b)'))
dev.off()

## violin plot of C, S and R   #####
library(plyr)
mean_c <- as.data.frame(ddply(haplotypeM, "group", summarise, grp.mean=mean(C.)))
head(mean_c)

median_c <- as.data.frame(ddply(haplotypeM, "group", summarise, grp.median=median(C.)))
head(median_c)

median_s <- as.data.frame(ddply(haplotypeM, "group", summarise, grp.median=median(S.)))
head(median_s)

median_r <- as.data.frame(ddply(haplotypeM, "group", summarise, grp.median=median(R.)))
head(median_r)


ggplot(haplotypeM, aes(x=C., color=group)) +
  geom_density()+
  geom_vline(data=mean_c, aes(xintercept=grp.mean, color=group),
             linetype="dashed") +
  scale_color_brewer(palette="Dark2") + theme_minimal()

mean_s <- ddply(haplotypeM, "group", summarise, grp.mean=mean(S.))
head(mean_s)

ggplot(haplotypeM, aes(x=S., color=group)) +
  geom_density()+
  geom_vline(data=mean_s, aes(xintercept=grp.mean, color=group),
             linetype="dashed") +
  scale_color_brewer(palette="Dark2") + theme_minimal()

mean_r <- ddply(haplotypeM, "group", summarise, grp.mean=mean(R.))
head(mean_r)


raincloud_theme <- theme(
  text = element_text(size = 10),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  axis.text = element_text(size = 14),
  axis.text.x = element_text(angle = 0, vjust = 0.5),
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 16),
  legend.position = c(0.9, 0.9),
  plot.title = element_text(lineheight = .8, face = "bold", size = 16),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
  axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"))

library(gghalves)
haplotypeM$group <- as.factor(haplotypeM$group)

plot_violin_C <- ggplot(haplotypeM, aes(group, C., fill=group))  +
  geom_violin() + 
  geom_boxplot(width=.1, outlier.shape = NA) +
  scale_fill_manual(values=c( "#d95f02",  "#7570b3", "#e7298a"))+
  # scale_fill_brewer(palette = "Dark2") +
  ylab("C-score (%)") + 
  stat_summary(fun=mean, geom="crossbar", width = 0.09,size = 0.3, color="white") +
  theme_bw() +
  raincloud_theme + theme(legend.position = "none")

plot_violin_S <- ggplot(haplotypeM, aes(group, S., fill=group))  +
  geom_violin() + #draw_quantiles= c(0.25, 0.75),stat = "half_ydensity", 
  geom_boxplot(width=.1, outlier.shape = NA) +
  scale_fill_manual(values=c( "#d95f02",  "#7570b3", "#e7298a"))+
  # scale_fill_brewer(palette = "Dark2") +
  ylab("S-score (%)") + 
  stat_summary(fun=mean, geom="crossbar", width = 0.09,size = 0.3, color="white") +
  theme_bw() +
  raincloud_theme + theme(legend.position = "none")

plot_violin_R <- ggplot(haplotypeM, aes(group, R., fill=group))  +
  geom_violin() +
  geom_boxplot(width=.1, outlier.shape = NA) +
  scale_fill_manual(values=c( "#d95f02",  "#7570b3", "#e7298a"))+
  # scale_fill_brewer(palette = "Dark2") +
  ylab("R-score (%)") + 
  stat_summary(fun=mean, geom="crossbar", width = 0.09,size = 0.3, color="white") +
  theme_bw() +
  raincloud_theme + theme(legend.position = "none")

pdf("Fig.3 c_S_3groups_violin.pdf",  useDingbats=FALSE, width=4.5, height=6)

plot_grid(plot_violin_C, plot_violin_S, nrow=2, labels=c('(a)', '(b)'))
dev.off()

pdf("Fig.S3 R_3groups_violin.pdf",  useDingbats=FALSE, width=4.5, height=4)
plot_violin_R
dev.off()

### end here   

library(rcompanion)
PT_S_3group = pairwisePermutationTest(S. ~ group, data = haplotypeM,
                                                    method="fdr")

PT_S_3group
## display letters
cldList(p.adjust ~ Comparison,
        data = PT_S_3group,
        threshold  = 0.05)


PT_C_3group = pairwisePermutationTest(C. ~ group, data = haplotypeM,
                                      method="fdr")

PT_C_3group
## display letters
cldList(p.adjust ~ Comparison,
        data = PT_C_3group,
        threshold  = 0.05)

PT_R_3group = pairwisePermutationTest(R. ~ group, data = haplotypeM,
                                      method="fdr")

PT_R_3group
## display letters
cldList(p.adjust ~ Comparison,
        data = PT_C_3group,
        threshold  = 0.05)

### non-parametric test as the data does not allow use ANOVA
kruskal.test(C. ~ group,
             data = haplotypeM)

kruskal.test(S. ~ group,
             data = haplotypeM)

library(rcompanion)
epsilonSquared(x = haplotypeM$group, g = haplotypeM$S.)
# output
epsilon.squared 

library(FSA)
dunnTest(S. ~ group,
         data = haplotypeM, method = "bonferroni")

pairwise.wilcox.test(haplotypeM$S., haplotypeM$group,
                     p.adjust.method = "BH")
pairwise.wilcox.test(haplotypeM$C., haplotypeM$group,
                     p.adjust.method = "BH")

kruskal.test(R. ~ group,
             data = haplotypeM)

ks.test(as.factor(haplotypeM_only$group),haplotypeM_only$C., 
        alternative = c("two.sided", "less", "greater"),
        exact = NULL, simulate.p.value = FALSE, B = 2000)

wilcox.test(haplotypeM_only$group,haplotypeM_only$S., 
       # alternative = c("two.sided", "less", "greater"),
        simulate.p.value = T, B = 20000)

### using library coin to run the imputation text
library(coin)
set.seed(20211206)
## c
sig_test_C_3group <- coin::oneway_test(C. ~ group, data = haplotypeM,
                                       distribution=approximate(nresample =100000))

print(sig_test_C_3group)


sig_test_S_3group <- coin::oneway_test(S. ~ group, data = haplotypeM,
                   distribution=approximate(nresample =100000))

print(sig_test_S_3group)


# using independence_test, can obtain the post hoc test results, suitable for multiple groups
sig_test_S_3group2 <- independence_test(S. ~ group, data = haplotypeM,
                                            distribution = approximate(nresample = 100000),
                                        xtrafo = mcp_trafo(group = "Tukey"))
pvalue(sig_test_S_3group2, method = "global")
pvalue(sig_test_S_3group2, method = "step-down", distribution = c("joint", "marginal"), 
       type = c("Bonferroni")) ##type = c("Sidak")) has the same results, but much slower

## R
sig_test_R_3group <- coin::oneway_test(R. ~ group, data = haplotypeM,
                                       distribution=approximate(nresample =100000))

print(sig_test_R_3group)



#### plot the collection points as supporting Fig. S1 ####
library(maps)
library(ggplot2)

world_map <- map_data("world")

#Creat a base plot with gpplot2
p <- ggplot() + coord_fixed() +
  xlab("") + ylab("")


#Add map to base plot
base_world_messy <- p + geom_polygon(data=world_map, aes(x=long, y=lat, group=group),
                                     colour="grey", fill= NA)

base_world_messy

base_world_messy +
  geom_point(data=phenotypic.group.mean.final,
             aes(x=phenotypic.group.mean.final$Longitude2, y= phenotypic.group.mean.final$Latitude2,
                 colour=group, shape = group), size=2) +  
  scale_color_manual(values=c("#1b9e77", "#d95f02", "#66a61e", "#7570b3", "#e7298a", "#e6ab02"))+
  scale_shape_manual(values = c(16, 20, 6, 7, 18, 10))+
  scale_y_continuous(limits=c(-58,84)) +
  # scale_x_continuous(limits = c(-180, 190)) +
  theme(plot.title = element_text(face=quote(bold)))  +
  guides(colour = guide_legend(override.aes = list(size=10)))+
  theme_bw() +
  labs( title="Phylogeographic groups",x="", y="") +
  theme(legend.position = c(0.1, 0.35))

pdf("fig.biogeographic_regions.pdf",  useDingbats=FALSE, width=12, height=8)

dev.off()


### add additional mixed model as one of the reviewers suggested
AIC(lm(C. ~ Cx +biocru_z_s1 + biocru_z_s2 + biocru_z_s3, phenotypic.group.mean.final))
plot(lm(S. ~ Cx +biocru_z_s1 + biocru_z_s2 + biocru_z_s3, phenotypic.group.mean.final))

library(lmerTest)
## C
C_mixedM <- glmer.nb(C. ~ Cx +biocru_z_s1 + biocru_z_s2 + biocru_z_s3 + (1|group), 
                     phenotypic.group.mean.final,verbose = T,
                    glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

summary(C_mixedM)

library(car)
Anova(C_mixedM)
plot(C_mixedM)
library(effects)
plot(allEffects(C_mixedM, typical=median))
## S
S_mixedM <- glmer.nb(S. ~ Cx +biocru_z_s1 + biocru_z_s2 + biocru_z_s3 + (1|group), 
                     phenotypic.group.mean.final,verbose = T,
                     glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1000000)))

summary(S_mixedM)
Anova(S_mixedM)
plot(S_mixedM)
plot(allEffects(S_mixedM, typical=median))

## R
R_mixedM <- glmer.nb(R. ~ Cx +biocru_z_s1 + biocru_z_s2 + biocru_z_s3 + (1|group), 
                     phenotypic.group.mean.final,verbose = T,
                     glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

summary(R_mixedM)
Anova(R_mixedM)
plot(R_mixedM)
plot(allEffects(R_mixedM, typical=median))
