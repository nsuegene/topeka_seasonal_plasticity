## Analysis of seasonal cold tolerance (Topeka) °degree symbol for figures (just copy-paste)
getwd()
setwd("~/Dropbox/manuscripts/topeka/EcoEvo_revision/") #
# final version 

library(ggplot2)
library(plyr)
library(gridExtra)
library(reshape2)
library(scales)
library(car)
library(MASS)
library(lme4)
library(e1071)

fly_manuscript <- theme_classic() + theme(axis.line.x = element_line(color="black", size = .7),axis.line.y = element_line(color="black", size = .7),text=element_text(size=14))

#### degree day data for all collection dates across CCR and RCH experiments ####
# http://www.degreedays.net/#
dd <- as.data.frame(cbind(hdd25=c(5.1,39.2,139.9, 25.3,34.9,108.5, 29.7,26.1,27.9,75.6,126.3,200.1, 51.6,39.1,33.5,51.6),hdd18=c(0,7.4,62.8, 1.9,5,41.5, 0,3.3,1.7,23.3,48.7,106, 6.8,1.1,5.4,6.8),cdd25=c(75.3,30.9,1.1, 37.8,36.8,5, 19.4,31,26.6,15.4,3.9,0, 19.6,17.6,25.7,19.6),cdd18=c(168.3,97.1,22.2, 112.4,1.5,36.1, 87.8,106,98.2,61,24.5,4.7, 73.1,77.7,95.4,73.1)))
dd$time <- c(1,2,3, 4,5,6, 7,8,9,10,11,12, 13,14,15,16)
dd$label <- c("'12/Jul ","'12/Sep ","'12/Oct ","'13/Jul ","'13/Sep ","'13/Oct ","'14/Jun ","'14/Jul ","'14/Aug ","'14/Sep ","'14/Oct ", "'14/Nov ","'15/Jun ","'15/Jul ","'15/Sep ","'15/Oct ")
dd$year <- c("2012","2012","2012","2013","2013","2013","2014","2014","2014","2014","2014","2014","2015","2015","2015","2015")
#write.table(dd[,c(5:7,1:4)], file="degreedays_data.txt", quote=F, row.names=F, col.names=T)

skewness(dd$hdd25) #1.192
skewness(dd$hdd18) #1.588
skewness(dd$cdd25) #1.153
skewness(dd$cdd18) #0.132

dd.long <- melt(dd, id.vars=c("time","year","label"), variable.name="metric")

ggplot(subset(dd.long, metric=="cdd18"), aes(y=value, x=year, shape=year)) + geom_boxplot(show.legend=F, width=0.5) + geom_point(aes(fill=year), position = position_nudge(x=0.3), show.legend=F) + scale_shape_manual(values=c(21, 22, 24, 23)) + fly_manuscript + scale_fill_manual(values=c("#f7f7f7","#cccccc","#969696","#525252")) + ylab("Cumulative heat exposure over "*18~degree*C) + xlab("Collection year")

# figure 1
postscript(file=paste("fig_1_degree_days_bw",format(Sys.time(),"%Y%m%d"),"eps",sep="."),onefile=F, width=5, height=4)
ggplot(subset(dd.long, metric=="cdd18"), aes(y=value, x=year, shape=year)) + geom_boxplot(show.legend=F, width=0.5) + geom_point(aes(fill=year), position = position_nudge(x=0.3), show.legend=F) + scale_shape_manual(values=c(21, 22, 24, 23)) + fly_manuscript + scale_fill_manual(values=c("#f7f7f7","#cccccc","#969696","#525252")) + ylab("Cumulative heat exposure over "*18~degree*C) + xlab("Collection year")
dev.off()


#### graphical exploration of CCR data ####
# read processed data
all.ccr <- read.table("ccr_cleaned_data.20170420.txt", h=T)
str(all.ccr)  
all.ccr$year <- factor(all.ccr$year)
all.ccr$temp <- factor(all.ccr$temp, levels=c("25","18"))
all.ccr$group <- factor(all.ccr$group, levels=c("12jul","12sep","12oct","13jul","13sep","13oct","14jun","14jul","14aug","14sep"))
all.ccr$cat <- factor(all.ccr$cat, levels=c("low","mid","high"))

# line means taken for ease of exploration
all.filtered.fig <- ddply(all.ccr, .(line,sex,temp,group,year, cdd18, cat), summarise, mean=mean(time),sd=sd(time))
summary(all.filtered.fig)

# checking data, response against fixed predictors
ggplot(all.filtered.fig, aes(x=temp, y=mean, colour=line, group=line)) + geom_line(show.legend=F) + facet_grid(.~cat) + geom_point(shape=21,fill="white", show.legend = F) + fly_manuscript
ggplot(all.filtered.fig, aes(x=cat, y=mean, colour=line, group=line)) + geom_line(show.legend=F) + facet_grid(.~temp) + geom_point(shape=21,fill="white", show.legend = F) + fly_manuscript

# mean and variance are correlated in a quasipoisson pattern
groupvar <- with(all.ccr, tapply(time, list(gentemp), var))
groupmean <- with(all.ccr, tapply(time, list(gentemp), mean))
ggplot(data.frame(groupmean, groupvar), aes(x=groupmean, y=groupvar)) + geom_point() + geom_smooth(method="lm", formula=y~x-1) + fly_manuscript

# as expected for quasipoisson, log removes relationship
groupvar <- with(all.ccr, tapply(log(time), list(gentemp), var))
groupmean <- with(all.ccr, tapply(log(time), list(gentemp), mean))
ggplot(data.frame(groupmean, groupvar), aes(x=groupmean, y=groupvar)) + geom_point() + geom_smooth(method="lm") + fly_manuscript


#### analyze CCR data ####
## test two hypotheses:
# 1 - does cdd18 (seasonal weather) affect time25 (all years), controling for sex in line in year (random)
# 2 - how does cdd18 affect delta-time25-18 (2012 vs. 2014)
# quasipoisson
pql.ccr <- glmmPQL(time ~ cdd18 * temp, ~1|year/line/sex, family=quasipoisson(link="log"), data=all.ccr)

plot(ranef(pql.ccr))
plot(fitted(pql.ccr),residuals(pql.ccr),  lines(smooth.spline(fitted(pql.ccr), residuals(pql.ccr)),col="blue",lty=2,lwd=2))
qqnorm(residuals(pql.ccr))
qqline(residuals(pql.ccr))
plot(all.ccr$year, residuals(pql.ccr))
plot(all.ccr$group, residuals(pql.ccr))
plot(all.ccr$sex, residuals(pql.ccr))
plot(all.ccr$temp, residuals(pql.ccr))
plot(all.ccr$cat, residuals(pql.ccr))

summary(pql.ccr) 
# cooling degree days (reference 18 C) are representative of the cumulative hot weather or weather-related natural selection experienced by wild flies
# this represents selection pressure because tested flies are common garden offspring
# while the effect of temp is the degree of developmental plasticity
# and temp:cdd18 is the GxE component
# cdd18 adds to CCR wakeup time, 0.001 +/- <0.001 per unit (t(287) = 7.092, p < 0.001)
# temp18 decreases wakeup time, -0.163 +/- 0.011 (t(38589) = -14.963, p < 0.001)
# cdd18 and raising flies at 18 interact, so that per unit cdd18?, wakeup is shorter -0.001 +/- <0.001 (t(38589) = -8.830, p < 0.001)

# random variance between sexes is about 0.069 StdDev
# random variance between lines is about 0.102
# random variance between years is about 0.077
# residual variance is 1.751

fixef(pql.ccr)
slope.ccr <- fixef(pql.ccr)[2]
intercept.ccr <- fixef(pql.ccr)[1]
temp18 <- fixef(pql.ccr)[3]
cdd.temp18 <- fixef(pql.ccr)[4]

ranef(pql.ccr)["year"]
ccr.year2012 <- ranef(pql.ccr)$year[1,]
ccr.year2013 <- ranef(pql.ccr)$year[2,]
ccr.year2014 <- ranef(pql.ccr)$year[3,]

abline.data <- data.frame(int=c(intercept.ccr+ccr.year2012, intercept.ccr+ccr.year2014, intercept.ccr+ccr.year2012+temp18, intercept.ccr+ccr.year2014+temp18), slo=c(slope.ccr, slope.ccr, slope.ccr+cdd.temp18, slope.ccr+cdd.temp18), year=c("2012","2014","2012","2014"), temp=c("25", "25", "18", "18"))

temp_names <- c(
  `25` = "(a) Development at 25 °C",
  `18` = "(b) Development at 18 °C"
)

year_names <- c(
  `2012` = "(c) 2012",
  `2014` = "(d) 2014"
)

plot.ccr.1 <- ggplot(subset(all.filtered.fig, year!=2013), aes(x=cdd18, y=mean, fill=year, shape=year)) + geom_jitter(size=1, width=2.5) + fly_manuscript + facet_grid(.~temp, labeller=as_labeller(temp_names)) + scale_y_continuous(trans="log", breaks=seq(1,45,3)) + scale_shape_manual(values=c(21, 24)) + theme(strip.text = element_text(size=16), strip.background = element_blank()) + scale_fill_manual(values=c("#f7f7f7","#969696")) + geom_abline(aes(intercept=int, slope=slo, colour=year, linetype=year), abline.data) + xlab(expression("Heat above 18 °C experienced by line founder")) + ylab("Mean chill coma recovery")+ scale_colour_manual(values=c("#CCCCCC","#333333"))

plot.ccr.2 <- ggplot(subset(all.filtered.fig, year!=2013), aes(x=temp, y=mean, fill=year)) + geom_boxplot() +facet_grid(.~year, labeller=as_labeller(year_names)) + fly_manuscript + scale_x_discrete(labels=c("25 °C","18 °C")) + scale_y_continuous(trans="log", breaks=seq(1,45,3)) + theme(strip.text = element_text(size=16), strip.background = element_blank()) + scale_fill_manual(values=c("#f7f7f7","#969696")) + xlab("Developmental temperature") + ylab("Mean chill coma recovery")

# figure 2: glm model fit, panels by sex and year and dev temp
postscript(file=paste("fig_2_ccr_glm_model_fit_bw",format(Sys.time(),"%Y%m%d"),"eps",sep="."),onefile=F, width=8, height=5)
grid.arrange(plot.ccr.1, plot.ccr.2)
dev.off()

# a version that includes 2013
abline.data.1 <- data.frame(int=c(intercept.ccr+ccr.year2012 , intercept.ccr+ccr.year2013, intercept.ccr+ccr.year2014, intercept.ccr+ccr.year2012+temp18, intercept.ccr+ccr.year2013+temp18, intercept.ccr+ccr.year2014+temp18), slo=c(slope.ccr, slope.ccr, slope.ccr, slope.ccr+cdd.temp18, slope.ccr+cdd.temp18, slope.ccr+cdd.temp18), year=c("2012","2013","2014","2012","2013","2014"), temp=c("25", "25", "25", "18", "18", "18"))

year_names.1 <- c(
  `2012` = "(c) 2012",
  `2013` = "(d) 2013",
  `2014` = "(e) 2014"
)

plot.ccr.1a <- ggplot(all.filtered.fig, aes(x=cdd18, y=mean, fill=year, shape=year)) + geom_jitter(size=1, width=2.5) + fly_manuscript + facet_grid(.~temp, labeller=as_labeller(temp_names)) + scale_y_continuous(trans="log", breaks=seq(1,45,3)) + scale_shape_manual(values=c(21, 22, 24)) + theme(strip.text = element_text(size=16), strip.background = element_blank()) + scale_fill_manual(values=c("#f7f7f7","#cccccc","#969696")) + geom_abline(aes(intercept=int, slope=slo, colour=year, linetype=year), abline.data.1) + xlab("Heat above 18 C experienced by line founder") + ylab("Mean chill coma recovery") + scale_colour_manual(values=c("#CCCCCC","#ffffff","#333333"))

plot.ccr.2a <- ggplot(all.filtered.fig, aes(x=temp, y=mean, fill=year)) + geom_boxplot() +facet_grid(.~year, labeller=as_labeller(year_names.1)) + fly_manuscript + scale_x_discrete(labels=c("25 C","18 C")) + scale_y_continuous(trans="log", breaks=seq(1,45,3)) + theme(strip.text = element_text(size=16), strip.background = element_blank()) + scale_fill_manual(values=c("#f7f7f7","#cccccc","#969696")) + xlab("Developmental temperature") + ylab("Mean chill coma recovery")

grid.arrange(plot.ccr.1a, plot.ccr.2a)


# residual plots
pql.ccr.ranef.year <- data.frame(ranef(pql.ccr)["year"])
pql.ccr.ranef.year$year <- factor(rownames(pql.ccr.ranef.year), levels=c("2012", "2013", "2014"))
colnames(pql.ccr.ranef.year)[1] <- "intercept"

ccr.ranef <- ggplot(data.frame(ranef(pql.ccr)["line"]), aes(x=X.Intercept.)) + geom_histogram(binwidth=0.025) + geom_vline(data=pql.ccr.ranef.year, aes(xintercept=intercept, linetype=year), show.legend=T) + fly_manuscript + xlab("Random effect of line in year") + ylab("Frequency") + scale_linetype_manual(values=c(1,2,3,1,2,3,1,2,3,4), name="Year") + ggtitle("(a)")
ccr.ranef

all.ccr$residuals <- residuals(pql.ccr)
ccr.cdd <- ggplot(all.ccr, aes(x=cat, y=residuals)) + geom_boxplot() + fly_manuscript+ xlab("CDD18 category") + scale_x_discrete(labels=c("low", "mid", "high")) + ylab("Residuals") + ggtitle("(b)")
ccr.temp <- ggplot(all.ccr, aes(x=temp, y=residuals)) + geom_boxplot() + fly_manuscript + xlab("Developmental temperature") + scale_x_discrete(labels=c("25°","18°")) + ylab("Residuals") + ggtitle("(c)")
ccr.sex <- ggplot(all.ccr, aes(x=sex, y=residuals)) + geom_boxplot() + fly_manuscript + scale_x_discrete(labels=c("female","male")) + xlab("Sex") + ylab("Residuals")
ccr.group <- ggplot(all.ccr, aes(x=group, y=residuals)) + geom_boxplot() + fly_manuscript+ xlab("Founder collection year/month")+scale_x_discrete(labels=c("'12/Jul ","'12/Sep ","'12/Oct ","'14/Jul ","'14/Sep ","'14/Oct ","'14/Jun ","'14/Jul ","'14/Aug ","'14/Sep ")) + theme(axis.text.x=element_text(angle=45, hjust=1)) + ylab("Residuals") + ggtitle("(d)")


# figure s1: distribution of residuals across fixed and random effects
postscript(file=paste("fig_a1_ccr_ranef_resid",format(Sys.time(),"%Y%m%d"),"eps",sep="."),onefile=F, width=8, height=5)
grid.arrange(ccr.ranef, ccr.cdd, ccr.temp, ccr.group, ncol=2, nrow=2)
dev.off()


#### graphical exploration of RCH data #### 
# read processed rch data
all.rch <- read.table("rch_cleaned_data.20170420.txt", h=T)
str(all.rch)
all.rch$year <- factor(all.rch$year)
all.rch$time <- factor(all.rch$time, levels=c("2014july","2014aug","2014sept","2014oct","2014nov","2015june","2015july","2015sept","2015oct"))
all.rch$trt <- factor(all.rch$trt, levels=c("NON","ACC"))
all.rch$cat <- factor(all.rch$cat, levels=c("low","mid","high")) 

# line means taken for ease of exploration
rch.filtered.fig <- ddply(all.rch, .(cage,sex,trt,time,year, cdd18, cat), summarise, mean=mean(alive/(dead+1)),sd=sd(alive/(dead+1)))
summary(rch.filtered.fig)

# checking data, response against fixed predictors
ggplot(rch.filtered.fig, aes(x=trt, y=mean, colour=cage, group=cage)) + geom_line(show.legend=F) + facet_grid(.~cat) + geom_point(shape=21,fill="white", show.legend = F) + fly_manuscript
ggplot(rch.filtered.fig, aes(x=cat, y=mean, colour=cage, group=cage)) + geom_line(show.legend=F) + facet_grid(.~trt) + geom_point(shape=21,fill="white", show.legend = F) + fly_manuscript


groupvar <- with(all.rch, tapply(alive/(dead+1), list(gentrt), var))
groupmean <- with(all.rch, tapply(alive/(dead+1), list(gentrt), mean))
ggplot(data.frame(groupmean, groupvar), aes(x=groupmean, y=groupvar)) + geom_point() + geom_smooth(method="lm") + fly_manuscript
# mean and variance are correlated, quasipoisson pattern, so maybe analyze with quasibinomial?

groupvar <- with(all.rch, tapply(log(alive/(dead+1)), list(gentrt), var))
groupmean <- with(all.rch, tapply(log(alive/(dead+1)), list(gentrt), mean))
ggplot(data.frame(groupmean, groupvar), aes(x=groupmean, y=groupvar)) + geom_point() + geom_smooth(method="lm") + fly_manuscript
# log removes much but not all of the correlation between mean and variance

groupvar <- with(all.rch, tapply(surv, list(gentrt), var))
groupmean <- with(all.rch, tapply(surv, list(gentrt), mean))
ggplot(data.frame(groupmean, groupvar), aes(x=groupmean, y=groupvar)) + geom_point() + geom_smooth(method="lm") + fly_manuscript
# this is a much less workable relationship


## seeing if data fit binomial
qqp(all.rch$alive/(all.rch$dead+1), "norm") # no
qqp(c(all.rch$alive,all.rch$dead), "binom", size=round(mean(all.rch$total)), prob=mean(all.rch$surv))


#### analyze RCH data ####
# quasi families cannot be fit with glmer
## test two hypotheses
# 1 - does cdd/hdd (seasonal weather) affect c(alive,dead) (all years), controling for sex in cage in year (random)
# 2 - how does cdd/hdd affect delta-survivorship between no treatment and acclimation treatment 
glme.rch.full <- glmer(cbind(alive,dead) ~ cdd18*trt + (1|year/cage/sex), data=all.rch, family=binomial, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))) # does not converge
glme.rch <- glmer(cbind(alive,dead) ~ cdd18 + trt + (1|year/cage/sex), data=all.rch, family=binomial, control=glmerControl(optimizer="bobyqa"))
anova(glme.rch.full, glme.rch) #interaction term cdd18:trt is likely overfit, and model does not converge
relgrad <- with(glme.rch.full@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
rm(glme.rch.full)


plot(ranef(glme.rch))[1]
plot(fitted(glme.rch),residuals(glme.rch),  lines(smooth.spline(fitted(glme.rch), residuals(glme.rch)),col="blue",lty=2,lwd=2))
qqnorm(residuals(glme.rch))
qqline(residuals(glme.rch))
plot(all.rch$year, residuals(glme.rch))
plot(all.rch$time, residuals(glme.rch))
plot(all.rch$sex, residuals(glme.rch))
plot(all.rch$trt, residuals(glme.rch))
plot(all.rch$cat, residuals(glme.rch))

summary(glme.rch) 
glme.rch.se <- sqrt(diag(vcov(glme.rch))) 
glme.tab <- cbind(Est=fixef(glme.rch), LL=fixef(glme.rch)-1.96*glme.rch.se, UL=fixef(glme.rch)+1.96*glme.rch.se)
exp(glme.tab) # odds ratio, no longer in logit scale

# cooling degree days (reference 18 C) are representative of the cumulative hot weather or weather-related natural selection experienced by wild flies
# for hypothesis 2, this represent selection pressure because tested flies are common garden offspring
# while the effect of trt is the degree of short term plasticity
# and trt:cdd18 is the GxE component -> unable to test for rch
# cdd18 does not increase significantly survivorship (p = 0.161) (odds ratio 1.004 (0.998 - 1.010))
# trt increases survivorship, 2.001 +/- 0.045 (p < 0.001) (odds ratio 7.395 (6.764 - 8.085))
# random variance between sexes is about 0.169 (sd 0.411)
# random variance between cages is about 0.316 (sd 0.563)
# random variance between years is about 1.479 (sd 1.216)


print(glme.rch, corr=FALSE)
ranef(glme.rch)
rch.year2015 <- ranef(glme.rch)$year[2,]
rch.year2014 <- ranef(glme.rch)$year[1,]

fixef(glme.rch)
intercept.rch <- fixef(glme.rch)[1]
slope.rch <- fixef(glme.rch)[2]
trt <- fixef(glme.rch)[3]

abline.data.2 <- data.frame(int=c( intercept.rch+rch.year2014, intercept.rch+rch.year2014+trt, intercept.rch+rch.year2015, intercept.rch+rch.year2015+trt), year=c("2014","2014","2015","2015"), trt=c("NON","ACC","NON","ACC"))

trt_names <- c(
  `NON` = "(a) No acclimation",
  `ACC` = "(b) Short-term acclimation"
)

year_names.2 <- c(
  `2014` = "(c) 2014",
  `2015` = "(d) 2015"  
)


plot.rch.1 <- ggplot(rch.filtered.fig, aes(x=cdd18, y=mean, fill=year, shape=year)) + geom_jitter(size=1, width=2.5) + facet_grid(.~trt, labeller=as_labeller(trt_names)) + fly_manuscript + scale_y_continuous(trans="log", breaks=c(0,0.01,0.1,1,2,4,8,16)) + theme(strip.text = element_text(size=16), strip.background = element_blank()) + scale_fill_manual(values=c("#969696","#525252")) + scale_shape_manual(values=c(23, 24)) + xlab("Heat above 18 °C experienced by line founder") + ylab("Odds ratio of survival")+ scale_colour_manual(values=c("#CCCCCC","#333333")) + geom_abline(aes(intercept=int, slope=slope.rch, colour=year, linetype=year), abline.data.2) 

plot.rch.2 <- ggplot(rch.filtered.fig, aes(x=trt, y=mean, fill=year)) + geom_boxplot() +facet_grid(.~year, labeller=as_labeller(year_names.2)) + fly_manuscript + scale_x_discrete(labels=c("None","Short-term")) + scale_y_continuous(trans="log", breaks=c(0,0.01,0.1,1,2,4,8,16)) + theme(strip.text = element_text(size=16), strip.background = element_blank()) + scale_fill_manual(values=c("#969696","#525252")) + xlab("Acclimation treatment") + ylab("Odds ratio of survival")


# result figure 3: glm model fit, panels by sex and year and dev temp
postscript(file=paste("fig_3_rch_glm_model_fit_bw",format(Sys.time(),"%Y%m%d"),"eps",sep="."),onefile=F, width=8, height=5)
grid.arrange(plot.rch.1, plot.rch.2)
dev.off()


all.rch$residuals <- residuals(glme.rch)

glme.rch.ranef <- data.frame(ranef(glme.rch)["year"])
glme.rch.ranef$year <- factor(rownames(glme.rch.ranef), levels=c("2014","2015"))
colnames(glme.rch.ranef)[1] <- "intercept"
rch.ranef <- ggplot(data.frame(ranef(glme.rch)["cage:year"]), aes(x=X.Intercept.)) + geom_histogram(binwidth=0.2) + geom_vline(data=glme.rch.ranef, aes(xintercept=intercept, linetype=year, show.legend=T)) + fly_manuscript + xlab("Random effect of cage in year") + ylab("Frequency") + scale_linetype_manual(values=c(2,1,2,3,1,2,3,1,2,3,4), name="Year") + ggtitle("(a)")

rch.cdd <- ggplot(all.rch, aes(x=cat, y=residuals)) + geom_boxplot() + fly_manuscript + xlab("CDD18 category") + scale_x_discrete(labels=c("low", "mid", "high")) + ylab("Residuals") + ggtitle("(b)")
rch.trt <- ggplot(all.rch, aes(x=trt, y=residuals)) + geom_boxplot() + fly_manuscript + scale_x_discrete(labels=c("none","short-term")) +xlab("Acclimation treatment") + ylab("Residuals") + ggtitle("(c)")
rch.group <- ggplot(all.rch, aes(x=time, y=residuals)) + geom_boxplot() + fly_manuscript + xlab("Founder collection year/month")+scale_x_discrete(labels=c("'14/Jul ","'14/Aug ","'14/Sep ","'14/Oct ","'14/Nov ","'15/Jun ","'15/Jul ","'15/Sep ","'15/Oct ")) + theme(axis.text.x=element_text(angle=45, hjust=1)) + ylab("Residuals") + ggtitle("(d)")

# results figure s2: distribution of residuals across fixed and random effects
postscript(file=paste("fig_a2_rch_ranef_resid",format(Sys.time(),"%Y%m%d"),"eps",sep="."),onefile=F, width=8, height=5)
grid.arrange(rch.ranef, rch.cdd, rch.trt, rch.group, ncol=2, nrow=2)
dev.off()

## Essentially for survivorship data, there's a huge difference between 2014 and 2015 but not much between different collection times. Unfortunately cdd18 and year are confounded. The model is only able to pick up that treatment makes a difference but not much else. 


#### summary analyses of both experiments ####
all.filtered.fig2 <- ddply(subset(all.ccr, year!=2013), .(line,sex,temp,group,year, cdd18,cat), summarise, mean=mean(time),sd=sd(time))
temp <- dcast(all.filtered.fig2, year+group+line+sex+cdd18+cat ~ temp, value.var="mean")
colnames(temp)[7:8] <- c("ccr25","ccr18")
temp$plasticity <- temp$ccr25 - temp$ccr18
ccr.fig <- droplevels(temp)
ggplot(ccr.fig, aes(x=ccr25, y=plasticity, group=group)) + geom_point() + facet_grid(year~.) + scale_x_reverse()


rch.filtered.fig2 <- ddply(all.rch, .(cage,sex,trt,time,year, cdd18,cat), summarise, mean=mean(surv),sd=sd(surv))
temp <- dcast(rch.filtered.fig2, year+time+cage+sex+cdd18+cat ~ trt, value.var="mean")
temp$acclimation <- temp$ACC - temp$NON
rch.fig <- droplevels(temp)
ggplot(rch.fig, aes(x=NON, y=acclimation, group=time)) + geom_point() + facet_grid(year~.)


summary(ccr.fig$cdd18)
table(ccr.fig$cdd18)
summary(rch.fig$cdd18)
table(rch.fig$cdd18)

#cut(cdd18, c(4,59,114,169)) # this can't work because there's only one time in low, and only one time in high. All else are mid
# just using previous grouping: cat

cat_names <- c(
  `low` = "(a) Low heat",
  `mid` = "(b) Mid heat",
  `high` = "(c) High heat"  
)

cat_names2 <- c(
  `low` = "(d) Low heat",
  `mid` = "(e) Mid heat",
  `high` = "(f) High heat"  
)


cat_text1 <- data.frame(cat=c("low","mid","high"),r=c(round(cor(subset(ccr.fig, cat=="low")[,c(7,9)])[2,1], digit=2),round(cor(subset(ccr.fig, cat=="mid")[,c(7,9)])[2,1], digit=2),round(cor(subset(ccr.fig, cat=="high")[,c(7,9)])[2,1], digit=2)))
cat_text1$lab <- paste0("r==",cat_text1$r)

ccr_summary <- ggplot(ccr.fig, aes(x=ccr25, y=plasticity, group=group, fill=year, shape=year)) + geom_point() + scale_x_reverse() + facet_grid(.~cat, labeller=as_labeller(cat_names)) + fly_manuscript + scale_fill_manual(values=c("#f7f7f7","#969696")) + scale_shape_manual(values=c(21, 24)) + ylab("Developmental plasticity") + xlab("Basal cold tolerance (chill coma recovery)") + geom_text(aes(label=lab, x=20, y=20, size=14), parse=T, inherit.aes=F, data=cat_text1, show.legend=F) + theme(strip.text = element_text(size=16), strip.background = element_blank()) 


cat_text2 <- data.frame(cat=c("low","mid","high"),r=c(round(cor(subset(rch.fig, cat=="low")[,c(7,9)])[2,1], digit=2),round(cor(subset(rch.fig, cat=="mid")[,c(7,9)])[2,1], digit=2),round(cor(subset(rch.fig, cat=="high")[,c(7,9)])[2,1], digit=2)))
cat_text2$lab <- paste0("r==",cat_text2$r)

rch_summary <- ggplot(rch.fig, aes(x=NON, y=acclimation, group=time, fill=year, shape=year)) + geom_point() + facet_grid(.~cat, labeller=as_labeller(cat_names2)) + fly_manuscript + scale_fill_manual(values=c("#969696","#525252")) + scale_shape_manual(values=c(24,23)) + ylab("Acclimation score") + xlab("Basal cold tolerance (cold stress surval)") + geom_text(aes(label=lab, x=.5, y=.8, size=14), parse=T, inherit.aes=F, data=cat_text2, show.legend=F) + theme(strip.text = element_text(size=16), strip.background = element_blank()) + scale_x_continuous(labels = c(0,.25,.5,.75,1))

# figure 4 - summary of plasticity and basal tolerance
postscript(file=paste("fig_4_plasticity_ccr_rch",format(Sys.time(),"%Y%m%d"),"eps",sep="."),onefile=F, width=8, height=6.5)
grid.arrange(ccr_summary, rch_summary, nrow=2)
dev.off()


## formal test of different R across cdd18 categories
# checking to make sure factors are specified properly
str(ccr.fig)
str(rch.fig)

# checking distributions of response variables
hist(ccr.fig$plasticity)
hist(rch.fig$acclimation)

# fitting linear models and checking significance of factors
ccr.aov <- lm(plasticity ~ cat*ccr25, data=ccr.fig)
anova(ccr.aov)
# the interaction is not significant F=0.746, P=0.475
summary(ccr.aov)$fstatistic

rch.aov <- lm(acclimation ~ cat*NON, data=rch.fig)
anova(rch.aov)
# the interaction is significant F=3.746, P=0.027
summary(rch.aov)$fstatistic

sessionInfo()


