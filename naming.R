require(lme4)
require(lmerTest)
require(ggplot2)
require(dplyr)
require(ggsignif)
require(ggpubr)
require(report)
require(multcomp)
require(tidyverse)
require(emmeans)
require(marginaleffects)
require(effects)
require(patchwork)
require(ggeffects)
require(broom.helpers)
require(ggstats)
require(gtsummary)
require(stringr)



## clean up data
data <- read.csv('/Users/agnesgao/Dropbox/UCD/Recording_data/l2_1_v2.csv', header=T)
df <- data.frame(data)
df <-filter(df, keypress == 'Z' | keypress =='M'|keypress =='X')

write.csv(df,"/Users/agnesgao/Dropbox/UCD/Recording_data/l2_1_clean_v2.csv", row.names = T)

##analysis

data <- read.csv('/Users/agnesgao/Desktop/pilot_NEW 2.csv', header=T)


data <- subset(data, data$RT>200)

# drop non words
data1 <- subset(data, bin!=555)
write.csv(data1,"/Users/agnesgao/Desktop/filtered.csv", row.names = T)

data1 <- read.csv('/Users/agnesgao/Desktop/filtered2.csv', header=T)

# drop wrong ld trials
data2 <- subset(data1, LD_ACC!=0)

#mark relatedness and pred results 1-relate 2-unrelate
data2$relate <- NA
data2$relate[which(data2$bin < 200)] <- 'unrelated'
data2$relate[which(data2$bin >= 200)] <- 'related'

#1-NAMED #2-UNNAMED
data2$naming <- NA 
data2$naming[which(data2$resp == 1)] <- 'Named Related'
data2$naming[which(data2$resp == 2)] <- 'Unnamed Related'
data2$naming[which(data2$resp == 3)] <- 'Unnamed Unrelated'
data2$naming[which(data2$resp == 4)] <- 'No Response'
data2$naming[which(data2$resp == 5)] <- 'Non-word'


data3 <- subset(data2, resp!=5)


##LD GRAPH
graph.data2 <- data3 %>% 
  group_by(naming) %>% summarise(
    mean = mean(RT, na.rm = T),
    sd = sd(RT, na.rm = T),
    N = length(RT)) %>% mutate(se = sd/sqrt(N))

ggplot(data = graph.data2, aes(x = naming, y = mean,fill=naming)) + 
  geom_bar(stat="identity") + 
  geom_errorbar(width = .3, aes(ymin = mean - se, ymax = mean + se)) + 
  #stat_compare_means(comparisons = my_comparisons, label.y = c(540, 560))+
  theme_classic() + 
  labs(y = "mean rt by condition", x = "condition", 
       title = "Main effects of condition")

##noresponse only
data_nr <- subset(data2, resp == 4)

##LD GRAPH
graph.data_nr <- data_nr %>% 
  group_by(relate) %>% summarise(
    mean = mean(RT, na.rm = T),
    sd = sd(RT, na.rm = T),
    N = length(RT)) %>% mutate(se = sd/sqrt(N))

ggplot(data = graph.data_nr, aes(x = relate, y = mean,fill=relate)) + 
  geom_bar(stat="identity") + 
  geom_errorbar(width = .3, aes(ymin = mean - se, ymax = mean + se)) + 
  #stat_compare_means(comparisons = my_comparisons, label.y = c(540, 560))+
  theme_classic() + 
  labs(y = "mean rt by condition", x = "condition", 
       title = "Main effects of condition")


## ONSET RT GRAPH



data3$onset<-as.numeric(data3$onset_og)

data_onset <- subset(data3, condition!="zero_response")

summary(data_onset$onset)
sum(data3$acc == "2")


graph.data3 <- data_onset %>% 
  group_by(condition) %>% summarise(
    mean = mean(onset, na.rm = T),
    sd = sd(onset, na.rm = T),
    N = length(onset)) %>% mutate(se = sd/sqrt(N))

ggplot(data = graph.data3, aes(x = condition, y = mean,fill=condition)) + 
  geom_bar(stat="identity") + 
  geom_errorbar(width = .3, aes(ymin = mean - se, ymax = mean + se)) + 
  #stat_compare_means(comparisons = my_comparisons, label.y = c(540, 560))+
  theme_classic() + 
  labs(y = "mean naming onset by condition", x = "condition", 
       title = "Main effects of condition")

data_onset_fil <- subset(data_onset, RT<=1200)

ggplot(data = data_onset_fil, mapping = aes(x = onset, y = RT)) +
  geom_point(alpha=0.1,aes(color = condition))+
  geom_smooth(aes(color = condition, fill = condition),method = "lm") +
  facet_wrap(facets =  vars(condition))



##lmer
#all condition

data4 <- subset(data3, resp!=4&5)
model <- lmer(log(RT) ~ naming +
                (1 | Sub)+
                #(0+naming_acc|Sub)+
                #(0+relatedness|Sub)+
                #(0+naming_acc:relatedness|Sub)+
                (1 | item),
                #(0+naming_acc|item)+
                #(0+relatedness|item),
                #(0+relatedness:naming_acc|item),
              REML=F, data=data4, control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 200000)))
summary(model)

drop1(model, .~.,test="Chisq")



m.em<-emmeans(model, list(pairwise ~ naming))
m.em
m.em.tran<-emmeans(model, list(pairwise ~ naming),type="response")
m.em.tran
m.rt.tran<-pairs(regrid(m.em),type = "response")



#plot main  - back log
plot1<-as.data.frame(m.em.tran)
plot1 <- plot1[-c(4,5,6),]

ggplot(data = plot1, aes(x = naming, y = response,color=naming)) + 
  geom_point(size=4) + 
  geom_errorbar(width = 0,aes(ymin = asymp.LCL, ymax = asymp.UCL)) + 
  scale_y_continuous(limits=c(300,500))+
  theme_classic(base_size = 25) + 
  theme(legend.position="none")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  labs(y = "Mean Response Time", x = "Condition")

#only look at relatedness
data_relate <- subset(data_onset, condition !="related_named")

model_r <- lmer(RT ~  condition +
                (condition | Sub),
              REML=F, data=data_relate, control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 200000)))


drop1(model_r, .~.,test="Chisq") #not significant
summary(model_r)



pairwise.t.test(data4$RT, data4$naming, p.adj = "none")

#no response - relatedness effect is significant
model2 <- lmer(log(RT) ~ relate+
                (relate | Sub),
              REML=F, data=data_nr, control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 200000)))
summary(model2)


drop1(model2, .~.,test="Chisq")


nr.em<-emmeans(model2, list(pairwise ~ relate))
nr.em
nr.em.tran<-emmeans(model2, list(pairwise ~ relate),type="response")
nr.em.tran
nr.rt.tran<-pairs(regrid(nr.em),type = "response")



#plot main  - back log
plot_nr<-as.data.frame(nr.em.tran)
plot_nr <- plot_nr[-c(3),-2]

ggplot(data = plot_nr, aes(x = relate, y = response,color=relate)) + 
  geom_point(size=4) + 
  geom_errorbar(width = 0, aes(ymin = lower.CL, ymax = upper.CL)) + 
  scale_y_continuous(limits=c(300,650))+
  theme_classic(base_size = 25) + 
  theme(legend.position="none")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  labs(y = "Mean Response Time", x = "Semantic Relatedness")




data4$naming<-as.factor(data4$naming)
data4 <- subset(data4, LEN2 < 11)

graph.data7 <- data4 %>% 
  group_by(LEN2,naming) %>% summarise(
    mean = mean(rt, na.rm = T),
    sd = sd(rt, na.rm = T),
    N = length(rt)) %>% mutate(se = sd/sqrt(N))


ggplot(data = graph.data7, aes(x = LEN2, y = mean, fill=naming, palette = "jco")) + 
  geom_point(aes(color = LEN2))+
  scale_x_continuous(limits=c(2, 10))+
  scale_y_continuous(limits=c(300,630))+
  geom_smooth(method ="lm")+
  geom_errorbar(width = .3, aes(ymin = mean - se, ymax = mean + se),position = position_dodge(), colour="sky blue") + 
  facet_wrap(~naming)




model3 <- lmer(onset ~ condition+
               (condition | Sub),
               REML=F, data=data_onset, control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 200000)))
summary(model3)
drop1(model3, .~.,test="Chisq")

pairwise.t.test(data_onset$onset, data_onset$condition, p.adj = "none")

### LM for RT as a function of response onset

data_relate_name <- subset(data_onset_fil, condition=="related_named")
data_relate_unname <- subset(data_onset_fil, condition=="related_unnamed")
data_unrelate_unname <- subset(data_onset_fil, condition=="unrelated_unnamed")

lm1 = lm(RT ~ onset, data=data_relate_name)
summary(lm1)
plot(RT ~ onset, data=data_relate_name)
abline(lm1)

lm2 = lm(RT ~ onset, data=data_relate_unname)
summary(lm2)
plot(RT ~ onset, data=data_relate_unname)
abline(lm2)


lm3 = lm(RT ~ onset, data=data_unrelate_unname)
summary(lm3)
plot(RT ~ onset, data=data_unrelate_unname)
abline(lm3)

data_relate <- subset(data3, relatedness!="unrelated")
data_incorrect <- subset(data3, naming_acc!="correct")



data_RP <- subset(data3, naming=='correctly named')


######lmer freq

data4<-subset(data3, naming != "No Response")
data4<-subset(data4, FREQ != "#N/A")
data4$naming<-as.factor(data4$naming)
data4$FREQ<-as.factor(data4$FREQ)

model_fr <- lmer(log(RT) ~  naming+FREQ+naming:FREQ+
                   (1 | Sub)+
                   (1 | item),
                 REML=F, data=data4, control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 200000)))
summary(model_fr)



fr.em<-emmeans(model_fr, list(pairwise ~ FREQ | naming,pairwise ~ naming | FREQ))
fr.em
fr.em.tran<-emmeans(model_fr, list(pairwise ~ FREQ | naming,pairwise ~ naming | FREQ),type="response")
fr.em.tran
fr.rt.tran<-pairs(regrid(fr.em),type = "response")
fr.rt.tran


#plot fr  - back log
plot_fr<-as.data.frame(fr.em.tran)
plot_fr <- plot_fr[-c(7:21),]


ggplot(data = plot_fr, aes(x = FREQ, y = response,color=naming)) + 
  geom_point(size=4) + 
  geom_errorbar(width = 0, aes(ymin = asymp.LCL, ymax = asymp.UCL)) + 
  scale_y_continuous(limits=c(300,550))+
  theme_classic(base_size = 25) + 
  theme(legend.position="none")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  labs(y = "Mean Response Time", x = "Frequency")+
  facet_wrap(~naming)

###conc
data4$CONC<-as.factor(data4$CONC)

model_co <- lmer(log(RT) ~  naming+CONC+naming:CONC+
                   (1 | Sub)+
                   (1 | item),
                 REML=F, data=data4, control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 200000)))
summary(model_co)



co.em<-emmeans(model_co, list(pairwise ~ CONC | naming,pairwise ~ naming | CONC))
co.em
co.em.tran<-emmeans(model_co, list(pairwise ~ CONC | naming,pairwise ~ naming | CONC),type="response")
co.em.tran
co.rt.tran<-pairs(regrid(co.em),type = "response")
co.rt.tran


#plot co  - back log
plot_co<-as.data.frame(co.em.tran)
plot_co <- plot_co[-c(7:21),]


ggplot(data = plot_co, aes(x = CONC, y = response,color=naming)) + 
  geom_point(size=4) + 
  geom_errorbar(width = 0, aes(ymin = asymp.LCL, ymax = asymp.UCL)) + 
  scale_y_continuous(limits=c(300,500))+
  theme_classic(base_size = 25) + 
  theme(legend.position="none")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  labs(y = "Mean Response Time", x = "Semantic Concreteness")+
  facet_wrap(~naming)




###Or
data4$ORTHO_N.1<-as.factor(data4$ORTHO_N.1)

model_or <- lmer(log(RT) ~  naming+ORTHO_N.1+naming:ORTHO_N.1+
                   (1 | Sub)+
                   (1 | item),
                 REML=F, data=data4, control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 200000)))
summary(model_or)



or.em<-emmeans(model_or, list(pairwise ~ ORTHO_N.1 | naming,pairwise ~ naming | ORTHO_N.1))
or.em
or.em.tran<-emmeans(model_or, list(pairwise ~ ORTHO_N.1 | naming,pairwise ~ naming | ORTHO_N.1),type="response")
or.em.tran
or.rt.tran<-pairs(regrid(or.em),type = "response")
or.rt.tran


#plot or  - back log
plot_or<-as.data.frame(or.em.tran)
plot_or <- plot_or[-c(7:21),]


ggplot(data = plot_or, aes(x = ORTHO_N.1, y = response,color=naming)) + 
  geom_point(size=4) + 
  geom_errorbar(width = 0, aes(ymin = asymp.LCL, ymax = asymp.UCL)) + 
  scale_y_continuous(limits=c(300,500))+
  theme_classic(base_size = 25) + 
  theme(legend.position="none")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  labs(y = "Mean Response Time", x = "Orthograhpic Neighborhood Size")+
  facet_wrap(~naming)


###pho
data4$PHONO_N.1<-as.factor(data4$PHONO_N.1)

model_pho <- lmer(log(RT) ~  naming+PHONO_N.1+naming:PHONO_N.1+
                   (1 | Sub)+
                   (1 | item),
                 REML=F, data=data4, control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 200000)))
summary(model_pho)



pho.em<-emmeans(model_pho, list(pairwise ~ PHONO_N.1 | naming,pairwise ~ naming | PHONO_N.1))
pho.em
pho.em.tran<-emmeans(model_pho, list(pairwise ~ PHONO_N.1 | naming,pairwise ~ naming | PHONO_N.1),type="response")
pho.em.tran
pho.rt.tran<-pairs(regrid(pho.em),type = "response")
pho.rt.tran


#plot pho  - back log
plot_pho<-as.data.frame(pho.em.tran)
plot_pho <- plot_pho[-c(7:21),]


ggplot(data = plot_pho, aes(x = PHONO_N.1, y = response,color=naming)) + 
  geom_point(size=4) + 
  geom_errorbar(width = 0, aes(ymin = asymp.LCL, ymax = asymp.UCL)) + 
  scale_y_continuous(limits=c(300,500))+
  theme_classic(base_size = 25) + 
  theme(legend.position="none")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  labs(y = "Mean Response Time", x = "Phonological Neighborhood Size")+
  facet_wrap(~naming)


###len
data4 <- subset(data4, LEN2<12)
model_len <- lmer(log(RT) ~  naming+LEN+naming:LEN+
                    (1 | Sub)+
                    (1 | item),
                  REML=F, data=data4, control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 200000)))
summary(model_len)


l.em.tran<-emmeans(model_len, list( pairwise ~ LEN|naming))
l.em.tran
l.rt.tran<-pairs(regrid(l.em.tran),type = "response")
l.rt.tran

#emmeans (model_len, specs = pairwise ~ LEN2, type = 'response')

#estimated slopes for length
emt<-emtrends(model_len, pairwise ~ naming, var = "LEN2",type="response")
test(emt)

emmip(model_len, naming ~ LEN2, cov.reduce = range,type="response", CIs = TRUE)+
  xlab("Word Length") +
  ylab("Linear Prediction")+
  geom_line(size=2)+
  theme(text = element_text(size = 20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")
  )

##plot ld rt by length
graph.data7 <- data4 %>% 
  group_by(LEN2,naming) %>% summarise(
    mean = mean(RT, na.rm = T),
    sd = sd(RT, na.rm = T),
    N = length(RT)) %>% mutate(se = sd/sqrt(N))
##NAMING RT BY LENGTH
ggplot(data = graph.data7, aes(x = LEN2, y = mean, color=naming)) + 
  geom_point(aes(shape=naming),size = 5)+
  #geom_smooth(method ="lm")+
  #geom_errorbar(width = .1, aes(ymin = mean - se, ymax = mean + se),position = position_dodge(), colour="black") + 
  labs(y = "Mean Response Time") + 
  theme_classic(base_size = 25)+
  scale_x_continuous(breaks = seq(2,11,1))+
  scale_y_continuous(breaks = seq(400,800,50))+
  labs(y = "Mean Response Time", x = "Word Length")+
  theme(legend.position = c(.3,.98))+
  theme(legend.title=element_blank())


data_RU <- subset(data3, naming=='incorrectly named related')
model_RU <- lmer(RT ~ FREQ  +
                   (1 | Sub)+
                   (1 | item),
                 REML=F, data=data_RU, control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 200000)))
summary(model_RU)
drop1(model_RU, .~.,test="Chisq")

data_UU <- subset(data3, naming=='incorrectly named unreltaed')
model_UU <- lmer(RT ~  FREQ +
                   (1 | Sub)+
                   (1 | item),
                 REML=F, data=data_UU, control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 200000)))
summary(model_UU)
drop1(model_UU, .~.,test="Chisq")


model_NR <- lmer(RT ~  FREQ +
                   (1 | Sub)+
                   (1 | item),
                 REML=F, data=data_nr, control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 200000)))
summary(model_NR)
drop1(model_NR, .~.,test="Chisq")

data3 <- subset(data3, FREQ!="#N/A")
graph.data_freq <- data3 %>% 
  group_by(FREQ,naming) %>% summarise(
    mean = mean(RT, na.rm = T),
    sd = sd(RT, na.rm = T),
    N = length(RT)) %>% mutate(se = sd/sqrt(N))

ggplot(data = graph.data_freq, aes(x = FREQ, y = mean, fill=naming)) + 
  geom_bar(stat="identity") + 
  geom_errorbar(width = .3, aes(ymin = mean - se, ymax = mean + se)) + 
  theme_classic() + 
  #geom_errorbar(width = .3, aes(ymin = mean - se, ymax = mean + se),position = position_dodge(), colour="black") + 
  labs(y = "mean LD rt", 
       title = "Frequency effect")+
  facet_wrap(~naming)

model_freqr <- lmer(RT ~  relatedness+FREQ+relatedness:FREQ+
                   (1 | Sub)+
                   (1 | item),
                 REML=F, data=data_incorrect, control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 200000)))
summary(model_freqr)

drop1(model_freqr, .~.,test="Chisq")


model_freqn <- lmer(RT ~  naming_acc+FREQ+naming_acc:FREQ+
                     (1 | Sub)+
                     (1 | item),
                   REML=F, data=data_relate, control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 200000)))
summary(model_freqn)

drop1(model_freqn, .~.,test="Chisq")

model_freqn <- lmer(RT ~  naming_acc+FREQ+naming_acc:FREQ+
                      (1 | Sub)+
                      (1 | item),
                    REML=F, data=data_relate, control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 200000)))
summary(model_freqn)

drop1(model_freqn, .~.,test="Chisq")

##########len

model_lenr <- lmer(RT ~  relatedness+LEN+relatedness:LEN+
                      (1 | Sub)+
                      (1 | item),
                    REML=F, data=data_incorrect, control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 200000)))
summary(model_lenr)

drop1(model_lenr, .~.,test="Chisq")


model_lenn <- lmer(RT ~  naming_acc+LEN+naming_acc:LEN+
                      (1 | Sub)+
                      (1 | item),
                    REML=F, data=data_relate, control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 200000)))
summary(model_lenn)

drop1(model_lenn, .~.,test="Chisq")


############phono
data3 <- subset(data3, PHONO_N.1!="#N/A")
graph.data_phono <- data3 %>% 
  group_by(PHONO_N.1,naming) %>% summarise(
    mean = mean(RT, na.rm = T),
    sd = sd(RT, na.rm = T),
    N = length(RT)) %>% mutate(se = sd/sqrt(N))

ggplot(data = graph.data_phono, aes(x = PHONO_N.1, y = mean, fill=naming)) + 
  geom_bar(stat="identity") + 
  geom_errorbar(width = .3, aes(ymin = mean - se, ymax = mean + se)) + 
  theme_classic() + 
  #geom_errorbar(width = .3, aes(ymin = mean - se, ymax = mean + se),position = position_dodge(), colour="black") + 
  labs(y = "mean LD rt", 
       title = "PHONO effect")+
  facet_wrap(~naming)


data_RP <- subset(data3, naming=='correctly named')
model_RP <- lmer(RT ~  PHONO_N.1 +
                   (1 | Sub)+
                   (1 | item),
                 REML=F, data=data_RP, control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 200000)))
summary(model_RP)
drop1(model_RP, .~.,test="Chisq")


data_RU <- subset(data3, naming=='incorrectly named related')
model_RU <- lmer(RT ~ PHONO_N.1  +
                   (1 | Sub)+
                   (1 | item),
                 REML=F, data=data_RU, control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 200000)))
summary(model_RU)
drop1(model_RU, .~.,test="Chisq")

data_UU <- subset(data3, naming=='incorrectly named unrelated')
model_UU <- lmer(RT ~  PHONO_N.1 +
                   (1 | Sub)+
                   (1 | item),
                 REML=F, data=data_UU, control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 200000)))
summary(model_UU)
drop1(model_UU, .~.,test="Chisq")


model_phonor <- lmer(RT ~  relatedness+PHONO_N.1+relatedness:PHONO_N.1+
                     (1 | Sub)+
                     (1 | item),
                   REML=F, data=data_incorrect, control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 200000)))
summary(model_phonor)

drop1(model_phonor, .~.,test="Chisq")


model_phono_nr <- lmer(RT ~  PHONO_N.1+
                     (1 | Sub)+
                     (1 | item),
                   REML=F, data=data_nr, control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 200000)))
summary(model_phono_nr)

drop1(model_phono_nr, .~.,test="Chisq")


##############CONC
data3 <- subset(data3, CONC!="#N/A")
graph.data_conc <- data3 %>% 
  group_by(CONC,naming) %>% summarise(
    mean = mean(RT, na.rm = T),
    sd = sd(RT, na.rm = T),
    N = length(RT)) %>% mutate(se = sd/sqrt(N))

ggplot(data = graph.data_conc, aes(x = CONC, y = mean, fill=naming)) + 
  geom_bar(stat="identity") + 
  geom_errorbar(width = .3, aes(ymin = mean - se, ymax = mean + se)) + 
  theme_classic() + 
  #geom_errorbar(width = .3, aes(ymin = mean - se, ymax = mean + se),position = position_dodge(), colour="black") + 
  labs(y = "mean LD rt", 
       title = "CONC effect")+
  facet_wrap(~naming)


data_RP <- subset(data3, naming=='correctly named')
model_RP <- lmer(RT ~  CONC +
                   (1 | Sub)+
                   (1 | item),
                 REML=F, data=data_RP, control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 200000)))
summary(model_RP)
drop1(model_RP, .~.,test="Chisq")


data_RU <- subset(data3, naming=='incorrectly named related')
model_RU <- lmer(RT ~ CONC  +
                   (1 | Sub)+
                   (1 | item),
                 REML=F, data=data_RU, control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 200000)))
summary(model_RU)
drop1(model_RU, .~.,test="Chisq")

data_UU <- subset(data3, naming=='incorrectly named unrelated')
model_UU <- lmer(RT ~  CONC +
                   (1 | Sub)+
                   (1 | item),
                 REML=F, data=data_UU, control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 200000)))
summary(model_UU)
drop1(model_UU, .~.,test="Chisq")

model_CONC_nr <- lmer(RT ~  CONC+
                         (1 | Sub)+
                         (1 | item),
                       REML=F, data=data_nr, control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 200000)))
summary(model_CONC_nr)

drop1(model_CONC_nr, .~.,test="Chisq")

model_CONCr <- lmer(RT ~  relatedness+CONC+relatedness:CONC+
                       (1 | Sub)+
                       (1 | item),
                     REML=F, data=data_incorrect, control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 200000)))
summary(model_CONCr)

drop1(model_CONCr, .~.,test="Chisq")


model_CONCi <- lmer(RT ~  naming_acc+CONC+naming_acc:CONC+
                      (1 | Sub)+
                      (1 | item),
                    REML=F, data=data_relate, control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 200000)))
summary(model_CONCi)

drop1(model_CONCi, .~.,test="Chisq")