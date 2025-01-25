
##################################################

################                  ##################

####################################################

library(ggpubr)
library(readr)
library(dplyr)
library(readr)
library(ggpmisc)  
library(tidyverse)
library(ggplot2)
library(ggprism)


dat <- read_csv("data/data_recurrent_all.csv")
dat = dat[!is.na(dat$date_visit),]
dat = dat[!is.na(dat$date_bl),]
dat$pid = stringr::str_sub(dat$id, start= -3)
dat = dat %>% group_by(pid) %>% mutate(date_visit=as.Date(date_visit,tryFormats = c("%m/%d/%Y")),date_bl=as.Date(date_bl,tryFormats = c("%m/%d/%Y")),duration=as.numeric(date_visit-date_bl),maxd = max(duration,na.rm=T))
sort(dat$maxd[dat$duration==0])
dat$pid=factor(dat$pid, levels=unique(dat$pid[order(dat$maxd,decreasing = TRUE)]), ordered=TRUE)

dat$pvcopies_ul=dat$pvcopy_ul
dat$pvs25_copiesul=dat$pvs25_copies_ul

##### classification of recurrence
dat=dat%>%mutate(rec_ins=factor(case_when(infcat_status=="Recruitment"~0,recurrenceinstance==1~1,recurrenceinstance==2~2,recurrenceinstance%in%c(3:4)~3)))
dat$rec_ins = factor(dat$rec_ins, labels=c("Baseline","First","Second","Third & above"))
####### only recurrence instance
dat=dat%>%mutate(rec_eps=factor(case_when(infcat_status=="Recruitment"~0,recurrenceinstance==1~1,recurrenceinstance==2~2,recurrenceinstance==3~3,recurrenceinstance %in%c(4:5)~4)))
dat$rec_eps = factor(dat$rec_eps, labels=c("Baseline","First","Second","Third","Fourth"))

table(dat$clinical_status)
symp=dat%>%filter(clinical_status%in%c("Recruitment","symptomatic"))
dim(symp)

asymp=dat%>%filter(clinical_status%in%c("Recruitment","asymptomatic"))
dim(asymp)

######## summary statistics for symptomatic individual by recurrence instance 
group_by(symp, rec_ins) %>%
  summarise(
    count = n(),
    median = median(pvcopies_ul, na.rm = TRUE),q1 = quantile(pvcopies_ul, 0.25,na.rm = TRUE),
    q3 = quantile(pvcopies_ul, 0.75,na.rm = TRUE),IQR = IQR(pvcopies_ul, na.rm = TRUE))
### median and IQR for pvs25 copies by recurrence instance 
group_by(symp, rec_ins) %>%
  summarise(
    count = n(),
    median = median(pvs25_copiesul, na.rm = TRUE),q1 = quantile(pvs25_copiesul, 0.25,na.rm = TRUE),
    q3 = quantile(pvs25_copiesul, 0.75,na.rm = TRUE),IQR = IQR(pvs25_copiesul, na.rm = TRUE))


######## summary statistics for asymptomatic individual by recurrence instance
group_by(asymp, rec_ins) %>%
  summarise(
    count = n(),
    median = median(pvcopies_ul, na.rm = TRUE),q1 = quantile(pvcopies_ul, 0.25,na.rm = TRUE),
    q3 = quantile(pvcopies_ul, 0.75,na.rm = TRUE),IQR = IQR(pvcopies_ul, na.rm = TRUE))
### median and IQR for pvs25 copies by recurrence instance  

group_by(asymp, rec_ins) %>%
  summarise(
    count = n(),
    median = median(pvs25_copiesul, na.rm = TRUE),q1 = quantile(pvs25_copiesul, 0.25,na.rm = TRUE),
    q3 = quantile(pvs25_copiesul, 0.75,na.rm = TRUE),IQR = IQR(pvs25_copiesul, na.rm = TRUE))


cols=c("#D55E00","#0072B2","#92C5DE","#009E73")

##### creating new variable recurrence instance by combing persistant individuals over follow up period to see distribution of the data point

dat1=dat%>%filter(clinical_status%in%c("asymptomatic","Recruitment","symptomatic"))

library(scales)
options(scipen=999)

pv18s=ggplot(data=filter(dat1,total_recurrence!=0,rec_ins!="NA"),aes(x = factor(rec_ins), y = pvcopies_ul,color=clinical_status))+ 
  
  geom_jitter(alpha=0.2, size=2, position=position_jitterdodge(jitter.width=0.2, dodge.width=0.75)) +
  geom_boxplot(aes(fill=clinical_status), alpha=0.5, position=position_dodge(width=0.75)) +
  theme_prism() +
  
  labs(x = "Recurrent episode", y = "Pv18S copies/μL") +
  #ggtitle("") +
  scale_y_log10(label =comma_format(big.mark = ""))+
  theme_classic() +theme(axis.text.x = element_text(face = "bold"),legend.text = element_text(size = 10, family = "Times New Roman"))+
  
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),limits=c(1,100000000),
                labels = trans_format("log10", math_format(10^.x)))+
  theme(legend.position = "Top" )+
  scale_color_manual(values=cols,labels=c("Asymp","Baseline", "Symp"))+
  scale_fill_manual(values=cols,labels=c("Asymp","Baseline", "Symp"))+
  
  guides(fill=F)+theme_prism()

##pvs25 copies
pvs25=ggplot(data=filter(dat1,total_recurrence!=0,rec_ins!="NA"),aes(x = factor(rec_ins), y =pvs25_copiesul ,color=clinical_status),size=2)+ 
  geom_jitter(alpha=0.2, size=2, position=position_jitterdodge(jitter.width=0.2, dodge.width=0.75)) +
  geom_boxplot(aes(fill=clinical_status), alpha=0.5, position=position_dodge(width=0.75)) +
  theme_prism() +
  labs(x = "Recurrent episode", y = "Pvs25 transcripts copy/µL") +
  scale_y_log10(label =comma_format(big.mark = ""))+
  theme_classic() +theme(axis.text.x = element_text(face = "bold"),legend.text = element_text(size = 10, family = "Times New Roman"))+
  
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),limits=c(1,100000000),
                labels = trans_format("log10", math_format(10^.x)))+
  theme(legend.position = "Top" )+
  scale_color_manual(values=cols,labels=c("Asymp","Baseline", "Symp"))+
  scale_fill_manual(values=cols,labels=c("Asymp","Baseline", "Symp"))+
  
  guides(fill=F)+theme_prism()

pvs25


# 

dat=dat%>%mutate(clinical_status1=factor(case_when(clinical_status=="Recruitment"~0,clinical_status=="symptomatic"~1,clinical_status=="asymptomatic"~2,clinical_status=="Asym_persisting"~3,clinical_status=="Symp_persisting"~3)))
table(dat$clinical_status1)

dat$clinical_status1 = factor(dat$clinical_status1, labels=c("Recruitment","Symptomatic recurrence","Asymptomatic recurrence","Persisting"))

table(dat$clinical_status1)

dat=dat%>%mutate(clinical_status2=factor(case_when(clinical_status=="Recruitment"~0,clinical_status=="symptomatic"~1,clinical_status=="asymptomatic"~2)))
table(dat$clinical_status2)

dat$clinical_status2 = factor(dat$clinical_status2, labels=c("Recruitment","Symptomatic recurrence","Asymptomatic recurrence"))




cols=c("#0072B2","#92C5DE","#D55E00","#009E73")



means <- dat %>%filter(clinical_status2!="NA") %>%group_by(clinical_status2) %>% 
  summarise(mean_pv18s = 10^mean(log10(pvcopies_ul+0.1), na.rm=TRUE),
            mean_pvs25 = 10^mean(log10(pvs25_copiesul+0.1), na.rm=TRUE))


p1 <- ggplot() +
  geom_density(data = filter(dat,clinical_status2!="NA"), aes(x = pvcopies_ul, fill = clinical_status2, color = clinical_status2),
               alpha = 0.2, linewidath = 0.8) +
  geom_vline(data = means, aes(xintercept = mean_pv18s, color = clinical_status2),
             linewidath = 1.3) +theme(axis.text.x = element_text(face = "bold"),legend.text = element_text(size = 10, family = "Times New Roman"))+
  scale_x_log10(label =comma_format(big.mark = ""))+
  scale_fill_manual("Clinical status ", values = cols) +
  scale_color_manual("Clinical status", values = cols) +
  labs(x = "Pv18S copies/uL", y = "Density") +
  ggtitle("") +theme_classic()+ylim(0,0.6)+
  guides(fill = FALSE)+theme(legend.position = "bottom")+theme_prism()
### pvs25
p2 <- ggplot() +
  geom_density(data = filter(dat,clinical_status2!="NA"), aes(x = pvs25_copiesul, fill = clinical_status2, color = clinical_status2),
               alpha = 0.2, linewidath = 0.8) +
  geom_vline(data = means, aes(xintercept = mean_pvs25, color = clinical_status2),
             linewidath = 1.3) +theme(axis.text.x = element_text(face = "bold"),legend.text = element_text(size = 10, family = "Times New Roman"))+
  scale_x_log10(label =comma_format(big.mark = ""))+
  scale_fill_manual("clinical status ", values = cols) +
  scale_color_manual("clinical status  ", values = cols) +
  labs(x = "Pvs25 transcripts copy/µL", y = "Density") +
  ggtitle("") +theme_classic()+scale_y_continuous(breaks=seq(0,0.5,0.1))+ylim(0,0.5)+
  guides(fill = FALSE)+theme(legend.position = "bottom")+theme_prism()

#### mean and 95% CI for parasite and gametosite density 

means <-dat %>% group_by(clinical_status2) %>% 
  summarise(n=n(),mean_pv18s = 10^mean(log10(pv18s_copy_ul+0.1), na.rm=TRUE),
            pv18s_lower = 10^t.test(log10(pv18s_copy_ul+0.1))$conf.int[1],
            pv18s_upper = 10^t.test(log10(pv18s_copy_ul+0.1))$conf.int[2])
## for gam dens
means <-dat %>% group_by(clinical_status2) %>% 
  summarise(n=n(),mean_pv25s = 10^mean(log10(pvs25_copiesul+0.1), na.rm=TRUE),
            pv25s_lower = 10^t.test(log10(pvs25_copiesul+0.1))$conf.int[1],
            pv25s_upper = 10^t.test(log10(pvs25_copiesul+0.1))$conf.int[2])


ggarrange(p1, p2,pv18s,pvs25, common.legend = TRUE,labels = c("(A)", "(C)","(B)","(D)"), legend = "bottom",nrow=2,ncol=2)
ggsave("plots/Figure_4AD_updated.tiff", device="tiff",dpi=600, height=10, width=12, units="in", compression = "lzw")

ggsave("plots/Figure_4AD_updated.png",device="png", width=11, height=8, dpi=600)


#### corelation  plot based on clinical status 
cols=c("#0072B2","#92C5DE","#D55E00","#009E73")


dat$parl=log10(dat$pvcopies_ul+0.1)
dat$gaml=log10(dat$pvs25_copiesul+0.1)

table(dat$clinical_status1)
dat_cl=dat%>%filter(clinical_status%in%c("Recruitment","asymptomatic","symptomatic","persisting"))


###########
###################

dat$gam=dat$pvs25_copiesul
dat$par=dat$pvcopies_ul

cor_clin=dat%>%filter(infcat_status!="NA",clinical_status1!="NA")%>%
  ggplot()+(aes(pvcopies_ul+0.1,pvs25_copiesul+0.1,color=clinical_status1)) +
  geom_smooth(method = gam,se=T,color="black")+
  geom_point(aes(x = par, y = gam, color = clinical_status1)) +
  stat_correlation(method = "pearson",mapping = use_label(c("r", "P", "n")),size=4)+
  scale_y_log10()+
  scale_x_log10()+
  labs(x ="Pv18S copies/μL", y ="Pvs25 transcripts copy/µL") +xlim(0,8)+
  scale_color_manual("Clinical status", values = cols) +
  theme_prism()+guides(shape=F)+
  scale_y_log10(limits=c(-0.1, 10000000), breaks=c(1, 10, 100, 1000, 10000, 100000, 1000000,10000000),labels = trans_format("log10", math_format(10^.x))) +
  
  scale_x_log10(limits=c(-0.1, 10000000), breaks=c(1, 10, 100, 1000, 10000, 100000, 1000000,10000000),labels = trans_format("log10", math_format(10^.x))) +
  
    theme(legend.position = "bottom")




