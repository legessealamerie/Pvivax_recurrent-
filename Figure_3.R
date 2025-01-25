
##########################################################

############### Author: Dr. Jordache Ramjith, Radboud University Medical Center


#############################################################


# Load necessary libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(ggprism)
library(ggpubr)
library(survival)
library(MASS)
library(pammtools)
library(mgcv)
library(flextable)

dtfull <- read_excel("data_recurrent.xlsx")
dtfull <- dtfull %>% 
  mutate(id = factor(id), date = as.Date(date)) %>%
  group_by(id) %>% 
  mutate(firstdate = first(date),
         days=as.numeric(difftime(date, firstdate, units="days")),
         status = ifelse(type_infection %in% "symptomatic recruitment","Recruitment",ifelse(type_infection %in% c("asymptomatic recurrence","submicroscopic recurrence"),"asymptomatic",ifelse(type_infection %in% "symptomatic recurrence","symptomatic","follow-up")))) %>%
  select(-firstdate)

dtfull = dtfull[!is.na(dtfull$id),]

dtinf <- na.omit(dtfull %>% select(id, days, trt_given))
dtsub <- na.omit(dtfull %>% select(id,date, age, gender, sick_malaria_1_year_bl,sick_malaria_3_months_bl)) %>%
  mutate(age_cat1 = factor(ifelse(age<15,"<=15 years","Above 15 years")))


dtpar <- dtfull %>% select(id,days, status, pvcopy_ul, pvs25_copies_ul) %>%
  mutate(pvcopy_ul = ifelse(is.na(pvcopy_ul) | status=="follow-up", NA, pvcopy_ul), pvs25_copies_ul = ifelse(is.na(pvs25_copies_ul) | status=="follow-up", NA, pvs25_copies_ul))

dtpar0 <- dtpar %>% group_by(id) %>% mutate(keep = ifelse(sum(!is.na(pvcopy_ul))==1 ,1,0)) %>% filter(keep==1) %>% select(-keep) %>%  na.omit()
dtpar<- dtpar %>% group_by(id) %>% mutate(keep = ifelse(sum(!is.na(pvcopy_ul))==1,1,0)) %>% filter(keep==0) %>% select(-keep)  %>% mutate(keep = ifelse(!is.na(pvcopy_ul) | days==last(days),1,0)) %>% filter(keep==1) %>% select(-keep) %>% mutate(pvcopy_ul=lag(pvcopy_ul),pvs25_copies_ul=lag(pvs25_copies_ul))

dtpar = rbind(dtpar0,dtpar) %>% arrange(id, days)


# Load the dataset (update with the correct path)
dt <- dtfull %>% select(id, days, status)

dt = ((dt %>% left_join(dtinf, by=c("id","days"))) %>% left_join(dtsub, by="id")) %>% left_join(dtpar, by=c("id","days","status"))

dt$trt_given = ifelse(is.na(dt$trt_given) | dt$trt_given == "follow-up","None",dt$trt_given)

dt=data.frame(dt)


#For those ids who have only Recruitment and follow-ups we want two generate two rows of data. Row 1 where the strata is defined as recruitment->symptomatic with status=0 and row 2 where the strata is defined as recruitment->asymptomatic with status=0 and days in both should be last followup date. All other variables as the original row of recruitment.


#For those ids with recruitment and then a subsequent asymptomatic infection, we want to generate two rows of data. Row 1 where the strata is defined as recruitment->symptomatic with status=0 and days= days of the infection. Row 2 where the strata is defined as recruitment->asymptomatic with status=1 and days=days of the infection. All other variables as the original row of recruitment.

#For those ids with recruitment and then a subsequent symptomatic infection, we want to generate two rows of data. Row 1 where the strata is defined as recruitment->symptomatic with status=1 and days= days of the infection. Row 2 where the strata is defined as recruitment->asymptomatic with status=0 and days=days of the infection. All other variables as the original row of recruitment.

#In addition to the above rows, For those ids with an asymptomatic infection, we want to generate two rows of data for each asymptomatic infection. Row1 where the strata is defined as asymptomatic->symptomatic and the status is 1 if there is a symptomatic infection afterwards and 0 if not. The day should either be the day of the subsequent symptomatic infection (status=1) or if not (status=0) then the day of the subsequent asymptomatic infection or if not then the last day of follow up. All other variables as the original row of the asymptomatic infection. Row 2 where the strata is defined as asymptomatic->asymptomatic and the status is 1 if there is an asymptomatic infection afterwards and 0 if not. The day should either be the day of the subsequent asymptomatic infection (status=1) or if not (status=0) then the day of the subsequent symptomatic infection or if not then the last day of follow up. All other variables as the original row of the asymptomatic infection.


#In addition to the above rows, For those ids with a symptomatic infection, we want to generate two rows of data for each symptomatic infection. Row1 where the strata is defined as symptomatic->asymptomatic and the status is 1 if there is an asymptomatic infection afterwards and 0 if not. The day should either be the day of the subsequent asymptomatic infection (status=1) or if not (status=0) then the day of the subsequent symptomatic infection or if not then the last day of follow up. All other variables as the original row of the symptomatic infection. Row 2 where the strata is defined as symptomatic->symptomatic and the status is 1 if there is a symptomatic infection afterwards and 0 if not. The day should either be the day of the subsequent symptomatic infection (status=1) or if not (status=0) then the day of the subsequent asymptomatic infection or if not then the last day of follow up. All other variables as the original row of the symptomatic infection.


dt0 = dt %>% group_by(id) %>% mutate(anyevents = sum(I(status=="asymptomatic"|status=="symptomatic")), days = last(days), gap=days, rows=1:n()) %>% ungroup() %>%filter(anyevents==0, rows==1) %>% select(-rows,-anyevents)

dt0 = rbind(
  dt0 %>% mutate(strata = "symptomatic->symptomatic", delta=0),
  dt0 %>% mutate(strata = "symptomatic->asymptomatic", delta=0)
)

dt1 = dt %>% group_by(id) %>% mutate(rows=1:n(),keep = ifelse(rows==n()|status=="asymptomatic"|status=="symptomatic",1,0)) %>% ungroup() %>%filter(keep==1) %>% select(-rows,-keep) %>% group_by(id) %>% mutate(gap=days-lag(days)) %>% mutate(rows = sum(1:n())) %>% filter(rows>1) %>% select(-rows) %>% mutate(prevstat = lag(status), trt_given=lag(trt_given))  %>% ungroup() %>% filter(!is.na(gap)) 


dt1 = rbind(  
  dt1 %>% filter(status!="follow-up") %>% mutate(strata=paste0(prevstat,"->",status),delta=1),
  dt1 %>% filter(status!="follow-up") %>% mutate(strata=paste0(prevstat,"->",ifelse(status=="symptomatic","asymptomatic","symptomatic")),delta=0),
  dt1 %>% filter(status=="follow-up") %>% mutate(strata=paste0(prevstat,"->","symptomatic"),delta=0),
  dt1 %>% filter(status=="follow-up") %>% mutate(strata=paste0(prevstat,"->","asymptomatic"),delta=0)) %>% select(-prevstat)


dt2 = dt %>% group_by(id) %>% mutate(keep = ifelse(status=="Recruitment" | status=="asymptomatic" | status=="symptomatic",1,0)) %>%filter(keep==1) %>% mutate(rows=1:n(), gap=days-lag(days), trt_given=lag(trt_given)) %>% filter(rows==2) %>% select(-rows,-keep) 

dt2 = rbind(  
  dt2 %>%  mutate(strata=paste0("symptomatic","->",status),delta=1),
  dt2 %>%  mutate(strata=paste0("symptomatic","->",ifelse(status=="symptomatic","asymptomatic","symptomatic")),delta=0)
) 

dtmsm = rbind(dt0,dt1,dt2) %>% arrange(id, days) %>% select(-status)
dtmsm$gap = ifelse(dtmsm$gap==0,0.5,dtmsm$gap)


bsl <- read_excel("INDIE_baseline_data_factor.xlsx")

dtmsm = dtmsm %>% left_join(bsl, by="id")

xlsx::write.xlsx(data.frame(dtmsm), "survmsm_dt.xlsx", row.names=F)

dtmsm$indiv = as.factor(dtmsm$id)
dtmsm$id = 1:nrow(dtmsm)

# Collapse the dataset by removing duplicates (ignoring id and strata)
dtsurv <- dtmsm %>%
  # Transform 'strata' based on the conditions provided
  mutate(
    strata = case_when(
      strata %in% c("asymptomatic->asymptomatic", "asymptomatic->symptomatic") ~ "Asymptomatic recurrence",
      strata %in% c("symptomatic->asymptomatic", "symptomatic->symptomatic") ~ "Symptomatic recurrence",
      TRUE ~ strata  # Keep original value if it doesn't match any condition
    )
  ) %>%
  # Select all columns except 'id' and 'strata'
  select(-id) %>%
  # Group by all remaining columns except 'delta'
  group_by(across(-delta)) %>%
  # For duplicates, keep the row where delta is 1 (if it exists)
  summarise(delta = max(delta), .groups = "drop")%>%
  # Reorganize columns so that 'id', 'date', 'days', 'gap', 'delta' come first
  relocate(indiv, date, days, gap, delta,strata) %>%
  arrange(indiv, gap) %>%
  group_by(indiv) %>%
  mutate(tstart = lag(days, default=0))%>%
  ungroup() %>%
  mutate(id = 1:n()) %>%
relocate(id, indiv, date,tstart, days, gap, delta, strata) %>%
  mutate(strata=ifelse(tstart==0, "Symptomatic recruitment", strata))
  







table(dtsurv$strata,dtsurv$delta)

dtsurv$strata= factor(dtsurv$strata,
                      levels = c("Symptomatic recruitment", "Symptomatic recurrence", "Asymptomatic recurrence")
)

library(ggsurvfit)

km_fit <- survfit(Surv(gap, delta) ~ strata, data = dtsurv)

logrank_pval <- round(survdiff(Surv(gap, delta) ~ strata, data = dtsurv)$pvalue,4)
logrank_pval <- ifelse(logrank_pval<0.0001,"p<0.0001",paste0("p=",logrank_pval))



gg_km <- ggsurvfit(km_fit, lwd=1) +
  labs(
    x = "Time (days) between infections",
    y = "Recurrence-free probability"
  ) +
 # scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(labels = scales::number_format(scale = 100), breaks = seq(0, 1, by = 0.25)) +
  
  theme_prism()+
  theme(legend.position = "top")+
  scale_fill_manual(name="",values=c("#0072B2", "#92C5DE", "#D55E00"),labels=c("Symptomatic recruitment", "Symptomatic recurrence", "Asymptomatic recurrence"))+ 
  scale_color_manual(name="",values=c("#0072B2", "#92C5DE", "#D55E00"), labels=c("Symptomatic recruitment", "Symptomatic recurrence", "Asymptomatic recurrence"))+
  xlab("Time (days) between infections")+
  add_confidence_interval(type="ribbon")+
  guides(fill="none")+
  geom_text(aes(x=350,y=0.9,label=logrank_pval))+
  scale_x_continuous(breaks=seq(0,360,30))

gg_km


gg_km+
  add_risktable()


library(survminer)

# Define the main plot
gg_km <- ggsurvfit(km_fit, lwd=1) +
  labs(
    x = "gap times (days) between infections",
    y = "Recurrence-free probability"
  ) +
  scale_x_continuous(labels = scales::number_format(scale = 100), breaks = seq(0, 1, by = 0.25)) +
  theme_prism() +
  theme(legend.position = "top") +
  scale_fill_manual(name="", values=c("#0072B2", "#92C5DE", "#D55E00"), labels=c("Symptomatic recruitment", "Symptomatic recurrence", "Asymptomatic recurrence")) + 
  scale_color_manual(name="", values=c("#0072B2", "#92C5DE", "#D55E00"), labels=c("Symptomatic recruitment", "Symptomatic recurrence", "Asymptomatic recurrence")) +
  xlab("Time (days) between infections") +
  add_confidence_interval(type="ribbon") +
  guides(fill="none") +
  geom_text(aes(x=350, y=0.9, label=logrank_pval)) +
  scale_x_continuous(breaks=seq(0,360,30))

# Add the risk table with customized font size and removal of strata labels
gg_km <- gg_km + 
  add_risktable(fontsize = 3.5, risk.table.y.text = FALSE) +  # Adjust the fontsize as needed and remove y-axis text
  theme(
    axis.text.y = element_blank(),  # Remove the y-axis text in the main plot
    axis.title.y = element_blank()   # Remove the y-axis title in the main plot
  )

gg_km




