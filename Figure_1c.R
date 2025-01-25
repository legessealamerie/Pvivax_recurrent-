

#### package
library(ggpubr)
library(prismatic)
library(ggprism)
library(reshape2)
library(tidyverse)
library(scales)
library(ggplot2)
library(dplyr)

###################

data=read_csv("data/Figure_1c.csv")


data$Year <- as.factor(data$Year)

data$Year <- as.factor(data$Year)
data$months <- factor(data$months, 
                           levels = c('January', 'February', 'March', 'April', 'May', 'June', 
                                      'July', 'August', 'September', 'October', 'November', 'December'))

# Create a custom date for x-axis
data$date <- as.Date(paste(data$Year, data$months, "01", sep="-"), format="%Y-%B-%d")

# 
p1 <- ggplot(data, aes(x = date)) +
  geom_line(aes(y = Pf, color = 'Pf')) +
  geom_point(aes(y = Pf, color = 'Pf')) +
  geom_line(aes(y = Pv, color = 'Pv')) +
  geom_point(aes(y = Pv, color = 'Pv')) +
  geom_line(aes(y = Recruitment, color = 'Recruitment')) +
  geom_point(aes(y = Recruitment, color = 'Recruitment')) +
  geom_bar(aes(y = Rainfall * (700/130), fill = 'Rainfall'), stat = 'identity', alpha = 0.5) +
  scale_y_continuous(
    name = 'Total cases',
    limits = c(0, 700),breaks = c(0,200,300,400,500,600,700,100),
    sec.axis = sec_axis(~. * (130/700), name = 'Rainfall', breaks = seq(0, 130, 20))
  ) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b", limits = c(as.Date("2020-01-01"), as.Date("2022-12-28"))) +
  labs(x = '', color = 'Legend', fill = 'Legend') +
  theme_classic() +  scale_color_manual("",values=c("#FFA07A", "#66CDAA",  "grey80","skyblue","red3","yellow","grey67"))+
  scale_fill_manual("",values=c("skyblue","#FFA07A", "#66CDAA",  "grey80","red3","yellow","grey67"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position = "bottom")+
  annotate("text", x = as.Date("2020-06-15"), y = -50, label = "2020", size = 4, vjust = 1.5) +
  annotate("text", x = as.Date("2021-06-15"), y = -50, label = "2021", size = 4, vjust = 1.5) +
  annotate("text", x = as.Date("2022-06-15"), y = -50, label = "2022", size = 4, vjust = 1.5)


print(p1)



