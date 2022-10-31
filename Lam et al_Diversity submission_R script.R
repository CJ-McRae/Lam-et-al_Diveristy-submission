###Diversity Manuscript Submission: R script###


#Long-term consistent reproduction and completion of a brooding coral life cycle through ex situ culture
#Lam et al. 
#Submission date to Diversity: November 2022 


#Pocillopora acuta planulae were collected and counted each day at ~9am 
#Experiment duration: Lunar Feb 2021- lunar January 2022 [~ Gregorian calendar dates: March 2021 - February 2022]
#FTS = Flow-through system
#RAS = Recirculating aquaculture system

#R citations:
#citation()
#citation(package = "tidyverse")
#citation(package = "car")
#citation(package = "lme4")
#citation(package = "lmerTest")
#citation(package = "emmeans")
#citation(package = "AER")
#citation(package = "glmmTMB")
#citation(package = "DHARMa")


#clear workspace
rm(list = ls()) 

#set language to English
Sys.setenv(LANG = "en")



####________________________####
####<TANK TEMPERATURE>####


#clear workspace
rm(list = ls()) 

####1) Libraries####

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(car)

####2) Data####

tank_temp_avg <- read_csv("tank_temp_avg.csv")


#class check

summary(tank_temp_avg)

tank_temp_avg$Gregorian_date <- as.Date(tank_temp_avg$Gregorian_date, format = "%Y%m%d")
tank_temp_avg$system <- as.factor(tank_temp_avg$system) 
tank_temp_avg$lunar_day<- as.factor(tank_temp_avg$lunar_day) 
tank_temp_avg$lunar_month <- as.factor(tank_temp_avg$lunar_month) 
tank_temp_avg$Gregorian_month <- as.factor(tank_temp_avg$Gregorian_month) 
tank_temp_avg$season <- as.factor(tank_temp_avg$season) 

summary(tank_temp_avg)
str(tank_temp_avg)
head(tank_temp_avg)



####3) Plot####

tank_temp_plot  <- 
  ggplot(tank_temp_avg, aes(x=Gregorian_date, y=temp, colour=system)) +
  geom_line(size=1)+
  geom_smooth()+
  scale_color_manual(values=c("steelblue", "seagreen"))+
  theme_classic()+
  theme(text=element_text(size=14,  family="sans"))+
  scale_x_date(date_labels="%m",date_breaks  ="1 month")+
  scale_y_continuous(limits=c(22, 30), breaks=c(22, 24, 26, 28, 30))+
  labs(title="Tank Temperature", y="Temperature (⁰C)", x="Gregorian Month")

tank_temp_plot

tank_temp_plot + theme(legend.position = "none")


#calculate tank mean and sd for each season

#remove NA rows to allow for mean calculation
season_avg_no.NA <- tank_temp_avg %>% filter(!is.na(temp))

season_avg <-
  season_avg_no.NA%>%
  group_by(season, system)%>%
  summarise(mean.temp = mean(temp), sd.temp =sd(temp))

season_avg 


####4) Analysis####

temp.mm <- lmer(temp~system*season + (1|lunar_month), data=tank_temp_avg)
summary(temp.mm)

#check assumptions 
plot(fitted(temp.mm),residuals(temp.mm))
hist(residuals(temp.mm))
qqnorm(residuals(temp.mm))
vif(temp.mm)
acf(resid(temp.mm))

#between system seasonal post hoc
temp_posthoc1 <- emmeans(temp.mm, pairwise~system |season)
temp_posthoc1

#within system seasonal post hoc
temp_posthoc2 <-emmeans(temp.mm, pairwise ~ season | system)
temp_posthoc2 



#Just for comparison check lm (note: no inclusion of random effect results in high ACF)
temp.lm <- lm(temp~system*season, data=tank_temp_avg)
summary(temp.lm)

#check assumptions 
plot(fitted(temp.lm),residuals(temp.lm))
hist(residuals(temp.lm))
qqnorm(residuals(temp.lm))
vif(temp.lm)
acf(resid(temp.lm))


####________________________####
####<WATER CHEMISTRY>#####

#clear workspace
rm(list = ls())

####1) Libraries####

library(tidyverse)
library(lme4)
library(lmerTest)
library(car)
library(AER)

####2) Data####

water_chem_WIDE <- read_csv("water_chem.csv")


#change from wide to long format

water_chem_LONG <-
  water_chem_WIDE%>%      
  pivot_longer(cols=c('Ca', 'Mg', 'NH3', 'NO2', 'NO3', 
                      'PO4', 'Alkalinity', 'pH', 'Salinity'),
               names_to='parameter',
               values_to='value')

#rename
water_chem <- water_chem_LONG 


#class check
summary(water_chem)

water_chem$lunar_month <- as.factor(water_chem$lunar_month)
water_chem$Gregorian_month <- as.factor(water_chem$Gregorian_month)
water_chem$system <- as.factor(water_chem$system)
water_chem$tank_id <- as.factor(water_chem$tank_id)
water_chem$time <- as.numeric(water_chem$time)
water_chem$parameter <- as.factor(water_chem$parameter)
water_chem$value <- as.numeric(water_chem$value)

summary(water_chem)
head(water_chem)
str(water_chem)


# system subsets

FTS_chem <-
  water_chem%>%
  filter(system =='FTS')

RAS_chem <-
  water_chem%>%
  filter(system =='RAS')



#Parameter subsets

Ca_subset <-
  water_chem%>%
  filter(parameter =="Ca")

Mg_subset <-
  water_chem%>%
  filter(parameter =="Mg")

NH3_subset <-
  water_chem%>%
  filter(parameter =="NH3")

NO2_subset <-
  water_chem%>%
  filter(parameter =="NO2")

NO3_subset <-
  water_chem%>%
  filter(parameter =="NO3")

PO4_subset <-
  water_chem%>%
  filter(parameter =="PO4")

Alkalinity_subset <-
  water_chem%>%
  filter(parameter =="Alkalinity")

pH_subset <-
  water_chem%>%
  filter(parameter =="pH")

salinity_subset <-
  water_chem%>%
  filter(parameter =="Salinity")



####3) Plots####

#system summary
water_chem_sum <-
  water_chem%>%
  group_by(system, parameter, time)%>%
  summarise(mean = mean(value), n = n(), sd = sd(value), se = sd / sqrt(n))


water_chem_sum$parameter <- factor(water_chem_sum$parameter, 
                                   levels=c("Ca", "Mg", "NH3", "NO2", "NO3", 
                                            "PO4", "Alkalinity", "pH", "Salinity"))




#system comparison plot
chem_plot <-
  ggplot(water_chem_sum, aes(x=time, y=mean, group = system, color=system))+ 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, 
                position=position_dodge(0.05)) +
  geom_line(aes(linetype=system)) + 
  geom_point(aes(shape=system), size =2)+
  scale_x_continuous(limits=c(1, 12), breaks=c(1,2,3,4,5,6,7,8,9,10,11,12))+
  labs(title="",x="", y = "")+
  theme_bw()+
  facet_wrap(~parameter, scales="free")

chem_plot

chem_plot + scale_color_manual(values=c('steelblue','seagreen'))

chem_plot + scale_color_manual(values=c('steelblue','seagreen'))+
  theme(legend.position = "none") +
  theme(text=element_text(size=11,  family="sans"))



#tank summary (for table)
water_chem_tank_sum <-
  water_chem%>%
  group_by(system, parameter)%>%
  summarise(mean = mean(value), n = n(), sd = sd(value), se = sd / sqrt(n))

water_chem_tank_sum 



####4) Analysis####


####*Ca####

plot(Ca_subset$value)
hist(Ca_subset$value)

#lmer mixed model (2 random effects)
Ca_mm.2r <- lmer(value~system + (1|tank_id) + (1|lunar_month), data= Ca_subset)
summary(Ca_mm.2r)
#singularity warming message

rand(Ca_mm.2r)
#tank_id was not significant, so dropped as a random effect from model
#lunar_month was significant; kept in model due to repeated measures

#lmer mixed model (1 random effect)
Ca_mm.1r <- lmer(value~system + (1|lunar_month), data= Ca_subset)
summary(Ca_mm.1r)

#compare models
anova(Ca_mm.2r, Ca_mm.1r)

#check assumptions of model with lowest AIC
plot(fitted(Ca_mm.1r),residuals(Ca_mm.1r))
hist(residuals(Ca_mm.1r))
qqnorm(residuals(Ca_mm.1r))



####*Mg####

plot(Mg_subset$value)
hist(Mg_subset$value)

#lmer mixed model (2 random effects)
Mg_mm.2r <- lmer(value~system + (1|tank_id) + (1|lunar_month), data= Mg_subset)
summary(Mg_mm.2r)
#singularity warming message

rand(Mg_mm.2r)
#tank_id was not significant, so dropped as a random effect from model
#lunar_month was not significant; but kept in model due to repeated measures

#lmer mixed model (1 random effect)
Mg_mm.1r <- lmer(value~system + (1|lunar_month), data= Mg_subset)
summary(Mg_mm.1r)

#compare models
anova(Mg_mm.2r, Mg_mm.1r)

#check assumptions of model with lowest AIC
plot(fitted(Mg_mm.1r),residuals(Mg_mm.1r))
hist(residuals(Mg_mm.1r))
qqnorm(residuals(Mg_mm.1r))


####*NH3####

plot(NH3_subset$value)
hist(NH3_subset$value)

#use a tobit model approach because of zero inflation 

NH3.tobit <- tobit(formula = value~system, data = NH3_subset)
summary(NH3.tobit)

plot(fitted(NH3.tobit),residuals(NH3.tobit))
hist(residuals(NH3.tobit))
qqnorm(residuals(NH3.tobit))


####*NO2####

plot(NO2_subset$value)
hist(NO2_subset$value)

#no analysis undertaken since only one non-zero value in data set


####*NO3####

plot(NO3_subset$value)
hist(NO3_subset$value)

#use a tobit model approach because of zero inflation 

NO3.tobit <- tobit(formula = value~system, data = NO3_subset)
summary(NO3.tobit)

plot(fitted(NO3.tobit),residuals(NO3.tobit))
hist(residuals(NO3.tobit))
qqnorm(residuals(NO3.tobit))



####*PO4####

plot(PO4_subset$value)
hist(PO4_subset$value)

#use a tobit model approach because of zero inflation 

PO4.tobit <- tobit(formula = value~system, data = PO4_subset)
summary(PO4.tobit)

plot(fitted(PO4.tobit),residuals(PO4.tobit))
hist(residuals(PO4.tobit))
qqnorm(residuals(PO4.tobit))

#But no non-zero data for the RAS, so this analysis was deemed not appropriate


####*alkalinity####

plot(Alkalinity_subset$value)
hist(Alkalinity_subset$value)


#lmer mixed model (2 random effects)
Alkalinity_mm.2r <- lmer(value~system + (1|tank_id) + (1|lunar_month), data= Alkalinity_subset)
summary(Alkalinity_mm.2r)
#singularity warming message

rand(Alkalinity_mm.2r)
#tank_id was not significant, so dropped as a random effect from model
#lunar_month was not significant, but kept in model due to repeated measures

#lmer mixed model (1 random effect)
Alkalinity_mm.1r <- lmer(value~system + (1|lunar_month), data= Alkalinity_subset)
summary(Alkalinity_mm.1r)

#compare models
anova(kH_mm.2r, kH_mm.1r)

#check assumptions of model with lowest AIC
plot(fitted(Alkalinity_mm.1r),residuals(Alkalinity_mm.1r))
hist(residuals(Alkalinity_mm.1r))
qqnorm(residuals(Alkalinity_mm.1r))




####*pH####

plot(pH_subset$value)
hist(pH_subset$value)

#lmer mixed model (2 random effects)
pH_mm.2r <- lmer(value~system + (1|tank_id) + (1|lunar_month), data= pH_subset)
summary(pH_mm.2r)
#singularity warming message

rand(pH_mm.2r)
#tank_id was not significant, so dropped as a random effect from model
#lunar_month was significant; kept in model due to repeated measures

#lmer mixed model (1 random effect)
pH_mm.1r <- lmer(value~system + (1|lunar_month), data= pH_subset)
summary(pH_mm.1r)

#compare models
anova(pH_mm.2r, pH_mm.1r)

#check assumptions of model with lowest AIC
plot(fitted(pH_mm.1r),residuals(pH_mm.1r))
hist(residuals(pH_mm.1r))
qqnorm(residuals(pH_mm.1r))



####*salinity####

plot(salinity_subset$value)
hist(salinity_subset$value)

#lmer mixed model (2 random effects)
salinity_mm.2r <- lmer(value~system + (1|tank_id) + (1|lunar_month), data= salinity_subset)
summary(salinity_mm.2r)
#singularity warming message

rand(salinity_mm.2r)
#tank_id was not significant, so dropped as a random effect from model
#lunar_month was significant; kept in model due to repeated measures

#lmer mixed model (1 random effect)
salinity_mm.1r <- lmer(log(value)~system + (1|lunar_month), data= salinity_subset)
summary(salinity_mm.1r)

#compare models
anova(salinity_mm.2r, salinity_mm.1r)

#check assumptions of model with lowest AIC
plot(fitted(salinity_mm.1r),residuals(salinity_mm.1r))
hist(residuals(salinity_mm.1r))
qqnorm(residuals(salinity_mm.1r))




####________________________####
####<REPRODUCTIVE OUTPUT>####


#clear workspace
rm(list = ls()) 

####1) Libraries####

library(tidyverse)
library(emmeans)
library(glmmTMB)
library(DHARMa)


####2) Data####

RO_data <- read_csv("RO_data.csv")

####***class check####

#RO_data
summary(RO_data)

RO_data$system <- as.factor(RO_data$system)
RO_data$colony_id <- as.factor(RO_data$colony_id)
RO_data$tank_id <- as.factor(RO_data$tank_id)
RO_data$reef_site <- as.factor(RO_data$reef_site)
RO_data$study_month <- as.numeric(RO_data$study_month)
RO_data$season <- as.factor(RO_data$season)
RO_data$lunar_month_year <- as.factor(RO_data$lunar_month_year)
RO_data$gregorian_month <- as.factor(RO_data$gregorian_month)
RO_data$Lday1_date <- as.Date(RO_data$Lday1_date)
RO_data$RO <- as.integer(RO_data$RO)
RO_data$mean_monthly_temp <- as.numeric(RO_data$mean_monthly_temp)
RO_data$colony_size <- as.numeric(RO_data$colony_size)
RO_data$number_lunar_days <- as.numeric(RO_data$number_lunar_days)


summary(RO_data)
head(RO_data)
str(RO_data)

####***system subsets####

#FTS
FTS_subset <-
  RO_data%>%
  filter(system =="FTS")

summary(FTS_subset)


#RAS
RAS_subset <-
  RO_data%>%
  filter(system =="RAS")

summary(RAS_subset)




####3) Plots####

####*FTS - monthly barplot####

#set order for x axis
FTS_subset$lunar_month_year <- factor(FTS_subset$lunar_month_year, 
                                      levels=c("Feb_2021", "Mar_2021", "Apr_2021", "May_2021", "Jun_2021", "Jul_2021", "Aug_2021", "Sep_2021","Oct_2021", "Nov_2021", "Dec_2021", "Jan_2022"))
#set colony_id order
FTS_subset$colony_id <- factor(FTS_subset$colony_id, 
                               levels=c("FTS_colony1", "FTS_colony2", "FTS_colony3", "FTS_colony4", "FTS_colony5", "FTS_colony6", 
                                        "FTS_colony7", "FTS_colony8", "FTS_colony9", "FTS_colony10", "FTS_colony11", "FTS_colony12")) 


#Total RO by month
FTS_barplot <- 
  ggplot(FTS_subset, aes(fill=colony_id, y=RO, x=lunar_month_year)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = c("#cee9f9", "#aedaf5", "#8ecbf1", "#6dbced", "#4dade9", "#2d9ee4", "#1b8cd2", "#1677b2", "#126192", "#0e4b71", "#0a3651", "#062031"))+
  scale_y_continuous(limits=c(0, 40000), breaks=c(0, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000))+
  labs(y = "# of planulae released", x="Lunar month", title ="Reproductive output: Flow-through system")+
  theme_classic()+
  theme(legend.position="none")+
  theme(text=element_text(size=12,  family="sans"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.2))

FTS_barplot  


####*FTS - facet barplot####
#Colony total RO by month (facet wrap)
FTS_barplot_facet <- 
  ggplot(FTS_subset, aes(fill=colony_id, y=RO, x=lunar_month_year)) + 
  geom_bar(position="stack", stat="identity", colour="black") +
  scale_fill_manual(values = c("#cee9f9", "#aedaf5", "#8ecbf1", "#6dbced", "#4dade9", "#2d9ee4", "#1b8cd2", "#1677b2", "#126192", "#0e4b71", "#0a3651", "#062031"))+
  scale_y_continuous(limits=c(0, 7000), breaks=c(0, 2000, 4000, 6000))+
  labs(y = "# of planulae released", x="Lunar month", title ="Reproductive output: Flow-through system")+
  theme_bw()+
  facet_wrap(~colony_id)

FTS_barplot_facet + theme(axis.text.x = element_text(angle = 90, vjust = 0.4)) + 
  theme(legend.position = "none")+
  theme(text=element_text(size=11,  family="sans")) 


####*RAS - monthly barplot####

#set order for x axis
RAS_subset$lunar_month_year <- factor(RAS_subset$lunar_month_year, 
                                      levels=c("Feb_2021", "Mar_2021", "Apr_2021", "May_2021", "Jun_2021", "Jul_2021", "Aug_2021", "Sep_2021","Oct_2021", "Nov_2021", "Dec_2021", "Jan_2022"))


RAS_subset$colony_id <- factor(RAS_subset$colony_id, 
                               levels=c("RAS_colony1", "RAS_colony2", "RAS_colony3", "RAS_colony4", "RAS_colony5", "RAS_colony6", 
                                        "RAS_colony7", "RAS_colony8", "RAS_colony9", "RAS_colony10", "RAS_colony11", "RAS_colony12", 
                                        "RAS_colony13", "RAS_colony14", "RAS_colony15"))


#Total RO by month
RAS_barplot <- 
  ggplot(RAS_subset, aes(fill=colony_id, y=RO, x=lunar_month_year)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = c("#b8fedc", "#41fda1", "#7cfebe", "#53fdaa", "#02e878", "#02d46d", "#02c063", "#02ac59", "#02974e", "#018344", "#016f39", "#015b2f", "#014724", "#01321a", "#001e10"))+
  scale_y_continuous(limits=c(0, 40000), breaks=c(0, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000))+
  labs(y = "# of planulae released", x="Lunar month", title ="Reproductive output: Recirculating system") +
  theme_classic()+
  theme(legend.position="none")+
  theme(text=element_text(size=12,  family="sans"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.2))

RAS_barplot



####*RAS - facet barplot####
#Colony total RO by month (facet wrap)
RAS_barplot_facet <- 
  ggplot(RAS_subset, aes(fill=colony_id, y=RO, x=lunar_month_year)) + 
  geom_bar(position="stack", stat="identity", colour="black") +
  scale_fill_manual(values = c("#b8fedc", "#41fda1", "#7cfebe", "#53fdaa", "#02e878", "#02d46d", "#02c063", "#02ac59", "#02974e", "#018344", "#016f39", "#015b2f", "#014724", "#01321a", "#001e10"))+
  scale_y_continuous(limits=c(0, 7000), breaks=c(0, 2000, 4000, 6000))+
  labs(y = "# of planulae released", x="Lunar month", title ="Reproductive output: Recirculating system")+
  theme_bw()+
  facet_wrap(~colony_id)

RAS_barplot_facet + theme(axis.text.x = element_text(angle = 90, vjust = 0.4)) + 
  theme(legend.position = "none")+
  theme(text=element_text(size=10,  family="sans"))

#note: y axis scaling omits outlier in RAS_colony12 in lunar May 2021 (RO=11848); 
#this is indicated in the figure caption



####*Both systems - boxplot####

#set order for x axis
RO_data$lunar_month_year <- factor(RO_data$lunar_month_year, 
                                   levels=c("Feb_2021", "Mar_2021", "Apr_2021", "May_2021", "Jun_2021", "Jul_2021", "Aug_2021", "Sep_2021","Oct_2021", "Nov_2021", "Dec_2021", "Jan_2022"))

#***by month####
#plot is just for reference; not included in manuscript
#note: 'season' was used in the glmmTMB model because the use of month did not allow for model convergence

system_boxplot <- 
  ggplot(data=RO_data, aes(x=lunar_month_year, y=RO, fill=system)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), 
             size=1.5) +
  scale_y_continuous(limits=c(0, 12000), breaks=c(0, 2000, 4000, 6000, 8000, 10000, 12000))+
  scale_fill_manual(breaks = c("FTS", "RAS"), 
                    values=c("steelblue", "seagreen"))+
  labs(y = "# of planulae released", x="Lunar Month")+
  theme_classic()

system_boxplot

system_boxplot + theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=0))+
  theme(text=element_text(size=12,  family="sans"))



#***by season####
#note: seasons were defined based on the lunar month's approximation on the Gregorian calender
#e.g., On the Gregorian calender we define 'spring' as March, April, May, therefore in the context 
#of lunar months lunar February, March, and April were defined as 'spring'. The Gregorian calender date
#for 'lunar day1' for each month is indicated in the RO_data.csv file in the column 'Lday1_date.


#set order for x axis
RO_data$season <- factor(RO_data$season, 
                         levels=c("spring", "summer", "fall", "winter"))


system_boxplot_season <- 
  ggplot(data=RO_data, aes(x=season, y=RO, fill=system)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), 
             size=1.5) +
  scale_y_continuous(limits=c(0, 12000), breaks=c(0, 2000, 4000, 6000, 8000, 10000, 12000))+
  scale_fill_manual(breaks = c("FTS", "RAS"), 
                    values=c("steelblue", "seagreen"))+
  labs(y = "# of planulae released", x="Season", title="Seasonal Reproductive Output")+
  theme_classic()

system_boxplot_season

system_boxplot_season + theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=0))+
  theme(text=element_text(size=12,  family="sans"))



####*Both systems - scatterplots####

#***colony_size####
colony_size_line <- 
  ggplot(RO_data, aes(x=colony_size, y=RO, colour=system)) +
  geom_jitter(width = 0.1, height = 0.1)+
  geom_smooth(method=lm , se=TRUE) +
  theme_bw()+
  scale_colour_manual(breaks = c("FTS", "RAS"), 
                      values=c("steelblue", "seagreen"))+
  labs(y = "# of planulae released", x="Total linear extension", title = "Reproductive output by colony size")+
  facet_wrap(~system)+
  theme(text=element_text(size=12,  family="sans"))

colony_size_line + theme(legend.position = "none")



#***temperature####

temp_line <- 
  ggplot(RO_data, aes(x=mean_monthly_temp, y=RO, colour=system)) +
  geom_jitter(width = 0.1, height = 0.1)+
  geom_smooth(method=lm , se=TRUE) +
  theme_bw()+
  scale_colour_manual(breaks = c("FTS", "RAS"), 
                      values=c("steelblue", "seagreen"))+
  labs(y = "# of planulae released", x="Temperature", title = "Reproductive output by mean monthly temperature")+
  facet_wrap(~system)+
  theme(text=element_text(size=12,  family="sans"))

temp_line + theme(legend.position = "none")



#***study month####

study_month_line <- 
  ggplot(RO_data, aes(x=study_month, y=RO, colour=system)) +
  geom_jitter(width = 0.1, height = 0.1)+
  geom_smooth(method=lm , se=TRUE) +
  theme_bw()+
  scale_x_continuous(limits=c(1, 12), breaks=c(1,2,3,4,5,6,7,8,9,10,11,12))+
  scale_colour_manual(breaks = c("FTS", "RAS"), 
                      values=c("steelblue", "seagreen"))+
  labs(y = "# of planulae released", x="Study month", title = "Reproductive output by study month")+
  facet_wrap(~system)+
  theme(text=element_text(size=12,  family="sans"))

study_month_line 

study_month_line  + theme(legend.position = "none")


####4) Analysis ####

####*Seasonal system summary####

RO_data <- RO_data %>% filter(!is.na(RO))

season_sum <-
  RO_data%>%
  group_by(system, season)%>%
  summarise(mean=mean(RO), sd=sd(RO))

season_sum


####*Main RO model (m1)####

summary(RO_data)
head(RO_data)
str(RO_data)

#take a look at reproductive output
plot(RO_data$RO)
hist(RO_data$RO)


#m1 (zero-inflated negative binomial mixed-effects model)

m1 <- glmmTMB(RO ~ system + colony_size + season +
                (system*colony_size) + (system*season) + 
                (1|tank_id) + (1|colony_id) + 
                ar1(lunar_month_year + 0|colony_id),
              zi = ~1,
              data = RO_data, 
              family = nbinom2())


summary(m1)
testDispersion(m1)
simulationOutput.m1 <- simulateResiduals(fittedModel = m1, plot = T)

#remove NAs
RO_data <- RO_data %>% filter(!is.na(system) & 
                                !is.na(season) &
                                !is.na(colony_size) &
                                !is.na(colony_id) &
                                !is.na(tank_id))



#compare residuals with each predictor
plotResiduals(simulationOutput.m1, form = RO_data$system) #ok
plotResiduals(simulationOutput.m1, form = RO_data$season) #ok
plotResiduals(simulationOutput.m1, form = RO_data$colony_size) #ok
plotResiduals(simulationOutput.m1, form = RO_data$colony_id) #ok
plotResiduals(simulationOutput.m1, form = RO_data$tank_id) #OK



RO_posthoc1 <- emmeans(m1, pairwise~system |season)
RO_posthoc1


RO_posthoc2 <-emmeans(m1, pairwise ~ season | system)
RO_posthoc2 



####*Temp RO model (m2)####

#clear workspace and load data again so removal of NAs in m1 analyses doesn't influence m2 model
rm(list = ls()) 

RO_data <- read_csv("RO_data.csv")

#RO_data
summary(RO_data)

RO_data$system <- as.factor(RO_data$system)
RO_data$colony_id <- as.factor(RO_data$colony_id)
RO_data$tank_id <- as.factor(RO_data$tank_id)
RO_data$reef_site <- as.factor(RO_data$reef_site)
RO_data$study_month <- as.numeric(RO_data$study_month)
RO_data$season <- as.factor(RO_data$season)
RO_data$lunar_month_year <- as.factor(RO_data$lunar_month_year)
RO_data$gregorian_month <- as.factor(RO_data$gregorian_month)
RO_data$Lday1_date <- as.Date(RO_data$Lday1_date)
RO_data$RO <- as.integer(RO_data$RO)
RO_data$mean_monthly_temp <- as.numeric(RO_data$mean_monthly_temp)
RO_data$colony_size <- as.numeric(RO_data$colony_size)
RO_data$number_lunar_days <- as.numeric(RO_data$number_lunar_days)


summary(RO_data)
head(RO_data)
str(RO_data)

#create scaled temperature column
RO_data <-
  RO_data %>%
  mutate(mean_monthly_temp_scaled = scale(mean_monthly_temp))

plot(RO_data$mean_monthly_temp)
plot(RO_data$mean_monthly_temp_scaled)


m2<- glmmTMB(RO ~ mean_monthly_temp_scaled * system + 
               (1|tank_id) + (1|colony_id) + 
               ar1(lunar_month_year + 0|colony_id),
             zi = ~1,
             data = RO_data, 
             family = nbinom2())


summary(m2)
testDispersion(m2)
simulationOutput.m2 <- simulateResiduals(fittedModel = m2, plot = T)


#remove NAs
RO_data <- RO_data %>% filter(!is.na(mean_monthly_temp_scaled) &
                                !is.na(system) &
                                !is.na(colony_id) &
                                !is.na(tank_id))

#compare residuals with each predictor
plotResiduals(simulationOutput.m2, form = RO_data$mean_monthly_temp) #not perfect, but OK
plotResiduals(simulationOutput.m2, form = RO_data$system) #OK
plotResiduals(simulationOutput.m2, form = RO_data$colony_id) #OK
plotResiduals(simulationOutput.m2, form = RO_data$tank_id) #OK



####________________________####
####<REPRODUCTIVE TIMING>####


#clear workspace
rm(list = ls()) 

####1) Libraries####

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(car)
library(ggridges)
library(Hmisc)

#library(hrbrthemes)
#library(viridis)
#library(writexl)


####2) Data####

RT_data <- read_csv("RT_data.csv")

#class check

summary(RT_data)

RT_data$Gregorian_date <- as.Date(RT_data$Gregorian_date, format = "%Y%m%d") 
RT_data$lunar_day <- as.numeric(RT_data$lunar_day) 
RT_data$number_lunar_days <- as.numeric(RT_data$number_lunar_days) 
RT_data$lunar_month <- as.factor(RT_data$lunar_month) 
RT_data$Gregorian_month <- as.factor(RT_data$Gregorian_month) 
RT_data$season <- as.factor(RT_data$season) 
RT_data$mean_daily_temp <- as.numeric(RT_data$mean_daily_temp) 
RT_data$system <- as.factor(RT_data$system) 
RT_data$tank_id <- as.factor(RT_data$tank_id) 
RT_data$colony_id <- as.factor(RT_data$colony_id) 
RT_data$RO<- as.integer(RT_data$RO) 
RT_data$RT_angle<- as.numeric(RT_data$RT_angle)


summary(RT_data)
str(RT_data)
head(RT_data)

#create column with mean monthly temp
RT_data <- 
  RT_data%>%
  group_by(tank_id, lunar_month) %>%
  mutate(month.temp = mean(mean_daily_temp , na.rm= T))

#remove NA rows
RT_data_no.NA <- RT_data %>% filter(!is.na(RT_angle))

summary(RT_data_no.NA)

#calculate lunar day
RT_data_no.NA <-
  RT_data_no.NA %>%
  mutate(LD = number_lunar_days*RT_angle/360)


####*weighted mean and sd####

#add in monthly sum to assess RT weight based on RO
RT_data_no.NA_sum <- 
  RT_data_no.NA%>%
  group_by(colony_id, lunar_month) %>%
  mutate(month.sum = sum(RO , na.rm= T))

#calculate RT weight
RT_data_no.NA_sum <- 
  RT_data_no.NA_sum%>%
  mutate(RT_weight = RO/month.sum)

#calculate weighted mean of lunar day
WT.MLD <- 
  RT_data_no.NA_sum%>%
  group_by(colony_id, system, lunar_month, tank_id, season, month.temp)%>%
  summarise(W.MLD = weighted.mean(LD,RT_weight), MLD = mean(LD))

#weighted mean and sd per season for each system

WT.MLD.SD <- 
  RT_data_no.NA_sum%>%
  group_by(system,season)%>%
  summarise(W.MLD = weighted.mean(LD,RT_weight), MLD = mean(LD), wt.var = wtd.var(LD, RT_weight), wt.sd = sqrt(wt.var))

WT.MLD.SD 

####*system subsets####

#FTS
FTS_RT_subset <-
  RT_data_no.NA_sum%>%
  filter(system =="FTS")

summary(FTS_RT_subset)

#RAS
RAS_RT_subset <-
  RT_data_no.NA_sum%>%
  filter(system =="RAS")

summary(RAS_RT_subset)


####3) Plots####

####*FTS: weighted LD####

#set lunar_month order
FTS_RT_subset$lunar_month <- factor(FTS_RT_subset$lunar_month, 
                                    levels=c("Jan_2022", "Dec_2021", "Nov_2021", "Oct_2021", "Sep_2021", "Aug_2021", "Jul_2021","Jun_2021", "May_2021", "Apr_2021", "Mar_2021", "Feb_2021"))

#set colony_id order
FTS_RT_subset$colony_id <- factor(FTS_RT_subset$colony_id, 
                                  levels=c("FTS_colony1", "FTS_colony2", "FTS_colony3", "FTS_colony4", "FTS_colony5", "FTS_colony6", 
                                           "FTS_colony7", "FTS_colony8", "FTS_colony9", "FTS_colony10", "FTS_colony11", "FTS_colony12")) 

FTS_W.LD_Rplot <- 
  ggplot(data=FTS_RT_subset, aes(y=lunar_month, x=LD)) +
  geom_density_ridges(aes(height=..density..,
                          weight=RT_weight),    
                      scale= 0.95,
                      stat="density", fill="steelblue")+
  scale_x_continuous(limits=c(1, 30), breaks=c(1,5,10,15,20,25,30))+
  geom_vline(xintercept=15, linetype="dashed", 
             color = "black", size=0.5)+
  theme_classic()+
  theme(text=element_text(size=12,  family="sans"))

FTS_W.LD_Rplot

ggsave("FTS_W.LD_Rplot.png")


####*FTS: facet weighted LD####

#Independent colony plots
FTS_colony_RT_Rplot <- 
  FTS_W.LD_Rplot +   theme_bw()+ facet_wrap(~colony_id) +  theme(text=element_text(size=9,  family="sans"))

FTS_colony_RT_Rplot

ggsave("FTS_colony_RT_Rplot.png")


####*RAS: weighted LD####

#set lunar_month order
RAS_RT_subset$lunar_month <- factor(RAS_RT_subset$lunar_month, 
                                    levels=c("Jan_2022", "Dec_2021", "Nov_2021", "Oct_2021", "Sep_2021", "Aug_2021", "Jul_2021","Jun_2021", "May_2021", "Apr_2021", "Mar_2021", "Feb_2021"))
#set colony_id order
RAS_RT_subset$colony_id <- factor(RAS_RT_subset$colony_id, 
                                  levels=c("RAS_colony1", "RAS_colony2", "RAS_colony3", "RAS_colony4", "RAS_colony5", "RAS_colony6", 
                                           "RAS_colony7", "RAS_colony8", "RAS_colony9", "RAS_colony10", "RAS_colony11", "RAS_colony12", 
                                           "RAS_colony13", "RAS_colony14", "RAS_colony15"))

RAS_W.LD_Rplot <- 
  ggplot(data=RAS_RT_subset, aes(y=lunar_month, x=LD)) +
  geom_density_ridges(aes(height=..density..,
                          weight=RT_weight),    
                      scale= 0.95,
                      stat="density", fill="seagreen")+
  scale_x_continuous(limits=c(1, 30), breaks=c(1,5,10,15,20,25,30))+
  geom_vline(xintercept=15, linetype="dashed", 
             color = "black", size=0.5)+
  theme_classic()+
  theme(text=element_text(size=12,  family="sans"))

RAS_W.LD_Rplot

ggsave("RAS_W.LD_Rplot.png")


####*RAS: facet weighted LD####
RAS_colony_RT_Rplot <- 
  RAS_W.LD_Rplot + theme_bw()+  facet_wrap(~colony_id) +  theme(text=element_text(size=9,  family="sans"))

RAS_colony_RT_Rplot

ggsave("RAS_colony_RT_Rplot.png")




####*Temp vs MLD####

temp_WT.MLD <- 
  ggplot(WT.MLD, aes(x=month.temp, y=W.MLD, colour=system)) +
  geom_jitter(width = 0.1, height = 0.1)+
  geom_smooth(method=lm , se=TRUE) +
  theme_bw()+
  scale_colour_manual(breaks = c("FTS", "RAS"), 
                      values=c("steelblue", "seagreen"))+
  labs(y = "Weighted MLD", x="Temperature", title = "MLD by mean monthly temperature")+
  facet_wrap(~system)+
  theme(legend.position = "none")+
  theme(text=element_text(size=12,  family="sans"))

temp_WT.MLD 

ggsave("temp_WT.MLD.png")



####*Reproductive days####

RT_Rday_count <-
  RT_data_no.NA %>%
  group_by(colony_id, lunar_month, system, number_lunar_days, season, tank_id)%>%
  summarise(n = n())

summary(RT_Rday_count)


#system subsets
FTS_Rday_count <-
  RT_Rday_count%>%
  filter(system =="FTS")

RAS_Rday_count <-
  RT_Rday_count%>%
  filter(system =="RAS")



#to calculate percentage of colonies reproducing each month
Pcol <- 
  RT_Rday_count%>%
  group_by(lunar_month, system)%>%
  summarise(n.reproducing = n())

Pcol$lunar_month <- factor(Pcol$lunar_month, 
                           levels=c("Feb_2021", "Mar_2021", "Apr_2021", "May_2021", "Jun_2021", "Jul_2021", "Aug_2021", "Sep_2021","Oct_2021", "Nov_2021", "Dec_2021", "Jan_2022"))

Pcol <- 
  Pcol%>%
  mutate(percentage = n.reproducing/12*100)


Pcol <- 
  Pcol%>%
  group_by(system)%>%
  summarise(months.mean = mean(percentage), sd = sd (percentage))

Pcol
#FTS           79.9  20.2
#RAS           70.8  26.5


#FTS
#colony mean reproductive days +/-se across months

FTS_Rdays <-
  FTS_Rday_count%>%
  group_by(lunar_month)%>%
  summarise(Rday.mean = mean (n), Rday.sd = sd(n), sample.size = n(), Rday.se = Rday.sd / sqrt(sample.size))

#overall mean reproductive days
mean(FTS_Rdays$Rday.mean)
#8.023782 
sd(FTS_Rdays$Rday.mean)
#3.076652



#All FTS colonies

#set order for x axis
FTS_Rdays$lunar_month <- factor(FTS_Rdays$lunar_month, 
                                levels=c("Feb_2021", "Mar_2021", "Apr_2021", "May_2021", "Jun_2021", "Jul_2021", "Aug_2021", "Sep_2021","Oct_2021", "Nov_2021", "Dec_2021", "Jan_2022"))


FTS_Rdays_plot <-
  ggplot(FTS_Rdays) +
  geom_bar( aes(x=lunar_month, y=Rday.mean), stat="identity", fill="steelblue", alpha=0.7, color="black")+ 
  geom_errorbar( aes(x=lunar_month, ymin=Rday.mean, ymax=Rday.mean+Rday.se), width=0.2, colour="black", alpha=0.9, size=0.5)+
  labs(y = "Number of reproductive days)", x="Lunar month", title ="Number of Reproductive Days: Flow-through system")+
  scale_y_continuous(limits=c(0, 16), breaks=c(0, 5, 10, 15))+
  theme_classic()+
  theme(text=element_text(size=12,  family="sans"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.2))

FTS_Rdays_plot



#RAS

RAS_Rdays <-
  RAS_Rday_count %>%
  group_by(lunar_month)%>%
  summarise(Rday.mean = mean (n), Rday.sd = sd(n), sample.size = n(), Rday.se = Rday.sd / sqrt(sample.size))


#overall mean reproductive days
mean(RAS_Rdays$Rday.mean)
#6.580913
sd(RAS_Rdays$Rday.mean)
#3.075451


#set order for x axis
RAS_Rdays$lunar_month <- factor(RAS_Rdays$lunar_month, 
                                levels=c("Feb_2021", "Mar_2021", "Apr_2021", "May_2021", "Jun_2021", "Jul_2021", "Aug_2021", "Sep_2021","Oct_2021", "Nov_2021", "Dec_2021", "Jan_2022"))

#All FTS colonies
RAS_Rdays_plot <-
  ggplot(RAS_Rdays) +
  geom_bar( aes(x=lunar_month, y=Rday.mean), stat="identity", fill="seagreen", alpha=0.7, color="black")+ 
  geom_errorbar( aes(x=lunar_month, ymin=Rday.mean, ymax=Rday.mean+Rday.se), width=0.2, colour="black", alpha=0.9, size=0.5)+
  labs(y = "Number of reproductive days)", x="Lunar month", title ="Number of Reproductive Days: Recirculating aquaculture system")+
  scale_y_continuous(limits=c(0, 16), breaks=c(0, 5, 10, 15))+
  theme_classic()+
  theme(text=element_text(size=12,  family="sans"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.2))

RAS_Rdays_plot



####4) Analysis####

###*Weighted MLD####

WT.MLD.mm <- lmer(W.MLD~system*season + (1|tank_id) + (1|colony_id), data=WT.MLD)
summary(WT.MLD.mm)

#check assumptions
plot(fitted(WT.MLD.mm),residuals(WT.MLD.mm))
hist(residuals(WT.MLD.mm))
qqnorm(residuals(WT.MLD.mm))
vif(WT.MLD.mm)
acf(resid(WT.MLD.mm))

WT.MLD_posthoc1 <- emmeans(WT.MLD.mm, pairwise~system |season)
WT.MLD_posthoc1


WT.MLD_posthoc2 <-emmeans(WT.MLD.mm, pairwise ~ season | system)
WT.MLD_posthoc2 




####*Weighted MLD and temp####

#FTS
FTS_MLD_subset <-
  WT.MLD %>%
  filter(system =="FTS")

summary(FTS_MLD_subset)

#RAS
RAS_MLD_subset <-
  WT.MLD %>%
  filter(system =="RAS")

summary(RAS_MLD_subset)



#FTS
FTS.temp.MLD.lm <- lm(W.MLD~month.temp, data=FTS_MLD_subset)
summary(FTS.temp.MLD.lm)

#check assumptions
plot(fitted(FTS.temp.MLD.lm),residuals(FTS.temp.MLD.lm))
hist(residuals(FTS.temp.MLD.lm))
qqnorm(residuals(FTS.temp.MLD.lm))



#RAS
RAS.temp.MLD.lm <- lm(W.MLD~month.temp, data=RAS_MLD_subset)
summary(RAS.temp.MLD.lm)

#check assumptions
plot(fitted(RAS.temp.MLD.lm),residuals(RAS.temp.MLD.lm))
hist(residuals(RAS.temp.MLD.lm))
qqnorm(residuals(RAS.temp.MLD.lm))




####*Reproductive days####

summary(RT_Rday_count)

#for monthly mean sd / system
Rday_analysis <-
  RT_Rday_count %>%
  group_by(lunar_month, system)%>%
  summarise(Rday.mean = mean (n), Rday.sd = sd(n), sample.size = n(), Rday.se = Rday.sd / sqrt(sample.size))

Rday_analysis


Rdays.mm <- lmer(n~system*season + (1|tank_id) + (1|colony_id), data=RT_Rday_count)
summary(Rdays.mm)
#singularity warning message likely due to some colonies in RAS being replaced, i.e., not reported each month
#(due to mortality of other colonies)

#check assumptions
plot(fitted(Rdays.mm),residuals(Rdays.mm))
hist(residuals(Rdays.mm))
qqnorm(residuals(Rdays.mm))
vif(Rdays.mm)
acf(resid(Rdays.mm))

RT_posthoc1 <- emmeans(Rdays.mm, pairwise~system |season)
RT_posthoc1


RT_posthoc2 <-emmeans(Rdays.mm, pairwise ~ season | system)
RT_posthoc2 


####________________________####
####<DIURNAL REPRODUCTIVE TIMING>####


#clear workspace
rm(list = ls()) 

####1) Libraries####

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(car)

#library(gamlss)
#library(glmmTMB)
#library(DHARMa)

####2) Data####

RT_diurnal <- read_csv("RT_diurnal.csv")

#class check

summary(RT_diurnal)

RT_diurnal$light.dark <- as.factor(RT_diurnal$light.dark)
RT_diurnal$study_time <- as.numeric(RT_diurnal$study_time)
RT_diurnal$colony_id <- as.factor(RT_diurnal$colony_id)
RT_diurnal$tank_id <- as.factor(RT_diurnal$tank_id)
RT_diurnal$system <- as.factor(RT_diurnal$system)
RT_diurnal$RO <- as.integer(RT_diurnal$RO)

summary(RT_diurnal)
str(RT_diurnal)
head(RT_diurnal)


#system subsets
FTS_di <-
  RT_diurnal%>%
  filter(system =="FTS")

summary(FTS_di)

RAS_di <-
  RT_diurnal%>%
  filter(system =="RAS")

summary(RAS_di)


####3) Plots####

#FTS
#Total RO per hour (n =9 colonies; 3/tank)
FTS_barplot.hour <- 
  ggplot(FTS_di, aes(fill=colony_id, y=RO, x=study_time)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = c("#aedaf5", "#6dbced", "#4dade9", "#1b8cd2", "#1677b2", "#126192", "#0e4b71", "#0a3651", "#062031"))+
  scale_y_continuous(limits=c(0, 200), breaks=c(0, 50, 100, 150, 200))+
  scale_x_continuous(limits=c(0, 48), breaks=c(1, 6,12,18,24,30,36,42,48))+
  labs(y = "# of planulae released", x="Hour", title ="Diurnal reproductive output: Flow-through system")+
  theme_classic()+
  theme(legend.position="none")+
  theme(text=element_text(size=12,  family="sans"))

FTS_barplot.hour 


#RAS
#Total RO per hour (n = 9 colonies; 3/tank)
RAS_barplot.hour <- 
  ggplot(RAS_di, aes(fill=colony_id, y=RO, x=study_time)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = c("#b8fedc", "#41fda1", "#7cfebe", "#02e878", "#02d46d", "#02ac59", "#015b2f", "#01321a", "#001e10"))+
  scale_y_continuous(limits=c(0, 200), breaks=c(0, 50, 100, 150, 200))+
  scale_x_continuous(limits=c(1, 48), breaks=c(1, 6,12,18,24,30,36,42,48))+
  labs(y = "# of planulae released", x="Hour", title ="Diurnal reproductive output: Recirculating system") +
  theme_classic()+
  theme(legend.position="none")+
  theme(text=element_text(size=12,  family="sans"))

RAS_barplot.hour 


####4) Analysis####

RT_di_sum <-
  RT_diurnal%>%
  group_by(colony_id, light.dark, tank_id, system)%>%
  summarise(total = sum(per.release))

#for mean percent of light.dark release over the 48h period

mean48 <-
  RT_di_sum%>%
  group_by(system, light.dark)%>%
  summarise(system.mean = mean(total))


#mixed model
plot(RT_di_sum$total)
hist(RT_di_sum$total)


RT_di_mm <- lmer(total ~ system + light.dark +
                   (system*light.dark) + 
                   (1|tank_id) + (1|colony_id), data=RT_di_sum)



summary(RT_di_mm)


#check assumptions 
plot(fitted(RT_di_mm),residuals(RT_di_mm))
hist(residuals(RT_di_mm))
qqnorm(residuals(RT_di_mm))
vif(RT_di_mm)
acf(resid(RT_di_mm))

#between system light.dark post hoc
di_posthoc1 <- emmeans(RT_di_mm, pairwise~system |light.dark)
di_posthoc1

#within system light.dark post hoc
di_posthoc2 <- emmeans(RT_di_mm, pairwise~light.dark|system)
di_posthoc2





####________________________####
####<F2>#####

#clear workspace
rm(list = ls())

####1) Libraries####

library(tidyverse)
library(car)
library(lme4)
library(lmerTest)


####2) Data####

F2_data <- read_csv("F2_data.csv")


#class check

summary(F2_data)

F2_data$F2_id <- as.factor(F2_data$F2_id)
F2_data$tank_id <- as.factor(F2_data$tank_id)
F2_data$system <- as.factor(F2_data$system)
F2_data$RO <- as.integer(F2_data$RO)

summary(F2_data)
str(F2_data)
head(F2_data)

####3) Plots####


####*F2 initial size####

F2_RO_start.size <- 
  ggplot(data=F2_data, aes(x=system, y=initial_diameter, fill=system)) + 
  geom_boxplot()+
  geom_point(position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.01), 
             size=1.5) +
  scale_y_continuous(limits=c(1,3), breaks=c(0, 1, 2, 3))+
  scale_fill_manual(breaks = c("FTS", "RAS"), 
                    values=c("steelblue", "seagreen"))+
  labs(y = "Diamter (cm)", x="System")+
  theme_classic()

F2_RO_start.size +theme(text=element_text(size=12,  family="sans")) + theme(legend.position = "none")



####*F2 RO####

F2_RO <- 
  ggplot(data=F2_data, aes(x=system, y=RO, fill=system)) + 
  geom_violin()+
  geom_point(position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.01), 
             size=1.5) +
  scale_y_continuous(limits=c(0, 400), breaks=c(0, 100, 200, 300, 400))+
  scale_fill_manual(breaks = c("FTS", "RAS"), 
                    values=c("steelblue", "seagreen"))+
  labs(y = "# of planulae released", x="System")+
  theme_classic()

F2_RO +theme(text=element_text(size=12,  family="sans")) + theme(legend.position = "none")

ggsave("F2_RO.png")


#boxplot option
F2_RO_box <- 
  ggplot(data=F2_data, aes(x=system, y=RO, fill=system)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), 
             size=1.5) +
  scale_y_continuous(limits=c(0, 400), breaks=c(0, 100, 200, 300, 400))+
  scale_fill_manual(breaks = c("FTS", "RAS"), 
                    values=c("steelblue", "seagreen"))+
  labs(y = "# of planulae released", x="System")+
  theme_bw()

F2_RO_box +theme(text=element_text(size=12,  family="sans")) + theme(legend.position = "none")

#stack option
F2_RO_stack <- 
  ggplot(data=F2_data, aes(x=system, y=RO, fill=F2_id)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = c("#cee9f9", "#aedaf5", "#8ecbf1", "#6dbced", "#4dade9", "#2d9ee4", "#1b8cd2", "#1677b2", "#126192", "#0e4b71", "#0a3651", "#062031", "steelblue", "darkblue", "blue","#b8fedc", "#41fda1", "#7cfebe", "#53fdaa", "#02e878", "#02d46d", "#02c063", "#02ac59", "#02974e", "#018344", "#016f39", "#015b2f", "#014724", "#01321a", "#001e10"))+
  theme_classic()+
  theme(legend.position="none")

F2_RO_stack 

####*F2 size####
F2_TLE <- 
  ggplot(data=F2_data, aes(x=system, y=TLE, fill=system)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), 
             size=1.5) +
  scale_y_continuous(limits=c(5, 15), breaks=c(5,7.5, 10, 15, 12.5))+
  scale_fill_manual(breaks = c("FTS", "RAS"), 
                    values=c("steelblue", "seagreen"))+
  labs(y = "Total Linear Extension", x="System")+
  theme_classic()

F2_TLE + theme(text=element_text(size=12,  family="sans")) + theme(legend.position = "none")


####*F2 RO vs. TLE####
F2_RO_TLE <- 
  ggplot(F2_data, aes(x=TLE, y=RO, colour=system)) +
  geom_jitter(width = 0.1, height = 0.1)+
  geom_smooth(method=lm , se=TRUE) +
  theme_bw()+
  scale_colour_manual(breaks = c("FTS", "RAS"), 
                      values=c("steelblue", "seagreen"))+
  labs(y = "# of planulae released", x="Total linear extension", title = "Reproductive output by size")+
  facet_wrap(~system)+
  theme(text=element_text(size=12,  family="sans"))

F2_RO_TLE + theme(text=element_text(size=12,  family="sans")) + theme(legend.position = "none")

#No facet
F2_RO_TLE2 <- 
  ggplot(F2_data, aes(x=TLE, y=RO, colour=system)) +
  geom_jitter(width = 0.1, height = 0.1)+
  geom_smooth(method=lm , se=TRUE) +
  theme_classic()+
  scale_y_continuous(limits=c(-100, 400), breaks=c(0, 100, 200, 300, 400))+
  scale_x_continuous(limits=c(5, 15), breaks=c(5,7.5,10,12.5,15))+
  scale_colour_manual(breaks = c("FTS", "RAS"), 
                      values=c("steelblue", "seagreen"))+
  labs(y = "# of planulae released", x="Total linear extension", title = "Reproductive output by size")+
  theme(text=element_text(size=12,  family="sans"))

F2_RO_TLE2 + theme(text=element_text(size=12,  family="sans")) + theme(legend.position = "none")


ggsave("F2_RO_TLE2.png")

####4) Analysis####

FTS_F2 <-
  F2_data %>%
  filter(system =="FTS")

RAS_F2 <-
  F2_data %>%
  filter(system =="RAS")


####*Initial size####

#start mean, sd
start_sum <-
  F2_data%>%
  group_by(system)%>%
  summarise(mean.start = mean(initial_diameter), sd.start = sd(initial_diameter))

start_sum 

#test for normality and equal variance
qqPlot(FTS_F2$initial_diameter)
qqPlot(RAS_F2$initial_diameter)


#check for equal variance
leveneTest(initial_diameter~system, data=F2_data)


#t-test to check initial size
t.test(initial_diameter~system, data=F2_data, var.equal=TRUE)



####*size: system comparison####

F2_size <- lmer(TLE~system + (1|tank_id), data=F2_data)
summary(F2_size)  

#check assumptions 
plot(fitted(F2_size),residuals(F2_size))
hist(residuals(F2_size))
qqnorm(residuals(F2_size))


####*RO: system comparison####

F2_RO <- lmer(sqrt(RO)~system + (1|tank_id), data=F2_data)
summary(F2_RO)  

#check assumptions 
plot(fitted(F2_RO),residuals(F2_RO))
hist(residuals(F2_RO))
qqnorm(residuals(F2_RO))



####*RO and size####


#FTS
FTS_F2_mm <- lm(sqrt(RO)~TLE, data=FTS_F2)
summary(FTS_F2_mm) 

#check assumptions 
plot(fitted(FTS_F2_mm),residuals(FTS_F2_mm))
hist(residuals(FTS_F2_mm))
qqnorm(residuals(FTS_F2_mm))



#RAS
RAS_F2_mm <- lm(sqrt(RO)~TLE, data=RAS_F2)
summary(RAS_F2_mm) 

#check assumptions 
plot(fitted(RAS_F2_mm),residuals(RAS_F2_mm))
hist(residuals(RAS_F2_mm))
qqnorm(residuals(RAS_F2_mm))


####________________________####
