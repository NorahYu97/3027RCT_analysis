# Downloaded baseline char from EDC  
```{r summarized table}
library(data.table)
library(dplyr)
setwd("/DATA2/xuanjingyu/RQ3027 IIT")
sum_data <- fread("./Summarized information of vaccinated individuals with name.csv") %>% as.data.frame()
## Record COVID-19 confirmed cases since intention-to-treat population enrolled
sum_data[,8:14] <- "no"
### Post V1 infection
post_1 <- c("FB147","FB211","FB299","FB277","FB361")
sum_data[match(post_1,sum_data$ID),8] <- "yes"
### Post V2 infection
post_2 <- c("FB320","FB342","FB336","FB335","FB117")
sum_data[match(post_2,sum_data$ID),9] <- "yes"
### Post V3 infection
post_3 <- c("FB198","FB206","FB189","FB295","FB169","FB374","FB349","FB369")
sum_data[match(post_3,sum_data$ID),10] <- "yes"
### Post V4 infection
post_4 <- c("FB066","FB106","FB027","FB222","FB064","FB070","FB113","FB116",
          "FB120","FB126","FB135","FB150","FB151","FB175","FB105","FB273",
          "FB291","FB247","FB372","FB314","FB331","FB329","FB332")
sum_data[match(post_4,sum_data$ID),11] <- "yes"
### Post V5 infection
pos_5 <- c("FB023","FB028","FB048","FB036","FB007","FB082","FB090","FB019",
           "FB015","FB033","FB005","FB204","FB127","FB079","FB093","FB122",
           "FB196","FB348","FB308","FB230")
sum_data[match(pos_5,sum_data$ID),12] <- "yes"
## COVID-19 cases count
sum_data$Cohort <- "coh1"
sum_data[188:nrow(sum_data),"Cohort"] <- "coh2"
COV19_count <- subset(sum_data,PostV3BTI=="yes" | PostV4BTI=="yes" | PostV5BTI=="yes",select=c(V1,Cohort))
COV19_count_vac <- table(COV19_count$V1,COV19_count$Cohort)
knitr::kable(COV19_count_vac)
## Actual enrolled subjects
actual_vac <- data.frame(coh1 = c(NA, NA,NA), 
                         coh2 = c(NA, NA,NA), 
                         row.names = c("RQ3013", "RQ3025","RQ3027"))
actual_vac$coh1 <- c(63,123,NA)
actual_vac$coh2 <- c(63,NA,125)
actual_vac$sum <- c(126,123,125)
knitr::kable(actual_vac)
## Divided by actual
divided <- round(COV19_count_vac/actual_vac,4)
knitr::kable(divided)


### day 14 ninfec count
sum_data_B4_infec <- subset(sum_data, ID %in% 
                        c("FB147","FB211","FB299","FB277", "FB361",
                             "FB117" ,"FB320" , "FB342", "FB336", "FB335" ,
                             "FB374",
                             "FB198" ,"FB206","FB189", "FB295",
                             "FB169","FB349","FB369"))
table(sum_data_B4_infec$V1,sum_data_B4_infec$Cohort)
```

# Table1 Baseline char of participants in the ITT cohort
```{r Demographics and Baseline Characteristic}
library(data.table)
library(dplyr)
setwd("/DATA2/xuanjingyu/RQ3027 IIT/")
tb1 <- sum_data
infec <- fread("./infection.csv") %>% as.data.frame()
## Match "previous infection"
colnames(infec) <- c("Name","ID","Infection","Time")
tb1_1 <- full_join(tb1,infec,by="ID")
## Rename tb1_1 colnames
tb1_2 <- tb1_1 %>% dplyr::rename(Sex=Gender)
tb1_2 %>% filter(is.na(Infection))
tb1_2$`1207` <- ifelse(tb1_2$`1207` == "#N/A", "no", "yes")
tb1_2$Infection <- ifelse(tb1_2$Infection %in% c("A.是，记得确切日期____________",                                      "B.是，只记得模糊日期，请填写月份，例如2022-09____________"), "yes",
                          ifelse(tb1_2$Infection == "C.否", "no", tb1_2$Infection))
write_excel_csv(tb1_2,file = "tb1_2.csv")
## 0.1 digits
options(digits=10)
## Age mean/SD/median/IQR
tapply(tb1_2$Age,tb1_2$V1,function(x) quantile(x,c(.25,.75)))
## Sex n%
sex_n <- xtabs(~ V1+Sex, data = tb1_2)
prop.table(sex_n, 1)*100
## Count and classify BMI
BMI_c <- tb1_2 %>% 
  mutate(BMI = round(Weight/(Height*0.01)^2,1))
BMI_c$cla[BMI_c$BMI < 18.5] <- "underweight"
BMI_c$cla[BMI_c$BMI <= 24.9 & BMI_c$BMI >= 18.5] <- "healthy weight"
BMI_c$cla[BMI_c$BMI <= 29.9 & BMI_c$BMI >=25.0] <- "overweight"
BMI_c$cla[BMI_c$BMI >= 30.0] <- "obesity"
## BMI n%
BMI_cla_n <- xtabs(~ V1+cla, data = BMI_c)
## BMI cla mean
BMI_cla_n_mean <- BMI_c %>%
  group_by(V1, cla) %>%
  summarise(BMI_cla_n_mean = mean(BMI, na.rm = TRUE))
## Prior COV19 infection
pri_cov <- xtabs(~V1+Infection, data = tb1_2)


```

# Arrange 14 days AE data from EDC
```{r Arrange 14 days AE data}
##Download 14 Days AE diaries from EDC
library(tidyverse)
library(dplyr)
a <- dir("./edcAE/")
orign_AE <- fread(paste0("./edcAE/",a[1])) %>% as.data.frame()
for(i in 2: length(a)){
  data2 <- fread(paste0("./edcAE/",a[i])) %>% as.data.frame()
  orign_AE <- rbind(orign_AE,data2)
}
##Input unique ID
combined_AE <- fread("./conbined_AE.csv") %>% as.data.frame()
combined_AE %>% arrange(ID)
##Transfer Grade 5 <- 0
combined_AE[,5:15][combined_AE[,5:15] == 5] <- 0
##Select each sympotoms GradeMAX within 14 days
max_AE <- aggregate(combined_AE[,5:15],by = list(combined_AE$ID), FUN=max)
## Vlookup
max_14AE <- merge(max_AE,sum_data, by.x="Group.1", by.y="ID")
se_max_14AE <- max_14AE %>% select(Group.1,Pain,`Scleroma/Swelling`,Redness,`Acute allergic reaction`,`Nausea/Vomitting`,`Joint pain`,`Muscular Pain`,Headache,Chill,Fatigue,Fever, V1)
## Divide into Cohort1 and Cohort2
co1_AE <- se_max_14AE %>% head(n=187)
co2_AE <- se_max_14AE %>% tail(n=189)
## Convert into "longer"
co1_AE_longer <- co1_AE %>%
                    pivot_longer(2:12,
                             names_to="reactions",
                             values_to = "grade") %>% 
                    filter(grade != 0)
co2_AE_longer <- co2_AE %>% 
                   pivot_longer(2:12,
                                names_to = "reactions",
                                values_to = "grade") %>% 
                   filter(grade != 0)
```


# Fig.1.1 Cohort1 Adverse reactions within 14 days after vaccination in the ITT cohort
```{r Fig Cohort1 Adverse reactions within 14 days after vaccination}
library(ggplot2)
library(data.table)
library(tidyverse)
# Cohort 1 safety analysis
## Add new col "V1_total"
co1_AE_longer$V1_total <- co1_AE_longer$V1
co1_AE_longer[which(co1_AE_longer$V1 =="RQ3013"),5] <- "63"
co1_AE_longer[which(co1_AE_longer$V1 == "RQ3025"),5] <- "124"
## Count by V1,reactions and grade
nco1_AE_longer <- co1_AE_longer %>%
   group_by(V1, reactions, grade,V1_total) %>%
   count()  %>% 
   as.data.frame()
nco1_AE_longer$percent <- nco1_AE_longer$freq / as.numeric(nco1_AE_longer$V1_total)
## Set order
nco1_AE_longer$V1 <- factor(nco1_AE_longer$V1,levels = c("RQ3025","RQ3013"))
nco1_AE_longer$grade <- ifelse(nco1_AE_longer$grade == 1,"Grade 1",
                         ifelse(nco1_AE_longer$grade == 2,"Grade 2",
                         ifelse(nco1_AE_longer$grade == 3,"Grade 3",
                         nco1_AE_longer$grade)))
nco1_AE_longer$grade <- factor(nco1_AE_longer$grade,
                               levels = c("Grade 3","Grade 2","Grade 1"))
## Reactions str wrap and set order
nco1_AE_longer$reactions <- str_wrap(nco1_AE_longer$reactions,width = 15)
nco1_AE_longer$reactions <- factor(nco1_AE_longer$reactions,
                                   levels = c("Pain",
                                              "Scleroma/Swelling",
                                              "Redness",
                                              "Acute allergic\nreaction",
                                              "Nausea/Vomitting",
                                              "Joint pain",
                                              "Muscular Pain",
                                              "Headache",
                                              "Chill",
                                              "Fatigue",
                                              "Fever"))
## Using ggplot2
p_co1_AE <- ggplot(nco1_AE_longer, 
                   aes(V1, percent*100, 
                       fill=grade)) + geom_bar(stat="identity",
                                               width = 0.7)+
       scale_fill_manual(values = c('Grade 1'='#3A5D95',
                                    'Grade 2'='#35A77B',
                                    'Grade 3'='#E14133'),
                         breaks = c("Grade 3","Grade 2","Grade 1"))+
  facet_grid(.~reactions) 
p_co1_AE <- p_co1_AE + labs(y="Percentage (%)",
                         x="",
                         fill="")
p_co1_AE <- p_co1_AE + scale_y_continuous(expand = c(0,0))+ 
        theme(axis.text.x = element_text(size=10,
                                         angle = 30,
                                         vjust=0.5,
                                         color = "black"),
              axis.text.y = element_text(size=10,
                                         color="black",),
              axis.title = element_text(size=10),
              legend.position = 'top',
              legend.text=element_text(size=10),
              legend.title=element_text(size=10))
p_co1_AE<-p_co1_AE+theme(panel.grid.major =element_blank(),
             panel.grid.minor = element_blank(),
             panel.background = element_blank(),
             axis.line = element_line(colour = "black"),
             axis.ticks.y = element_blank())
p_co1_AE<-p_co1_AE+theme(strip.text=element_text(size=8),
                         panel.spacing.x = unit(0, "lines")
             )
p_co1_AE
```

# Fig 1.2 Cohort2 Adverse reactions within 14 days after vaccination in the ITT cohort 
```{r Fig Cohort1 Adverse reactions within 14 days after vaccination}

# Cohort 2 safety analysis
## Manually modify Grade 4
co2_AE_longer[39,4] <- 1
co2_AE_longer[68,4] <- 1
## Add new col "V1_total"
co2_AE_longer$V1_total <- co2_AE_longer$V1
co2_AE_longer[which(co2_AE_longer$V1 =="RQ3013"),5] <- "63"
co2_AE_longer[which(co2_AE_longer$V1 == "RQ3027"),5] <- "126"
## Count by V1,reactions and grade
nco2_AE_longer <- co2_AE_longer %>%
   group_by(V1, reactions, grade,V1_total) %>%
   count()  %>% 
   as.data.frame()
nco2_AE_longer$percent <- nco2_AE_longer$freq / as.numeric(nco2_AE_longer$V1_total)
## Set order
nco2_AE_longer$V1 <- factor(nco2_AE_longer$V1,levels = c("RQ3027","RQ3013"))
nco2_AE_longer$grade <- ifelse(nco2_AE_longer$grade == 1,"Grade 1",
                         ifelse(nco2_AE_longer$grade == 2,"Grade 2",
                         ifelse(nco2_AE_longer$grade == 3,"Grade 3",
                         nco2_AE_longer$grade)))
nco2_AE_longer$grade <- factor(nco2_AE_longer$grade,
                               levels = c("Grade 3","Grade 2","Grade 1"))
## Reactions str wrap and set order
nco2_AE_longer$reactions <- str_wrap(nco2_AE_longer$reactions,width = 15)
nco2_AE_longer$reactions <- factor(nco2_AE_longer$reactions,
                                   levels = c("Pain",
                                              "Scleroma/Swelling",
                                              "Redness",
                                              "Acute allergic\nreaction",
                                              "Nausea/Vomitting",
                                              "Joint pain",
                                              "Muscular Pain",
                                              "Headache",
                                              "Chill",
                                              "Fatigue",
                                              "Fever"))
## Using ggplot2
p_co2_AE <- ggplot(nco2_AE_longer, 
                   aes(V1, percent*100, 
                       fill=grade)) + geom_bar(stat="identity",
                                               width = 0.7)+
       scale_fill_manual(values = c('Grade 1'='#3A5D95',
                                    'Grade 2'='#35A77B',
                                    'Grade 3'='#E14133'),
                         breaks = c("Grade 3","Grade 2","Grade 1"))+
  facet_grid(.~reactions) 
p_co2_AE <- p_co2_AE + labs(y="Percentage (%)",
                         x="",
                         fill="")
p_co2_AE <- p_co2_AE + scale_y_continuous(expand = c(0,0))+ 
        theme(axis.text.x = element_text(size=10,
                                         angle = 30,
                                         vjust=0.5,
                                         color = "black"),
              axis.text.y = element_text(size=10,
                                         color="black",),
              axis.title = element_text(size=10),
              legend.position = 'top',
              legend.text=element_text(size=10),
              legend.title=element_text(size=10))
p_co2_AE<-p_co2_AE+theme(panel.grid.major =element_blank(),
             panel.grid.minor = element_blank(),
             panel.background = element_blank(),
             axis.line = element_line(colour = "black"),
             axis.ticks.y = element_blank())
p_co2_AE<-p_co2_AE+theme(strip.text=element_text(size=8),
                         panel.spacing.x = unit(0, "lines")
             )
p_co2_AE
```

# Fig.2 Immunogenicity of RQ3027, RQ3025 and RQ3013
## Overall immunogenicity analysis
``` {r Overall immunogenicity analysis}
library(patchwork)
library(Rmisc)
library(rstatix)
library(ggbeeswarm)
library(scales)
setwd("/DATA2/xuanjingyu/RQ3027 IIT/GMT")
## Using B1B4  data
B1 <- fread("./GMT/B1_0722.csv") %>% as.data.frame()
B4 <- fread("./GMT/B4_0722.csv") %>% as.data.frame()
B1B4 <- rbind(B1,B4)
## Del infected row
# mer_B1B4_sum <- inner_join(B1B4, sum_data, by = "ID")
# B1B4_ninfec <- mer_B1B4_sum %>% filter(PostV1BTI != "yes" & 
#                                          PostV2BTI != "yes" & 
#                                          PostV3BTI != "yes")
B1B4_ninfec <- B1B4[!B1B4$ID %in% c("FB147","FB211","FB299","FB277", "FB361",
                             "FB117" ,"FB320" , "FB342", "FB336", "FB335" ,
                             "FB374",
                             "FB198" ,"FB206","FB189", "FB295",
                             "FB169","FB349","FB369"),]
## Del non-relevant column
B1B4_clean <- B1B4_ninfec[,-c(6:23)]
## Rename
B1B4_clean$Patch <- ifelse(B1B4_clean$Patch == "B1","Day0",
                               ifelse(B1B4_clean$Patch == "B4", "Day14",
                                      B1B4_clean$Patch))
colnames(B1B4_clean)[colnames(B1B4_clean) == "SARS-CoV-2"] <- "VOC"
## Del V4 drop-out row
ov_B1B4 <- B1B4_clean %>% filter(ID != "FB030" &
                                          ID != "FB188"&
                                          ID != "FB291")

## Select XBB.1.16, BQ.1.1, BF.7,WT,BA.45 data
ov_B1B4_0703 <- subset(ov_B1B4, VOC %in% c("XBB.1.16","BQ.1.1","BF.7",
                                               "Prototype","BA.4/5"))
## GMT
ov_B1B4_0703_log <- ov_B1B4_0703
ov_B1B4_0703_log$Levels <- log(ov_B1B4_0703_log$Levels)
GMT_ov_B1B4<-summarySE(ov_B1B4_0703_log, measurevar="Levels",
                    groupvars=c("Patch","Group","VOC"),
                    na.rm = TRUE)
GMT_ov_B1B4$"GMT95%CI.n" <- paste0(round(exp(GMT_ov_B1B4$Levels),2)," (",paste(round(exp(GMT_ov_B1B4$Levels-GMT_ov_B1B4$ci),2),"-",paste(round(exp(GMT_ov_B1B4$Levels+GMT_ov_B1B4$ci),2)),";","n=",GMT_ov_B1B4$N),")") 
GMT_ov_B1B4$GMT <- exp(GMT_ov_B1B4$Levels)
GMT_ov_B1B4$lowerCI <- round(exp(GMT_ov_B1B4$Levels-GMT_ov_B1B4$ci),2)
GMT_ov_B1B4$upperCI <- round(exp(GMT_ov_B1B4$Levels+GMT_ov_B1B4$ci),2)

## GMFR
ov_B1B4_0703_di <- ov_B1B4_0703 %>% 
  pivot_wider(names_from = "Patch",
              values_from = "Levels") %>% 
  group_by(VOC,Group) %>% 
  mutate(GMFR = Day14/Day0) %>% 
  ungroup()
ov_B1B4_0703_di$GMFR <- log(ov_B1B4_0703_di$GMFR)
GMFR_ov_B1B4<-summarySE(ov_B1B4_0703_di, measurevar="GMFR",
                    groupvars=c("Group","VOC"),
                    na.rm = TRUE)
GMFR_ov_B1B4$"GMFR95%CI.n" <- paste0(round(exp(GMFR_ov_B1B4$GMFR),1),
                                    " (",
                                    paste(round(exp(GMFR_ov_B1B4$GMFR-GMFR_ov_B1B4$ci),1),
                                          "-",
                                          paste(round(exp(GMFR_ov_B1B4$GMFR+GMFR_ov_B1B4$ci),1)),
                                          ";",
                                          "n=",
                                          GMFR_ov_B1B4$N),
                                    ")") 
## GMR ov_B4_Prototype
ov_B4_Prototype <- subset(ov_B1B4_0703,ov_B1B4_0703$VOC=="Prototype" )
ov_B4_Prototype$log <- log(ov_B4_Prototype$Levels)
ov_B4_Prototype_wider <- reshape::cast(ov_B4_Prototype[,-c(4,5)],
                                     ID+Group~Patch)
fit_ov_B4_Prototype <- lm(Day14~Day0+Group,data=ov_B4_Prototype_wider)
summary(fit_ov_B4_Prototype)
round(exp(confint(fit_ov_B4_Prototype,'GroupRQ3025',level = 0.95)),2)
round(exp(confint(fit_ov_B4_Prototype,'GroupRQ3027',level = 0.95)),2)

## GMR ov_B4_BA.4/5
ov_B4_BA.45 <- subset(ov_B1B4_0703,ov_B1B4_0703$VOC=="BA.4/5" )
ov_B4_BA.45$log <- log(ov_B4_BA.45$Levels)
ov_B4_BA.45_wider <- reshape::cast(ov_B4_BA.45[,-c(4,5)],
                                     ID+Group~Patch)
fit_ov_B4_BA.45 <- lm(Day14~Day0+Group,data=ov_B4_BA.45_wider)
summary(fit_ov_B4_BA.45)
round(exp(confint(fit_ov_B4_BA.45,'GroupRQ3025',level = 0.95)),2)
round(exp(confint(fit_ov_B4_BA.45,'GroupRQ3027',level = 0.95)),2)

## GMR ov_B4_BF.7
ov_B4_BF.7 <- subset(ov_B1B4_0703,ov_B1B4_0703$VOC=="BF.7" )
ov_B4_BF.7$log <- log(ov_B4_BF.7$Levels)
ov_B4_BF.7_wider <- reshape::cast(ov_B4_BF.7[,-c(4,5)],
                                     ID+Group~Patch)
fit_ov_B4_BF.7 <- lm(Day14~Day0+Group,data=ov_B4_BF.7_wider)
summary(fit_ov_B4_BF.7)
round(exp(confint(fit_ov_B4_BF.7,'GroupRQ3025',level = 0.95)),1)
round(exp(confint(fit_ov_B4_BF.7,'GroupRQ3027',level = 0.95)),1)
## GMR ov_B4_BQ.1.1
ov_B4_BQ.1.1 <- subset(ov_B1B4_0703,ov_B1B4_0703$VOC=="BQ.1.1" )
ov_B4_BQ.1.1$log <- log(ov_B4_BQ.1.1$Levels)
ov_B4_BQ.1.1_wider <- reshape::cast(ov_B4_BQ.1.1[,-c(4,5)],
                                     ID+Group~Patch)
fit_ov_B4_BQ.1.1 <- lm(Day14~Day0+Group,data=ov_B4_BQ.1.1_wider)
summary(fit_ov_B4_BQ.1.1)
round(exp(confint(fit_ov_B4_BQ.1.1,'GroupRQ3025',level = 0.95)),1)
round(exp(confint(fit_ov_B4_BQ.1.1,'GroupRQ3027',level = 0.95)),1)

## GMR ov_B4_XBB.1.16
ov_B4_XBB.1.16 <- subset(ov_B1B4_0703,ov_B1B4_0703$VOC=="XBB.1.16" )
ov_B4_XBB.1.16$log <- log(ov_B4_XBB.1.16$Levels)
ov_B4_XBB.1.16_wider <- reshape::cast(ov_B4_XBB.1.16[,-c(4,5)],
                                     ID+Group~Patch)
fit_ov_B4_XBB.1.16 <- lm(Day14~Day0+Group,data=ov_B4_XBB.1.16_wider)
summary(fit_ov_B4_XBB.1.16)
round(exp(confint(fit_ov_B4_XBB.1.16,'GroupRQ3025',level = 0.95)),1)
round(exp(confint(fit_ov_B4_XBB.1.16,'GroupRQ3027',level = 0.95)),1)

## GMR 2725_B4_XBB.1.16
s_B4_XBB.1.16 <- subset(ov_B4_XBB.1.16,Group %in% c("RQ3027",
                                                    "RQ3025"))
s_B4_XBB.1.16$log <- log(s_B4_XBB.1.16$Levels)
s_B4_XBB.1.16_wider <- reshape::cast(s_B4_XBB.1.16[,-c(4,5)],
                                     ID+Group~Patch)
fit_s_B4_XBB.1.16 <- lm(Day14~Day0+Group,data=s_B4_XBB.1.16_wider)
summary(fit_s_B4_XBB.1.16)
round(exp(confint(fit_s_B4_XBB.1.16,'GroupRQ3027',level = 0.95)),1)

```

``` {r Overall B1B4 Immunogenicity plot}
# Levels
GMT_ov_B1B4$Patch <- factor(GMT_ov_B1B4$Patch,levels = c("Day0","Day14"))
GMT_ov_B1B4$Group <- factor(GMT_ov_B1B4$Group,levels = c("RQ3013",
                                                   "RQ3025",
                                                   "RQ3027"))
GMT_ov_B1B4$VOC <- factor(GMT_ov_B1B4$VOC,levels = c("Prototype",
                                               "BA.4/5",
                                               "BF.7",
                                               "BQ.1.1",
                                               "XBB.1.16"))
ov_B1B4_0703$VOC <- factor(ov_B1B4_0703$VOC, levels = c("Prototype",
                                               "BA.4/5",
                                               "BF.7",
                                               "BQ.1.1",
                                               "XBB.1.16"))
# ggplot
p_co_B1B4 <- ggplot(GMT_ov_B1B4, aes(x=Patch,y=GMT,color=VOC)) +
  geom_bar(stat="identity",position = position_dodge(),
           width=0.9,
           size=0.5,
           fill="white") +
  facet_grid(~Group) 
p_co_B1B4 <- p_co_B1B4 + geom_errorbar(aes(ymin=lowerCI,
                                               ymax=upperCI),
                                           width=.2,
                                           size=0.5,
                                           position = position_dodge(0.9))
p_co_B1B4 <- p_co_B1B4 + geom_beeswarm(data = ov_B1B4_0703,
                                           dodge.width = 0.9,
                                           aes(y = Levels,
                                               x= Patch,
                                               color=VOC),
                                           cex = 0.7,
                                           size = 1.4,
                                           alpha=.3,
                                           shape=16,
                                           show.legend = FALSE,)
p_co_B1B4 <- p_co_B1B4 + theme_classic() +
  scale_color_manual(values = c("Prototype" = "#09357A",
                                "BA.4/5" = "#02915A",
                                "BF.7" = "#740D91",
                                "BQ.1.1" = "#FFB240",
                                "XBB.1.16" ="#DA1100")) +
  theme(legend.position = 'top',
        legend.text = element_text(size = 12)) +
  theme(axis.text = element_text(size = 12, color="black"),
        axis.line.y = element_line(color = 'black'),
        axis.title.y = element_text(size = 12, color="black")) +
  theme(strip.text.x = element_text(size = 12),
        strip.switch.pad.grid = unit(10,"lines"),
        strip.background.x = element_rect(color = "white"),
        panel.spacing.x = unit(0.3,"cm"),
        panel.grid = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(margin = margin(t=5,
                                                   r=0,
                                                   b=0,
                                                   l=0)))
p_co_B1B4 <- p_co_B1B4 + scale_y_continuous(
                     trans = log10_trans(),
                     breaks =trans_breaks("log10", function(x) 10^x),
                     labels =trans_format("log10", math_format(10^.x)),
                     expand = c(0,0),
                     limits = c(10^0,10^6)) + 
  labs(title = "",x = "",y = "Neutralizing Antibody ID50 \nGMT (95%CI)")
p_co_B1B4

```


## Data management of immunogenicity of Cohort1
``` {r Immunogenicity analysis Coh1}
library(tidyverse)
library(ggplot2)
library(dplyr)
library(patchwork)
library(rstatix)
setwd("/DATA2/xuanjingyu/RQ3027 IIT/GMT")
## Using B1B4  data
B1 <- fread("./GMT/B1_0717.csv") %>% as.data.frame()
B4 <- fread("./GMT/B4_0717.csv") %>% as.data.frame()
B1B4 <- rbind(B1,B4)
## Del infected row
mer_B1B4_sum <- inner_join(B1B4, sum_data, by = "ID")
B1B4_ninfec <- mer_B1B4_sum %>% filter(PostV1BTI != "yes" & 
                                         PostV2BTI != "yes" & 
                                         PostV3BTI != "yes")
# B1B4_infec <- mer_B1B4_sum %>% filter(PostV1BTI == "yes", 
#                                          PostV2BTI == "yes", 
#                                          PostV3BTI == "yes")
## Del non-relevant column
B1B4_clean <- B1B4_ninfec[,-c(6:23)]
## Rename
B1B4_clean$Patch <- ifelse(B1B4_clean$Patch == "B1","Day0",
                               ifelse(B1B4_clean$Patch == "B4", "Day14",
                                      B1B4_clean$Patch))
colnames(B1B4_clean)[colnames(B1B4_clean) == "SARS-CoV-2"] <- "VOC"
## Select Cohort1 FB001-FB187
coh1_B1B4_wdrop <- B1B4_clean[B1B4_clean$ID %in% 
                           paste0("FB",sprintf("%03d",1:187)),]
## Del V4 drop-out row
coh1_B1B4 <- coh1_B1B4_wdrop %>% filter(ID != "FB030")


## Select XBB.1.16, BQ.1.1, BF.7,WT,BA.45 data
coh1_B1B4_0703 <- subset(coh1_B1B4, VOC %in% c("XBB.1.16","BQ.1.1","BF.7",
                                               "Prototype","BA.4/5"))


```

## Cohort1 GMT
```{r GMT & GMR}
library(epiDisplay)
library(Rmisc)
library(reshape2)
## GMT
coh1_B1B4_0703 -> B.data
B.data$Levels <- log(B.data$Levels)
GMT.Wild<-summarySE(B.data, measurevar="Levels",
                    groupvars=c("Patch","Group","VOC"),
                    na.rm = TRUE)
GMT.Wild$"GMT95%CI.n" <- paste0(round(exp(GMT.Wild$Levels),2)," (",paste(round(exp(GMT.Wild$Levels-GMT.Wild$ci),2),"-",paste(round(exp(GMT.Wild$Levels+GMT.Wild$ci),2)),";","n=",GMT.Wild$N),")") 
GMT.Wild$GMT <- exp(GMT.Wild$Levels)
GMT.Wild$lowerCI <- round(exp(GMT.Wild$Levels-GMT.Wild$ci),2)
GMT.Wild$upperCI <- round(exp(GMT.Wild$Levels+GMT.Wild$ci),2)

## GMR co1_B4_BQ.1.1
co1_B1B4_BQ.1.1 <- subset(coh1_B1B4_0703,coh1_B1B4_0703$VOC=="BQ.1.1" )
co1_B1B4_BQ.1.1$log <- log(co1_B1B4_BQ.1.1$Levels)
co1_B1B4_BQ.1.1_wider <- reshape::cast(co1_B1B4_BQ.1.1[,-c(4,5)],
                                     ID+Group~Patch)
fit_co1_B4_BQ.1.1 <- lm(Day14~Day0+Group,data=co1_B1B4_BQ.1.1_wider)
summary(fit_co1_B4_BQ.1.1)
round(exp(confint(fit_co1_B4_BQ.1.1,'GroupRQ3025',level = 0.95)),2)
## GMR co1_B4_XBB.1.16
co1_B1B4_XBB.1.16 <- subset(coh1_B1B4_0703,coh1_B1B4_0703$VOC=="XBB.1.16" )
co1_B1B4_XBB.1.16$log <- log(co1_B1B4_XBB.1.16$Levels)
co1_B1B4_XBB.1.16_wider <- reshape::cast(co1_B1B4_XBB.1.16[,-c(4,5)],
                                     ID+Group~Patch)
fit_co1_B4_XBB.1.16 <- lm(Day14~Day0+Group,data=co1_B1B4_XBB.1.16_wider)
summary(fit_co1_B4_XBB.1.16)
round(exp(confint(fit_co1_B4_XBB.1.16,'GroupRQ3025',level = 0.95)),2)
## GMR co1_B4_BF.7
co1_B1B4_BF.7 <- subset(coh1_B1B4_0703,coh1_B1B4_0703$VOC=="BF.7" )
co1_B1B4_BF.7$log <- log(co1_B1B4_BF.7$Levels)
co1_B1B4_BF.7_wider <- reshape::cast(co1_B1B4_BF.7[,-c(4,5)],
                                     ID+Group~Patch)
fit_co1_B4_BF.7 <- lm(Day14~Day0+Group,data=co1_B1B4_BF.7_wider)
summary(fit_co1_B4_BF.7)
round(exp(confint(fit_co1_B4_BF.7,'GroupRQ3025',level = 0.95)),2)
```

## Barplot of Cohort1 immunogenicity 
```{r barplot}

library(ggplot2)
library(ggbeeswarm)
library(scales)
# Levels
GMT.Wild$Patch <- factor(GMT.Wild$Patch,levels = c("Day0","Day14"))
GMT.Wild$Group <- factor(GMT.Wild$Group,levels = c("RQ3013",
                                                   "RQ3025"))
GMT.Wild$VOC <- factor(GMT.Wild$VOC,levels = c("Prototype",
                                               "BA.4/5",
                                               "BF.7",
                                               "BQ.1.1",
                                               "XBB.1.16"))
coh1_B1B4_0703$VOC <- factor(coh1_B1B4_0703$VOC, levels = c("Prototype",
                                                            "BA.4/5",
                                                            "BF.7",
                                                            "BQ.1.1",
                                                            "XBB.1.16"))
# ggplot
p_coh1_B1B4 <- ggplot(GMT.Wild, aes(x=Patch,y=GMT,color=VOC)) +
  geom_bar(stat="identity",position = position_dodge(),
           width=0.9,
           size=0.5,
           fill="white") +
  facet_grid(~Group) 
p_coh1_B1B4 <- p_coh1_B1B4 + geom_errorbar(aes(ymin=lowerCI,
                                               ymax=upperCI),
                                           width=.2,
                                           size=0.5,
                                           position = position_dodge(0.9))
p_coh1_B1B4 <- p_coh1_B1B4 + geom_beeswarm(data = coh1_B1B4_0703,
                                           dodge.width = 0.9,
                                           aes(y = Levels,
                                               x = Patch,
                                               color = VOC),
                                           cex = 0.8,
                                           size = 1.4,
                                           alpha = .3,
                                           shape = 16,
                                           show.legend = FALSE)
p_coh1_B1B4 <- p_coh1_B1B4 + theme_classic() +
  scale_color_manual(values = c("Prototype" = "#09357A",
                                "BA.4/5" = "#02915A",
                                "BF.7" = "#740D91",
                                "BQ.1.1" = "#FFB240",
                                "XBB.1.16" ="#DA1100")) +
  theme(legend.position = 'top',
        legend.text = element_text(size = 12),
        legend.spacing.x = unit(0.3,"cm")) +
  theme(axis.text = element_text(size = 12, color="black"),
        axis.line.y = element_line(color = 'black'),
        axis.title.y = element_text(size = 12, color="black")) +
  theme(strip.text.x = element_text(size = 12),
        strip.switch.pad.grid = unit(10,"lines"),
        strip.background.x = element_rect(color = "white"),
        panel.spacing.x = unit(0.5,"cm"),
        panel.grid = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(margin = margin(t=5,
                                                   r=0,
                                                   b=0,
                                                   l=0)))
p_coh1_B1B4 <- p_coh1_B1B4 + scale_y_continuous(
                     trans = log10_trans(),
                     breaks =trans_breaks("log10", function(x) 10^x),
                     labels =trans_format("log10", math_format(10^.x)),
                     expand = c(0,0),
                     limits = c(10^0,10^5)) + 
  labs(title = "",x = "",y = "Neutralizing Antibody ID50 \nGMT (95%CI)")
p_coh1_B1B4


```

## Data management of immunogenicity of Cohort2
``` {r Immunogenicity analysis Coh2}
library(tidyverse)
library(ggplot2)
library(dplyr)
library(patchwork)
library(rstatix)
setwd("/DATA2/xuanjingyu/RQ3027 IIT/GMT")

## Select Cohort2 FB188-FB376
coh2_B1B4_wdrop <- B1B4_clean[B1B4_clean$ID %in% 
                           paste0("FB",sprintf("%03d",188:376)),]
## Del V4 drop-out row
coh2_B1B4 <- coh2_B1B4_wdrop %>% filter(ID != "FB188"&
                                          ID != "FB291")


## Select XBB.1.16, BQ.1.1, BF.7 data
coh2_B1B4_0703 <- subset(coh2_B1B4, VOC %in% c("XBB.1.16","BQ.1.1","BF.7",
                                               "Prototype","BA.4/5"))

```
## Cohort2 GMT
```{r }
library(epiDisplay)
library(Rmisc)
library(reshape2)
## GMT
coh2_B1B4_0703_log <- coh2_B1B4_0703
coh2_B1B4_0703_log$Levels <- log(coh2_B1B4_0703_log$Levels)
GMT_coh2_B1B4<-summarySE(coh2_B1B4_0703_log, measurevar="Levels",
                    groupvars=c("Patch","Group","VOC"),
                    na.rm = TRUE)
GMT_coh2_B1B4$"GMT95%CI.n" <- paste0(round(exp(GMT_coh2_B1B4$Levels),2)," (",paste(round(exp(GMT_coh2_B1B4$Levels-GMT_coh2_B1B4$ci),2),"-",paste(round(exp(GMT_coh2_B1B4$Levels+GMT_coh2_B1B4$ci),2)),";","n=",GMT_coh2_B1B4$N),")") 
GMT_coh2_B1B4$GMT <- exp(GMT_coh2_B1B4$Levels)
GMT_coh2_B1B4$lowerCI <- round(exp(GMT_coh2_B1B4$Levels-GMT_coh2_B1B4$ci),2)
GMT_coh2_B1B4$upperCI <- round(exp(GMT_coh2_B1B4$Levels+GMT_coh2_B1B4$ci),2)

## GMR co1_B4_BQ.1.1
co2_B1B4_BQ.1.1 <- subset(coh2_B1B4_0703,coh2_B1B4_0703$VOC=="BQ.1.1" )
co2_B1B4_BQ.1.1$log <- log(co2_B1B4_BQ.1.1$Levels)
co2_B1B4_BQ.1.1_wider <- reshape::cast(co2_B1B4_BQ.1.1[,-c(4,5)],
                                     ID+Group~Patch)
fit_co2_B4_BQ.1.1 <- lm(Day14~Day0+Group,data=co2_B1B4_BQ.1.1_wider)
summary(fit_co2_B4_BQ.1.1)
round(exp(confint(fit_co2_B4_BQ.1.1,'GroupRQ3027',level = 0.95)),2)
## GMR co1_B4_XBB.1.16
co2_B1B4_XBB.1.16 <- subset(coh2_B1B4_0703,coh2_B1B4_0703$VOC=="XBB.1.16" )
co2_B1B4_XBB.1.16$log <- log(co2_B1B4_XBB.1.16$Levels)
co2_B1B4_XBB.1.16_wider <- reshape::cast(co2_B1B4_XBB.1.16[,-c(4,5)],
                                     ID+Group~Patch)
fit_co2_B4_XBB.1.16 <- lm(Day14~Day0+Group,data=co2_B1B4_XBB.1.16_wider)
summary(fit_co2_B4_XBB.1.16)
round(exp(confint(fit_co2_B4_XBB.1.16,'GroupRQ3027',level = 0.95)),2)
## GMR co1_B4_BF.7
co2_B1B4_BF.7 <- subset(coh2_B1B4_0703,coh2_B1B4_0703$VOC=="BF.7" )
co2_B1B4_BF.7$log <- log(co2_B1B4_BF.7$Levels)
co2_B1B4_BF.7_wider <- reshape::cast(co2_B1B4_BF.7[,-c(4,5)],
                                     ID+Group~Patch)
fit_co2_B4_BF.7 <- lm(Day14~Day0+Group,data=co2_B1B4_BF.7_wider)
summary(fit_co2_B4_BF.7)
round(exp(confint(fit_co2_B4_BF.7,'GroupRQ3027',level = 0.95)),2)

```
## Barplot of Cohort2 immunogenicity 
```{r barplot}

library(ggplot2)
library(ggbeeswarm)
library(scales)

# Levels
GMT_coh2_B1B4$Patch <- factor(GMT_coh2_B1B4$Patch,levels = c("Day0","Day14"))
GMT_coh2_B1B4$Group <- factor(GMT_coh2_B1B4$Group,levels = c("RQ3013",
                                                   "RQ3027"))
GMT_coh2_B1B4$VOC <- factor(GMT_coh2_B1B4$VOC,levels = c("Prototype",
                                               "BA.4/5",
                                               "BF.7",
                                               "BQ.1.1",
                                               "XBB.1.16"))
coh2_B1B4_0703$VOC <- factor(coh2_B1B4_0703$VOC, levels = c("Prototype",
                                               "BA.4/5",
                                               "BF.7",
                                               "BQ.1.1",
                                               "XBB.1.16"))
# ggplot
p_coh2_B1B4 <- ggplot(GMT_coh2_B1B4, aes(x=Patch,y=GMT,color=VOC)) +
  geom_bar(stat="identity",position = position_dodge(),
           width=0.9,
           size=0.5,
           fill="white") +
  facet_grid(~Group) 
p_coh2_B1B4 <- p_coh2_B1B4 + geom_errorbar(aes(ymin=lowerCI,
                                               ymax=upperCI),
                                           width=.2,
                                           size=0.5,
                                           position = position_dodge(0.9))
p_coh2_B1B4 <- p_coh2_B1B4 + geom_beeswarm(data = coh2_B1B4_0703,
                                           dodge.width = 0.9,
                                           aes(y = Levels,
                                               x= Patch,
                                               color=VOC),
                                           cex = 0.8,
                                           size = 1.4,
                                           alpha=.3,
                                           shape=16,
                                           show.legend = FALSE,)
p_coh2_B1B4 <- p_coh2_B1B4 + theme_classic() +
  scale_color_manual(values = c("Prototype" = "#09357A",
                                "BA.4/5" = "#02915A",
                                "BF.7" = "#740D91",
                                "BQ.1.1" = "#FFB240",
                                "XBB.1.16" ="#DA1100")) +
  theme(legend.position = 'top',
        legend.text = element_text(size = 12)) +
  theme(axis.text = element_text(size = 12, color="black"),
        axis.line.y = element_line(color = 'black'),
        axis.title.y = element_text(size = 12, color="black")) +
  theme(strip.text.x = element_text(size = 12),
        strip.switch.pad.grid = unit(10,"lines"),
        strip.background.x = element_rect(color = "white"),
        panel.spacing.x = unit(0.3,"cm"),
        panel.grid = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(margin = margin(t=5,
                                                   r=0,
                                                   b=0,
                                                   l=0)))
p_coh2_B1B4 <- p_coh2_B1B4 + scale_y_continuous(
                     trans = log10_trans(),
                     breaks =trans_breaks("log10", function(x) 10^x),
                     labels =trans_format("log10", math_format(10^.x)),
                     expand = c(0,0),
                     limits = c(10^0,10^5)) + 
  labs(title = "",x = "",y = "Neutralizing Antibody ID50 \nGMT (95%CI)")
p_coh2_B1B4


```



## Cumulative event rate of COVID-19
```{r COVID-19 cases}
library(ggplot2)
library(survival)
library(survminer)
## Manage confirmed COVID-19 cases
cul <- sum_data[,c("ID", "V1", "PostV1BTI", "PostV2BTI", 
                   "PostV3BTI", "PostV4BTI", "PostV5BTI", 
                   "PostV6BTI", "PostV7BTI", "Cohort")]
## Infection day
cul$inf_day <- 70
cul$inf_day[cul$ID == "FB147" | cul$ID == "FB211"] <- c(0,0)
cul$inf_day[cul$ID == "FB299" | cul$ID == "FB277" |
              cul$ID == "FB361"] <- c(1,1,1)
cul$inf_day[cul$ID == "FB117" | cul$ID == "FB320" |
              cul$ID == "FB342" | cul$ID == "FB336" |
              cul$ID == "FB335"] <- 
  c(7,7,7,7,7)
cul$inf_day[cul$ID == "FB374"] <- 13
cul$inf_day[cul$ID == "FB198" | cul$ID == "FB206" |
              cul$ID == "FB189" | cul$ID == "FB295" | 
              cul$ID == "FB169" | cul$ID == "FB349" |
              cul$ID == "FB369"] <- 
  c(14,14,14,14,14,14,14)
cul$inf_day[ cul$ID == "FB027" | cul$ID == "FB222"] <- c(16,16)
cul$inf_day[cul$ID == "FB066"] <- 18
cul$inf_day[cul$ID == "FB106" |  cul$ID == "FB247"] <- c(20,20)
cul$inf_day[cul$ID== "FB105"] <- 22
cul$inf_day[cul$ID == "FB372"] <- 23
cul$inf_day[cul$ID == "FB064" | cul$ID == "FB070" | 
              cul$ID == "FB113" | cul$ID == "FB116" | 
              cul$ID == "FB120" | cul$ID == "FB126" | 
              cul$ID == "FB135" | cul$ID == "FB150" | 
              cul$ID == "FB151" | cul$ID == "FB175" | 
              cul$ID == "FB273" | cul$ID == "FB291" | 
              cul$ID == "FB314" | cul$ID == "FB331" | 
              cul$ID == "FB329" | cul$ID == "FB332" ] <-
  c(28,28,28,28,28,28,28,28,28,28,28,28,28,28,28)
cul$inf_day[cul$ID == "FB023" | cul$ID == "FB028" |
              cul$ID == "FB048" | cul$ID == "FB036" | 
              cul$ID == "FB033"] <-
  c(35,35,35,35,35)
cul$inf_day[cul$ID == "FB019" | cul$ID == "FB007" |
              cul$ID == "FB082" | cul$ID == "FB090" ] <- 
  c(36,36,36,36)
cul$inf_day[cul$ID == "FB204"] <- 37
cul$inf_day[cul$ID == "FB015" |  cul$ID == "FB348"] <- c(38,38)
cul$inf_day[cul$ID == "FB005"] <- 47
cul$inf_day[cul$ID == "FB127"] <- 49
cul$inf_day[cul$ID == "FB079" | cul$ID == "FB093" |
              cul$ID == "FB122" | cul$ID == "FB196"] <-
  c(51,51,51,51)
cul$inf_day[cul$ID == "FB308"] <- 61
cul$inf_day[cul$ID == "FB230"] <- 70
## Status
cul$status <- 0
infected <- c("FB198","FB206","FB189","FB295","FB169","FB374",
              "FB349","FB369","FB066","FB106","FB027","FB222","FB064","FB070",
              "FB113","FB116","FB120","FB126","FB135","FB150","FB151","FB175",
              "FB105","FB273","FB291","FB247","FB372","FB314","FB331","FB329",
              "FB332","FB023","FB028","FB048","FB036","FB007","FB082","FB090",
              "FB019","FB015","FB033","FB005","FB204","FB127","FB079","FB093",
              "FB122","FB196","FB348","FB308")
cul[match(infected,cul$ID),12] <- 1
## Del protocol deviation
cul_0714 <- cul %>% filter(ID != "FB147" &
                             ID != "FB211")
## plot
f <- surv_fit(Surv(inf_day,status)~V1,data=cul_0714)
summary(f)
p_cul_0714 <- ggsurvplot(f, fun = "event",
           # conf.int = TRUE, 
           # conf.int.style="step",
           # conf.int.alpha=0.45, 
           surv.plot.height= 0.7,
           # linetype = "strata",
           legend.title = "Group",
           legend.labs = c("RQ3013", "RQ3025","RQ3027"),
           legend=c(0.2,0.8),
           xlab="Days",
           ylab="Cumulative Event Rate (%)",
           xlim = c(0,70), ylim = c(0,0.2),
           risk.table = TRUE,
           # tables.height = 0.15,
           # #tables.theme = theme_cleantable(),
           font.x = c(12, "plain", "black"),
           break.x.by = 7,
           font.y = c(12, "plain", "black"),
           font.tickslab = c(12, "plain", "black"),
           # risk.table.col="strata",
           # risk.table.height=0.2,
           palette = c("#09357A","#02915A","#740D91")
          )
p_cul_0714
```

## ZZJ GMT check
```{r GMT 95%CI ZZJ check}
geom_mean <- function(x, CI = 0.95){
  x=x[!is.na(x)]
  n = length(x)
  X = log(x)
  geomM <- exp(mean( X ) )
  a = (1-CI)/2
  Int = abs( qt(a, n-1)*sd(X)/sqrt(n)  )
  re = c( round( geomM , digits = 3), 
          paste( round( exp(mean(X)- Int) , 
          digits = 3), 
          round( exp(mean(X) + Int ) ,                                    
                 digits = 3), sep = "," ) )
  names(re) = c("GeometricalMean","CI_interval")
  return(re)
}

ov_B4_BQ <- ov_B4_BQ.1.1 %>%
  filter(ov_B4_BQ.1.1$Patch == "Day14")
GMT_ov_B4_BQ <- tapply(ov_B4_BQ$Levels,
                       ov_B4_BQ$Group,
                       geom_mean)

```