```{r N-IgG}
library(ggplot2)
library(tidyverse)

N_IgG <- fread("/DATA2/xuanjingyu/RQ3027 IIT/Anti_N_IgG.csv") %>% as.data.frame()
p_N_IgG <- ggplot(N_IgG, aes(x=Concentration)) +
  geom_histogram()
  # geom_histogram(mapping = aes(y = after_stat(density))) +
  # geom_density(color = "red", linewidth = 1)

p_ECDF <- ggplot(N_IgG, aes(x = Concentration,
                            y =after_stat(y)))  +
  stat_ecdf(geom = "step")

p_vac_NIGG <- ggplot(N_IgG, aes(x = Concentration,
                                color = baseline,
                                fill = baseline)) +
  geom_density(alpha = .4)


```

```{r Correlation between Day 14 WT NAbs and N-IgG}

B4_WT <- subset(ov_B4_Prototype, ov_B4_Prototype$Patch=="Day14")
merge_N_WT_B4 <- merge(N_IgG, B4_WT, by = "ID")
## Correlation test
N_WT_B4 <- merge_N_WT_B4[,c(2,4,7)]
cor(N_WT_B4)
cor.test(N_WT_B4[,1],N_WT_B4[,2])
## plot
library(ggpubr)
p1 <- ggplot(N_WT_B4, aes(x = Concentration,
                          y = Levels,
                          color = baseline)) +
  geom_point(aes(color = baseline)) +
  geom_smooth(method = 'lm', formula = y~x, se = T) +
  stat_cor(data = N_WT_B4, method = "pearson")


```


```{r Correlation between Day 0 WT NAbs and N-IgG}

B1_WT <- subset(ov_B4_Prototype, ov_B4_Prototype$Patch=="Day0")
merge_N_WT_B1 <- merge(N_IgG, B1_WT, by = "ID")
## Correlation test
N_WT_B1 <- merge_N_WT_B1[,c(2,4,7)]

## plot
library(ggpubr)
p2 <- ggplot(N_WT_B1, aes(x = Concentration,
                          y = Levels,
                          color = baseline)) +
  geom_point(aes(color = baseline)) +
  geom_smooth(method = 'lm', formula = y~x, se = T) +
  stat_cor(data = N_WT_B1, method = "pearson")


```