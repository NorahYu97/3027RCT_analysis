```{r Error in function? tidyverse column_to_rownames}
library(Racmacs)
library(dplyr)
library(tidyverse)
library(data.table)
map_14 <- subset(cul,inf_day %in% c("1","7","13","14"))
B4 <- fread("/DATA2/xuanjingyu/RQ3027 IIT/GMT/B4_0722.csv") %>% 
  as.data.frame()
colnames(B4)[colnames(B4) == "SARS-CoV-2"] <- "VOC"
sub_map14 <- subset(B4, ID %in% map_14$ID)
NAb_map14 <- subset(sub_map14, VOC %in% c("XBB.1.16",
                        "BQ.1.1",
                        "BF.7",
                        "Prototype",
                        "BA.4/5"))
NAb_map14_w <- NAb_map14 %>% pivot_wider(id_cols = "VOC",
                                         names_from = "ID",
                                         values_from = "Levels")
# NAb_map14c <- column_to_rownames(NAb_map14_w, var = "VOC")
# dt_14 <- apply(NAb_map14c, 2, function(x) max(x) - x)

map <- acmap(titer_table = NAb_map14c)
map1 <- optimizeMap(map = map,
                   number_of_dimensions = 2,
                   number_of_optimizations = 500,
                   minimum_column_basis = "none")
view(map1)


map3d <- make.acmap(titer_table = NAb_map14c,
                    number_of_dimensions = 3,
                    number_of_optimizations = 500,
                    minimum_column_basis = 'none')
view(map3d)
```

```{r }
library(Racmacs)
library(dplyr)
library(data.table)
test <- fread("/home/xuanjingyu/2.csv") %>% as.data.frame()
test <- test[,-1]
voc <- c("BQ.1.1",
         "XBB.1.16",
         "BF.7",
         "BA.4/5",
         "Prototype")
rownames(test)<-voc

map <- acmap(titer_table = test)
map <- optimizeMap(map=map,
                   number_of_dimensions = 2,
                   number_of_optimizations = 500,
                   minimum_column_basis = "none",)
view(map)


map3d <- make.acmap(titer_table = test,
                    number_of_dimensions = 3,
                    number_of_optimizations = 500,
                    minimum_column_basis = 'none')
view(map3d)
```