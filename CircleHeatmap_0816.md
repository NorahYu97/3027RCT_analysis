---
title: "Untitled"
author: "Norah Yu"
date: '2023-08-16'
output: html_document
---

```{r Circle Heatmap of B1 NAbs}
library(tidyverse)
library(data.table)
library(dplyr)
library('ComplexHeatmap')
library('circlize')
library("RColorBrewer")
library(dendextend)
library(gridBase)

B1 <- fread("/DATA2/xuanjingyu/RQ3027 IIT/GMT/B1_0816.csv") %>% as.data.frame()
B4 <- fread("/DATA2/xuanjingyu/RQ3027 IIT/GMT/B4_0814.csv") %>% as.data.frame()

B1B4 <- rbind(B1,B4)

## Del infected row
B1B4_ninfec <- B1B4[!B1B4$ID %in% c("FB147","FB211","FB299","FB277", "FB361",
                             "FB117" ,"FB320" , "FB342", "FB336", "FB335" ,
                             "FB374",
                             "FB198" ,"FB206","FB189", "FB295",
                             "FB169","FB349","FB369"),]

## Rename
B1B4_ninfec$Patch <- ifelse(B1B4_ninfec$Patch == "B1","Day0",
                               ifelse(B1B4_ninfec$Patch == "B4", "Day14",
                                      B1B4_ninfec$Patch))
colnames(B1B4_ninfec)[colnames(B1B4_ninfec) == "SARS-CoV-2"] <- "VOC"

## Del V4 drop-out row
B1B4_ndrop <- B1B4_ninfec %>% filter(ID != "FB030" &
                                          ID != "FB188"&
                                          ID != "FB291") 

## Select XBB.1.16, BQ.1.1, BF.7,WT,BA.45 data
B1B4_variants <- subset(B1B4_ndrop, VOC %in% c("XBB.1.16",
                                           "XBB.1.9.1",
                                           "XBB.1.5",
                                           "BQ.1.1",
                                           "BF.7",
                                           "BA.4/5",
                                           "Beta",
                                           "Alpha",
                                           "Prototype"))

B1B4_wider <- B1B4_variants %>% pivot_wider(id_cols = c("ID",
                                                    "Group",
                                                    "Patch"),
                                        names_from = "VOC",
                                        values_from = "Levels")
B1B4_wider$newID <- paste0(B1B4_wider$ID,"_",B1B4_wider$Patch)

B1B4_wider_new <- column_to_rownames(B1B4_wider, var = "ID")
B1B4_w <- B1B4_wider_new[,-c(1)]

# g_rowB1 <- B1_wider[,3]
# g_rowB1 <- as.matrix(g_rowB1)

cir_B1B4 <- as.matrix(B1B4_w)
voc_order <- c("XBB.1.9.1",
               "XBB.1.16",
               "XBB.1.5",
               "BQ.1.1",
               "BF.7",
               "BA.4/5",
               "Beta",
               "Alpha",
               "Prototype")

## Normalize data
cir_B1B4 <- (scale(log2(B1B4_wider)))
cir_B1 <- cir_B1[,order(factor(colnames(cir_B1),
                               levels = voc_order))]

range(cir_B1)

## Color
mycol=colorRamp2(c(-2, 0, 2),
                 c("#58A4FD","white","#5E1D9D"))
 
mycol2=colorRamp2(c(-5, 0, 5),
                 c("blue","white","red"))

p1 <- Heatmap(cir_B1,
              row_names_gp = gpar(fontsize = 4),
        column_names_gp= gpar(fontsize = 8),
        col= mycol2,
        name="legend",
        cluster_columns = F,
        cluster_rows = F)

library(pheatmap)
pheatmap(cir_B1)




```

# B1 3027 circle heatmap NAbs
```{r Circle Heatmap of B1 NAbs}
library(tidyverse)
library(data.table)
library(dplyr)
library('ComplexHeatmap')
library('circlize')
library("RColorBrewer")
library(dendextend)
library(gridBase)

B1 <- fread("/DATA2/xuanjingyu/RQ3027 IIT/GMT/B1_0816.csv") %>% as.data.frame()

## Del infected row
B1_ninfec <- B1[!B1$ID %in% c("FB147","FB211","FB299","FB277", "FB361",
                             "FB117" ,"FB320" , "FB342", "FB336", "FB335" ,
                             "FB374",
                             "FB198" ,"FB206","FB189", "FB295",
                             "FB169","FB349","FB369"),]

## Rename
B1_ninfec$Patch <- ifelse(B1_ninfec$Patch == "B1","Day0",
                               ifelse(B1_ninfec$Patch == "B4", "Day14",
                                      B1_ninfec$Patch))
colnames(B1_ninfec)[colnames(B1_ninfec) == "SARS-CoV-2"] <- "VOC"

## Del V4 drop-out row
B1_ndrop <- B1_ninfec %>% filter(ID != "FB030" &
                                          ID != "FB188"&
                                          ID != "FB291") 

## Select XBB.1.16, BQ.1.1, BF.7,WT,BA.45 data
B1_variants <- subset(B1_ndrop, VOC %in% c("Prototype",
               "Alpha",
               "Beta",
               "BA.4/5",
               "BF.7",
               "BQ.1.1",
               "XBB.1.5",
               "XBB.1.16",
               "XBB.1.9.1"))
# B1_variants$Levels <- log(B1_variants$Levels)
B1_3027 <- subset(B1_variants, Group == "RQ3027")

B1_27w <- B1_3027 %>% pivot_wider(id_cols = c("ID"),
                                        names_from = "VOC",
                                        values_from = "Levels")
B1_27w_new <- column_to_rownames(B1_27w, var = "ID")


ma_B1_27 <- as.matrix(B1_27w_new)
voc_order <- c("XBB.1.9.1",
               "XBB.1.16",
               "XBB.1.5",
               "BQ.1.1",
               "BF.7",
               "BA.4/5",
               "Beta",
               "Alpha",
               "Prototype")

## Normalize data
cir_B1_27 <- t(scale(t(ma_B1_27)))
pheatmap(cir_B1_27)

cir_B1_27 <- cir_B1_27[,order(factor(colnames(cir_B1_27), 
                                   levels = voc_order))]

range(cir_B1_27)

## Color
mycol=colorRamp2(c(-2, 0, 2),
                 c("#58A4FD","white","#5E1D9D"))
mycol2=colorRamp2(c(-2, 0, 2),
                 c("blue","white","red")) 


Heatmap(cir_B1_27,
        row_names_gp = gpar(fontsize = 7),
        column_names_gp= gpar(fontsize = 8),
        col= mycol2,
        name="legend",
        cluster_columns = F)

## Circos
circos.par(gap.after=c(250)) 
circos.heatmap(cir_B1_27,
               col=mycol,
               dend.side="inside",
               rownames.side="outside",
               track.height = 0.28, 
               rownames.col="black",
               bg.border="black",
               rownames.cex = 0.35)
circos.clear()

```


# B1 3025 circle heatmap NAbs
```{r Circle Heatmap of B1 NAbs}
library(tidyverse)
library(data.table)
library(dplyr)
library('ComplexHeatmap')
library('circlize')
library("RColorBrewer")
library(dendextend)
library(gridBase)

B1 <- fread("/DATA2/xuanjingyu/RQ3027 IIT/GMT/B1_0816.csv") %>% as.data.frame()

## Del infected row
B1_ninfec <- B1[!B1$ID %in% c("FB147","FB211","FB299","FB277", "FB361",
                             "FB117" ,"FB320" , "FB342", "FB336", "FB335" ,
                             "FB374",
                             "FB198" ,"FB206","FB189", "FB295",
                             "FB169","FB349","FB369"),]

## Rename
B1_ninfec$Patch <- ifelse(B1_ninfec$Patch == "B1","Day0",
                               ifelse(B1_ninfec$Patch == "B4", "Day14",
                                      B1_ninfec$Patch))
colnames(B1_ninfec)[colnames(B1_ninfec) == "SARS-CoV-2"] <- "VOC"

## Del V4 drop-out row
B1_ndrop <- B1_ninfec %>% filter(ID != "FB030" &
                                          ID != "FB188"&
                                          ID != "FB291") 

## Select XBB.1.16, BQ.1.1, BF.7,WT,BA.45 data
B1_variants <- subset(B1_ndrop, VOC %in% c("XBB.1.9.1",
               "XBB.1.16",
               "XBB.1.5",
               "BQ.1.1",
               "BF.7",
               "BA.4/5",
               "Beta",
               "Alpha",
               "Prototype"))
# B1_variants$Levels <- log(B1_variants$Levels)
B1_3025 <- subset(B1_variants, Group == "RQ3025")

B1_25w <- B1_3025 %>% pivot_wider(id_cols = c("ID"),
                                        names_from = "VOC",
                                        values_from = "Levels")
B1_25w_new <- column_to_rownames(B1_25w, var = "ID")


ma_B1_25 <- as.matrix(B1_25w_new)
voc_order <- c("XBB.1.9.1",
               "XBB.1.16",
               "XBB.1.5",
               "BQ.1.1",
               "BF.7",
               "BA.4/5",
               "Beta",
               "Alpha",
               "Prototype")

## Normalize data
cir_B1_25 <- t(scale(t(ma_B1_25)))

cir_B1_25 <- cir_B1_25[,order(factor(colnames(cir_B1_25), 
                                   levels = voc_order))]

range(cir_B1_25)

## Color
mycol=colorRamp2(c(-2, 0, 2),
                 c("#58A4FD","white","#5E1D9D"))
 
Heatmap(cir_B1_25,
        row_names_gp = gpar(fontsize = 7),
        column_names_gp= gpar(fontsize = 8),
        col= mycol2,
        name="legend",
        cluster_columns = F)
## Circos
circos.par(gap.after=c(250)) 
circos.heatmap(cir_B1_25,
               col=mycol,
               dend.side="inside",
               rownames.side="outside",
               track.height = 0.28, 
               rownames.col="black",
               bg.border="black",
               rownames.cex = 0.35)

circos.clear()

```

# B1 3013 circle heatmap NAbs
```{r Circle Heatmap of B1 NAbs}
library(tidyverse)
library(data.table)
library(dplyr)
library('ComplexHeatmap')
library('circlize')
library("RColorBrewer")
library(dendextend)
library(gridBase)

B1 <- fread("/DATA2/xuanjingyu/RQ3027 IIT/GMT/B1_0816.csv") %>% as.data.frame()

## Del infected row
B1_ninfec <- B1[!B1$ID %in% c("FB147","FB211","FB299","FB277", "FB361",
                             "FB117" ,"FB320" , "FB342", "FB336", "FB335" ,
                             "FB374",
                             "FB198" ,"FB206","FB189", "FB295",
                             "FB169","FB349","FB369"),]

## Rename
B1_ninfec$Patch <- ifelse(B1_ninfec$Patch == "B1","Day0",
                               ifelse(B1_ninfec$Patch == "B4", "Day14",
                                      B1_ninfec$Patch))
colnames(B1_ninfec)[colnames(B1_ninfec) == "SARS-CoV-2"] <- "VOC"

## Del V4 drop-out row
B1_ndrop <- B1_ninfec %>% filter(ID != "FB030" &
                                          ID != "FB188"&
                                          ID != "FB291") 

## Select XBB.1.16, BQ.1.1, BF.7,WT,BA.45 data
B1_variants <- subset(B1_ndrop, VOC %in% c("XBB.1.9.1",
               "XBB.1.16",
               "XBB.1.5",
               "BQ.1.1",
               "BF.7",
               "BA.4/5",
               "Beta",
               "Alpha",
               "Prototype"))
# B1_variants$Levels <- log(B1_variants$Levels)
B1_3013 <- subset(B1_variants, Group == "RQ3013")

B1_13w <- B1_3013 %>% pivot_wider(id_cols = c("ID"),
                                        names_from = "VOC",
                                        values_from = "Levels")
B1_13w_new <- column_to_rownames(B1_13w, var = "ID")


ma_B1_13 <- as.matrix(B1_13w_new)
voc_order <- c("XBB.1.9.1",
               "XBB.1.16",
               "XBB.1.5",
               "BQ.1.1",
               "BF.7",
               "BA.4/5",
               "Beta",
               "Alpha",
               "Prototype")

## Normalize data
cir_B1_13 <- t(scale(t(ma_B1_13)))

cir_B1_13 <- cir_B1_13[,order(factor(colnames(cir_B1_13), 
                                   levels = voc_order))]

range(cir_B1_13)

## Color
mycol=colorRamp2(c(-2, 0, 2),
                 c("#58A4FD","white","#5E1D9D"))
 
Heatmap(cir_B1_13,
        row_names_gp = gpar(fontsize = 7),
        column_names_gp= gpar(fontsize = 8),
        col= mycol2,
        name="legend",
        cluster_columns = F)
## Circos
circos.par(gap.after=c(250)) 
circos.heatmap(cir_B1_13,
               col=mycol,
               dend.side="inside",
               rownames.side="outside",
               track.height = 0.28, 
               rownames.col="black",
               bg.border="black",
               rownames.cex = 0.35)

leg_B1=Legend(title="value",
              col_fun=mycol,
              direction = c("vertical"))
grid.draw(leg_B1)

circos.clear()

```
