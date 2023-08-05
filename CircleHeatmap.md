# B1 circle heatmap NAbs
```{r Circle Heatmap of B1 NAbs}
library(tidyverse)
library(data.table)
library(dplyr)
library('ComplexHeatmap')
library('circlize')
library("RColorBrewer")
library(dendextend)
library(gridBase)

B1 <- fread("/DATA2/xuanjingyu/RQ3027 IIT/GMT/B1_0722.csv") %>% as.data.frame()

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
B1_variants <- subset(B1_ndrop, VOC %in% c("XBB.1.16","BQ.1.1","BF.7",
                                               "Prototype","BA.4/5"))
# B1_variants$Levels <- log(B1_variants$Levels)


B1_wider <- B1_variants %>% pivot_wider(id_cols = c("ID",
                                                    "Group"),
                                        names_from = "VOC",
                                        values_from = "Levels")
B1_wider_new <- column_to_rownames(B1_wider, var = "ID")
B1_w <- B1_wider_new[,-c(1)]

g_rowB1 <- B1_wider[,3]
g_rowB1 <- as.matrix(g_rowB1)

cir_B1 <- as.matrix(B1_w)
voc_order <- c("XBB.1.16",
               "BQ.1.1",
               "BF.7",
               "BA.4/5",
               "Prototype")

## Normalize data
cir_B1 <- t(scale(t(B1_w)))
# cir_B1 <- cir_B1[,order(factor(colnames(cir_B1), 
#                                levels = voc_order))]

range(cir_B1)

## Color
mycol=colorRamp2(c(-2, 0, 2),
                 c("blue","white","red"))
 


Heatmap(cir_B1,row_names_gp = gpar(fontsize = 3),
 
        column_names_gp= gpar(fontsize = 8),
 
        col= mycol,
 
        name="legend")

## Circos
circos.par(gap.after=c(25)) 
circos.heatmap(cir_B1,
               col=mycol,
               dend.side="inside",
               rownames.side="outside",#rownames.side：控制矩阵行名的方向,与dend.side相同；但注意二者不能在同一侧，必须一内一外
               track.height = 0.28, #轨道的高度，数值越大圆环越粗
               rownames.col="black",
               bg.border="black", #背景边缘颜色
               # split = g_rowB1,#用行注释分裂热图
               # show.sector.labels = T,
               # rownames.cex=0.4,#字体大小
               # rownames.font=1,#字体粗细
               # cluster=TRUE,#cluster=TRUE为对行聚类，cluster=FALSE则不显示聚类
               # # dend.track.height=0.18,#调整行聚类树的高度
               # dend.callback=function(dend,m,si) {#dend.callback：用于聚类树的回调，当需要对聚类树进行重新排序，或者添加颜色时使用包含的三个参数：dend：当前扇区的树状图；m：当前扇区对应的子矩阵；si：当前扇区的名称
               #   color_branches(dend,k=10,col=1:10)#color_branches():修改聚类树颜色#聚类树颜色改为1，即单色/黑色
               # }
)


leg_B1=Legend(title="value",col_fun=mycol,
 
          direction = c("vertical"),
 
          #title_position= c('topcenter')，
 
)
 
grid.draw(leg_B1)
 
#draw(lg, x = unit(0.9,"npc"), y = unit(0.5,"npc"), just = c("right","center"))#画在右边
 
#添加列名：
 
circos.track(track.index=get.current.track.index(),panel.fun=function(x,y){
  if(CELL_META$sector.numeric.index==1){
    cn=colnames(cir_B1)
    n=length(cn)
    circos.text(rep(CELL_META$cell.xlim[2],n) + 
                  convert_x(0.5,"mm"),#x坐标
                1:n+5,#调整y坐标,行距+距离中心距(1:n)*1.2+5,
                cn,
                cex=0.5,
                adj=c(0,1),
                facing="inside")}
  },bg.border=NA)

circos.clear()
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

B1 <- fread("/DATA2/xuanjingyu/RQ3027 IIT/GMT/B1_0722.csv") %>% as.data.frame()

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
B1_variants <- subset(B1_ndrop, VOC %in% c("XBB.1.16","BQ.1.1","BF.7",
                                               "Prototype","BA.4/5"))
# B1_variants$Levels <- log(B1_variants$Levels)
B1_3027 <- subset(B1_variants, Group == "RQ3027")

B1_27w <- B1_3027 %>% pivot_wider(id_cols = c("ID"),
                                        names_from = "VOC",
                                        values_from = "Levels")
B1_27w_new <- column_to_rownames(B1_27w, var = "ID")


ma_B1_27 <- as.matrix(B1_27w_new)
voc_order <- c("Prototype",
               "BA.4/5",
               "BF.7",
               "BQ.1.1",
               "XBB.1.16")

## Normalize data
cir_B1_27 <- t(scale(t(ma_B1_27)))

cir_B1_27 <- cir_B1_27[,order(factor(colnames(cir_B1_27), 
                                   levels = voc_order))]

range(cir_B1_27)

## Color
mycol=colorRamp2(c(-2, 0, 2),
                 c("#58A4FD","white","#5E1D9D"))
 


Heatmap(cir_B1_27,row_names_gp = gpar(fontsize = 3),
 
        column_names_gp= gpar(fontsize = 8),
 
        col= mycol,
 
        name="legend")

## Circos
circos.par(gap.after=c(250)) 
circos.heatmap(cir_B1_27,
               col=mycol,
               dend.side="inside",
               rownames.side="outside",
               track.height = 0.28, 
               rownames.col="black",
               bg.border="black", 

)

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

B1 <- fread("/DATA2/xuanjingyu/RQ3027 IIT/GMT/B1_0722.csv") %>% as.data.frame()

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
B1_variants <- subset(B1_ndrop, VOC %in% c("XBB.1.16","BQ.1.1","BF.7",
                                               "Prototype","BA.4/5"))
# B1_variants$Levels <- log(B1_variants$Levels)
B1_3025 <- subset(B1_variants, Group == "RQ3025")

B1_25w <- B1_3025 %>% pivot_wider(id_cols = c("ID"),
                                        names_from = "VOC",
                                        values_from = "Levels")
B1_25w_new <- column_to_rownames(B1_25w, var = "ID")


ma_B1_25 <- as.matrix(B1_25w_new)
voc_order <- c("Prototype",
               "BA.4/5",
               "BF.7",
               "BQ.1.1",
               "XBB.1.16")

## Normalize data
cir_B1_25 <- t(scale(t(ma_B1_25)))

cir_B1_25 <- cir_B1_25[,order(factor(colnames(cir_B1_25), 
                                   levels = voc_order))]

range(cir_B1_25)

## Color
mycol=colorRamp2(c(-2, 0, 2),
                 c("#58A4FD","white","#5E1D9D"))
 

## Circos
circos.par(gap.after=c(250)) 
circos.heatmap(cir_B1_25,
               col=mycol,
               dend.side="inside",
               rownames.side="outside",
               track.height = 0.28, 
               rownames.col="black",
               bg.border="black", 

)

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

B1 <- fread("/DATA2/xuanjingyu/RQ3027 IIT/GMT/B1_0722.csv") %>% as.data.frame()

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
B1_variants <- subset(B1_ndrop, VOC %in% c("XBB.1.16","BQ.1.1","BF.7",
                                               "Prototype","BA.4/5"))
# B1_variants$Levels <- log(B1_variants$Levels)
B1_3013 <- subset(B1_variants, Group == "RQ3013")

B1_13w <- B1_3013 %>% pivot_wider(id_cols = c("ID"),
                                        names_from = "VOC",
                                        values_from = "Levels")
B1_13w_new <- column_to_rownames(B1_13w, var = "ID")


ma_B1_13 <- as.matrix(B1_13w_new)
voc_order <- c("Prototype",
               "BA.4/5",
               "BF.7",
               "BQ.1.1",
               "XBB.1.16")

## Normalize data
cir_B1_13 <- t(scale(t(ma_B1_13)))

cir_B1_13 <- cir_B1_13[,order(factor(colnames(cir_B1_13), 
                                   levels = voc_order))]

range(cir_B1_13)

## Color
mycol=colorRamp2(c(-2, 0, 2),
                 c("#58A4FD","white","#5E1D9D"))
 

## Circos
circos.par(gap.after=c(250)) 
circos.heatmap(cir_B1_13,
               col=mycol,
               dend.side="inside",
               rownames.side="outside",
               track.height = 0.28, 
               rownames.col="black",
               bg.border="black")

leg_B1=Legend(title="value",
              col_fun=mycol,
              direction = c("vertical"))
grid.draw(leg_B1)

circos.clear()

```