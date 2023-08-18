```{r overall}

library(pheatmap)
library(tidyverse)

voc_order <- c("Prototype",
                                                     "Alpha",
                                                     "Beta",
                                                     "BA.4/5",
                                                     "BF.7",
                                                     "BQ.1.1",
                                                     "XBB.1.5",
                                                     "XBB.1.16",
                                                     "XBB.1.9.1")
B1 <- fread("/DATA2/xuanjingyu/RQ3027 IIT/GMT/B1_0816.csv") %>% as.data.frame()
B4 <- fread("/DATA2/xuanjingyu/RQ3027 IIT/GMT/B4_0814.csv") %>% as.data.frame()

B1_variants <- subset(B1, `SARS-CoV-2` %in% c("XBB.1.16",
                                           "XBB.1.9.1",
                                           "XBB.1.5",
                                           "BQ.1.1",
                                           "BF.7",
                                           "BA.4/5",
                                           "Beta",
                                           "Alpha",
                                           "Prototype"))

B4_variants <- subset(B4, `SARS-CoV-2` %in% c("XBB.1.16",
                                           "XBB.1.9.1",
                                           "XBB.1.5",
                                           "BQ.1.1",
                                           "BF.7",
                                           "BA.4/5",
                                           "Beta",
                                           "Alpha",
                                           "Prototype"))

B1_w <- pivot_wider(B1_variants,id_cols = c("ID","Group",
                                                    "Patch"),
                                        names_from = "SARS-CoV-2",
                                        values_from = "Levels")
B4_w <- pivot_wider(B4_variants,id_cols = c("ID","Group",
                                                    "Patch"),
                                        names_from = "SARS-CoV-2",
                                        values_from = "Levels")

B1_w <- B1_w[,c(colnames(B1_w)[1:3],voc_order)]
B4_w <- B4_w[,c(colnames(B4_w)[1:3],voc_order)]

B1B4_w <- merge(B1_w[,-c(2:3)],B4_w[,-c(2:3)],by="ID")
B1B4_ninfec <- B1B4_w[!B1B4_w$ID %in% c("FB147","FB211","FB299","FB277", "FB361",
                             "FB117" ,"FB320" , "FB342", "FB336", "FB335" ,
                             "FB374",
                             "FB198" ,"FB206","FB189", "FB295",
                             "FB169","FB349","FB369"),]
B1B4_ndrop <- B1B4_ninfec %>% filter(ID != "FB030" &
                                          ID != "FB188"&
                                          ID != "FB291") %>% 
  column_to_rownames(.,"ID")


Day <- factor(rep(c("Day0","Day14"),each=9),levels=c("Day0","Day14"))
VOC <- factor(rep(voc_order,2),levels = voc_order)
Group <- B1_w[!B1_w$ID %in% c("FB147","FB211","FB299","FB277", "FB361",
                             "FB117" ,"FB320" , "FB342", "FB336", "FB335" ,
                             "FB374",
                             "FB198" ,"FB206","FB189", "FB295",
                             "FB169","FB349","FB369",
                             "FB030","FB188","FB291"),]
Group_d <- Group[,c(1,2)] %>% column_to_rownames(.,"ID")

annotation_col=data.frame(Day=Day,VOC=VOC)
rownames(annotation_col)=colnames(B1B4_ndrop)
annotation_row=Group_d
rownames(annotation_row)=rownames(Group_d)

colors=list(Day = c(Day0="#919191",Day14="#FF9200"),
            VOC = c("Prototype" = "#09357A",
                                "Alpha" = "#7725C6",
                                "Beta" = "#5A39E8",
                                "BA.4/5" = "#02915A",
                                "BF.7" = "#5181FF",
                                "BQ.1.1" = "#74D3D0",
                                "XBB.1.5" = "#F1EC5B",
                                "XBB.1.16" = "#FF903E",
                                "XBB.1.9.1" = "#FF322C"),
            Group = c(RQ3013 = "#09357A",
                      RQ3025 = "#02915A",
                      RQ3027 = "#740D91"))
p1 <- pheatmap(log2(B1B4_ndrop+1),
               scale = "row",
         show_colnames = F,
         show_rownames = F,
         cluster_rows = T,
         cluster_cols = F,
         fontsize=8,
         fontsize_row=13,
         border=F,
         annotation_colors=colors,
         annotation_row = annotation_row,
         annotation_names_row=T,
         annotation_names_col=F,
         annotation_col = annotation_col,
         breaks = seq(-1,1,length.out=100),
         gaps_col = 9,
         cutree_rows = 3,
         color = colorRampPalette(c("#58A4FD","white","#5E1D9D"))(100))
         # ,filename = '../output/DESeq/GSE137354/GSE137354_heatmap.pdf',height = 8,width =14)

```


```{r RQ3027+RQ3025+RQ3013}

library(pheatmap)
library(tidyverse)

voc_order <- c("Prototype",
                                                     "Alpha",
                                                     "Beta",
                                                     "BA.4/5",
                                                     "BF.7",
                                                     "BQ.1.1",
                                                     "XBB.1.5",
                                                     "XBB.1.16",
                                                     "XBB.1.9.1")
B1 <- fread("/DATA2/xuanjingyu/RQ3027 IIT/GMT/B1_0816.csv") %>% as.data.frame()
B4 <- fread("/DATA2/xuanjingyu/RQ3027 IIT/GMT/B4_0814.csv") %>% as.data.frame()

B1_variants <- subset(B1, `SARS-CoV-2` %in% c("XBB.1.16",
                                           "XBB.1.9.1",
                                           "XBB.1.5",
                                           "BQ.1.1",
                                           "BF.7",
                                           "BA.4/5",
                                           "Beta",
                                           "Alpha",
                                           "Prototype"))

B4_variants <- subset(B4, `SARS-CoV-2` %in% c("XBB.1.16",
                                           "XBB.1.9.1",
                                           "XBB.1.5",
                                           "BQ.1.1",
                                           "BF.7",
                                           "BA.4/5",
                                           "Beta",
                                           "Alpha",
                                           "Prototype"))

B1_w <- pivot_wider(B1_variants,id_cols = c("ID","Group",
                                                    "Patch"),
                                        names_from = "SARS-CoV-2",
                                        values_from = "Levels")
B4_w <- pivot_wider(B4_variants,id_cols = c("ID","Group",
                                                    "Patch"),
                                        names_from = "SARS-CoV-2",
                                        values_from = "Levels")

B1_w <- B1_w[,c(colnames(B1_w)[1:3],voc_order)]
B4_w <- B4_w[,c(colnames(B4_w)[1:3],voc_order)]

B1B4_w <- merge(B1_w,B4_w,by="ID")
B1B4_ninfec <- B1B4_w[!B1B4_w$ID %in% c("FB147","FB211","FB299","FB277", "FB361",
                             "FB117" ,"FB320" , "FB342", "FB336", "FB335" ,
                             "FB374",
                             "FB198" ,"FB206","FB189", "FB295",
                             "FB169","FB349","FB369"),]
B1B4_ndrop <- B1B4_ninfec %>% filter(ID != "FB030" &
                                          ID != "FB188"&
                                          ID != "FB291")

## Subset RQ3027
B1B4_27 <- subset(B1B4_ndrop, Group.x == "RQ3027")
B1B4_p27 <- B1B4_27[,-c(2,3,13,14)]
rownames(B1B4_p27) <- B1B4_p27$ID
B1B4_p27 <- B1B4_p27[,-1]

## Subset RQ3025
B1B4_25 <- subset(B1B4_ndrop, Group.x == "RQ3025")
B1B4_p25 <- B1B4_27[,-c(2,3,13,14)]
rownames(B1B4_p25) <- B1B4_p25$ID
B1B4_p25 <- B1B4_p25[,-1]

## Subset RQ3013
B1B4_13 <- subset(B1B4_ndrop, Group.x == "RQ3013")
B1B4_p13 <- B1B4_13[,-c(2,3,13,14)]
rownames(B1B4_p13) <- B1B4_p13$ID
B1B4_p13 <- B1B4_p13[,-1]

# Pheatmap
Day <- factor(rep(c("Day0","Day14"),each=9),levels=c("Day0","Day14"))
VOC <- factor(rep(voc_order,2),levels = voc_order)
# Group <- B1_w[!B1_w$ID %in% c("FB147","FB211","FB299","FB277", "FB361",
#                              "FB117" ,"FB320" , "FB342", "FB336", "FB335" ,
#                              "FB374",
#                              "FB198" ,"FB206","FB189", "FB295",
#                              "FB169","FB349","FB369",
#                              "FB030","FB188","FB291"),]
# Group_d <- Group[,c(1,2)] %>% column_to_rownames(.,"ID")

annotation_col=data.frame(Day=Day,VOC=VOC)
rownames(annotation_col)=colnames(B1B4_p27)
# annotation_row=Group_d
# rownames(annotation_row)=rownames(Group_d)

colors=list(Day = c(Day0="#919191",Day14="#FF9200"),
            VOC = c("Prototype" = "#09357A",
                                "Alpha" = "#7725C6",
                                "Beta" = "#5A39E8",
                                "BA.4/5" = "#02915A",
                                "BF.7" = "#5181FF",
                                "BQ.1.1" = "#74D3D0",
                                "XBB.1.5" = "#F1EC5B",
                                "XBB.1.16" = "#FF903E",
                                "XBB.1.9.1" = "#FF322C"))
            # Group = c(RQ3013 = "#09357A",
            #           RQ3025 = "#02915A",
            #           RQ3027 = "#740D91"))
pheat_27 <- pheatmap(log2(B1B4_p27+1),
               scale = "row",
         show_colnames = F,
         show_rownames = F,
         cluster_rows = T,
         cluster_cols = F,
         fontsize=8,
         fontsize_row=13,
         border=F,
         annotation_colors=colors,
         annotation_names_row=T,
         annotation_names_col=F,
         annotation_col = annotation_col,
         breaks = seq(-1,1,length.out=100),
         gaps_col = 9,
         color = colorRampPalette(c("#58A4FD","white","#5E1D9D"))(100))
         # ,filename = '../output/DESeq/GSE137354/GSE137354_heatmap.pdf',height = 8,width =14)

pheat_25 <- pheatmap(log2(B1B4_p25+1),
               scale = "row",
         show_colnames = F,
         show_rownames = F,
         cluster_rows = T,
         cluster_cols = F,
         fontsize=8,
         fontsize_row=13,
         border=F,
         annotation_colors=colors,
         annotation_names_row=T,
         annotation_names_col=F,
         annotation_col = annotation_col,
         breaks = seq(-1,1,length.out=100),
         gaps_col = 9,
         color = colorRampPalette(c("#58A4FD","white","#5E1D9D"))(100))
         # ,filename = '../output/DESeq/GSE137354/GSE137354_heatmap.pdf',height = 8,width =14)

pheat_13 <- pheatmap(log2(B1B4_p13+1),
               scale = "row",
         show_colnames = F,
         show_rownames = F,
         cluster_rows = T,
         cluster_cols = F,
         fontsize=8,
         fontsize_row=13,
         border=F,
         annotation_colors=colors,
         annotation_names_row=T,
         annotation_names_col=F,
         annotation_col = annotation_col,
         breaks = seq(-1,1,length.out=100),
         gaps_col = 9,
         color = colorRampPalette(c("#58A4FD","white","#5E1D9D"))(100))
         # ,filename = '../output/DESeq/GSE137354/GSE137354_heatmap.pdf',height = 8,width =14)
```