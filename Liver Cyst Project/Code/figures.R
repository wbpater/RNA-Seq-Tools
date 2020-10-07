### Fig A. ####################################################################################################
# Total heatmap on batch corrected, vst normalised counts
salmon_load(make_tx2gene = T, annotation = "/scratch/wpater/genomes_and_annotations/gencode.v29.annotation.gtf")

salmon_load(data_table = "/scratch/wpater/fig_a_table.txt", salmon_dir = "/scratch/wpater/analysis/liver_organoids/salmon_output/",
            DESeq_design = "Mutation", control = "None", returnMode = "vst", batch_correction = T, batch_group = "batch", k=2) # Load salmon quant file, batch correct counts, vst normalise
vst=ensembl_id2symbol(vst) # set rownames to gene symbols
final_genes = c("PRKCSH", "SEC63","GANAB","SEC61B", "PKD1", "PKD2","ALG8","HNF1B", "LRP5") # list with final genes
column_order = c("MDL_1-1", "MDL_1-2","MDL_3.3", "MDL_3.4","MDL_3.5","MDL_3.6","EDG_n2","EDG_n3","EDG_n4","EDG_n10",
                 "MDL_1-5", "EDG_9p6","EDG_16p6","MDL_2-6","MDL_3.8","MDL_2-7","EDG_15p6","MDL_2-8","MDL_3.7",
                 "MDL_1-6", "EDG_8p6", "MDL_1-7", "EDG_10p4","MDL_1-8","EDG_12p5")
vst= vst[final_genes,column_order] # select rows with final_genes as rownames
library(pheatmap)
df = as.data.frame(colData(vst)[,c("Selection","Mutation", "LOH", "Library_prep","Duplicates", "Sex","heatmap_lab")]) # make dataframe for heatmap annotation
colnames(df)[colnames(df) == "Library_prep"] = "Library prep"
ann_colors = list(Mutation=c(None="#999999",PRKCSH="#00ff00",Unknown="#ff0000",SEC63="#ff00ff",PKD1="#0000ff"),
                  "Library prep" = c("PolyA selection"="#009292","Ribo-Zero run 1"="#b66dff","Ribo-Zero run 2"="#db6d00"),
                  Duplicates=c("Unique"= "#dbdbdb","Duplicate 1"= "#ff0000","Duplicate 2"= "#ffd500", "Duplicate 3"= "#00fff3", "Duplicate 4"="#0051ff" ),
                  Sex=c("Male"= "steelblue", "Female"= "plum"),
                  LOH=c("Yes"="royalblue","No"="grey100"),
                  Selection=c("Not selected"="grey", "PLD"="grey100", "Control"="grey100")) # Colours for the heatmap annotation
breakslist = c(-2.5,-2,1.9,-1.8,-1.6,-1.4,-1.2,-1.0,-0.8,-0.5,0,0.5,0.75,1,1.25,1.5,2,2.5)
pheatmap(assay(vst), cluster_rows=F, show_rownames=T, cluster_cols=F, breaks = breakslist, show_colnames = F,
         annotation_col=df[,1:6], annotation_colors = ann_colors, scale = "row", color = colorRampPalette(c("mediumblue", "white","red"))(length(breakslist)+2),
         main = paste0("Cyst related genes in liver organoids"), na_col = "white", labels_col = as.character(df$heatmap_lab), display_numbers = F,
         angle_col = 315)
# Save as 8.5 *6.5 inch pdf

### Fig B. ####################################################################################################
## PCA of selected controls and PLD samples of batch corrected, vst normalised counts
# using an annotation file where genes mapping to chrY are removed
# bash code: grep -v "chrY" my.gtf > no_y.gtf
salmon_load(make_tx2gene = T, annotation = "/scratch/wpater/genomes_and_annotations/noY_gencode.v29.annotation.gtf")

salmon_load(data_table = "/scratch/wpater/fig_b_table.txt", salmon_dir = "/scratch/wpater/analysis/liver_organoids/salmon_output/",
            DESeq_design = "Selection", control = "Control", returnMode = "vst", batch_correction = T, batch_group = "batch", k=2)
vst=ensembl_id2symbol(vst)
final_samples = c("MDL_1-5","MDL_2-6","MDL_2-7","MDL_2-8","MDL_3.3","MDL_3.4","MDL_3.5","MDL_3.6","EDG_16p6","EDG_n2","EDG_n3","EDG_n4")
vst= vst[,final_samples]
cbPalette = c(Control="#0000ff",PLD="#ff0000",na="#999999")
PCAplotR(vst, shape = "Mutation", colour = "Selection" ,PCX=1, PCY=2, ntop = 500, palette = cbPalette)+
  ggtitle("Selected controls and PLD samples", subtitle = "Top 500 highest variance genes of batch corrected, vst normalised counts")
# Save as 6.5 *5 inch pdf

## PC that shows the effect of sex
cbPalette = c(Male="#0000ff",Female="#ff0000",na="#999999")
for (i in 1:11) {print(PCAplotR(vst, shape = "Mutation", colour = "Sex" ,PCX=i, PCY=i+1, ntop = 500, palette = cbPalette)+
  ggtitle("The effect of sex in controls and PLD samples ", subtitle = "Top 500 highest variance genes of batch corrected, vst normalised counts"))}

### Table A. ####################################################################################################
## Standard DEG analysis
salmon_load(data_table = "/scratch/wpater/fig_b_table.txt", salmon_dir = "/scratch/wpater/analysis/liver_organoids/salmon_output/",
            DESeq_design = "Selection", control = "Control", returnMode = "dds", batch_correction = T, batch_group = "batch", k=2)
sample1 = "PLD"
sample2 = "Control"
res = results(dds, contrast=c("Selection", sample1, sample2))                         # compute DEG between sample 1 and sample 2
res=ensembl_id2symbol(res)
res[is.na(res$padj),] = 1
res = res[res$padj < 0.05,]
length(res$log2FoldChange[which(res$log2FoldChange > 0)]) #number upregulated
length(res$log2FoldChange[which(res$log2FoldChange < 0)]) # number downregulated
write.csv(res, file = "/scratch/wpater/PLD_v_Control.csv")
# Import pathway analysis table from Enrichr

### Fig C. Panel C ####################################################################################################
# PRKCSH vs control comparison
salmon_load(data_table = "/scratch/wpater/fig_b_table.txt", salmon_dir = "/scratch/wpater/analysis/liver_organoids/salmon_output/",
            DESeq_design = "Targeted", control = "Control", returnMode = c("vst","dds"), batch_correction = T, batch_group = "batch", k=2)

sample1 = "PRKCSH"
sample2 = "Control"
res = results(dds, contrast=c("Targeted", sample1, sample2))                         # compute DEG between sample 1 and sample 2
res=ensembl_id2symbol(res)
res[is.na(res$padj),] = 1
res = res[res$padj < 0.05,]
length(res$log2FoldChange[which(res$log2FoldChange > 0)]) #number upregulated
length(res$log2FoldChange[which(res$log2FoldChange < 0)]) # number downregulated
write.csv(res, file = "/scratch/wpater/PRKCSH_v_Control.csv")

vst=ensembl_id2symbol(vst)
final_samples = c("MDL_3.3","MDL_3.4","MDL_3.5","MDL_3.6","MDL_1-5","EDG_16p6", "EDG_n2", "EDG_n3", "EDG_n4")
final_genes = c("PDIA3","HSPA5","PDIA6","PDIA4","OS9","CALR")
vst= vst[final_genes,final_samples] # select rows with final_genes as rownames
library(pheatmap)
df = as.data.frame(colData(vst)[,c("Mutation","Sex", "Library_prep","name")]) # make dataframe for heatmap annotation
colnames(df)[colnames(df) == "Library_prep"] = "Library prep"
ann_colors = list(Mutation=c(None="#999999","PRKCSH LOH"="#00ff00"),
                  "Library prep" = c("PolyA selection"="#009292","Ribo-Zero run 1"="#b66dff","Ribo-Zero run 2"="#db6d00"),
                  Sex=c("Male"= "steelblue", "Female"= "plum")) # Colours for the heatmap annotation
breakslist = seq(-2.1,2.0, by = 0.10) # scaling of the heatmap colours
set.seed(10)
pheatmap(assay(vst), cluster_rows=T, show_rownames=T, cluster_cols=T, breaks = breakslist,
         annotation_col=df[,1:2], annotation_colors = ann_colors, scale = "row", color = colorRampPalette(c("blue", "white","red"))(length(breakslist)),
         main = paste0("ER protein processing genes in PLD organoids"), na_col = "white", labels_col = as.character(df$name), 
         display_numbers = F, show_colnames = F)
# Save as 8.5 *6.5 inch pdf

### Fig C. Panel A ####################################################################################################
salmon_load(data_table = "/scratch/wpater/fig_b_table.txt", salmon_dir = "/scratch/wpater/analysis/liver_organoids/salmon_output/",
            DESeq_design = "Targeted", control = "Control", returnMode = c("vst","dds"), batch_correction = T, batch_group = "batch", k=2)
sample1 = "PRKCSH"
sample2 = "Control"
res = results(dds, contrast=c("Targeted", sample1, sample2))                         # compute DEG between sample 1 and sample 2
res=ensembl_id2symbol(res)
res[is.na(res$padj),] = 1
res = res[res$padj < 0.05,]
row.names(res)

vst=ensembl_id2symbol(vst)
final_samples = c("MDL_3.3","MDL_3.4","MDL_3.5","MDL_3.6","MDL_1-5","EDG_16p6", "EDG_n2", "EDG_n3", "EDG_n4")
final_genes = row.names(res)
# Note: removing non coding genes 
final_genes = row.names(res)[!grepl(pattern = "ENSG", x = row.names(res))]
vst= vst[final_genes,final_samples] # select rows with final_genes as rownames
library(pheatmap)
df = as.data.frame(colData(vst)[,c("Mutation","Sex", "Library_prep","name")]) # make dataframe for heatmap annotation
colnames(df)[colnames(df) == "Library_prep"] = "Library prep"
ann_colors = list(Mutation=c(None="#999999","PRKCSH LOH"="#00ff00"),
                  "Library prep" = c("PolyA selection"="#009292","Ribo-Zero run 1"="#b66dff","Ribo-Zero run 2"="#db6d00"),
                  Sex=c("Male"= "steelblue", "Female"= "plum")) # Colours for the heatmap annotation
breakslist = seq(-2.1,2.0, by = 0.10) # scaling of the heatmap colours
set.seed(10)
pheatmap(assay(vst), cluster_rows=T, show_rownames=T, cluster_cols=T, breaks = breakslist,
         annotation_col=df[,1:2], annotation_colors = ann_colors, scale = "row", color = colorRampPalette(c("blue", "white","red"))(length(breakslist)),
         main = paste0("Significantly different protein coding genes between control and PRKCSH LOH"), na_col = "white", labels_col = as.character(df$name), 
         display_numbers = F, show_colnames = F)
# Save as 8.5 *6.5 inch pdf

### Fig C. Panel B ####################################################################################################
salmon_load(data_table = "/scratch/wpater/fig_b_table.txt", salmon_dir = "/scratch/wpater/analysis/liver_organoids/salmon_output/",
            DESeq_design = "Targeted", control = "Control", returnMode = c("vst","dds"), batch_correction = T, batch_group = "batch", k=2)
sample1 = "PRKCSH"
sample2 = "Control"
res = results(dds, contrast=c("Targeted", sample1, sample2))                         # compute DEG between sample 1 and sample 2
res=ensembl_id2symbol(res)
res[is.na(res$padj),] = 1
res = res[res$padj < 0.05,]
row.names(res)

vst=ensembl_id2symbol(vst)
final_samples = c("MDL_3.3","MDL_3.4","MDL_3.5","MDL_3.6","MDL_1-5","EDG_16p6", "EDG_n2", "EDG_n3", "EDG_n4")
final_genes = row.names(res)[grepl(pattern = "ENSG", x = row.names(res))]

vst= vst[final_genes,final_samples] # select rows with final_genes as rownames
library(pheatmap)
df = as.data.frame(colData(vst)[,c("Mutation","Sex", "Library_prep","name")]) # make dataframe for heatmap annotation
colnames(df)[colnames(df) == "Library_prep"] = "Library prep"
ann_colors = list(Mutation=c(None="#999999","PRKCSH LOH"="#00ff00"),
                  "Library prep" = c("PolyA selection"="#009292","Ribo-Zero run 1"="#b66dff","Ribo-Zero run 2"="#db6d00"),
                  Sex=c("Male"= "steelblue", "Female"= "plum")) # Colours for the heatmap annotation
breakslist = seq(-2.1,2.0, by = 0.10) # scaling of the heatmap colours
set.seed(10)
pheatmap(assay(vst), cluster_rows=T, show_rownames=T, cluster_cols=T, breaks = breakslist,
         annotation_col=df[,1:2], annotation_colors = ann_colors, scale = "row", color = colorRampPalette(c("blue", "white","red"))(length(breakslist)),
         main = paste0("Significantly different non-coding genes between control and PRKCSH LOH"), na_col = "white", labels_col = as.character(df$name), 
         display_numbers = F, show_colnames = F)
# Save as 8.5 *6.5 inch pdf




## Wide view heatmap
salmon_load(data_table = "/scratch/wpater/fig_b_table.txt", salmon_dir = "/scratch/wpater/analysis/liver_organoids/salmon_output/",
            DESeq_design = "Targeted", control = "Control", returnMode = c("vst","dds"), batch_correction = T, batch_group = "batch", k=2)
vst=ensembl_id2symbol(vst)
final_genes = c("PDIA3","HSPA5","PDIA6","PDIA4","OS9","CALR")
final_samples = c( "MDL_1-5" , "MDL_1-6" , "MDL_1-7",  "MDL_1-8" , "MDL_2-6" , "MDL_2-7" , "MDL_2-8" , "MDL_3.3" ,"MDL_3.4" , "MDL_3.5" , "MDL_3.6" , "MDL_3.7" , "MDL_3.8",
                   "EDG_n2" ,  "EDG_n3" , "EDG_n4" ,  "EDG_8p6" , "EDG_9p6" , "EDG_10p4" ,"EDG_12p5", "EDG_15p6", "EDG_16p6")
vst= vst[final_genes,final_samples] # select rows with final_genes as rownames
library(pheatmap)
df = as.data.frame(colData(vst)[,c("Mutation","Sex", "Library_prep","name")]) # make dataframe for heatmap annotation
colnames(df)[colnames(df) == "Library_prep"] = "Library prep"
ann_colors = list(Mutation=c(None="#999999",PRKCSH="#006400","PRKCSH LOH"="#00ff00",Unknown="#ff0000","SEC63 LOH"="#ff00ff", "PKD1"="#0000ff"),
                  "Library prep" = c("PolyA selection"="#009292","Ribo-Zero run 1"="#b66dff","Ribo-Zero run 2"="#db6d00"),
                  Sex=c("Male"= "steelblue", "Female"= "plum")) # Colours for the heatmap annotation
breakslist = seq(-2.1,2.0, by = 0.10) # scaling of the heatmap colours
set.seed(10)
pheatmap(assay(vst), cluster_rows=T, show_rownames=T, cluster_cols=T, breaks = breakslist,
         annotation_col=df[,1:2], annotation_colors = ann_colors, scale = "row", color = colorRampPalette(c("blue", "white","red"))(length(breakslist)),
         main = paste0("ER protein processing genes in PLD organoids"), na_col = "white", labels_col = as.character(df$name), 
         display_numbers = F, show_colnames = T)

#### Fig D-E Pathway overview ####################################################################################################
library(org.Hs.eg.db)
kegg = org.Hs.egPATH2EG
mapped = mappedkeys(kegg)
kegg2 = as.list(kegg[mapped])
genes = as.character(kegg2$`04141`)
library("AnnotationDbi")
genes = as.data.frame(mapIds(org.Hs.eg.db, keys=genes, column="SYMBOL",                 # Make an extra column with the gene sybmol matching the ENSMBL id
       keytype="ENTREZID", multiVals="first")  )[,1]                           # 
salmon_load(data_table = "/scratch/wpater/fig_b_table.txt", salmon_dir = "/scratch/wpater/analysis/liver_organoids/salmon_output/",
            DESeq_design = "Targeted", control = "Control", returnMode = "dds", batch_correction = T, batch_group = "batch", k=2)
sample1 = "PRKCSH"
sample2 = "Control"
res = results(dds, contrast=c("Targeted", sample1, sample2))                         # compute DEG between sample 1 and sample 2
res = ensembl_id2symbol(res)
res[is.na(res$padj),] = 1
res = res[as.character(genes[as.character(genes) %in% rownames(res)]),] # filter out genes not in vst dataset
write.csv(res, file = "/scratch/wpater/ER_pathway_genes.csv")

# the ones not in pathway
undetected = as.character(genes[!(genes %in% rownames(res))])


### Supplemental B ############################################################################################
## Batch correction rationale and "validation"

# Not corrected version PC 1+2
salmon_load(data_table = "/scratch/wpater/fig_sa_table.txt", salmon_dir = "/scratch/wpater/analysis/liver_organoids/salmon_output/",
            DESeq_design = "Mutation", control = "None", returnMode = "vst")
cbPalette = c(Control="#0000ff",PLD="#ff0000","Not selected" ="#999999")
PCAplotR(vst, shape = "Library_prep", colour = "Selection" ,PCX=1, PCY=2, ntop = 500, palette = cbPalette, label = "Duplicates")+
  ggtitle("All controls and PLD samples", subtitle = "Top 500 highest variance genes of vst normalised counts")+
  scale_shape_discrete(name = "Library prep")
# Save as 6.5 *5 inch pdf

# Batch corrected PCA with all samples
salmon_load(data_table = "/scratch/wpater/fig_sa_table.txt", salmon_dir = "/scratch/wpater/analysis/liver_organoids/salmon_output/",
            DESeq_design = "Mutation", control = "None", returnMode = "vst", batch_correction = T, batch_group = "batch", k=2)
cbPalette = c(Control="#0000ff",PLD="#ff0000","Not selected" ="#999999")
PCAplotR(vst, shape = "Library_prep", colour = "Selection" ,PCX=1, PCY=2, ntop = 500, palette = cbPalette, label = "Duplicates")+
  ggtitle("All controls and PLD samples", subtitle = "Top 500 highest variance genes of batch corrected, vst normalised counts")+
  scale_shape_discrete(name = "Library prep")
# Save as 6.5 *5 inch pdf

### Supplemental C ####################################################
## PCA Controls vs. PRKCSH LOH
salmon_load(data_table = "/scratch/wpater/fig_sb_table.txt", salmon_dir = "/scratch/wpater/analysis/liver_organoids/salmon_output/",
            DESeq_design = "Selection", control = "Control", returnMode = "vst", batch_correction = T, batch_group = "batch", k=2)
vst=ensembl_id2symbol(vst)
final_samples = c("MDL_1-5","MDL_3.3","MDL_3.4","MDL_3.5","MDL_3.6","EDG_16p6","EDG_n2","EDG_n3","EDG_n4")
vst= vst[,final_samples]
cbPalette = c(Control="#0000ff","PRKCSH LOH"="#00ff00",na="#999999")
PCAplotR(vst, colour = "Selection" ,PCX=1, PCY=2, ntop = 500, palette = cbPalette)+
  ggtitle("Selected controls and PRKCSH LOH samples ", subtitle = "Top 500 highest variance genes of batch corrected, vst normalised counts")
# Save as 6.5 *5 inch pdf


#### Supplemental D SEC63 LOH DEG ####################################################################################################
salmon_load(data_table = "/scratch/wpater/fig_g_table.txt", salmon_dir = "/scratch/wpater/analysis/liver_organoids/salmon_output/",
            DESeq_design = "Targeted", control = "Control", returnMode = c("vst","dds"), batch_correction = T, batch_group = "batch", k=2)
sample1 = "SEC63 LOH"
sample2 = "Control"
res = results(dds, contrast=c("Targeted", sample1, sample2))                         # compute DEG between sample 1 and sample 2
res=ensembl_id2symbol(res)
res[is.na(res$padj),] = 1
res = res[res$padj < 0.05,]
length(res$log2FoldChange[which(res$log2FoldChange > 0)]) #number upregulated
length(res$log2FoldChange[which(res$log2FoldChange < 0)]) # number downregulated
write.csv(res, file = "/scratch/wpater/SEC63_v_Control.csv")

vst=ensembl_id2symbol(vst)
final_samples = c("MDL_3.3","MDL_3.4","MDL_3.5","MDL_3.6", "EDG_n2", "EDG_n3", "EDG_n4","MDL_2-7")
vst= vst[,final_samples]
cbPalette = c(Control="#0000ff","SEC63 LOH"="#ff00ff",na="#999999")
PCAplotR(vst,  colour = "Targeted" ,PCX=1, PCY=2, ntop = 500, palette = cbPalette)+
  ggtitle("Selected controls and SEC63 LOH sample ", subtitle = "Top 500 highest variance genes of batch corrected, vst normalised counts")
# Save as 6.5 *5 inch pdf

# final_genes = c("HIST1H2AM","HIST1H2BO","HIST1H3J","HIST1H2AH","HIST1H2BI","HIST1H2AJ","HIST1H3A","HIST2H3A","HIST1H3F","HIST1H2BF","HIST1H2AG","HIST1H3B","HIST1H3C","HIST2H3C","HLA-DQB1") # lupus
# final_genes=c("HIST1H2AM","HIST1H2BO","HIST1H3J","HDAC1","HIST1H2AH","HIST1H2BI","HIST1H2AJ","HIST1H3A","HIST2H3A","ADORA2B","HIST1H3F","HIST1H2BF","HIST1H2AG","HIST1H3B","HIST1H3C","HIST2H3C") # alcoholism
# vst= vst[final_genes,final_samples] # select rows with final_genes as rownames
# library(pheatmap)
# df = as.data.frame(colData(vst)[,c("Mutation","Sex", "Library_prep","name")]) # make dataframe for heatmap annotation
# colnames(df)[colnames(df) == "Library_prep"] = "Library prep"
# ann_colors = list(Mutation=c(None="#999999","SEC63 LOH"="#ff00ff"),
#                   "Library prep" = c("PolyA selection"="#009292","Ribo-Zero run 1"="#b66dff","Ribo-Zero run 2"="#db6d00"),
#                   Sex=c("Male"= "steelblue", "Female"= "plum")) # Colours for the heatmap annotation
# breakslist = seq(-2.1,2.0, by = 0.10) # scaling of the heatmap colours
# set.seed(10)
# pheatmap(assay(vst), cluster_rows=T, show_rownames=T, cluster_cols=T, breaks = breakslist,
#          annotation_col=df[,1:2], annotation_colors = ann_colors, scale = "row", color = colorRampPalette(c("blue", "white","red"))(length(breakslist)),
#          main = paste0("SEC63 enriched genes in PLD organoids"), na_col = "white", labels_col = as.character(df$name), display_numbers = F, show_colnames = F)
# ## Wide view heatmap
# salmon_load(data_table = "/scratch/wpater/fig_g_table.txt", salmon_dir = "/scratch/wpater/analysis/liver_organoids/salmon_output/",
#             DESeq_design = "Targeted", control = "Control", returnMode = c("vst","dds"), batch_correction = T, batch_group = "batch", k=2)
# vst=ensembl_id2symbol(vst)
# final_samples = c( "MDL_1-5" , "MDL_1-6" , "MDL_1-7",  "MDL_1-8" , "MDL_2-6" , "MDL_2-7" , "MDL_2-8" , "MDL_3.3" ,"MDL_3.4" , "MDL_3.5" , "MDL_3.6" , "MDL_3.7" , "MDL_3.8",
#                    "EDG_n2" ,  "EDG_n3" , "EDG_n4" ,  "EDG_8p6" , "EDG_9p6" , "EDG_10p4" ,"EDG_12p5", "EDG_15p6", "EDG_16p6")
# final_genes = c("HIST1H2AM","HIST1H2BO","HIST1H3J","HIST1H2AH","HIST1H2BI","HIST1H2AJ","HIST1H3A","HIST2H3A","HIST1H3F","HIST1H2BF","HIST1H2AG","HIST1H3B","HIST1H3C","HIST2H3C","HLA-DQB1") # lupus
# final_genes=c("HIST1H2AM","HIST1H2BO","HIST1H3J","HDAC1","HIST1H2AH","HIST1H2BI","HIST1H2AJ","HIST1H3A","HIST2H3A","ADORA2B","HIST1H3F","HIST1H2BF","HIST1H2AG","HIST1H3B","HIST1H3C","HIST2H3C") # alcoholism
# vst= vst[final_genes,final_samples] # select rows with final_genes as rownames
# library(pheatmap)
# df = as.data.frame(colData(vst)[,c("Mutation","Sex", "Library_prep","name")]) # make dataframe for heatmap annotation
# colnames(df)[colnames(df) == "Library_prep"] = "Library prep"
# ann_colors = list(Mutation=c(None="#999999",PRKCSH="#006400","PRKCSH LOH"="#00ff00",Unknown="#ff0000","SEC63 LOH"="#ff00ff", "PKD1"="#0000ff"),
#                   "Library prep" = c("PolyA selection"="#009292","Ribo-Zero run 1"="#b66dff","Ribo-Zero run 2"="#db6d00"),
#                   Sex=c("Male"= "steelblue", "Female"= "plum")) # Colours for the heatmap annotation
# breakslist = seq(-2.1,2.0, by = 0.10) # scaling of the heatmap colours
# set.seed(10)
# pheatmap(assay(vst), cluster_rows=T, show_rownames=T, cluster_cols=T, breaks = breakslist,
#          annotation_col=df[,1:2], annotation_colors = ann_colors, scale = "row", color = colorRampPalette(c("blue", "white","red"))(length(breakslist)),
#          main = paste0("SEC63 enriched genes in PLD organoids"), na_col = "white", labels_col = as.character(df$name), display_numbers = F, show_colnames = T)



#### Fig H PRKCSH hetrozygous DEG ####################################################################################################
salmon_load(data_table = "/scratch/wpater/fig_h_table.txt", salmon_dir = "/scratch/wpater/analysis/liver_organoids/salmon_output/",
            DESeq_design = "Targeted", control = "Control", returnMode = c("vst","dds"), batch_correction = T, batch_group = "batch", k=2)
sample1 = "PRKCSH"
sample2 = "Control"
res = results(dds, contrast=c("Targeted", sample1, sample2))                         # compute DEG between sample 1 and sample 2
res=ensembl_id2symbol(res)
res[is.na(res$padj),] = 1
res = res[res$padj < 0.05,]
length(res$log2FoldChange[which(res$log2FoldChange > 0)]) #number upregulated
length(res$log2FoldChange[which(res$log2FoldChange < 0)]) # number downregulated
write.csv(res, file = "/scratch/wpater/PRKCSH_HET_v_Control.csv")
# no enriched pathways

vst=ensembl_id2symbol(vst)
final_samples = c("MDL_3.3","MDL_3.4","MDL_3.5","MDL_3.6", "EDG_n2", "EDG_n3", "EDG_n4","MDL_2-6")
vst= vst[,final_samples]
cbPalette = c(Control="#0000ff","PRKCSH"="#006400",na="#999999")
PCAplotR(vst,  colour = "Targeted" ,PCX=1, PCY=2, ntop = 500, palette = cbPalette)+
  ggtitle("Selected controls and heterozygous PRKCSH sample ", subtitle = "Top 500 highest variance genes of batch corrected, vst normalised counts")
# Save as 6.5 *5 inch pdf


#### Fig I Unknown mutation DEG ####################################################################################################
salmon_load(data_table = "/scratch/wpater/fig_i_table.txt", salmon_dir = "/scratch/wpater/analysis/liver_organoids/salmon_output/",
            DESeq_design = "Targeted", control = "Control", returnMode = c("vst","dds"), batch_correction = T, batch_group = "batch", k=2)
sample1 = "Unknown"
sample2 = "Control"
res = results(dds, contrast=c("Targeted", sample1, sample2))                         # compute DEG between sample 1 and sample 2
res=ensembl_id2symbol(res)
res[is.na(res$padj),] = 1
res = res[res$padj < 0.05,]
length(res$log2FoldChange[which(res$log2FoldChange > 0)]) #number upregulated
length(res$log2FoldChange[which(res$log2FoldChange < 0)]) # number downregulated
write.csv(res, file = "/scratch/wpater/Unknown_v_Control.csv")
# Some enriched pathways, but only 3 significant genes

vst=ensembl_id2symbol(vst)
final_samples = c("MDL_3.3","MDL_3.4","MDL_3.5","MDL_3.6", "EDG_n2", "EDG_n3", "EDG_n4","MDL_2-8")
vst= vst[,final_samples]
cbPalette = c(Control="#0000ff","Unknown"="#ff0000",na="#999999")
PCAplotR(vst,  colour = "Targeted" ,PCX=1, PCY=2, ntop = 500, palette = cbPalette)+
  ggtitle("Selected controls and unknown mutation sample ", subtitle = "Top 500 highest variance genes of batch corrected, vst normalised counts")
# Save as 6.5 *5 inch pdf

#final_genes = c() #
#vst= vst[final_genes,final_samples] # select rows with final_genes as rownames
#library(pheatmap)
#df = as.data.frame(colData(vst)[,c("Mutation","Sex", "Library_prep","name")]) # make dataframe for heatmap annotation
#colnames(df)[colnames(df) == "Library_prep"] = "Library prep"
#ann_colors = list(Mutation=c(None="#999999",PRKCSH="#006400"),
#                  "Library prep" = c("PolyA selection"="#009292","Ribo-Zero run 1"="#b66dff","Ribo-Zero run 2"="#db6d00"),
#                  Sex=c("Male"= "steelblue", "Female"= "plum")) # Colours for the heatmap annotation
#breakslist = seq(-2.1,2.0, by = 0.10) # scaling of the heatmap colours
#set.seed(10)
#pheatmap(assay(vst), cluster_rows=T, show_rownames=T, cluster_cols=T, breaks = breakslist,
#         annotation_col=df[,1:2], annotation_colors = ann_colors, scale = "row", color = colorRampPalette(c("blue", "white","red"))(length(breakslist)),
#         main = paste0("Unknown mutation enriched genes in PLD organoids"), na_col = "white", labels_col = as.character(df$name), display_numbers = F, show_colnames = F)




### Fig SF. ####################################################################################################
salmon_load(make_tx2gene = T, annotation = "/scratch/wpater/genomes_and_annotations/noY_gencode.v29.annotation.gtf")

salmon_load(data_table = "/scratch/wpater/fig_sc_table.txt", salmon_dir = "/scratch/wpater/analysis/liver_organoids/salmon_output/",
            DESeq_design = "Selection", control = "Control", returnMode = c("vst","dds"), batch_correction = T, batch_group = "batch", k=2)
sample1 = "PLD"
sample2 = "Control"
res = results(dds, contrast=c("Selection", sample1, sample2))                         # compute DEG between sample 1 and sample 2
res=ensembl_id2symbol(res)
res[is.na(res$padj),] = 1
res = res[res$padj < 0.05,]
length(res$log2FoldChange[which(res$log2FoldChange > 0)]) #number upregulated
length(res$log2FoldChange[which(res$log2FoldChange < 0)]) # number downregulated
write.csv(res, file = "/scratch/wpater/OVERALL_PLD_v_Control.csv")

vst=ensembl_id2symbol(vst)
final_samples = c("MDL_1-5","MDL_1-6","MDL_1-7","MDL_1-8","MDL_2-6","MDL_2-7","MDL_2-8","MDL_3.3","MDL_3.4","MDL_3.5","MDL_3.6","MDL_3.7","MDL_3.8",
                  "EDG_n2","EDG_n3","EDG_n4","EDG_8p6","EDG_9p6","EDG_10p4","EDG_12p5","EDG_16p6")
vst= vst[,final_samples]
cbPalette = c(Control="#0000ff",PLD="#ff0000",na="#999999")
PCAplotR(vst, shape = "Mutation", colour = "Selection" ,PCX=1, PCY=2, ntop = 500, palette = cbPalette)+
  ggtitle("All PLD samples and controls", subtitle = "Top 500 highest variance genes of batch corrected, vst normalised counts")
# Save as 6.5 *5 inch pdf

vst= vst[row.names(res),final_samples] # select rows with final_genes as rownames
library(pheatmap)
df = as.data.frame(colData(vst)[,c("Mutation","Sex", "Library_prep","name")]) # make dataframe for heatmap annotation
colnames(df)[colnames(df) == "Library_prep"] = "Library prep"
ann_colors = list(Mutation=c(None="#999999",PRKCSH="#006400","PRKCSH LOH"="#00ff00",Unknown="#ff0000","SEC63 LOH"="#ff00ff"),
                  "Library prep" = c("PolyA selection"="#009292","Ribo-Zero run 1"="#b66dff","Ribo-Zero run 2"="#db6d00"),
                  Sex=c("Male"= "steelblue", "Female"= "plum")) # Colours for the heatmap annotation
breakslist = seq(-2.1,2.0, by = 0.10) # scaling of the heatmap colours
set.seed(10)
pheatmap(assay(vst), cluster_rows=T, show_rownames=F, cluster_cols=T, breaks = breakslist,
         annotation_col=df[,1:2], annotation_colors = ann_colors, scale = "row", color = colorRampPalette(c("blue", "white","red"))(length(breakslist)),
         main = paste0("Differential genes between all PLD samples and controls"), na_col = "white", labels_col = as.character(df$name), 
         display_numbers = F, show_colnames = F)
# Save as 8.5 *6.5 inch pdf