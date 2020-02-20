#starting with the counts matrix derived in snRNAhackathon_UWdata_filtered.R:
library(Seurat)
library(scater)
library(batchelor)

counts <- dat[rowSums(dat != 0) >= 500,]
dim(counts)
write.csv(as.matrix(counts),file='/Users/relyanow/Desktop/mount2/users/relyanow/test_data/UW_sn_AD.csv')

dat <- CreateSeuratObject(counts = as.matrix(counts), project = "HighLow", min.cells = 1, min.features = 3)
counts<-0
dat

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
dat[["percent.mt"]] <- PercentageFeatureSet(dat, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(dat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(dat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

dat <- subset(dat, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 15)

dat <- NormalizeData(dat, normalization.method = "LogNormalize", scale.factor = 10000)

str(dat)

dat <- FindVariableFeatures(object = dat,nfeatures = 10000)

str(dat@assays$RNA@var.features)

dat <- ScaleData(dat,features=dat@assays$RNA@var.features)

dat <- RunPCA(dat,npcs=50)

# Examine and visualize PCA results a few different ways
print(dat[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(dat, dims = 1:2, reduction = "pca")

DimPlot(dat, reduction = "pca")

ElbowPlot(dat,50)

dat <- FindNeighbors(dat, dims = 1:25)
dat <- FindClusters(dat, resolution = 0.25)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
dat <- RunUMAP(dat, dims = 1:25)
dat <- RunTSNE(dat, dims = 1:25)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(dat, reduction = "umap")
DimPlot(dat, reduction = "tsne")

#saveRDS(dat, file = "UW_single_cell.rds")
dat<-readRDS(file = "UW_single_cell.rds")

length(rownames(dat))
counts <- dat@assays$RNA@counts

length(Labels$SEX)
head(Labels)



library(scran)

sce <- SingleCellExperiment(list(counts=counts))
clusters <- quickCluster(sce, min.size=100)
sce <- computeSumFactors(sce, cluster=clusters)

sce <- normalize(sce)
head(logcounts(sce[0:20,0:20]))
head(rownames(logcounts(sce)))

dim(logcounts(sce))
length(rownames(Labels))
head(Labels)
unique(Labels$Clinical.DX)

sce$projid = Labels$ids
sce$sex = Labels$SEX
sce$sex[Labels$SEX=='M']='male'
sce$sex[Labels$SEX!='M']='female'
sce$simpleDiagnosis[sce$projid==6672]='Control'
sce$simpleDiagnosis[sce$projid!=6672]='AD'
sce$pmi = Labels$PMI
sce$projid = as.character(Labels$ids)
sce$dataset = Labels$SEX
sce$dataset[sce$dataset=='male']='UW'
sce$dataset[sce$dataset!='female']='UW'

gene_names <-readRDS(file = "gene_names.rds")
length(gene_names)
genes=c()
for (gene in gene_names){
    genes = c(genes,which(rownames(logcounts(sce))==gene)[1])
}
length(genes)

saveRDS(sce, file = "UW_scran_normalized.rds")

sce = readRDS("UW_scran_normalized.rds")
Labels = data.frame(sce$projid)
Labels$projid = sce$projid
Labels$sex = sce$sex 
Labels$simpleDiagnosis = sce$simpleDiagnosis
Labels$pmi = sce$pmi 
Labels$dataset = sce$dataset
rownames(Labels) = colnames(logcounts(sce))

library(scater)
library(batchelor)
library(monocle3)

counts2 = logcounts(sce)
#counts2 = batchelor::fastMNN(counts2,batch=sce$projid)

dim(counts2)
gene_short_name <- data.frame(rownames(counts2))
rownames(gene_short_name) <- rownames(counts2)
gene_short_name$gene_short_name <- rownames(counts2)
head(gene_short_name)
length(gene_short_name)

cds <- new_cell_data_set(counts2,
                     cell_metadata = Labels,
                     gene_metadata = gene_short_name)


sum(is.na(cds$pmi))

cds$projid = as.numeric(cds$projid)
cds = preprocess_cds(cds, num_dim = 30,method="PCA",residual_model_formula_str="~projid+pmi",norm_method="size_only")
cds = reduce_dimension(cds)


plot_pc_variance_explained(cds)

cds$projid = as.character(cds$projid)
plot_cells(cds, color_cells_by="projid",cell_size=.1,label_cell_groups=0,show_trajectory_graph=FALSE)


plot_cells(cds, color_cells_by="sex",cell_size=.1,label_cell_groups=0,show_trajectory_graph=FALSE)

plot_cells(cds, color_cells_by="simpleDiagnosis",cell_size=.1,label_cell_groups=0,show_trajectory_graph=FALSE)

cds = cluster_cells(cds)
plot_cells(cds, color_cells_by="partition",cell_size=.1,show_trajectory_graph=FALSE)

plot_cells(cds, genes="CD74",cell_size=.1,show_trajectory_graph=FALSE)

counts<-0
counts2<-0
dat<-0
saveRDS(cds, file = "../data/UW_10X_RNA/UW_monocle_cds.rds")
gene_correlation <- cor(as.matrix(counts(cds)), method = c("pearson"))
dim(gene_correlation)

# chi square test for significant difference female/male ratio

cds_subset = cds[,partitions(cds)==8]
test = data.frame(clusters(cds_subset))
test$cluster = cds_subset$sex
test$cluster[inds] = 'Mic1'
test$cluster[test$cluster!='Mic1'] = 'Other'
test$sex = cds_subset$sex
test$cluster = as.factor(test$cluster)
test$sex = as.factor(test$sex)
tbl = table(test$cluster, test$sex)
library(MASS)
chisq.test(tbl) 

k=softConnectivity(datE=as.matrix(counts(cds),power=6)
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")

cds_subset = cds[,Labels$sex!="asdfmale"] # mic3 subcluster
cds_subset = cds_subset[,partitions(cds_subset)==9] # oligo
cds_subset = cds_subset[rowSums(exprs(cds_subset) != 0) >= 10,] # only keep genes non-zero in at least 10 cells
l = c("FTL","CAV1"
,"COX7A2"
,"RPS7"
,"B2M"
,"RPS27A"
,"HLA-A"
,"EEF1A1"
,"RPS16"
,"RPS2"
,"RPLP0"
,"RPL21"
,"BTF3"
,"RPS27"
,"RPL41"
,"PLEKHA6"
,"GSTP1"
,"RPL19"
,"RPS18"
,"RPS5"
,"HLA-C"
,"ZNF90"
,"H2AFZ"
,"DAD1"
,"POMP"
,"CD9"
,"RPS17"
,"UBXN1"
,"BAIAP2L1"
,"POLR2J"
,"PHGDH"
,"UBB")
inds = c()
for (gene in l){
    if ((gene %in% rownames(cds_subset))){
    inds = c(inds,which(rownames(cds_subset)==gene))
    }
}
inds
cds_subset = cds_subset[inds,]
cds_subset = preprocess_cds(cds_subset, num_dim = 5,method="PCA",norm_method="size_only",residual_model_formula_str="~projid")
cds_subset = reduce_dimension(cds_subset)
cds_subset = cluster_cells(cds_subset)
plot_cells(cds_subset, color_cells_by = "cluster", 
           label_groups_by_cluster=FALSE,label_leaves=FALSE,
           label_branch_points=FALSE,cell_size=2,label_cell_groups=0)
# cds_subset$test = as.numeric(clusters(cds_subset))
# cds_subset$test[clusters(cds_subset)==10]=1
# cds_subset$test[clusters(cds_subset)!=10]=0
# cds_subset$test = as.factor(cds_subset$test)
# cds_subset$test = t
plot_cells(cds_subset,label_cell_groups=0,genes=c("FTL","APOE","TPT1","RPL13"),cell_size=1,show_trajectory_graph=FALSE)
plot_cells(cds_subset,label_cell_groups=0,color_cells_by="partition",cell_size=1,show_trajectory_graph=FALSE)

# cds_subset@clusters$UMAP$partitions[cds_subset@clusters$UMAP$partitions == "2"] <- "1"
# cds_subset@clusters$UMAP$partitions[cds_subset@clusters$UMAP$partitions == "3"] <- "1"
# cds_subset <- learn_graph(cds_subset,use_partition=F)
# cds_subset$ftl = as.numeric(counts(cds_subset)[which(rownames(cds_subset)=="FTL"),])
# plot_cells(cds_subset, color_cells_by = "sex", 
#            label_groups_by_cluster=FALSE,label_leaves=FALSE,
#            label_branch_points=FALSE,cell_size=2,label_cell_groups=0)
# plot_cells(cds_subset, color_cells_by = "projid", 
#            label_groups_by_cluster=FALSE,label_leaves=FALSE,
#            label_branch_points=FALSE,cell_size=2,label_cell_groups=0)

# cds_subset = order_cells(cds_subset,root_cells=row.names(colData(cds_subset))[1])
# plot_cells(cds_subset, color_cells_by = "pseudotime", 
#            label_groups_by_cluster=FALSE,label_leaves=FALSE,
#            label_branch_points=FALSE,cell_size=2,label_cell_groups=0)

# gene_fits = fit_models(cds_subset, model_formula_str = "~test")
# fit_coefs = coefficient_table(gene_fits)
# fit_coefs <- subset(fit_coefs, term != "(Intercept)")

# fit_coefs <- subset(fit_coefs, q_value < 0.05)
# fit_coefs[order(fit_coefs$q_value),]
# write.csv(fit_coefs[order(fit_coefs$q_value),],file='/Users/relyanow/Documents/Sage/notebooks/figures/sn_Mathys/oli_UW.csv')


t = cds_subset$test

library(limma)
library(edgeR)
head(Labels$SEX)
cds$projid = sce$projid


l = c("FTL","SPP1","RPLP1","TMEM163","SLC11A1","RPLP2","RPL28","RPS20","C1QC","FTH1","RPL19","RPS15","RPL35",
"C1QB","RPL13","RPS19","RPS11","APOE","RPS27A","ACSL1","C1QA","RPL13A","PLEKHA7","RPS24","TMSB4X","RPL32","RPS2",
"VSIG4","RPS6","RPS8","RPL31","RPL27A","RPS16","RPS3","RPL10","TPT1","CST3","TMSB10","HLA-DRA","CD74","CD14",
"TXNRD1","CTSB","HLA-DRB1","ACTB","RPS9","GPX1","HIF1A","CYBA","AMBRA1","TYROBP","SLC2A5","LAPTM5","HLA-B","SIPA1L1",
"MS4A6A","TBC1D14","DPYD","TUBA1B","HCLS1","PDE4B","SRGN","EEF1A1","HSP90AA1","SAT1","EPB41L3","ADGRG1","PSAP",
"DENND3",
"RB1",
"GAPDH",
"PHC2",
"HSPA1A",
"RNF149",
"HCK",
"SH3TC1",
"NUMB")
inds = c()
for (gene in l){
    if ((gene %in% rownames(cds))){
    inds = c(inds,which(rownames(cds)==gene))
    }
}
inds
cds_subset = cds[inds,] # mic3 subcluster
cds_subset = cds_subset[,partitions(cds_subset)==8] # mic3 subcluster
#'6774' '6802' '6829' '6845' '6874'
#cds_subset = cds_subset[,c(which(cds_subset$projid=='6774'),which(cds_subset$projid=='6829'),which(cds_subset$projid=='6874'))]
#cds_subset = cds_subset[rowSums(exprs(cds_subset) != 0) >= 20,] # only keep genes non-zero in at least 10 cells
cds_subset = preprocess_cds(cds_subset, num_dim = 20,method="PCA",residual_model_formula_str="~projid+sex",norm_method="size_only")
cds_subset = reduce_dimension(cds_subset)
cds_subset = cluster_cells(cds_subset)
#cds_subset = cds_subset[,cds_subset$sex=="female"] # mic3 subcluster
plot_pc_variance_explained(cds_subset)
plot_cells(cds_subset, color_cells_by = "cluster", 
           label_groups_by_cluster=FALSE,label_leaves=FALSE,
           label_branch_points=FALSE,cell_size=2,label_cell_groups=1)
pdf("/Users/relyanow/Documents/Sage/ROSMAP_sn_writeup/figures/UW_mic_tsne_using_mic1_genes.pdf") 
plot_cells(cds_subset, genes = c("FTL","APOE","TPT1","RPL13"), 
           label_groups_by_cluster=FALSE,label_leaves=FALSE,
           label_branch_points=FALSE,cell_size=1,label_cell_groups=0)
dev.off()
plot_cells(cds_subset, color_cells_by = "projid", 
           label_groups_by_cluster=FALSE,label_leaves=FALSE,
           label_branch_points=FALSE,cell_size=2,label_cell_groups=0)
plot_cells(cds_subset, color_cells_by = "sex", 
           label_groups_by_cluster=FALSE,label_leaves=FALSE,
           label_branch_points=FALSE,cell_size=2,label_cell_groups=0)

dput(as.numeric(c(which(clusters(cds_subset)==4),which(clusters(cds_subset)==11),
                which(clusters(cds_subset)==3),which(clusters(cds_subset)==15))))




library(limma)
library(edgeR)
head(Labels$SEX)
cds$projid = sce$projid

inds= c(125, 186, 231, 267, 268, 272, 342, 369, 385, 391, 431, 437, 
469, 574, 692, 819, 859, 1194, 1202, 1209, 1214, 1235, 1243, 
1245, 1252, 1267, 1283, 1298, 1307, 1310, 1312, 1316, 1318, 1328, 
1333, 1338, 1342, 1347, 1361, 1368, 1370, 1385, 1393, 1440, 1445, 
1453, 1454, 1458, 1467, 1476, 1479, 1482, 1503, 1518, 1527, 1546, 
1555, 1564, 1585, 1601, 1605, 1651, 1655, 1683, 1688, 1695, 1707, 
1718, 1755, 1861, 1864, 1870, 1872, 1873, 1889, 1893, 1897, 1906, 
2013, 2046, 2052, 2076, 2132, 2138, 3, 6, 20, 33, 42, 56, 61, 
65, 66, 67, 76, 80, 85, 134, 156, 162, 172, 173, 177, 191, 235, 
241, 253, 255, 292, 296, 297, 300, 303, 304, 313, 319, 320, 321, 
324, 330, 338, 341, 348, 349, 359, 427, 458, 460, 462, 478, 486, 
487, 493, 496, 591, 635, 638, 640, 674, 709, 855, 875, 887, 905, 
918, 938, 998, 1006, 1027, 1052, 1065, 1110, 1115, 1120, 1125, 
1138, 1142, 1148, 1165, 1174, 1180, 1182, 1188, 1191, 1197, 1203, 
1212, 1217, 1227, 1231, 1233, 1236, 1248, 1254, 1263, 1265, 1275, 
1297, 1306, 1309, 1317, 1346, 1358, 1362, 1366, 1369, 1376, 1379, 
1380, 1394, 1401, 1402, 1403, 1405, 1409, 1415, 1419, 1421, 1426, 
1438, 1444, 1450, 1470, 1473, 1475, 1481, 1506, 1508, 1513, 1515, 
1530, 1537, 1541, 1542, 1568, 1569, 1571, 1574, 1575, 1581, 1593, 
1600, 1602, 1604, 1611, 1617, 1622, 1631, 1656, 1661, 1668, 1671, 
1690, 1704, 1709, 1715, 1716, 1720, 1728, 1740, 1786, 1793, 1818, 
1821, 1822, 1829, 1830, 1834, 1842, 1849, 1856, 1866, 1867, 1877, 
1878, 1880, 1883, 1886, 1896, 1902, 1916, 1919, 1923, 1924, 1954, 
1959, 1961, 1965, 1975, 1995, 2000, 2014, 2020, 2029, 2035, 2038, 
2047, 2067, 2072, 2074, 2077, 2086, 2099, 2107, 2109, 2118, 2124, 
2125, 2130, 2136, 2137, 11, 21, 29, 60, 81, 82, 90, 94, 127, 
128, 129, 130, 140, 160, 176, 180, 194, 201, 237, 254, 262, 263, 
275, 280, 283, 289, 299, 308, 311, 323, 367, 372, 390, 399, 408, 
434, 465, 473, 481, 482, 485, 489, 494, 509, 510, 545, 553, 571, 
586, 611, 624, 639, 648, 649, 659, 700, 736, 740, 766, 767, 778, 
816, 838, 842, 850, 899, 903, 906, 914, 922, 925, 929, 930, 940, 
947, 967, 971, 973, 995, 1004, 1008, 1010, 1011, 1034, 1037, 
1049, 1057, 1063, 1085, 1089, 1094, 1112, 1123, 1129, 1154, 1331, 
1422, 1424, 1429, 1437, 1509, 1512, 1514, 1532, 1547, 1550, 1629, 
1632, 1640, 1645, 1678, 1702, 1706, 1710, 1722, 1729, 1734, 1742, 
1753, 1754, 1769, 1779, 1792, 1795, 1824, 1836, 1840, 1844, 1875, 
1891, 1928, 1949, 1957, 1991, 1994, 2007, 2025, 2045, 2063, 2116, 
2123, 2133, 26, 45, 68, 143, 171, 228, 232, 239, 302, 375, 384, 
392, 407, 423, 457, 497, 583, 598, 604, 615, 637, 650, 694, 696, 
724, 725, 745, 757, 772, 807, 820, 828, 841, 848, 849, 852, 932, 
950, 958, 961, 1012, 1015, 1053, 1069, 1079, 1093, 1096, 1097, 
1124, 1139, 1152, 1155, 1156, 1157, 1159, 1199, 1350, 1351, 1374, 
1375, 1384, 1567, 1594, 1616, 1652, 1665, 1666, 1771, 1797, 1814, 
1815, 1904, 1911, 1918, 1936, 1941, 2005, 2048, 2065, 2078, 2082, 
2093, 2098, 2114, 2127)



cds_subset = cds[,partitions(cds)==8] 
'all'
length(colnames(cds)[cds$sex=="male"] )/length(colnames(cds) )
length(colnames(cds)[cds$sex=="female"])/length(colnames(cds) )
'mic'
length(colnames(cds_subset) )
length(colnames(cds_subset)[cds_subset$sex=="male"] )/length(colnames(cds_subset) )
length(colnames(cds_subset)[cds_subset$sex=="female"])/length(colnames(cds_subset) )

'mic1'
cds_subset = cds_subset[,inds]
length(colnames(cds_subset) )
length(colnames(cds_subset)[cds_subset$sex=="male"] )/length(colnames(cds_subset) )
length(colnames(cds_subset)[cds_subset$sex=="female"])/length(colnames(cds_subset) )

cds_subset = cds[,] 
cds_subset = cds_subset[,partitions(cds_subset)==8] 
cds_subset = cds_subset[rowSums(exprs(cds_subset) != 0) >= 20,] # only keep genes non-zero in at least 10 cells
#cds_subset = preprocess_cds(cds_subset, num_dim = 10,method="PCA",residual_model_formula_str="~projid",norm_method="size_only")

cds_subset$mic1 = cds_subset$projid
cds_subset$mic1[inds] = 'Mic1'
cds_subset$mic1[cds_subset$mic1!='Mic1']='Mic0'

#cds_subset = cds_subset[,cds_subset$sex=="female"] 
#cds_subset = preprocess_cds(cds_subset, num_dim = 20,method="PCA",residual_model_formula_str="~pmi",norm_method="size_only")
cds_subset = reduce_dimension(cds_subset)
cds_subset = cluster_cells(cds_subset)


gene_fits = fit_models(cds_subset, model_formula_str = "~mic1+pmi")
fit_coefs = coefficient_table(gene_fits)
fit_coefs[which(fit_coefs$gene_short_name=='FTL'),]
fit_coefs[which(fit_coefs$gene_short_name=='APOE'),]
fit_coefs[which(fit_coefs$gene_short_name=='NACA'),]
fit_coefs[which(fit_coefs$gene_short_name=='HLA-B'),]
fit_coefs[which(fit_coefs$gene_short_name=='SPP1'),]
fit_coefs <- subset(fit_coefs, term != "(Intercept)")
fit_coefs <- subset(fit_coefs, term != "projid")
fit_coefs <- subset(fit_coefs, term != "pmi")
fit_coefs <- subset(fit_coefs, status == "OK")

fit_coefs2 <- subset(fit_coefs, q_value < 0.05)
fit_coefs2[order(fit_coefs2$q_value),]
#write.csv(fit_coefs2[order(fit_coefs2$q_value),],file='/Users/relyanow/Documents/Sage/notebooks/figures/sn_Mathys/mic0-1_subcluster_UW.csv')

library(EnhancedVolcano)
fit_coefs$pvalue=fit_coefs$q_value
fit_coefs$log2FoldChange=fit_coefs$normalized_effect
#pdf("/Users/relyanow/Documents/Sage/ROSMAP_sn_writeup/figures/volcano_UW.pdf") 
EnhancedVolcano(fit_coefs,
    lab = fit_coefs$gene_short_name,
    x = 'log2FoldChange',
    y = 'pvalue',
    xlim = c(-4,4),
    xlab = bquote(~Log[2]~ 'normalized effect'),
    shape = c(17, 35, 17, 18),
    pCutoff = .05,
    FCcutoff = 1,
    colAlpha = 1,
    transcriptPointSize = 3.0,
    transcriptLabSize = 4.0,
    transcriptLabCol = 'black',
    transcriptLabFace = 'bold',
    boxedlabels = FALSE,
    legend=c('NS','Log (base 2) normalized effect','P value',
      'P value & Log (base 2) normalized effect'),
    legendPosition = 'right',
    legendLabSize = 14,
    legendIconSize = 5.0)
#dev.off()

# library
library(ggplot2)
library(viridis)

 
# create a dataset
specie <- c(rep("Human \n(Mathys et al.)" , 2) , rep("Human" , 2) ,rep("Mouse \n(Figerio et al.)",2) )#"Mouse (Figerio et al.)"
Sex <- rep(c("Female" , "Male") , 3)
Percent <- c(59,41,78,22,54,46)
data <- data.frame(specie,Sex,Percent)
#pdf("/Users/relyanow/Documents/Sage/ROSMAP_sn_writeup/figures/mic1_percentage.pdf") 
# Grouped
ggplot(data, aes(fill=Sex, y=Percent, x=specie)) + 
    geom_bar(position="dodge", stat="identity")+ theme_classic()+theme(text=element_text(size=21))
#dev.off()

plot_cells(cds_subset, color_cells_by = "cluster", 
           label_groups_by_cluster=FALSE,label_leaves=FALSE,
           label_branch_points=FALSE,cell_size=2,label_cell_groups=0)

plot_cells(cds, genes=c("CD74"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,reduction_method="UMAP",cell_size=.5)
plot_cells(cds_subset, genes=c("APOE","CD74","FTL","ACTB","CST3","RPL13","IFI44L","SAT1","LYN"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,reduction_method="UMAP",cell_size=1)
unique(cds_subset$projid)
cds_subset$test = cds_subset$projid
cds_subset$test[cds_subset$test!='asdf']='Other'
cds_subset$test[cds_subset$projid=='6672']='test'
plot_cells(cds_subset, color_cells_by = "test",
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,reduction_method="UMAP",cell_size=1)
cds_subset$test = cds_subset$projid
cds_subset$test[cds_subset$test!='asdf']='Other'
cds_subset$test[cds_subset$projid=='6687']='test'
plot_cells(cds_subset, color_cells_by = "test",
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,reduction_method="UMAP",cell_size=1)
cds_subset$test = cds_subset$projid
cds_subset$test[cds_subset$test!='asdf']='Other'
cds_subset$test[cds_subset$projid=='6726']='test'
plot_cells(cds_subset, color_cells_by = "test",
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,reduction_method="UMAP",cell_size=1)
cds_subset$test = cds_subset$projid
cds_subset$test[cds_subset$test!='asdf']='Other'
cds_subset$test[cds_subset$projid=='6774']='test'
plot_cells(cds_subset, color_cells_by = "test",
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,reduction_method="UMAP",cell_size=1)
#'6774' '6802' '6829' '6845' '6874'
#'6672' '6687' '6726' '6774' '6802' '6829' '6845' '6874'
cds_subset$test = cds_subset$projid
cds_subset$test[cds_subset$test!='asdf']='Other'
cds_subset$test[cds_subset$projid=='6802']='test'
plot_cells(cds_subset, color_cells_by = "test",
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,reduction_method="UMAP",cell_size=1)
cds_subset$test = cds_subset$projid
cds_subset$test[cds_subset$test!='asdf']='Other'
cds_subset$test[cds_subset$projid=='6829']='test'
plot_cells(cds_subset, color_cells_by = "test",
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,reduction_method="UMAP",cell_size=1)
cds_subset$test = cds_subset$projid
cds_subset$test[cds_subset$test!='asdf']='Other'
cds_subset$test[cds_subset$projid=='6845']='test'
plot_cells(cds_subset, color_cells_by = "test",
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,reduction_method="UMAP",cell_size=1)
cds_subset$test = cds_subset$projid
cds_subset$test[cds_subset$test!='asdf']='Other'
cds_subset$test[cds_subset$projid=='6874']='test'
plot_cells(cds_subset, color_cells_by = "test",
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,reduction_method="UMAP",cell_size=1)


# plot_cells(cds_subset, color_cells_by="partition",cell_size=1,show_trajectory_graph=FALSE)
# cds_subset$partition = as.factor(partitions(cds_subset))
# gene_fits = fit_models(cds_subset, model_formula_str = "~partition")
# fit_coefs = coefficient_table(gene_fits)
# fit_coefs$q_value[which(fit_coefs$gene_short_name=='SLC11A1')]
# fit_coefs <- subset(fit_coefs, q_value < 0.05)
# fit_coefs <- subset(fit_coefs, term != "(Intercept)")
# fit_coefs <- subset(fit_coefs, term != "projid")
# fit_coefs <- subset(fit_coefs, term != "educ")
# fit_coefs <- subset(fit_coefs, term != "apoe_genotype44")
fit_coefs[order(fit_coefs$q_value)[0:100],]

fit_coefs[which(fit_coefs$gene_short_name=='APOE')[1],]
fit_coefs[which(fit_coefs$gene_short_name=='FTL')[1],]
fit_coefs[which(fit_coefs$gene_short_name=='SLC11A1')[1],]
fit_coefs[which(fit_coefs$gene_short_name=='C1QC')[1],]
fit_coefs[which(fit_coefs$gene_short_name=='RPS20')[1],]



match(markers,dat@assays$RNA@var.features)
#markers=c('PDGFRA', 'COL9A1', 'ST18' ,'TGFBR2' ,'GFAP' ,'FGFR3', 'CX3CR1', 'C3', 'CD74' ,'GAD2', 'GRIP2', 'SLC17A7' ,'SATB2' ,'NRGN' ,'GABRA1' ,'PROM1' )


# find all markers of cluster 0,1
#cluster1.markers <- FindMarkers(dat, ident.1 = c(0,1), min.pct = 0.25)
#head(cluster1.markers, n = 20)

markers = c('PDGFRA','COL9A1','TMEM63A','ST18','TGFBR2',
            'FGFR3','CX3CR1','C3','CD74','GAD2','GRIP2','SLC17A7','SATB2','NRGN','GABRG2','GABRA1','PROM1','GABRG1','NEUROD6')
length(markers)
names=c('OPC','OPC','Olig','Olig','Endo',
                    'Astro','Mglia','Mglia','Mglia','GABA','GABA',
                    'Glut','Glut','Glut','Neurons','Neurons','BSC','Neurons','Neurons')
length(names)
clusts=c()
for (marker in rownames(head(cluster1.markers,n=100))){
    if (length(which(marker %in% markers))>0){
        clusts = c(clusts,names[which(marker %in% markers)])
    }
}
clusts

new.cluster.ids <- c()
markers
for (i in 0:20){
    #####
    num=i
    a.markers <- FindMarkers(dat, ident.1 = c(num), min.pct = 0.25,features=markers)
    a.markers<-a.markers[a.markers$avg_logFC>0,]
    clusts=c()
    print(num)
    print(head(a.markers,n=2))
    for (marker in rownames(head(a.markers,n=1))){
        if (length(which(markers == marker))>0){
            clusts = c(clusts,names[which(markers==marker)])
        }
    }
    new.cluster.ids <- c(new.cluster.ids,clusts)

}
new.cluster.ids
names(new.cluster.ids) <- levels(dat)
dat <- RenameIdents(dat, new.cluster.ids)
DimPlot(dat, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


subset <- SubsetData(object = dat, cells = colnames(dat)[Labels$ids == 6672])
types=c('OPC','Olig','Endo','Astro','Mglia','GABA','Glut','Neurons','BSC')
for (t in types){
    print(c(t,as.character(length(subset@active.ident[subset@active.ident == t]) / length(subset@active.ident))))
}
DimPlot(subset, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


subset <- SubsetData(object = dat, cells = colnames(dat)[Labels$ids == 6687])
types=c('OPC','Olig','Endo','Astro','Mglia','GABA','Glut','Neurons','BSC')
for (t in types){
    print(c(t,as.character(length(subset@active.ident[subset@active.ident == t]) / length(subset@active.ident))))
}
DimPlot(subset, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


subset <- SubsetData(object = dat, cells = colnames(dat)[Labels$ids == 6726])
types=c('OPC','Olig','Endo','Astro','Mglia','GABA','Glut','Neurons','BSC')
for (t in types){
    print(c(t,as.character(length(subset@active.ident[subset@active.ident == t]) / length(subset@active.ident))))
}
DimPlot(subset, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


subset <- SubsetData(object = dat, cells = colnames(dat)[Labels$ids == 6774])
types=c('OPC','Olig','Endo','Astro','Mglia','GABA','Glut','Neurons','BSC')
for (t in types){
    print(c(t,as.character(length(subset@active.ident[subset@active.ident == t]) / length(subset@active.ident))))
}
DimPlot(subset, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


subset <- SubsetData(object = dat, cells = colnames(dat)[Labels$ids == 6802])
types=c('OPC','Olig','Endo','Astro','Mglia','GABA','Glut','Neurons','BSC')
for (t in types){
    print(c(t,as.character(length(subset@active.ident[subset@active.ident == t]) / length(subset@active.ident))))
}
DimPlot(subset, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


subset <- SubsetData(object = dat, cells = colnames(dat)[Labels$ids == 6829])
types=c('OPC','Olig','Endo','Astro','Mglia','GABA','Glut','Neurons','BSC')
for (t in types){
    print(c(t,as.character(length(subset@active.ident[subset@active.ident == t]) / length(subset@active.ident))))
}
DimPlot(subset, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


subset <- SubsetData(object = dat, cells = colnames(dat)[Labels$ids == 6845])
types=c('OPC','Olig','Endo','Astro','Mglia','GABA','Glut','Neurons','BSC')
for (t in types){
    print(c(t,as.character(length(subset@active.ident[subset@active.ident == t]) / length(subset@active.ident))))
}
DimPlot(subset, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


subset <- SubsetData(object = dat, cells = colnames(dat)[Labels$ids == 6874])
types=c('OPC','Olig','Endo','Astro','Mglia','GABA','Glut','Neurons','BSC')
for (t in types){
    print(c(t,as.character(length(subset@active.ident[subset@active.ident == t]) / length(subset@active.ident))))
}
DimPlot(subset, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

library(monocle3)
dim(dat)
#Extract data, phenotype data, and feature data from the SeuratObject
inds <- match(dat@assays$RNA@var.features,rownames(dat))
counts <- as(as.matrix(dat@assays$RNA@data[inds,]), 'sparseMatrix')
gene_short_name<-as.matrix(rownames(dat)[inds])
nrow(counts)
nrow(gene_short_name)
rownames(counts)<-as.vector(gene_short_name)
gene_short_name <- data.frame(rownames(counts))
rownames(gene_short_name) <- rownames(counts)
gene_short_name$gene_short_name <- rownames(counts)
#counts <- sweep(as.matrix(counts2), MARGIN=2, colSums(as.matrix(counts)), `/`)*median(colSums(as.matrix(counts2)))
cds <- new_cell_data_set(counts2,
                     cell_metadata = Labels,
                     gene_metadata = gene_short_name)
head(colSums(counts))

colData(cds)$ident = as.vector(dat@active.ident)

cds = preprocess_cds(cds, num_dim = 50,method="PCA")
cds = reduce_dimension(cds)


plot_pc_variance_explained(cds)

plot_cells(cds, color_cells_by="ePRS",cell_size=.1,label_cell_groups=0,show_trajectory_graph=FALSE)

plot_cells(cds, color_cells_by="ident",cell_size=.1,label_cell_groups=0,show_trajectory_graph=FALSE)

colData(cds)$ids = as.character(colData(cds)$ids)
plot_cells(cds, color_cells_by="ids",cell_size=.1,label_cell_groups=0,show_trajectory_graph=FALSE)

str(dat)

plot_cells(cds, color_cells_by="SEX",cell_size=.1,label_cell_groups=0,show_trajectory_graph=FALSE)

plot_cells(cds, color_cells_by="Clinical.DX",cell_size=.1,label_cell_groups=0,show_trajectory_graph=FALSE)

plot_cells(cds, color_cells_by="PMI",cell_size=.1,label_cell_groups=0,show_trajectory_graph=FALSE)

plot_cells(cds, color_cells_by="Size_Factor",cell_size=.1,label_cell_groups=0,show_trajectory_graph=FALSE)
