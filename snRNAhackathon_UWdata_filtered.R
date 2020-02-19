library(readxl)
library(Matrix)
library(synapser)
library(rqdatatable)
#to login to synapse manually: synLogin("username", "password")
synLogin()

system( paste0('tar -xvf ', synapser::synGet('syn21614190')$path) )
system( 'gunzip home/dnanexus/6672/outs/filtered_feature_bc_matrix/matrix.mtx.gz' )
system( 'gunzip home/dnanexus/6672/outs/filtered_feature_bc_matrix/features.tsv.gz')
system( 'gunzip home/dnanexus/6672/outs/filtered_feature_bc_matrix/barcodes.tsv.gz')
dat_6672 <- readMM('home/dnanexus/6672/outs/filtered_feature_bc_matrix/matrix.mtx')
tags_6672  <- read.delim('home/dnanexus/6672/outs/filtered_feature_bc_matrix/barcodes.tsv', header=FALSE)
#the gene names file is identical for all patients
gene_metadata <- read.delim('home/dnanexus/6672/outs/filtered_feature_bc_matrix/features.tsv', header=FALSE)


system( paste0('tar -xvf ', synapser::synGet('syn21614188')$path) )
system( 'gunzip home/dnanexus/6687/outs/filtered_feature_bc_matrix/matrix.mtx.gz' )
system( 'gunzip home/dnanexus/6687/outs/filtered_feature_bc_matrix/barcodes.tsv.gz')
dat_6687 <- readMM('home/dnanexus/6687/outs/filtered_feature_bc_matrix/matrix.mtx')
tags_6687  <- read.delim('home/dnanexus/6687/outs/filtered_feature_bc_matrix/barcodes.tsv', header=FALSE)

system( paste0('tar -xvf ', synapser::synGet('syn21614184')$path) )
system( 'gunzip home/dnanexus/6726/outs/filtered_feature_bc_matrix/matrix.mtx.gz' )
system( 'gunzip home/dnanexus/6726/outs/filtered_feature_bc_matrix/barcodes.tsv.gz')
dat_6726 <- readMM('home/dnanexus/6726/outs/filtered_feature_bc_matrix/matrix.mtx')
tags_6726  <- read.delim('home/dnanexus/6726/outs/filtered_feature_bc_matrix/barcodes.tsv', header=FALSE)

system( paste0('tar -xvf ', synapser::synGet('syn21614182')$path) )
system( 'gunzip home/dnanexus/6774/outs/filtered_feature_bc_matrix/matrix.mtx.gz' )
system( 'gunzip home/dnanexus/6774/outs/filtered_feature_bc_matrix/barcodes.tsv.gz')
dat_6774 <- readMM('home/dnanexus/6774/outs/filtered_feature_bc_matrix/matrix.mtx')
tags_6774  <- read.delim('home/dnanexus/6774/outs/filtered_feature_bc_matrix/barcodes.tsv', header=FALSE)

system( paste0('tar -xvf ', synapser::synGet('syn21614179')$path) )
system( 'gunzip home/dnanexus/6802/outs/filtered_feature_bc_matrix/matrix.mtx.gz' )
system( 'gunzip home/dnanexus/6802/outs/filtered_feature_bc_matrix/barcodes.tsv.gz')
dat_6802 <- readMM('home/dnanexus/6802/outs/filtered_feature_bc_matrix/matrix.mtx')
tags_6802  <- read.delim('home/dnanexus/6802/outs/filtered_feature_bc_matrix/barcodes.tsv', header=FALSE)

system( paste0('tar -xvf ', synapser::synGet('syn21614178')$path) )
system( 'gunzip home/dnanexus/6829/outs/filtered_feature_bc_matrix/matrix.mtx.gz' )
system( 'gunzip home/dnanexus/6829/outs/filtered_feature_bc_matrix/barcodes.tsv.gz')
dat_6829 <- readMM('home/dnanexus/6829/outs/filtered_feature_bc_matrix/matrix.mtx')
tags_6829  <- read.delim('home/dnanexus/6829/outs/filtered_feature_bc_matrix/barcodes.tsv', header=FALSE)

system( paste0('tar -xvf ', synapser::synGet('syn21614192')$path) )
system( 'gunzip home/dnanexus/6845/outs/filtered_feature_bc_matrix/matrix.mtx.gz' )
system( 'gunzip home/dnanexus/6845/outs/filtered_feature_bc_matrix/barcodes.tsv.gz')
dat_6845 <- readMM('home/dnanexus/6845/outs/filtered_feature_bc_matrix/matrix.mtx')
tags_6845  <- read.delim('home/dnanexus/6845/outs/filtered_feature_bc_matrix/barcodes.tsv', header=FALSE)

system( paste0('tar -xvf ', synapser::synGet('syn21614191')$path) )
system( 'gunzip home/dnanexus/6874/outs/filtered_feature_bc_matrix/matrix.mtx.gz' )
dat_6874 <- readMM('home/dnanexus/6874/outs/filtered_feature_bc_matrix/matrix.mtx')
system( 'gunzip home/dnanexus/6874/outs/filtered_feature_bc_matrix/barcodes.tsv.gz')
tags_6874  <- read.delim('home/dnanexus/6874/outs/filtered_feature_bc_matrix/barcodes.tsv', header=FALSE)

#need to add column names and row names to each matrix; row names are the gene names, cols are barcodes + patientID
#only need column 2 of gene_meta1 (gene names)
#delete extra column in gene metadata:
gene_metadata <- subset(gene_metadata, select=c(1,2))
colnames(gene_metadata)[2] <- "gene_short_name"
#need to make duplicated gene short names unique
gene_metadata$gene_short_name <- make.names(gene_metadata$gene_short_name, unique=TRUE)

#each rowname is a gene, from the features file:
rownames(dat_6672) <- gene_metadata[,2]
rownames(dat_6687) <- gene_metadata[,2]
rownames(dat_6726) <- gene_metadata[,2]
rownames(dat_6774) <- gene_metadata[,2]
rownames(dat_6802) <- gene_metadata[,2]
rownames(dat_6829) <- gene_metadata[,2]
rownames(dat_6845) <- gene_metadata[,2]
rownames(dat_6874) <- gene_metadata[,2]

#each column is one of the barcodes, or cell tags:
#colnames(barcodes1)[1] <- "tags"
colnames(dat_6672) <- tags_6672[,1]
colnames(dat_6687) <- tags_6687[,1]
colnames(dat_6726) <- tags_6726[,1]
colnames(dat_6774) <- tags_6774[,1]
colnames(dat_6802) <- tags_6802[,1]
colnames(dat_6829) <- tags_6829[,1]
colnames(dat_6845) <- tags_6845[,1]
colnames(dat_6874) <- tags_6874[,1]

colnames(dat_6672) <- paste(colnames(dat_6672), "6672", sep="_")
colnames(dat_6687) <- paste(colnames(dat_6687), "6687", sep="_")
colnames(dat_6726) <- paste(colnames(dat_6726), "6726", sep="_")
colnames(dat_6774) <- paste(colnames(dat_6774), "6774", sep="_")
colnames(dat_6802) <- paste(colnames(dat_6802), "6802", sep="_")
colnames(dat_6829) <- paste(colnames(dat_6829), "6829", sep="_")
colnames(dat_6845) <- paste(colnames(dat_6845), "6845", sep="_")
colnames(dat_6874) <- paste(colnames(dat_6874), "6874", sep="_")

#join the individual patient matrices into one
dat <- cbind2(dat_6672,dat_6687)
dat <- cbind2(dat,dat_6726)
dat <- cbind2(dat,dat_6774)
dat <- cbind2(dat,dat_6802)
dat <- cbind2(dat,dat_6829)
dat <- cbind2(dat,dat_6845)
dat <- cbind2(dat,dat_6874)
dim(dat)
head(colnames(dat))
head(rownames(dat))

#read in metadata
syn21598855 <- synapser::synGet(entity='syn21598855')
metadata <- read_excel(syn21598855$path)
head(metadata)
metadata <- as.data.frame(metadata)
names(metadata)[names(metadata) == "Sample ID"] <- "ids"
names(metadata)[names(metadata) == "Clinical DX"] <- "Clinical.DX"
names(metadata)[names(metadata) == "path DX"] <- "path.DX"
names(metadata)[names(metadata) == "Apo E"] <- "ApoE"
metadata$ids <- as.factor(metadata$ids)
head(metadata)

#this count threshold can be changed
counts <- dat[rowSums(dat != 0) >= 250,]
dim(counts)

gene_short_name <- data.frame(rownames(counts))
rownames(gene_short_name) <- rownames(counts)
gene_short_name$gene_short_name <- rownames(counts)
head(gene_short_name)
dim(gene_short_name)

Labels = data.frame(colnames(counts))
Labels$ids <- ifelse(grepl("6672", Labels$colnames.counts.)==T, "6672",
                     ifelse(grepl("6687", Labels$colnames.counts.)==T, "6687",
                            ifelse(grepl("6726", Labels$colnames.counts.)==T, "6726",
                                   ifelse(grepl("6774", Labels$colnames.counts.)==T, "6774",
                                          ifelse(grepl("6802", Labels$colnames.counts.)==T, "6802",
                                                 ifelse(grepl("6829", Labels$colnames.counts.)==T, "6829",
                                                        ifelse(grepl("6845", Labels$colnames.counts.)==T, "6845", "6874")))))))


Labels <- natural_join(Labels, metadata, 
                       by = "ids",
                       jointype = "FULL")
rownames(Labels) = colnames(counts)
head(Labels)
dim(Labels)


#preprocessing and creating a monocle object
# must first unload synapser because causes multiple definitions of S4Vectors
detach("package:synapser", unload=TRUE)
unloadNamespace("PythonEmbedInR") 

library(monocle3)
library(scater)
library(ggplot2)
library(gridExtra)

#create monocle object:
cds_uw <- new_cell_data_set(counts,
                            cell_metadata = Labels,
                            gene_metadata = gene_short_name)

### preprocessing and reduce dimensionality
### preprocessing includes library size normalization and regressing out pmi
cds_uw = preprocess_cds(cds_uw, num_dim = 30,method="PCA", norm_method="log", residual_model_formula_str="~PMI")
cds_uw = reduce_dimension(cds_uw)
cds_uw = cluster_cells(cds_uw)
plot_pc_variance_explained(cds_uw)

cds_uw$Diagnosis = cds_uw$Clinical.DX
cds_uw$Sex = cds_uw$Sex
cds_uw$Samples = cds_uw$ids

cds_uw$Sex = Labels$SEX
head(cds_uw$Sex)

p1<-plot_cells(cds_uw, color_cells_by="partition",cell_size=.001,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 10))+theme(legend.position = "none")
p2<-plot_cells(cds_uw, color_cells_by="Sex",cell_size=.001,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 10))
p3<-plot_cells(cds_uw, color_cells_by="Diagnosis",cell_size=.001,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 8))
p4<-plot_cells(cds_uw, color_cells_by="Samples",cell_size=.001,label_cell_groups=0,show_trajectory_graph=FALSE)
p6<-plot_cells(cds_uw, color_cells_by="cluster",cell_size=.1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm") )+theme(legend.position = "none")
#pdf(paste0("/Users/relyanow/Documents/Sage/ROSMAP_sn_writeup/figures/UW_sn_summary.pdf"))
grid.arrange(arrangeGrob(p1,p6, ncol=2),arrangeGrob(p3,p2,ncol=2),p4, heights=c(2,2,4), ncol=1)
#dev.off()



#read in marker genes from mathys et al (supplemental table 8)
mathy_marker_genes <- read.csv('~/mathys_marker_genes.csv', header=TRUE)
mathy_marker_genes$adj.pvals = as.numeric(as.character(mathy_marker_genes$adj.pvals))
#mathy_marker_genes <- mathy_marker_genes[mathy_marker_genes$adj.pvals<1e-20,]
head(mathy_marker_genes)

genes<-c()
for (gene in unique(c(as.vector(mathy_marker_genes$gene.name),c("SYT1","SNAP25","GRIN1","GAD1","GAD2","SLC17A7","CAMK2A","NRGN","AQP4",
                                                                "GFAP","MBP","MOBP","PLP1","PDGFRA","VCAN","CD74","CSF1R","C3","FLT1","CLDN5")))){
  if (gene %in% rownames(cds_uw)){
    genes <- c(genes,which(rownames(cds_uw)==gene))
  }
}
length(genes)

cds_subset = cds_uw[genes,]
cds_subset = preprocess_cds(cds_subset, num_dim = 30,method="PCA",residual_model_formula_str="~PMI",norm_method="size_only")
cds_subset = reduce_dimension(cds_subset)
cds_subset = cluster_cells(cds_subset)

p1<-plot_cells(cds_subset, color_cells_by="partition",cell_size=.001,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 10))+theme(legend.position = "none")
p2<-plot_cells(cds_subset, color_cells_by="Sex",cell_size=.001,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 10))
p3<-plot_cells(cds_subset, color_cells_by="Diagnosis",cell_size=.001,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 8))
p4<-plot_cells(cds_subset, color_cells_by="Samples",cell_size=.001,label_cell_groups=0,show_trajectory_graph=FALSE)
p6<-plot_cells(cds_subset, color_cells_by="cluster",cell_size=.1,label_cell_groups=0,show_trajectory_graph=FALSE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm") )+theme(legend.position = "none")
#pdf(paste0("/Users/relyanow/Documents/Sage/ROSMAP_sn_writeup/figures/UW_sn_summary_marker_genes.pdf"))
grid.arrange(arrangeGrob(p1,p6, ncol=2),arrangeGrob(p3,p2,ncol=2),p4, heights=c(2,2,4), ncol=1)
#dev.off()

saveRDS(cds_subset, file = "UW_monocle_preprocessed_cds.rds")
cds_uw <- readRDS(file = "UW_monocle_preprocessed_cds.rds")

#plot the broad marker genes
#png(paste0("/Users/relyanow/Documents/Sage/ROSMAP_sn_writeup/figures/UW_sn_diffexp_markers.png"),width = 150, height = 150, units='mm', res = 300)
plot_cells(cds_subset, genes=c("SYT1","SNAP25","GRIN1","GAD1","GAD2","SLC17A7","CAMK2A","NRGN","AQP4",
                               "GFAP","MBP","MOBP","PLP1","PDGFRA","VCAN","CD74","CSF1R","C3","FLT1","CLDN5"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,reduction_method="UMAP",cell_size=.1)
#dev.off()

plot_cells(cds_subset, color_cells_by="partition",cell_size=.001,label_cell_groups=1)

cds_subset$broad.cell.type = cds_subset$Sex
for (partition in unique(partitions(cds_subset))){
  interneuron <- mean(counts(cds_subset)[rownames(cds_subset)=="GAD1",partitions(cds_subset)==partition])
  interneuron <- mean(c(interneuron,mean(counts(cds_subset)[rownames(cds_subset)=="GAD2",partitions(cds_subset)==partition])))
  Ex <- mean(counts(cds_subset)[rownames(cds_subset)=="CAMK2A",partitions(cds_subset)==partition])
  Ex <- mean(c(Ex,mean(counts(cds_subset)[rownames(cds_subset)=="NRGN",partitions(cds_subset)==partition])))
  Ast <- mean(counts(cds_subset)[rownames(cds_subset)=="AQP4",partitions(cds_subset)==partition])
  Ast <- mean(c(Ast,mean(counts(cds_subset)[rownames(cds_subset)=="GFAP",partitions(cds_subset)==partition])))/2
  Oli <- mean(counts(cds_subset)[rownames(cds_subset)=="MBP",partitions(cds_subset)==partition])
  Oli <- mean(c(Oli,mean(counts(cds_subset)[rownames(cds_subset)=="PLP1",partitions(cds_subset)==partition])))/3
  OPC <- mean(counts(cds_subset)[rownames(cds_subset)=="PDGFRA",partitions(cds_subset)==partition])
  OPC <- mean(c(OPC,mean(counts(cds_subset)[rownames(cds_subset)=="VCAN",partitions(cds_subset)==partition])))
  Mic <- mean(counts(cds_subset)[rownames(cds_subset)=="CD74",partitions(cds_subset)==partition])
  Mic <- mean(c(Mic,mean(counts(cds_subset)[rownames(cds_subset)=="CSF1R",partitions(cds_subset)==partition])))
  End <- mean(counts(cds_subset)[rownames(cds_subset)=="FLT1",partitions(cds_subset)==partition])
  End <- mean(c(End,mean(counts(cds_subset)[rownames(cds_subset)=="CLDN5",partitions(cds_subset)==partition])))
  names <- c('In','Ex','Ast','Oli','Opc','Mic','End')
  means <- c(interneuron, Ex, Ast, Oli, OPC, Mic, End)
  best_name <- names[which(means == max(means))]
  cds_subset$broad.cell.type[partitions(cds_subset)==partition] = best_name
}
plot_cells(cds_subset, color_cells_by="broad.cell.type",cell_size=.001,label_cell_groups=0,show_trajectory_graph=FALSE)

for (celltype in unique(cds_subset$broad.cell.type)){
  print(celltype)
  print(length(cds_subset$broad.cell.type[cds_subset$broad.cell.type==celltype])/length(colnames(cds_subset)))
}

#recode diagnosis to control vs ad (pt 6672 was the only control)
cds_subset$dementia = cds_subset$Diagnosis
cds_subset$Diagnosis = cds_subset$dementia
cds_subset$Diagnosis[cds_subset$ids==6672]='Control'
cds_subset$Diagnosis[cds_subset$ids!=6672]='AD'


#####UNTESTED CODE FROM REBECCA--STILL WORKING OUT ISSUES######

## Identify Mic1 subcluster in UW data
#l = c(l,as.vector(mathy_marker_genes$gene.name[mathy_marker_genes$subpopulation=='Mic0']))
l = as.vector(mathy_marker_genes$gene.name[mathy_marker_genes$subpopulation=='Mic1'])
#l = c(l,as.vector(mathy_marker_genes$gene.name[mathy_marker_genes$subpopulation=='Mic2']))
l = unique(l)
length(l)
inds = c()
for (gene in l){
  if ((gene %in% rownames(cds_uw))){
    inds = c(inds,which(rownames(cds_uw)==gene))
  }
}
cds_subset <- cds_uw[inds,cds_uw$broad.cell.type=='Mic']
cds_subset <- preprocess_cds(cds_subset, num_dim = 30,residual_model_formula_str="~PMI")
cds_subset = reduce_dimension(cds_subset)
cds_subset = cluster_cells(cds_subset)

#cs<-cs[,cs$apoe>0]
cs = cds_uw[,cds_uw$broad.cell.type=='Mic']
cs$Mic1 = cs$Sex
cs$Mic1[clusters(cds_subset)==7]='Mic1'
cs$Mic1[clusters(cds_subset)==11]='Mic1'
cs$Mic1[clusters(cds_subset)==5]='Mic1'
cs$Mic1[log(as.vector(counts(cs)[rownames(cs)=='FTL',]))>1.5]='Mic1'
cs$Mic1[cs$Mic1!='Mic1']='Mic0'
cds_subset$Mic1 = cs$Mic1
plot_cells(cds_subset, genes=c("FTL","APOE","TPT1","RPL13","SLC5A11","NLGN1"),cell_size=1,label_cell_groups=0,show_trajectory_graph=FALSE)
p1<-plot_cells(cds_subset, color_cells_by="ePRS",cell_size=1,label_cell_groups=0,show_trajectory_graph=TRUE,label_leaves=FALSE,label_branch_points=FALSE)
p2<-plot_cells(cds_subset, color_cells_by="Diagnosis",cell_size=.1,label_cell_groups=0,show_trajectory_graph=TRUE,label_leaves=FALSE,label_branch_points=FALSE)
p3<-plot_cells(cds_subset, color_cells_by='Mic1',cell_size=1,label_cell_groups=1,show_trajectory_graph=TRUE,label_leaves=FALSE,label_branch_points=FALSE)
p4<-plot_cells(cds_subset, color_cells_by="ids",cell_size=1,label_cell_groups=0,show_trajectory_graph=TRUE,label_leaves=FALSE,label_branch_points=FALSE)
grid.arrange(p1,p2,p3,p4,ncol=2)

cs2 <- cs[,cs$Sex=='M']
gene_fits = fit_models(cs2, model_formula_str = "~Mic1+PMI")
fit_coefs = coefficient_table(gene_fits)
fit_coefs  <- subset(fit_coefs, term != "(Intercept)")
fit_coefs2 <- subset(fit_coefs, term == "Mic1Mic1")
fit_coefs2 <- subset(fit_coefs2, status == 'OK')
fit_coefs2 <- fit_coefs2[order(fit_coefs2$q_value),]
fit_coefs2 <- subset(fit_coefs2, q_value < .05)
fit_coefs2
fit_coefs[which(fit_coefs$gene_short_name=='CST3'),]
keeps = c("gene_short_name","test_val","p_value","normalized_effect","q_value")
fit_coefs2 = fit_coefs2[keeps]
#write.csv(fit_coefs2[order(-fit_coefs2$normalized_effect),],file='/Users/relyanow/Documents/Sage/notebooks/figures/sn_Mathys/mic1.csv')