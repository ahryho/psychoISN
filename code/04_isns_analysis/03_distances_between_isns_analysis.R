rm(list = ls())

library(dplyr)
library(data.table)
library(rdist)
library(gdsfmt)

library(ggplot2)
library(dendextend)
library(factoextra)

# source("~/kul/dex-stim-human-array-isns/code/00_functions/load_gds_files.R")
# dir_pre <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/"
# dnam_pca_dir <- paste0("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/output/data/methylation/02_dmp/")

source("~/bio/code/kul/dex-stim-human-array-isns/code/00_functions/load_gds_files.R")
dir_pre <- "~/bio/code/kul/dex-stim-human-array-isns/"

rslt_dir     <- paste0(dir_pre, "results/5_fold_cv/all/")
treatment    <- "veh"

##

pheno           <- fread(paste0(dir_pre, "input/pheno/pheno_full_for_kimono.csv"), dec = ",")[Include == T][Dex == ifelse(treatment == "veh", 0, 1)]

dist_veh_df     <- fread(paste0(rslt_dir, treatment, "/smccnet_isns_euclidean_dist_mtrx.csv"))
dist_mtrx       <- as.dist(dist_veh_df)
hclust_avg      <- hclust(dist_mtrx, method = "ward.D")

plot(hclust_avg)
nr_cluster      <- 4
# rect.hclust(hclust_avg , k = nr_cluster, border = 2:6)

cut_avg         <- cutree(hclust_avg, k = nr_cluster)

avg_dend_obj <- as.dendrogram(hclust_avg)
avg_col_dend <- color_branches(avg_dend_obj, k = nr_cluster)
plot(avg_col_dend)

dist_veh_df      <- mutate(dist_veh_df, DNA_ID = colnames(dist_veh_df), cluster = as.factor(cut_avg))
count(dist_veh_df,cluster)

dist_veh_df <- dist_veh_df[match(pheno$DNA_ID, dist_veh_df$DNA_ID), ]

fviz_cluster(list(data = dist_veh_df[, 1:196], 
                  cluster = dist_veh_df$cluster),
             geom = "point") + 
  scale_shape_manual(values = as.factor(pheno$Status)) +
  theme_minimal()



## 

# dnam_pca_df <- readRDS(paste0(dnam_pca_dir, "dex_methyl_pca.rds"))
# nr_pcs   <- ncol(dnam_pca_df$rotation)
# 
# dnam_pcs <- data.frame(dnam_pca_df$rotation)
# dnam_pcs <- dnam_pcs[match(pheno$DNAm_ID, rownames(dnam_pcs)), ]
# 
# colnames(dnam_pcs) <- paste("DNAm", colnames(dnam_pcs), sep = "_")
# 
# plt_df <- data.frame(scale(dnam_pcs), DNA_ID = pheno$DNA_ID, 
#                      Treatment = as.factor(pheno$Dex), Sex = as.factor(pheno$Sex), MDD = as.factor(pheno$Status),
#                      pheno[, .(Age, BMI_D1)])
# plt_df <- inner_join(plt_df, dist_df[, .(DNA_ID, cluster)])
# 
# ggplot(data = plt_df, aes(DNAm_PC4, DNAm_PC6, col = cluster)) + #, label = DNA_ID)) +
#   geom_vline(xintercept = 0, linetype = "dashed", size = 0.4) +
#   geom_hline(yintercept = 0, linetype = "dashed", size = 0.4) +
#   geom_point(size = 2, alpha = 1) +
# #  geom_text() +
# #  xlab(paste0("DNAm PC1 (", signif(propvar[1], 2), "%)")) +
# #  ylab(paste0("DNAm PC2 (", signif(propvar[2], 2), "%)")) +
#   theme( panel.background = element_blank(),
#          plot.title = element_text(size = 10),
#          axis.title = element_text(size = 10),
#          axis.text.x = element_text(angle = 0, hjust = 1),
#          legend.position = "none")
# 
## PCA on multiomics data

dnam_gds_fn <- paste0(dir_pre, "input/dnam/mad_filtered/gds/methyl_beta_mtrx_corrected_for_cov_mad80_filtered_", treatment, ".gds")
snps_gds_fn <- paste0(dir_pre, "input/snps/ld_pruned/gds/dex_geno_imputed_maf_ld_pruned_from_gen.gds")

global_net  <- readRDS(paste0(rslt_dir, treatment, "/networks/smccnet_global_chr_all.rds"))

dnam_mtrx <- LoadMethyl(dnam_gds_fn, is_mad = F)
snps_mtrx <- LoadGenotype(snps_gds_fn, is_ld = F)

dnam_mtrx <- dnam_mtrx[, colnames(dnam_mtrx) %in% colnames(global_net)]
snps_mtrx <- snps_mtrx[, colnames(snps_mtrx) %in% colnames(global_net)]

pca_mtrx <- cbind(dnam_mtrx, snps_mtrx)

pca_res <- prcomp(pca_mtrx, scale = T, center = T)

fviz_eig(pca_res)

fviz_pca_var(pca_res,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

pca_eigenvec <- as.data.frame(pca_res$x)
# pca_eigenvec <- mutate(pca_eigenvec, DNA_ID = rownames(pca_eigenvec))
pca_eigenvec <- pca_eigenvec[match(pheno$DNA_ID, rownames(pca_eigenvec)), ]

plt_df <- data.frame((pca_eigenvec), DNA_ID = pheno$DNA_ID, 
                     Treatment = as.factor(pheno$Dex), Sex = as.factor(pheno$Sex), MDD = as.factor(pheno$Status),
                     pheno[, .(Age, mAge_Hovath, BMI_D1)])
plt_df <- inner_join(plt_df, dist_df[, .(DNA_ID, cluster)])

fviz_pca_ind(pca_res,axes = c(1, 2), 
             label="var", col.var = "black", #setas
             geom = "point", pointsize = 2, col.ind=plt_df$cluster, 
             addEllipses = TRUE, ellipse.type ="confidence", palette = "aaas",
             legend.title = "Cluster") + theme_minimal()

ggplot(data = plt_df, aes(PC1, PC2, col = cluster)) + #, label = DNA_ID)) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.4) +
  geom_point(size = 2, alpha = 1) +
  #  geom_text() +
  #  xlab(paste0("DNAm PC1 (", signif(propvar[1], 2), "%)")) +
  #  ylab(paste0("DNAm PC2 (", signif(propvar[2], 2), "%)")) +
  theme( panel.background = element_blank(),
         plot.title = element_text(size = 10),
         axis.title = element_text(size = 10),
         axis.text.x = element_text(angle = 0, hjust = 1),
         legend.position = "none") +
  scale_color_brewer(palette = "Set2")


## Average DNAm beta value per sample

dnam_avg_beta_df <- data.frame(DNA_ID = rownames(dnam_mtrx), avg_dnam_beta = apply(dnam_mtrx, 1, median))
plt_df           <- left_join(plt_df, dnam_avg_beta_df)

plt.size <- 10
ggplot(plt_df, aes(x = cluster, y = BMI_D1)) +
  geom_violin(aes(color = cluster, fill = cluster)) +
  geom_boxplot(width = 0.05, fill = "white", color = "black") +
  theme(legend.position = "none", # c(0.9, 0.9), 
        legend.title = element_blank(), 
        legend.text = element_text(size = plt.size),
        panel.grid.major = element_blank(), 
        panel.background = element_blank(),
        plot.title = element_text(size = plt.size), 
        axis.title = element_text(size = plt.size),
        axis.text = element_text(size = plt.size, colour = "black")) + 
  labs(title = "",
       y = "") +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2")
