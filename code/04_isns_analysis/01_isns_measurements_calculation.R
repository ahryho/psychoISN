library(dplyr)
library(tidyverse)
library(data.table)
library(gdsfmt)

library(igraph)

library(tidyverse)

source("~/kul/dex-stim-human-array-isns/code/00_functions/load_gds_files.R")
source("~/kul/dex-stim-human-array-isns/code/00_functions/save_isns_as_single_object.R")
source("~/kul/dex-stim-human-array-isns/code/00_functions/get_chord_diagramm.R")

rslt_dir <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/results/5_fold_cv/all/"
setwd(rslt_dir)

# 1. Load global network

treatment <- "dex"
dex_global_net  <- readRDS(paste0(treatment, "/networks/smccnet_global_chr_all.rds"))
dex_isns_fn_lst <- list.files(paste0(treatment, "/networks"), pattern = "individual", full.names = T)

dex_isns_lst <- save_isns_as_single_object(dex_isns_fn_lst[1:2], treatment, 
                                           ofile = paste0(rslt_dir, treatment, "/networks/smccnet_", treatment, "_isns_chr_all.rds"))


treatment <- "veh"
veh_global_net <- readRDS(paste0(treatment, "/networks/smccnet_global_chr_all.rds"))

# 2. Get the overlap between dex and veh features

dex_features <- rownames(dex_global_net) # 5'388
veh_features <- rownames(veh_global_net) # 4'519

olap_feat <- intersect(dex_features, veh_features) # 758
olap_cpgs <- olap[grepl("cg", olap, ignore.case = T)] # 566
olap_snps <- olap[grepl("rs", olap, ignore.case = T)] # 192

# 3. Calculate net measurements

## Vertex Degree — give the number of edges for each vertex
## Betweenness Centrality — fraction of shortest paths that pass through the vertex
## Vertex Correlation Similarity - correlation similarity between actors

## Create network
net <- veh_global_net[olap_feat, olap_feat]
summary(as.numeric(net))
net_trsh <- 0.02

net_long <- net %>% 
  as.matrix() %>% as.data.frame() %>%
  rownames_to_column %>%
  gather(key = 'key', value = 'value', -rowname) %>% 
  subset(value > net_trsh) %>% 
  mutate(value = signif(value, 3))


psycho_graph <- graph_from_data_frame(net_long, directed = F)
# plot(psycho_graph, edge.arrow.size = 0)


## Vertex degree

vertex_degree <- degree(psycho_graph)
network_distr <- degree_distribution(psycho_graph)

# Betweenness Centrality
betweenness_vertex <- betweenness(psycho_graph, v = V(psycho_graph), directed = F)
