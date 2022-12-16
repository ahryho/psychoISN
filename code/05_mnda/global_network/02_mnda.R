library(data.table)

library(mnda)
library(keras)
# install_keras()
library(aggregation)

# rslt_dir   <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/results/5_fold_cv/all/"
# rslt_dir    <- "~/bio/code/kul/dex-stim-human-array-isns/results/data/5_fold_cv/all/"
# input_graph <- paste0(rslt_dir, "mnda/global_network_mnda_input.csv")
# out_fn      <- paste0(rslt_dir, "mnda/global_netqork_mnda_distances.rds")

args        <- commandArgs(T)
input_graph <- as.character(args[1]) 
out_fn      <- as.character(args[2]) 

## Load data

graph_data <- fread(input_graph, sep = ";")

## Create embedding space

embedding_space_list <- mnda_embedding_2layer(graph_data, train.rep = 50, walk.rep = 100,
                                              embedding.size = 10, random.walk = T, null.perm = TRUE)

## Calculate the node-pair distances and p-values

mnda_dist <- mnda_node_detection_2layer(embedding_space_list, p.adjust.method = "none")

saveRDS(mnda_dist, out_fn)