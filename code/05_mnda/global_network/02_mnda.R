library(data.table)

library(mnda)
library(keras)
# install_keras()
library(aggregation)

library(tictoc)

# rslt_dir   <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/results/5_fold_cv/all/"

# args        <- commandArgs(T)
# input_graph <- as.character(args[1]) 
# out_fn      <- as.character(args[2]) 

rslt_dir    <- "~/bio/code/kul/dex-stim-human-array-isns/results/data/5_fold_cv/all/"
input_graph <- paste0(rslt_dir, "mnda/global_network_mnda_input_2K.csv")
out_fn      <- paste0(rslt_dir, "mnda/global_netqork_mnda_distances_2K.rds")

## Load data

graph_data <- fread(input_graph, sep = ";")

## Create embedding space

print("Start creating embedding space ...", quote = F)
tic("Create embedding space")

embedding_space_list <- mnda_embedding_2layer(graph_data, train.rep = 100, 
                                              walk.rep = 100,
                                              n.steps = 10,
                                              embedding.size = 10, 
                                              batch.size = 32,
                                              epochs = 50,
                                              random.walk = T, null.perm = TRUE)

toc()

## Calculate the node-pair distances and p-values

print("Start calculation of the node-pair distances and p-values ...", quote = F)
tic("Calculate the node-pair distances and p-values")

mnda_dist <- mnda_node_detection_2layer(embedding_space_list, p.adjust.method = "BH")

toc()

## Save results

tic("Save results")
saveRDS(list(embedding_space = embedding_space_list, mnda_dist = mnda_dist), out_fn)
toc()