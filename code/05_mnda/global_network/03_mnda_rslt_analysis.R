library(data.table)

library(mnda)
library(keras)
# install_keras()
library(aggregation)

# rslt_dir <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/results/5_fold_cv/all/"
rslt_dir     <- "~/bio/code/kul/dex-stim-human-array-isns/results/data/5_fold_cv/all/"
input_graph  <- paste0(rslt_dir, "mnda/global_network_mnda_input_2K.csv")
mnda_rslt_fn <- paste0(rslt_dir, "mnda/global_netqork_mnda_distances_2K.rds")

graph_data <- fread(input_graph, sep = ";")
mnda_rslt  <- readRDS(mnda_rslt_fn)
mnda_dist  <- mnda_rslt$mnda_dist

nodes <- mnda_dist$high_var_nodes

graph_to_plot = cbind(graph_data[,1:2], W = graph_data[,3] - graph_data[,4])
G = mnda::as.igraph(graph_to_plot)

nodes <- mnda_dist$significant_nodes[1:15]
subgraph_plot(G, nodes)


mnda_embedding_space <- mnda_rslt$embedding_space
