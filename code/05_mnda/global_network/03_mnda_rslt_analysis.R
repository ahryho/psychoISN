library(data.table)

library(mnda)
library(keras)
# install_keras()
library(aggregation)

# rslt_dir <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/results/5_fold_cv/all/"
rslt_dir     <- "~/bio/code/kul/dex-stim-human-array-isns/results/data/5_fold_cv/all/"
input_graph  <- paste0(rslt_dir, "mnda/global_network_mnda_input.csv")
mnda_rslt_fn <- paste0(rslt_dir, "mnda/global_netqork_mnda_distances.rds")

mnda_dist <- fread(mnda_rslt_fn)

nodes <- mnda_dist$high_ranked_nodes

graph_to_plot = cbind(graph_data[300:350,1:2], W = graph_data[300:350,3] - graph_data[300:350,4])
G = mnda::as.igraph(graph_to_plot)

subgraph_plot(G, nodes)
