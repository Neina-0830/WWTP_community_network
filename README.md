# WWTP_community_network
Code for "Keystone taxa mediate the trade-off between microbial community stability and performance in activated sludges"

## System Requirements

R (code tested on v 4.0.3)

packages: RMThreshold, igraph, psych, tnet, dplyr, reshape2, ggplot2, statnet, circlize, ggalluvial, ggExtra, bnlearn, simba


## Installation

R https://www.r-project.org/

rstudio https://rstudio.com/

(time required to install software and packages <1hr)

## Instruction to Use

data/

contains the data (in .csv format) used for the analysis.


code/

contains code to replicate results. 

RMT_graph-logAbun.R -> construct the cooccurrence network and calculate basic topology parameters.

Within-Among-module connectivity.R -> calculate within-module connectivities and among-module connectivities of the network.

Power-law-distribution.R -> fit the degreee distribution of AS network and random network.

network_natural_connectivity.R -> network robustness test.

network_node_centrality.R -> calculate node centrality.

get_node_summary.R -> get the node information of AS network.

get_edge_summary.R -> get the edge information of AS network.

Compare_keystone_others.R -> compare the keystone and others taxa in AS network.

Compare_keystone_samples.R -> compare the keystone-positive samples and keystone-negative samples.

Compare_keystone_hubs.R -> compare the keystone and hubs taxa in AS network.

network_robustness_different_targeted_deletion-compare.R -> compare impact of different targeted taxa on AS network stability.

keystone_presence_analysis.R -> select environmental factors favoring the presence of keystone taxa.

consumer_resource_model_analysis.R -> analysis code of consumer-resource model.# WWTP_community_prediction
