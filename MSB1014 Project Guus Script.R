##########################################################################
##                                                                      ##
##  Guus Wilmink, i6156211                                              ##
##  MSB1014 Network Biology project                                     ##
##  Identifying key genes, protein complexes, and biological processes  ##
##  shared between the five most common forms of cancer                 ##
##                                                                      ##
##########################################################################


# Setup --------------------------------------------------

# Libraries
if(!"BiocManager" %in% installed.packages()){
  install.packages("BiocManager")
}
library(BiocManager)
if(!require(rstudioapi)){
  BiocManager::install("rstudioapi")
  library(rstudioapi)
}
if(!require(RCy3)){
  BiocManager::install("RCy3")
  library(RCy3)
}
if(!require(igraph)){
  BiocManager::install("igraph")
  library(igraph)
}
if(!require(clusterProfiler)){
  BiocManager::install("clusterProfiler")
  library(clusterProfiler)
}
if(!require(org.Hs.eg.db)){
  BiocManager::install("org.Hs.eg.db")
  library(org.Hs.eg.db)
}
if(!require(STRINGdb)){
  BiocManager::install("STRINGdb")
  library(STRINGdb)
}
if(!require(dplyr)){
  BiocManager::install("dplyr")
  library(dplyr)
}

# Set working directory to current folder
setwd(dirname(getActiveDocumentContext()$path)) 

# Set seed for reproducible results
set.seed(123)

# Check connection to Cytoscape
cytoscapePing ()


# STRING graph of CGC PPIs ------------------------------------------------

# Load the data
COSMIC <- read.csv("COSMIC.csv")
write.table(COSMIC$Gene.Symbol, file='genelist.txt', row.names=F, col.names=F, quote=F)
# Set up STRING query
string.cmd <- paste('string protein query taxonID=9606 cutoff=0.8 query=', 
                    paste(COSMIC$Gene.Symbol, collapse=","),'"', sep="")
# Send network to Cytoscape
commandsGET(string.cmd)
# Analyze the network
analyzeNetwork()


# Hubs and bottlenecks ------------------------------------------

# Plot degree and betweenness distributions
Centralities <- getTableColumns(columns=c("display name", "Degree", "BetweennessCentrality"))
deg_dis <- data.frame(deg=sort(Centralities$Degree, decreasing=T),
                      node=1:nrow(Centralities)
                      )
plot(deg_dis, xlab="Degree", ylab="Node")
title("Degree distribution")

bet_dis <- data.frame(bet=sort(Centralities$BetweennessCentrality, decreasing=T),
                      node=1:nrow(Centralities)
                      )
plot(bet_dis, xlab="Betweenness", ylab="Node")
title("Betweenness distribution")

# Based on the degree distribution, we select degree >= 80 for hub nodes
# For betweenness, we select 0.03
deg_thres <- 80
bet_thres <- 0.03
# Hubs
hubs <- createColumnFilter("Hubs", "Degree", deg_thres, "GREATER_THAN_OR_EQUAL")
# Bottlenecks
bottlenecks <- createColumnFilter("Bottlenecks", "BetweennessCentrality", bet_thres, "GREATER_THAN_OR_EQUAL")
# Create encodings for hubs and bottlenecks
string_mapped <- getTableColumns()
hub_bottle_df <- data.frame("name" = string_mapped$'name',
                            "Hubs" = rep("0", nrow(string_mapped)),
                            "Bottlenecks" = rep("0", nrow(string_mapped)),
                            "Hub_bottle_both" = rep("0", nrow(string_mapped)))
# Binary hub vector
hub_bottle_df$Hubs[hub_bottle_df$'name' %in% hubs$nodes] <- "hub"
# Binary bottleneck vector
hub_bottle_df$Bottlenecks[hub_bottle_df$'name' %in% bottlenecks$nodes] <- "bottleneck"
# Vector with hub = 1, bottleneck = 2, hub and bottleneck simultaneously = 3
hub_bottle_df$"Hub_bottle_both" <- hub_bottle_df$Hubs
hub_bottle_df$"Hub_bottle_both"[hub_bottle_df$'name' %in% bottlenecks$nodes] <- "bottleneck"
hub_bottle_df$"Hub_bottle_both"[hub_bottle_df$'name' %in% hubs$nodes & 
                                hub_bottle_df$'name' %in% bottlenecks$nodes] <- "both"
# Add new vectors to node table
loadTableData(hub_bottle_df, "name")
# Assign hub nodes as red, bottlenecks as blue, and hubs+bottlenecks as green
# Make all other nodes gray
setNodeColorMapping("Hub_bottle_both", table.column.values=c("hub", "bottleneck", "both"), 
                    colors=c("#ff4d50", "#5267ff", "#91ff87"),  mapping.type='d',
                    style="STRING style v1.5")
setNodeColorMapping("Hub_bottle_both", table.column.values=c("hub", "bottleneck", "both"), 
                    colors=c("#aeaeae", "#aeaeae", "#aeaeae"), default.color="#aeaeae",  mapping.type='d',
                    style="STRING style v1.5")
setEdgeColorDefault("#9867c5", style='STRING style v1.5')

# Clustering: largest connected component --------------------------------------------------------

# Largest cluster
networkSuid <- getNetworkSuid()
setCurrentNetwork(network=networkSuid)
renameNetwork("Filtered Network (largest component)", network=as.numeric(networkSuid))
clearSelection(network=networkSuid)
# Filter network
clustermaker <- paste("cluster connectedcomponents createGroups=FALSE network=SUID:",networkSuid, sep="")
res <- commandsGET(clustermaker)
createColumnFilter("Largest component", "__ccCluster", 1, "IS_NOT", network=networkSuid)
applyFilter("Largest component", network=networkSuid)
deleteSelectedNodes(network=networkSuid)

# Clustering: Run MCODE ----------------------------------------------

mcode <- paste("mcode cluster haircut=TRUE degreeCutoff=3 kCore=2 network=SUID:",
               networkSuid, sep=" ")
res_mcode <- commandsGET(mcode)


