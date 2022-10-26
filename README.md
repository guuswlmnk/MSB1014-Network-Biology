# MSB1014-Network-Biology
Scripts and cytoscape session MSB1014 project

Required data for R scripts (gene list) can be downloaded from: https://cancer.sanger.ac.uk/census

The Cytoscape session can be readily used without downloading the data, or alternatively the network and clusters can be generated using the gene list and R scripts.

R script instructions:
1. Name the downloaded csv file with gene names "COSMIC.csv" and place in same directory as scripts
2. Open Cytoscape (3.9.1).
3. Run "MSB1014 Project Guus Script.R"
4. In Cytoscape, create networks from each MCODE cluster.
5. For each of the cluster networks, run EnrichmentTable app from Cytoscape. Set display name as gene ID column and species as homo sapiens. Set correction method to FDR with corrected p-value < 0.05.
6. Save the table as "cluster x enrichment.csv" for each network, substituting x with the cluster number.
7. Run the "enrichment styling.R" script. Make sure you have the network you want to visualize enrichments for opened in Cytoscape. In line 1, specified the correct enrichment table csv file "cluster x enrichment.csv". This adds node columns in Cytoscape with the frequencies of genes occurring in GO molecular function terms or WikiPathway terms that can be used for visualization.
