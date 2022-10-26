enrichment <- read.csv('cluster 5 enrichment.csv')

ind <- which(enrichment$source %in% c("WikiPathways"))
enrichment_WP <- enrichment[ind,]

ind <- which(enrichment$source %in% c("Gene Ontology Molecular Function"))
enrichment_GO <- enrichment[ind,]

genecount <- getTableColumns(columns='display name')
genecount$countGO <- c(rep(0, nrow(genecount)))
genecount$countWP <- c(rep(0, nrow(genecount)))

iter <- 0
for (gene in genecount$`display name`) {
  iter <- iter+1
  for (row in 1:nrow(enrichment_GO)) {
    if (gene %in% unlist(strsplit(enrichment_GO$intersecting.genes[row], '|', fixed=T))) {
      genecount$countGO[iter] <- genecount$countGO[iter] + 1
    }
    if (gene %in% unlist(strsplit(enrichment_WP$intersecting.genes[row], '|', fixed=T))) {
      genecount$countWP[iter] <- genecount$countWP[iter] + 1
    }
  }
}
rownames(genecount) <- genecount$`display name`
loadTableData(genecount, table.key.column='display name')  
