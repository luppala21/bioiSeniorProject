BiocManager::install(c("GOFunction","org.Hs.eg.db","pathview","ReactomePA"))

library(org.Hs.eg.db)

library(pathview)

library(ReactomePA)

print(getwd())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

print(list.files())

up <- read.csv(file="RESULTS-IBS.csv")

#convert the eigth column(symbol) to Entrez Gene
upsymbol2eg <- id2eg(as.character(up[,1]),category="symbol",org="Hs")
print(head(upsymbol2eg))
dim(upsymbol2eg)

# save the second column of upsymbol2eg map into a vector called entrez_up
entrez_up <- upsymbol2eg[,2]
entrez_up <- as.numeric(unique(entrez_up[!is.na(entrez_up)]))

# read all the genes for background
all <- read.csv(file="RESULTS-diseases-compiled.csv")
allsymbol2eg <- id2eg(as.character(all[,1]),category="symbol",org="Hs")
print(head(allsymbol2eg))

# save the second column of upsymbol2eg map into a vector called entrez_all
entrez_all <- allsymbol2eg[,2]
entrez_all <- as.numeric(unique(entrez_all[!is.na(entrez_all)]))

p.values <- up$microbiome_pvalue
#assign entrez ids to pvalues
names(p.values) <-upsymbol2eg[,2]
#remove entries with no entrez ids
p.values <- p.values[!is.na(names(p.values))]

pv.out <- pathview(gene.data = -log10(p.values), pathway.id = "hsa05321", species = "hsa", out.suffix = "kegg_pathway")

x <- enrichPathway(gene = entrez_up, pvalueCutoff = 0.2,qvalueCutoff = 0.2, readable=T)

head(data.frame(x))

# save the result into a file **reactome.up.csv**
write.table(data.frame(x), file="reactome.up.csv", quote=F, row.names=F, col.names=T)
# plot the top 10 pathways enriched (Click on svZoom to see everything in the figure)
barplot(x, showCategory = 10)

#Open this file to complete this exercise. Make sure you see two files (1) GO Plot and (2) CSV file is seen in your working directory
source("gofunction-Deliverable4.R")