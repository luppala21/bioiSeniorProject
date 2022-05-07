BiocManager::install(c("GOFunction","org.Hs.eg.db","pathview","ReactomePA"))

library(org.Hs.eg.db)
library(pathview)
library(ReactomePA)

print(getwd())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

print(list.files())

diffExpress <- read.csv(file="RESULTS-IBS.csv")

#convert gene symbols to Entrez ids
convert_entrez <- id2eg(as.character(diffExpress[,1]),category="symbol",org="Hs")

# save entrez_ids separately for genes being examined
entrez_ids <- convert_entrez[,2]
entrez_ids <- as.numeric(unique(entrez_ids[!is.na(entrez_ids)]))

# read all the genes for background
background <- read.csv(file="RESULTS-diseases-compiled.csv")
background_entrez <- id2eg(as.character(background[,1]),category="symbol",org="Hs")
print(head(background_entrez))

# save entrez_all separately for background genes
entrez_all <- background_entrez[,2]
entrez_all <- as.numeric(unique(entrez_all[!is.na(entrez_all)]))

p.values <- diffExpress$microbiome_pvalue
names(p.values) <-convert_entrez[,2] #assign entrez ids to pvalues
p.values <- p.values[!is.na(names(p.values))] #remove entries with no entrez ids

# all functional annotation is based on homo sapiens pathways

#----KEGG----
# hsa05321 -> IBD
# hsa05212 -> pancreatic malignancies
pv.out <- pathview(gene.data = -log10(p.values), pathway.id = "hsa05321", species = "hsa", out.suffix = "kegg_pathway")

x <- enrichPathway(gene = entrez_ids, pvalueCutoff = 0.2,qvalueCutoff = 0.2, readable=T)

head(data.frame(x))

#----REACTOME----
write.table(data.frame(x), file="reactome.up.csv", quote=F, row.names=F, col.names=T)
barplot(x, showCategory = 10) # only top 10 enriched pathways are shown

#----GOFunction----
sigUpTermCC <- GOFunction(entrez_up,entrez_all,organism="org.Hs.eg.db",
                          ontology="CC",fdrmethod="BY",bmpSize = 4000,
                          filename="PDAC-GO-MF",fdrth=0.1)
