BiocManager::install(c("GOFunction","org.Hs.eg.db","pathview","ReactomePA"))

library(org.Hs.eg.db)
library(pathview)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("include.R")

up <- read.csv(file="RESULTS-IBD.csv")

#convert the eigth column(symbol) to Entrez Gene
upsymbol2eg <- id2eg(as.character(up[,1]),category="symbol",org="Hs")
print(head(upsymbol2eg))
dim(upsymbol2eg)

# save the second column of upsymbol2eg map into a vector called entrez_up
entrez_up <- upsymbol2eg[,3]
entrez_up <- as.numeric(unique(entrez_up[!is.na(entrez_up)]))

# read all the genes for background
all <- read.csv(file="RESULTS-diseases-compiled.csv")
allsymbol2eg <- id2eg(as.character(all[,1]),category="symbol",org="Hs")
print(head(allsymbol2eg))

# save the second column of upsymbol2eg map into a vector called entrez_all
entrez_all <- allsymbol2eg[,2]
entrez_all <- as.numeric(unique(entrez_all[!is.na(entrez_all)]))

# Perform gene ontology cellular component enrichment for up-regulated genes
sigUpTermCC <- GOFunction(entrez_up,entrez_all,organism="org.Hs.eg.db",
                          ontology="CC",fdrmethod="BY",bmpSize = 4000,
                          filename="PDAC-GO-MF",fdrth=0.1)
# There will be two output files with prefix 'upSigTermCC'

sigTerm <- GOFunction(entrez_up,entrez_all,organism="org.Hs.eg.db",
                      ontology="CC",fdrmethod="BH",filename="PDAC-GO-CC",
                      fdrth=0.1,bmpSize = 4000)
# There will be two output files with prefix 'upSigTermBP'



