# Sources:
#   NCBI GEO2R and examples from Genetic Sequence Analysis
#   https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html

library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE63379", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL17586", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
group_membership <- "1111111111111111111111111111111100000000000000000000000000000000000"
sml <- strsplit(group_membership, split="")[[1]]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("IBS","healthy"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","SPOT_ID"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

# summarize test results by expression
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)

# volcano plot (log P-value vs log fold change)
colnames(fit2) # contrast names
ct <- 1        # contrast of interest
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))
