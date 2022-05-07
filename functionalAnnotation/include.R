require("DBI")
createGODAG <-function(sigNodes, ontology = "BP") {
  
  nodeLabel <- new.env(hash = T, parent = emptyenv())
  
  isNodeInDAG <- function(node) {
    return(exists(node, envir = nodeLabel, mode = 'logical', inherits = FALSE))
  }
  setNodeInDAG <- function(node) {
    assign(node, TRUE, envir = nodeLabel)
  }
  
  GOParents <- get(paste('GO', ontology, 'PARENTS', sep = ''))
  
  ROOT <- as.character(revmap(GOParents)$all)
  
  adjList <- as.list(GOParents)  
  edgeEnv <- new.env(hash = T, parent = emptyenv())  
  envAddEdge <- function(u, v, type) {
    assign(v, switch(type, is_a = 0, part_of = 1, -1), envir = get(u, envir = edgeEnv))
  }
  
  createNodesRelationship <- function(node) {
    if(isNodeInDAG(node))
      return(1)
    
    ## we put the node in the graph and we get his parents
    setNodeInDAG(node)    # we visit the node
    assign(node, new.env(hash = T, parent = emptyenv()), envir = edgeEnv) # adj list
    
    if(node == ROOT) 
      return(2)
    
    adjNodes <- adjList[[node]]
    
    if(length(adjNodes) == 0)
      cat('\n There are no adj nodes for node: ', node, '\n')
    
    for(i in 1:length(adjNodes)) {
      x <- as.character(adjNodes[i])
      envAddEdge(node, x, names(adjNodes[i]))
      createNodesRelationship(x)
    }
    
    return(0)
  }
  
  ## we start from the most specific nodes
  lapply(sigNodes, createNodesRelationship)
  
  .graphNodes <- ls(edgeEnv)
  .edgeList <- eapply(edgeEnv,
                      function(adjEnv) {
                        aux <- as.list(adjEnv)
                        return(list(edges = match(names(aux), .graphNodes),
                                    weights = as.numeric(aux)))
                      })
  
  ## now we can build the graphNEL object
  GOgraph.topo <- new('graphNEL',
                      nodes = .graphNodes,
                      edgeL = .edgeList,
                      edgemode = 'directed')
  
  require("SparseM") || stop("package SparseM is required")
  
  GOgraph.topo <- sparseM2Graph(t(graph2SparseM(GOgraph.topo, TRUE)),
                                .graphNodes, edgemode = "directed")
  return(GOgraph.topo)
}

localRedundancy <-function(sigTerm, generalAnn, sigTermRelation, annRef, annInterest, ppth, pcth)
  
{
  annRef <- unique(annRef[,1])
  allRefnum <- length(annRef)
  annInterest <- unique(annInterest[,1])
  allInterestnum <- length(annInterest)
  
  sigTermRelationRe <- sigTermRelation
  sigLabel <- array(0,dim=c(nrow(sigTerm),1))
  sigTerm$Label <- sigLabel
  sigTerm$SeLabel <- sigLabel
  La <- sigTerm[,1]
  noRelationTerm <- setdiff(La,union(sigTermRelation[,1],sigTermRelation[,2]))
  sigTerm[sigTerm[,1] %in% noRelationTerm,7] <- 1
  La <- setdiff(La,noRelationTerm)
  
  while (length(La) > 0) {
    
    sigTermRelationRe <- sigTermRelationRe[sigTermRelationRe[,2] %in% La,]
    if (nrow(sigTermRelationRe) == 0)
      leafnode <- La
    else
      leafnode <- setdiff(sigTermRelationRe[,2],sigTermRelationRe[,1])
    
    La <- setdiff(La,leafnode)
    
    
    for (j in c(1:length(leafnode))){
      node <- leafnode[j];
      genes <- generalAnn[generalAnn[,2]==node,1]
      genes <- intersect(genes,annRef)
      sgenes <- intersect(genes,annInterest)
      childnode <- sigTermRelation[sigTermRelation[,1]==node,2]
      #cat("Panode:",node,"\n")
      #cat("Chnode:",childnode,"\n")
      if (length(childnode)==0) {
        sigTerm[sigTerm[,1]==node,7] <- 1
      }
      else {         
        activeChild <- sigTerm[(sigTerm[,1] %in% childnode) & sigTerm[,7]==1,1]
        if (length(activeChild)==0) {
          sigTerm[sigTerm[,1]==node,7] <- 1
        }
        else {
          allcgenes <- generalAnn[generalAnn[,2]==activeChild[1],1]
          if (length(activeChild)>1) {
            for (k in c(2:length(activeChild))) {
              cgenes <- generalAnn[generalAnn[,2]==activeChild[k],1]
              allcgenes <- union(allcgenes,cgenes)
            }
          }
          allcgenes <- intersect(allcgenes,annRef)
          allcsiggenes <- intersect(allcgenes,annInterest)
          
          if(length(setdiff(allcgenes,genes))>0){
            # cat("Parent Node:",node," annotated with ",length(genes)," genes while child Nodes:",as.character(activeChild), "annotated with ",length(allcgenes), " genes.\n")
            #cat("The number of different genes between parent node and child nodes is", length(setdiff(allcgenes,genes)),".\n")
            #cat("This may be the problem from annoation package GOFunction used!\n")
          }else{
            
            extragenes <- setdiff(genes,allcgenes)
            extrasiggenes <- intersect(extragenes,annInterest)
            if (length(extragenes)!=0) {
              fp <- length(extrasiggenes)/length(extragenes)
              fc <- length(allcsiggenes)/length(allcgenes)
              #cat("FP:",fp," FC:",fc," length(extrasiggenes):",length(extrasiggenes)," length(extragenes):",length(extragenes)," length(allcsiggenes):", length(allcsiggenes), " length(sgenes):", length(sgenes), " length(genes):",length(genes)," length(allcgenes):",length(allcgenes),"\n\n\n\n")
              
              p <- 1-phyper(length(extrasiggenes)-1,allInterestnum,allRefnum-allInterestnum,length(extragenes),lower.tail = TRUE,log.p= FALSE)
              pc <- 1-phyper(length(allcsiggenes)-1,length(sgenes),length(genes)-length(sgenes),length(allcgenes),lower.tail = TRUE,log.p= FALSE)
              if (fp>=fc | p<=ppth) {
                sigTerm[sigTerm[,1]==node,7] <- 1
                if (pc>pcth)
                  sigTerm[sigTerm[,1] %in% activeChild,8] <- sigTerm[sigTerm[,1] %in% activeChild,8]+1
              }
            }
          }
        }
      }
    }
  }
  sigTerm[,7] <- sigTerm[,7]+sigTerm[,8]
  sigTerm[sigTerm[,7]>1,7] <- 0;
  sigTermRedun <- sigTerm[sigTerm[,7]==1,c(1:6)]
  return(sigTermRedun)
}

enrichmentFunction <-function(annRef, annInterest, method, fdrth)
  
{
  allRefnum <- length(unique(annRef[,1]))
  allInterestnum <- length(unique(annInterest[,1]))
  
  allAnnterm <- unique(annRef[,2])
  allAnntermL <- length(allAnnterm)
  
  
  refCount <- tapply(annRef[,1],annRef[,2],length)
  refTerm <- levels(factor(annRef[,2]))
  refTermCount <- data.frame(goid=refTerm,refnum=array(refCount,length(refCount)))
  interestCount <- tapply(annInterest[,1],annInterest[,2],length)
  interestTerm <- levels(factor(annInterest[,2]))
  interestTermCount <- data.frame(goid=interestTerm,interestnum=array(interestCount,length(interestCount)))
  
  ref_interest_TermCount <- refTermCount;
  ref_interest_TermCount$interestnum = array(0,dim=c(length(ref_interest_TermCount$goid),1))
  ref_interest_TermCount[ref_interest_TermCount$goid %in% interestTermCount[,1],3]=interestTermCount$interestnum;
  
  
  n <- nrow(ref_interest_TermCount)
  pv <- array(0,dim=c(n,1))
  for (i in c(1:n)){
    p <- 1-phyper(ref_interest_TermCount[i,3]-1,allInterestnum,allRefnum-allInterestnum,ref_interest_TermCount[i,2],lower.tail = TRUE,log.p= FALSE)
    pv[i,1] <- p
  }
  ref_interest_TermCount$pvalue <- pv
  
  ref_interest_TermCount <- ref_interest_TermCount[order(ref_interest_TermCount[,4]),]
  
  if (length(grep(method,"bonferroni"))){
    adjustp <- p.adjust(ref_interest_TermCount$pvalue,method ="bonferroni");
    ref_interest_TermCount$adjustp <- adjustp
    sigTerm <- ref_interest_TermCount[ref_interest_TermCount[,5]<=fdrth,];
  }
  else {
    if (length(grep(method,"BH"))){
      adjustp <- p.adjust(ref_interest_TermCount$pvalue,method ="BH");
      ref_interest_TermCount$adjustp <- adjustp
      sigTerm <- ref_interest_TermCount[ref_interest_TermCount[,5]<=fdrth,];
    }
    else {
      if (length(grep(method,"BY"))){
        adjustp <- p.adjust(ref_interest_TermCount$pvalue,method ="BY");
        ref_interest_TermCount$adjustp <- adjustp
        sigTerm <- ref_interest_TermCount[ref_interest_TermCount[,5]<=fdrth,];
      }
    }
  }
  rownames(sigTerm) <- NULL
  return(list(sigTerm=sigTerm,allTerm=ref_interest_TermCount))       
}

globalRedundancy <-function(generalAnn, sigTermRelation, annRef, annInterest, sigTermRedun, poth, peth)
{
  annRef <- unique(annRef[,1])
  allRefnum <- length(annRef)
  annInterest <- unique(annInterest[,1])
  allInterestnum <- length(annInterest)
  sigTermRedun$overlap = array(0,dim=c(nrow(sigTermRedun),1));
  sigTermenv <- new.env(hash=T,parent=emptyenv())
  assign("sigTerm",sigTermRedun,envir=sigTermenv)
  
  calculateEachTerm <- function (term1) {
    sigTermRedun <- get("sigTerm",sigTermenv)
    gene1 <- generalAnn[generalAnn[,2]==term1,1]
    gene1 <- intersect(gene1, annRef)
    siggene1 <- intersect(gene1, annInterest)
    extrterm <- setdiff(sigTermRedun[,1], term1);
    calculateExtraTerm <- function(term2) {
      gene2 <- generalAnn[generalAnn[,2]==term2,1];
      gene2 <- intersect(gene2, annRef);
      siggene2 <- intersect(gene2, annInterest)
      po <- sigTermRelation[(sigTermRelation[,1]==term1 & sigTermRelation[,2]==term2) | (sigTermRelation[,1]==term2 & sigTermRelation[,2]==term1),]
      if (nrow(po)==0){
        refov <- intersect(gene1,gene2);
        if (length(refov)>0) {
          sigov <- intersect(siggene1,siggene2)
          extra1 <- setdiff(gene1,refov)
          extrasig1 <- intersect(extra1, annInterest)
          extra2 <- setdiff(gene2,refov)
          extrasig2 <- intersect(extra2, annInterest)
          if(length(extra2)==0){
            return(0)
          }
          else{
            pex2 <- 1-phyper(length(extrasig2)-1,allInterestnum,allRefnum-allInterestnum,length(extra2),lower.tail = TRUE,log.p= FALSE)
            po2 <- 1-phyper(length(sigov)-1,length(siggene2),length(gene2)-length(siggene2),length(refov),lower.tail = TRUE,log.p= FALSE)
            if(length(extra1)==0){
              if ((po2>poth) | (po2<=poth & pex2<=peth)){
                sigTermRedun[sigTermRedun[,1]==term1,7] <- 1
                assign("sigTerm",sigTermRedun,envir=sigTermenv)
              }
            }
            else{            
              pex1 <- 1-phyper(length(extrasig1)-1,allInterestnum,allRefnum-allInterestnum,length(extra1),lower.tail = TRUE,log.p= FALSE)
              po1 <- 1-phyper(length(sigov)-1,length(siggene1),length(gene1)-length(siggene1),length(refov),lower.tail = TRUE,log.p= FALSE)
              if((po1<=poth) & (pex1>peth)){
                if ((po2>poth) | (po2<=poth & pex2<=peth)){
                  sigTermRedun[sigTermRedun[,1]==term1,7] <- 1
                  assign("sigTerm",sigTermRedun,envir=sigTermenv)
                }
              }
            }
          }
        }
      }
    }
    lapply(extrterm,calculateExtraTerm)
  }
  lapply(sigTermRedun[,1],calculateEachTerm)
  sigTermRedun <- get("sigTerm",sigTermenv)
  sigTermRedun <- sigTermRedun[sigTermRedun[,7]==0,c(1:6)]
  return(sigTermRedun);
}

GOFunction <-function(interestGenes, refGenes, organism = "org.Hs.eg.db", ontology = "BP", fdrmethod = "BY", fdrth = 0.05, ppth = 0.05, pcth = 0.05, poth = 0.05, peth = 0.05, bmpSize = 2000, filename = "sigTerm")
{
  
  ## Loading the GO annotation package
  require(organism, character.only = TRUE) || stop(paste("package", organism, "is required", sep=" "))
  
  ## Extracting the gene annotation data
  .sql <-  paste("select distinct t1.gene_id,t2.go_id",
                 " from genes as t1 inner join",paste("go", tolower(ontology), "all", sep = "_"),
                 " as t2 on t1._id=t2._id",seq="")
  organism <- strsplit(organism,".db")
  organism <- organism[[1]]
  conn <- get(paste(organism, "_dbconn", sep = ""))()
  generalAnn <- dbGetQuery(conn, .sql)
  annRef <- generalAnn[generalAnn[,1] %in% refGenes,]
  annInterest <- generalAnn[generalAnn[,1] %in% interestGenes,]
  
  ## Calculating the statistically significant terms for interestGenes
  cat("Finding statistically significant terms...\n") 
  termInfo <- enrichmentFunction(annRef, annInterest,fdrmethod,fdrth)
  sigTerm <- termInfo$sigTerm
  if(nrow(sigTerm)==0){
    warning("There is no significant term! \n")
    return(NULL)
  }
  allTerm <- termInfo$allTerm
  
  ## Loading the GO structure data
  require("GO.db") || stop("package GO.db is required")
  conn <- get("GO_dbconn")()
  
  ## Finding the GO term name for statistically siginificant terms
  .sql <- paste("select distinct go_id goid,term name from go_term where ontology='",
                toupper(ontology), "'", sep="")
  allTermName <-  dbGetQuery(conn,.sql)
  sigTermName <- allTermName[allTermName[,1] %in% sigTerm[,1],]
  sigTermName <- sigTermName[order(sigTermName[,1]),]
  sigTerm <- sigTerm[order(sigTerm[,1]),]
  sigTerm$name <- sigTermName[,2]
  sigTerm <- sigTerm[,c(1,6,2,3,4,5)]
  
  
  
  ## Finding the relationship between statistically significant terms
  .sql <- paste("select distinct t1.go_id parentid,t2.go_id childid from ", paste("go", tolower(ontology),   "offspring", sep = "_"), " as t3 inner join  go_term as t1 on t1._id=t3._id inner join go_term as t2", " on t2._id=t3._offspring_id", sep="")
  allTermRelation <-  dbGetQuery(conn,.sql)
  sigTermRelation <- allTermRelation[(allTermRelation[,1] %in% sigTerm[,1]) & (allTermRelation[,2] %in% sigTerm[,1]),]
  rm(allTermRelation)
  
  ## Finding the offspring terms for each statistically significant terms
  
  
  
  ## Local redundance
  cat("Treating for local redundant terms...\n")
  sigTerm_LocalRedun <- localRedundancy(sigTerm, generalAnn, sigTermRelation, annRef, annInterest, ppth, pcth)
  
  ##Global redundance
  
  cat("Treating for global redundant terms...\n")
  sigTerm_GlobalRedun <- globalRedundancy(generalAnn, sigTermRelation, annRef, annInterest, sigTerm_LocalRedun, poth,  peth)
  
  
  ## Display the GO DAG plot for the significant terms
  
  cat("Visualizing the GO DAG...\n")
  require("graph") || stop("package graph is required")
  sigDAG <- createGODAG(as.character(sigTerm[,1]),ontology)
  allDAGTerm <- allTerm[allTerm[,1] %in% nodes(sigDAG),]
  
  dagTermName <- allTermName[allTermName[,1] %in% allDAGTerm[,1],]
  dagTermName <- dagTermName[order(dagTermName[,1]),]
  allDAGTerm <- allDAGTerm[order(allDAGTerm[,1]),]
  allDAGTerm$name <- dagTermName[,2]
  allDAGTerm <- allDAGTerm[,c(1,6,2,3,4,5)]
  sigTermID <- as.character(sigTerm[,1])
  sigTerm_LocalRedunID <- as.character(sigTerm_LocalRedun[,1])
  sigTerm_GlobalRedunID <- as.character(sigTerm_GlobalRedun[,1])
  
  showSigNodes(sigDAG, sigTermID, sigTerm_LocalRedunID, sigTerm_GlobalRedunID, allDAGTerm, bmpSize, filename)
  
  ## Saving the table of statistically significant terms
  label <- array("",dim=c(nrow(sigTerm),1))
  rmsigLocalTerm <- setdiff(sigTermID, sigTerm_LocalRedunID)
  label[sigTermID %in% rmsigLocalTerm,1] <- "Local"
  rmsigGlobalTerm <- setdiff(sigTerm_LocalRedunID, sigTerm_GlobalRedunID)
  label[sigTermID %in% rmsigGlobalTerm,1] <- "Global"
  label[sigTermID %in% sigTerm_GlobalRedunID,1] <- "Final"
  
  sigTerm$FinalResult <- label
  tablename <- paste(filename,".csv",sep="")
  write.csv(sigTerm, tablename, row.names=F)
  cat("\n********************  Results  ********************\n")
  cat("The number of annotated interesting genes:", length(unique(annInterest[,1])), "\n")
  cat("The number of annotated reference genes:", length(unique(annRef[,1])), "\n");
  cat("The number of statistically significant terms:", nrow(sigTerm), "\n");
  cat("The number of terms after treating local redundancy:", nrow(sigTerm_LocalRedun), "\n")
  cat("The number of terms after treating global redundancy:", nrow(sigTerm_GlobalRedun), "\n")
  note = paste("Please see details about the significant terms in the files ",filename,".csv and ",filename,".bmp!\n",sep="")
  cat(note)
  return(sigTerm_GlobalRedun)
}

showSigNodes <-function(DAG, sigTerm, sigTerm_Local, sigTerm_Global, dagTermInfo, bmpSize, filename) {
  
  require('Rgraphviz') || stop('package Rgraphviz is required')
  
  graphAttrs <- getDefaultAttrs(layoutType = 'dot')
  graphAttrs$cluster <- NULL
  graphAttrs$node$shape <- 'ellipse'
  graphAttrs$node$fontsize <- '20'
  
  nodeAttrs <- list()
  edgeAttrs <- list()
  
  allTerm <- as.character(dagTermInfo[,1])
  nodeAttrs$label[allTerm] <- allTerm
  
  
  rmLocalTerm <- setdiff(sigTerm, sigTerm_Local)
  nodeAttrs$color[rmLocalTerm] <- rep('red', length(rmLocalTerm))
  nodeAttrs$shape[rmLocalTerm] <- rep('circle', length(rmLocalTerm))
  
  rmGlobalTerm <- setdiff(sigTerm_Local, sigTerm_Global)
  nodeAttrs$color[rmGlobalTerm] <- rep('red', length(rmGlobalTerm))
  nodeAttrs$shape[rmGlobalTerm] <- rep('box', length(rmGlobalTerm))
  nodeAttrs$height[rmGlobalTerm] <- rep('0.7', length(rmGlobalTerm))
  nodeAttrs$width[rmGlobalTerm] <- rep('0.7', length(rmGlobalTerm))
  
  nodeAttrs$color[sigTerm_Global] <- rep('red', length(sigTerm_Global))
  nodeAttrs$shape[sigTerm_Global] <- rep('rectangle', length(sigTerm_Global))
  nodeAttrs$height[sigTerm_Global] <- rep('0.7', length(sigTerm_Global))
  nodeAttrs$width[sigTerm_Global] <- rep('1.1', length(sigTerm_Global))
  
  
  
  dagTermInfo[dagTermInfo[,5]<2.2E-16,5] <- 2.2E-16;
  dagTermInfo[dagTermInfo[,6]<2.2E-16,6] <- 2.2E-16;
  dagTermInfo$colorran <- round(log10(dagTermInfo[,6])-range(log10(dagTermInfo[,6]))[1] + 1)
  mm <- max(dagTermInfo$colorran)
  colorMap <- heat.colors(mm)
  nodeAttrs$fillcolor[allTerm] <- unlist(lapply(dagTermInfo$colorran, function(x) return(colorMap[x])))
  
  weightsList <- edgeWeights(DAG)
  to <- lapply(weightsList, names)
  from <- nodes(DAG)
  edge.names <- paste(rep(from, listLen(to)), unlist(to), sep = "~")
  edge.weights <- unlist(weightsList)
  names(edge.weights) <- edge.names
  ##    0 for a is_a relation,  1 for a part_of relation
  edgeAttrs$color <- ifelse(edge.weights == 0, 'black', 'red')
  filename <- paste(filename,".bmp",sep="")
  bmp(filename, width = bmpSize, height = bmpSize, res = 300, antialias = "none");
  plot(DAG, attrs = graphAttrs, nodeAttrs = nodeAttrs, edgeAttrs = edgeAttrs)  
  dev.off()
}