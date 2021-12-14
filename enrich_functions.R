suppressMessages({
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(GO.db)
  library(KEGG.db)
  library(reactome.db)
  library(AnnotationDbi)
  library(plyr)
  library(stats)
  library(jsonlite)
  library(gprofiler2)
  
})

#' Get enrichment.
#'
#' Submit a list of genes, specify the ogranism, and annotation type to get enrichment.
#'
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @importFrom AnnotationDbi select keys
#' @importFrom plyr ddply . summarise
#' @importFrom stats fisher.test
#'
#' @param genelist List of gene identifiers for which you want enriched GO.
#' @param annType Type of GO annotation 'GO', 'KEGG', or 'REACTOME' default:'GO'.
#' @param organism Organism 'hs' or 'mm' default:'hs'.
#' @param adjMethod Pvalue adjustment method 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none'. deafult:'BH'
#' @return Data frame containing enriched GO with their corresponding pValue, EASE score and list of genes representing GO
#' @examples
#' \dontrun{
#' annotation_enrichment(genelist=list_of_genes.df, annType="GO", organism="hs")
#' annotation_enrichment(genelist=list_of_genes.df)
#' }
#' @keywords internal
#' @export
annotation_enrichment <- function(genelist, keyType="SYMBOL", annType="GO", organism="hs", adjMethod="BH", background=NULL,pvalueCutoff=0.05){
  
  #Get organism annotation library
  orgLibs <- c("hsa"="org.Hs.eg.db", "mmu"="org.Mm.eg.db")
  orgDB <- orgLibs[organism]
  cat("Organism libary: ", orgDB, "\n")
  orgDB <- get(orgDB)
  
  #Get annotation library
  if(annType=="REACTOME"){
    annDB <- "reactome.db"
    cat("Annotation libary: ", annDB, "\n")
    annDB <- get(annDB)
  }else{
    cat("Annotation libary: ", orgLibs[organism], "\n")
    annDB <- orgDB
  }
  
  #Check gene identifier mapping to organism library
  idx <- which(genelist %in% keys(orgDB, keytype=keyType))
  if(length(idx)==0){
    stop("FAILED TO MAP ANY IDENTIFIER TO ORGANISM LIBRARY!")
  }else if(length(idx)<length(genelist)){
    cat("Unable to map ", length(genelist)-length(idx), " to organism library!")
    unmapped_gIDs <- genelist[-idx]
    genelist <- genelist[idx]
  }
  
  if(annType=="REACTOME" && keyType!="ENTREZID"){
    cat("Converting ", keyType, " to ENTREZID for mapping to reactome.db annotation library...\n")
    selectDF <- AnnotationDbi::select(orgDB, keys=genelist, columns="ENTREZID", keytype=keyType)
    selectDF <- selectDF[,c(1,2)]
    idx <- which(is.na(selectDF[,2]))
    if(length(idx)>0){
      cat("Remove NA. ", length(idx), " unmapped!")
      selectDF <- selectDF[-idx,]
    }
    genelist <- selectDF[,2]
    keyType <- "ENTREZID"
  }
  
  genelist.len <- length(genelist)
  annTypeMap <- c("GO"="GO", "KEGG"="PATH", "REACTOME"="PATHID")
  annType <- annTypeMap[annType]
  
  ##print("Creating Res DF...")
  selectDF <- AnnotationDbi::select(annDB, keys=genelist, columns=annType, keytype=keyType)
  ##print("head(selectDF):")
  ##print(head(selectDF))
  
  if(annType=="GO"){
    selectDF <- selectDF[selectDF$EVIDENCE != "ND",]
  }
  selectDF <- selectDF[,c(1,2)]
  
  if(length(which(is.na(selectDF[,2]))) > 0){
    ##print("Remove NA. Unmapped to annotation...")
    selectDF <- selectDF[-which(is.na(selectDF[,2])),]
  }
  
  colnames(selectDF) <- c("gID", "annID")
  ##print("head(selectDF):")
  ##print(head(selectDF))
  selectDF <- unique(selectDF)
  gID <- c()
  annID <- c()
  resDF <- plyr::ddply(selectDF, plyr::.(annID), plyr::summarise, paste(gID, collapse=","), length(gID))
  colnames(resDF) <- c("annID", "gID", "User_Genes")
  resDF$User_Genes_Not_In <- genelist.len-resDF$User_Genes
  ##print("Getting All Genes for mapped IDs...")
  ##print(class(resDF$annID))
  ##print(head(resDF$annID))
  selectAnnDF <- AnnotationDbi::select(annDB, keys=resDF$annID, columns="ENTREZID", keytype=annType)
  #selectAnnDF <- AnnotationDbi::select(annDB, keys=resDF$annID, columns="SYMBOL", keytype=annType)
  colnames(selectAnnDF)[1] <- "annID"
  ##print("Getting Gene count for each ID...")
  resDF$ALL_Genes <- plyr::ddply(selectAnnDF,.(annID),nrow)[,2]
  
  ##print("Getting Total Gene Count...")
  #genome.len <- length(unique(keys(annDB, keytype="ENTREZID")))
  
  if(!is.null(background)){
    genome.len <- length(unique(refMapDF[,"ENTREZID"]))
  }else {
    refMapDF <- select(annDB, keys=keys(annDB, annType), column="ENTREZID", keytype=annType)
    genome.len <- length(unique(refMapDF[,"ENTREZID"]))
    cat("Genome Length: ", genome.len, "\n")
    resDF$ALL_Genes_Not_In <- genome.len-resDF$ALL_Genes
  }
  
  ##print("Calculating PValue...")
  ##print("head(resDF):")
  ##print(head(resDF))
  resDF$pValue <- apply(resDF, 1, function(x){
    xx <- as.numeric(x[3:6])
    mat <- matrix(c(xx[1],xx[2],xx[3],xx[4]), nrow=2, dimnames=list(c("In", "Not_In"), c("User", "All")))
    testRes <- fisher.test(mat)
    testRes$p.value
  })
  ##print("Calculating Adjusted Pvalue...")
  resDF$pValueAdj <- p.adjust(resDF$pValue, method=adjMethod)
  ##print("Calculating EASE...")
  resDF$ease_score <- apply(resDF, 1, function(x){
    xx <- as.numeric(x[3:6])
    genesIn <- xx[1]-1
    mat <- matrix(c(genesIn,xx[2],xx[3],xx[4]), nrow=2, dimnames=list(c("In", "Not_In"), c("User", "All")))
    testRes <- stats::fisher.test(mat)
    testRes$p.value
  })
  
  if(annType=="GO"){
    selectAnnDF <- AnnotationDbi::select(GO.db, keys=as.vector(resDF$annID), columns=c("TERM", "ONTOLOGY"), keytype="GOID")
    resDF <- cbind(selectAnnDF, resDF[,c(2:ncol(resDF))])
  }
  ##print("Completed. Returning...")
  resDF = resDF[resDF$pValueAdj <= pvalueCutoff,]
  resDF
}

get_kegg_hierarchy <- function(){
  kegg_hierarchy <- data.frame(Level1=character(), Level2=character(), Pathway=character(), ID=numeric())
  keggJSON <- jsonlite::fromJSON("https://www.genome.jp/kegg-bin/download_htext?htext=br08901.keg&format=json&filedir=")
  for(i in 1:length(keggJSON$children$name)){
    l1 <- keggJSON$children$name[i]
    l2_DF <- keggJSON$children$children[[i]]
    for(j in 1:nrow(l2_DF)){
      l2 <- l2_DF[j,1]
      ###print(l2)
      pathTxt <- l2_DF[j,2][[1]][,1]
      ###print(length(pathTxt))
      pathTxtSplit <- strsplit(pathTxt, "  ")
      pathTxtDF <- as.data.frame(t(as.data.frame(pathTxtSplit)),row.names=c(1:length(pathTxt)))
      ###print(head(pathTxtDF))
      ###print(dim(pathTxtDF))
      colnames(pathTxtDF) <- c("ID", "Pathway")
      tmp <- data.frame(Level1=l1, Level2=l2, Pathway=pathTxtDF$Pathway, ID=pathTxtDF$ID, row.names=c(1:length(pathTxt)))
      ###print(dim(tmp))
      kegg_hierarchy <- rbind(kegg_hierarchy, tmp)
    }       
  }
  return(kegg_hierarchy)
}

enrich = function(x, type, org, pval, adjust_method,sig = TRUE, mis = 0, only_annotated = TRUE){
  if(only_annotated){
    domain_size = "annotated"
  }else{
    domain_size = "known"
  }
  
  print("before gprofiler")
  
  out = gprofiler2::gost(query = as.character(x[,1]),sources = type,organism = org,domain_scope = domain_size,
                          user_threshold = pval,correction_method = adjust_method,significant = sig, evcodes = TRUE)
  out = out$result
  
  if(is.null(out)) out = data.frame(query= character(), significant= integer(), p_value= integer(), term_size = integer(),            
                                    query_size =  integer(), intersection_size =  integer(),precision =  integer(),recall =  integer(),               
                                    term_id = character(),source = character(),term_name = character(),effective_domain_size = integer(),
                                    source_order = integer(),parents = character(),evidence_codes = character(),intersection = character())
  
  print(out)
  
  ### UPDATE gprofiler version
  # out = gProfileR::gprofiler(query = as.character(x[,1]),src_filter=type,organism=org,domain_size = domain_size,
  #                            max_p_value = pval, correction_method = adjust_method,significant = sig,min_isect_size = mis)

  # out=out[,c("term.id","intersection","p.value","p.value","term.name")]
  out=out[,c("term_id","intersection","p_value","p_value","term_name")]
  
  colnames(out) = c( "annID","gID","pValue","pValueAdj","Description")
  
  print(out)
  # if(nrow(out)>0){
  #   out$annID= gsub("KEGG:","",out$annID)
  #   if(org=="hsapiens")
  #     out$annID= gsub("REAC:","R-HSA-",out$annID)
  #   else{
  #     out$annID= gsub("REAC:","R-MMU-",out$annID)
  #   }
  
  if(nrow(out)>0)
    out$annID= gsub("KEGG:|REAC:","",out$annID)
  
  return(out)
}

convert_genes = function(organism = "hsapiens", GList, annType = "SYMBOL"){
  
  #save(GList,organism, annType,file = "inside_convert_genes.RData")
  library(org.Rn.eg.db)
  orgLibs <- list("Human"=org.Hs.eg.db, "Mouse"=org.Mm.eg.db, "Rat" = org.Rn.eg.db)
  orgDB <- orgLibs[[organism]]
  
  if(annType == "SYMBOL"){
    tmp = GList[[1]]
    tmp <- AnnotationDbi::select(orgDB, keys=as.character(tmp[,1]), columns="SYMBOL", keytype=annType)
    return(GList)
  }
  
  for(i in 1:length(GList)){
    M=GList[[i]]
    print("accessing M[,1]")
    print(class(M))
    print(dim(M))
    genes = M[,1]
    print("selected genes")
    print(genes)
    
    selectAnnDF <- AnnotationDbi::select(orgDB, keys=genes, columns="SYMBOL", keytype=annType)
    
    # if(annType=="SYMBOL"){
    #   selectAnnDF = cbind(selectAnnDF,selectAnnDF)
    #   rownames(selectAnnDF) = selectAnnDF$SYMBOL
    # }
    
    print(head(selectAnnDF))
    
    toRem = which(is.na(selectAnnDF[,2]))
    if(length(toRem)>0){
      selectAnnDF = selectAnnDF[-toRem,]
    }
    
    
    #remove eventual duplicates 
    selectAnnDF = selectAnnDF[!duplicated(selectAnnDF$SYMBOL),]
    
    M = M[M[,1] %in% selectAnnDF[,1],]
    
    print("----- > head M")
    print(head(M))
    
    M = as.matrix(M)
    print("----- > head M after dataframe")
    print(head(M))
    
    # macke sure they match by the key type
    matches = match(selectAnnDF[,1],M[,1])
    print("matches")
    print(matches)
    
    #all(M$ID[matches] == selectAnnDF[,1])
    M[matches,1] = selectAnnDF[,2]
    rownames(M) = M[,1]
    
    #update gene symbols
    GList[[i]] = M
    
    # M = M[which(M[,1] %in% selectAnnDF[,1]),]
    # GList[[i]] = M
  }
  
  return(GList)
}


# all_GO_BP = lapply(all_data,enrich,"GO:BP")
# all_GO_CC = lapply(all_data,enrich,"GO:CC")
# all_GO_MF = lapply(all_data,enrich,"GO:MF")
# 
# all_KEGG= lapply(all_data,enrich,"KEGG")
# all_REACT = lapply(all_data,enrich,"REAC")


