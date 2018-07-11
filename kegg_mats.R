library(clusterProfiler)
library(gtools)
library(readxl)

read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

#Funzione che restituisce liste di dataframe
#funzione che prende li liste di dataframe e si crea la matrice dei pathways

# gene_sig is the list with the genes for each treatment. It has a vector of genes for each treatment
# key_type is the type of genes taken in inpyt #SYMBOL, ENTREZGENE, ecc
# annType is the type of enrichment we want to perform (eg. GO, KEGG, REACTOME)
build_dataframe_enrichment <- function(gene_sig, organism = 'hsa', pvalueCutoff = 0.05,pAdjustMethod = "fdr",keyType="SYMBOL", annType="GO") {
  conmp_names = names(gene_sig)
  #List that will contains the results of the enrichment
  EnrichDatList = list()
  for (i in 1:length(gene_sig)){
    kegg <- NA
    entrez_sym <- as.character(gene_sig[[i]])
    entrez_sym[entrez_sym %in% ""] = NA
    entrez_sym <- entrez_sym[complete.cases(entrez_sym)]
    if(!is.na(entrez_sym) && length(entrez_sym)>0){
      try(kegg <- annotation_enrichment(genelist =na.omit(entrez_sym) , keyType=keyType, annType=annType, organism=organism, adjMethod=pAdjustMethod,pvalueCutoff = pvalueCutoff))
      if(!gtools::invalid(kegg)){
        EnrichDatList[[conmp_names[[i]]]] = kegg
      }else{
        EnrichDatList[[conmp_names[[i]]]] = NULL
      }
    }
  }
  return(EnrichDatList)
}

#This function filter the GO enrichment results only for the ontology specified by the user
filterGO = function(EnrichDatList,go_type="BP"){
  goTerm =c()
  for(i in 1:length(EnrichDatList)){
    df = EnrichDatList[[i]]
    # idx = which(df$ONTOLOGY %in% go_type)
    # if(length(idx)>0){
    #   EnrichDatList[[i]] = df[idx,]
    #   goTerm = union(goTerm,df[idx,"GOID"])
    # }else{
    #   EnrichDatList[[i]] = NULL
    # }
    goTerm = union(goTerm,df[,"annID"])
   
  }
  
  XX = select(x = GO.db,columns = c("GOID","ONTOLOGY"),keys = goTerm)
  toRem = which(is.na(XX[,2]))
  if(length(toRem)>0){
    toRemGo = XX[toRem,1]
    for(i in 1:length(EnrichDatList)){
      df = EnrichDatList[[i]]
      toRem2 = which(df$annID %in% toRemGo)
      if(length(toRem2)>0){
        EnrichDatList[[i]] = EnrichDatList[[i]][-toRem2,]
      }
    }
  }
  
  if(length(toRem)>0){
    goTerm = XX[-toRem,1]
    
  }
  
  return(list(EnrichDatList=EnrichDatList,goTerm=goTerm))
}


# it takes in input a table with samples on the columns and genes in the rows and return a list long N (with N the number of samples)
# In every i-th position of the list there are the gene of sample i-th
create_list_of_genes_from_data = function(DAT, isHeader = TRUE, file_type){
  ##print(dim(DAT))
  
  if(isHeader){
    #print("The file has colnames")
    name_contrasts = colnames(DAT)
  }else{
    #print("The file has no colnames")
    name_contrasts = paste("V",1:ncol(DAT),sep="")
  }
  
  if(file_type=="GenesOnly"){
    #list with ensemble gene id for each contrast
    LIST = list()
    for(i in 1:ncol(DAT)){
      genes = DAT[,i]
      LIST[[name_contrasts[i]]] = genes
    }
    return(LIST)
  }else{
    LIST1 = list()
    LIST2 = list()
    
    for(i in 1:ncol(DAT)){
      if(i %% 2 == 1){
        genes = DAT[,i]
        LIST1[[name_contrasts[i]]] = genes
      }else{
        LIST2[[name_contrasts[i]]] = DAT[is.na(DAT[,i])==FALSE,i]
      }
    }
    names(LIST2) = names(LIST1)
    return(list(LIST1 = LIST1, LIST2 = LIST2))
  }
}

# EnrichDatList is the list of dataframe with the results of the enrichment for every sample
kegg_mat_p <- function(EnrichDatList,hierarchy) {
  cat("Inside kegg_mat_p")
  conmp_names = names(EnrichDatList)
  cat("Inside kegg_mat_p GO")
  names_un = unique(hierarchy$ID)
  kegg_mat_cell <- matrix(data = NA, nrow = length(conmp_names),ncol = length(names_un), dimnames = list(conmp_names,names_un))
  
  for (i in 1:length(EnrichDatList)){
    kegg = EnrichDatList[[i]]
    if(!gtools::invalid(kegg)){
      if (nrow(kegg)>0)
        for(kg in 1:nrow(kegg)){
          if(kegg$annID[kg] %in% colnames(kegg_mat_cell) ){
            kegg_mat_cell[conmp_names[[i]],kegg$annID[kg]] <- kegg$pValueAdj[kg]
          }else{
            #print(paste(kegg$annID[kg],"not found in the provided kegg hierarchy"))
          }
        }
    }
  }  
  
  kegg_mat_cell <- kegg_mat_cell[,colSums(!is.na(kegg_mat_cell)) >0]
  
  #kegg_mat_cell <- kegg_mat_cell[,colSums(!is.na(kegg_mat_cell)) < nrow(kegg_mat_cell)]
  return(kegg_mat_cell)
}

kegg_mat_fc <- function(EnrichDatList,hierarchy,GList, summ_fun=median) {
  cat("Inside kegg_mat_fc")
  conmp_names = names(EnrichDatList)
  names_un = unique(hierarchy$ID)
  kegg_mat_cell <- matrix(data = NA, nrow = length(conmp_names),ncol = length(names_un), dimnames = list(conmp_names,names_un))
  
  ##print(head(kegg_mat_cell))
  
  for (i in 1:length(EnrichDatList)){
    kegg = EnrichDatList[[i]]
    if(!gtools::invalid(kegg)){
      if (nrow(kegg)>0)
        for(kg in 1:nrow(kegg)){
          if(kegg$annID[kg] %in% colnames(kegg_mat_cell) ){
            genes_in_path <- unlist(strsplit(as.character(kegg$gID[kg]),","))
            MM = GList[[i]]
            summFC <- summ_fun(MM[tolower(MM[,1]) %in% tolower(genes_in_path),2])
            kegg_mat_cell[conmp_names[[i]],kegg$annID[kg]] <- summFC
          }else{
            #print(paste(kegg$annID[kg],"not found in the provided kegg hierarchy"))
          }
        }
    }
  }  
  kegg_mat_cell <- kegg_mat_cell[,colSums(!is.na(kegg_mat_cell)) >0]
  
  #kegg_mat_cell <- kegg_mat_cell[,colSums(kegg_mat_cell,na.rm = T)>0]
  return(kegg_mat_cell)
}

#plots all 3 collapsed leayers of a kegg matrix and the sub matrices obtained by splitting at level 1
plot_kegg_mat <- function(kegg_hierarchy, kegg_mat_cell, group_col,discrete=T, pre_title="",square_colors=c(),color_leg=c(),path_text_size=10,treat_text_size=10) {
  
  #collapse at level 1 using the median to summarize vaulues
  kegg_nano_1 <- collapse_paths(kegg_hierarchy = kegg_hierarchy,kegg_mat_cell = kegg_mat_cell, collapse_level = 1)
  #extract collapsed matrix and collapsed hierarachy
  mat <- kegg_nano_1[[1]]
  hier <- kegg_nano_1[[2]]
  #plot the collapsed matrix
  plot_grid(path_mat = mat,path_hier = hier,experiment_ann = group_col,discrete =  discrete,level_col = 1,title = paste(pre_title,"level1"),square_colors,color_leg,path_text_size = path_text_size,treat_text_size = treat_text_size)
  
  #collapse at level 1 using the median to summarize vaulues
  kegg_nano_2 <- collapse_paths(kegg_hierarchy = kegg_hierarchy,kegg_mat_cell = kegg_mat_cell,collapse_level = 2)
  mat <- kegg_nano_2[[1]]
  hier <- kegg_nano_2[[2]]
  plot_grid(path_mat = mat,path_hier = hier,experiment_ann = group_col,discrete =  discrete,level_col = 1,title = paste(pre_title,"level2"),square_colors,color_leg,path_text_size = path_text_size,treat_text_size = treat_text_size)
  
  #collapse at level 3 using the median to summarize vaulues
  #this call is on the last levele hence no summarizatrion is done
  #we exploit the side effect of reordering and filtering of the hierarchy
  kegg_nano_3 <- collapse_paths(kegg_hierarchy = kegg_hierarchy,kegg_mat_cell = kegg_mat_cell, collapse_level = 3)
  mat <- kegg_nano_3[[1]]
  hier <- kegg_nano_3[[2]]
  plot_grid(path_mat = mat,path_hier = hier,experiment_ann = group_col, discrete =  discrete,level_col = 1, title = paste(pre_title,"level3"),square_colors,color_leg,path_text_size = path_text_size,treat_text_size = treat_text_size)
  
  #split the matrix with levels defined by level 1 
  path_by_lev_list <- paths_bylev(kegg_hierarchy = kegg_hierarchy, kegg_mat_cell = kegg_mat_cell, split_level = 1)
  # plot sub matrices one at time grouping and coloring at level 2
  for (i in 1:length(path_by_lev_list)){
    mat <- path_by_lev_list[[i]][[1]]
    hier <- path_by_lev_list[[i]][[2]]
    plot_grid(path_mat = mat,path_hier = hier,experiment_ann = group_col,discrete =  discrete,level_col = 2,title = paste(pre_title,names(path_by_lev_list)[i]),square_colors,color_leg,path_text_size = path_text_size,treat_text_size = treat_text_size)
  }
}



# # EnrichDatList is the list of dataframe with the results of the enrichment for every sample
# kegg_mat_p <- function(EnrichDatList,kegg_hierarchy,mm_reactome_hierarchy,mouse_map,hm_reactome_hierarchy,human_map, go_hierarchy, org = "mm",annType="GO",go_type = "BP") {
#   cat("Inside kegg_mat_p")
#   conmp_names = names(EnrichDatList)
#   
#   if(annType=="GO"){
#     cat("Inside kegg_mat_p GO")
#     names_un = unique(go_hierarchy[,3])
#     kegg_mat_cell <- matrix(data = NA, nrow = length(conmp_names),ncol = length(names_un), dimnames = list(conmp_names,names_un))
#     
#     for (i in 1:length(EnrichDatList)){
#       kegg = EnrichDatList[[i]]
#       if(!gtools::invalid(kegg)){
#         if (nrow(kegg)>0)
#           for(kg in 1:nrow(kegg)){
#             if(kegg$TERM[kg] %in% colnames(kegg_mat_cell) ){
#               kegg_mat_cell[conmp_names[[i]],kegg$TERM[kg]] <- kegg$pValueAdj[kg]
#             }else{
#               #print(paste(kegg$TERM[kg],"not found in the provided kegg hierarchy"))
#             }
#           }
#       }
#     }
#       
#   }
#   if(annType == "KEGG"){
#     cat("Inside kegg_mat_p KEGG")
#     
#     names_un = unique(kegg_hierarchy$Pathway)
#     kegg_mat_cell <- matrix(data = NA, nrow = length(conmp_names),ncol = length(names_un), dimnames = list(conmp_names,names_un))
#     
#     for (i in 1:length(EnrichDatList)){
#       kegg = EnrichDatList[[i]]
#       if(!gtools::invalid(kegg)){
#         
#         kegg$Description = kegg_hierarchy$Pathway[kegg_hierarchy$ID %in% kegg$annID]
#         
#         if (nrow(kegg)>0)
#           for(kg in 1:nrow(kegg)){
#             if(kegg$Description[kg] %in% colnames(kegg_mat_cell) ){
#               kegg_mat_cell[conmp_names[[i]],kegg$Description[kg]] <- kegg$pValueAdj[kg]
#             }else{
#               #print(paste(kegg$Description[kg],"not found in the provided kegg hierarchy"))
#             }
#           }
#       }
#     }
#   }
#   if(annType == "REACTOME"){
#     cat("Inside kegg_mat_p REACTOME")
#     
#     if(org == "mm"){
#       reactome_hierarchy = mm_reactome_hierarchy
#       names_un = mouse_map[unlist(reactome_hierarchy$Pathway),2]
#       
#     }else{
#       reactome_hierarchy = hm_reactome_hierarchy
#       names_un = human_map[unlist(reactome_hierarchy$Pathway),2]
#       
#     }
#     
#     #names_un = unlist(unique(reactome_hierarchy$Pathway))
#     kegg_mat_cell <- matrix(data = NA, nrow = length(conmp_names),ncol = length(names_un), dimnames = list(conmp_names,names_un))
#     
#     for (i in 1:length(EnrichDatList)){
#       kegg = EnrichDatList[[i]]
#       #print(i)
#       if(!gtools::invalid(kegg)){
#         
#         kegg$Description = mouse_map[kegg$annID,2]#reactome_hierarchy$Pathway[kegg_hierarchy$ID %in% kegg$annID]
#         
#         if (nrow(kegg)>0)
#           for(kg in 1:nrow(kegg)){
#             if(kegg$Description[kg] %in% colnames(kegg_mat_cell) ){
#               kegg_mat_cell[conmp_names[[i]],kegg$Description[kg]] <- kegg$pValueAdj[kg]
#             }else{
#               #print(paste(kegg$Description[kg],"not found in the provided kegg hierarchy"))
#             }
#           }
#       }
#     }
#   
#   }
#   
#   kegg_mat_cell <- kegg_mat_cell[,colSums(kegg_mat_cell,na.rm = T)>0]
#   return(kegg_mat_cell)
# }



# kegg_mat_p <- function(conmp_names, gene_sig,kegg_hierarchy, organism = 'hsa', pvalueCutoff = 0.05,pAdjustMethod = "fdr",keyType="SYMBOL", annType="GO") {
#   names_un = unique(kegg_hierarchy$Pathway)
#   kegg_mat_cell <- matrix(data = NA, nrow = length(conmp_names),ncol = length(names_un), dimnames = list(conmp_names,names_un))
# 
#   for (i in 1:length(gene_sig)){
#     entrez_sym <- NA
#     kegg <- NA
#     
#     entrez_sym <- as.character(gene_sig[[i]])
#     entrez_sym[entrez_sym %in% ""] = NA
#     entrez_sym <- entrez_sym[complete.cases(entrez_sym)]
#     
#     if(!is.na(entrez_sym) && length(entrez_sym)>0){
#       #try(kegg2 <- enrichKEGG(gene = na.omit(entrez_sym), organism = organism, pvalueCutoff = pvalueCutoff,pAdjustMethod = pAdjustMethod))
#       try(kegg <- annotation_enrichment(genelist =na.omit(entrez_sym) , keyType=keyType, annType=annType, organism=organism, adjMethod=pAdjustMethod,pvalueCutoff = pvalueCutoff))
#       if(!gtools::invalid(kegg)){
#         #kegg <- kegg@result
#         kegg$Description = kegg_hierarchy$Pathway[kegg_hierarchy$ID %in% kegg$annID]
#         
#         if (nrow(kegg)>0)
#           for(kg in 1:nrow(kegg)){
#             if(kegg$Description[kg] %in% colnames(kegg_mat_cell) ){
#               kegg_mat_cell[conmp_names[[i]],kegg$Description[kg]] <- kegg$pValueAdj[kg]
#             }else{
#               #print(paste(kegg$Description[kg],"not found in the provided kegg hierarchy"))
#             }
#           }
#       }
#     }
#   }
#   
#   kegg_mat_cell <- kegg_mat_cell[,colSums(kegg_mat_cell,na.rm = T)>0]
#   return(kegg_mat_cell)
# }
# 
# kegg_mat_fc <- function(conmp_names, gene_sig, gene_sig_fc, kegg_hierarchy, discr=T, organism = 'hsa', pvalueCutoff = 0.05,pAdjustMethod = "fdr", summ_fun=median,keyType="SYMBOL", annType="GO") {
#   names_un = unique(kegg_hierarchy$Pathway)
#   kegg_mat_FC <- matrix(data = NA, nrow = length(conmp_names),ncol = length(names_un), dimnames = list(conmp_names,names_un))
#   
#   for (i in 1:length(gene_sig)){
#     kegg <- NA
#     genes_in_path <- NA
#     summFC <- NA
# 
#     if(!all(is.na(gene_sig[[i]])) && length(gene_sig[[i]])>0){
#       #try(kegg <- enrichKEGG(gene = na.omit(gene_sig[[i]]), organism = organism, pvalueCutoff = pvalueCutoff,pAdjustMethod = pAdjustMethod))
#       try(kegg <- annotation_enrichment(genelist =na.omit(entrez_sym) , keyType=keyType, annType=annType, organism=organism, adjMethod=pAdjustMethod,pvalueCutoff = pvalueCutoff))
#       
#       if(!gtools::invalid(kegg)){
#         #kegg <- kegg@result
#         kegg$Description = kegg_hierarchy$Pathway[kegg_hierarchy$ID %in% kegg$annID]
#         
#         if (nrow(kegg)>0)
#           for(kg in 1:nrow(kegg)){
#             if(kegg$Description[kg] %in% colnames(kegg_mat_FC) ){
#               genes_in_path <- unlist(strsplit(as.character(kegg$gID[kg]),","))
#               summFC <- summ_fun(gene_sig_fc[[i]][gene_sig[[i]] %in% genes_in_path])
#               kegg_mat_FC[conmp_names[[i]],kegg$Description[kg]] <- summFC
#             }else{
#               #print(paste(kegg$Description[kg],"not found in the provided kegg hierarchy"))
#             }
#           }
#       }
#     }
#   }
#   
#   kegg_mat_FC <- kegg_mat_FC[,colSums(!is.na(kegg_mat_FC))>0]
#   
#   if (discr){
#     kegg_mat_FC[kegg_mat_FC>0] <- 1
#     kegg_mat_FC[kegg_mat_FC<0] <- -1
#   }
#   return(kegg_mat_FC)
# }