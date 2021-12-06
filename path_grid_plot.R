# utility function to cut a string to the desired number of characters
shorten <- function(string) {
  #substr(string, 1, 0)
  string
}

#splits a matrix of pathways X treatments based on the provided hierarchy and the deired level of separation
#return a list of (matrix, subHierarchy) one element for each distinct level

paths_bylev <- function (kegg_hierarchy,kegg_mat_cell,split_level=1) {
  levels_mat <- list()
  #select pathways that are present in the matrix
  hierarchy_sub <- kegg_hierarchy[kegg_hierarchy$Pathway %in% colnames(kegg_mat_cell),]
  #order the matrix columns following the provided hyerarchy
  mat_cell_ord <- kegg_mat_cell[,as.character(hierarchy_sub$Pathway)]

  #split the input matrix using pathways groups defined by the provided split_level in the hierarchy
  for(level in unique(hierarchy_sub[,split_level])){
    #take sub-hierarchy rooted in level
    lev_hierarchy <- hierarchy_sub[hierarchy_sub[,split_level]==level,]
    #extract columns with pathways in the defined sub-hierarchy
    kegg_lev <- mat_cell_ord[,colnames(mat_cell_ord) %in% lev_hierarchy$Pathway]
    if (sum(colnames(mat_zcell_ord) %in% lev_hierarchy$Pathway)==1){
      kegg_lev <- as.matrix(kegg_lev)
      colnames(kegg_lev) <- lev_hierarchy$Pathway
      }
    #create a nemed list element containing sub matrix and sub pathway
    levels_mat[[length(levels_mat)+1]] <- list(kegg_lev,lev_hierarchy)
    names(levels_mat)[length(levels_mat)] <- level

  }
  return(levels_mat)
}

update_hierarchy=function(kegg_hierarchy,lev1_content,lev2_content,lev3_content){
   # if(is.na(lev1_content) == F) print("selected lev1 ", lev1_content ,"\n")
   # if(is.na(lev2_content) == F) print("selected lev2 ", lev2_content ,"\n")
   # if(is.na(lev3_content) == F) print("selected lev3 ", lev3_content ,"\n")

  print("I'm updating the hierarchy...")
  if("All" %in% lev1_content){
    print("All in lev1")
    idx1 = 1:nrow(kegg_hierarchy)
  }else{
    idx1= which(kegg_hierarchy[,1] %in% lev1_content)
  }
  if("All" %in% lev2_content){
    print("All in lev2")
    idx2 = 1:nrow(kegg_hierarchy)
  }else{
    idx2= which(kegg_hierarchy[,2] %in% lev2_content)
  }
  if("All" %in% lev3_content){
    print("All in lev3")
    idx3 = 1:nrow(kegg_hierarchy)
  }else{
    idx3= which(kegg_hierarchy[,3] %in% lev3_content)
  }

  idx = intersect(idx1,intersect(idx2,idx3))
  #print("Selected idx: ")
  #print(idx)

  return(kegg_hierarchy[idx,])

}

#collapse a matrix of pathways X treatments based on the provided hierarchy and the deired level
#return a list of (matrix, subHierarchy) one element for each distinct collapsing level
# the value is computed summarising row elements with the provided function col_fun (default is median)
collapse_paths <- function (kegg_hierarchy,kegg_mat_cell, collapse_level=1,col_fun=function(x){median(x,na.rm = T)}) {
  #save(kegg_hierarchy,kegg_mat_cell,collapse_level,col_fun, file="demo/collapse_path.RData")

  #select pathways that are present in the matrix
  hierarchy_sub <- kegg_hierarchy[kegg_hierarchy$ID %in% colnames(kegg_mat_cell),]
  #order the matrix columns following the provided hyerarchy
  mat_cell_ord <- kegg_mat_cell[,as.character(hierarchy_sub$ID)]

  #get the list of categories on which we will collapse our pathways
  #keg_lev <- unique(hierarchy_sub[,collapse_level])
  if(collapse_level==1){
    keg_lev <- unique(hierarchy_sub[,collapse_level])
  }else{
    if(collapse_level==3){
      keg_lev <- hierarchy_sub[,collapse_level]
    }
    else{
      keg_lev <- unique(hierarchy_sub[,1:collapse_level])
      keg_lev = keg_lev[,2]
    }
  }


  #create a new matrix with columns corresponding to collapse categories
  if(is.null(nrow(mat_cell_ord))){
    mat_cell_ord = matrix(mat_cell_ord,ncol=1,nrow=length(mat_cell_ord),dimnames = list(as.character(names(mat_cell_ord)),hierarchy_sub$ID))
  }
  lev_mat <- matrix(NA,ncol = length(keg_lev),nrow = nrow(mat_cell_ord),dimnames = list(rownames(mat_cell_ord),keg_lev))

  #for each category take all pathways in that category and summaryze their values for each row (sample)
  for (lev in keg_lev){ #controllare che il padre e' arricchito. se e' arricchito ci metti il padre (ci metti asterisco), altrimenti fai il summary
    #take sub-hierarchy rooted in lev
    path_lev <- unique(hierarchy_sub[hierarchy_sub[,collapse_level] == lev, ]$ID)
    #extract columns with pathways in the defined sub-hierarchy
    mat_cell_sub <- as.matrix(mat_cell_ord[,colnames(mat_cell_ord) %in% path_lev])
    #summarize vaules over rows using the provided col_fun function
    summary <- unlist(apply(mat_cell_sub,1,FUN = col_fun))
    #store summarized column for culumn lev
    ## if duplicate column names in mat it was assigning only the first
    # lev_mat[,lev] <- summary

    lev_mat[,colnames(lev_mat) %in% lev] <- summary
  }
  #create a nemed list element containing collapsed matrix and collapsed hierarchy as well
  lev_hierarchy <- unique(as.data.frame(hierarchy_sub[,1:collapse_level]))
  colnames(lev_hierarchy) <- colnames(hierarchy_sub)[1:collapse_level]

  return(list(lev_mat,lev_hierarchy))
}



#plots a matrix having colums for wich a hierachy is provided
# path_mat is matrix to be plotted
# path_hier is the hierarchy defined over columns
# experiment_ann is a vector of the same length as # of rows of path_mat defining grouping for samples
# discrete tells the function if the value are to be plotted using a countinuos scale or a discrete scale
# square_colors if a discrete scale is chosen, than the colours for each possible value have to be provided
# color_leg     if a discrete scale is chosen, than the colour legend for each possible value has to be provided
# level_col   level (column number) from the hierarchy used to group columns (pathways)

plot_grid <- function(path_mat,path_hier, title="", experiment_ann=c(),discrete=F,square_colors=c(),color_leg=c(),level_col=1,treat_text_size=8,path_text_size=6, asRatio = FALSE) {
  #save(path_mat, path_hier,experiment_ann,title,discrete,square_colors,color_leg,
  #     level_col,treat_text_size,path_text_size,file="demo/demo_plot.RData")

  #path_mat = path_mat[rownames(experiment_ann),]
  #define the groups from the hierarchy and the chosen level
  path_col <- factor(path_hier[,level_col], levels = unique(path_hier[,level_col]))
  #prepare a set of colors to assign to each group
  #darkcols <- c(brewer.pal(8, "Dark2"),brewer.pal(8, "Accent")[-4],brewer.pal(12, "Paired")[-11])

  #using ggplot2 we need to melt the input matrix
  kegg_melt <- melt(path_mat)
  colnames(kegg_melt)[1:2] = c("Var1","Var2")
  #if discrete create value-color mapping using input values or this default scale
  if (discrete){
    if(length(square_colors)==0){
      colors <- c("-1"="darkgreen","0"="white","1"="red")
      color_leg <- c("-1"= "negative","0"="neutral","1"="positive")
    }else {colors <- square_colors}
    #if discrete we need to convert values in the matrix in factors to esure ggplot will plot on a discrete scale
    kegg_melt$value <- as.factor(kegg_melt$value)
  }
  # if a grouping for the sample is provided assign it to the melted dataframe rows (ncol by ncol)
  #this will be used in faceting
  if(length(experiment_ann[,1])>0){
    #print(path_mat)
    #print(kegg_melt)
    #print(experiment_ann)
    #print(ncol(path_mat))
    x = experiment_ann[which(experiment_ann[,2] %in% rownames(path_mat)),1]
    #print(x)
    kegg_melt$experiment <- rep(x,ncol(path_mat))
  }else{
    #otherwise assign a default group to all samples
    kegg_melt$experiment <- as.factor("treatment")
  }

  #assign group fromf rom the hierarchy and the chosen level (path_col) to the melted dataframe rows (nrow by nrow)
  #this will be used in faceting
  kegg_melt$path_group <- factor(as.character(rep(path_col,each=nrow(path_mat))),levels=unique(as.character(rep(path_col,each=nrow(path_mat)))))

  #prepare a set of colors to assign to each group
  darkcols <- randomColor(length(unique(kegg_melt$path_group)),luminosity = "dark")


  # prepare the ggplot instance -- to render groups we will use faceting
  #using tile plot to render thje matrix as a grid of squares



  #print(kegg_melt$experiment)
  print(head(kegg_melt))
  ukm = unique(kegg_melt$experiment)
  print(ukm)

  for(i in 1:nrow(kegg_melt)){
    kegg_melt[i,"experiment"] = which(ukm %in% kegg_melt[i,"experiment"] )
  }
  print(head(kegg_melt))

  #Ordinamento per mantenere ordinamento i gruppi ordinati secondo ordine numerico
  labelx = sort(unique(as.numeric(as.character(kegg_melt$experiment))))
  print(labelx)
  kegg_melt$experiment = factor(kegg_melt$experiment,levels = labelx)

  #Ordiniamo i campioni con lo stesso ordine delle foglie dell'albero gerarchico
  kegg_melt$Var1 = factor(kegg_melt$Var1,levels = unique(kegg_melt$Var1))

  if(discrete == FALSE){
    #the values in kegg_melt must be numeric, not factors
    kegg_melt$value = as.numeric(as.vector(kegg_melt$value))
  }


  p <- ggplot(kegg_melt, aes(Var1, Var2, fill = value)) +
    facet_grid(path_group~experiment,scales="free",space="free",labeller = labeller(path_group=shorten)) +
    # theme(axis.text.x = element_text(colour = as.numeric(as.factor(experiment_ann))),
    #       axis.text.y = element_text(colour = darkcols[as.numeric(path_col)]))   +
    geom_tile(colour = "black",size=0.1) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    #coord_equal() +
    theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=0, size=treat_text_size),
          axis.text.y = element_text(size=path_text_size),
          #plot.margin = margin(1, 1, 1, 1, "cm"),
          strip.text.y = element_text(angle = 0,size = 9,face = "bold")) +
          scale_x_discrete(position = "top") +
          labs(x = "", y="")

    if(asRatio) p = p+theme(aspect.ratio = 1)


  if (discrete){
    p <- p + scale_fill_manual(values=colors,labels=color_leg,na.value = 'gray50')
  }else{
    #limits=c(min(kegg_melt$value[is.na(kegg_melt$value)==FALSE]), max(kegg_melt$value[is.na(kegg_melt$value)==FALSE]))
    p <- p + scale_fill_gradient(low = "darkgreen", high = "red",na.value = 'gray50')
  }

  #afrer faceting we need to open the ggplot object to put different row/column label colors to different facets
  png(file.path(tempdir(),"veer.png"))
  gplot <- ggplotGrob( p )
  dev.off()
  nms <- lapply( gplot$grobs , function(x) names( x[]$children ) )
  #we search for axis.line.x and axis.line.y containers
  grbs_x_id <- which( sapply( lapply( nms , function(x) grepl( "axis.line.x" , x ) ) , any ) == 1 )
  grbs_y_id <- which( sapply( lapply( nms , function(x) grepl( "axis.line.y" , x ) ) , any ) == 1 )
  #grob variable iterates over blocks of colums of the ggplot
  #for (grob in (1:length(grbs_x_id))){
  #  gplot$grobs[[grbs_x_id[grob]]]$children$axis$grobs[[1]]$children[[1]]$gp$col=grob
  #}
  #grob variable iterates over blocks of rows of the ggplot
  #for (grob in (1:length(grbs_y_id))){
  #  gplot$grobs[[grbs_y_id[grob]]]$children$axis$grobs[[1]]$children[[1]]$gp$col=darkcols[grob]
  #}
     #grob variable iterates over blocks of colums of the ggplot
  for (grob in (1:length(grbs_x_id))){
    to_select <- grbs_x_id[grob]
    if (length(to_select) >= 1) { # Check added because there was a 0 length vector. Luca
      item_at_grob <- gplot$grobs[[to_select]]
      if (! is.null(item_at_grob)) { # Check added because there where nulls. Luca
        item_at_grob$children$axis$grobs[[1]]$children[[1]]$gp$col=grob
      }
    }
  }
  #grob variable iterates over blocks of rows of the ggplot
  for (grob in (1:length(grbs_y_id))){
    to_select <- grbs_y_id[grob]
    if (length(to_select) >= 1) {
      item_at_grob <- gplot$grobs[[to_select]]
      if (! is.null(item_at_grob)) {
        item_at_grob$children$axis$grobs[[1]]$children[[1]]$gp$col=darkcols[grob]
      }
    }
  }


#  save(list = ls(all.names = TRUE),file = "ggplotting")

  #plot the final result
  #plot(gplot) #Commented by Veer
  return(gplot)
}


#plots all 3 collapsed leayers of a kegg matrix and the sub matrices obtained by splitting at level 1
plot_kegg_mat <- function(kegg_hierarchy, kegg_mat_cell, group_col,discrete=T, pre_title="",square_colors=c(),color_leg=c(),path_text_size=12,treat_text_size=12,lev_3_cex_row=10,lev_3_cex_col=10) {

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
  plot_grid(path_mat = mat,path_hier = hier,experiment_ann = group_col, discrete =  discrete,level_col = 1, title = paste(pre_title,"level3"),square_colors,color_leg,path_text_size = lev_3_cex_row,treat_text_size = lev_3_cex_col)

  #split the matrix with levels defined by level 1
  path_by_lev_list <- paths_bylev(kegg_hierarchy = kegg_hierarchy, kegg_mat_cell = kegg_mat_cell, split_level = 1)
  # plot sub matrices one at time grouping and coloring at level 2
  for (i in 1:length(path_by_lev_list)){
    mat <- path_by_lev_list[[i]][[1]]
    hier <- path_by_lev_list[[i]][[2]]
    plot_grid(path_mat = mat,path_hier = hier,experiment_ann = group_col,discrete =  discrete,level_col = 2,title = paste(pre_title,names(path_by_lev_list)[i]),square_colors,color_leg,path_text_size = path_text_size,treat_text_size = treat_text_size)
  }
}

#Plot grid map for genes of a selected pathway
plot_grid_genes <- function(path_mat, title="", experiment_ann=c(), gene_group=NULL, discrete=F, square_colors=c(), color_leg=c(), treat_text_size=8, asRatio=TRUE){

  kegg_melt <- melt(path_mat)
  colnames(kegg_melt)[1:2] = c("Var1","Var2")

  #if discrete create value-color mapping using input values or this default scale
  if (discrete){
    if(length(square_colors)==0){
      colors <- c("-1"="darkgreen","0"="white","1"="red")
      color_leg <- c("-1"= "negative","0"="neutral","1"="positive")
    }else{
      colors <- square_colors
    }
    #if discrete we need to convert values in the matrix in factors to esure ggplot will plot on a discrete scale
    kegg_melt$value <- as.factor(kegg_melt$value)
  }

  # if a grouping for the sample is provided assign it to the melted dataframe rows (ncol by ncol)
  #this will be used in faceting

  print("experiment_ann")
  print(experiment_ann)
  if(length(experiment_ann[,1])>0){
    print("Experiment is present!")
    #print(path_mat)
    #print(kegg_melt)
    #print(experiment_ann)
    #print(ncol(path_mat))

    x = experiment_ann[which(experiment_ann[,2] %in% rownames(path_mat)),1]
    #print(x)

    print("dim(kegg_melt)")
    print(dim(kegg_melt))

    print("str(kegg_melt)")
    print(str(kegg_melt))

    print("str(rep(x,ncol(path_mat)))")
    print(str(rep(x,ncol(path_mat))))

    kegg_melt$experiment <- rep(x,ncol(path_mat))
  }else{
    print("Creating fake 'experiment'!")
    #otherwise assign a default group to all samples
    kegg_melt$experiment <- as.factor("treatment")
  }

  #using tile plot to render thje matrix as a grid of squares

  #print(kegg_melt$experiment)
  print("head(kegg_melt)")
  print(head(kegg_melt))
  ukm = unique(kegg_melt$experiment)
  print("Unique experiments:")
  print(ukm)

  for(i in 1:nrow(kegg_melt)){
    kegg_melt[i,"experiment"] = which(ukm %in% kegg_melt[i,"experiment"])
  }
  print("head(kegg_melt)")
  print(head(kegg_melt))

  #Ordinamento per mantenere ordinamento i gruppi ordinati secondo ordine numerico

  labelx = sort(unique(as.numeric(as.character(kegg_melt$experiment))))
  print("labelx")
  print(labelx)

  kegg_melt$experiment = factor(kegg_melt$experiment, levels=labelx)

  #Ordiniamo i campioni con lo stesso ordine delle foglie dell'albero gerarchico

  kegg_melt$Var1 = factor(kegg_melt$Var1, levels=unique(kegg_melt$Var1))

  if(discrete == FALSE){
    #the values in kegg_melt must be numeric, not factors
    kegg_melt$value = as.numeric(as.vector(kegg_melt$value))
  }

  if(!is.null(gene_group)){
    print("class(gene_group)")
    print(class(gene_group))

    print("str(gene_group)")
    print(str(gene_group))

    # print("gene_group")
    # print(gene_group)

    kegg_melt$gene_group = gene_group
  }

  #Make pathway melt
  kegg_melt_path <- dplyr::filter(kegg_melt, gene_group=="Pathway")

  print("str(kegg_melt_path)")
  print(str(kegg_melt_path))

  #Pathway map
  pPath <- ggplot(kegg_melt_path, aes(Var1, Var2, fill=value)) +
    facet_grid(.~experiment, scales="free", space="free") +
    geom_tile(colour="black", size=0.1) +
    ggtitle(title) +
    theme(plot.title=element_text(hjust=0.5)) +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=0, size=treat_text_size),
          axis.text.y = element_text(size=treat_text_size, color="black"),
          strip.text.y=element_text(angle=0, size=9, face="bold")) +
          scale_x_discrete(position="top") +
          labs(x="", y="")

    if(asRatio) pPath = pPath+theme(aspect.ratio = 1)

    if (discrete){
      pPath <- pPath + scale_fill_manual(values=colors,labels=color_leg,na.value = 'gray50')
    }else{
      pPath <- pPath + scale_fill_gradient(low = "darkgreen", high = "red",na.value = 'gray50')
    }

  #afrer faceting we need to open the ggplot object to put different row/column label colors to different facets
  png(file.path(tempdir(),"juen.png"))
  gplotPath <- ggplotGrob( pPath )
  dev.off()
  nms <- lapply( gplotPath$grobs , function(x) names( x[]$children ) )
  #we search for axis.line.x and axis.line.y containers
  grbs_x_id <- which( sapply( lapply( nms , function(x) grepl( "axis.line.x" , x ) ) , any ) == 1 )
  grbs_y_id <- which( sapply( lapply( nms , function(x) grepl( "axis.line.y" , x ) ) , any ) == 1 )
  #grob variable iterates over blocks of colums of the ggplot
  # for (grob in (1:length(grbs_x_id))){
  #   gplotPath$grobs[[grbs_x_id[grob]]]$children$axis$grobs[[1]]$children[[1]]$gp$col=grob
  # }
  
  for (grob in (1:length(grbs_x_id))){
    to_select <- grbs_x_id[grob]
    if (length(to_select) >= 1) {
      item_at_grob <- gplotPath$grobs[[to_select]]
      if (! is.null(item_at_grob)) {
        item_at_grob$children$axis$grobs[[1]]$children[[1]]$gp$col=grob
      }
    }
  }
  

  #Make genes melt
  kegg_melt_gene <- dplyr::filter(kegg_melt, gene_group=="Genes")

  print("str(kegg_melt_gene)")
  print(str(kegg_melt_gene))

  ##Genes map
  p <- ggplot(kegg_melt_gene, aes(Var1, Var2, fill=value)) +
    facet_grid(.~experiment, scales="free", space="free") +
    geom_tile(colour="black", size=0.1) +
    ggtitle(title) +
    theme(plot.title=element_text(hjust=0.5)) +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=0, size=treat_text_size),
          axis.text.y = element_text(size=treat_text_size, color="black"),
          strip.text.y=element_text(angle=0, size=9, face="bold")) +
          scale_x_discrete(position="top") +
          labs(x="", y="")

  ### This commented plot was used for faceted pathway and genes map ###
  # p <- ggplot(kegg_melt, aes(Var1, Var2, fill=value))
  #
  # if(is.null(gene_group)){
  #   p <- p + facet_grid(.~experiment, scales="free", space="free")
  # }else{
  #   p <- p + facet_grid(gene_group~experiment, scales="free", space="free")
  # }
  #
  # p <- p +
  #   geom_tile(colour="black", size=0.1) +
  #   ggtitle(title) +
  #   theme(plot.title=element_text(hjust=0.5)) +
  #   theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=0, size=treat_text_size),
  #         axis.text.y = element_text(size=treat_text_size, color="black"),
  #         strip.text.y=element_text(angle=0, size=9, face="bold")) +
  #         scale_x_discrete(position="top") +
  #         labs(x="", y="")

  if(asRatio) p = p+theme(aspect.ratio = 1)

  if (discrete){
    p <- p + scale_fill_manual(values=colors,labels=color_leg,na.value = 'gray50')
  }else{
    p <- p + scale_fill_gradient(low = "darkgreen", high = "red",na.value = 'gray50')
  }

  #afrer faceting we need to open the ggplot object to put different row/column label colors to different facets
  png(file.path(tempdir(),"juen.png"))
  gplot <- ggplotGrob( p )
  dev.off()
  nms <- lapply( gplot$grobs , function(x) names( x[]$children ) )
  #we search for axis.line.x and axis.line.y containers
  grbs_x_id <- which( sapply( lapply( nms , function(x) grepl( "axis.line.x" , x ) ) , any ) == 1 )
  grbs_y_id <- which( sapply( lapply( nms , function(x) grepl( "axis.line.y" , x ) ) , any ) == 1 )
  #grob variable iterates over blocks of colums of the ggplot
  # for (grob in (1:length(grbs_x_id))){
  #   gplot$grobs[[grbs_x_id[grob]]]$children$axis$grobs[[1]]$children[[1]]$gp$col=grob
  # }

  for (grob in (1:length(grbs_x_id))){
    to_select <- grbs_x_id[grob]
    if (length(to_select) >= 1) {
      item_at_grob <- gplot$grobs[[to_select]]
      if (! is.null(item_at_grob)) {
        item_at_grob$children$axis$grobs[[1]]$children[[1]]$gp$col=grob
      }
    }
  }
  
  # #grob variable iterates over blocks of rows of the ggplot
  # for (grob in (1:length(grbs_y_id))){
  #   gplot$grobs[[grbs_y_id[grob]]]$children$axis$grobs[[1]]$children[[1]]$gp$col=darkcols[grob]
  # }

  #return(p)
  #return(gplot)

  #Combine ggplot for pathway and genes
  #gplotComb <- arrangeGrob(gplotPath, gplot, nrow=2, ncol=1, widths=1)
  gplotComb <- rbind(gplotPath, gplot, size="first")
  print("str(gplotComb)")
  print(str(gplotComb))
  return(gplotComb)
}
