library(GO.db)
library(GOSim)
library(igraph)

# res_x$Ontology = as.factor(vapply(res_x$GO_ID,function(x){ if (x %in% keys(GOTERM)) res=GOTERM[[x]]@Ontology else res=""; return(res) },c("")))
# 
# 
# res_x$Ontology = as.factor(vapply(res_x$GO_ID,function(x){ if (x %in% keys(GOTERM)) res=GOTERM[[x]]@Ontology else res=""; return(res) },c("")))
# 
# term_bp=unique(res_x$GO_ID[res_x$Ontology == "BP"])
# 


createGoHier <- function(go_terms,type="BP") {
  
  if(type=="BP"){
    root="GO:0008150"
    anc_fun=GOBPANCESTOR
    child_fun=GOBPCHILDREN
    
  }else if (type=="CC"){
    root="GO:0005575"
    anc_fun=GOCCANCESTOR
    child_fun=GOCCCHILDREN
  }else if(type=="MF"){
    root="GO:0003674"
    anc_fun=GOMFANCESTOR
    child_fun=GOMFCHILDREN
  }
  
  a=tibble(level1="",level2="",level3="")
  a=a[-1,]
  
  ancestors=unique(c(unlist(as.list(anc_fun[go_terms])),go_terms))
  
  bp1=intersect(unlist(as.list(child_fun[root])),ancestors)
  pb = txtProgressBar(min =1,max = length(go_terms),style = 3)
  
  for(i in 1:length(go_terms)){
    lev1_anc = intersect(unlist(as.list(anc_fun[go_terms[i]])),bp1)
    if(length(lev1_anc)>0)
      for (j in 1:length(lev1_anc)){
        lev2_anc = unlist(as.list(child_fun[lev1_anc[j]]))
        lev2_anc = lev2_anc[grep("is_a|part_of",names(lev2_anc))]
        lev2_anc=intersect(lev2_anc,ancestors)
        for (l in 1:length(lev2_anc)){}
         # a=add_row(a,level1=lev1_anc[j],level2=lev2_anc[l],level3=go_terms[i])
      }
    setTxtProgressBar(pb,i)
  }
  close(pb)
  return(a)
}

# hier=createGoHier(go_terms = term_bp, type = "BP")
# 
# hier_names=hier
# hier_names$level1 = vapply( hier$level1, function(x){ if (x %in% keys(GOTERM)) res=GOTERM[[x]]@Term else res=""; return(res) },c(""))
# hier_names$level2 = vapply( hier$level2, function(x){ if (x %in% keys(GOTERM)) res=GOTERM[[x]]@Term else res=""; return(res) },c(""))
# hier_names$level3 = vapply( hier$level3, function(x){ if (x %in% keys(GOTERM)) res=GOTERM[[x]]@Term else res=""; return(res) },c(""))

