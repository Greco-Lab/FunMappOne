#Universal Bioconductor package installation function
install.bioc <- function(pkg){
	vers <- getRversion()
	if (vers >= "3.6"){
		if (!requireNamespace("BiocManager", quietly = TRUE)) {
			install.packages("BiocManager")
		}
		BiocManager::install(pkg)
	}else{
		if (!requireNamespace("BiocInstaller", quietly = TRUE)){
			source("https://bioconductor.org/biocLite.R")
			biocLite(pkg, suppressUpdates=TRUE)
		}else{
			BiocInstaller::biocLite(pkg, suppressUpdates=TRUE)
		}
	}
}

#Install CRAN dependencies
cran_pkgs <- c(
	"DT", "RColorBrewer", "V8", "cellranger", "devtools", "gProfileR", "ggplot2",
	"ggplotify", "grid", "gridExtra", "gtable", "gtools", "igraph", "randomcoloR",
	"readxl", "reshape", "scales", "shiny", "shinyBS", "shinycssloaders",
	"shinyjs", "tibble", "tidyverse", "tidyverse", "xlsx"
)
cran_pkgs.inst <- cran_pkgs[!(cran_pkgs %in% rownames(installed.packages()))]
if(length(cran_pkgs.inst)>0){
	print(paste0("Missing ", length(cran_pkgs.inst), " CRAN Packages:"))
	for(pkg in cran_pkgs.inst){
		print(paste0("Installing Package:'", pkg, "'..."))
		install.packages(
			pkg, repo="http://cran.rstudio.org", dependencies=TRUE
		)
		print("Installed!!!")
	}
}

#Install Bioconductor dependencies
bioc_pkgs <- c(
	"GOSim", "KEGG.db", "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db",
	"reactome.db"
)
bioc_pkgs.inst <- bioc_pkgs[!(bioc_pkgs %in% rownames(installed.packages()))]
if(length(bioc_pkgs.inst)>0){
	print(paste0("Missing ", length(bioc_pkgs.inst), " Bioconductor Packages:"))
	for(pkg in bioc_pkgs.inst){
		print(paste0("Installing Package:'", pkg, "'..."))
		install.bioc(pkg)
		print("Installed!!!")
	}
}

