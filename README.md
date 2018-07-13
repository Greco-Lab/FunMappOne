# FunMappOne

### Linux system library dependencies

```BASH
libv8-3.14-dev
libxml2-dev 
libssl-dev
```

### Install dependencies

```R
#Install Bioconductor dependencies
source("http://bioconductor.org/biocLite.R")
bioc_pkgs <- c("org.Hs.eg.db", "org.Mm.eg.db", "KEGG.db", "reactome.db", "GOSim")
bioc_pkgs.inst <- bioc_pkgs[!(bioc_pkgs %in% rownames(installed.packages()))]
if(length(bioc_pkgs.inst)>0){
  print(paste0("Missing ", length(bioc_pkgs.inst), " Bioconductor Packages:"))
  for(pkg in bioc_pkgs.inst){
    print(paste0("Installing Package:'", pkg, "'..."))
    biocLite(pkg, suppressUpdates=TRUE)
    print("Installed!!!")
  }
}

#Install CRAN dependencies
cran_pkgs <- c("ggplotify", "RColorBrewer", "reshape", "ggplot2", "shiny", "shinyjs", "tibble", 
               "gProfileR", "DT", "randomcoloR", "readxl", "cellranger", "devtools", "scales", "gtools")
cran_pkgs.inst <- cran_pkgs[!(cran_pkgs %in% rownames(installed.packages()))]
if(length(cran_pkgs.inst)>0){
  print(paste0("Missing ", length(cran_pkgs.inst), " CRAN Packages:"))
  for(pkg in cran_pkgs.inst){
    print(paste0("Installing Package:'", pkg, "'..."))
    install.packages(pkg, repo="http://cran.rstudio.org", dependencies=TRUE)
    print("Installed!!!")
  }
}
```
### Run FunMappOne from GitHub
```R
# Load 'shiny' library
library(shiny)
library(shinyjs)
# Using runGitHub
runGitHub("FunMappOne", "Greco-Lab")
```

#### How to run locally
```R
  # Clone the git repository
  git clone https://github.com/Greco-Lab/FunMappOne FunMappOne

  # Start R session and run by using runApp()
  library(shiny)
  runApp("FunMappOne/")
```
