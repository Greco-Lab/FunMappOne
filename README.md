# FunMappOne

## Running the FunMappOne Docker image

A FunMappOne docker image is available at https://hub.docker.com/r/grecolab/funmappone/

All you need to do is to downalad the image with docker "docker pull grecolab/funmappone"
and run it by mapping the Docker http port 3838 on the host port 8787.

## Installing FunMappOne your local system 

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

### Run FunMappOne 
#### From GitHub
```R
# Load 'shiny' library
library(shiny)
library(shinyjs)
# Using runGitHub
runGitHub("FunMappOne", "Greco-Lab")
```

#### From a local copy 
```R
  # Clone the git repository
  git clone https://github.com/Greco-Lab/FunMappOne FunMappOne
  # Start R session and run by using runApp()
  library(shiny)
  library(shinyjs)
  runApp("FunMappOne/")
```

## Using FunMappOne 

Once your FunMappOne instance is istalled and runnig just open your browser and visit:
http://localhost:8787/FunMappOne

If everything has been installed correctly the FunMappOne welcome page should be now available.
