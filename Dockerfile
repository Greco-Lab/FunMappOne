FROM rocker/r-base:latest

RUN apt-get update && apt-get install -y \
	openjdk-8-jdk \
	libcurl4-openssl-dev \
	libxml2-dev \
	libv8-dev 

RUN ["R", "CMD", "javareconf"]

USER docker
COPY install_dependencies.R /home/docker/install_dependencies.R
RUN ["Rscript", "/home/docker/install_dependencies.R"]

ENTRYPOINT ["Rscript", "/home/docker/shiny_app/shiny_start.R"]
