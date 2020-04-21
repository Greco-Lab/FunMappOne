FROM rocker/shiny:3.5.3

RUN apt-get update && apt-get install -y \
	openjdk-8-jdk \
	libcurl4-openssl-dev \
	libxml2-dev \
	libv8-dev \
	libssl-dev \
	liblzma-dev \
	libbz2-dev 

RUN ["R", "CMD", "javareconf"]

RUN useradd -m docker
USER docker
RUN mkdir -p "/home/docker/.R/libs"
RUN echo '.libPaths("/home/docker/.R/libs")' >> /home/docker/.Rprofile
COPY install_dependencies.R /home/docker/install_dependencies.R
RUN ["Rscript", "/home/docker/install_dependencies.R"]
RUN mkdir /home/docker/shiny_app
COPY ./*.R /home/docker/shiny_app/
COPY ./*.RData /home/docker/shiny_app/
COPY ./www /home/docker/shiny_app/

COPY ./Rprofile.site /usr/lib/R/etc/

EXPOSE 3838

ENTRYPOINT ["Rscript", "/home/docker/shiny_app/shiny_start.R"]
