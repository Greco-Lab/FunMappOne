library(shiny)
library(shinyjs)

runApp(
	appDir="/home/docker/shiny_app", 
	port=3838, 
	launch.browser=FALSE, 
	host="0.0.0.0"
)
