#this script generates a shiny app for the Suillus genome set to visualize / download
#cultures morphologies under different conditions 
#ITS sequences
#growth rates at diferent temperatures


#load libraries
library(shiny)
#library(data.table)
library(tidyverse)
library(shinyBS)
library(ggtree)
library(ape)
library(plotly)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(mailtoR)
library(lubridate) #to work with dates
library(xts) #to work with dates
library(data.table)

##set seed for reproducibility
set.seed(666)

#set wd
setwd("~/Desktop/SuilluScope")

#DB<-read.delim("<X>.txt", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE, check.names = TRUE)

#Metadata
DB<-read.delim("Suillus_app_metadata.txt", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE, check.names = TRUE)

#growth data
GD<-read.delim("MOCK_temperature_assay_2023.csv", header = TRUE, sep = ",", fill = TRUE, strip.white = TRUE, check.names = FALSE)
#fix names 

#remove "contam" rows and NA rows
GD_clean <- GD %>%
  filter(!if_any(everything(), ~str_detect(., "contam")))

#remove "replicate" col
GD_clean$replicate <- NULL

input<- GD_clean
#prep data by reformatting to long form. NOTE unhash if you want to average over the mean)
prep_data_stats<- function(input) {
  #convert to numeric
  cols<- ncol(input)-2
  input[3:cols] <- lapply(input[3:cols], as.numeric)
  keys <- colnames(input)[!grepl('2023',colnames(input))]
  X <- as.data.table(input)
  #input_clean<- X[,lapply(.SD,mean),keys] #YOU ARE HERE CALCULATE SE. 
  
  #make long
  input_clean_long <- X %>%
    pivot_longer(
      cols = -c("Temp", "Strain"),
      names_to = "date",
      values_to = "area",
      values_transform = list(area = as.numeric, Strain = as.character)
    )
  
  #change date format to read as a date in lubridate and xts
  input_clean_long$date2<- lubridate::parse_date_time(input_clean_long$date, "dmy")
  
  #change to number of days 
  days <- yday(input_clean_long$date2) - 165 #first day was Jun 14 = day 0 (the day we inoculated, wich is 165 days into the year)
  input_clean_long$n_days<- days
  
  #clean up unnecessary cols
  input_clean_long$date <- NULL
  input_clean_long$date2 <- NULL
  
  #return output object
  return(input_clean_long)
}

#run function to make input data
prep_stats_test<-prep_data_stats(GD_clean)

#Separate by Temperature treatment
C10<- prep_stats_test[prep_stats_test$Temp == 10,]
C20<- prep_stats_test[prep_stats_test$Temp == 20,]
C24<- prep_stats_test[prep_stats_test$Temp == 24,]
C27<- prep_stats_test[prep_stats_test$Temp == 27,]
C30<- prep_stats_test[prep_stats_test$Temp == 30,]
C34<- prep_stats_test[prep_stats_test$Temp == 34,]
C37<- prep_stats_test[prep_stats_test$Temp == 37,]


#calculate mean and se for each temp
#Function to get the mean and the mean standard error for each group
#+++++++++++++++++++++++++
# data = a data frame
# varname  = the name of a column containing the variable to be summarized
# groupnames = vector of column names to be used as grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  library(plotrix)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      se = std.error(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


#run data_summary function
C10 <- data_summary(C10, varname="area", 
                    groupnames=c("Strain", "n_days"))
C20 <- data_summary(C20, varname="area", 
                    groupnames=c("Strain", "n_days"))
C24 <- data_summary(C24, varname="area", 
                    groupnames=c("Strain", "n_days"))
C27 <- data_summary(C27, varname="area", 
                    groupnames=c("Strain", "n_days"))
C30 <- data_summary(C30, varname="area", 
                    groupnames=c("Strain", "n_days"))
C34 <- data_summary(C34, varname="area", 
                    groupnames=c("Strain", "n_days"))
C37 <- data_summary(C37, varname="area", 
                    groupnames=c("Strain", "n_days"))

datafiles <- list(
  C10 = C10,
  C20 = C20,
  C24 = C24,
  C27 = C27,
  C30 = C30,
  C34 = C34,
  C37 = C37
)

#run stats
#test %>% rstatix::shapiro_test(area)
#if p val non-significant = data has a pretty normal distribution


#Test graph before adding it

#test_plot<- ggplot(C10, aes(x = n_days, y = area, group = Strain, color = Strain)) + 
#  geom_point(aes(color = Strain), size = 2) +
#  geom_line(aes(color = Strain), size = .8) +
#  scale_color_manual(values = c("#003f5c", 
#                                "#668eaa",
#                                "#c2e7ff",
#                                "#2f4b7c",
#                                "#8293bc",
#                                "#d6e2ff",
#                                "#665191", 
#                                "#a794c7",
#                                "#ebdcff",
#                                "#a05195",
#                                "#ce94c4",
#                                "#fcd8f5",
#                                "#d45087",
#                                "#ec95b6",
#                                "#ffd5e5",
#                                "#f95d6a",
#                                "#ff9d9e",
#                                "#ffd6d5",
#                                "#ff7c43",
#                               "#ffac82",
#                                "#ffd9c6",
#                                "#ffa600",
#                                "#ffc171",
#                                "#fbddbe",
#                                "#962B09")) +
#  theme_minimal()+
#  xlab("dpi")+ ylab("colony area (mm^2)")+
#  ylab(bquote("Colony area "~(mm^2)))+
#  ggtitle("S. hirtellus EM16")+
#  geom_errorbar(aes(ymin=area-se, ymax=area+se, alpha=0.6), width=2, show.legend = FALSE,
#                position=position_dodge(1)) 



##############################

##-----------create ui object-----------
ui <- fluidPage(
  fluidRow( 
    column(3, img(src='Suilluscope_logo.png', align = "left", height="42%", width="42%", res = 128), style = "padding-top:6px; padding-bottom:6px; padding-left:25px; padding-right:0px; margin-right:-75px"),
    column(9, tags$head(tags$style(HTML("
  .tabbable > .nav {
    padding-top: 60px;
  }
  .tabbable > .nav > li[class=active] > a {
    color: white;
    background-color: #666666;
    width: 110px;
    height: 35px;
    border-radius: 0px;
    font-family: Arial;
    font-size: 15px;
    padding-bottom: 6px;
    display: flex;
    align-items: center;
    justify-content: center;
  }
  .tabbable > .nav > li > a {
    color: black;
    background-color: white;
    width: 110px;
    height: 35px;
    border-radius: 0px;
    font-family: Arial;
    font-size: 15px;
    padding-bottom: 6px;
    display: flex;
    align-items: center;
    justify-content: center;
    white-space: nowrap;
  }
  .tabbable > .nav > li > a:hover {
    color: black;
    background-color: #E6E6E6;
    width: 110px;
    height: 35px;
    border-radius: 0px;
    font-family: Arial;
    font-size: 15px;
    padding-bottom: 6px;
    display: flex;
    align-items: center;
    justify-content: center;
  }
  .tabbable > .nav > li {
    padding-right: 48px;
  }
")))
    ),

    
    tabsetPanel(type = "tabs",
#                tabPanel("HOME", fluidRow(
#                  column(12, img(src='App_home.png', align = "left", height="42%", width="42%", res = 128), style = "padding-top:6px; padding-bottom:6px; padding-left:25px; padding-right:0px; margin-right:-75px")), htmlOutput(outputId ="home", style = "background-color:transparent; padding:20px; margin-top:10px")),
tags$head(
  tags$style(
    HTML("
        #home_bottom2 a, #home_right a {
          color: #DDC87B;
          text-decoration: none;
        }
        .background-image-container {
          background-size: cover;
          background-position: center;
          background-repeat: no-repeat;
          position: relative;
          margin-top: -15px; 
        }
      ")
  )
),

tabPanel("HOME",
         tags$div(
           class = "background-image-container",
           style = "background-image: url('App_home.png'); height: 300px;",
           fluidRow(
             column(
               width = 6, 
               htmlOutput(outputId ="home_left", style = "background-color: rgba(0, 0, 0, 0.5); color: white; padding:20px; margin-top:54px; margin-bottom:35px; margin-right:355px; margin-left:-282px;")
             ),
             column(
               width = 3, 
               htmlOutput(outputId ="home_right", style = "background-color: rgba(0, 0, 0, 0.5); color: white; padding:20px; margin-top:54px; margin-bottom:35px; margin-left:-270px; margin-right:-36px;")
             )
           ),
           fluidRow(
             column(
               width = 12,
               div(
                 style = "position: relative;", 
                 img(src='lower_background.png', align = "center", height="60%", width="101%", res = 128),
                 htmlOutput(outputId ="home_bottom1", style = "position: absolute; top: 10px; left: 10px; color: #4D4D4D; margin: 0; margin-top:5px; margin-left:450px; margin-right:450px; padding: 15px; background-color: rgba(0, 0, 0, .5); background-color: rgba(0, 0, 0, 0); border: 5px solid gray; border-radius: 20px")
               )
             )
           ),
           fluidRow(
             column(
               width = 12,
               div(
                 style = "position: relative;", 
                 htmlOutput(outputId ="home_bottom2", style = "position: absolute; top: 10px; left: 10px; color: black; margin: 0; margin-top:-150px; margin-left:30px; padding: 0")
               )
             )
           )
         )
),           
                tabPanel(
                  "MORPHOLOGY",
                  fluidRow(

                  ),
                  fluidRow(
                    column( 
                      1, 
                      style = "height:700px; width:50px; background-color:rgba(119,120,55,0.5); padding-left: 0px; padding-right: 20px"), #little buffer col.
                    column( 
                      2, 
                      style = "height:700px; width:250px; background-color:rgba(119,120,55,0.5); padding-left: 0px; padding-right: 20px",
                      br(), 
                      #offset = 1,
                      selectizeInput(
                        inputId = "strain_image",
                        label = "SELECT STRAIN",
                        multiple = FALSE,
                        choices = c("", DB$Full.Name),
                        options = list(
                          openOnFocus = FALSE,
                          maxItems = 1,
                          render = I(
                            '{
            option: function(item, escape) {
              if (item.label === "") {
                return "<div data-value=\'" + escape(item.value) + "\'></div>";
              } else {
                var words = item.label.split(" ");
                var label = "<div>" + words.slice(0, -1).map(function(word) { return "<em>" + escape(word) + "</em>"; }).join(" ") + " " + escape(words.slice(-1)[0]) + "</div>";
                return "<div data-value=\'" + escape(item.value) + "\'>" + label + "</div>";
              }
            },
            item: function(item, escape) {
              if (item.label === "") {
                return "<div data-value=\'" + escape(item.value) + "\'></div>";
              } else {
                var words = item.label.split(" ");
                var label = "<div>" + words.slice(0, -1).map(function(word) { return "<em>" + escape(word) + "</em>"; }).join(" ") + " " + escape(words.slice(-1)[0]) + "</div>";
                return "<div data-value=\'" + escape(item.value) + "\'>" + label + "</div>";
              }
            }
          }'
                          ),
                          onChange = I(
                            'function(value) {
            if (value !== "") {
              var lastWord = value.split(" ").pop();
              var imageName = lastWord + ".png";
              var imagePath = imageName;
              var img = new Image();
              img.src = imagePath;
              img.setAttribute("width", "150%");
              var targetDiv = document.getElementById("plate_img");
              targetDiv.innerHTML = "";
              targetDiv.appendChild(img);
            } else {
              var targetDiv = document.getElementById("plate_img");
              targetDiv.innerHTML = "";
            }
          }'
                          )
                        )
                      ),
                      #radioButtons(inputId = "switch_input", "Choose a dataset", choices = c("top of plate", "bottom of plate")),
                      #style = "padding-right:40px"
                    ),
                    column(
                      9, 
                      div(id = "plate_img"),
                      plotlyOutput("plate_img", height = "715px", width = "165px"), 
                      style = "background-color:transparent; padding-left:100px; padding-right:400px; padding-top:0px; padding-bottom:200px; margin-top:25px" #note margin-top changes how close the .png is to the top of the row
                    )
                  ),
                  br(),
                  htmlOutput(outputId = "morphology", width = "340px", style = "background-color:#E6E6E6; padding:18px")
                ),

                tabPanel("GROWTH RATE", fluidRow(
                  column(1, 
                         style = "height:700px; width:50px; background-color:rgba(119,120,55,0.5); padding-left: 0px; padding-right: 20px"), # little buffer col
                  column(2, 
                         style = "height:700px; width:250px; background-color:rgba(119,120,55,0.5); padding-left: 0px; padding-right: 20px",
                         br(),
                         selectInput(inputId = "datafiles",
                                     label = "SELECT TEMPERATURE (°C):",
                                     choices = c("C10" = "10", "C20" = "20", "C24" = "24", "C27" = "27", "C30" = "30", "C34" = "34", "C37" = "37"))),
                  column(8, 
                         plotlyOutput("growth_plot", height = "600px", width = "700px"),
                         htmlOutput(outputId = "growth", style = "background-color:transparent; padding-right:5px; padding-left:10px; margin-top:10px"),
                         offset = 1),  # Adjust the offset value as needed
                  
                  )), 
                tabPanel("METADATA", fluidRow(
                  column(12, style = "height:2px; background-color:#000000")), htmlOutput(outputId ="metadata", style = "background-color:transparent; padding:20px; margin-top:10px")),
                
                
                )
    )
  )



##-----------create server object-----------
######################

server<- function(input, output){ 
  #define 'HOME' page
  output$home_left <- renderText(paste(tags$b("WELCOME TO", style = "font-size: 24px;"),
                                       tags$b("SuilluScope v1.0", style = "font-size: 24px; color: #D3AD0D"),
                                       tags$br(),
                                       "SuilluScope is an interactive, open access database for <i>Suillus</i> fungi, a model genus for ectomycorrhizal ecology.",
                                       "There are extensive genomic resources available for", tags$i("Suillus"), "including more annotated full-genome assemblies 
                                       than any other ectomycorrhizal group.",
                                       "This database characterizes the phenotypic traits and responses of these genome-sequenced isolates, aiming to help researchers identify optimal 
                                       culture conditions, predict and compare trait responses across diverse species within the genus, and empower further 
                                       research linking genotypes to phenotypes.",
                                       tags$br()
  ))
  
  
  output$home_right<- renderText(paste(tags$b("ISOLATE COLLECTION", style = "font-size: 24px; color: #D3AD0D"),
                                       tags$br(),
                                       "All isolates used in the database are from the", tags$i("Suillus"), "genome strain collection, curated at Duke University.", 
                                       "For more information on the <i>Suillus</i> system, please visit",
                                       tags$a(href = "http://www2.hawaii.edu/~nn33/suillus/", "The International ", tags$i("Suillus"), " Consortium.", target = "_blank"),
                                       "Genomic resources for these strains are available on the MycoCosm ",tags$i("Suillus"), 
                                       tags$a(href="https://mycocosm.jgi.doe.gov/Suillus/Suillus.info.html", "web portal", target="_blank"), 
                                       "opperated by DOE Joint Genome Institute.",
                                       "Isolates amenable to cryopreservation are publicly available as part of the", 
                                       tags$a(href="https://nrrl.ncaur.usda.gov/", "Agricultural Research Culture Collection (NRRL).", target="_blank"), 
                                       "Isolates which cannot be cryopreserved are available by contacting us", mailtoR(email = "LotusLofgren@gmail.com", text = "here.")
  ))
  
  output$home_bottom1  <- renderText(paste(tags$b("HOW TO CITE", style = "font-size: 24px; color: #666666"), tags$br(), 
                                           "If you find SuilluScope useful in your own research, please cite the assocaited publication: <>", 
                                           tags$br(), 
                                           "If you use data on a specific genome strain in your own research, please cite SuilluScope as well as the paper that originally published that strain (listed in the Metadata tab)." 
  )) 
  output$home_bottom2 <- renderText(paste("This is version v1.0 of the database, released on 25.July.2023",
                                         tags$br(),
                                         "SuilluScope was built using the open source programming language", 
                                         tags$a(href="https://www.r-project.org/about.html", "R", target="_blank"), 
                                         "with reactive programming via", 
                                         tags$a(href="https://shiny.rstudio.com/", "R shiny", target ="_blank"),
                                         tags$br(),
                                         "All of the code necessary to run the program is publicly available at", tags$a(href="https://github.com/MycoPunk/SUILLUSAPP", "the SuilluScope GitHub.", target= "_blank"), 
                                         "Please report issues to the git issues page",
                                         tags$br(),
                                         "If you have feature requests or a published dataset that you would like us to consider adding to this site, please contact the author", mailtoR(email = "LotusLofgren@gmail.com", text = "here.")

  ))
  
  
  
#define 'GROWTH' page displaying growth at different temperatures
# use datafiles list of DFs to make the DF choice reactive (STATIC)
  
#  datasetInput <- reactive({
#    switch(input$datafiles,
#           "10" = C10,
#           "20" = C20,
#           "24" = C24, 
#           "27" = C27,
#           "30" = C30,
#           "34" = C34, 
#           "37" = C37)
#  })
  
  
  
  # #plot growth data (STATIC)
  # output$growth_plot<- renderPlotly({
  #   ggplot(datasetInput(), aes(x = n_days, y = area, group = Strain, color = Strain)) + 
  #     geom_point(aes(color = Strain), size = 2) +
  #     geom_line(aes(color = Strain), size = .8) +
  #     scale_color_manual(values = c("#003f5c", 
  #                                   "#668eaa",
  #                                   "#c2e7ff",
  #                                   "#2f4b7c",
  #                                   "#8293bc",
  #                                   "#d6e2ff",
  #                                   "#665191", 
  #                                   "#a794c7",
  #                                   "#ebdcff",
  #                                   "#a05195",
  #                                   "#ce94c4",
  #                                   "#fcd8f5",
  #                                   "#d45087",
  #                                   "#ec95b6",
  #                                   "#ffd5e5",
  #                                   "#f95d6a",
  #                                   "#ff9d9e",
  #                                   "#ffd6d5",
  #                                   "#ff7c43",
  #                                   "#ffac82",
  #                                   "#ffd9c6",
  #                                   "#ffa600",
  #                                   "#ffc171",
  #                                   "#fbddbe",
  #                                   "#962B09")) +
  #     theme_minimal()+
  #     xlab("dpi")+ ylab("colony area (mm^2)")+
  #     ylab(bquote("Colony area "~(mm^2)))+
  #     #ggtitle("S. <> <>")+
  #     geom_errorbar(aes(ymin=area-se, ymax=area+se, alpha=0.6), width=2, show.legend = FALSE,
  #                   position=position_dodge(1)
  #     )})}


  

#use datafiles list of DFs to make the DF choice reactive and filter the data by temperature
    get_filtered_data <- function(temp) {
      switch(temp,
             "10" = C10,
             "20" = C20,
             "24" = C24,
             "27" = C27,
             "30" = C30,
             "34" = C34,
             "37" = C37)
    }
    
    
#plot that shit
    output$growth_plot <- renderPlotly({
      dataset <- get_filtered_data(input$datafiles)
      
      plot_ly(
        data = dataset,
        x = ~n_days,
        y = ~area,
        color = ~Strain,
        colors = c("#003f5c", "#668eaa", "#c2e7ff", "#2f4b7c", "#8293bc", "#d6e2ff", "#665191", "#a794c7", "#ebdcff",
                   "#a05195", "#ce94c4", "#fcd8f5", "#d45087", "#ec95b6", "#ffd5e5", "#f95d6a", "#ff9d9e", "#ffd6d5",
                   "#ff7c43", "#ffac82", "#ffd9c6", "#ffa600", "#ffc171", "#fbddbe", "#962B09")) %>%
        add_lines(size = 1, line = list(width = 2), hoverinfo = "text", text = ~Strain) %>%
        layout(
          title = "",
          xaxis = list(title = "days post inoculation"),
          yaxis = list(title = HTML("colony area (mm<sup>2</sup>)")),
          showlegend = TRUE,
          legend = list(title = "Strain"),
          hovermode = "closest"
        )
    })
    
#click the plot to update line visibility
    observeEvent(event_data("plotly_click", source = "growth_plot"), {
      selected_strain <- event_data("plotly_click", source = "growth_plot")$y
      dataset <- get_filtered_data(input$datafiles)
      
#change visibility for the selected strain
      if (is.null(selected_strain)) {
        visible <- TRUE
      } else {
        visible <- ifelse(dataset$Strain %in% selected_strain, FALSE, TRUE)
      }
      
      plotlyProxy("growth_plot") %>%
        plotlyProxyInvoke("restyle", list(visible = visible), list(2)) %>%
        plotlyProxyInvoke("restyle", list(visible = visible), list(3))
    })
  

}




###run as app
shinyApp(ui = ui, server = server)
