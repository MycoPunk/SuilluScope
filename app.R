####
#This script holds all of the code necessary to run the R shiny app, SuilluScope. 
# this is version 1.0, released in August of 2023. 
####

#load libraries
library(shiny)
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
library(DT) #for working with interactive datatables

##set seed for reproducibility
set.seed(666)

#set wd
setwd("~/Desktop/SuilluScope")

#Metadata
DB<-read.delim("Suillus_app_metadata.txt", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE, check.names = TRUE)

#metadata for the app
DB_slim<- DB
DB_slim$Assembly.Length<- NULL
DB_slim$n.Genes<- NULL
DB_slim$Sequencing.depth<- NULL 
DB_slim$n.Gaps<- NULL
DB_slim$Scaffold.N50<- NULL
DB_slim$Scaffold.L50.Mbp<- NULL
DB_slim$n.contigs<- NULL
DB_slim$n.scaffolds<- NULL
DB_slim$seq.platform<- NULL
DB_slim$Host.speices.for.assays<- NULL
DB_slim$Range_and_region<- NULL

#growth data
#GD<-read.delim("MOCK_temperature_assay_2023.csv", header = TRUE, sep = ",", fill = TRUE, strip.white = TRUE, check.names = FALSE)
GD<-read.delim("temperature_assay_2023_20s.csv", header = TRUE, sep = ",", fill = TRUE, strip.white = TRUE, check.names = FALSE)

#fix names (from strain codes to the full species names + strain codes)
lookup_table <- DB %>%
  select(Culture.code, Full.Name)
GD$Culture.code <- lookup_table$Full.Name[match(GD$Culture.code, lookup_table$Culture.code)]

#remove "contam" rows and NA rows
GD_clean <- GD %>%
  filter(!if_any(everything(), ~str_detect(., "contam")))

#make a copy for the growth rate data
GD_clean_with_rep<- GD_clean

#remove "replicate" col for parser
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
      cols = -c("Temp", "Culture.code"),
      names_to = "date",
      values_to = "area",
      values_transform = list(area = as.numeric, Culture.code = as.character)
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
                    groupnames=c("Culture.code", "n_days"))
C20 <- data_summary(C20, varname="area", 
                    groupnames=c("Culture.code", "n_days"))
C24 <- data_summary(C24, varname="area", 
                    groupnames=c("Culture.code", "n_days"))
C27 <- data_summary(C27, varname="area", 
                    groupnames=c("Culture.code", "n_days"))
C30 <- data_summary(C30, varname="area", 
                    groupnames=c("Culture.code", "n_days"))
C34 <- data_summary(C34, varname="area", 
                    groupnames=c("Culture.code", "n_days"))
C37 <- data_summary(C37, varname="area", 
                    groupnames=c("Culture.code", "n_days"))

datafiles <- list(
  C10 = C10,
  C20 = C20,
  C24 = C24,
  C27 = C27,
  C30 = C30,
  C34 = C34,
  C37 = C37
)


###prep growth rate data for 20 degree control
#get only 20 degree
C20_all_data<- GD_clean_with_rep[GD_clean_with_rep$Temp == 20,]

#make longer (keep rep codes)
#make long
input_clean_long_reps <- C20_all_data %>%
  pivot_longer(
    cols = -c("Temp", "Culture.code", "replicate"),
    names_to = "date",
    values_to = "area",
    values_transform = list(area = as.numeric, Culture.code = as.character)
  )

#make data.table
input_clean_long_reps<- as.data.table(input_clean_long_reps)

#change date format to read as a date in lubridate and xts
input_clean_long_reps$date2<- lubridate::parse_date_time(input_clean_long_reps$date, "dmy")

#change to number of days 
days <- yday(input_clean_long_reps$date2) - 165 #first day was Jun 14 = day 0 (the day we inoculated, wich is 165 days into the year)
input_clean_long_reps$n_days<- days

#clean up unnecessary cols
input_clean_long_reps$date <- NULL
input_clean_long_reps$date2 <- NULL
input_clean_long_reps$Temp <- NULL


#calculate growth rate, group by each unique combination of Sp and replicate, sort by date for each group, and use lag to calculate the difference in area between that date and the one before it.
growth_rates <- input_clean_long_reps %>%
  group_by(Culture.code, replicate) %>%
  group_modify(~mutate(., increase_in_area = area - lag(area, default = first(area))))

#Remove the grouping
growth_rates <- ungroup(growth_rates)

#get the average increase_in_area and se for each Sp and n_days combination
df_average <- growth_rates %>%
  group_by(Culture.code, n_days) %>%
  do(data.frame(avg_increase = mean(.$increase_in_area),
                std_error = sd(.$increase_in_area) / sqrt(length(.$increase_in_area)))) %>%
  ungroup()







#run stats
#test %>% rstatix::shapiro_test(area)
#if p val non-significant = data has a pretty normal distribution


#Test graph before adding it

#test_plot<- ggplot(C10, aes(x = n_days, y = area, group = Culture.code, color = Culture.code)) + 
#  geom_point(aes(color = Culture.code), size = 2) +
#  geom_line(aes(color = Culture.code), size = .8) +
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
    width: 130px;
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
    width: 130px;
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
    width: 130px;
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
                      selectizeInput(
                        inputId = "Culture.code_image",
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
                      htmlOutput(outputId ="textbox_morpholgy", style = "background-color: rgba(0, 0, 0, 0.5); color: white; padding:5px; margin-top:200px; margin-bottom:0px; margin-left:0px; margin-right:0px;")
                      
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

                tabPanel("TEMPERATURE", 
                         fluidRow(
                  column(1, 
                         style = "height:700px; width:50px; background-color:rgba(119,120,55,0.5); padding-left: 0px; padding-right: 20px"), # little buffer col
                  column(2, 
                         style = "height:700px; width:250px; background-color:rgba(119,120,55,0.5); padding-left: 0px; padding-right: 20px",
                         br(),
                         selectInput(inputId = "datafiles",
                                     label = "SELECT TEMPERATURE:",
                                     choices = c("10°C" = "10", "20°C" = "20", "24°C" = "24", "27°C" = "27", "30°C" = "30", "34°C" = "34", "37°C" = "37")),
                         htmlOutput(outputId ="textbox_temperature", style = "background-color: rgba(0, 0, 0, 0.5); color: white; padding:5px; margin-top:200px; margin-bottom:0px; margin-left:0px; margin-right:0px;")
                  ),
                  column(8, 
                         plotlyOutput("temp_plot", height = "600px", width = "700px"),
                         htmlOutput(outputId = "growth", style = "background-color:transparent; padding-right:5px; padding-left:10px; margin-top:10px"),
                         offset = 1),
                  
                  )),

tabPanel("GROWTH RATES", 
         fluidRow(
           column(11, 
                  plotlyOutput("growth_rate_plot", height = "600px", width = "700px"),
                  htmlOutput(outputId = "growth_rate", style = "background-color:transparent; padding-right:5px; padding-left:10px; margin-top:10px"),
                  offset = 1),
         )),
                 
                 
tabPanel(
  "METADATA",
  fluidRow(
    column(12, style = "height:2px; background-color:#000000")
  ),
  htmlOutput(outputId ="metadata", style = "background-color:transparent; padding:20px; margin-top:10px"),
  titlePanel("Interactive Spreadsheet"),
  sidebarLayout(
    sidebarPanel(
      # Checkbox group for column selection
      checkboxGroupInput(
        "columns",
        "Select Columns:",
        choices = colnames(DB_slim),
        selected = colnames(DB_slim)  # By default, all columns are selected
      )
    ),
    mainPanel(
      # Data table output
      DTOutput("data_table")  # Add this line to display the data table
    )
  )
),


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
                                       "All isolates used in the database are from the", tags$i("Suillus"), "genome Culture.code collection, curated at Duke University.", 
                                       "For more information on the <i>Suillus</i> system, please visit",
                                       tags$a(href = "http://www2.hawaii.edu/~nn33/suillus/", "The International ", tags$i("Suillus"), " Consortium.", target = "_blank"),
                                       "Genomic resources for these Culture.codes are available on the MycoCosm ",tags$i("Suillus"), 
                                       tags$a(href="https://mycocosm.jgi.doe.gov/Suillus/Suillus.info.html", "web portal", target="_blank"), 
                                       "opperated by DOE Joint Genome Institute.",
                                       "Isolates amenable to cryopreservation are publicly available as part of the", 
                                       tags$a(href="https://nrrl.ncaur.usda.gov/", "Agricultural Research Culture Collection (NRRL).", target="_blank"), 
                                       "Isolates which cannot be cryopreserved are available by contacting us", mailtoR(email = "LotusLofgren@gmail.com", text = "here.")
  ))
  
  output$home_bottom1  <- renderText(paste(tags$b("HOW TO CITE", style = "font-size: 24px; color: #666666"), tags$br(), 
                                           "If you find SuilluScope useful in your own research, please cite the assocaited publication: <>", 
                                           tags$br(), 
                                           "If you use data on a specific genome Culture.code in your own research, please cite SuilluScope as well as the paper that originally published that Culture.code (listed in the Metadata tab)." 
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
  
  
  output$textbox_morpholgy <- renderText(paste("This is a test"
  ))
  
  output$textbox_temperature <- renderText(paste("This is a test"
  ))
  


  

### 
#TEMPERATURE PLOT
###
  
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
    output$temp_plot <- renderPlotly({
      dataset <- get_filtered_data(input$datafiles)
      
#Function to format the legend text with italics for genus / specific ep.
      format_legend_text <- function(culture_code) {
        words <- strsplit(culture_code, " ")[[1]]
        if (length(words) >= 3) {
          return(paste0("<i>", words[1], " ", words[2], "</i>", " ", paste(words[-c(1, 2)], collapse = " ")))
        } else if (length(words) == 2) {
          return(paste0("<i>", words[1], " ", words[2], "</i>"))
        } else {
          return(culture_code)
        }
      }
      
#Create a new column in the with the formatted legend text
      dataset$Formatted_Legend <- sapply(dataset$Culture.code, format_legend_text)
   
#Define the color palette
species_colors <- c("#003f5c", "#668eaa", "#c2e7ff", "#2f4b7c", "#8293bc", "#d6e2ff", "#665191",
                          "#a794c7", "#ebdcff", "#a05195", "#ce94c4", "#fcd8f5", "#d45087", "#ec95b6",
                          "#ffd5e5", "#f95d6a", "#ff9d9e", "#ffd6d5", "#ff7c43", "#ffac82", "#ffd9c6",
                          "#ffa600", "#ffc171", "#fbddbe", "#962B09")
      
         
#Plot
      plot_ly(data = dataset, x = ~n_days, y = ~area, line = list(width = 2),
              legendgroup = ~Formatted_Legend, name = ~Formatted_Legend,
              color = ~Formatted_Legend, colors = species_colors,
              hoverinfo = "text", text = ~Culture.code,
              showlegend = TRUE) %>%
        layout(
          title = "",
          xaxis = list(title = "days post inoculation"),
          yaxis = list(title = HTML("colony area (mm<sup>2</sup>)")),
          legend = list(title = "Culture.code"),
          hovermode = "closest"
        )
    })
    
    

#click the plot to update line visibility
    observeEvent(event_data("plotly_click", source = "temp_plot"), {
      selected_Culture.code <- event_data("plotly_click", source = "temp_plot")$y
      dataset <- get_filtered_data(input$datafiles)
      
#change visibility for the selected Culture.code
      if (is.null(selected_Culture.code)) {
        visible <- TRUE
      } else {
        visible <- ifelse(dataset$Culture.code %in% selected_Culture.code, FALSE, TRUE)
      }
      
      plotlyProxy("temp_plot") %>%
        plotlyProxyInvoke("restyle", list(visible = visible), list(2)) %>%
        plotlyProxyInvoke("restyle", list(visible = visible), list(3))
    })
  


#
#output$growth_rate_plot <- renderPlot({
#  ggplot(df_average, aes(x = as.factor(n_days), y = avg_increase, color = Culture.code)) +
#    geom_point(position = position_jitter(width = 0.2), size = 3) +
#    geom_errorbar(aes(ymin = avg_increase - std_error, ymax = avg_increase + std_error),
#                  width = 0.2, position = position_dodge(width = 0.2)) +
#    geom_smooth(method = "loess", se = FALSE, size = 1, aes(group = Culture.code, color = Culture.code)) +
#    labs(x = "n_days", y = "Average Increase in Area", title = "Average Increase in Area by n_days and Sp") +
#    theme_minimal() +
#    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better readability
#    scale_color_manual(values = species_colors)
#})


###
#GROWTH RATE PLOT
###
    
    output$growth_rate_plot <- renderPlotly({
      species_colors <- c("#003f5c", "#668eaa", "#c2e7ff", "#2f4b7c", "#8293bc", "#d6e2ff", "#665191",
                          "#a794c7", "#ebdcff", "#a05195", "#ce94c4", "#fcd8f5", "#d45087", "#ec95b6",
                          "#ffd5e5", "#f95d6a", "#ff9d9e", "#ffd6d5", "#ff7c43", "#ffac82", "#ffd9c6",
                          "#ffa600", "#ffc171", "#fbddbe", "#962B09")
      
      #Function to format the legend text with italics for genus / specific ep.(again)
      format_legend_text <- function(culture_code) {
        words <- strsplit(culture_code, " ")[[1]]
        if (length(words) >= 3) {
          return(paste0("<i>", words[1], " ", words[2], "</i>", " ", paste(words[-c(1, 2)], collapse = " ")))
        } else if (length(words) == 2) {
          return(paste0("<i>", words[1], " ", words[2], "</i>"))
        } else {
          return(culture_code)
        }
      } 
      
      #fun function for format legend
      df_average$Formatted_Legend <- sapply(df_average$Culture.code, format_legend_text)
      
      
      # Create the main plot with markers (dots) only
      plot_ly(data = df_average, x = ~as.factor(n_days), y = ~avg_increase, 
              legendgroup = ~Formatted_Legend, name = ~Formatted_Legend,
              color = ~Formatted_Legend, colors = species_colors,
              hoverinfo = "text", text = ~Culture.code,
              showlegend = FALSE) %>%
        add_trace(type = "scatter", mode = "markers", 
                  marker = list(size = 4)) %>%
        add_ribbons(ymin = ~(avg_increase - std_error), ymax = ~(avg_increase + std_error),
                    fill = "rgba(0,100,80,0.2)", line = list(color = 'transparent'),
                    showlegend = FALSE) %>%
        add_lines(showlegend = TRUE) %>%
        layout(title = "",
               xaxis = list(title = "days post-inoculation"),
               yaxis = list(title = "average growth rate"),
               legend = list(title = "Culture.code"),
               hovermode = "closest")
      
    })


###METADATA TAB
# Reactive function to filter selected columns from the DB_slim dataframe
selected_columns <- reactive({
  if (length(input$columns) > 0) {
    # Filter the selected columns
    DB_slim %>%
      select(all_of(input$columns))
  } else {
    # If no columns selected, return an empty dataframe
    data.frame()
  }
})

# Render the datatable based on the selected columns
output$data_table <- renderDT({
  req(input$columns) # Ensure that a column is selected
  selected_columns <- c("Culture.code", input$columns)
  datatable(DB_slim[, selected_columns, drop = FALSE], options = list(pageLength = 10))
})



}

###run as app
shinyApp(ui = ui, server = server)
