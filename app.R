##########################################
#This script contains all of the code necessary to run the R shiny app, SuilluScope. 
#This is version 1.0, released in August of 2023. 
#Version 1.0 highlights 25 strains of Suillus, using two different assays: culture phenotype on four different media types
#and a replicated (n=4) growth rate assay conducted across 7 temperatures.
##########################################

#load libraries
library(shiny)
library(tidyverse)
library(shinyBS)
library(plotly)
library(ggplot2)
library(ggpubr)
library(mailtoR)
library(lubridate) #to work with dates
library(xts) #to work with dates
library(data.table)
library(DT) #for working with interactive datatables
library(shinyjs)
library(plotrix)
library(plyr)
library(dplyr)


##set seed for reproducibility
set.seed(666)

#set wd
#setwd("~/Desktop/SuilluScope")

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
DB_slim$Range.and.region<- NULL



#growth data
GD<-read.delim("temperature_assay_28Jul2023.csv", header = TRUE, sep = ",", fill = TRUE, strip.white = TRUE, check.names = FALSE)

#fix names (from strain codes to the full species names + strain codes)
lookup_table <- DB %>%
  select(Culture.code, Full.Name)
GD$Culture.code <- lookup_table$Full.Name[match(GD$Culture.code, lookup_table$Culture.code)]

#remove "contam" rows and NA rows
GD_clean <- GD %>%
  dplyr::filter(!if_any(everything(), ~str_detect(., "contam")))

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


#function to calculate mean and se for each temp
data_summary <- function(data, varname, groupnames){
  #require(plyr)
  #library(plotrix)
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


###prep growth rate data for growth rate calculation
#get only 20 degree
C10_all_data<- GD_clean_with_rep[GD_clean_with_rep$Temp == 10,]
C20_all_data<- GD_clean_with_rep[GD_clean_with_rep$Temp == 20,]
C24_all_data<- GD_clean_with_rep[GD_clean_with_rep$Temp == 24,]
C27_all_data<- GD_clean_with_rep[GD_clean_with_rep$Temp == 27,]
C30_all_data<- GD_clean_with_rep[GD_clean_with_rep$Temp == 30,]
C34_all_data<- GD_clean_with_rep[GD_clean_with_rep$Temp == 34,]
C37_all_data<- GD_clean_with_rep[GD_clean_with_rep$Temp == 37,]


#function to make longer (while keeping A-D rep codes)
input_clean_long_reps<- function(data) {
  cleaned_data <- data %>%
    pivot_longer(
      cols = -c("Temp", "Culture.code", "replicate"),
      names_to = "date",
      values_to = "area",
      values_transform = list(area = as.numeric, Culture.code = as.character)
    )
  
  return(cleaned_data)
}

#run function
C10_all_data_long<- as.data.table(input_clean_long_reps(C10_all_data))
C20_all_data_long<- as.data.table(input_clean_long_reps(C20_all_data))
C24_all_data_long<- as.data.table(input_clean_long_reps(C24_all_data))
C27_all_data_long<- as.data.table(input_clean_long_reps(C27_all_data))
C30_all_data_long<- as.data.table(input_clean_long_reps(C30_all_data))
C34_all_data_long<- as.data.table(input_clean_long_reps(C34_all_data))
C37_all_data_long<- as.data.table(input_clean_long_reps(C37_all_data))

#make data.table

#change date format to read as a date in lubridate and xts
#input_clean_long_reps$date2<- lubridate::parse_date_time(input_clean_long_reps$date, "dmy")
C10_all_data_long$date2<- lubridate::parse_date_time(C10_all_data_long$date, "dmy")
C20_all_data_long$date2<- lubridate::parse_date_time(C20_all_data_long$date, "dmy")
C24_all_data_long$date2<- lubridate::parse_date_time(C24_all_data_long$date, "dmy")
C27_all_data_long$date2<- lubridate::parse_date_time(C27_all_data_long$date, "dmy")
C30_all_data_long$date2<- lubridate::parse_date_time(C30_all_data_long$date, "dmy")
C34_all_data_long$date2<- lubridate::parse_date_time(C34_all_data_long$date, "dmy")
C37_all_data_long$date2<- lubridate::parse_date_time(C37_all_data_long$date, "dmy")


#change to number of days 
#days <- yday(input_clean_long_reps$date2) - 165 #first day was Jun 14 = day 0 (the day we inoculated, which is 165 days into the year)
C10_all_data_long$n_days <- yday(C10_all_data_long$date2) - 165 #first day was Jun 14 = day 0 (the day we inoculated, which is 165 days into the year)
C20_all_data_long$n_days <- yday(C20_all_data_long$date2) - 165
C24_all_data_long$n_days <- yday(C24_all_data_long$date2) - 165
C27_all_data_long$n_days <- yday(C27_all_data_long$date2) - 165
C30_all_data_long$n_days <- yday(C30_all_data_long$date2) - 165
C34_all_data_long$n_days <- yday(C34_all_data_long$date2) - 165
C37_all_data_long$n_days <- yday(C37_all_data_long$date2) - 165


#clean up unnecessary cols
C10_all_data_long <- C10_all_data_long %>% select(-date, -date2, -Temp)
C20_all_data_long <- C20_all_data_long %>% select(-date, -date2, -Temp)
C24_all_data_long <- C24_all_data_long %>% select(-date, -date2, -Temp)
C27_all_data_long <- C27_all_data_long %>% select(-date, -date2, -Temp)
C30_all_data_long <- C30_all_data_long %>% select(-date, -date2, -Temp)
C34_all_data_long <- C34_all_data_long %>% select(-date, -date2, -Temp)
C37_all_data_long <- C37_all_data_long %>% select(-date, -date2, -Temp)

#calculate growth rate, group by each unique combination of Sp and replicate, sort by date for each group, and use lag to calculate the difference in area between that date and the one before it.
growth_rates10 <- C10_all_data_long %>%
  group_by(Culture.code, replicate) %>%
  group_modify(~mutate(., increase_in_area = area - dplyr::lag(area, default = first(area))))
growth_rates20 <- C20_all_data_long %>%
  group_by(Culture.code, replicate) %>%
  group_modify(~mutate(., increase_in_area = area - dplyr::lag(area, default = first(area))))
growth_rates24 <- C24_all_data_long %>%
  group_by(Culture.code, replicate) %>%
  group_modify(~mutate(., increase_in_area = area - dplyr::lag(area, default = first(area))))
growth_rates27 <- C27_all_data_long %>%
  group_by(Culture.code, replicate) %>%
  group_modify(~mutate(., increase_in_area = area - dplyr::lag(area, default = first(area))))
growth_rates30 <- C30_all_data_long %>%
  group_by(Culture.code, replicate) %>%
  group_modify(~mutate(., increase_in_area = area - dplyr::lag(area, default = first(area))))
growth_rates34 <- C34_all_data_long %>%
  group_by(Culture.code, replicate) %>%
  group_modify(~mutate(., increase_in_area = area - dplyr::lag(area, default = first(area))))
growth_rates37 <- C37_all_data_long %>%
  group_by(Culture.code, replicate) %>%
  group_modify(~mutate(., increase_in_area = area - dplyr::lag(area, default = first(area))))


#Remove the grouping
growth_rates10<-ungroup(growth_rates10)
growth_rates20<-ungroup(growth_rates20)
growth_rates24<-ungroup(growth_rates24)
growth_rates27<-ungroup(growth_rates27)
growth_rates30<-ungroup(growth_rates30)
growth_rates34<-ungroup(growth_rates34)
growth_rates37<-ungroup(growth_rates37)


####
#Calculate increase in area per day at room temp. 
###

DF_for_paper<- growth_rates20

#get all unique numbers of days between sampling points
unique_n_days <- unique(DF_for_paper$n_days)

#Function to calculate the number of days between n_days and the next lowest value
get_next_lowest_n_days <- function(n_days) {
  next_val_index <- which(unique_n_days < n_days)
  if (length(next_val_index) == 0) {
    next_lowest_val <- n_days[1]
  } else {
    next_lowest_val <- max(unique_n_days[next_val_index])
  }
  n_days - next_lowest_val
}

#run the function to add the new column 
DF_for_paper <- DF_for_paper %>%
  mutate(next_lowest_n_days = sapply(n_days, get_next_lowest_n_days))

#calculate average increase per day for each replicate
DF_for_paper$average_growth_per_day<- DF_for_paper$increase_in_area / DF_for_paper$next_lowest_n_days

#remove NAs 
DF_for_paper <- DF_for_paper[complete.cases(DF_for_paper[, "average_growth_per_day"]), ]

#per strain, average over all replicates for each time point 
group.means<-ddply(DF_for_paper,c("Culture.code","n_days"),summarise,ave=mean(average_growth_per_day))

#calculate overall average for each strain 
group.means.means<- aggregate(ave ~ Culture.code, group.means, mean)

#get the row index with the highest value in the 'mean' column and print it
#row_index_highest_mean <- which.max(group.means.means$ave)
#print(group.means.means[row_index_highest_mean, ]) #18 Suillus quiescens FC197 90.26757

#get the row index with the lowest value in the 'mean' column and print it
#row_index_lowest_mean <- which.min(group.means.means$ave)
#print(group.means.means[row_index_lowest_mean, ]) #21 Suillus spraguei EM44 17.26775

#get the highest value at any time point. 
#highest_overall <- group.means %>%
#  group_by(Culture.code) %>%
#  slice(which.max(ave))
#View(highest_overall)



####


# 
# #
# #get the average increase_in_area and se for each Sp and n_days combination
# df_average10 <- growth_rates10 %>%
#   group_by(Culture.code, n_days) %>%
#   do(data.frame(avg_increase = mean(.$increase_in_area),
#                 std_error = sd(.$increase_in_area) / sqrt(length(.$increase_in_area)))) %>%
#   ungroup()
# 
# df_average20 <- growth_rates20 %>%
#   group_by(Culture.code, n_days) %>%
#   do(data.frame(avg_increase = mean(.$increase_in_area),
#                 std_error = sd(.$increase_in_area) / sqrt(length(.$increase_in_area)))) %>%
#   ungroup()
# 
# df_average24 <- growth_rates24 %>%
#   group_by(Culture.code, n_days) %>%
#   do(data.frame(avg_increase = mean(.$increase_in_area),
#                 std_error = sd(.$increase_in_area) / sqrt(length(.$increase_in_area)))) %>%
#   ungroup()
# 
# df_average27 <- growth_rates27 %>%
#   group_by(Culture.code, n_days) %>%
#   do(data.frame(avg_increase = mean(.$increase_in_area),
#                 std_error = sd(.$increase_in_area) / sqrt(length(.$increase_in_area)))) %>%
#   ungroup()
# 
# df_average30 <- growth_rates30 %>%
#   group_by(Culture.code, n_days) %>%
#   do(data.frame(avg_increase = mean(.$increase_in_area),
#                 std_error = sd(.$increase_in_area) / sqrt(length(.$increase_in_area)))) %>%
#   ungroup()
# 
# df_average34 <- growth_rates34 %>%
#   group_by(Culture.code, n_days) %>%
#   do(data.frame(avg_increase = mean(.$increase_in_area),
#                 std_error = sd(.$increase_in_area) / sqrt(length(.$increase_in_area)))) %>%
#   ungroup()
# 
# df_average37 <- growth_rates37 %>%
#   group_by(Culture.code, n_days) %>%
#   do(data.frame(avg_increase = mean(.$increase_in_area),
#                 std_error = sd(.$increase_in_area) / sqrt(length(.$increase_in_area)))) %>%
#   ungroup()




##########################################
##-----------create ui object-----------##
##########################################
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
                               htmlOutput(outputId ="home_left", style = "background-color: rgba(0, 0, 0, 0.6); color: white; padding:20px; margin-top:54px; margin-bottom:35px; margin-right:355px; margin-left:-282px;")
                             ),
                             column(
                               width = 3, 
                               htmlOutput(outputId ="home_right", style = "background-color: rgba(0, 0, 0, 0.6); color: white; padding:20px; margin-top:54px; margin-bottom:35px; margin-left:-270px; margin-right:-36px;")
                             )
                           ),
                           
                           
                           useShinyjs(),
                           
                           fluidRow(
                             column(
                               width = 12,
                               div(
                                 style = "position: relative;", 
                                 img(src='lower_background.png', align = "center", height="60%", width="101%", res = 128),
                                 div(
                                   id = "citation_info",
                                   class = "well well-lg",
                                   style = "position: absolute; top: 10px; left: 10px; color: #4D4D4D; margin: 0; margin-top:5px; margin-left:450px; margin-right:450px; padding: 15px; background-color: rgba(0, 0, 0, .5); background-color: rgba(0, 0, 0, 0); border: 5px solid gray; border-radius: 20px",
                                   tags$b("HOW TO CITE", style = "font-size: 24px; color: #666666"),
                                   tags$br(),
                                   "If you find SuilluScope useful in your own research, please cite the associated publication",
                                   tags$a("HERE", href = "#", 
                                          id = "tooltip_link",
                                          style = "color: black; text-decoration: underline;"
                                   ),"as well as the paper(s) that originally published the strains you used in your work (listed in the Metadata tab)."
                                 )
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
                      style = "height:700px; width:33px; background-color:rgba(119,120,55,0.5); padding-left: 0px; padding-right: 20px"), #little buffer col.
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
                                  style = "height:700px; width:33px; background-color:rgba(119,120,55,0.5); padding-left: 0px; padding-right: 20px"), # little buffer col
                           column(2, 
                                  style = "height:700px; width:250px; background-color:rgba(119,120,55,0.5); padding-left: 0px; padding-right: 20px",
                                  br(),
                                  selectInput(inputId = "datafiles",
                                              label = "SELECT TEMPERATURE",
                                              choices = c("10°C" = "10", "20°C" = "20", "24°C" = "24", "27°C" = "27", "30°C" = "30", "34°C" = "34", "37°C" = "37")),
                                  htmlOutput(outputId ="textbox_temperature", style = "background-color: rgba(0, 0, 0, 0.5); color: white; padding:5px; margin-top:200px; margin-bottom:0px; margin-left:0px; margin-right:0px;")
                           ),
                           column(6, 
                                  plotlyOutput("temp_plot", height = "600px", width = "900px"),
                                  htmlOutput(outputId = "growth", style = "background-color:transparent; padding-right:5px; padding-left:10px; margin-top:10px; margin-left:-100px"),),
                           column(3,
                                  plotlyOutput("growth_rate_plot", height = "600px", width = "900px"),
                                  htmlOutput(outputId = "growth_rate", style = "background-color:transparent; padding-right:5px; padding-left:10px; margin-top:10px")),
                           
                         )),
                
                
                tabPanel(
                  "METADATA",
                  htmlOutput(outputId ="metadata", style = "background-color:transparent; padding:20px; margin-top:-40px"),
                  #titlePanel("Interactive Spreadsheet"),
                  sidebarLayout(
                    sidebarPanel(
                      # Checkbox group for column selection
                      style = "width: 300px;", # Adjust the width here (e.g., width: 150px;)
                      checkboxGroupInput(
                        "columns",
                        "SELECT COLUMNS TO DISPLAY",
                        choices = colnames(DB_slim),
                        selected = colnames(DB_slim)
                      )
                    ),
                    mainPanel(
                      style = "margin-left:-150px",
                      # Data table output
                      DTOutput("data_table")
                    )
                  )
                ),
                
    )
  )
)




##########################################
##---------create server object---------##
##########################################
server<- function(input, output, session){ 
  #define 'HOME' page
  output$home_left <- renderText(paste(tags$b("WELCOME TO", style = "font-size: 24px;"),
                                       tags$b("SuilluScope v1.0", style = "font-size: 24px; color: #D3AD0D"),
                                       tags$br(),
                                       "SuilluScope is an interactive, open access database for <i>Suillus</i> fungi, a model genus for ectomycorrhizal ecology and evolution.",
                                       "There are extensive genomic resources available for", tags$i("Suillus"), "including more annotated genome assemblies 
                                       than for any other ectomycorrhizal group.",
                                       "This platform proves data on the phenotypic traits and responses of these genome-sequenced isolates, aiming to help researchers identify optimal 
                                       culture conditions, predict and compare trait responses across species, and empower further 
                                       research linking genotypes to phenotypes.",
                                       tags$br()
  ))
  
  
  output$home_right<- renderText(paste(tags$b("ISOLATE COLLECTION", style = "font-size: 24px; color: #D3AD0D"),
                                       tags$br(),
                                       "All isolates used in the database are from the", tags$i("Suillus"), "genome culture collection, curated at Duke University.", 
                                       "For more information on the <i>Suillus</i> system, please visit",
                                       tags$a(href = "http://www2.hawaii.edu/~nn33/suillus/", "The International ", tags$i("Suillus"), " Consortium.", target = "_blank"),
                                       "Genomic resources for these cultures are available on the MycoCosm ",tags$i("Suillus"), 
                                       tags$a(href="https://mycocosm.jgi.doe.gov/Suillus/Suillus.info.html", "web portal", target="_blank"), 
                                       "opperated by DOE Joint Genome Institute.",
                                       "Isolates amenable to cryopreservation will soon be avalible as part of the", 
                                       tags$a(href="https://nrrl.ncaur.usda.gov/", "Agricultural Research Culture Collection (NRRL).", target="_blank"), 
                                       "Isolates which cannot be cryopreserved are available from the authors", mailtoR(email = "LotusLofgren@gmail.com", text = "here.")
  ))
  
  #render shinyjs tool for citation pop up 
  shinyjs::runjs('
    $("#tooltip_link").click(function() {
      if ($("#custom_tooltip").length > 0) {
        // If tooltip is already visible, hide it
        $("#custom_tooltip").remove();
      } else {
        // If tooltip is hidden, show it to the right of the link
        var tooltipContent = "<div style=\'width: 350px; word-wrap: break-word;\'> <b><i>Suillus</i>: an emerging model for the study of ectomycorrhizal ecology and evolution.</b> <br>Lofgren L, Nguyen NH, Kennedy P, Pérez-Pazo4 E, Fletcher J, Liao H-L, Wang H, Zhang K, Ruytinx J, Smith A, Ke Y-H, Cotter H.V., Engwall E, Hameed K.M., Vilgalys R, Branco S. <br><i>In review</i></div>";
        var tooltipHtml = \'<div id="custom_tooltip" style="position: absolute; top: 50%; left: 100%; transform: translate(-5px, -50%); background-color: white; border: 1px solid black; padding: 10px; max-width: 400px;">\' + tooltipContent + \'</div>\';
        $(this).parent().append(tooltipHtml);
      }
    });
    
    // Hide tooltip when clicking outside of it
    $(document).on("click", function(event) {
      if (!$(event.target).closest("#tooltip_link, #custom_tooltip").length) {
        $("#custom_tooltip").remove();
      }
    });
  ')
  
  output$home_bottom2 <- renderText(paste("This is version v1.0 of the database, released in beta on 28.July.2023",
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
  
  
  output$textbox_morpholgy <- renderText(paste(tags$b("METHODS:"), 
                                               tags$br(), 
                                               "Cultures were grown on four media types including Modified Melin-Norkrans (MMN), 
                                               Modified Fries Media (Fries), Modified Hagem’s Agar (Hagem’s), and Pachlewski’s Media (Pachlewski’s, or Px). 
                                               All media types were prepared at their full respective carbon concentrations and adjusted to pH 6 prior to autoclaving. 
                                               Cultures were grown for 28 days, at room temperature, in the dark, prior to being photographed."
  ))
  
  output$textbox_temperature <- renderText(paste(tags$b("METHODS:"), 
                                                 tags$br(), 
                                                 "Cultures were started by placing 3mm agar plugs on Modified Melin-Norkrans (MMN) media, 
                                                 adjusted to pH 6 prior to autoclaving, and grown in temperature adjustable incubators, 
                                                 in the dark, for a total of 33 days (at n=4 replicates per species per temperature treatment). 
                                                 Starting on day 8, colony area was recorded twice per week over the course of the assay by 
                                                 marking the colony margin on the back of each petri dish with a fine-tip sharpie. 
                                                 After 33 days of growth, the back of the petri dishes were imaged using a flatbed scanner, 
                                                 and colony area was calculated for each time point using the program imageJ. "
  ))
  
  
  
  
  
  #### 
  #TEMPERATURE PLOT
  ####
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
  
  
  #set rendering variables for formatting
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
    
    #create a new column in the with the formatted legend text
    dataset$Formatted_Legend <- sapply(dataset$Culture.code, format_legend_text)
    
    #cefine the color palette
    species_colors <- c("#003f5c", "#668eaa", "#c2e7ff", "#2f4b7c", "#8293bc", "#d6e2ff", "#665191",
                        "#a794c7", "#ebdcff", "#a05195", "#ce94c4", "#fcd8f5", "#d45087", "#ec95b6",
                        "#ffd5e5", "#f95d6a", "#ff9d9e", "#ffd6d5", "#ff7c43", "#ffac82", "#ffd9c6",
                        "#ffa600", "#ffc171", "#fbddbe", "#962B09")
    
    #plot that shit
    plot_ly(data = dataset, x = ~n_days, y = ~area, type = "scatter", mode="lines+markers", line = list(width = 2),
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
  
  
  
  ###
  ##METADATA TAB
  ###
  
  #define the choices and format with with spaces instead of periods for pretty rendering
  choices <- setdiff(colnames(DB_slim), "Culture.code")
  choices <- gsub("\\.", " ", choices) # Replace periods with spaces
  
  #update checkboxs for group choices when the app starts
  observe({
    updateCheckboxGroupInput(session, inputId = "columns", choices = choices)
  })
  
  #render the datatable based on the selected columns
  output$data_table <- renderDT({
    req(input$columns) # Ensure that a column is selected
    
    #modify column names of the DB_slim dataframe to match the pretty format
    colnames(DB_slim) <- gsub("\\.", " ", colnames(DB_slim)) # Replace periods with spaces
    
    selected_columns <- c("Culture.code", input$columns)
    # Filter the selected columns, excluding "Culture.code"
    selected <- setdiff(selected_columns, "Culture.code")
    
    datatable(DB_slim[, selected, drop = FALSE], options = list(pageLength = 10))
  })
  
  
  
}


##########################################
##-----------run as app-----------
##########################################
shinyApp(ui = ui, server = server)
