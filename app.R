##########################################
#This script contains all of the code necessary to run the R shiny app, SuilluScope. 
#This is version 1.0b (beta), released in August of 2023. 
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
GD<-read.delim("temperature_assay_21Aug2023.csv", header = TRUE, sep = ",", fill = TRUE, strip.white = TRUE, check.names = FALSE)


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
ui = navbarPage(
  fluidRow(
    column(3, img(src='Suilluscope_logo.png', align = "left", height="110px", width="110px", res = 128), style = "padding-top:0px; padding-bottom:20px; padding-left:10px; padding-right:0px; margin-right:-75px", class = "logo-column"),
    column(9,
           div(id = "home_left"),
           div(id = "home_right"),
           tags$head(tags$style(HTML("
             .navbar, .navbar-collapse {
               background-color: white;
               border: none;
               margin-bottom: 12px; /* note, this is how you adjust the space between the nav row and the rest of the page */
             }

             .navbar .nav {
               padding-left: 200px; 
               background-color: white;
               margin-top: 40px;
               margin-bottom: 25px;
             }
             .tabbable > .nav {
               padding-top: 60px;
               background-color: white;
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
             .logo-column {
               margin-bottom: 20px; 
             }
             #home_left a, #home_right a{
               color: #DDC87B !important; /* Set the link color to yellow */
               text-decoration: none;
             }
            #home_bottom2 {
               color: black !important; 
               text-decoration: none;
             }
              #home_bottom2 a {
               color: #DDC87B !important; /* Set the link color to yellow */
               text-decoration: underline;
             }
             @media (max-width: 767px) { /* Apply styles only mobile < than 768px width */
               .navbar {
                 margin-bottom: 100px; 
               }
               .logo-column {
                 margin-bottom: 30px; 
               }
               .background-image-container {
                 height: auto;
                 margin-bottom: 50px;
               }
               .tab-content .tab-pane {
                 padding: 0; /* Remove padding for tab content */
               }
               .html-widget {
                 padding: 10px; 
               }
            @media (max-width: 768px) {
              .image-column {
              width: 100%; /* Take up full width on mobile */
              }
              }
            @media (max-width: 768px) {
              .green-column {
              width: 100%; /* Take up full width on mobile*/
              }
              }
             #data_table_container table {
             width: 100%; /* so the table takes up 100% width of its container */
             max-width: none; 
             }
             @media (max-width: 768px) {
             #data_table_container table {
             font-size: 12px; /* Adjust font size for mobile */
             margin-left: 150px; /* Adjust margin for mobile, since we start with it nudged to the left for desktop */
             }
             }
}
           ")), 
           )
    )
  )
  
  
  ,
  
  tabPanel("HOME",
           tags$div(
             class = "background-image-container-top",
             fluidRow(
               column(
                 width = 6,
                 htmlOutput(outputId ="home_left", class = "html-widget", style = "background-color: rgba(0, 0, 0, 0.6); color: white; padding:20px; margin-top:35px; margin-bottom:35px; margin-right:60px; margin-left:60px;")
               ),
               column(
                 width = 6,
                 htmlOutput(outputId ="home_right", class = "html-widget", style = "background-color: rgba(0, 0, 0, 0.6); color: white; padding:20px; margin-top:35px; margin-bottom:35px; margin-left:60px; margin-right:60px;")
               )
             )
           ),
           useShinyjs(),
           tags$div(
             class = "background-image-container-bottom",
             fluidRow(
               column(
                 width = 12,
                 div(
                   id = "citation_info",
                   class = "well well-lg",
                   style = "position: relative; color: #4D4D4D; margin: 0 auto; max-width: 600px; padding: 15px; background-color: rgba(0, 0, 0, .1); border: 5px solid gray; border-radius: 20px; margin-top:35px; margin-bottom:135px", 
                   tags$b("HOW TO CITE", style = "font-size: 24px; color: #666666"),
                   tags$br(),
                   "If you find SuilluScope useful in your own research, please cite the associated publication",
                   tags$a("HERE", href = "#",
                          id = "tooltip_link",
                          style = "color: black; text-decoration: underline;"
                   ),"as well as the paper(s) that originally published the strains you used in your work (listed in the Metadata tab)."
                 )
               ),
               column(
                 width = 12,
                 div(
                   style = "position: relative;", 
                   htmlOutput(outputId ="home_bottom2", style = "color: black; margin-left: 30px; margin-right: 30px; margin-bottom:15px")
                 )
               )
             )
           ),
           tags$head(tags$style(HTML("
    .background-image-container-top {
      background-image: url('background_top.png');
      background-size: cover;
    }
    .background-image-container-bottom {
      background-image: url('background_bottom.png');
      background-size: cover;
      margin-top: -1px; /* to remove the white line between the two background images */
    }
  ")))
  ),
  
  
  
  tabPanel(
    "MORPHOLOGY",
    fluidRow(  tags$style(
      HTML(
        "@media (max-width: 768px) {
         .image-column {
           width: 100%;
           padding-left: 0;
           padding-right: 0;
         }}

         #plate_img {
           padding-left: 0;
           padding-right: 15px;
           min-width: 350px;
           margin-left: 25px
         }
         
         .desktop-only {
            display: none;
          }
          
          .mobile-only {
            display: block;
          }
      }"
      )
    ),
    column( 
      3, 
      class = "image-column col-md-3 col-12",
      style = "height: 700px; background-color: rgba(119,120,55,0.5); padding-left: 25px; padding-right: 20px; margin-bottom: -600px; margin-left: 15px; margin-right: -15px",
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
      img.setAttribute("width", "100%");  // Set width to 100%
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
    ),useShinyjs(),
    column(
      9,
      class = "image-column",
      div(id = "plate_img"),
    ),
    div(
      id = "textbox_morpholgy",
      class = "well well-lg",
      style = "position: absolute; bottom: 0; color: #4D4D4D; margin: 0 auto; width: auto; padding: 15px; background-color: #A8A78A; border: 0px solid gray; border-radius: 20px; margin-bottom: 0px; margin-top: -250; margin-left: 34px; max-width: 320px",
      tags$a("Methods", href = "#",
             id = "textbox_morphology_link",
             style = "color: black; text-decoration: underline; font-weight: bold;"
      ),
      div(
        id = "textbox_morpholgy_content",
        style = "display: none; max-height:300px; width: auto; overflow-y: auto"
      )
    )
    ),
  )
  ,
  
  
  
  tabPanel(
    "TEMPERATURE",
    fluidRow(  tags$style(
      HTML(
        "@media (max-width: 768px) {
         .green-column {
           width: 100%;
           padding-left: 0;
           padding-right: 0;
         }
         
         #temp_img {
           padding-left: 0;
           padding-right: 0;
           min-width: 350px;
         }
         
         .desktop-only {
            display: none;
         }
         
         .mobile-only {
            display: block;
         }
         
         #temp_plot {
            max-height: 300px; 
         }
      }"
      )
    ),
    column( 
      3, 
      class = "image-column col-md-3 col-12",
      style = "height: 700px; background-color: rgba(119,120,55,0.5); padding-left: 25px; padding-right: 20px; margin-bottom: -600px; margin-left: 15px; margin-right: -15px",
      br(),
      selectInput(
        inputId = "datafiles",
        label = "SELECT TEMPERATURE",
        choices = c("10°C" = "10", "20°C" = "20", "24°C" = "24", "27°C" = "27", "30°C" = "30", "34°C" = "34", "37°C" = "37")
      )
      ,
    ),useShinyjs(),
    column(
      9,
      class = "green-column",
      div(id = "temp_img"),
      plotlyOutput("temp_plot", height = "500px")
    ),
    div(
      id = "textbox_temperature",
      class = "well well-lg",
      style = "position: absolute; bottom: 0; color: #4D4D4D; margin: 0 auto; width: auto; padding: 15px; background-color: #A8A78A; border: 0px solid gray; border-radius: 20px; margin-bottom: 0px; margin-top: -250; margin-left: 34px; max-width: 320px",
      tags$a("Methods", href = "#",
             id = "textbox_temperature_link",
             style = "color: black; text-decoration: underline; font-weight: bold;"
      ),
      div(
        id = "textbox_temperature_content",
        style = "display: none; max-height:300px; width: auto; overflow-y: auto"
      )
    )
    ),
  )
  ,
  
  
  
  tabPanel(
    "METADATA",
    htmlOutput(outputId ="metadata", style = "background-color:transparent; padding:20px; margin-top:-40px"),
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
      #render data table
      mainPanel(
        div(
          style = "margin-left:-150px",
          id = "main_panel_content",
          div(DTOutput("data_table"), id = "data_table_container"))
      )
    )
  )
  
  , id = "navBar", collapsible = TRUE # we can change to see when we shrink the windows: TRUE will pop up the hamburger menu icon 
)

# Set the title of the browser tab using JavaScript
js_code <- '
  document.title = "SuilluScope";
'

ui <- tagList(
  tags$head(tags$script(HTML(js_code))),
  ui
)



##########################################
##---------create server object---------##
##########################################
server<- function(input, output, session){ 
  #define 'HOME' page
  output$home_left <- renderText(paste(tags$b("WELCOME TO", style = "font-size: 24px;"),
                                       tags$b("SuilluScope", style = "font-size: 24px; color: #D3AD0D"),
                                       tags$br(),
                                       "SuilluScope is an interactive, open access database for <i>Suillus</i> fungi, a model genus for ectomycorrhizal ecology and evolution.",
                                       "There are extensive genomic resources available for", tags$i("Suillus"), "including more annotated genome assemblies 
                                       than for any other ectomycorrhizal group.",
                                       "This platform provides data on the phenotypic traits and responses of these genome-sequenced isolates, aiming to help researchers identify optimal 
                                       culture conditions, predict and compare trait responses across species, and empower further 
                                       research linking genotypes to phenotypes.",
                                       tags$br()
  ))
  
  
  output$home_right<- renderText(paste(tags$b("ISOLATE COLLECTION", style = "font-size: 24px; color: #D3AD0D"),
                                       tags$br(),
                                       "All isolates used in the database are part of the", tags$i("Suillus"), "Genome Strain Culture Collection.", 
                                       "Isolates amenable to cryopreservation have been integrated into the", 
                                       tags$a(href="https://nrrl.ncaur.usda.gov/", "Agricultural Research Service Culture Collection (NRRL).", target="_blank"), 
                                       "Please see the Metadata tab for NRRL accession numbers.",
                                       "Isolates which cannot be cryopreserved are available by contacting the author", mailtoR(email = "LotusLofgren@gmail.com", text = "here."),
                                       "Genomic resources for these cultures are available on the MycoCosm ",tags$i("Suillus"), 
                                       tags$a(href="https://mycocosm.jgi.doe.gov/Suillus/Suillus.info.html", "web portal", target="_blank"), 
                                       "operated by DOE Joint Genome Institute.",
                                       "For more information on the <i>Suillus</i> system, please visit",
                                       tags$a(href = "http://www2.hawaii.edu/~nn33/suillus/", "The International ", tags$i("Suillus"), " Consortium.", target = "_blank")
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
  
  
  #render shinyjs tool for methods popup (morphology)
  shinyjs::runjs('
  $(document).ready(function() {
    $("#textbox_morphology_link").click(function() {
      var tooltipContent = "<div style=\'width: 340px; word-wrap: break-word;\'>"
                          + "Cultures were grown on four media types<br>" 
                          + "including Modified Melin-Norkrans (MMN),<br>"
                          + " Modified Fries Media (Fries), Modified <br>"
                          + "Hagem’s Agar (Hagem’s), and Pachlewski’s<br>"
                          + "Media (Pachlewski’s, or Px). All media types<br>"
                          + "were prepared at their full respective carbon concentrations and adjusted to pH 6 prior to<br> autoclaving."
                          + " Cultures were grown for 28 days,<br> at room temperature, in the dark, prior to<br>"
                          + "being photographed.</div>";

      var tooltipHtml = \'<div id="custom_tooltip_content">\' + tooltipContent + \'</div>\';
      $("#textbox_morpholgy_content").html(tooltipHtml).show();
    });

    // Hide tooltip when clicking outside of it
    $(document).on("click", function(event) {
      if (!$(event.target).closest("#textbox_morphology_link, #custom_tooltip_content").length) {
        $("#textbox_morpholgy_content").hide();
      }
    });
  });
')
  
  
  
  #render shinyjs tool for methods popup (temperature)
  shinyjs::runjs('
  $(document).ready(function() {
    $("#textbox_temperature_link").click(function() {
      var tooltipContent = "<div style=\'width: 340px; word-wrap: break-word;\'>"
                          + "Cultures were started by placing 3mm agar<br>"
                          + "plugs on Modified Melin-Norkrans (MMN)<br>"
                          + "media adjusted to pH 6 prior to autoclaving,<br>"
                          + "and grown in temperature adjustable<br>" 
                          + "incubators, in the dark, for a total of 33 days<br>"
                          + "(at n=4 replicates per species per temperature treatment)." 
                          + " Starting on day 8, colony area was recorded twice per week over the course<br>"
                          + "of the assay by marking the colony margin on<br>"
                          + "the back of each petri dish with a fine-tip<br>"
                          + "sharpie. After 33 days of growth, the back of<br>"
                          + "the petri dishes were imaged using a flatbed<br>" 
                          + "scanner, and colony area was calculated for<br>" 
                          + "each time point using the program imageJ.</div>";

      var tooltipHtml = \'<div id="custom_tooltip_content">\' + tooltipContent + \'</div>\';
      $("#textbox_temperature_content").html(tooltipHtml).show();
    });

    // Hide tooltip when clicking outside of it
    $(document).on("click", function(event) {
      if (!$(event.target).closest("#textbox_temperature_link, #custom_tooltip_content").length) {
        $("#textbox_temperature_content").hide();
      }
    });
  });
')
  
  
 
##Add version in formation and SuilluScope logo to homepage
  output$home_bottom2 <- renderUI({
    text <- "This is version v1.0Beta of the database, released on 28.July.2023"
    text <- paste(text, tags$br(),
                  "SuilluScope was built using the open source programming language",
                  tags$a(href = "https://www.r-project.org/about.html", "R", target = "_blank"),
                  "with reactive programming via",
                  tags$a(href = "https://shiny.rstudio.com/", "R shiny", target = "_blank"),
                  tags$br(),
                  "All of the code necessary to run the program is publicly available at",
                  tags$a(href = "https://github.com/MycoPunk/SuilluScope", "the SuilluScope GitHub.", target = "_blank"),
                  "Please report issues to the git issues page.",
                  tags$br(),
                  "If you have feature requests or a published dataset that you would like us to consider adding to this site, please contact the author", mailtoR(email = "LotusLofgren@gmail.com", text = "here.")
    )
    
    #Add SuilluScope Logo
    image_path <- "ISC_logo_2023.png"
    image_tag <- tags$a(href = "http://www2.hawaii.edu/~nn33/suillus/",
                        tags$img(src = image_path, style = "max-width: 100px; height: auto;"))

    
    #Combine text and logo within div elements
    image_div <- tags$div(style = "text-align: center; margin-top: -110px; margin-bottom: 75px;", image_tag)
    text_div <- tags$div(style = "text-align: left; margin-top: -50px;", HTML(text))
    
    # Return the UI content
    return(tagList(image_div, text_div))
  })
  
  
  

  
  
  
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
  
  ###
  ##CLOSE NAVBAR HAMBURGER AFTER CLICK
  ###
  observeEvent(input$navBar, {
    runjs('
      var elem = document.getElementsByClassName("navbar-collapse")[0]
      elem.setAttribute("aria-expanded", "false");
      elem.setAttribute("class", "navbar-collapse collapse");
    ')
  })
  
  
  
}


##########################################
##-----------run as app-----------
##########################################
shinyApp(ui = ui, server = server)