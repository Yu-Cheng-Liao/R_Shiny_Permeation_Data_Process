library(shiny)
library(readxl)
library(stringr)
library(tidyverse)
library(xlsx)
library(mmand)
library(shinyalert)


ui <- 
  fluidPage(
    titlePanel("Please Upload Permeation Data"),
    sidebarLayout(
      sidebarPanel(
        tags$a(href="https://drive.google.com/file/d/1Yb-na48zLb2EYytdHHG7EWZKDWoFzH93/view?usp=sharing",
               tags$b("Input Data Template")),
        tags$br(),
        tags$a(href="https://drive.google.com/file/d/1mJOW2jvNLw-im6stY1vuXewYHiFkNPRa/view?usp=sharing",
               tags$b("Input Data Sample")),
        h4("Please click the link above, download excel file, fill out and upload in below box."),
        h4("After uploading the file, Press `Download Xlsx`"),
        fileInput("file", "Choose Xlsx File to Upload",
                  accept = c(".xlsx")),
        useShinyalert(),
        #actionButton("Run", "Run!"),
        downloadButton("Download", "Download Xlsx"),
        width = 8),
      mainPanel(
        tableOutput("contents"),
      )
    ),
    fluidRow(
      column(8,
             h4("Standard Curve Formula : "),
             textOutput("Standard_Conc")
      )
    )
  )


server <- function(input, output) {
  state <- reactiveValues()
  final_out <- reactiveValues()
  
  observe({
    req(input$file)
    inFile <- input$file
    
    state$inFile <- inFile
    
    if(is.null(inFile))
      return(NULL)
    file.rename(inFile$datapath,
                paste(inFile$datapath, ".xlsx", sep=""))
    
    state$df1 <- 
      read_excel(paste(inFile$datapath, ".xlsx", sep=""), sheet = 1)
    state$df2 <- 
      read_excel(paste(inFile$datapath, ".xlsx", sep=""), sheet = 2)
    state$df3 <- 
      read_excel(paste(inFile$datapath, ".xlsx", sep=""), sheet = 3)
    state$df4 <- 
      read_excel(paste(inFile$datapath, ".xlsx", sep=""), sheet = 4)
    
    Formulation <- state$df1
    Area <- state$df2
    Condition <- state$df3
    Curve <- state$df4
    
    Area_long <- 
      Area %>% 
      pivot_longer(-No, names_to = "Time", values_to = "Peak_Area") %>% 
      drop_na()
    
    Area_long$Time <- as.numeric(Area_long$Time)
    Area_long$Peak_Area <- as.numeric(Area_long$Peak_Area)
    
    #Curve
    #R-squared < 0.998 error
    std_model <- lm(Peak_Area_Average ~ Conc, data = Curve)
    std_model_coef <- coef(std_model) %>% as.data.frame()
    summary(std_model)
    
    #R-squared
    Rsquared <- 
      summary(std_model)$r.squared
    Slope <-
      coef(std_model)["Conc"] %>%
      round(digits = 2) %>%
      as.character()
    Intercept <-
      coef(std_model)["(Intercept)"] %>%
      round(digits = 2) %>%
      as.character()
    
    png(file = "C:/Permeation/png/curve.png", width = 450, height = 300)
    plot(Peak_Area_Average ~ Conc,
         data = Curve,
         main = paste0("Average peak Area = ",
                       Slope,
                       " * Concentration + ",
                       Intercept,
                       "\n",
                       "R.squared = ",
                       Rsquared))
    abline(std_model, col="red")
    dev.off()
    
    #
    Time_object <- 
      Area_long$Time %>% 
      unique() %>% 
      as.numeric() %>% 
      sort()
    #
    Area_long_Conc <- 
      Area_long %>% 
      drop_na() %>% 
      mutate(
        Conc_ug_ml = 
          threshold((Peak_Area - std_model_coef[1,]) / std_model_coef[2,], 0,
                    binarise = F),
        Conc_ug = Conc_ug_ml * 14)
    
    Qt_ug <- list()
    for (i in seq_along(Area_long_Conc$Conc_ug)) {
      # well_number <- Area_long_Conc[Area_long_Conc$Conc_ug == Area_long_Conc$Conc_ug[i], 1]
      # Time_number <- Area_long_Conc[Area_long_Conc$Conc_ug == Area_long_Conc$Conc_ug[i], 2]
      well_number <- Area_long_Conc[i, "No"]
      Time_number <- Area_long_Conc[i, "Time"]
      previous_Conc <- ifelse(which(Time_object == as.numeric(Time_number)) == 1, 0, 
                              Area_long_Conc[Area_long_Conc$No == as.numeric(well_number) &  
                                               Area_long_Conc$Time == Time_object[which(Time_object == as.numeric(Time_number)) - 1], 
                                             "Conc_ug_ml"]) %>% 
        as.numeric()
      previous_total <- previous_Conc * (14- Condition$Sampling_Volume)
      Qt_ug[i] <- Area_long_Conc$Conc_ug[i] - previous_total
    }
    
    Qt_ug_tibble <- 
      as_tibble_col(Qt_ug, column_name = "Qt_ug")
    
    Area_long_Conc <- 
      cbind(Area_long_Conc, Qt_ug_tibble) 
    
    Area_long_Conc$Qt_ug <-  Area_long_Conc$Qt_ug %>% as.numeric()
    
    #Q(0-t)/cm2(μg/cm2)
    Area_long_Conc <- 
      Area_long_Conc %>% 
      group_by(No) %>% 
      mutate(
        Qt_ug_cm2 = cumsum(Qt_ug) / Condition$Surface_of_Skin
      )
    
    for(j in seq_along(Area_long_Conc$Qt_ug)) {
      this_time <- Area_long_Conc[j, "Time"] %>% as.numeric()
      Area_long_Conc$Time_diff[j] <- 
        this_time - ifelse(which(Time_object == this_time) == 1, 0, Time_object[(which(Time_object == this_time)-1)])
    }
    
    Area_long_Conc$Flux <- 
      Area_long_Conc$Qt_ug / Condition$Surface_of_Skin / Area_long_Conc$Time_diff
    
    #test Average data & plot ####
    T_F <- Formulation[, c(1:3)]
    
    plot_avg_qt <- list()
    plot_avg_qt_sk <- list()
    plot_avg_fx <- list() 
    plot_avg_fx_sk <- list()
    Average_Amount <- list()
    Average_Amount_by_skin <- list()
    Average_Flux <- list()
    Average_Flux_by_skin <- list()
    for (k in names(T_F[-c(1, 3)])) {
      
      ##### 
      #plot_avg_qt_v2
      dtf <- T_F[-c(1, 3)]
      dtn <- names(dtf[k])
      #qt by variables 
      dataset <- 
        T_F %>% 
        left_join(Area_long_Conc) %>%
        select(-c("No")) %>% 
        group_by_(dtn, "Time") %>% 
        summarise(Average = mean(Qt_ug_cm2, na.rm = T),
                  sd = sd(Qt_ug_cm2, na.rm = T)) %>%
        drop_na()
      
      p <- 
        ggplot(dataset,
               aes_string(
                 x = "Time",
                 y = "Average",
                 color = dtn,
                 group = dtn
               )
        ) +
        geom_line(stat = "identity", size = 1) +
        geom_point() +
        geom_errorbar(aes(ymin=Average-sd, ymax=Average+sd), width=.2) +
        labs(title = paste0(
          Condition$Project, " " , Condition$Compound,
          "\n Cumulative Amount by ", names(dtf[k])),
          color = paste0(names(dtf[k]))) +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_x_continuous(name = "Sampling Time")
      name <- paste(names(dtf[k]), sep = "_")
      list_p <- list(p)
      plot_avg_qt[[name]] <- list_p
      
      #qt by variables & skin
      dataset <- 
        T_F %>%
        left_join(Area_long_Conc) %>%
        select(-c("No")) %>%
        group_by_(dtn, "Skin", "Time") %>%
        summarise(Average = mean(Qt_ug_cm2, na.rm = T)) %>%
        unite(interact, as.name(dtn), Skin, sep = "-") %>% 
        drop_na()
      
      p <- 
        ggplot(
          dataset,
          aes_string(
            x = "Time",
            y = "Average",
            # group = dtn,
            color = "interact",
            group = "interact"
            # color = interaction(dtn, "Skin")
          )
        ) + 
        geom_line(stat = "identity", size = 1) +
        geom_point() +
        labs(title = paste0(
          Condition$Project, " " ,Condition$Compound,
          "\n Cumulative Amount by ", names(dtf[k]), " by Skin"),
          color = paste0(names(dtf[k]))) +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_x_continuous(name = "Sampling Time")
      name <- paste(names(dtf[k]), sep = "_")
      list_p <- list(p)
      plot_avg_qt_sk[[name]] <- list_p
      
      #plot_avg_fx_v2
      #fx by variable
      dataset <- 
        T_F %>%
        left_join(Area_long_Conc) %>%
        select(-c("No")) %>%
        group_by_(dtn, "Time") %>%
        summarise(Average = mean(Flux, na.rm = T),
                  sd = sd(Flux, na.rm = T)) %>%
        drop_na()
      p <- 
        ggplot(
          dataset,
          aes_string(
            x = "Time",
            y = "Average",
            group = dtn,
            color = dtn
          )
        ) + 
        geom_line(stat = "identity", size = 1) +
        geom_point() +
        geom_errorbar(aes(ymin=Average-sd, ymax=Average+sd), width=.2) +
        labs(title = paste0(
          Condition$Project," ",Condition$Compound,
          "\n Flux by ", names(dtf[k])),
          color = paste0(names(dtf[k]))) +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_x_continuous(name = "Sampling Time")
      name <- paste(names(dtf[k]), sep = "_")
      list_p <- list(p)
      plot_avg_fx[[name]] <- list_p
      
      
      #fx by variable & skin
      dataset <- 
        T_F %>%
        left_join(Area_long_Conc) %>%
        select(-c("No")) %>%
        group_by_(dtn, "Skin", "Time") %>%
        summarise(Average = mean(Flux, na.rm = T)) %>%
        unite(interact, as.name(dtn), Skin, sep = "-") %>%
        drop_na()
      p <- 
        ggplot(
          dataset,
          aes_string(
            x = "Time",
            y = "Average",
            group = "interact",
            color = "interact"
          )
        ) + 
        geom_line(stat = "identity", size = 1) +
        geom_point() +
        labs(title = paste0(
          Condition$Project," ",Condition$Compound,
          "\n Flux by ", names(dtf[k]), " by Skin"),
          color = paste0(names(dtf[k]))) +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_x_continuous(name = "Sampling Time")
      name <- paste(names(dtf[k]), sep = "_")
      list_p <- list(p)
      plot_avg_fx_sk[[name]] <- list_p
      
      #
      SKY_plot <- list(plot_avg_fx_sk = plot_avg_fx_sk,
                       plot_avg_qt_sk = plot_avg_qt_sk)
      SKN_plot <- list(plot_avg_fx = plot_avg_fx,
                       plot_avg_qt = plot_avg_qt)
      Avg_plot <- list(SKY_plot = SKY_plot,
                       SKN_plot = SKN_plot)
      
      #Average Qt_ug_cm2 by variable  
      Average_Amount[[dtn]] <-
        T_F %>%
        left_join(Area_long_Conc) %>%
        group_by_(dtn, "Time") %>%
        drop_na() %>%
        summarise(Average = mean(Qt_ug_cm2, na.rm = T)) %>%
        pivot_wider(id_cols =dtn,
                    names_from = Time,
                    values_from = Average) %>% 
        mutate_if(is.numeric, round, digits = 2)
      
      #Average Qt_ug_cm2 by variable & skin  
      Average_Amount_by_skin[[dtn]] <-
        T_F %>%
        left_join(Area_long_Conc) %>%
        group_by_(dtn, "Skin", "Time") %>%
        drop_na() %>%
        summarise(Average = mean(Qt_ug_cm2, na.rm = T)) %>%
        pivot_wider(id_cols = c(dtn, Skin),
                    names_from = Time,
                    values_from = Average) %>% 
        mutate_if(is.numeric, round, digits = 2)
      #Average Flux by variable
      Average_Flux[[dtn]] <-
        T_F %>%
        left_join(Area_long_Conc) %>% 
        group_by_(dtn, "Time") %>%
        drop_na() %>% 
        summarise(Average = mean(Flux, na.rm = T)) %>%
        pivot_wider(id_cols = dtn,
                    names_from = Time,
                    values_from = Average) %>%
        mutate_if(is.numeric, round, digits = 2) 
      #Average Flux by variable & skin
      Average_Flux_by_skin[[dtn]] <-
        T_F %>%
        left_join(Area_long_Conc) %>% 
        group_by_(dtn, "Skin", "Time") %>%
        drop_na() %>% 
        summarise(Average = mean(Flux, na.rm = T)) %>%
        pivot_wider(id_cols = c(dtn, Skin),
                    names_from = Time,
                    values_from = Average) %>%
        mutate_if(is.numeric, round, digits = 2) 
      
    }
    SKY <- list(Average_Flux_by_skin = Average_Flux_by_skin,
                Average_Amount_by_skin = Average_Amount_by_skin)
    
    SKN <- list(Average_Flux = Average_Flux,
                Average_Amount = Average_Amount)
    Avg_data <- list(SKY = SKY,
                     SKN = SKN)
    
    End_data <- list()
    for(m in names(Area_long_Conc[c(4:7, 9)])) {
      End_data[[m]] <-   
        Area_long_Conc %>% 
        pivot_wider(
          id_cols = "No",
          names_from = "Time",
          values_from = m) %>%
        right_join(Formulation[, c(1:3)]) %>% 
        select("No", "Formulation_No", "Skin", everything()) %>% 
        arrange(No) %>% 
        mutate_if(is.numeric, round, digits = 2)
    }
    state$Area_long_Conc <- Area_long_Conc
    state$End_data <- End_data
    state$T_F <- T_F
    state$Avg_data <- Avg_data
    state$Avg_plot <- Avg_plot
    state$SKY_plot <- SKY_plot
    state$SKN_plot <- SKN_plot
    
    #final_out$out <- "finish"
    
  })
  
  output$contents <- renderTable(
    {
      print(state$df4)
    }
  )
  
  output$Standard_Conc <- renderText({
    
    req(input$file)
    
    std_model <-
      lm(Peak_Area_Average ~ Conc, data = state$df4)
    rsquared <- summary(std_model)$r.squared
    Slope <-
      coef(std_model)["Conc"] %>%
      round(digits = 2) %>%
      as.character()
    Intercept <-
      coef(std_model)["(Intercept)"] %>%
      round(digits = 2) %>%
      as.character()
    Standard_Conc <- paste0("Average peak Area = ",
                            Slope,
                            " * Concentration + ",
                            Intercept,
                            " ; ",
                            "R.squared = ",
                            rsquared)
    Standard_Conc
  })
  
  output$Download <- downloadHandler(
    
    filename = function(){
      paste0(Sys.Date(), "_", state$df3$Project, "_", "Permeation.xlsx")
    },
    content = function(file){
      #wb <- createWorkbook()
      wb <- loadWorkbook(paste(state$inFile$datapath, ".xlsx", sep=""))
      ws <- createSheet(wb, sheetName = "Permeation") #calculation process
      wa <- createSheet(wb, sheetName = "Avg_plot") # Avg & plot
      #Standard Curve df & plot save to xlsx----
      title_style <- CellStyle(wb) +
        Font(wb, heightInPoints = 16,
             isBold = TRUE)
      
      addDataFrame(as.data.frame(state$df4),
                   sheet = ws,
                   row.names = F,
                   startRow = 2)
      rows <- createRow(ws, rowIndex = 1)
      sheetTitle <- createCell(rows, colIndex = 1)
      setCellValue(sheetTitle[[1, 1]], paste0(state$df3$Project,
                                              " ",
                                              state$df3$Compound,
                                              " ",
                                              "Standard Curve"))
      setCellStyle(sheetTitle[[1, 1]], title_style)
      
      addPicture(file = "C:/Permeation/png/curve.png",
                 sheet = ws, 
                 startColumn = 8)
      
      #calculation process####
      startRow = 20
      startRow_h = 18
      chart_unit = 22
      data_name <- names(state$Area_long_Conc[c(4:7, 9)])
      title_style <- CellStyle(wb) +
        Font(wb, heightInPoints = 16,
             isBold = TRUE)
      for (n in 1:length(state$End_data)) {
        addDataFrame(as.data.frame(state$End_data[[n]]),
                     sheet = ws, 
                     row.names = F,
                     startRow = startRow)
        startRow = startRow + chart_unit
        
        rows <- createRow(ws, rowIndex = startRow_h + (chart_unit * (n - 1)))
        sheetTitle <- createCell(rows, colIndex = 1)
        setCellValue(sheetTitle[[1, 1]], paste0(state$df3$Project,
                                                " ",
                                                state$df3$Compound,
                                                " ",
                                                data_name[n]))
        setCellStyle(sheetTitle[[1, 1]], title_style)
        
      }
      
      #Avg & plot####
      startRow = 3
      startRow_h = 1
      chart_unit = 22
      var_name <- names(state$T_F[-c(1, 3)])
      title_style <- CellStyle(wb) +
        Font(wb, heightInPoints = 16,
             isBold = TRUE)
      for (z in c("SKY", "SKN")) {
        for (y in seq_along(var_name)) {
          for (x in 1:length(state$Avg_data)) {
            ad <- state$Avg_data[[z]][[x]][[y]]
            addDataFrame(as.data.frame(ad),
                         sheet = wa,
                         row.names = F,
                         startRow = startRow,
                         startColumn = 13)
            startRow = startRow + chart_unit
            
            rows <- createRow(wa, rowIndex = startRow_h)
            startRow_h = startRow_h + chart_unit
            sheetTitle <- createCell(rows, colIndex = 13)
            setCellValue(sheetTitle[[1, 1]], paste0(state$df3$Project,
                                                    " ",
                                                    state$df3$Compound,
                                                    " ",
                                                    names(state$Avg_data[[z]][x])
                                                    # names(state$Avg_data[[z]][[x]][y]),
                                                    # " ",
                                                    #names(state$Avg_data[z])
            ))
            setCellStyle(sheetTitle[[1, 1]], title_style)
          }
        }
      }  
      
      # Export Pic -----------
      pic_path <- list()
      for (sk in c("SKY_plot", "SKN_plot")) {
        for (kn in seq_along(names(state$T_F[-c(1, 3)]))) {
          for (fq in 1:length(state$Avg_plot)) {
            path_t <- paste0(names(state$Avg_plot[sk]),
                             "_",
                             names(state$Avg_plot[[sk]][fq]),
                             "_",
                             kn,
                             ".png")
            pic_path[[paste0(sk, "_", fq, "_", kn)]] <- path_t
            
            png(filename = path_t,  height=1000, width=1500, res=250, pointsize=8)
            print(state$Avg_plot[[sk]][[fq]][[kn]])
            dev.off()
            
          }
        }
      }
      
      startRow = 1
      startRow_h = 1
      chart_unit = 22
      for (path_n in seq_along(pic_path)) {
        
        addPicture(file = pic_path[[path_n]], sheet = wa, startRow = startRow)
        startRow = startRow + chart_unit
      }
      
      saveWorkbook(wb, file = file)
    }
  )
  
}

shinyApp(ui=ui, server=server)
