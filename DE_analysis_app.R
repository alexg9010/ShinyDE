# Load packages -----------------------------------------------------
library(shiny)
library(ggplot2)
library(dplyr)
library(tools)
library(stringr)
library(DT)
# 
# # Load data ---------------------------------------------------------
# load("data/DEtable.RData")

# Define UI ---------------------------------------------------------
ui <- fluidPage(
  
  # App title
  titlePanel("Volcano Plot of DESeq2 analysis"),
  
  # Sidebar layout with a input and output definitions
  sidebarLayout(
    
    # Inputs: Select variables to plot
    sidebarPanel(
      
      fileInput("file1", "Choose CSV File",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      
      # Set adjusted p-value threshold
      sliderInput(inputId = "p.adj.sign",
                  label = "Maximum Adjusted P-value:",
                  min = 0, max = 1,
                  value = 0.05),


      # Set adjusted p-value threshold
      sliderInput(inputId = "fc.sign",
                  label = "Minimum log2 Fold Change :",
                  min = 0, max = 10,
                  value = 1),
      
      # Select variable for color
      selectInput(inputId = "col", 
                  label = "Coloring",
                  choices = c("Significance" = "sig", 
                              "Group" = "group"), 
                  selected = "Significance"),
      
      
      # Show data table
      checkboxInput(inputId = "show_data",
                    label = "Show data table",
                    value = TRUE),
      
      # Horizontal line for visual separation
      hr(),
      
      # Select which types to plot
      checkboxGroupInput(inputId = "selected_type",
                         label = "Select type(s):",
                         choices = NULL,
                         selected = NULL)
    ),
    
    # Output:
    mainPanel(
      
      # Show scatterplot
      plotOutput(outputId = "scatterplot"),
      br(),    # a little bit of visual separation
      
      # Print number of obs plotted
      uiOutput(outputId = "n"),
      br(),    # a little bit of visual separation
      
      # Print plot description
      uiOutput(outputId = "plot_description_text"),
      br(),    # a little bit of visual separation

      # Show data table
      DT::dataTableOutput(outputId = "de_table")
    )
  )
)

# Define server function --------------------------------------------
server <- function(input, output,session) {
  
  data <- reactive({
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    infile <- input$file1
    if(is.null(infile)) {return(NULL)}
    read.csv(infile$datapath,row.names = 1)
  })
  
  groups <- reactive({
    # groups are derived from input file with column groups
    req(data())
    
    levels(data()$group)
    
  })
  
  observe({

    if(length(groups())) { 
      updateCheckboxGroupInput(session , inputId = "selected_type",  choices = groups(),
                             selected = groups()[1])
    }
  })
  
  # Create a subset of data filtering for selected types
  DE_subset <- reactive({
    req(input$selected_type) # ensure availability of value before proceeding
    results <- req(data())
    cbind(geneID=rownames(results),results) %>%
      dplyr::filter(group %in% input$selected_type)
  })

  # Create scatterplot object the plotOutput function is expecting
  output$scatterplot <- renderPlot({
    
    res <- DE_subset()
    DE_table = as.data.frame(dplyr::mutate(as.data.frame(res), 
                                          sig=ifelse(res$padj<input$p.adj.sign, 
                                                     ifelse(abs(res$log2FoldChange)>=input$fc.sign,
                                                            paste0("FDR<",input$p.adj.sign," & \nlog2(FC)>=",input$fc.sign), 
                                                            paste0("FDR<",input$p.adj.sign)), 
                                                     "Not Sig"))
                            , row.names=rownames(res))
    
    #DE_table <- na.omit(DE_table)
    
    p <- ggplot2::ggplot(DE_table, 
                    ggplot2::aes(x = log2FoldChange,y=  -log10(pvalue),
                                 label= rownames(DE_table) )) +
      ggplot2::geom_point(ggplot2::aes_string(col = input$col))  +
      ggplot2::ggtitle("Volcano Plot of DESeq2 analysis")
    
    s = input$de_table_rows_selected

    if (length(s)) p <- p + ggrepel::geom_text_repel(data=DE_table[s,],
                                      ggplot2::aes(label=DE_table[s,"geneID"]))
    p
    
    # sig.DE_table <- DE_table$sig==paste0("FDR<",p.adj.sign," & \nlog2(FC)>=",fc.sign)
    # p + ggrepel::geom_text_repel(data=DE_table[sig.DE_table,][1:10,],
    #                                   ggplot2::aes(label=rownames(DE_table[sig.DE_table,][1:10,])))
  })
  
  # Print number of genes plotted
  output$n <- renderUI({
    types <- DE_subset()$group %>% 
      factor(levels = input$selected_type) 
    counts <- table(types)
    
    significant <- filter(DE_subset(),
                          ( !is.na(padj) & padj <= input$p.adj.sign),
                          abs(log2FoldChange) >= input$fc.sign)

    counts.sign <- table(factor(significant$group,levels = input$selected_type))
    
    HTML(paste("There are", counts, input$selected_type, "genes in this dataset,
               with ",counts.sign," significant ones. <br>"))
  })
  
  # Print plot description
  # With eventReactive() only update plot description when 
  # action button is clicked
  output$plot_description_text <- eventReactive(
    eventExpr = input$update_plot_description,
    valueExpr = { HTML(input$plot_description) }
    )
  
  # Print data table if checked
  output$de_table <- DT::renderDataTable(
    if(input$show_data){
      
      DE_table <- DE_subset()
      DE_table <- DE_table[( !is.na(DE_table$padj) & DE_table$padj<=input$p.adj.sign) &
                             abs(DE_table$log2FoldChange)>=input$fc.sign,]
      DT::datatable(data = DE_table,   
                    extensions = c('Buttons', 'FixedColumns'), 
                    options = list(fixedColumns = TRUE, 
                                   scrollX = TRUE,
                                   dom = 'Bfrtip',
                                   buttons = list(c('csv','excel'),I('colvis'))),
                    rownames = FALSE,
                    filter = 'bottom'
      )
    }
  )
}

# Create the Shiny app object ---------------------------------------
shinyApp(ui, server
         ,options = list(launch.browser=TRUE)
         )

