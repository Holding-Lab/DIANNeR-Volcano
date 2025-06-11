library(shiny)
library(ggplot2)
library(plotly)
library(dplyr)
library(stringr)


#rsconnect::deployApp('.', appName = "GR-volcano-plot") 

DEResults<-read.csv("DE_results.csv")

contrasts <- colnames(DEResults)
contrastsLFC <- contrasts[grepl("log2\\.fold\\.change", contrasts)]
contrastsPval <- contrasts[grepl("p\\.val", contrasts)]

get_base_contrast <- function(x) {
  gsub("_log2\\.fold\\.change", "", x)
}

base_contrasts <- unique(sapply(contrastsLFC, get_base_contrast))


#Shiny

ui <- fluidPage(
  titlePanel("Interactive Volcano Plot Viewer"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("contrast", "Select Contrast", choices = base_contrasts),
      textInput("search", "Search Gene", placeholder = "e.g., UBE2L3"),
      helpText("Hover over points to view details.")
    ),
    
    mainPanel(
      plotlyOutput("volcanoPlot", height = "600px")
    )
  )
)

# Server
server <- function(input, output, session) {
  contrast_data <- reactive({
    lfc_col <- paste0(input$contrast, "_log2.fold.change")
    pval_col <- paste0(input$contrast, "_p.val")
    
    df <- DEResults %>%
      select(Protein.ID, Gene.Name, !!lfc_col, !!pval_col) %>%
      rename(
        log2FC = !!lfc_col,
        pval = !!pval_col
      ) %>%
      mutate(
        negLog10P = -log10(pval),
        highlight = ifelse(str_detect(Gene.Name, regex(input$search, ignore_case = TRUE)), "yes", "no")
      )
    
    df
  })
  
  output$volcanoPlot <- renderPlotly({
    df <- contrast_data()
    
    p <- ggplot(df, aes(x = log2FC, y = negLog10P,
                        text = paste("Gene:", Gene.Name, "<br>log2FC:", round(log2FC, 2), "<br>p-val:", signif(pval, 3)),
                        color = highlight)) +
      geom_point(alpha = 0.7) +
      scale_color_manual(values = c("yes" = "red", "no" = "grey")) +
      theme_minimal() +
      labs(
        x = "log2 Fold Change",
        y = "-log10(p-value)",
        title = paste("Volcano Plot for", input$contrast)
      )
    
    ggplotly(p, tooltip = "text")
  })
}

# Run app
shinyApp(ui, server)