library(shiny)
library(ggplot2)
library(plotly)
library(dplyr)
library(stringr)
library(ggpubr)

#rsconnect::deployApp('.', appName = "GR-volcano-plot") 

DEResults<-read.csv("DE_results.csv")

contrasts <- colnames(DEResults)
contrastsLFC <- contrasts[grepl("log2\\.fold\\.change", contrasts)]
contrastsPval <- contrasts[grepl("p\\.adj", contrasts)]

get_base_contrast <- function(x) {
  gsub("_log2\\.fold\\.change", "", x)
}

base_contrasts <- unique(sapply(contrastsLFC, get_base_contrast))






#Shiny

ui <- fluidPage(
  titlePanel(NULL),  # Remove default title
  
  # Custom Title and Metadata Block
  tags$div(
    style = "text-align: center; margin-bottom: 30px;",
    tags$h2("Interactive Figure: Cross-Tissue Analysis of Glucocorticoid Receptor Interactome by label-free DIA-NN-RIME (DIANNeR)"),
    tags$img(src = "Figures/Abstract.png", height = "300px", style = "margin-top: 10px;"),
    tags$p(
      tags$strong("Weiye Zhao, Susanna Rose, Thomas F Grimes, et al.; Andrew Holding*")
    ),
    tags$p(
      em("York Biomedical Research Institute, Department of Biology, University of York. Additional affiliations in manuscript."),
      style = "font-size: 14px;"
    ),
    tags$p(
      "Correspondence: ",
      tags$a(href = "mailto:andrew.holding@york.ac.uk", "andrew.holding@york.ac.uk"),
      style = "font-size: 14px;"
    )
  ),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("contrast", "Select Contrast", choices = base_contrasts,
      selected  = "PrimaryBreastEpis_IgG_vs_PrimaryBreastEpis_GR"),
      textInput("search", "Search by protein (HGNC names) or partial matches",
                value       = "NRIP1",   
                placeholder = "e.g., UBE2L3"),
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
    pval_col <- paste0(input$contrast, "_p.adj")
    
    df <- DEResults %>%
      select(Protein.ID, Gene.Name, !!lfc_col, !!pval_col) %>%
      rename(
        log2FC = !!lfc_col,
        pval = !!pval_col
      ) %>%
      mutate(
        log2FC = -log2FC,
        negLog10P = -log10(pval),
        highlight = ifelse(str_detect(Gene.Name, regex(input$search, ignore_case = TRUE)), "Yes", "No")
      )
    
    df
  })
  
  output$volcanoPlot <- renderPlotly({
    df <- contrast_data()
    
    p <- ggplot(df, aes(x = log2FC, y = negLog10P,
                        text = paste("Gene:", Gene.Name, "<br>log2FC:", round(log2FC, 2), "<br>p-val:", signif(pval, 3)),
                        color = highlight)) +
      geom_point(alpha = 0.7) +
      scale_color_manual(values = c("Yes" = "magenta", "No" = "grey70")) +
      geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "grey30") +
      geom_hline(yintercept = 2, linetype = "dashed", color = "grey30") +
      ggpubr::theme_pubr() +
      labs(
        x = "log2 Fold Change",
        y = "-log10(adjust-p)",
        color = "Protein match",
        title = paste0(
          "Volcano Plot comparing ",
          gsub("_vs_", " (left) and ", input$contrast),
          " (right)"
        )
      ) + theme(legend.position = "none")
    ggplotly(p, tooltip = "text")  %>%
      layout(margin = list(t = 80))
  })
}

# Run app
shinyApp(ui, server)