library(shiny)
library(shinydashboard)
library(plotly)
library(readr)

sidebar <- dashboardSidebar(
  hr(),
  sidebarMenu(id="tabs",
              menuItem("EDA", tabName = "EDA", icon=icon("table"),
                       menuSubItem("Basic info", tabName = "Basic", icon = icon("angle-right")),
                       menuSubItem("Heterogeneity", tabName = "Heterogeneity", icon = icon("angle-right"))),
              menuItem("Analyses", tabName="Analyses", icon=icon("line-chart"), selected=TRUE),
              menuItem("About", tabName = "About", icon = icon("laugh-squint"))
  )
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "Basic",
            fluidRow(box( title ="ScRNA-seq data", width = NULL, status = "primary",solidHeader = TRUE,
                          DT::dataTableOutput("scRNA"))),
            fluidRow(box(  title ="Summary of scRNA-seq data ", width = NULL,
                           background = "navy", solidHeader = TRUE,  
                           tabsetPanel(
                             tabPanel(verbatimTextOutput("summaryR"))
                           ))
                     
            ),
            fluidRow(box( title ="An example of DNA methylation data", width = NULL, status = "primary",solidHeader = TRUE,
                          DT::dataTableOutput("DNAm"))),
            fluidRow(box( title ="Summary of DNA methylation data", width = NULL,
                          background = "navy", solidHeader = TRUE,  
                          tabsetPanel(
                            tabPanel(verbatimTextOutput("summaryM"))
                          )))
            
    ),
    tabItem(tabName = "Heterogeneity",
            fluidRow(
              column(width = 6,
                     box( title ="Heatmap of scRNA-seq data", width = NULL,
                          background = "navy", solidHeader = TRUE,
                          plotlyOutput("Plotheat")
                     )),
              column(width = 6,
                     box(  title ="Means and standard deviations of ESCs", width = NULL,
                           background = "navy", solidHeader = TRUE,
                           plotlyOutput("Plotmsd")
                     ))
            ),
            fluidRow(
              column(width = 6,
                     box( title ="Zero-inflated effects", width = NULL,
                          background = "navy", solidHeader = TRUE,
                          plotlyOutput("Plotcol")
                     )),
              column(width = 6,
                     box(  title ="Zero-inflated effects", width = NULL,
                           background = "navy", solidHeader = TRUE,
                           plotlyOutput("Plotrow")
                     ))
            )
    ),
    tabItem(tabName = "Analyses",
            fluidRow(
              column(width = 5,
                     box(title ="BIC", width = NULL, status = "primary",solidHeader = TRUE,
                         tabPanel(h5("BIC"),
                                  sliderInput("K", "Values of clustering parameter:", value=4, min=1, max = 4, step=1)
                         )
                     ),
                     box( title ="Clusters", width = NULL, status = "primary",solidHeader = TRUE,
                          sliderInput("clu", "Number of clusters:",
                                      min = 2, max = 5, value = 2, step= 1)
                     )),
              column(width = 7, 
                     box( title ="BIC values", width = NULL, background = "navy", solidHeader = TRUE,
                          plotlyOutput("PlotBIC"))
              )),
            fluidRow(
              column(width = 6,
                     box( title ="Genomic feature indicators", width = NULL,
                          solidHeader = TRUE, background = "navy",
                          plotlyOutput("PlotInd"))),
              column(width = 6,
                     box(  title ="The clustering result", width = NULL,
                           background = "navy", solidHeader = TRUE,
                           plotlyOutput("Plotres"))
              ))
            
    ),
    tabItem(tabName = "About",
            includeMarkdown("About.Rmd")
    )
  )
)

ui <- dashboardPage(
  dashboardHeader(title = "Integrative single-cell clustering analysis", titleWidth = 500),
  sidebar,
  body
)
