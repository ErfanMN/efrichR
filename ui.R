library(shiny)
library(DT)

options(shiny.sanitize.errors = TRUE)


# Define UI for miles per gallon app ----
ui <- pageWithSidebar(
  
  # App title ----
  headerPanel("eFrichR"),
  
  # Sidebar panel for inputs ----
  sidebarPanel(
    textAreaInput("TextArea" , "Your Gene Set:" , height = 400 ),
    br(),
    actionButton("goButton", "Go!" , icon = icon("bar-chart-o") ,  width = 250  ),
    p("Click the button to update the value displayed in the main panel.")
    
  ),
  
  # Main panel for displaying outputs ----
  mainPanel(
    plotOutput("plot"),
    DT::dataTableOutput("mytable") , 
    htmlOutput("myImage")
  )
)



