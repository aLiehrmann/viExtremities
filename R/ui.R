ui <- function(){
  shiny::fluidPage(
    shiny::titlePanel("Visualization of Extremities"),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::fileInput(
          'file', 
          'Choose BAM file', 
          accept = c('.bam')
        ),
        shiny::fileInput(
          'annot', 
          'Choose annotation file', 
          accept = c('.gff3','.gtf','.gff')
        ),
        shiny::selectInput(
          "chr", 
          "Chromosme", 
          choices = c(), 
          selected=c()
        ),
        shiny::checkboxGroupInput(
          'strand', 
          'Strand', 
          choices = c('Forward' = '+', 'Reverse' = '-'), 
          selected = c('+', '-')
        ),
        shiny::numericInput(
          'prob_start', 
          'Start position',
          value = 1
        ),
        shiny::numericInput(
          'prob_end', 
          'End position', 
          value = 1000
        ),
        shiny::numericInput(
          'max_range', 
          'Maximum range', 
          value = 100
        ),
        shiny::selectInput(
          "annot_lab", 
          "Annotation label", 
          choices = c(), 
          selected=c()
        ),
        shiny::selectInput(
          "annot_type", 
          "Annotation type", 
          choices = c(), 
          selected=c()
        ),
        shiny::downloadButton(
          "download", 
          "Download Plot"
        )
      ),
      shiny::mainPanel(
        shiny::plotOutput("scatterplot", height = "800px")
      )
    )
  )
}