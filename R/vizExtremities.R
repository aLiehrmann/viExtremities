#' vizExtremities
#' @description Shiny application launcher
#' @param maxMemory Maximum memory dedicated to the Shiny application (Go)
#' @export
vizExtremities <- function(maxMemory=5) {
  options(shiny.maxRequestSize = maxMemory*1000*1024^2)
  shiny::shinyApp(ui = ui(), server = server)
}