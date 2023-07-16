server <- function(input, output, session) {
  data <- shiny::reactiveValues()
  
  #- Retrieve reads from the user-specified BAM file --------------------------#
  shiny::observeEvent(
    input$file, {
      shiny::req(input$file)
      bam_file       <- input$file$datapath
      scan_param     <- Rsamtools::ScanBamParam()
      data$all_reads <- as(GenomicAlignments::readGAlignments(
        file  = bam_file, 
        param = scan_param), "GRanges")
      shiny::updateSelectInput(
        session, 
        "chr", 
        choices  = levels(GenomicRanges::seqnames(data$all_reads)), 
        selected = levels(GenomicRanges::seqnames(data$all_reads))[[1]]
      )
    }
  )
  
  #- Retrieve annotations from the user-specified gff, gff3 or gtf file -------#
  shiny::observeEvent(
    input$annot, {
      shiny::req(input$annot)
      annot_file       <- input$annot$datapath
      data$annotations <- rtracklayer::import(annot_file)
      shiny::updateSelectInput(
        session, 
        "annot_type", 
        choices  = levels(data$annotations$type), 
        selected = levels(data$annotations$type)[[1]]
      )
      shiny::updateSelectInput(
        session, 
        "annot_lab", 
        choices  = names(GenomicRanges::mcols(data$annotations)), 
        selected = names(GenomicRanges::mcols(data$annotations))[[1]]
      )
    }
  )
  
  output$scatterplot <- shiny::renderPlot({
    shiny::req(data$all_reads)
    
    #- (3',5') counts ---------------------------------------------------------#
    sub_reads <- data$all_reads[
      as.character(GenomicRanges::seqnames(data$all_reads))==input$chr &
      GenomicRanges::start(data$all_reads)<input$prob_end & 
      GenomicRanges::start(data$all_reads)>input$prob_start-input$max_range &  
      GenomicRanges::end(data$all_reads)>input$prob_start & 
      GenomicRanges::end(data$all_reads)<input$prob_end+input$max_range & 
      as.character(GenomicRanges::strand(data$all_reads))%in%input$strand,
    ]
    
    ext_counts <- table(paste0(
      GenomicRanges::start(sub_reads),
      "_",
      GenomicRanges::end(sub_reads)
    ))
    ext_id     <- names(ext_counts)
    ext_start  <- as.integer(stringr::str_extract(ext_id, "\\d+(?=_)"))
    ext_end    <- as.integer(stringr::str_extract(ext_id, "(?<=_)\\d+"))
    ext_df     <- data.frame(
      ext_id, 
      ext_start, 
      ext_end, 
      counts = as.integer(ext_counts)
    )
    
    #- Plot -------------------------------------------------------------------#
    if(is.null(input$annot)){
      #- Annotations not provided by the user ---------------------------------#
      ggplot2::ggplot(
        ext_df, 
        ggplot2::aes(
          x     = ext_end, 
          y     = ext_start, 
          size  = counts, 
          alpha = counts
        )
      ) +
        ggplot2::geom_abline(
          slope     = 1, 
          intercept = 0, 
          color     = "grey70"
        )+
        ggplot2::geom_point() +
        ggplot2::xlab(ifelse(input$strand=="+", "3'", "5'"))+
        ggplot2::ylab(ifelse(input$strand=="+", "5'", "3'"))+
        ggplot2::coord_cartesian(
          xlim = c(min(ext_df$ext_end)-100,max(ext_df$ext_end)+100),         
          ylim = c(min(ext_df$ext_start)-100,max(ext_df$ext_start)+100)
        )+
        ggplot2::theme_bw()+
        ggplot2::theme(text = ggplot2::element_text(size=20))
    } else {
      #- Retrieve overlapped annotations --------------------------------------#
      sub_annotations <- data$annotations[
        as.character(GenomicRanges::seqnames(data$annotations))==input$chr &
        GenomicRanges::start(data$annotations)<max(GenomicRanges::end(sub_reads)) & 
        GenomicRanges::end(data$annotations)>min(GenomicRanges::start(sub_reads)) & 
        as.character(GenomicRanges::strand(data$annotations)) %in% input$strand &
        as.character(data$annotations$type)==input$annot_type,
      ]
      ggplot2::ggplot(
        ext_df, 
        ggplot2::aes(
          x     = ext_end, 
          y     = ext_start, 
          size  = counts, 
          alpha = counts
        )
      ) +
        ggplot2::geom_rect(
          ggplot2::aes(
            xmin = start, 
            xmax = end, 
            ymin = start, 
            ymax = end, 
            fill = .data[[input$annot_lab]], 
            x    = NULL, 
            y    = NULL
          ), 
          data  = data.frame(sub_annotations), 
          alpha = 0.3,
          size  = NULL
        )+
        ggplot2::geom_abline(
          slope     = 1, 
          intercept = 0, 
          color     = "grey70")+
        ggplot2::geom_point() +
        ggplot2::xlab(ifelse(input$strand=="+", "3'", "5'"))+
        ggplot2::ylab(ifelse(input$strand=="+", "5'", "3'"))+
        ggplot2::coord_cartesian(
          xlim=c(min(ext_df$ext_end)-100,max(ext_df$ext_end)+100),         
          ylim=c(min(ext_df$ext_start)-100,max(ext_df$ext_start)+100)
        )+
        ggplot2::theme_bw()+
        ggplot2::theme(text=ggplot2::element_text(size=20))
    }
  })
  
  output$download <- shiny::downloadHandler(
    filename = function() {
      paste("plot_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file, width=14, height=12)
      
      #- !!!! COPY OF RENDERPLOT OTHERWISE ERROR MESSAGE !!!! -----------------#
      print({
        shiny::req(data$all_reads)
        
        sub_reads <- data$all_reads[
          as.character(GenomicRanges::seqnames(data$all_reads))==input$chr &
          GenomicRanges::start(data$all_reads)<input$prob_end & 
          GenomicRanges::start(data$all_reads)>input$prob_start-input$max_range &  
          GenomicRanges::end(data$all_reads)>input$prob_start & 
          GenomicRanges::end(data$all_reads)<input$prob_end+input$max_range & 
          as.character(GenomicRanges::strand(data$all_reads)) %in% input$strand,
        ]
        
        ext_counts <- table(paste0(
          GenomicRanges::start(sub_reads),
          "_",
          GenomicRanges::end(sub_reads)
        ))
        ext_id     <- names(ext_counts)
        ext_start  <- as.integer(stringr::str_extract(ext_id, "\\d+(?=_)"))
        ext_end    <- as.integer(stringr::str_extract(ext_id, "(?<=_)\\d+"))
        ext_df     <- data.frame(
          ext_id, 
          ext_start, 
          ext_end, 
          counts=as.integer(ext_counts)
        )
        
        if(is.null(input$annot)){
          ggplot2::ggplot(
            ext_df, 
            ggplot2::aes(
              x     = ext_end, 
              y     = ext_start, 
              size  = counts, 
              alpha = counts)
          ) +
            ggplot2::geom_abline(
              slope     = 1, 
              intercept = 0, 
              color     = "grey70"
            )+
            ggplot2::geom_point() +
            ggplot2::xlab(ifelse(input$strand=="+", "3'", "5'"))+
            ggplot2::ylab(ifelse(input$strand=="+", "5'", "3'"))+
            ggplot2::coord_cartesian(
              xlim=c(min(ext_df$ext_end)-100,max(ext_df$ext_end)+100),         
              ylim=c(min(ext_df$ext_start)-100,max(ext_df$ext_start)+100)
            )+
            ggplot2::theme_bw()+
            ggplot2::theme(text=ggplot2::element_text(size=20))
        } else {
          sub_annotations <- data$annotations[
            as.character(GenomicRanges::seqnames(data$annotations))==input$chr &
            GenomicRanges::start(data$annotations)<max(GenomicRanges::end(sub_reads)) & 
            GenomicRanges::end(data$annotations)>min(GenomicRanges::start(sub_reads)) & 
            as.character(GenomicRanges::strand(data$annotations)) %in% input$strand &
            as.character(data$annotations$type)==input$annot_type,
          ]
          
          ggplot2::ggplot(
            ext_df, 
            ggplot2::aes(
              x     = ext_end, 
              y     = ext_start, 
              size  = counts, 
              alpha = counts)
            ) +
            ggplot2::geom_rect(
              ggplot2::aes(
                xmin = start, 
                xmax = end, 
                ymin = start, 
                ymax = end, 
                fill = .data[[input$annot_lab]], 
                x    = NULL, 
                y    = NULL
              ), 
              data  = data.frame(sub_annotations), 
              alpha = 0.3,
              size  = NULL
            )+
            ggplot2::geom_abline(
              slope     = 1, 
              intercept = 0, 
              color     = "grey70"
            )+
            ggplot2::geom_point() +
            ggplot2::xlab(ifelse(input$strand=="+", "3'", "5'"))+
            ggplot2::ylab(ifelse(input$strand=="+", "5'", "3'"))+
            ggplot2::coord_cartesian(
              xlim=c(min(ext_df$ext_end)-100,max(ext_df$ext_end)+100),         
              ylim=c(min(ext_df$ext_start)-100,max(ext_df$ext_start)+100)
            )+
            ggplot2::theme_bw()+
            ggplot2::theme(text=ggplot2::element_text(size=20))
        }
      })
      dev.off()
    }
  )
}