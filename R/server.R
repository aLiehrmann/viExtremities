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
    #browser()
    #- (3',5') counts ---------------------------------------------------------#
    ext_df <- do.call(rbind, lapply(c("+","-"), function(strand){
      sub_reads <- data$all_reads[
        as.character(GenomicRanges::seqnames(data$all_reads))==input$chr &
          GenomicRanges::start(data$all_reads)<input$prob_end & 
          GenomicRanges::start(data$all_reads)>input$prob_start-input$max_range &  
          GenomicRanges::end(data$all_reads)>input$prob_start & 
          GenomicRanges::end(data$all_reads)<input$prob_end+input$max_range &
          as.character(GenomicRanges::strand(data$all_reads))==strand,
      ]
      
      #- GenomicRanges::coverage ? --------------------------------------------#
      ext_counts <- table(paste0(
        GenomicRanges::start(sub_reads),
        "_",
        GenomicRanges::end(sub_reads)
      ))
      ext_id <- names(ext_counts)
      if (strand=="+") {
        ext_start <- as.integer(stringr::str_extract(ext_id, "\\d+(?=_)"))
        ext_end   <- as.integer(stringr::str_extract(ext_id, "(?<=_)\\d+"))
      } else {
        ext_end   <- as.integer(stringr::str_extract(ext_id, "\\d+(?=_)"))
        ext_start <- as.integer(stringr::str_extract(ext_id, "(?<=_)\\d+"))
      } 
      ext_df <- data.frame(
        ext_id, 
        ext_start,
        ext_end, 
        counts = as.integer(ext_counts)
      )
    }))
    
    #- NA ?! remove it for now ------------------------------------------------# 
    ext_df <- ext_df[!is.na(ext_df$ext_start),]
    no_ext <- nrow(ext_df)==0
    no_sub_annotations <- TRUE
    #browser()
    if (!(is.null(input$annot) | no_ext)) {
      sub_annotations <- data$annotations[
        as.character(GenomicRanges::seqnames(data$annotations))==input$chr &
          GenomicRanges::start(data$annotations)<max(ext_df$ext_end) & 
          GenomicRanges::end(data$annotations)>min(ext_df$ext_start) &
          as.character(data$annotations$type)==input$annot_type,
      ]
      if (length(sub_annotations)>0) {
        sub_annotations <- data.frame(sub_annotations)
        sub_annotations <- do.call(rbind, lapply(1:nrow(sub_annotations), function(i){
          start <- sub_annotations$start[[i]]
          end <- sub_annotations$end[[i]]
          if (sub_annotations$strand[[i]]=="+"){
            data.frame(
              x     = c(start,end),
              xend  = c(end,end),
              y     = c(start,start),
              yend  = c(start,end),
              color = factor(
                rep(sub_annotations[[input$annot_lab]][[i]],2),
                levels=unique(sub_annotations[[input$annot_lab]])
              )
            )
          } else{
            data.frame(
              x     = c(start,start),
              xend  = c(start,end),
              y     = c(start,end),
              yend  = c(end,end),
              color = factor(
                rep(sub_annotations[[input$annot_lab]][[i]],2),
                levels=unique(sub_annotations[[input$annot_lab]])
              )
            )
          }
        }))
        no_sub_annotations <- FALSE
      } 
    }
    #browser()
    if (!(no_ext | no_sub_annotations)) {
      ggplot2::ggplot(
        data = ext_df, 
        ggplot2::aes(
          x     = ext_end, 
          y     = ext_start, 
          size  = counts, 
          alpha = counts
        )
      ) +
        ggplot2::geom_point() +
        ggplot2::geom_segment(
          inherit.aes = FALSE,
          ggplot2::aes(
            x     = x, 
            xend  = xend, 
            y     = y, 
            yend  = yend, 
            color = color
          ), 
          data       = sub_annotations, 
          linewidth  = 1
        )+
        ggplot2::scale_color_discrete(name = "annotation")+
        ggplot2::geom_abline(
          slope     = 1, 
          intercept = 0, 
          color     = "grey70")+
        ggplot2::xlab("3'")+
        ggplot2::ylab("5'")+
        ggplot2::coord_cartesian(
          xlim=c(min(ext_df$ext_end)-100,max(ext_df$ext_end)+100),         
          ylim=c(min(ext_df$ext_start)-100,max(ext_df$ext_start)+100)
        )+
        ggplot2::theme_bw()+
        ggplot2::theme(text=ggplot2::element_text(size=20))
    } else if (!no_ext & no_sub_annotations) {
      ggplot2::ggplot(
        data = ext_df, 
        ggplot2::aes(
          x     = ext_end, 
          y     = ext_start, 
          size  = counts, 
          alpha = counts
        )
      ) +
        ggplot2::geom_point() +
        ggplot2::geom_abline(
          slope     = 1, 
          intercept = 0, 
          color     = "grey70")+
        ggplot2::xlab("3'")+
        ggplot2::ylab("5'")+
        ggplot2::coord_cartesian(
          xlim=c(min(ext_df$ext_end)-100,max(ext_df$ext_end)+100),         
          ylim=c(min(ext_df$ext_start)-100,max(ext_df$ext_start)+100)
        )+
        ggplot2::theme_bw()+
        ggplot2::theme(text=ggplot2::element_text(size=20))
    } else {
      ggplot2::ggplot()+
        ggplot2::xlab("3'")+
        ggplot2::ylab("5'")+
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
        #browser()
        #- (3',5') counts ---------------------------------------------------------#
        ext_df <- do.call(rbind, lapply(c("+","-"), function(strand){
          sub_reads <- data$all_reads[
            as.character(GenomicRanges::seqnames(data$all_reads))==input$chr &
              GenomicRanges::start(data$all_reads)<input$prob_end & 
              GenomicRanges::start(data$all_reads)>input$prob_start-input$max_range &  
              GenomicRanges::end(data$all_reads)>input$prob_start & 
              GenomicRanges::end(data$all_reads)<input$prob_end+input$max_range &
              as.character(GenomicRanges::strand(data$all_reads))==strand,
          ]
          
          #- GenomicRanges::coverage ? --------------------------------------------#
          ext_counts <- table(paste0(
            GenomicRanges::start(sub_reads),
            "_",
            GenomicRanges::end(sub_reads)
          ))
          ext_id <- names(ext_counts)
          if (strand=="+") {
            ext_start <- as.integer(stringr::str_extract(ext_id, "\\d+(?=_)"))
            ext_end   <- as.integer(stringr::str_extract(ext_id, "(?<=_)\\d+"))
          } else {
            ext_end   <- as.integer(stringr::str_extract(ext_id, "\\d+(?=_)"))
            ext_start <- as.integer(stringr::str_extract(ext_id, "(?<=_)\\d+"))
          } 
          ext_df <- data.frame(
            ext_id, 
            ext_start,
            ext_end, 
            counts = as.integer(ext_counts)
          )
        }))
        
        #- NA ?! remove it for now ------------------------------------------------# 
        ext_df <- ext_df[!is.na(ext_df$ext_start),]
        no_ext <- nrow(ext_df)==0
        no_sub_annotations <- TRUE
        #browser()
        if (!(is.null(input$annot) | no_ext)) {
          sub_annotations <- data$annotations[
            as.character(GenomicRanges::seqnames(data$annotations))==input$chr &
              GenomicRanges::start(data$annotations)<max(ext_df$ext_end) & 
              GenomicRanges::end(data$annotations)>min(ext_df$ext_start) &
              as.character(data$annotations$type)==input$annot_type,
          ]
          if (length(sub_annotations)>0) {
            sub_annotations <- data.frame(sub_annotations)
            sub_annotations <- do.call(rbind, lapply(1:nrow(sub_annotations), function(i){
              start <- sub_annotations$start[[i]]
              end <- sub_annotations$end[[i]]
              if (sub_annotations$strand[[i]]=="+"){
                data.frame(
                  x     = c(start,end),
                  xend  = c(end,end),
                  y     = c(start,start),
                  yend  = c(start,end),
                  color = factor(
                    rep(sub_annotations[[input$annot_lab]][[i]],2),
                    levels=unique(sub_annotations[[input$annot_lab]])
                  )
                )
              } else{
                data.frame(
                  x     = c(start,start),
                  xend  = c(start,end),
                  y     = c(start,end),
                  yend  = c(end,end),
                  color = factor(
                    rep(sub_annotations[[input$annot_lab]][[i]],2),
                    levels=unique(sub_annotations[[input$annot_lab]])
                  )
                )
              }
            }))
            no_sub_annotations <- FALSE
          } 
        }
        #browser()
        if (!(no_ext | no_sub_annotations)) {
          ggplot2::ggplot(
            data = ext_df, 
            ggplot2::aes(
              x     = ext_end, 
              y     = ext_start, 
              size  = counts, 
              alpha = counts
            )
          ) +
            ggplot2::geom_point() +
            ggplot2::geom_segment(
              inherit.aes = FALSE,
              ggplot2::aes(
                x     = x, 
                xend  = xend, 
                y     = y, 
                yend  = yend, 
                color = color
              ), 
              data       = sub_annotations, 
              linewidth  = 1
            )+
            ggplot2::scale_color_discrete(name = "annotation")+
            ggplot2::geom_abline(
              slope     = 1, 
              intercept = 0, 
              color     = "grey70")+
            ggplot2::xlab("3'")+
            ggplot2::ylab("5'")+
            ggplot2::coord_cartesian(
              xlim=c(min(ext_df$ext_end)-100,max(ext_df$ext_end)+100),         
              ylim=c(min(ext_df$ext_start)-100,max(ext_df$ext_start)+100)
            )+
            ggplot2::theme_bw()+
            ggplot2::theme(text=ggplot2::element_text(size=20))
        } else if (!no_ext & no_sub_annotations) {
          ggplot2::ggplot(
            data = ext_df, 
            ggplot2::aes(
              x     = ext_end, 
              y     = ext_start, 
              size  = counts, 
              alpha = counts
            )
          ) +
            ggplot2::geom_point() +
            ggplot2::geom_abline(
              slope     = 1, 
              intercept = 0, 
              color     = "grey70")+
            ggplot2::xlab("3'")+
            ggplot2::ylab("5'")+
            ggplot2::coord_cartesian(
              xlim=c(min(ext_df$ext_end)-100,max(ext_df$ext_end)+100),         
              ylim=c(min(ext_df$ext_start)-100,max(ext_df$ext_start)+100)
            )+
            ggplot2::theme_bw()+
            ggplot2::theme(text=ggplot2::element_text(size=20))
        } else {
          ggplot2::ggplot()+
            ggplot2::xlab("3'")+
            ggplot2::ylab("5'")+
            ggplot2::theme_bw()+
            ggplot2::theme(text=ggplot2::element_text(size=20))
        }
      })
      dev.off()
    }
  )
}