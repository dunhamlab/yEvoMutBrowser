gene_pro_view_ui <- function(id) {
  ns <- NS(id)
  prefix <- ns("")

  tabPanel(
    "Gene View2", div("", style = "height: 10px;"),
    plotlyOutput(NS(id, "geneViewPlot"), width = "600px"), verbatimTextOutput(NS(id, "gene")),
    selectInput(NS(id, "geneSelectDropDown"), "Gene", choices = NULL),
    selectInput(NS(id, "uniprotid"), "Gene", choices = c("Q02486", "P37898", "P38631")),

    tags$div(
      id = "molstar-parent",
      style = "position: relative; width: 600px; height: 600px;",
      tags$canvas(
        id    = "molstar-canvas",
        style = "position: absolute; top: 0; left: 0; width: 100%; height: 100%;"
      ),
    )
    ,

    tags$script(HTML(sprintf("
      window.MY_MODULE_NS = '%s';
    ", prefix))),

    # ⬇️  place the script LAST so Shiny is ready
    tags$script(src = "static/molstar-custom.js"),
    verbatimTextOutput(NS(id, "resiinfo")),
    plotlyOutput(NS(id, "mutplot"), width = "600px", height = "150px"),
    div(style = "position: relative; top: 100px;",),
    plotlyOutput(NS(id, "domainplot"), width = "600px", height = "150px"),

    includeCSS("R/www/styling.css"),

    tags$div(class = "container",
      # Inner div with class "header"
      tags$div(class = "protein-cont",
        tags$div(class = "my-verbatim-text",  # <-- your custom class
          textOutput(NS(id, "content"))
        )
      ),
      # Content section
      tags$div(class = "content_div",
        tags$div(class = 'url_div',
          tags$p("Protein: ", uiOutput(NS(id, "url"), inline = TRUE)),
        ),      
        tags$div(class = "description_div",
          tags$p("Description: ", textOutput(NS(id, "description"), inline = TRUE))
        ),

        tags$div(class = "pathway_div",
          tags$p("Pathways: ", uiOutput(NS(id, "pathway"), inline = TRUE))
        ),


        
        )
),

    # uiOutput(NS(id, "url"))
   )
}

gene_pro_view_server <- function(id, total_spaces, filtered_data, genes_info, link, gene_info_link_function, color_vector) {
  moduleServer(id, function(input, output, session) {
    #gene view draopdown menu
    print(rainbow(3))
    observe({
      choices <- filtered_data() %>% pull(GENE)

      # Filter choices to include only those present in genes_info
      choices <- choices[choices %in% genes_info$GENE]

      choices <- sort(choices)

      updateSelectizeInput(session, "geneSelectDropDown",
                           choices = as.character(choices),
                           server = TRUE,
                           options = list(maxOptions = length(choices))
      )
    })

    gene <- reactive(input$geneSelectDropDown)


    cur_gene <- reactive({
      # merge & filter gene info and mutation info and filter for selected gene
      fd <- filtered_data() %>%
        merge(genes_info, by = intersect(names(filtered_data()), names(genes_info))) %>%
        arrange(GENE)

      filter(fd, GENE == gene())
    })

    output$content <- renderText("Protein Info Card")

    output$description <- renderText(paste0(first(cur_gene()$DESCRIPTION)))


    all_annotations <- c(
      "missense", "nonsense", "5'-upstream",
      "indel-frameshift", "indel-inframe", "synonymous")

    domains_info <- read.csv(ORGANISM_DOMAIN_INFO_PATH)
    pathway_info <- read.csv(ORGANISM_PATHWAY_INFO_PATH)


    output$pathway <- renderUI({
      validate(
        need(gene(), "LOADING"),
        need(filtered_data(), "Loading data...")
      )

      sgd <- first(cur_gene()$SGDID)
      row <- pathway_info[pathway_info$SGDID == sgd, ]
      BASE_URL <- "https://pathway.yeastgenome.org/YEAST/new-image?type=PATHWAY&object="

      if (nrow(row) == 0) {
                           tags$p("N/A")} else {

        links <- lapply(seq_len(nrow(row)), function(i) {
          url <- paste0(BASE_URL, row$PathwayID[i])

          tags$div(
            tags$a(href = url, target = "_blank", row$Pathway[i])
          )
        })

    # return a container with all links
    do.call(tags$div, links)
      }


  #     ifelse(nrow(row) == 0,
  #     paste0("N/Aa"),
  #     paste0(row$Pathway)
  #  )

    })

    genedatatable <- function(cur_gene) {
      pattern <- "(?<=\\d)([A-Za-z]|\\*|indel)$|([A-Za-z]|\\*)$"

      count_proteins <- cur_gene %>%
        mutate(indel = nchar(ALT) - nchar(REF)) %>% # Calculate indel difference
        group_by(POS, PROTEIN) %>%
        summarize(
          GENE = first(GENE),
          PROTEIN = first(PROTEIN),
          ANNOTATION = first(ANNOTATION),
          COUNTS = n(),
          Letter1 = substr(PROTEIN, 1, 1), # Extract the first character
          # Extract the numbers
          Numbers = as.numeric(str_extract(PROTEIN, "[0-9]+")),
          Letter2 = str_extract(PROTEIN, pattern),
          indel = first(indel)
        ) %>%
        ungroup()

      # print(count_proteins)
      count_proteins_same <- count_proteins %>%
        group_by(Numbers) %>%
        summarize(
          GENE = first(GENE),
          PROTEIN = list(PROTEIN),
          ANNOTATION = first(ANNOTATION),
          Counts_diff_mutation = list(COUNTS),
          Counts_tot = sum(COUNTS),
          indel = first(indel)
        ) %>%
        ungroup()

      all_annotations <- c(
        "missense", "nonsense", "5'-upstream",
        "indel-frameshift", "indel-inframe", "synonymous"
      )

      count_proteins_same <- count_proteins_same %>%
        mutate(
          combined = map2_chr(
            PROTEIN, Counts_diff_mutation,
            function(prot, counts) {
              prot_list <- str_split(prot, ", ") %>% unlist()
              counts_list <- str_split(counts, ", ") %>% unlist()

              combined_strings <- map2_chr(
                prot_list, counts_list,
                function(p, c) {
                  letter1 <- substr(p, 1, 1)
                  # numbers <- str_extract(p, "[0-9]+") %>% as.numeric()
                  letter2 <- str_extract(p, pattern)
                  paste("Count", letter1, "->", letter2, ":", c, "\n",
                    sep = " "
                  )
                }
              )
              paste(combined_strings, collapse = "")
            }
          )
        ) %>%
        mutate(
        PROTEIN = as.character(PROTEIN),
        # Extract the first character Amino Acid Wild Type
        AA_WT = substr(PROTEIN, 1, 1),
        AA_POS = if_else(ANNOTATION == "5'-upstream", -15,
            # Extract Amino Acid Position
            as.numeric(str_extract(PROTEIN, "[0-9]+"))
        ),
        # Amino Acid Mutation
        AA_M = substr(PROTEIN, nchar(PROTEIN), nchar(PROTEIN)),
        ANNOTATION = factor(ANNOTATION, levels = all_annotations)
        )


      count_proteins_same
    }


    # Learn about Gene button within gene viewer
    if(link != "NONE") {
    output$url <- renderUI({
      url <- a(gene(),
               href = gene_info_link_function(genes_info, input$geneSelectDropDown),
                target = "_blank"
      )
      url
    })
    }

  #   observeEvent(input$uniprotid, {
  #     uniprot_id <- input$uniprotid
  #     print(paste("Uniprot ", uniprot_id))
  #     session$sendCustomMessage("initMolstar", uniprot_id)
  # })

    output$resiinfo <- renderText({
      req(input$resiinfo)
      paste("Hovered residue:", input$resiinfo)
    })


    observeEvent(input$geneSelectDropDown, {
      gene_name <- input$geneSelectDropDown
      print(colnames(cur_gene()))
      uniprotid <- genes_info$UniprotID[genes_info$GENE == gene_name]
      print(paste("uniprot ", uniprotid))
      session$sendCustomMessage("initMolstar", uniprotid)
  })


    pfam_rectangles <- reactive({
      cg <- cur_gene()
      # Domains found based on UniparcID of the selected gene
      pd <- domains_info[domains_info$UniparcID == cg$UniparcID, ]
      pf   <- pd[pd$Database == "Pfam", ]
      if (nrow(pf) == 0) {
        return(tibble(xmin=double(), xmax=double(), 
                      ymin=double(), ymax=double(), 
                      text=character(),id=character()))
      }
      rects <- tibble(
        xmin = pf$Start,
        xmax = pf$End,
        ymin = 0,
        ymax = 1,
        text = paste0(pf$Start, "-", pf$End, "\n", pf$Description)
      ) 
        
      rects$id <- paste0(first(cg$GENE), "_", seq_len(nrow(rects)), "_", rainbow(nrow(rects)))
      rects$fill_color <- rainbow(nrow(rects))
      rects

    })

      annotation_colors <- set_names(color_vector, all_annotations)



    output$geneViewPlot <- renderPlotly({
      ranges <- reactiveValues(x = NULL, y = NULL)

      # Showing please select gene message
      # Construct the string with spaces on each side of the loading message
      geneview_message <- "Please select a gene below"
      # Calculate the number of empty spaces needed on each side
      gene_message_length <- nchar(geneview_message)
      gene_spaces_on_each_side <- floor(
        (total_spaces - gene_message_length) / 2
      )
      select_gene_message <- sprintf(
        "%s%s%s%s",
        "\n\n\n\n\n\n\n\n",
        paste(rep(" ", gene_spaces_on_each_side), collapse = ""),
        geneview_message,
        paste(rep(" ", gene_spaces_on_each_side), collapse = "")
      )

      # Render gene view plot

      validate(
        need(gene(), select_gene_message),
        need(filtered_data(), "Loading data...")
      )

      cur_gene <- cur_gene()
      count_proteins_same <- genedatatable(cur_gene)
      # Pattern for extracting the second part of protein


      annotation_colors <- set_names(color_vector, all_annotations)

      # Generating ranges for the plot size
      xmax <- genes_info %>%
        filter(GENE == gene()) %>%
        pull(PROTEIN_LENGTH) %>%
        unique() %>%
        as.numeric()
      ranges$x <- c(-50, xmax + 50)

      ranges$y <- c(0, max(count_proteins_same$Counts_tot) +
        ifelse(max(count_proteins_same$Counts_tot) < 5, 4, 1))

      p <- count_proteins_same %>%
        ggplot(
          aes(
            x = AA_POS,
            y = Counts_tot,
          )
        ) +
        # dummy geom_point to get legend showing all_annotations regardless of
        # what data is displayed
        geom_point(data = data.frame(
          PROTEIN = rep(NA_character_, length(all_annotations)),
          AA_WT = rep(NA_character_, length(all_annotations)),
          AA_POS = rep(NA_real_, length(all_annotations)),
          AA_M = rep(NA_character_, length(all_annotations)),
          ANNOTATION = factor(all_annotations, levels = all_annotations),
          Counts_diff_mutation = I(rep(
            list(numeric(0)),
            length(all_annotations)
          )),
          Counts_tot = rep(NA_integer_, length(all_annotations)),
          combined = rep(NA_character_, length(all_annotations))
        ), aes(y = 0, color = ANNOTATION, text = NULL), size = 2) +
        geom_hline(yintercept = 0, linetype = 2, alpha = .2, aes(text = NULL)) +
        geom_segment(aes(x = 0, xend = xmax, y = 0, yend = 0, text = NULL),
                     size = 15, color = "cornflowerblue"
        ) +
        geom_segment(aes(
          x = AA_POS, xend = AA_POS, y = 0, yend = Counts_tot,
          text = NULL
        ), color = "gray") +
        geom_point(aes(
          x = AA_POS,
          y = Counts_tot,
          color = ANNOTATION,
          text = ifelse(
            grepl("indel", ANNOTATION), # Check if ANNOTATION contains "indel"
            paste0(
              "Indel\n", abs(indel), " base ", ifelse(indel > 0,
                                                      "insertion", "deletion"
              ), "\nCount: ", Counts_diff_mutation,
              "\nPosition: ", AA_POS
            ), # Text for indel annotations
            ifelse(
              is.na(PROTEIN),
              paste0(
                ANNOTATION, "\nCount: ", Counts_diff_mutation,
                "\nPosition: ", -abs(unique(cur_gene$POS) -
                                       unique(cur_gene$START))
              ),
              paste0(combined, "\nPosition: ", AA_POS)
            )
          )
        ), size = 2) +
        ggtitle(as.character(gene())) +
        theme_classic() +
        theme(
          axis.text.x = element_text(),
          axis.title.y = element_text(),
          axis.text.y = element_text(),
          axis.ticks.y = element_line(),
          plot.title = element_text(hjust = 0.5),
          plot.margin = margin(20, 20, 0, 20)
        ) +
        xlab("Amino acid position") +
        ylab("Mutation Count") +
        scale_y_continuous(breaks = function(x) {
          seq(floor(min(x)),
              ceiling(max(x)),
              by = 1
          )
        }) +
        scale_color_manual(values = annotation_colors) + # Manual color scale
        coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE) +
        annotate("text",
                 x = 1, y = Inf,
                 label = "Drag over mutations to see more", hjust = 0, vjust = 2,
                 color = "black", size = 5
        ) +
        guides(color = guide_legend(title = "Annotation"))

      ggplotly(p, tooltip = "text")
    })

    # print(input$geneSelectDropDown)

    output$mutplot <- renderPlotly({

      cur_gene <- cur_gene()
      count_proteins_same <- genedatatable(cur_gene)
      cg <- cur_gene()
      print(first(cur_gene$PROTEIN_LENGTH))


      # print(count_proteins_same$combined)
      amino_acid_loc <- count_proteins_same$Numbers
      data_df <- data.frame(value = amino_acid_loc, dummy_y = 0) 

      annotation_colors <- set_names(color_vector, all_annotations)

      new_theme_empty <- theme_bw()
      new_theme_empty$line <- element_blank()
      new_theme_empty$rect <- element_blank()
      new_theme_empty$strip.text <- element_blank()
      new_theme_empty$axis.text.y <- element_blank()
      new_theme_empty$axis.title <- element_blank()
      new_theme_empty$plot.margin <- structure(c(0, 0, -1, -1),
                                              unit = "lines",
                                              valid.unit = 3L,
                                              class = "unit")

      
      # Create a scatter plot with ggplot2
      gg <- ggplot(data_df, aes(x = value, y = 0, key=count_proteins_same$ANNOTATION, color = count_proteins_same$ANNOTATION,
        text = ifelse(
          grepl("indel", count_proteins_same$ANNOTATION), # Check if ANNOTATION contains "indel"
          paste0(
            "Indel\n", abs(count_proteins_same$indel), " base ", ifelse(count_proteins_same$indel > 0,
                                                    "insertion", "deletion"
            ), "\nCount: ", count_proteins_same$Counts_diff_mutation,
            "\nPosition: ", count_proteins_same$AA_POS
          ), # Text for indel annotations
          ifelse(
            is.na(count_proteins_same$PROTEIN),
            paste0(
              count_proteins_same$ANNOTATION, "\nCount: ", count_proteins_same$Counts_diff_mutation,
              "\nPosition: ", -abs(unique(cur_gene$POS) -
                                      unique(cur_gene$START))
            ),
            paste0(count_proteins_same$combined, "\nPosition: ", count_proteins_same$AA_POS)
            )
          )
      )) +
        scale_color_manual(values = annotation_colors) + geom_point(size = 2.2, alpha = 0.3 ) +
        # geom_point(color = "#5b91e2", size = 2.2) + # Add points with specified color and size
        labs(
             x = "X-axis Variable",             # Add x-axis label
             y = "Y-axis Variable") +  xlim(0, first(cur_gene()$PROTEIN_LENGTH)) +

        new_theme_empty + labs(x = "My X-axis Label") + guides(fill = FALSE) + theme(legend.position = "none")



      gg <- ggplotly(gg, tooltip = "text", source = "mutplot", dynamicTicks = TRUE) %>%
      layout(
          xaxis = list(
            rangemode = "nonnegative",

            range = c(0, first(cg$PROTEIN_LENGTH)),   # <— force initial range [0, protein_length]
            autorange = FALSE         # <— prevent Plotly from auto-expanding it
          ) 
        ) %>%
      
        config(scrollZoom = TRUE)

      gg
    })


  observe({
    ed <- event_data("plotly_click", source="mutplot", priority = 'event')
    req(ed)
    print("DSDSD")
    ann <- as.character(ed$key[[1]])
    print(ann)
    hex <- annotation_colors[ann]
    print(hex)
    print(typeof(hex))
    hex_val <- unname(hex)[1]
    session$sendCustomMessage("highlightResidueWithSphere",
                              list(positions = ed$x,
                              colorHex=hex_val))
  })

    output$domainplot <- renderPlotly({

      new_theme_empty <- theme_bw()
      new_theme_empty$line <- element_blank()
      new_theme_empty$rect <- element_blank()
      new_theme_empty$strip.text <- element_blank()
      new_theme_empty$axis.text.y <- element_blank()
      new_theme_empty$axis.title <- element_blank()
      new_theme_empty$plot.margin <- structure(c(0, 0, -1, -1),
                                              unit = "lines",
                                              valid.unit = 3L,
                                              class = "unit")


      rects <- pfam_rectangles()
      print(rects$id)
      print(rects)
      cg <- cur_gene()

      if (nrow(rects) == 0) {
        return(plotly_empty(type="scatter", mode="markers") %>%
                 layout(title = "No Pfam domains"))
      }


      p <- ggplot(rects) +
        geom_rect(aes(xmin = xmin, xmax = xmax,
                    ymin = ymin, ymax = ymax,
                    text = text,
                    key = id,
                    fill = fill_color,
                    customdata = fill_color
                  ),

                  color = "black", 
                  alpha = 0.6) + scale_fill_identity() +
        labs(title = paste("Pfam Domains for", gene()), x = "Residue", y = "") +
        theme_minimal() + scale_x_continuous(limits = c(0, first(cg$PROTEIN_LENGTH))) + new_theme_empty +
        guides(fill = FALSE)


  gg <- ggplotly(p, tooltip = "text", dynamicTicks = TRUE, source = "domainplot") %>%
    layout(
      xaxis = list(
        rangemode = "nonnegative",

        range = c(0, first(cg$PROTEIN_LENGTH)),    # <— force initial range [0, protein_length]
        autorange = FALSE         # <— prevent Plotly from auto-expanding it
      ) 
    ) %>% config(scrollZoom = TRUE)

  gg
      
    })

    observe({
      # priortiy event will be triggered at every event/click
      ed <- event_data("plotly_click", source = "domainplot", priority = 'event')
      req(ed$key)

      # print(typeof(ed$key))

      str <- ((ed$key)[[1]])
      split_vector <- unlist(strsplit(str, "_"))
      hex <- tail(split_vector, 1)
      # req(ed$fill)
      # print(ed$fill)

      message("Clicked color: ", ed$customdata)
      domain <- pfam_rectangles()[pfam_rectangles()$id == ed$key, ]
      domain_start <- domain$xmin
      domain_end <- domain$xmax

      
      session$sendCustomMessage("highlightDomains",
      list(residueStart=domain_start, residueEnd=domain_end, colorHex=hex))
    })

  })
}