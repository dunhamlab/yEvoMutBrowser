gene_pro_view_ui <- function(id) {
  ns <- NS(id)
  prefix <- ns("")

  tabPanel(
    "Protein View",
    div(style = "margin-left: 20px; width: 600px;",

      div("", style = "height: 10px;"),
      tags$div(class = "info_div",
        verbatimTextOutput(NS(id, "info")),
      ),

      # Dropdown menu for gene selection
      selectInput(NS(id, "geneSelectDropDown"), "Gene", choices = NULL),
      actionButton(NS(id, "screenshot"), "Take Screenshot 📷", class = "button_style", width = "100%"),

      # Molstar viewer container
      tags$div(
        id = "molstar-parent",
        style = "position: relative; width: 600px; height: 600px;",
        tags$canvas(
          id    = "molstar-canvas",
          style = "position: absolute; top: 0; left: 0; width: 100%; height: 100%;"
        ),
      )
      ,

      # Creates namespace for the module to be used in JS
      tags$script(HTML(sprintf("
        window.MY_MODULE_NS = '%s';
      ", prefix))),

      tags$script(src = "static/molstar-custom.js"),

      # Residue info box
      tags$div(class = "resi_info_div",
        verbatimTextOutput(NS(id, "resiinfo")),
      ),

      uiOutput(NS(id, "mutation_legend")),

      # Mutation Display
      tags$div(class="annotation_div",
        tags$div(class="mut_button_div",
          actionButton(NS(id, 'mutations'), 'Mutations', class = "button_style"),
        ),

        tags$div(class="mut_plot_div",
          plotlyOutput(NS(id, "mutplot"), height = "140px"),
        ),
      ),

      # Domain Display
      tags$div(class="annotation_div",
        tags$div(class="domain_button_div",
          actionButton(NS(id, 'domain'), 'Pfam Domains', class = "button_style"),
        ),

        tags$div(class="domainplot_div",
          plotlyOutput(NS(id, "domainplot"), height = "65px"),
        ),
      ),

      # Motif Display
      tags$div(class="annotation_div",
        tags$div(class="mut_button_div",
          actionButton(NS(id, 'motif'), 'Motifs', class = "button_style"),  # Add class here
        ),

        tags$div(class="motif_plot_div",
          plotlyOutput(NS(id, "motifplot"), height = "65px"),
        ),
      ),

      includeCSS("R/www/styling.css"),


      # Protein Info Card
      tags$div(class = "container",
        tags$div(class = "protein-cont",

          tags$div(class = "protein_info_label",
            textOutput(NS(id, "content"))
          )
        ),

        tags$div(class = "content_div",
          tags$div(class = 'url_div',
            tags$p("Protein: ", uiOutput(NS(id, "url"), inline = TRUE)),
          ),
          tags$div(class = "description_div",
            tags$p("Description: ", textOutput(NS(id, "description"),
                                               inline = TRUE))
          ),

          tags$div(class = "pathway_div",
            tags$p("Pathways: ", uiOutput(NS(id, "pathway"), inline = TRUE))
          ),
        )

      ),
    )

  )
}

gene_pro_view_server <- function(id, total_spaces, filtered_data, genes_info, link, gene_info_link_function, color_vector) {
  moduleServer(id, function(input, output, session) {

    # Reactive values to track selected states
    mut_selected <- reactiveVal(FALSE)
    dom_selected <- reactiveVal(FALSE)
    motif_selected <- reactiveVal(FALSE)


    #gene view dropdown menu
    #creating stabilization so that the selected gene holds across different sample selections  
    stable_data <- debounce(filtered_data, 150)

    last_choices <- reactiveVal(NULL)

    #generates information based on the stable data (so that it doesn't flip back and forth during filtering)
    observeEvent(stable_data(), {

      new_choices <- tryCatch({
        stable_data() %>%
          dplyr::pull(GENE) %>%
          intersect(genes_info$GENE) %>%
          unique() %>%
          sort() %>%
          as.character()
      }, error = function(e) character(0))

      if (identical(new_choices, last_choices()))
        return()

      #holds onto selected gene
      current_sel <- isolate(input$geneSelectDropDown)

      selected_gene <- if (!is.null(current_sel) &&
                             nzchar(current_sel) &&
                             current_sel %in% new_choices) {
        current_sel
      } else if (length(new_choices) > 0) {
        new_choices[1]
      } else {
        NULL
      }

      updateSelectizeInput(
        session,
        "geneSelectDropDown",
        choices = new_choices,
        selected = selected_gene,
        server = TRUE,
        options = list(maxOptions = length(new_choices))
      )

      last_choices(new_choices)

    }, ignoreInit = FALSE)

    gene <- reactive({input$geneSelectDropDown})

    output$info <- renderText({
      req(gene())
      paste0("Note: The yeast protein structures are derived from AlphaFold prediction models and represent the canonical form of the protein."
      )
    })

    cur_gene <- reactive({
      # merge & filter gene info and mutation info and filter for selected gene, using the stablized data
      fd <- stable_data() %>%
        merge(genes_info, by = intersect(names(stable_data()), names(genes_info))) %>%
        arrange(GENE)

      filter(fd, GENE == gene())
    })

    output$content <- renderText("Protein Info Card")
    output$description <- renderText(paste0(first(cur_gene()$DESCRIPTION)))



    all_annotations <- c(
      "missense", "nonsense", "5'-upstream",
      "indel-frameshift", "indel-inframe", "synonymous", "transposon")

    letter_aa <- c("A", "R", "N", "D", "C", "E", "Q", "G", "H",
      "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"
    )

    three_letter_aa <- c("Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln",
                         "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe",
                         "Pro", "Ser", "Thr", "Trp", "Tyr", "Val")

    # Map one-letter to three-letter amino acid codes
    names(three_letter_aa) <- letter_aa

    # Read in data
    pfam <- read.csv(ORGANISM_PFAM_DOMAIN_INFO_PATH)
    pathway_info <- read.csv(ORGANISM_PATHWAY_INFO_PATH)
    motif_info <- read.csv(ORGANISM_MOTIF_INFO_PATH)

    output$pathway <- renderUI({
      validate(
        need(gene(), "LOADING"),
        need(filtered_data(), "Loading data...")
      )

      # Extract information for the selected gene
      sgd <- first(cur_gene()$SGDID)
      row <- pathway_info[pathway_info$SGDID == sgd, ]
      BASE_URL <- "https://pathway.yeastgenome.org/YEAST/new-image?type=PATHWAY&object="

      if (nrow(row) == 0) {
        tags$p("N/A")
      } else {

        links <- lapply(seq_len(nrow(row)), function(i) {
          url <- paste0(BASE_URL, row$PathwayID[i])

          tags$div(
            tags$a(href = url, target = "_blank", row$Pathway[i])
          )
        })

        # return a container with all links
        do.call(tags$div, links)
      }
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
          POS = first(POS),
          START = first(START),
          STOP = first(STOP),
          STRAND = first(STRAND),
          INFO = first(INFO), # Include INFO field for transposon type
          Letter1 = substr(PROTEIN, 1, 1), # Extract the first character
          # Extract the numbers - calculate protein position for transposons
          Numbers = if_else(ANNOTATION == "transposon",
                           # Calculate protein position from genomic position
                           if_else(STRAND == 1,
                                  # Forward strand: distance from start
                                  round((POS - START + 1) / 3),
                                  # Reverse strand: distance from stop
                                  round((STOP - POS + 1) / 3)),
                           as.numeric(str_extract(PROTEIN, "[0-9]+"))),
          Letter2 = str_extract(PROTEIN, pattern),
          indel = first(indel)
        ) %>%
        ungroup()

      count_proteins_same <- count_proteins %>%
        group_by(Numbers) %>%
        summarize(
          GENE = first(GENE),
          PROTEIN = list(PROTEIN),
          ANNOTATION = first(ANNOTATION),
          Counts_diff_mutation = list(COUNTS),
          Counts_tot = sum(COUNTS),
          indel = first(indel),
          INFO = first(INFO) # Propagate INFO field for transposons
        ) %>%
        ungroup()


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
                  letter2 <- str_extract(p, pattern)
                  paste("Count", three_letter_aa[letter1], "->", three_letter_aa[letter2], ":", c, "\n",
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
    if (link != "NONE") {
      output$url <- renderUI({
        url <- a(gene(),
            href = gene_info_link_function(genes_info, input$geneSelectDropDown),
            target = "_blank"
        )
        url
      })
    }

    hover_text <- function(data) {
      ifelse(
        data$ANNOTATION == "transposon", # Check if it's a transposon
        paste0(
          "Transposon insertion\n",
          "Type: ", data$INFO, "\n",
          "Count: ", data$Counts_diff_mutation,
          "\nPosition: ", data$AA_POS
        ),
        ifelse(
          grepl("indel", data$ANNOTATION), # Check if ANNOTATION contains "indel"
          paste0(
            "Indel\n", abs(data$indel), " base ", ifelse(data$indel > 0,
                                                    "insertion", "deletion"
            ), "\nCount: ", data$Counts_diff_mutation,
            "\nPosition: ", data$AA_POS
          ), # Text for indel annotations
          paste0(data$combined, "Position: ", data$AA_POS)
        )
      )
      }

    # Listens to hover events from JS and updates the text output
    output$resiinfo <- renderText({
      req(input$resi_aa, input$resi_num)
      new <- genedatatable(cur_gene())
      residue_num <- input$resi_num
      if (residue_num %in% new$Numbers) {
        # Filtered out rows with 'N/A' values
        row <- new[!is.na(new$Numbers) & new$Numbers == residue_num, ]
        hover_text(row, TRUE)
      } else {
        paste0("Position: ", input$resi_num, " Amino Acid: ", input$resi_aa)
      }
    })

    # Observe gene selection changes and loads in
    # AlphaFold Structure into Mol* Viewer
    observeEvent(input$geneSelectDropDown, {
      print(cur_gene())

      gene_name <- input$geneSelectDropDown
      # Reset selections
      mut_selected(FALSE)
      dom_selected(FALSE)
      motif_selected(FALSE)
      runjs(sprintf("$('#%s').removeClass('active');", session$ns("mutations")))
      runjs(sprintf("$('#%s').removeClass('active');", session$ns("domain")))
      runjs(sprintf("$('#%s').removeClass('active');", session$ns("motif")))
      uniprotid <- genes_info$UniprotID[genes_info$GENE == gene_name]
      # Renders AlphaFold Structure
      session$sendCustomMessage("initMolstar", uniprotid)
  })


    # General purpose function to create rectanglular data for plotting
    rect_data <- function(annotation_data, color_func, dtype="pfam") {
      cg <- cur_gene()
      pd <- annotation_data
      if (nrow(pd) == 0) {
        return(tibble(xmin = double(), xmax = double(),
                      ymin = double(), ymax = double(),
                      text = character(), id = character()))
      }
      rects <- tibble(
        xmin = pd$Start,
        xmax = pd$End,
        ymin = 0,
        ymax = 1,
        text = paste0(pd$Start, "-", pd$End, "\n", pd$Description)
      )

      rects$id <- paste0(dtype,first(cg$GENE), "_", seq_len(nrow(rects)), "_", color_func(nrow(rects)))
      rects$fill_color <- color_func(nrow(rects))
      rects

    }

    # Create a named vector of annotation colors
    annotation_colors <- set_names(color_vector, all_annotations)

    alpha_val <- 0.4
    transp_colors <- lapply(annotation_colors, function(col) {
      adjustcolor(col, alpha.f = alpha_val)
    })

    # Create a legend for the mutation types
    output$mutation_legend <- renderUI({
      tags$div(
        style = "display: flex; flex-direction: row;",
        lapply(names(transp_colors), function(mut) {
          tags$div(
            style = "display:flex; align-items:center; 
            margin-bottom:4px; margin-top:30px;",
            tags$div(
              style = sprintf("width:20px; height:20px;
              background:%s; margin-right:5px; border:1px 
              solid #000; margin-left:5px; margin-top: 4px;",
                              transp_colors[mut])
            ),
            tags$span(mut)
          )
        })
      )
    })

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


    output$mutplot <- renderPlotly({
      cur_gene <- cur_gene()
      count_proteins_same <- genedatatable(cur_gene)
      cg <- cur_gene()

      # Create a scatter plot with ggplot2
      gg <- ggplot(count_proteins_same, aes(x = count_proteins_same$Numbers,
        y = 0, key=count_proteins_same$ANNOTATION,
        color = count_proteins_same$ANNOTATION,
        text = hover_text(count_proteins_same)
      )) +
        scale_color_manual(values = annotation_colors) +
        geom_point(size = 3.2, alpha = 0.3, shape = 15) +
        xlim(0, first(cur_gene()$PROTEIN_LENGTH)) +

        new_theme_empty + labs(x = "My X-axis Label") + 
        guides(fill = FALSE) + theme(legend.position = "none")

      # Hover text
      gg <- ggplotly(gg, tooltip = "text", source = "mutplot", dynamicTicks = TRUE) %>%
        layout(
          title = list(text = NULL),
          margin = list(l = 6, r = 6, t = 6, b = 20, pad = 0),

          xaxis = list(
            rangemode = "nonnegative",
            # Force initial range [0, protein_length]
            range = c(0, first(cg$PROTEIN_LENGTH)),
            # Prevent Plotly from auto-expanding it
            autorange = FALSE
          )
        ) %>%

        config(
          scrollZoom = TRUE,
          toImageButtonOptions = list(
            format = 'svg',
            filename = 'mutation_plot',
            height = 400,
            width = 800
          )
        )

      gg
    })

  # Displays 'tracking tool' when hovering over mutation points
  observe({
    # Get the hover event data
    ed <- event_data("plotly_hover", source = "mutplot", priority = 'event')
    req(ed)
    # Extract annotation information
    ann <- as.character(ed$key[[1]])
    hex <- annotation_colors[ann]
    hex_val <- unname(hex)[1]
    positions <- c(ed$x)
    my_rectangle <- list(
      type = "rect",
      fillcolor = "#fbff00",
      opacity = 0.3,
      x0 = ed$x - 0.5,
      x1 = ed$x + 0.5,
      y0 = 0.1,      # Use relative coordinates (0-1)
      y1 = 0.9,      # Use relative coordinates (0-1)
      xref = "x",    # x-coordinates in data space
      yref = "paper" # y-coordinates as fraction of plot area (constant size)
    )

      # Update the plotlyProxy for the other plots
      plotlyProxy("mutplot", session) %>%
        plotlyProxyInvoke("relayout", list(shapes = list(my_rectangle)))

      plotlyProxy("domainplot", session) %>%
        plotlyProxyInvoke("relayout", list(shapes = list(my_rectangle)))

      plotlyProxy("motifplot", session) %>%
        plotlyProxyInvoke("relayout", list(shapes = list(my_rectangle)))
  })

    # Highlights and zooms into mutation points when clicked
    observe({
      # Get the click event data
      ed <- event_data("plotly_click", source = "mutplot", priority = "event")
      req(ed)
      ann <- as.character(ed$key[[1]])
      hex <- annotation_colors[ann]
      hex_val <- unname(hex)[1]
      positions <- c(ed$x)
      session$sendCustomMessage("zoomToResidue",
                                list(residueNumber = positions))
      session$sendCustomMessage("highlightResidueWithSphere",
                                list(positions = positions,
                                     colorHex = hex_val))
    })

    observeEvent(input$screenshot, {
      session$sendCustomMessage("takeScreenshot", list())
    })

    # Toggle mutation highlighting
    observeEvent(input$mutations, {
      mut_selected(!mut_selected())

      if (mut_selected()) {
        runjs(sprintf("$('#%s').addClass('active');", session$ns("mutations")))
      } else {
        runjs(sprintf("$('#%s').removeClass('active');", session$ns("mutations")))
      }

      if (mut_selected()) {
        gene_info <- genedatatable(cur_gene())

        for (i in 1:nrow(gene_info)) {
          row <- gene_info[i, ]
          residue <- (row["Numbers"][[1]])
          mut_type <- row["ANNOTATION"][[1]]
          hex <- annotation_colors[mut_type]
          hex_val <- unname(hex)[1]
          session$sendCustomMessage("highlightResidueWithSphere",
                                    list(positions = residue,
                                      colorHex = hex_val
                                    ))

        }
      } else {
        session$sendCustomMessage("clearSpheres", list())
      }
    })

  highlight_all <- function(rects) {
    for (i in 1:nrow(rects)) {
      row <- rects[i, ]
      domain_start <- row['xmin'][[1]]
      domain_end <- row["xmax"][[1]]
      hex <- row["fill_color"][[1]]
      hex <- unname(hex)
      len <- nchar(hex)
      hex <- substr(hex, 1, len - 2)
      session$sendCustomMessage("highlightDomains",
                                list(residueStart = domain_start,
                                     residueEnd = domain_end, colorHex = hex))
    }
  }

    # Toggle domain highlighting
    observeEvent(input$domain, {
      dom_selected(!dom_selected())
      if (dom_selected()) {
        runjs(sprintf("$('#%s').addClass('active');", session$ns("domain")))
        runjs(sprintf("$('#%s').removeClass('active');", session$ns("motif")))
      } else {
        runjs(sprintf("$('#%s').removeClass('active');", session$ns("domain")))
      }

      if (dom_selected()) {
        motif_selected(FALSE)
        cg <- cur_gene()
        # Fetches relevant Pfam domains
        pd <- pfam[pfam$UniparcID == first(cg$UniparcID), ]
        # Transforms the data into rectangles
        rects <- rect_data(pd, viridis)
        # Highlights each domain one by one
        session$sendCustomMessage("clearPaint", list())
        if (nrow(rects) != 0) {
          highlight_all(rects)
        }
      } else {
        session$sendCustomMessage("clearPaint", list())
      }
    })

    # Toggle motif highlighting
    observeEvent(input$motif, {
      motif_selected(!motif_selected())
      if (motif_selected()) {
        runjs(sprintf("$('#%s').addClass('active');", session$ns("motif")))
        runjs(sprintf("$('#%s').removeClass('active');", session$ns("domain")))
      } else {
        runjs(sprintf("$('#%s').removeClass('active');", session$ns("motif")))
      }


      if (motif_selected()) {
        dom_selected(FALSE)
        cg <- cur_gene()
        # Fetches relevant motif information
        pd <- motif_info[motif_info$UniprotID == first(cg$UniprotID), ]
        # Transforms the data into rectangles
        rects <- rect_data(pd, viridis)
        # Highlights each motif one by one
        if (nrow(rects) != 0) {
          highlight_all(rects)
        }
      } else {
        session$sendCustomMessage("clearPaint", list())
      }
    })


    # Plotting rectangular, continuous data
    plot_rect <- function(pd, cg, dtype = "pfam", plotc="domainplot", no_data_text="No Pfam domains") {
      rects <- rect_data(pd, viridis, dtype)

      if (nrow(rects) == 0) {
        return(plotly_empty(type="scatter", mode="markers") %>%
                 layout(title = no_data_text))
      }

      p <- ggplot(rects) +
        geom_rect(aes(xmin = xmin, xmax = xmax,
                    ymin = ymin, ymax = ymax,
                    text = text,
                    key = id,
                    fill = fill_color
                  ),

                  color = "black",
                  alpha = 0.6) + scale_fill_identity() +
        theme_minimal() +
        scale_x_continuous(limits = c(0, first(cg$PROTEIN_LENGTH))) +
        new_theme_empty +
        guides(fill = FALSE)


      gg <- ggplotly(p, tooltip = "text", 
                     dynamicTicks = TRUE, source = plotc) %>%
        layout(
          title = list(text = NULL),
          margin = list(l = 6, r = 6, t = 6, b = 20, pad = 0),

          xaxis = list(
            rangemode = "nonnegative",
            # Force initial range [0, protein_length]
            range = c(0, first(cg$PROTEIN_LENGTH)),
            # Prevent Plotly from auto-expanding it
            autorange = FALSE
          )
        ) %>%
        config(
          scrollZoom = TRUE,
          toImageButtonOptions = list(
            format = 'svg',
            filename = 'domain_plot',
            height = 200,
            width = 800
          )
        )

      gg

    }

    output$domainplot <- renderPlotly({
      cg <- cur_gene()
      # Get pfam domains for the current gene
      pd <- pfam[pfam$UniparcID == first(cg$UniparcID), ]
      plot_rect(pd, cg)
    })

    output$motifplot <- renderPlotly({
      cg <- cur_gene()
      # Get motif info for the current gene
      md <- motif_info[motif_info$UniprotID == first(cg$UniprotID), ]
      plot_rect(md, cg, dtype = "motif", plotc = "motifplot",
                no_data_text = "No Motif Data Available")
    })

    highlight_specific <- function(ed, annotation_data, dtype, uniprot_needed) {
      cg <- cur_gene()
      str <- ((ed$key)[[1]])
      split_vector <- unlist(strsplit(str, "_"))
      hex <- tail(split_vector, 1)
      len <- nchar(hex)
      hex <- substr(hex, 1, len - 2)
      if (uniprot_needed) {
        identifier <- first(cg$UniprotID)
        pd <- annotation_data[annotation_data$UniprotID == first(identifier), ]
      } else {
        identifier <- first(cg$UniparcID)
        pd <- annotation_data[annotation_data$UniparcID == first(identifier), ]
      }
      rects <- rect_data(pd, viridis, dtype=dtype)
      data <- rects[rects$id == ed$key, ]
      data_start <- data$xmin
      data_end <- data$xmax
      session$sendCustomMessage("zoomToResidue",
                                list(residueNumber = data_start))
      session$sendCustomMessage("highlightDomains",
                                list(residueStart = data_start,
                                     residueEnd = data_end, colorHex = hex))
    }

    observe({
      # priortiy event will be triggered at every click
      ed <- event_data("plotly_click", source = "motifplot", priority = 'event')
      req(ed$key)
      highlight_specific(ed, motif_info, dtype="motif", uniprot_needed = TRUE)
    })


    observe({
      # priortiy event will be triggered at every event/click
      ed <- event_data("plotly_click",
                       source = "domainplot", priority = "event")
      req(ed$key)
      highlight_specific(ed, pfam, dtype = "pfam", uniprot_needed = FALSE)
    })


    shared_zoom <- function(plot_id) {
      # Function to share zoom between mutplot and domainplot
      # This function uses a reactive value to store the zoom range
      # and applies it to both plots when either plot is zoomed.

      relayout_data <- event_data("plotly_relayout", source = plot_id)

      if (!is.null(relayout_data) && 
            "xaxis.range[0]" %in% names(relayout_data) && 
            "xaxis.range[1]" %in% names(relayout_data)) {

        new_range <- c(relayout_data[["xaxis.range[0]"]],
                       relayout_data[["xaxis.range[1]"]])

        plotlyProxy("domainplot") %>%
          plotlyProxyInvoke("relayout", list("xaxis.range" = new_range))

        plotlyProxy("mutplot") %>%
          plotlyProxyInvoke("relayout", list("xaxis.range" = new_range))

        plotlyProxy("motifplot") %>%
          plotlyProxyInvoke("relayout", list("xaxis.range" = new_range))

      }
    }

    observeEvent(event_data("plotly_relayout", source = "mutplot"), {
      shared_zoom("mutplot")
    })

    observeEvent(event_data("plotly_relayout", source = "domainplot"), {
      shared_zoom("domainplot")
    })

    observeEvent(event_data("plotly_relayout", source = "motifplot"), {
      shared_zoom("motifplot")
    })
  })
}