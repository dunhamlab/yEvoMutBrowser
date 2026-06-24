protein_prediction_ui <- function(id) {
  ns <- NS(id)
  tabPanel(
           "Protein Prediction",
           div(style = "margin-left: 20px; width: 600px;",
           div("", style = "height: 10px;"),
           tags$div(class = "info_div",
             verbatimTextOutput(NS(id, "info")),
           ),
            
            selectInput(NS(id, "geneSelectDropDown"), "Gene", choices = NULL),

            tags$div(style = "font-size: 12px; color: #666; margin-bottom: 8px;",
                     "Click mutations to select them (they turn black and are ",
                     "listed below); click again to remove. Available once a ",
                     "single isolate is selected (instructor, year and sample) ",
                     "or you upload your own data."),

            # Protein track: one amino acid per codon (3 bases), aligned to and
            # zoom-synced with the DNA track below. Drawn by sequence-track.js.
            tags$div(
              id = NS(id, "prottrack-wrap"),
              style = paste(
                "position: relative; width: 100%; height: 34px;",
                "border: 1px solid #ddd; border-radius: 4px;",
                "user-select: none;"
              ),
              tags$canvas(
                id = NS(id, "prottrack-canvas"),
                style = "width: 100%; height: 100%; display: block;"
              ),
              tags$div(
                id = NS(id, "prottrack-tip"),
                style = paste(
                  "position: absolute; display: none; pointer-events: none;",
                  "background: rgba(0,0,0,0.8); color: #fff; padding: 2px 6px;",
                  "border-radius: 3px; font: 12px monospace; white-space: nowrap;",
                  "z-index: 10;"
                )
              )
            ),

            div("", style = "height: 4px;"),

            # IGV-style sequence track: a canvas drawn by static/sequence-track.js.
            # Scroll to zoom, drag to pan; genomic coordinates show on hover.
            tags$div(
              id = NS(id, "seqtrack-wrap"),
              style = paste(
                "position: relative; width: 100%; height: 95px;",
                "border: 1px solid #ddd; border-radius: 4px;",
                "user-select: none;"
              ),
              tags$canvas(
                id = NS(id, "seqtrack-canvas"),
                style = "width: 100%; height: 100%; display: block;"
              ),
              tags$div(
                id = NS(id, "seqtrack-tip"),
                style = paste(
                  "position: absolute; display: none; pointer-events: none;",
                  "background: rgba(0,0,0,0.8); color: #fff; padding: 2px 6px;",
                  "border-radius: 3px; font: 12px monospace; white-space: nowrap;",
                  "z-index: 10;"
                )
              )
            ),

            div("", style = "height: 6px;"),

            # Mutations track: genomic coordinates synced with the sequence
            # track above. Markers are drawn by static/sequence-track.js.
            tags$div(
              id = NS(id, "muttrack-wrap"),
              style = paste(
                "position: relative; width:100%; height: 135px;",
                "border: 1px solid #ddd; border-radius: 4px;",
                "user-select: none;"
              ),
              tags$canvas(
                id = NS(id, "muttrack-canvas"),
                style = "width: 100%; height: 100%; display: block;"
              ),
              tags$div(
                id = NS(id, "muttrack-tip"),
                style = paste(
                  "position: absolute; display: none; pointer-events: none;",
                  "background: rgba(0,0,0,0.85); color: #fff; padding: 4px 7px;",
                  "border-radius: 3px; font: 12px monospace; white-space: pre;",
                  "z-index: 10;"
                )
              )
            ),

            uiOutput(NS(id, "mutation_legend")),

            # Box listing the mutations the user has clicked; filled by
            # static/sequence-track.js. Hidden until the first selection.
            tags$div(
              id = NS(id, "mut-selected"),
              style = paste(
                "display: none; margin-top: 12px; padding: 8px 10px;",
                "border: 1px solid #ccc; border-radius: 4px;",
                "max-width: 600px; font: 12px monospace;"
              )
            ),

            tags$script(src = "static/sequence-track.js"),

           ))
}

protein_prediction_server <- function(id, total_spaces, filtered_data, genes_info, link, gene_info_link_function, color_vector, alphafold_color_vector, form_complete = reactive(FALSE)) {
  moduleServer(id, function(input, output, session) {
    output$info <- renderText({
      "This section provides a visual representation of the predicted protein structure based on the gene sequence. The reference sequence is derived from S288c."
    })

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

    gene_sequence <- read.csv("data/gene_sequences.csv")

    # --- Translation (for the protein track) -------------------------------
    codon_table <- c(
      TTT="F", TTC="F", TTA="L", TTG="L", CTT="L", CTC="L", CTA="L", CTG="L",
      ATT="I", ATC="I", ATA="I", ATG="M", GTT="V", GTC="V", GTA="V", GTG="V",
      TCT="S", TCC="S", TCA="S", TCG="S", CCT="P", CCC="P", CCA="P", CCG="P",
      ACT="T", ACC="T", ACA="T", ACG="T", GCT="A", GCC="A", GCA="A", GCG="A",
      TAT="Y", TAC="Y", TAA="*", TAG="*", CAT="H", CAC="H", CAA="Q", CAG="Q",
      AAT="N", AAC="N", AAA="K", AAG="K", GAT="D", GAC="D", GAA="E", GAG="E",
      TGT="C", TGC="C", TGA="*", TGG="W", CGT="R", CGC="R", CGA="R", CGG="R",
      AGT="S", AGC="S", AGA="R", AGG="R", GGT="G", GGC="G", GGA="G", GGG="G")

    revcomp <- function(s) {
      comp <- chartr("ACGTN", "TGCAN", toupper(s))
      paste(rev(strsplit(comp, "", fixed = TRUE)[[1]]), collapse = "")
    }

    # Translate a coding (5'->3') sequence, residue 1 first ('*' = stop).
    translate_protein <- function(cds) {
      cds <- toupper(cds)
      ncod <- nchar(cds) %/% 3
      if (ncod < 1) return("")
      aa <- vapply(seq_len(ncod), function(j) {
        codon <- substr(cds, (j - 1) * 3 + 1, (j - 1) * 3 + 3)
        a <- codon_table[[codon]]
        if (is.null(a)) "X" else a
      }, character(1))
      paste(aa, collapse = "")
    }

    gene <- reactive({input$geneSelectDropDown})

    cur_gene <- reactive({
      # merge & filter gene info and mutation info and filter for selected gene, using the stablized data
      fd <- stable_data() %>%
        merge(genes_info, by = intersect(names(stable_data()), names(genes_info))) %>%
        arrange(GENE)

      filter(fd, GENE == gene())
    })

    # All mutations in view for this gene, one row per (POS, ALT). Always
    # shown; whether they can be *selected* is gated separately (see
    # `selectable` below).
    variant_rows <- reactive({
      cg <- cur_gene()
      if (is.null(cg) || nrow(cg) == 0) return(NULL)
      cg %>%
        group_by(POS, ALT) %>%
        summarize(
          REF = first(REF),
          ANNOTATION = first(ANNOTATION),
          proteins = paste(unique(na.omit(PROTEIN)), collapse = ", "),
          .groups = "drop"
        ) %>%
        arrange(POS)
    })

    # Push the reference sequence + this isolate's mutations to the canvas
    # renderer (static/sequence-track.js). Re-runs when the gene or the
    # upstream filter changes. The JS does all drawing and selection.
    observeEvent(list(gene(), stable_data(), form_complete()), {
      req(gene())

      group_id <- session$ns("proPred")
      seq_ids <- list(
        groupId  = group_id,
        canvasId = session$ns("seqtrack-canvas"),
        tipId    = session$ns("seqtrack-tip")
      )
      mut_ids <- list(
        groupId    = group_id,
        canvasId   = session$ns("muttrack-canvas"),
        tipId      = session$ns("muttrack-tip"),
        clickInput = session$ns("mut_click"),
        boxId      = session$ns("mut-selected")
      )
      prot_ids <- list(
        groupId  = group_id,
        canvasId = session$ns("prottrack-canvas"),
        tipId    = session$ns("prottrack-tip")
      )

      send_empty <- function(message) {
        session$sendCustomMessage("seqTrackData",
          c(seq_ids, list(empty = TRUE, message = message)))
        session$sendCustomMessage("protTrackData",
          c(prot_ids, list(empty = TRUE, message = "")))
        session$sendCustomMessage("mutTrackData",
          c(mut_ids, list(empty = TRUE, message = "")))
      }

      matched_row <- gene_sequence[gene_sequence$GENE == gene(), ]
      if (nrow(matched_row) < 1) {
        send_empty(paste0("No reference sequence available for ", gene(), "."))
        return()
      }

      matched_row <- matched_row[1, ]  # guard against duplicate rows
      fwd <- toupper(as.character(matched_row$sequence))
      gstart <- as.numeric(matched_row$start)  # genomic coord of forward base 0
      chrom <- as.character(matched_row$chrom)

      if (is.na(fwd) || !nzchar(fwd)) {
        send_empty(paste0("Sequence for ", gene(), " is empty."))
        return()
      }

      strand <- tryCatch(
        as.numeric(genes_info$STRAND[genes_info$GENE == gene()][1]),
        error = function(e) 1)
      if (is.na(strand)) strand <- 1

      # Present everything in the gene's 5'->3' (coding) orientation so residue
      # 1 is always on the left. For -1 strand the displayed sequence is the
      # reverse complement; `step` lets the JS map display index -> genomic
      # coordinate (which then counts down for -1 genes).
      L <- nchar(fwd)
      seq_str <- if (strand == -1) revcomp(fwd) else fwd
      disp_start <- if (strand == -1) gstart + L - 1 else gstart  # genomic @ idx 0
      step <- if (strand == -1) -1 else 1

      session$sendCustomMessage("seqTrackData", c(seq_ids, list(
        seq   = seq_str,
        start = disp_start,
        step  = step,
        chrom = chrom
      )))

      # Protein track: one amino acid per codon, aligned to the DNA above.
      session$sendCustomMessage("protTrackData", c(prot_ids, list(
        start     = disp_start,
        step      = step,
        chrom     = chrom,
        seqLength = L,
        aa        = translate_protein(seq_str)
      )))

      # Show all mutations for the gene; selecting (turning black + recording)
      # is only enabled when scoped to a single isolate.
      muts <- variant_rows()
      if (is.null(muts)) {
        session$sendCustomMessage("mutTrackData",
          c(mut_ids, list(empty = TRUE, message = "No mutations for this gene.")))
        return()
      }

      annotation_colors <- setNames(color_vector, all_annotations)
      mut_list <- lapply(seq_len(nrow(muts)), function(i) {
        ann <- as.character(muts$ANNOTATION[i])
        # Show the change on the displayed (coding) strand: reverse-complement
        # the alleles for -1 genes so they match the DNA shown.
        ref <- toupper(as.character(muts$REF[i]))
        alt <- toupper(as.character(muts$ALT[i]))
        if (strand == -1) { ref <- revcomp(ref); alt <- revcomp(alt) }
        change <- paste0(ref, "→", alt)
        base <- if (ann %in% c("indel-frameshift", "indel-inframe", "transposon") ||
                    nchar(alt) != 1 || !(alt %in% c("A", "C", "G", "T"))) {
          "other"
        } else {
          alt
        }
        list(
          pos        = muts$POS[i],     # genomic; JS maps to display index
          alt        = as.character(muts$ALT[i]),
          annotation = ann,
          change     = change,
          base       = base,
          color      = unname(annotation_colors[ann]),
          hover      = paste0(
            ann, "\n", change, "\n",
            if (nzchar(muts$proteins[i])) paste0(muts$proteins[i], "\n") else "",
            chrom, ":", muts$POS[i]
          )
        )
      })

      session$sendCustomMessage("mutTrackData", c(mut_ids, list(
        start      = disp_start,
        step       = step,
        chrom      = chrom,
        seqLength  = L,
        selectable = isTRUE(form_complete()),  # form filled -> clicks allowed
        mutations  = mut_list
      )))
    }, ignoreInit = FALSE)

    all_annotations <- c(
      "missense", "nonsense", "5'-upstream",
      "indel-frameshift", "indel-inframe", "synonymous", "transposon")

    three_letter_aa <- c("Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln",
                         "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe",
                         "Pro", "Ser", "Thr", "Trp", "Tyr", "Val")

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


})
}