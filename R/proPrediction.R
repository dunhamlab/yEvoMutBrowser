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
                     "Click the mutation you want to predict (turns black and is ",
                     "listed below); click again to remove."),

            # Protein track: one amino acid per codon (3 bases), aligned to and
            # zoom-synced with the DNA track below. Drawn by sequence-track.js.
            tags$div(
              id = NS(id, "prottrack-wrap"),
              style = paste(
                "position: relative; width: 100%; height: 95px;",
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

            # X-axis label for the wild-type DNA/protein tracks above.
            tags$div(
              "Position Number",
              style = paste(
                "text-align: center; font-size: 11px; color: #666;",
                "margin-top: 2px; margin-bottom: 4px;"
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

            # Shown by static/sequence-track.js only while a mutation is
            # selected; hidden otherwise (and on gene/sample changes).
            actionButton(NS(id, "mut_click"), "Generate Mutated Sequence",
                         style = "display: none; margin-top: 10px;"),

            # Mutant DNA + protein tracks (with the selected mutation applied),
            # in the same view group as the wild type so they zoom/pan together.
            # Shown by sequence-track.js after "Generate Mutated Sequence" is clicked.
            tags$div(
              id = NS(id, "mutant-wrap"),
              style = "display: none; margin-top: 14px;",
              tags$div("With selected mutation",
                       style = "font-weight: bold; font-size: 22px; margin-bottom: 4px;"),
              tags$div(
                id = NS(id, "mutprot-wrap"),
                style = paste(
                  "position: relative; width: 100%; height: 95px;",
                  "border: 1px solid #ddd; border-radius: 4px;",
                  "user-select: none;"
                ),
                tags$canvas(id = NS(id, "mutprot-canvas"),
                            style = "width: 100%; height: 100%; display: block;"),
                tags$div(id = NS(id, "mutprot-tip"), style = paste(
                  "position: absolute; display: none; pointer-events: none;",
                  "background: rgba(0,0,0,0.8); color: #fff; padding: 2px 6px;",
                  "border-radius: 3px; font: 12px monospace; white-space: nowrap;",
                  "z-index: 10;"))
              ),
              div("", style = "height: 4px;"),
              tags$div(
                id = NS(id, "mutdna-wrap"),
                style = paste(
                  "position: relative; width: 100%; height: 95px;",
                  "border: 1px solid #ddd; border-radius: 4px;",
                  "user-select: none;"
                ),
                tags$canvas(id = NS(id, "mutdna-canvas"),
                            style = "width: 100%; height: 100%; display: block;"),
                tags$div(id = NS(id, "mutdna-tip"), style = paste(
                  "position: absolute; display: none; pointer-events: none;",
                  "background: rgba(0,0,0,0.8); color: #fff; padding: 2px 6px;",
                  "border-radius: 3px; font: 12px monospace; white-space: nowrap;",
                  "z-index: 10;"))
              ),

              # X-axis label for the mutant DNA/protein tracks above.
              tags$div(
                "Position Number",
                style = paste(
                  "text-align: center; font-size: 11px; color: #666;",
                  "margin-top: 2px; margin-bottom: 4px;"
                )
              ),

              # Copies the mutant protein (without the stop "*"); wired in
              # static/sequence-track.js. Also opens AlphaFold Server in a new
              # tab so the user can immediately paste the copied sequence in
              # to predict the mutant structure.
              tags$button(id = NS(id, "copy-protein"),
                          class = "btn btn-default", style = "margin-top: 10px;",
                          onclick = "window.open('https://alphafoldserver.com/', '_blank');",
                          "Copy mutated protein")
            ),

            tags$script(src = "static/sequence-track.js"),

            # WT-vs-mutant structure comparison. The wild-type is auto-loaded
            # from AlphaFold DB; the mutant (predicted in AlphaFold from the
            # copied protein) is uploaded here, aligned to the WT orientation
            # and shown in its own viewer. The two viewers' cameras are linked
            # so rotating one rotates the other. The mutated residue is a black
            # strip on each cartoon.
            # Hidden until the user copies the mutated protein (that's the point
            # at which they have a sequence to fold in AlphaFold); revealed by
            # static/sequence-track.js on the copy click, re-hidden on any
            # gene/mutation change.
            tags$div(id = NS(id, "struct-compare"),
              style = "display: none; margin-top: 16px; width: 640px;",
              tags$div("Structure comparison — WT (grey) vs mutant (orange)",
                       style = "font-weight: bold; font-size: 22px; margin-bottom: 4px;"),
              fileInput(NS(id, "mutant_structure"),
                        "Upload AlphaFold mutant structure (.pdb / .cif)",
                        accept = c(".pdb", ".cif"),
                        width = "420px"),
              tags$div(
                style = "display: flex; flex-direction: column; gap: 12px;",
                tags$div(
                  tags$div("Wild type", style = "font-size: 12px; color: #666;"),
                  tags$div(
                    id = NS(id, "wt-parent"),
                    style = "position: relative; width: 600px; height: 460px;",
                    tags$canvas(
                      id = NS(id, "wt-canvas"),
                      style = paste("position: absolute; top: 0; left: 0;",
                                    "width: 100%; height: 100%;")
                    )
                  ),
                  # Residue under the cursor in the WT viewer.
                  verbatimTextOutput(NS(id, "wt_resiinfo"))
                ),
                tags$div(
                  tags$div("Mutant", style = "font-size: 12px; color: #666;"),
                  tags$div(
                    id = NS(id, "mut-parent"),
                    style = "position: relative; width: 600px; height: 460px;",
                    tags$canvas(
                      id = NS(id, "mut-canvas"),
                      style = paste("position: absolute; top: 0; left: 0;",
                                    "width: 100%; height: 100%;")
                    )
                  ),
                  # Residue under the cursor in the mutant viewer.
                  verbatimTextOutput(NS(id, "mut_resiinfo"))
                )
              ),
              # Zoom either viewer to the mutation residue; the linked camera
              # makes the other follow (useful when the two structures aren't
              # perfectly aligned and the mutation drifts off-screen).
              tags$div(
                style = "margin-top: 10px; display: flex; gap: 10px;",
                actionButton(NS(id, "zoom_wt"), "Hover to WT mutation"),
                actionButton(NS(id, "zoom_mut"), "Hover to Mutant mutation")
              ),
            ),
            # Step-by-step usage instructions. Styled identically to the
            # "info_div" box above (same class -> same border/background/font);
            # only the content differs: "Instructions:" is bold.
            tags$div(class = "info_div",
             uiOutput(NS(id, "instructions")),
           ),
           ))
}

protein_prediction_server <- function(id, total_spaces, filtered_data, genes_info, link, gene_info_link_function, color_vector, alphafold_color_vector, form_complete = reactive(FALSE)) {
  moduleServer(id, function(input, output, session) {

    output$info <- renderText({
      "This section provides a visual representation of the predicted protein structure based on the gene sequence. The reference sequence is derived from S288c."
    })

    # Same visual container as output$info (both live inside a class="info_div"
    # div), but built with renderUI instead of renderText/verbatimTextOutput so
    # "Instructions:" can be bold. Built as a single HTML() string (rather than
    # passing tags$b()/strings as separate tags$pre() arguments) because
    # htmltools pretty-prints separate tag arguments with its own indentation
    #/newlines, and since the CSS sets white-space: pre-wrap on this pre, that
    # extra formatting whitespace was rendering as visibly huge gaps between
    # lines. A single pre-built string has exactly the \n's we want and none
    # we don't.
    output$instructions <- renderUI({
      tags$pre(
        class = "shiny-text-output",
        HTML(paste0(
          "<b>Instructions:</b>\n",
          "1. Select a gene from the dropdown menu.\n",
          "2. Click on the mutation of interest in the mutation track to select it (it will turn black).\n",
          "3. Click 'Generate Mutated Sequence' to view the mutated DNA and protein tracks.\n",
          "4. Copy the mutated protein sequence and use AlphaFold to predict its structure.\n",
          "5. Upload the predicted mutant structure to compare it with the wild-type structure."
        ))
      )
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

    # Translate a coding (5'->3') sequence, residue 1 first. Stops at (and
    # includes) the first stop codon ('*'), so a premature stop introduced by a
    # nonsense mutation truncates the protein there.
    translate_protein <- function(cds) {
      cds <- toupper(cds)
      ncod <- nchar(cds) %/% 3
      if (ncod < 1) return("")
      aa <- character(ncod)
      n <- 0L
      for (j in seq_len(ncod)) {
        codon <- substr(cds, (j - 1) * 3 + 1, (j - 1) * 3 + 3)
        a <- codon_table[[codon]]
        if (is.null(a)) a <- "X"
        n <- n + 1L
        aa[n] <- a
        if (a == "*") break
      }
      paste(aa[seq_len(n)], collapse = "")
    }

    # Apply a selected coding-sequence mutation to the reference sequence.
    # The mark spans the mutant sequence positions [idx, idx + nchar(alt) - 1]
    # -- i.e. the full ALT allele, anchored at the mutation's genomic
    # position -- rather than only the subset of bases that differ
    # character-by-character from REF. That keeps the highlighted region
    # consistent with what's marked in the mutation track above: a deletion
    # like ATTG -> A marks just the single remaining anchor base, and an
    # insertion like T -> TTTT marks the full 4-base inserted run.
    apply_mutation_to_sequence <- function(seq, gl, sel) {
      if (is.null(gl) || is.null(sel) || is.null(sel$pos)) {
        return(list(seq = seq, markIndices = integer(0)))
      }

      idx <- coding_index(gl, sel$pos)
      if (is.na(idx) || idx < 0 || idx >= gl$L) {
        return(list(seq = seq, markIndices = integer(0)))
      }

      ref_genomic <- toupper(as.character(sel$ref))
      alt_genomic <- toupper(as.character(sel$alt))
      ref_coding <- if (gl$strand == -1) revcomp(ref_genomic) else ref_genomic
      alt_coding <- if (gl$strand == -1) revcomp(alt_genomic) else alt_genomic

      start_pos <- idx + 1L
      end_pos <- start_pos + nchar(ref_coding) - 1L

      mutated_seq <- paste0(
        substr(seq, 1L, start_pos - 1L),
        alt_coding,
        substr(seq, end_pos + 1L, nchar(seq))
      )

      mark_indices <- integer(0)
      ann <- if (!is.null(sel$annotation)) tolower(as.character(sel$annotation)) else ""
      if (nchar(ref_coding) == 1L && nchar(alt_coding) == 1L) {
        mark_indices <- idx
      } else if (identical(ann, "indel-inframe") && nchar(alt_coding) > 0L) {
        mark_indices <- seq.int(idx, length.out = nchar(alt_coding))
      }

      list(seq = mutated_seq, markIndices = mark_indices)
    }

    gene <- reactive({input$geneSelectDropDown})

    # ---- Structure comparison (Mol*): two synced side-by-side viewers ----
    wtc <- session$ns("wt-canvas")
    mutc <- session$ns("mut-canvas")

    # Load the wild-type AlphaFold structure into the WT viewer. Used both on
    # gene change and whenever the selected mutation changes (see below) --
    # in the latter case we clear + reload rather than leaving the existing
    # structure in place, since Mol* keeps a persistent residue
    # selection/highlight (e.g. from "Hover to WT mutation") that would
    # otherwise carry the previous mutation's highlight forward.
    load_wt_structure <- function() {
      uid <- genes_info$UniprotID[genes_info$GENE == gene()][1]
      if (length(uid) == 0 || is.na(uid) || !nzchar(uid)) return()
      session$sendCustomMessage("initMolstar", list(
        canvasId = wtc,
        parentId = session$ns("wt-parent"),
        numInput = session$ns("wt_resi_num"),
        aaInput  = session$ns("wt_resi_aa"),
        source   = list(kind = "alphafold", uniprotId = uid)
      ))
    }

    # Load the wild-type AlphaFold structure into the WT viewer on gene change.
    # Also wipe the mutant viewer + its file input so a mutant uploaded for the
    # previous gene doesn't carry over into the new one, and disable "Hover to
    # Mutant mutation" since there's now nothing loaded for it to zoom to.
    observeEvent(gene(), {
      req(gene())
      session$sendCustomMessage("clearViewer", list(canvasId = mutc))
      shinyjs::reset("mutant_structure")
      shinyjs::disable("zoom_mut")
      load_wt_structure()
    })

    # Upload a mutant structure -> load it in the mutant viewer aligned to the
    # WT orientation, link the two cameras, and mark the mutation residue
    # black. initAlignedMutant on the JS side is a long, un-awaited async
    # chain (parse the upload, fetch a second AlphaFold structure just to
    # compute the alignment, superpose, commit, delete the temp reference,
    # build the cartoon, highlight, then link cameras) -- the struct-compare
    # panel is already visible by this point, so a click on "Hover to Mutant
    # mutation" can land while that chain is still running and zoomToResidue()
    # finds an empty structure hierarchy and silently no-ops. Disable the
    # button here and only re-enable it once JS reports the chain is
    # genuinely finished (see the mutant_ready observer below).
    observeEvent(input$mutant_structure, {
      f <- input$mutant_structure
      req(f)
      shinyjs::disable("zoom_mut")
      uid <- genes_info$UniprotID[genes_info$GENE == gene()][1]
      txt <- paste(readLines(f$datapath, warn = FALSE), collapse = "\n")
      fmt <- if (grepl("\\.cif$", f$name, ignore.case = TRUE)) "mmcif" else "pdb"

      # Protein residue of the selected mutation, for the black cartoon strip.
      gl <- gene_layout()
      sel <- input$mut_selection
      res <- NULL
      if (!is.null(gl) && !is.null(sel) && !is.null(sel$pos)) {
        idx <- coding_index(gl, sel$pos)
        if (!is.na(idx) && idx >= 0) {
          aa_idx <- floor(idx / 3) + 1
          ann <- if (!is.null(sel$annotation)) tolower(as.character(sel$annotation)) else ""
          if (identical(ann, "nonsense")) aa_idx <- max(1L, aa_idx - 1L)
          res <- aa_idx
        }
      }

      session$sendCustomMessage("initAlignedMutant", list(
        canvasId     = mutc,
        parentId     = session$ns("mut-parent"),
        numInput     = session$ns("mut_resi_num"),
        aaInput      = session$ns("mut_resi_aa"),
        readyInput   = session$ns("mutant_ready"),
        colorHex     = "#d55e00",
        source       = list(kind = "data", format = fmt, data = txt),
        alignUniprot = if (!is.na(uid) && nzchar(uid)) uid else NULL,
        linkWith     = wtc,
        markResidue  = res   # NULL if no mutation selected
      ))
    })

    # Fired by initAlignedMutant on the JS side once the mutant structure is
    # fully loaded, aligned, highlighted, and camera-linked -- only then is
    # "Hover to Mutant mutation" safe to use.
    observeEvent(input$mutant_ready, {
      shinyjs::enable("zoom_mut")
    }, ignoreInit = TRUE)

    # The mutation-track click handler (static/sequence-track.js) hides the
    # mutant-wrap/struct-compare containers as soon as the selected mutation
    # changes, but it has no way to reach into the Mol* viewers themselves --
    # so both the previously-uploaded mutant structure AND the WT viewer's
    # residue highlight/selection (e.g. from "Hover to WT mutation") were
    # staying in place, just hidden. Clear and reload both viewers whenever
    # the selection changes (including deselection) so the next
    # "copy mutated protein" / upload starts from a blank slate, matching the
    # gene-change reset above. Also disable "Hover to Mutant mutation" again,
    # since the mutant viewer is now empty until the next upload finishes.
    observeEvent(input$mut_selection, {
      session$sendCustomMessage("clearViewer", list(canvasId = mutc))
      shinyjs::reset("mutant_structure")
      shinyjs::disable("zoom_mut")
      session$sendCustomMessage("clearViewer", list(canvasId = wtc))
      load_wt_structure()
    }, ignoreInit = TRUE, ignoreNULL = FALSE)

    # Protein residue number of the currently selected mutation (or NULL).
    # For nonsense mutations, the introduced stop codon truncates the mutant
    # protein at that site, so we highlight the preceding residue instead.
    mut_residue <- reactive({
      gl <- gene_layout()
      sel <- input$mut_selection
      if (is.null(gl) || is.null(sel) || is.null(sel$pos)) return(NULL)
      idx <- coding_index(gl, sel$pos)
      if (is.na(idx) || idx < 0) return(NULL)
      aa_idx <- floor(idx / 3) + 1
      ann <- if (!is.null(sel$annotation)) tolower(as.character(sel$annotation)) else ""
      if (identical(ann, "nonsense")) aa_idx <- max(1L, aa_idx - 1L)
      aa_idx
    })

    # Zoom a viewer to the mutation residue. The linked camera makes the other
    # viewer follow, but zooming each independently lets the user recenter on
    # its own residue when the structures aren't perfectly aligned.
    observeEvent(input$zoom_wt, {
      res <- mut_residue(); req(res)
      session$sendCustomMessage("zoomToResidue",
        list(canvasId = wtc, residueNumber = res, chainId = "A"))
    })
    observeEvent(input$zoom_mut, {
      res <- mut_residue(); req(res)
      session$sendCustomMessage("zoomToResidue",
        list(canvasId = mutc, residueNumber = res, chainId = "A"))
    })

    # Hovered-residue readouts for each viewer (fed by the viewers' hover ->
    # Shiny inputs), mirroring the Protein View tab's residue info box.
    output$wt_resiinfo <- renderText({
      req(input$wt_resi_num, input$wt_resi_aa)
      paste0("Position: ", input$wt_resi_num, " Amino Acid: ", input$wt_resi_aa)
    })
    output$mut_resiinfo <- renderText({
      req(input$mut_resi_num, input$mut_resi_aa)
      paste0("Position: ", input$mut_resi_num, " Amino Acid: ", input$mut_resi_aa)
    })
    # The comparison section is shown/hidden via raw JS (display:none), so Shiny
    # can't tell these become visible; keep them live so they update on hover.
    outputOptions(output, "wt_resiinfo", suspendWhenHidden = FALSE)
    outputOptions(output, "mut_resiinfo", suspendWhenHidden = FALSE)

    cur_gene <- reactive({
      # merge & filter gene info and mutation info and filter for selected gene, using the stablized data
      fd <- stable_data() %>%
        merge(genes_info, by = intersect(names(stable_data()), names(genes_info))) %>%
        arrange(GENE)

      filter(fd, GENE == gene())
    })

    # Coordinate layout for the current gene, in the gene's 5'->3' (coding)
    # orientation. Shared by the wild-type and mutant tracks so they map
    # display index -> genomic coordinate identically. NULL if unavailable.
    gene_layout <- reactive({
      g <- gene()
      if (is.null(g) || !nzchar(g)) return(NULL)
      mr <- gene_sequence[gene_sequence$GENE == g, ]
      if (nrow(mr) < 1) return(NULL)
      mr <- mr[1, ]
      fwd <- toupper(as.character(mr$sequence))
      if (is.na(fwd) || !nzchar(fwd)) return(NULL)
      gstart <- as.numeric(mr$start)
      L <- nchar(fwd)
      strand <- tryCatch(as.numeric(genes_info$STRAND[genes_info$GENE == g][1]),
                         error = function(e) 1)
      if (is.na(strand)) strand <- 1
      list(
        chrom      = as.character(mr$chrom),
        gstart     = gstart,
        L          = L,
        strand     = strand,
        coding     = if (strand == -1) revcomp(fwd) else fwd,
        disp_start = if (strand == -1) gstart + L - 1 else gstart,
        step       = if (strand == -1) -1 else 1
      )
    })

    # Coding-sequence display index (0-based) of a genomic position.
    coding_index <- function(gl, pos) {
      (as.numeric(pos) - gl$disp_start) * gl$step
    }

    # All coding-sequence mutations in view for this gene, one row per
    # (POS, ALT).
    variant_rows <- reactive({
      cg <- cur_gene()
      if (is.null(cg) || nrow(cg) == 0) return(NULL)
      # Drop mutation types that don't sit on the coding sequence.
      cg <- cg %>%
        filter(!(ANNOTATION %in% c("transposon", "synonymous", "5'-upstream")))
      if (nrow(cg) == 0) return(NULL)
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
    observeEvent(list(gene(), stable_data()), {
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
        clickInput = session$ns("mut_selection"),
        boxId      = session$ns("mut-selected"),
        buttonId   = session$ns("mut_click"),
        mutantWrapId = session$ns("mutant-wrap"),
        structCompareId = session$ns("struct-compare")
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

      gl <- gene_layout()
      if (is.null(gl)) {
        send_empty(paste0("No reference sequence available for ", gene(), "."))
        return()
      }
      # Everything is presented in the gene's 5'->3' (coding) orientation
      # (see gene_layout); `step` lets the JS map display index -> genomic.
      chrom <- gl$chrom
      L <- gl$L
      strand <- gl$strand
      seq_str <- gl$coding
      disp_start <- gl$disp_start
      step <- gl$step

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
        # missense and nonsense both land on the base they change to (ALT);
        # indels (non-single-base) go to the "indel" lane.
        base <- if (ann %in% c("indel-frameshift", "indel-inframe") ||
                    nchar(alt) != 1 || !(alt %in% c("A", "C", "G", "T"))) {
          "indel"
        } else {
          alt
        }
        list(
          pos        = muts$POS[i],     # genomic; JS maps to display index
          ref        = as.character(muts$REF[i]),
          alt        = as.character(muts$ALT[i]),
          annotation = ann,
          change     = change,
          base       = base,
          clickable  = ann %in% c("missense", "nonsense", "indel-inframe"),
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
        mutations  = mut_list
      )))
    }, ignoreInit = FALSE)

    # "Generate Sequence": apply the selected (missense/nonsense/indel-inframe)
    # mutation to the coding sequence and render mutant DNA + protein tracks in
    # the same view group as the wild type. Only the DNA track gets a black
    # highlight at the affected position(s); the protein track is left
    # unhighlighted (protMarkIndices is always empty) so amino acid residues
    # are never blacked out.
    observeEvent(input$mut_click, {
      gl <- gene_layout()
      sel <- input$mut_selection  # single object: list(pos, alt), or NULL
      if (is.null(gl) || is.null(sel) || is.null(sel$pos)) return()

      mutation_result <- apply_mutation_to_sequence(gl$coding, gl, sel)
      mut_seq <- mutation_result$seq
      mark_indices <- mutation_result$markIndices

      session$sendCustomMessage("mutantTracks", list(
        groupId      = session$ns("proPred"),
        wrapId       = session$ns("mutant-wrap"),
        dnaCanvasId  = session$ns("mutdna-canvas"),
        dnaTipId     = session$ns("mutdna-tip"),
        protCanvasId = session$ns("mutprot-canvas"),
        protTipId    = session$ns("mutprot-tip"),
        copyBtnId    = session$ns("copy-protein"),
        start        = gl$disp_start,
        step         = gl$step,
        chrom        = gl$chrom,
        seqLength    = gl$L,
        markIndex    = if (length(mark_indices) > 0) mark_indices[1] else NULL,
        markIndices  = mark_indices,
        dnaMarkIndices = mark_indices,
        protMarkIndices = integer(0),
        seq          = mut_seq,
        aa           = translate_protein(mut_seq)
      ))
    }, ignoreInit = TRUE)

    all_annotations <- c(
      "missense", "nonsense", "5'-upstream",
      "indel-frameshift", "indel-inframe", "synonymous", "transposon")

    # Mutation types shown on the track (those that sit on the coding sequence).
    track_annotations <- c(
      "missense", "nonsense", "indel-frameshift", "indel-inframe")

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
        style = "display: flex; flex-direction: row; justify-content: center;",
        lapply(track_annotations, function(mut) {
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
