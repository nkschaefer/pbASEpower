# ASE / ASC Power Analysis — Shiny App
# Power to detect allele-specific expression (scRNA-seq) or
# allele-specific chromatin accessibility (scATAC-seq)
#
# Statistical model: Beta-binomial with Wald test across pseudobulked individuals
# See the "Statistical Model" tab for full documentation.
#
# Dependencies: install.packages(c("shiny", "ggplot2", "bslib"))

library(shiny)
library(bslib)
library(ggplot2)

# ============================================================================
# Core Statistical Functions
# ============================================================================

#' Analytical power via normal approximation to the beta-binomial Wald test.
#'
#' @param depth  Mean allele-informative reads per feature per individual
#' @param n_ind  Number of individuals (all contribute)
#' @param pi_alt Allelic fraction under H1 (> 0.5)
#' @param phi    Beta-binomial concentration (alpha + beta)
#' @param alpha  Per-test significance level (after any correction)
#' @return Power in [0, 1]
calc_power <- function(depth, n_ind, pi_alt, phi, alpha) {
  if (depth <= 0 || n_ind < 1 || pi_alt <= 0.5 || pi_alt >= 1) return(0)
  phi <- max(phi, 0.01)

  var_h0 <- 0.25 * (depth + phi) / (n_ind * depth * (1 + phi))
  ncp <- (pi_alt - 0.5) / sqrt(var_h0)
  z_crit <- qnorm(1 - alpha / 2)

  power <- pnorm(ncp - z_crit) + pnorm(-ncp - z_crit)
  min(max(power, 0), 1)
}

#' Monte Carlo power estimate using the exact beta-binomial generative model.
simulate_power <- function(depth, n_ind, pi_alt, phi, alpha, n_sim = 10000) {
  d <- max(1, round(depth))
  K <- max(1, round(n_ind))
  phi <- max(phi, 0.01)

  a <- pi_alt * phi
  b <- (1 - pi_alt) * phi

  p_mat <- matrix(rbeta(K * n_sim, a, b), nrow = K)
  y_mat <- matrix(rbinom(K * n_sim, d, p_mat), nrow = K)
  pi_hat <- colMeans(y_mat / d)
  se_h0 <- sqrt(0.25 * (d + phi) / (K * d * (1 + phi)))
  z <- (pi_hat - 0.5) / se_h0

  mean(abs(z) > qnorm(1 - alpha / 2))
}

#' Solve for the effective per-test threshold under Benjamini-Hochberg FDR
#' control, given the power function at that threshold.
#'
#' Under BH at level q with m tests, of which m1 = pi1*m are truly alternative,
#' the data-adaptive threshold tau satisfies:
#'   tau = q * E[rejections] / m
#'   E[rejections] = m1 * Power(tau) + m0 * tau
#' We solve this fixed-point equation iteratively.
calc_bh_alpha <- function(q, n_features, frac_true, depth, n_ind, pi_alt, phi,
                          max_iter = 50, tol = 1e-12) {
  m  <- n_features
  m1 <- m * frac_true
  m0 <- m - m1
  tau <- q * frac_true  # initial guess
  for (i in seq_len(max_iter)) {
    pow <- calc_power(depth, n_ind, pi_alt, phi, tau)
    R   <- m1 * pow + m0 * tau
    tau_new <- q * R / m
    if (abs(tau_new - tau) < tol) break
    tau <- tau_new
  }
  max(tau, .Machine$double.xmin)
}

#' Vectorized power calculation — all arguments can be vectors.
calc_power_vec <- function(depth, n_ind, pi_alt, phi, alpha) {
  phi <- pmax(phi, 0.01)
  var_h0 <- 0.25 * (depth + phi) / (n_ind * depth * (1 + phi))
  ncp <- (pi_alt - 0.5) / sqrt(var_h0)
  z_crit <- qnorm(1 - alpha / 2)
  power <- pnorm(ncp - z_crit) + pnorm(-ncp - z_crit)
  power[depth <= 0 | n_ind < 1 | pi_alt <= 0.5 | pi_alt >= 1] <- 0
  pmin(pmax(power, 0), 1)
}

#' Vectorized BH threshold — solves the fixed-point equation for
#' vectors of (depth, n_ind, pi_alt, phi) simultaneously.
calc_bh_alpha_vec <- function(q, m, frac_true, depth, n_ind, pi_alt, phi,
                              max_iter = 50, tol = 1e-12) {
  m1 <- m * frac_true
  m0 <- m - m1
  tau <- rep(q * frac_true, length(depth))
  for (i in seq_len(max_iter)) {
    pow <- calc_power_vec(depth, n_ind, pi_alt, phi, tau)
    R <- m1 * pow + m0 * tau
    tau_new <- q * R / m
    if (max(abs(tau_new - tau)) < tol) break
    tau <- tau_new
  }
  pmax(tau, .Machine$double.xmin)
}

#' Map user inputs to model parameters.
#'
#' Pipeline:
#'   1. raw reads -> unique molecules via library saturation model
#'   2. unique molecules -> raw depth per feature per individual
#'   3. raw depth -> allele-informative depth via usable read fraction
derive_params <- function(n_cells, n_ind, frac_rarest, reads_per_cell,
                          library_complexity, read_length,
                          heterozygosity, n_features) {
  cells_per_ind <- n_cells * frac_rarest / n_ind
  unique_per_cell <- library_complexity *
    (1 - exp(-reads_per_cell / library_complexity))
  saturation <- unique_per_cell / library_complexity
  raw_depth <- unique_per_cell * cells_per_ind / n_features
  frac_usable <- 1 - (1 - heterozygosity)^read_length
  eff_depth <- raw_depth * frac_usable
  list(cells_per_ind  = cells_per_ind,
       unique_per_cell = unique_per_cell,
       saturation      = saturation,
       raw_depth       = raw_depth,
       frac_usable     = frac_usable,
       depth           = eff_depth)
}

# ============================================================================
# Shiny Module — Assay Panel (reused for RNA-seq and ATAC-seq)
# ============================================================================

assayUI <- function(id, defaults) {
  ns <- NS(id)

  sidebarLayout(
    sidebarPanel(
      width = 4,

      h4("Experimental Design"),
      sliderInput(ns("n_cells"), "Total number of cells",
                  min = 1000, max = 500000, value = defaults$n_cells, step = 1000),
      sliderInput(ns("n_ind"), "Number of individuals",
                  min = 2, max = 100, value = defaults$n_ind),
      sliderInput(ns("frac_rarest"), "Fraction of cells belonging to target type",
                  min = 0.01, max = 1.0, value = defaults$frac_rarest, step = 0.01),

      hr(),
      h4("Sequencing Parameters"),
      numericInput(ns("reads_per_cell"),
                   "Avg sequencing reads per cell",
                   value = defaults$reads_per_cell,
                   min = 1000, max = 200000, step = 1000),
      helpText(defaults$reads_help),
      numericInput(ns("library_complexity"),
                   "Library complexity (unique molecules/cell)",
                   value = defaults$library_complexity,
                   min = 1000, max = 100000, step = 1000),
      helpText(defaults$complexity_help),
      numericInput(ns("read_length"), "Read length (bp)",
                   value = defaults$read_length,
                   min = 25, max = 300, step = 25),

      hr(),
      h4("Genetic Parameters"),
      sliderInput(ns("heterozygosity"), "Per-base heterozygosity",
                  min = 0.0001, max = 0.01, value = 0.001, step = 0.0001),
      helpText("Average probability that any given base is heterozygous.",
               "~0.001 for humans (1 in 1000 bp). Together with read length,",
               "this determines the fraction of reads carrying allelic information."),
      sliderInput(ns("threshold"),
                  "Allele-specific fraction threshold",
                  min = 0.51, max = 0.9, value = 0.68, step = 0.01),
      helpText("Effect size: allelic fraction under the alternative hypothesis.",
               "0.6 means the over-represented allele carries 60% of reads."),

      hr(),
      h4("Model Parameters"),
      sliderInput(ns("phi"),
                  HTML("Dispersion &phi; = &alpha; + &beta;"),
                  min = 1, max = 200, value = defaults$phi, step = 1),
      helpText(defaults$phi_help),
      numericInput(ns("n_features"),
                   paste("Number of", defaults$feature_label),
                   value = defaults$n_features,
                   min = 100, max = 500000, step = 100),
      numericInput(ns("alpha"),
                   "Significance level / FDR target",
                   value = 0.05, min = 1e-6, max = 0.2, step = 0.005),
      radioButtons(ns("correction"), "Multiple testing correction",
                   choices = c("None" = "none",
                               "Bonferroni (FWER)" = "bonferroni",
                               "Benjamini-Hochberg (FDR)" = "bh"),
                   selected = "bh"),
      conditionalPanel(
        condition = sprintf("input['%s'] == 'bh'", ns("correction")),
        sliderInput(ns("frac_true_ase"),
                    "Expected fraction of features with true ASE/ASC",
                    min = 0.01, max = 0.5, value = 0.1, step = 0.01),
        helpText("Required for BH only. The BH threshold depends on the",
                 "proportion of truly non-null features, which must be",
                 "assumed for power analysis. Bonferroni does not use this",
                 "parameter. Typical values: 5-30%.")
      ),

      hr(),
      h4("Power Curve"),
      selectInput(ns("sweep_var"), "X-axis variable:",
                  choices = c("Number of cells" = "n_cells",
                              "Number of individuals" = "n_ind",
                              "Fraction rarest cell type" = "frac_rarest",
                              "Reads per cell" = "reads_per_cell",
                              "Library complexity" = "library_complexity",
                              "Read length" = "read_length",
                              "Per-base heterozygosity" = "heterozygosity",
                              "AS fraction threshold" = "threshold",
                              "Dispersion (phi)" = "phi")),
      hr(),
      actionButton(ns("sim_btn"), "Verify with simulation",
                   class = "btn-outline-secondary btn-sm"),
      br(), br(),
      uiOutput(ns("sim_result"))
    ),

    mainPanel(
      width = 8,

      # --- Derived quantities ---
      fluidRow(style = "margin-bottom:12px;",
        column(3, wellPanel(style = "text-align:center;",
          h6("Cells / ind. / cell type"),
          h3(textOutput(ns("out_cells"))),
          tags$small(HTML("&nbsp;"))
        )),
        column(3, wellPanel(style = "text-align:center;",
          h6("Unique molecules / cell"),
          h3(textOutput(ns("out_unique"))),
          tags$small(style = "color:#6c757d;",
                     textOutput(ns("out_saturation")))
        )),
        column(3, wellPanel(style = "text-align:center;",
          h6("Usable read fraction"),
          h3(textOutput(ns("out_frac_usable"))),
          tags$small(HTML("&nbsp;"))
        )),
        column(3, wellPanel(style = "text-align:center;",
          h6("Allele-info. depth / ind."),
          h3(textOutput(ns("out_eff_depth"))),
          tags$small(HTML("&nbsp;"))
        ))
      ),

      # --- Warnings ---
      uiOutput(ns("warnings")),

      # --- Power result ---
      wellPanel(style = "text-align:center;",
        h5("Statistical Power"),
        htmlOutput(ns("power_display")),
        p(style = "color:#6c757d;", textOutput(ns("alpha_info")))
      ),

      # --- Power by depth multiplier ---
      wellPanel(
        h5("Power across the feature depth distribution"),
        p(style = "color:#6c757d;",
          paste0(defaults$feature_cap,
                 " vary widely in ", defaults$signal_noun,
                 ". This table shows power at multiples of the ",
                 "average allele-informative depth.")),
        tableOutput(ns("power_table"))
      ),

      # --- Power curve ---
      wellPanel(
        plotOutput(ns("power_curve"), height = "420px")
      )
    )
  )
}

assayServer <- function(id) {
  moduleServer(id, function(input, output, session) {

    alpha_adj <- reactive({
      a <- input$alpha
      corr <- input$correction
      if (corr == "bonferroni") {
        a <- a / input$n_features
      } else if (corr == "bh") {
        d <- derived()
        a <- calc_bh_alpha(input$alpha, input$n_features,
                           input$frac_true_ase, d$depth,
                           input$n_ind, input$threshold, input$phi)
      }
      max(a, .Machine$double.xmin)
    })

    derived <- reactive({
      derive_params(input$n_cells, input$n_ind, input$frac_rarest,
                    input$reads_per_cell, input$library_complexity,
                    input$read_length, input$heterozygosity,
                    input$n_features)
    })

    power_val <- reactive({
      d <- derived()
      calc_power(d$depth, input$n_ind, input$threshold, input$phi, alpha_adj())
    })

    # -- Derived quantity outputs --
    output$out_cells <- renderText(sprintf("%.0f", derived()$cells_per_ind))
    output$out_unique <- renderText({
      u <- derived()$unique_per_cell
      if (u >= 100) sprintf("%.0f", u) else sprintf("%.1f", u)
    })
    output$out_saturation <- renderText(
      sprintf("%.0f%% of library saturated", derived()$saturation * 100))
    output$out_frac_usable <- renderText(
      sprintf("%.1f%%", derived()$frac_usable * 100))
    output$out_eff_depth <- renderText({
      d <- derived()$depth
      if (d >= 10) sprintf("%.1f", d) else sprintf("%.2f", d)
    })

    # -- Warnings --
    output$warnings <- renderUI({
      d <- derived()
      msgs <- character(0)
      if (d$depth < 1)
        msgs <- c(msgs, paste(
          "Very low allele-informative depth per feature (",
          round(d$depth, 3), "reads).",
          "Power will be near zero. Consider increasing cells, reads",
          "per cell, or reducing the number of features."))
      if (input$n_ind < 3)
        msgs <- c(msgs,
          "Fewer than 3 individuals. The normal approximation may be unreliable.")
      if (d$saturation > 0.95)
        msgs <- c(msgs,
          paste0("Library is ", round(d$saturation * 100),
                 "% saturated. Additional sequencing reads yield ",
                 "diminishing returns. Consider increasing library ",
                 "complexity or sequencing more cells instead."))
      if (length(msgs) == 0) return(NULL)
      tagList(lapply(msgs, function(m) {
        div(class = "alert alert-warning", role = "alert", m)
      }))
    })

    # -- Power display --
    output$power_display <- renderUI({
      p <- power_val()
      color <- if (p < 0.5) "#dc3545" else if (p < 0.8) "#fd7e14" else "#198754"
      tags$h1(sprintf("%.1f%%", p * 100),
              style = sprintf("color:%s; font-size:56px; font-weight:700;", color))
    })

    output$alpha_info <- renderText({
      corr <- input$correction
      if (corr == "bonferroni") {
        sprintf("Bonferroni-adjusted alpha = %.2e  (%g / %d features)",
                alpha_adj(), input$alpha, input$n_features)
      } else if (corr == "bh") {
        sprintf("BH effective threshold = %.2e  (FDR = %g, %.0f%% true ASE assumed)",
                alpha_adj(), input$alpha, input$frac_true_ase * 100)
      } else {
        sprintf("Nominal alpha = %g (no multiple testing correction)", input$alpha)
      }
    })

    # -- Depth multiplier table (vectorized) --
    output$power_table <- renderTable({
      d <- derived()
      mults <- c(0.1, 0.25, 0.5, 1, 2, 5, 10)
      depths <- d$depth * mults
      aa <- alpha_adj()
      pows <- calc_power_vec(depths, input$n_ind, input$threshold, input$phi, aa)
      data.frame(
        Multiplier = paste0(mults, "x"),
        `Informative reads` = round(depths, 2),
        Power = sprintf("%.1f%%", pows * 100),
        check.names = FALSE
      )
    }, align = "ccr", width = "100%", striped = TRUE, hover = TRUE)

    # -- Power curve (vectorized) --
    output$power_curve <- renderPlot({
      sv <- input$sweep_var
      n <- 100L

      ri <- switch(sv,
        n_cells = list(
          vals = seq(1000, max(300000, input$n_cells * 2.5), length.out = n),
          lab  = "Number of cells",        cur = input$n_cells),
        n_ind = list(
          vals = seq(2, max(80, input$n_ind * 3), by = 1),
          lab  = "Number of individuals",   cur = input$n_ind),
        frac_rarest = list(
          vals = seq(0.01, 1.0, length.out = n),
          lab  = "Fraction of cells belonging to target type", cur = input$frac_rarest),
        reads_per_cell = list(
          vals = seq(1000, max(100000, input$reads_per_cell * 3), length.out = n),
          lab  = "Sequencing reads per cell", cur = input$reads_per_cell),
        library_complexity = list(
          vals = seq(1000, max(50000, input$library_complexity * 3), length.out = n),
          lab  = "Library complexity",      cur = input$library_complexity),
        read_length = list(
          vals = seq(25, 300, by = 5),
          lab  = "Read length (bp)",        cur = input$read_length),
        heterozygosity = list(
          vals = seq(0.0001, 0.01, length.out = n),
          lab  = "Per-base heterozygosity",  cur = input$heterozygosity),
        threshold = list(
          vals = seq(0.51, 0.9, length.out = n),
          lab  = "AS fraction threshold",   cur = input$threshold),
        phi = list(
          vals = seq(1, 200, by = 1),
          lab  = "Dispersion (phi)",        cur = input$phi)
      )

      nv <- length(ri$vals)

      # Vectorized derive_params: pass swept variable as vector
      dp_args <- list(n_cells = input$n_cells, n_ind = input$n_ind,
                      frac_rarest = input$frac_rarest,
                      reads_per_cell = input$reads_per_cell,
                      library_complexity = input$library_complexity,
                      read_length = input$read_length,
                      heterozygosity = input$heterozygosity,
                      n_features = input$n_features)
      if (sv %in% names(dp_args)) dp_args[[sv]] <- ri$vals
      dd <- do.call(derive_params, dp_args)

      ni <- if (sv == "n_ind")       ri$vals else rep(input$n_ind, nv)
      th <- if (sv == "threshold")   ri$vals else rep(input$threshold, nv)
      ph <- if (sv == "phi")         ri$vals else rep(input$phi, nv)

      corr <- input$correction
      if (corr == "bonferroni") {
        aa <- rep(input$alpha / input$n_features, nv)
      } else if (corr == "bh") {
        aa <- calc_bh_alpha_vec(input$alpha, input$n_features,
                                input$frac_true_ase, dd$depth, ni, th, ph)
      } else {
        aa <- rep(input$alpha, nv)
      }

      pows <- calc_power_vec(dd$depth, ni, th, ph, aa)

      df <- data.frame(x = ri$vals, y = pows)

      ggplot(df, aes(x, y)) +
        geom_line(color = "#2c3e50", linewidth = 1.1) +
        geom_hline(yintercept = 0.8, linetype = "dashed",
                   color = "#198754", alpha = 0.6) +
        geom_vline(xintercept = ri$cur, linetype = "dashed",
                   color = "#0d6efd", alpha = 0.6) +
        annotate("text", x = ri$cur, y = 0.03, label = "current",
                 color = "#0d6efd", hjust = -0.15, size = 3.8) +
        annotate("text", x = min(ri$vals), y = 0.82, label = "80% power",
                 color = "#198754", hjust = 0, size = 3.8) +
        scale_y_continuous(limits = c(0, 1),
                           labels = function(x) paste0(x * 100, "%")) +
        labs(x = ri$lab, y = "Power",
             title = paste("Power vs.", ri$lab)) +
        theme_minimal(base_size = 14) +
        theme(plot.title = element_text(face = "bold"),
              panel.grid.minor = element_blank())
    })

    # -- Simulation verification --
    observeEvent(input$sim_btn, {
      d <- derived()
      withProgress(message = "Simulating (n = 10 000)...", {
        sim <- simulate_power(d$depth, input$n_ind, input$threshold,
                               input$phi, alpha_adj(), n_sim = 10000)
      })
      ana <- power_val()
      output$sim_result <- renderUI({
        tags$small(style = "color:#6c757d;",
          sprintf("Simulation: %.1f%%  |  Analytical: %.1f%%",
                  sim * 100, ana * 100))
      })
    })
  })
}

# ============================================================================
# Documentation Tab
# ============================================================================

model_doc_panel <- function() {
  withMathJax(fluidRow(column(8, offset = 2, style = "padding: 30px 15px;",

    h2("Statistical Model Documentation"),
    hr(),
    h3("0. Disclaimer"),
    p("This tool was generated using Claude code and reviewed/updated with human input."),
    h3("1. Overview"),
    p("This tool calculates the statistical power to detect",
      strong("allele-specific expression (ASE)"), "from scRNA-seq or",
      strong("allele-specific chromatin accessibility (ASC)"), "from scATAC-seq.",
      "Raw sequencing reads are first converted to unique molecules via a",
      "library saturation model, then filtered to allele-informative reads",
      "(those overlapping heterozygous sites). The resulting allele counts are",
      "modeled with a beta-binomial distribution that captures overdispersion",
      "beyond binomial sampling noise, and evidence is combined across",
      "pseudobulked individuals using a Wald test."),

    h3("2. Generative Model"),
    p("For individual \\(i\\) at feature \\(g\\) (gene or peak):"),
    p("$$p_{ig} \\sim \\text{Beta}\\!\\left(\\pi\\,\\phi,\\;",
      "(1-\\pi)\\,\\phi\\right)$$"),
    p("$$Y_{ig} \\mid p_{ig} \\sim",
      "\\text{Binomial}(d_{ig},\\; p_{ig})$$"),
    p("Marginalising over \\(p_{ig}\\) yields a beta-binomial:"),
    p("$$Y_{ig} \\sim \\text{BetaBinomial}\\!\\left(",
      "d_{ig},\\;\\pi\\phi,\\;(1-\\pi)\\phi\\right)$$"),
    tags$table(class = "table", style = "max-width: 700px;",
      tags$tbody(
        tags$tr(tags$td("\\(Y_{ig}\\)"),
                tags$td("Alternative allele read count (pseudobulked, at het sites)")),
        tags$tr(tags$td("\\(d_{ig}\\)"),
                tags$td("Allele-informative read depth")),
        tags$tr(tags$td("\\(\\pi\\)"),
                tags$td("Population mean allelic fraction (0.5 under \\(H_0\\))")),
        tags$tr(tags$td("\\(\\phi = \\alpha + \\beta\\)"),
                tags$td("Dispersion / concentration parameter")),
        tags$tr(tags$td("\\(\\rho = 1/(\\phi+1)\\)"),
                tags$td("Intra-class correlation (derived from \\(\\phi\\))"))
      )
    ),
    p("Moments of the marginal:"),
    p("$$E[Y_{ig}] = d_{ig}\\,\\pi$$"),
    p("$$\\text{Var}(Y_{ig}) = d_{ig}\\,\\pi(1-\\pi)",
      "\\cdot \\frac{d_{ig} + \\phi}{1 + \\phi}$$"),
    p("When \\(\\phi \\to \\infty\\), the variance reduces to",
      "the binomial \\(d\\,\\pi(1-\\pi)\\). The factor",
      "\\((d+\\phi)/(1+\\phi)\\) inflates the variance to capture biological",
      "variability in allelic ratios."),

    h3("3. Library Complexity and Unique Molecules"),
    p("Raw sequencing reads are converted to unique molecular observations",
      "using a saturation model:"),
    p("$$n_{\\text{unique}} = C \\times \\left(",
      "1 - e^{-R/C}\\right)$$"),
    p("where \\(R\\) is the number of raw sequencing reads per cell and",
      "\\(C\\) is the library complexity (total unique molecules available",
      "per cell). At low sequencing depth (\\(R \\ll C\\)),",
      "most reads are unique. As depth increases, duplicate reads",
      "accumulate and new information saturates."),
    p("The", strong("saturation fraction"), "\\(n_{\\text{unique}} / C\\)",
      "indicates how much of the library has been captured:"),
    tags$ul(
      tags$li("< 70%: additional sequencing will substantially increase unique molecules"),
      tags$li("70-90%: moderate returns from additional sequencing"),
      tags$li("> 90%: diminishing returns; more cells or higher-complexity",
              "libraries would be more cost-effective")
    ),
    p("ATAC-seq libraries typically have lower per-cell complexity than",
      "RNA-seq libraries (fewer Tn5 insertion events in accessible",
      "chromatin), so they saturate faster at equivalent sequencing depth.",
      "This is why increasing ATAC-seq sequencing depth yields rapidly",
      "diminishing returns compared to RNA-seq."),

    h3("4. Allele-Informative Read Fraction"),
    p("Only reads overlapping a heterozygous site carry allelic information.",
      "Given per-base heterozygosity \\(h\\) and read length \\(L\\) bp,",
      "the probability that a read overlaps at least one heterozygous",
      "position is:"),
    p("$$f_{\\text{usable}} = 1 - (1 - h)^L$$"),
    p("For humans (\\(h \\approx 0.001\\)) with 100 bp reads:"),
    p("$$f_{\\text{usable}} = 1 - 0.999^{100} \\approx 9.5\\%$$"),

    h3("5. Hypothesis Test"),
    p("We test for a population-level allelic imbalance at each feature:"),
    p("$$H_0:\\; \\pi = 0.5 \\qquad\\text{(no allelic imbalance)}$$"),
    p("$$H_1:\\; \\pi = \\pi_1 \\qquad\\text{(allelic fraction at threshold)}$$"),
    p("Per-individual allelic fraction estimate:",
      "\\(\\hat{\\pi}_i = Y_i / d_i\\)."),
    p("Combined estimator across \\(K = N_{\\text{ind}}\\) individuals with",
      "equal depth \\(d\\):"),
    p("$$\\bar{\\pi} = \\frac{1}{K}\\sum_{i=1}^{K} \\hat{\\pi}_i$$"),
    p("Wald test statistic:"),
    p("$$Z = \\frac{\\bar{\\pi} - 0.5}{",
      "\\sqrt{\\text{Var}_0(\\bar{\\pi})}}$$"),
    p("where the null variance is:"),
    p("$$\\text{Var}_0(\\bar{\\pi}) = \\frac{0.25\\,(d+\\phi)}{",
      "K\\,d\\,(1+\\phi)}$$"),
    p("We reject \\(H_0\\) at level \\(\\alpha\\) when \\(|Z| > z_{\\alpha/2}\\)."),

    h3("6. Power Formula"),
    p("Under \\(H_1\\) the non-centrality parameter is:"),
    p("$$\\delta = \\frac{\\pi_1 - 0.5}{",
      "\\sqrt{\\text{Var}_0(\\bar{\\pi})}}",
      "= (\\pi_1 - 0.5)\\sqrt{\\frac{K\\,d\\,(1+\\phi)}{",
      "0.25\\,(d+\\phi)}}$$"),
    p("Power (two-sided):"),
    p("$$\\text{Power} = \\Phi(\\delta - z_{\\alpha/2})",
      "+ \\Phi(-\\delta - z_{\\alpha/2})$$"),
    p("where \\(\\Phi\\) is the standard normal CDF."),

    h3("7. Full Derivation Pipeline"),
    p("User inputs map to the model through the following chain:"),
    p("$$\\text{cells per ind.} = \\frac{N_{\\text{cells}} \\times",
      "f_{\\text{rare}}}{N_{\\text{ind}}}$$"),
    p("$$n_{\\text{unique}} = C \\times",
      "\\left(1 - e^{-R/C}\\right)$$"),
    p("$$\\text{raw depth} = \\frac{n_{\\text{unique}} \\times",
      "\\text{cells per ind.}}{N_{\\text{features}}}$$"),
    p("$$f_{\\text{usable}} = 1 - (1 - h)^L$$"),
    p("$$d_{\\text{effective}} = \\text{raw depth} \\times f_{\\text{usable}}$$"),
    p("$$K = N_{\\text{ind}}$$"),

    h3("8. Interpreting the Dispersion Parameter"),
    p("The parameter \\(\\phi = \\alpha + \\beta\\) controls the concentration",
      "of the Beta distribution from which individual-level allelic fractions",
      "are drawn. It can be interpreted as the",
      strong("effective number of prior reads"), ": a Beta(\\(\\alpha, \\beta\\))",
      "prior on the allelic fraction \\(\\pi\\) carries the same information",
      "as having previously observed \\(\\phi = \\alpha + \\beta\\) reads."),
    tags$ul(
      tags$li("\\(\\phi \\to \\infty\\): binomial model; all variance comes from",
              "finite depth."),
      tags$li("\\(\\phi \\approx 50\\text{-}100\\): modest overdispersion,",
              "typical for highly expressed genes in bulk RNA-seq ASE",
              "(Castel et al. 2015)."),
      tags$li("\\(\\phi \\approx 21\\): moderate overdispersion.",
              "The default for this tool, corresponding to the geometric",
              "mean of the median-dispersion range reported in",
              "Petrova et al. 2025 (ASPEN, PLoS Comp Biol;",
              "doi:10.1371/journal.pcbi.1013837)."),
      tags$li("\\(\\phi < 10\\): substantial overdispersion; additional",
              "sequencing depth provides little benefit.")
    ),
    p("At \\(\\phi = 21\\) and \\(\\pi = 0.5\\), the Beta(10.5, 10.5) prior",
      "has SD \\(\\approx 0.10\\), so 95% of individual-level allelic fractions",
      "fall in approximately [0.29, 0.71]."),
    p("The critical insight: once allele-informative depth \\(d\\) substantially",
      "exceeds \\(\\phi\\), the variance \\(\\text{Var}(\\bar{\\pi})\\)",
      "\\(\\approx 0.25 / (K \\cdot \\phi)\\) and",
      strong("only more individuals reduce uncertainty.")),

    h3("9. Multiple Testing"),
    p("Three correction modes are available:"),

    tags$h5("None"),
    p("No correction is applied; power is reported at the nominal",
      "\\(\\alpha\\). Appropriate when testing a single pre-specified feature."),

    tags$h5("Bonferroni (FWER control)"),
    p("Sets \\(\\alpha_{\\text{adj}} = \\alpha / N_{\\text{features}}\\).",
      "Controls the family-wise error rate (probability of any false positive).",
      "Conservative; provides a", strong("lower bound"), "on achievable power.",
      "Does not depend on the fraction of truly non-null features."),

    tags$h5("Benjamini-Hochberg (FDR control)"),
    p("Controls the false discovery rate at level \\(q = \\alpha\\).",
      "Less conservative than Bonferroni and more commonly used in",
      "genome-wide ASE studies."),
    p("Unlike Bonferroni, BH power analysis requires an additional",
      "assumption:", strong("the expected fraction of features with true",
      "allelic imbalance"), "(\\(\\pi_1\\)).",
      "This is because the BH threshold is data-adaptive: it depends on",
      "the empirical p-value distribution, which in turn depends on how",
      "many features are truly non-null. Bonferroni does not use this",
      "parameter."),
    p("The effective per-test threshold \\(\\tau\\) is found by solving",
      "the fixed-point equation:"),
    p("$$\\tau = q \\cdot \\frac{m_1 \\cdot \\text{Power}(\\tau)",
      "+ m_0 \\cdot \\tau}{m}$$"),
    p("where \\(m_1 = \\pi_1 m\\) is the number of truly alternative features",
      "and \\(m_0 = m - m_1\\) is the number of true nulls. The equation is",
      "solved iteratively. When \\(\\pi_1\\) is small (few true positives),",
      "\\(\\tau\\) approaches the Bonferroni threshold; when \\(\\pi_1\\) is",
      "large, \\(\\tau\\) becomes much more liberal."),

    h3("10. Assumptions and Limitations"),
    tags$ol(
      tags$li(tags$b("Average depth:"),
              "The calculation uses the mean read depth across features.",
              "Expression/accessibility is highly skewed.",
              "The depth-multiplier table addresses this by showing power",
              "at different expression quantiles."),
      tags$li(tags$b("Equal depth across individuals:"),
              "Assumes similar cell counts per individual. Unequal sampling",
              "reduces effective sample size."),
      tags$li(tags$b("Independence:"),
              "Assumes individuals are unrelated. Population structure",
              "or cryptic relatedness would inflate false positives."),
      tags$li(tags$b("Single dispersion:"),
              "Uses one \\(\\phi\\) for all features. In practice,",
              "dispersion varies by expression level and genomic context."),
      tags$li(tags$b("Normal approximation:"),
              "The Wald test is accurate when \\(K \\times d\\) is",
              "moderate-to-large but may be anti-conservative at very",
              "low depths. Use the simulation button to check."),
      tags$li(tags$b("Population-level test:"),
              "Tests whether the", em("average"),
              "allelic fraction across individuals deviates from 0.5.",
              "Individual-level ASE tests have different power properties."),
      tags$li(tags$b("Uniform heterozygosity:"),
              "The usable-read fraction assumes het sites are uniformly",
              "distributed at rate \\(h\\) per base. In practice,",
              "heterozygosity varies across the genome."),
      tags$li(tags$b("Library saturation model:"),
              "The \\(C(1-e^{-R/C})\\) model assumes each unique molecule",
              "is sampled independently. This is a standard approximation",
              "but may not perfectly capture protocol-specific effects",
              "(e.g., GC bias, UMI collisions).")
    ),

    h3("11. Parameter Guidance"),
    tags$table(class = "table table-striped",
      tags$thead(tags$tr(
        tags$th("Parameter"), tags$th("scRNA-seq"), tags$th("scATAC-seq")
      )),
      tags$tbody(
        tags$tr(tags$td("Reads per cell"),
                tags$td("20,000-50,000"),
                tags$td("25,000-50,000")),
        tags$tr(tags$td("Library complexity"),
                tags$td(HTML("~15,000<br>(mRNA capture efficiency)")),
                tags$td(HTML("~8,000<br>(limited by Tn5 insertions)"))),
        tags$tr(tags$td("Features"),
                tags$td("10,000-20,000 genes"),
                tags$td("100,000-300,000 peaks")),
        tags$tr(tags$td(HTML("Dispersion &phi;")),
                tags$td(HTML("20-100<br>(higher for well-expressed genes)")),
                tags$td(HTML("10-50<br>(ATAC tends toward lower &phi;)"))),
        tags$tr(tags$td("Per-base heterozygosity"),
                tags$td("~0.001 (human)"),
                tags$td("~0.001 (human)")),
        tags$tr(tags$td("Read length"),
                tags$td("100-150 bp (10x)"),
                tags$td("50-150 bp"))
      )
    ),

    h3("12. Simulation Verification"),
    p("The", em("Verify with simulation"), "button runs 10,000 Monte Carlo",
      "iterations of the full generative model:"),
    tags$ol(
      tags$li("Draw per-individual allelic fractions:",
              "\\(p_i \\sim \\text{Beta}(\\pi_1 \\phi,\\;(1-\\pi_1)\\phi)\\)"),
      tags$li("Draw allele counts:",
              "\\(Y_i \\sim \\text{Binomial}(d, p_i)\\)"),
      tags$li("Compute the Wald \\(Z\\) statistic"),
      tags$li("Reject if \\(|Z| > z_{\\alpha/2}\\)")
    ),
    p("Agreement validates the normal approximation.",
      "Discrepancy at very low depth or few individuals",
      "indicates the approximation is inaccurate in that regime."),

    hr(),
    p(style = "color:#6c757d;",
      "Implementation uses a closed-form Wald test power formula.",
      "No simulation packages (e.g., splatter) are required.")
  )))
}

# ============================================================================
# App UI
# ============================================================================

ui <- navbarPage(
  title = "Pseudobulk single cell ASE / ASC Power Analysis",
  theme = bs_theme(version = 5, bootswatch = "yeti"),

  tabPanel("scRNA-seq",
    icon = icon("dna"),
    assayUI("rna", list(
      n_cells            = 50000,
      n_ind              = 10,
      frac_rarest        = 0.05,
      reads_per_cell     = 20000,
      library_complexity = 15000,
      read_length        = 100,
      phi                = 21,
      n_features         = 15000,
      feature_label      = "genes tested",
      feature_unit       = "gene",
      feature_cap        = "Genes",
      signal_noun        = "expression level",
      reads_help         = paste("Raw sequencing reads per cell (before deduplication).",
                                 "For 10x Chromium scRNA-seq, typically 20,000-50,000."),
      complexity_help    = paste("Total unique molecules per cell (library size).",
                                 "Determines how quickly additional sequencing",
                                 "saturates. RNA-seq: ~15,000 (mRNA capture limits)."),
      phi_help           = paste("Beta-binomial concentration (phi = alpha + beta).",
                                 "Interpretable as effective prior reads.",
                                 "Default 21 is the geometric mean of the median-",
                                 "dispersion range reported in Petrova et al. 2025",
                                 "(ASPEN, PLoS Comp Biol). Once allele-informative",
                                 "depth exceeds phi, only more individuals help.")
    ))
  ),

  tabPanel("scATAC-seq",
    icon = icon("layer-group"),
    assayUI("atac", list(
      n_cells            = 50000,
      n_ind              = 10,
      frac_rarest        = 0.05,
      reads_per_cell     = 25000,
      library_complexity = 8000,
      read_length        = 100,
      phi                = 21,
      n_features         = 150000,
      feature_label      = "peaks tested",
      feature_unit       = "peak",
      feature_cap        = "Peaks",
      signal_noun        = "accessibility",
      reads_help         = paste("Raw sequencing read pairs per cell.",
                                 "For 10x ATAC/Multiome, typically 25,000-50,000."),
      complexity_help    = paste("Total unique fragments per cell.",
                                 "ATAC-seq: ~8,000 (limited by Tn5 insertion",
                                 "events in accessible chromatin). Lower than",
                                 "RNA-seq, causing faster saturation."),
      phi_help           = paste("Beta-binomial concentration (phi = alpha + beta).",
                                 "Interpretable as effective prior reads.",
                                 "Default 21 (Petrova et al. 2025, ASPEN).",
                                 "ATAC-seq may warrant lower values (10-20)",
                                 "due to noisier allelic ratios in sparse",
                                 "chromatin accessibility data.")
    ))
  ),

  tabPanel("Statistical Model",
    icon = icon("book"),
    model_doc_panel()
  )
)

# ============================================================================
# App Server
# ============================================================================

server <- function(input, output, session) {
  # Uncomment to enable live theme picker during development:
  # bs_themer()
  assayServer("rna")
  assayServer("atac")
}

shinyApp(ui, server)
