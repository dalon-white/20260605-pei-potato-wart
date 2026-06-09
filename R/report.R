#' Report generation

#' Write Markdown report with equations, assumptions, and results.
#' @param cfg Configuration list.
#' @param annual_summary Annual summary table.
#' @param required_tbl Required n table.
#' @param oc_tbl OC curve table.
#' @param sensitivity_tbl Sensitivity table.
#' @param upstream_tbl Upstream filtering summary table.
#' @param output_path Output markdown path.
#' @return Output path.
write_report <- function(cfg, annual_summary, required_tbl, oc_tbl, sensitivity_tbl, upstream_tbl = NULL, output_path = "outputs/report.md") {
  upstream_lines <- if (is.null(upstream_tbl)) {
    character(0)
  } else {
    c(
      "",
      "## Upstream shipment filtering (scenario summary)",
      paste(capture.output(print(upstream_tbl)), collapse = "\n")
    )
  }

  lines <- c(
    "# Border Sampling Analysis Report",
    "",
    "## Context",
    "This repository simulates shipment-level and annual risk for detecting pest A in commodity B under upstream visual screening and border asymptomatic-only molecular sampling.",
    "",
    "## Assumptions",
    sprintf("- Risk target alpha: %.3f (confidence %.1f%%)", cfg$alpha, 100 * (1 - cfg$alpha)),
    sprintf("- Design prevalence p0: %.4f", cfg$p0),
    sprintf("- Effective border sensitivity Se_eff = Se_mol * Se_subsample = %.3f", Se_eff(cfg)),
    "- Zero-acceptance policy: any positive rejects shipment.",
    "",
    "## Equations",
    "- Binomial approximation: n >= log(alpha) / log(1 - p0)",
    "- Hypergeometric zero probability: P0 = choose(N-D, n) / choose(N, n), D = ceil(p0*N)",
    "- Imperfect sensitivity approximation: P0 = (1 - p0 * Se_eff)^n",
    "",
    "## Tiered sampling policy",
    sprintf("- Small shipments: test all when N_asym <= %s", cfg$tiers$N_small_threshold),
    "- Medium/Large shipments: use hypergeometric and/or binomial floor with sensitivity inflation.",
    "",
    "## Key annual summary",
    paste(capture.output(print(annual_summary)), collapse = "\n"),
    "",
    "## Required n vs N_asym (head)",
    paste(capture.output(print(utils::head(required_tbl, 20))), collapse = "\n"),
    "",
    "## OC curve (head)",
    paste(capture.output(print(utils::head(oc_tbl, 20))), collapse = "\n"),
    "",
    "## Sensitivity analysis (head)",
    paste(capture.output(print(utils::head(sensitivity_tbl, 20))), collapse = "\n"),
    upstream_lines,
    "",
    "## Figures",
    "- ![](figures/required_n_vs_N.png)",
    "- ![](figures/oc_curve.png)",
    "- ![](figures/annual_risk_sensitivity.png)",
    "",
    "## Interpretation",
    "Shipment-level assurance is controlled by design n at p0 and alpha. Annual assurance depends on shipment volume and residual false-negative risk after upstream conditioning.",
    ""
  )
  readr::write_lines(lines, output_path)
  output_path
}
