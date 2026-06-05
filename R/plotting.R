#' Plotting helpers

#' Plot required n versus border population size.
#' @param required_tbl Required-n summary table.
#' @param out_path Output PNG path.
#' @return ggplot object.
plot_required_n_vs_N <- function(required_tbl, out_path) {
  p <- ggplot2::ggplot(required_tbl, ggplot2::aes(x = N_asym, y = required_n, color = scenario)) +
    ggplot2::geom_step(alpha = 0.9) +
    ggplot2::labs(x = "Asymptomatic units at border (N_asym)", y = "Required sample size (n)", color = "Scenario") +
    ggplot2::theme_minimal()
  ggplot2::ggsave(filename = out_path, plot = p, width = 8, height = 5, dpi = 120)
  p
}

#' Plot OC curves.
#' @param oc_tbl OC table from \\code{oc_curve_table()}.
#' @param out_path Output PNG path.
#' @return ggplot object.
plot_oc_curve <- function(oc_tbl, out_path) {
  p <- ggplot2::ggplot(oc_tbl, ggplot2::aes(x = n, y = confidence, color = as.factor(p))) +
    ggplot2::geom_line() +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(x = "Sample size (n)", y = "Achieved confidence", color = "Prevalence p") +
    ggplot2::theme_minimal()
  ggplot2::ggsave(filename = out_path, plot = p, width = 8, height = 5, dpi = 120)
  p
}

#' Plot annual risk sensitivity.
#' @param sensitivity_tbl Sensitivity results table.
#' @param out_path Output PNG path.
#' @return ggplot object.
plot_annual_risk_sensitivity <- function(sensitivity_tbl, out_path) {
  p <- ggplot2::ggplot(sensitivity_tbl, ggplot2::aes(x = value, y = prob_any_infected_pass_undetected, color = parameter)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    ggplot2::labs(x = "Parameter value", y = "Annual undetected infected pass probability", color = "Parameter") +
    ggplot2::theme_minimal()
  ggplot2::ggsave(filename = out_path, plot = p, width = 8, height = 5, dpi = 120)
  p
}
