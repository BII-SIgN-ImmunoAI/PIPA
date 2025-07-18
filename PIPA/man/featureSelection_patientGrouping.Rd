% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/featureSelection_patientGrouping.r
\name{featureSelection_patientGrouping}
\alias{featureSelection_patientGrouping}
\title{Feature selection - patient grouping and prognostic significance}
\usage{
featureSelection_patientGrouping(
  cell_feature_data = NULL,
  feature_data = NULL,
  surv_data = NULL,
  log10_featureNms = NULL,
  neighbor_size_byFrac = 0.05,
  plot_bnm = "",
  output_dir = getwd(),
  seed = 999,
  detection_rate_cutoff = seq(0, 0.9, by = 0.1),
  min_cluster_size = 20
)
}
\arguments{
\item{cell_feature_data}{A data.frame containing tumor-average cellular feature values: tumors in rows, features in columns}

\item{feature_data}{A data.frame containing the feature selection of each feature (rows) in each Cox PH Lasso run (columns)}

\item{surv_data}{A data.frame: tumors in rows (tumor_ids in row.names), with two colums "time" and "cens"}

\item{log10_featureNms}{A vector of characters of feature names whose values to be log10-transformed prior to scaling (z-score)}

\item{neighbor_size_byFrac}{A fraction (0 to 1) indicating the fraction of total samples/tumors to be
used in PhenoGraph}

\item{plot_bnm}{A character string for output plot filename identifier}

\item{output_dir}{A file.path string of output directory}

\item{seed}{An integer of seed (for repeatability)}

\item{detection_rate_cutoff}{A vector of fractions (0 to 1) indicating the (prognostically selected)
feature detection rates, to be tested}

\item{min_cluster_size}{An integer indicating the minimum cluster size (i.e. no. of cases/tumors) for which
a tumor subtype containing tumors less than this number will be excluded in the survival analysis
in \code{surv_analysis_fn} and in the corresponding output plot \code{Viz_heatmap_bySubtype}}
}
\description{
Feature selection - patient grouping and prognostic significance:
using cell features obtained at each 10% increment of detection rates
(determined by \code{select_prognostic_feature}), tumors are divided into
groups which are evaluated for prognostic significance.
}
